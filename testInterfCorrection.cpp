// c++ -o testInterfCorrection `root-config --glibs --cflags` -lm Functions.cc testInterfCorrection.cpp


#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cassert>

#include "TH1.h"
#include "TFile.h"
#include "TF1.h"
#include "TLorentzVector.h"
// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "LHAPDF/LHAPDF.h"
#include "Functions.h"

using namespace ROOT::Math ;
using namespace std ;

#include "findInterferenceTools.h"


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t crystalBallLowHighRatio (Double_t * xx, Double_t * par)
{
  return crystalBallLowHigh (xx, par) / crystalBallLowHigh (xx, par + 7) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main (int argc, char ** argv) 
{
  if (argc < 2) 
    {
      cout << "no mass" << endl ;
      exit (1) ;
    }

  double mass = atof (argv[1]) ;
  cout << "mass " << mass << endl ;

  //PG input files
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  string filename_NLOsig = "" ;
  double XS_NLOsig = -1 ;

  if (mass == 500)
    {
      cout << "mass " << 500 << endl ;
      filename_NLOsig = "/Users/govoni/data/lvjj_samples/interference/powheg/5936/h_qqH_WW_LNU2J_500_46_r1.lhe" ;
      XS_NLOsig = 0.1561 * 0.08 * 3 ; // pb 500 GeV  H_XS * X_BR (e+ or u+) * (charge and flavour factors: 3/2 * 2)
    }
  else if (mass == 800)
    {
      cout << "mass " << 800 << endl ;
      filename_NLOsig = "/Users/govoni/data/lvjj_samples/interference/powheg/6060/h_qqH_WW_LNU2J_800_54_r1.lhe" ;
      XS_NLOsig = 1. ; // pb 500 GeV
    }
  else
    {
      cout << "wrong mass: " << mass << endl ;
      exit (1) ;
    }

  std::ifstream ifs_NLOsig (filename_NLOsig.c_str ()) ;
  LHEF::Reader reader_NLOsig (ifs_NLOsig) ;

  //PG output histograms
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TH1F * h_MWW     = new TH1F (TString ("h_MWW"), 
                               TString ("h_MWW"), 200., 200., 1800.) ;
                   
  TH1F * h_MWW_VBF = new TH1F (TString ("h_MWW_VBF"), 
                               TString ("h_MWW_VBF"), 200., 200., 1800.) ;

  //PG loop over events
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  double totalWeight = 0 ;

  //PG loop over input events
  while (reader_NLOsig.readEvent ()) 
    {
      if ( reader_NLOsig.outsideBlock.length() ) std::cout << reader_NLOsig.outsideBlock;
    
      vector<TLorentzVector> v_f_Ws ;
      vector<TLorentzVector> v_f_quarks ;
      vector<TLorentzVector> v_f_leptons ;
      vector<TLorentzVector> v_f_neutrinos ;
      vector<TLorentzVector> v_f_higgs ;
      
      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader_NLOsig.hepeup.IDUP.size (); ++iPart) 
        {
          TLorentzVector dummy = buildP (reader_NLOsig.hepeup, iPart) ;

          // intermediate particles          
          if (reader_NLOsig.hepeup.ISTUP.at (iPart) == 2)
            {
              if (abs (reader_NLOsig.hepeup.IDUP.at (iPart)) == 24) 
                {
                  v_f_Ws.push_back (dummy) ;
                }              
            } // intermediate particles
            
          // outgoing particles          
          if (reader_NLOsig.hepeup.ISTUP.at (iPart) == 1)
            {
           // quarks
           if (abs (reader_NLOsig.hepeup.IDUP.at (iPart)) < 7) 
             {
               v_f_quarks.push_back (dummy) ;        
             } // quarks
           else if (abs (reader_NLOsig.hepeup.IDUP.at (iPart)) == 11 ||
                    abs (reader_NLOsig.hepeup.IDUP.at (iPart)) == 13 ||
                    abs (reader_NLOsig.hepeup.IDUP.at (iPart)) == 15)
             {
               v_f_leptons.push_back (dummy) ;
             }
           else if (abs (reader_NLOsig.hepeup.IDUP.at (iPart)) == 12 ||
                    abs (reader_NLOsig.hepeup.IDUP.at (iPart)) == 14 ||
                    abs (reader_NLOsig.hepeup.IDUP.at (iPart)) == 16)
             {
               v_f_neutrinos.push_back (dummy) ;        
             }
           if (abs (reader_NLOsig.hepeup.IDUP.at (iPart)) == 25)
             {
               v_f_higgs.push_back (dummy) ;                       
             } 

         } // outgoing particles
       } // loop over particles in the event

//      if (totalCount < 10) cout << "PARTICLES " <<  v_f_leptons.size () << "\t" << v_f_neutrinos.size () << "\t" << v_f_quarks.size () << "\n" ;

//      pair<int, int> detaIndices = findPairWithLargestDeta (v_f_quarks) ;
//      if (v_f_quarks.at (detaIndices.second).Eta () - v_f_quarks.at (detaIndices.first).Eta () < 2) continue ;
//      TLorentzVector largestPair = v_f_quarks.at (detaIndices.second) + v_f_quarks.at (detaIndices.first) ;
//      if (largestPair.M () < 100) continue ; //PG selection applied in phantom

      //PG the first two are the VBF jets, the following ones the W jets
      sort (v_f_quarks.rbegin (), v_f_quarks.rend (), ptsort ()) ;  
      if (v_f_quarks.size () < 2) continue ;
      
      double mass = v_f_higgs.at (0).M () ;

      h_MWW->Fill (mass) ;

      //PG VBF cuts from Luca's analysis
      if (fabs (v_f_quarks.at (0).Eta () - v_f_quarks.at (1).Eta ()) < 2.5) continue ;
      TLorentzVector dijet = v_f_quarks.at (0) + v_f_quarks.at (1) ;
      if (dijet.M () < 250) continue ;

      h_MWW_VBF->Fill (mass) ;

    } //PG loop over input events

  //PG correction factor
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TF1 * f_ratio = new TF1 ("f_ratio", crystalBallLowHighRatio, 0, 2000, 14);
  f_ratio->SetLineWidth (3) ;
  f_ratio->SetLineColor (kGray + 2) ;
  double f_ratio_params[14] = {0.000673659, 502.995, -30.331, 1.3619, 2.60642, 0.833161, 2.59464,        // phantom
                               0.000698824, 499.888, -22.7945, -0.511835, 47.9117, -0.827871, 3.36354} ; // madgraph
  f_ratio->SetParameters (f_ratio_params) ;

  TH1F * h_MWW_corr     = (TH1F *) h_MWW->Clone ("h_MWW_corr") ;
  for (int i = 1 ; i < h_MWW_corr->GetNbinsX () ; ++i)
    {
      double CFactor = f_ratio->Eval (h_MWW_corr->GetBinCenter (i)) ;
      h_MWW_corr->SetBinContent (i, CFactor * h_MWW_corr->GetBinContent (i)) ;
    }

  TH1F * h_MWW_corr_VBF = (TH1F *) h_MWW_VBF->Clone ("h_MWW_corr_VBF") ;
  for (int i = 1 ; i < h_MWW_corr_VBF->GetNbinsX () ; ++i)
    {
      double CFactor = f_ratio->Eval (h_MWW_corr_VBF->GetBinCenter (i)) ;
      h_MWW_corr_VBF->SetBinContent (i, CFactor * h_MWW_corr_VBF->GetBinContent (i)) ;
    }

  TFile f ("testInterfCorrection.root", "recreate") ;
  h_MWW->Write () ;
  h_MWW_corr->Write () ;
  h_MWW_VBF->Write () ;
  h_MWW_corr_VBF->Write () ;
  f.Close () ;

  return 0 ;
}

/*
dove stanno i LHE file di segnale:
500 5936
800 6060

http://cms.cern.ch/iCMS/jsp/mcprod/admin/requestmanagement.jsp?dsn=VBF*HToWWToLAndTauNuQQ_M*&campid=Summer12

*/

