// c++ -o findInterference `root-config --glibs --cflags` `lhapdf-config --cppflags  --ldflags` -lm findInterference.cpp

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
#include "TLorentzVector.h"
// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "LHAPDF/LHAPDF.h"


using namespace ROOT::Math ;
using namespace std ;

#include "findInterferenceTools.h"


struct histos
{
  TH1F *  m_h_MWW ;
  TH1F *  m_h_scale ;
  TString m_name ;
  double  m_XS ;
  
  histos (TString name, double XS) : m_name (name), m_XS (XS)
    {
      m_h_MWW = new TH1F (TString ("h_MWW_") + name, 
                          TString ("h_MWW_") + name, 200., 200., 1800.) ;
      m_h_MWW->Sumw2 () ;
      m_h_scale = new TH1F (TString ("h_scale_") + name, 
                          TString ("h_scale_") + name, 100, 0., 1000.) ;
      m_h_scale->Sumw2 () ;
    }
 
  void norm (double total = 0)
    {
      double factor = m_XS / m_h_MWW->GetEntries () ;
      if (total != 0) factor = m_XS / total ;
      m_h_MWW->Scale (factor) ;
      m_h_scale->Scale (factor) ;
    }
  
  ~histos ()
    {
      delete m_h_MWW ;
      delete m_h_scale ;
    }  
    
  void save (TFile & outfile) 
    {
      outfile.cd () ;
      m_h_MWW->Write () ;
      m_h_scale->Write () ;
    }
  

} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


double 
fillHistos (LHEF::Reader & reader, histos & Histos, double XS, double referenceScale = 0, int max = -1)
{
  double totalCount = 0. ;
  int events = 0 ;
   
  //PG loop over input events
  while (reader.readEvent ()) 
    {
//      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    
      vector<TLorentzVector> v_f_Ws ;
      vector<TLorentzVector> v_f_quarks ;
      vector<TLorentzVector> v_f_leptons ;
      vector<TLorentzVector> v_f_neutrinos ;
      
      double x[2] = {0., 0.} ;
      int flavour[2] = {0, 0} ;

      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
          TLorentzVector dummy = buildP (reader.hepeup, iPart) ;

           // incoming particles        
           if (reader.hepeup.ISTUP.at (iPart) == -1) 
             {
               x[iPart] = dummy.P () / 4000. ;
               flavour[iPart] = reader.hepeup.IDUP.at (iPart) ;
             } // incoming particles         

          // intermediate particles          
          if (reader.hepeup.ISTUP.at (iPart) == 2)
            {
              if (abs (reader.hepeup.IDUP.at (iPart)) == 24) 
                {
                  v_f_Ws.push_back (dummy) ;
                }              
            } // intermediate particles
            
          // outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
           // quarks
           if (abs (reader.hepeup.IDUP.at (iPart)) < 7) 
             {
               v_f_quarks.push_back (dummy) ;        
             } // quarks
           else if (abs (reader.hepeup.IDUP.at (iPart)) == 11 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 13 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 15)
             {
               v_f_leptons.push_back (dummy) ;
             }
           else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 14 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 16)
             {
               v_f_neutrinos.push_back (dummy) ;        
             }
         } // outgoing particles
       } // loop over particles in the event

//      if (totalCount < 10) cout << "PARTICLES " <<  v_f_leptons.size () << "\t" << v_f_neutrinos.size () << "\t" << v_f_quarks.size () << "\n" ;

      int warnNum = 0 ;
      if (v_f_quarks.size () < 4)
        {
          cout << "warning, not enough quarks" << endl ;
          ++warnNum ;
        }
      if (v_f_leptons.size () < 1)
        {
          cout << "warning, not enough leptons" << endl ;
          ++warnNum ;
        }
      if (v_f_neutrinos.size () < 1)
        {
          cout << "warning, not enough neutrinos" << endl ;
          ++warnNum ;
        }
      if (warnNum > 0) continue ;


      pair<int, int> detaIndices = findPairWithLargestDeta (v_f_quarks) ;
      if (v_f_quarks.at (detaIndices.second).Eta () - v_f_quarks.at (detaIndices.first).Eta () < 2) continue ;
      TLorentzVector largestPair = v_f_quarks.at (detaIndices.second) + v_f_quarks.at (detaIndices.first) ;
      if (largestPair.M () < 100) continue ; //PG selection applied in phantom

      //PG do I need this cut?! FIXME
      int cont = 0 ;
      for (int iJ = 0 ; iJ < 4 ; ++iJ)
        for (int iJ2 = iJ + 1 ; iJ2 < 4 ; ++iJ2)
          if (v_f_quarks.at (iJ).DeltaR (v_f_quarks.at (iJ2)) < 0.4) cont = 1 ;
      if (cont == 1) continue ;

      //PG the first two are the VBF jets, the following ones the W jets
      sort (v_f_quarks.rbegin (), v_f_quarks.rend (), ptsort ()) ;  
      
//      pair<int, int> Wpair (2, 3) ;
      pair<int, int> Wpair = findPairWithWMass (v_f_quarks) ;

      if (Wpair.first > 3 || Wpair.second > 3)
        {
          cout << "warning, wrong quarks in W determination\n" ;
        }

      TLorentzVector total = (v_f_leptons.at (0) + v_f_neutrinos.at (0)) + 
                             (v_f_quarks.at (Wpair.first) + v_f_quarks.at (Wpair.second)) ;

      //PG the scale:
      float scale = reader.hepeup.SCALUP ;
      Histos.m_h_scale->Fill (scale) ;

      double weight = 1. ;
      if (referenceScale != 0 )
        weight = LHAPDF::xfx (x[0], referenceScale, flavour[0]) * LHAPDF::xfx (x[1], referenceScale, flavour[1]) /
                 (LHAPDF::xfx (x[0], scale, flavour[0]) * LHAPDF::xfx (x[1], scale, flavour[1])) ;

      Histos.m_h_MWW->Fill (total.M (), weight) ;
      totalCount += weight ;
      ++events ;
      if (max > 0 && max < events) 
        {
          cout << max << " events reached, exiting" << endl ;
          break ;
        }

    } //PG loop over input events

  Histos.norm (totalCount) ;

  return totalCount ;
  
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main (int argc, char ** argv) 
{

  const int SUBSET = 0 ;
  const string NAME = "cteq6ll" ; //"cteq6l1"

  LHAPDF::initPDFSet (NAME, LHAPDF::LHPDF, SUBSET) ;
  const int NUMBER = LHAPDF::numberPDF () ;

  LHAPDF::initPDF (0) ;

  if (argc < 2) 
    {
      cout << "mass is missing" << endl ;
      exit (1) ;
    }

  int maxEventsPerSample = -1 ;
  if (argc == 3)
    {
      int dummy = atoi (argv[2]) ;
      if (dummy > 0) maxEventsPerSample = dummy ;
    }
    
  double mass = atof (argv[1]) ;
  if (mass != 350 && mass != 500 && mass != 650 && mass != 800 && mass != 1000)
    {
      cout << "wrong mass: " << mass << endl ;
      exit (1) ;
    }

  //PG choose the samples
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

//  // Sandro's sample
//  string filename_phbkg = "/Users/govoni/data/lvjj_samples/interference/4jlv_h126/genh126/total.lhe" ;
//  double XS_phbkg = 0.07762748 * 2 ; // 7.76274847686845293E-002 // pb the factor 2 accounts for muons, electrons

  // pietro's sample
  string filename_phbkg = "/Users/govoni/data/lvjj_samples/interference/phantom/total.126.lhe" ;
  double XS_phbkg = 0.07756069 * 2 ; // 7.7560687011E-002 // pb the factor 2 accounts for muons, electrons
                    
  string filename_mg = "" ;
  double XS_mg = -1 ;
  string filename_phbkgsig = "" ;
  double XS_phbkgsig = -1 ;

  if (mass == 350)
    {
      filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/madgraph_qqHWWlv4j_mh350.lhe" ;
      XS_mg = 0.015413 ; // pb 350 GeV
    
      filename_phbkgsig = "/Users/govoni/data/lvjj_samples/interference/phantom/total.350.lhe" ;
      XS_phbkgsig =  0.087460 * 2 ; // 8.7459661306E-002 // pb 350 GeV
    }
  else if (mass == 500)
    {
      filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/madgraph_500GeV_4jlv.lhe" ;
      XS_mg = 0.009129 ; // pb 500 GeV
    
      filename_phbkgsig = "/Users/govoni/data/lvjj_samples/interference/phantom/" ;
//      XS_phbkgsig = 0.078904216 * 2 ; // 7.890421624985394E-002 // pb 500 GeV // sample from Sandro
      XS_phbkgsig = 0.078719  * 2 ; // 7.8718631364E-002 // pb 500 GeV // sample from Pietro
    }
  else if (mass == 650)
    {
      filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/madgraph_qqHWWlv4j_mh650.lhe" ;
      XS_mg = 0.019246 ; // pb 650 GeV
      
      filename_phbkgsig = "/Users/govoni/data/lvjj_samples/interference/phantom/total.650.lhe" ;
      XS_phbkgsig = 0.076244 * 2 ; // 7.62444816705E-002 // pb 650 GeV
    }
  else if (mass == 800)
    {
      filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/H800_lvl4jets.lhe" ;
      XS_mg = 0.0014136 ; // pb 800 GeV
      
      // XS_phbkgsig = 0.075067956 * 2 ; // 7.506795619825214E-002 // pb 800 GeV // sample from Sandro
      filename_phbkgsig = "/Users/govoni/data/lvjj_samples/interference/phantom/total.800.lhe" ;
      XS_phbkgsig = 0.075085 * 2 ; // 7.5085416255E-002 // pb 800 GeV // sample from Pietro
    }
  else if (mass == 1000)
    {
      filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/madgraph_qqHWWlv4j_mh1000.lhe" ;
      XS_mg = 0.0016269 ; // pb 1000 GeV
      
      filename_phbkgsig = "/Users/govoni/data/lvjj_samples/interference/phantom/total.1000.lhe" ;
      XS_phbkgsig =  0.074253 * 2 ; // 7.4253338203E-002 // pb 1000 GeV
    }

  //PG messages
  
  cout << "\nworking with mass : " << mass << endl ;
  cout << "signal  :\t" << filename_mg << endl ;
  cout << "bkg     :\t" << filename_phbkg << endl ;
  cout << "S + bkg :\t" << filename_phbkgsig << "\n\n" ;

  //PG ---- madgraph ---- signal only
  
  std::ifstream ifs_mg (filename_mg.c_str ()) ;
  LHEF::Reader reader_mg (ifs_mg) ;
  histos H_mg ("mg", XS_mg) ;
  double entries_mg = fillHistos (reader_mg, H_mg, XS_mg, mass, maxEventsPerSample) ;

  cout << "madgraph events : " << entries_mg << endl ;
  
  //PG ---- phantom ---- background only

  std::ifstream ifs_phbkg (filename_phbkg.c_str ()) ;
  LHEF::Reader reader_phbkg (ifs_phbkg) ;
  histos H_phbkg ("phbkg", XS_phbkg) ;
  double entries_phbkg = fillHistos (reader_phbkg, H_phbkg, XS_phbkg, mass, maxEventsPerSample) ;

  cout << "phantom bkg events : " << entries_phbkg << endl ;

  //PG ---- phantom ---- background and signal

  std::ifstream ifs_phbkgsig (filename_phbkgsig.c_str ()) ;
  LHEF::Reader reader_phbkgsig (ifs_phbkgsig) ;
  histos H_phbkgsig ("phbkgsig", XS_phbkgsig) ;
  double entries_phbkgsig = fillHistos (reader_phbkgsig, H_phbkgsig, XS_phbkgsig, mass, maxEventsPerSample) ;

  cout << "phantom bkg+sig events : " << entries_phbkgsig << endl ;

  //PG saving the histograms

  TString name = "findInterference." ;
  name += mass ; 
  name += ".root" ;
  TFile f (name, "recreate") ;
  H_phbkgsig.save (f) ;
  H_phbkg.save (f) ;
  H_mg.save (f) ;
  f.Close () ;

  cout << "madgraph signal (" << XS_mg << " pb):\n" << filename_mg << "\n" ;
  cout << "phantom signal + bkg (" << XS_phbkgsig << " pb):\n" << filename_phbkgsig << "\n" ;
  cout << "phantom bkg (" << XS_phbkg << " pb):\n" << filename_phbkg << "\n" ;
  cout << "\n produced " << name << endl ;


  return 0 ;
}
