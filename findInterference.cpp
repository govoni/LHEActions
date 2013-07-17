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

typedef LorentzVector<ROOT::Math::PxPyPzE4D<double> > lorentzVector ;



struct ptsort: public std::binary_function<TLorentzVector, TLorentzVector, bool>
{
  bool operator () (const TLorentzVector & x, const TLorentzVector & y)
    {
      return  (x.Perp () < y.Perp () ) ;
    }
} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----




struct etasort: public std::binary_function<TLorentzVector, TLorentzVector, bool>
{
  bool operator () (const TLorentzVector & x, const TLorentzVector & y)
    {
      return  (x.Eta () < y.Eta () ) ;
    }
} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


pair<int, int> findPairWithLargestDeta (const vector<TLorentzVector> & v_f_quarks)
{
  int minimum = min_element (v_f_quarks.begin (), v_f_quarks.end (), etasort ()) - v_f_quarks.begin () ;
  int maximum = max_element (v_f_quarks.begin (), v_f_quarks.end (), etasort ()) - v_f_quarks.begin () ;
  return pair<int, int> (minimum, maximum) ;
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int splitString (vector<string>& fields, const string & work, char delim, int rep = 0) {
    if (!fields.empty ()) fields.clear () ;  // empty vector if necessary
    string buf = "" ;
    int i = 0;
    while (i < work.length ()) {
        if (work[i] != delim)
            buf += work[i] ;
        else if (rep == 1) {
            fields.push_back (buf) ;
            buf = "";
        } else if (buf.length() > 0) {
            fields.push_back (buf) ;
            buf = "" ;
        }
        i++;
    }
    if (!buf.empty ())
        fields.push_back (buf) ;
    return fields.size () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


TLorentzVector buildP (const LHEF::HEPEUP & event, int iPart)
{
  TLorentzVector dummy ;
  dummy.SetPxPyPzE (
      event.PUP.at (iPart).at (0), // px
      event.PUP.at (iPart).at (1), // py
      event.PUP.at (iPart).at (2), // pz
      event.PUP.at (iPart).at (3) // E
    ) ;
  return dummy ;  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


lorentzVector buildLP (const LHEF::HEPEUP & event, int iPart)
{
  lorentzVector dummy ;
  dummy.SetPxPyPzE (
      event.PUP.at (iPart).at (0), // px
      event.PUP.at (iPart).at (1), // py
      event.PUP.at (iPart).at (2), // pz
      event.PUP.at (iPart).at (3) // E
    ) ;
  return dummy ;  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


struct histos
{
  TH1F *  m_h_MWW ;
  TH1F *  m_h_scale ;
  TString m_name ;
  double  m_XS ;
  
  histos (TString name, double XS) : m_name (name), m_XS (XS)
    {
      m_h_MWW = new TH1F (TString ("h_MWW_") + name, 
                          TString ("h_MWW_") + name, 50, 200., 2000.) ;
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
fillHistos (LHEF::Reader & reader, histos & Histos, double XS, double referenceScale = 0)
{
  double totalCount = 0. ;
   
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    
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

             } // incoming particle          

          // intermediate particles          
          if (reader.hepeup.ISTUP.at (iPart) == 2)
            {
              if (abs (reader.hepeup.IDUP.at (iPart)) == 24) 
                {
                  v_f_Ws.push_back (dummy) ;
                }              
            }
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

      pair<int, int> detaIndices = findPairWithLargestDeta (v_f_quarks) ;
      if (v_f_quarks.at (detaIndices.second).Eta () - v_f_quarks.at (detaIndices.first).Eta () < 2) continue ;
      TLorentzVector largestPair = v_f_quarks.at (detaIndices.second) + v_f_quarks.at (detaIndices.first) ;
      if (largestPair.M () < 100) continue ; //PG selection applied in phantom

      int cont = 0 ;
      for (int iJ = 0 ; iJ < 4 ; ++iJ)
        for (int iJ2 = iJ + 1 ; iJ2 < 4 ; ++iJ2)
          if (v_f_quarks.at (iJ).DeltaR (v_f_quarks.at (iJ2)) < 0.4) cont = 1 ;
      if (cont == 1) continue ;

      //PG the first two are the VBF jets, the following ones the W jets
      sort (v_f_quarks.rbegin (), v_f_quarks.rend (), ptsort ()) ;  
      TLorentzVector total = (v_f_leptons.at (0) + v_f_neutrinos.at (0)) + (v_f_quarks.at (2) + v_f_quarks.at (3)) ;

      //PG the scale:
      float scale = reader.hepeup.SCALUP ;
      Histos.m_h_scale->Fill (scale) ;

      double weight = 1. ;
      if (referenceScale != 0 )
        weight = LHAPDF::xfx (x[0], referenceScale, flavour[0]) * LHAPDF::xfx (x[1], referenceScale, flavour[1]) /
                 (LHAPDF::xfx (x[0], scale, flavour[0]) * LHAPDF::xfx (x[1], scale, flavour[1])) ;

      Histos.m_h_MWW->Fill (total.M ()) ;
      totalCount += weight ;

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


  //PG ---- madgraph ---- signal only
  
  string filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/madgraph_500GeV_4jlv.lhe" ;
  double XS_mg = 0.009129 ; // pb
  
  std::ifstream ifs_mg (filename_mg.c_str ()) ;
  LHEF::Reader reader_mg (ifs_mg) ;
  histos H_mg ("mg", XS_mg) ;
  double entries_mg = fillHistos (reader_mg, H_mg, XS_mg, 500.) ;

  cout << "madgraph events : " << entries_mg << endl ;
  
  //PG ---- phantom ---- background only

  string filename_phbkg = "/Users/govoni/data/lvjj_samples/interference/4jlv_h126/genh126/total.lhe" ;
  double XS_phbkg = 0.07762748 * 2 ; // 7.76274847686845293E-002 // pb the factor 2 accounts for muons, electrons
  
  std::ifstream ifs_phbkg (filename_phbkg.c_str ()) ;
  LHEF::Reader reader_phbkg (ifs_phbkg) ;
  histos H_phbkg ("phbkg", XS_phbkg) ;
  double entries_phbkg = fillHistos (reader_phbkg, H_phbkg, XS_phbkg, 500.) ;

  cout << "phantom bkg events : " << entries_phbkg << endl ;

  //PG ---- phantom ---- background and signal

  string filename_phbkgsig = "/Users/govoni/data/lvjj_samples/interference/4jlv/genh500/total.lhe" ;
//  string filename_phbkgsig = "/Users/govoni/data/lvjj_samples/interference/4jlv/genh800/total.lhe" ;
  double XS_phbkgsig = 0.078904216 * 2 ; // 7.890421624985394E-002 // pb
  
  std::ifstream ifs_phbkgsig (filename_phbkgsig.c_str ()) ;
  LHEF::Reader reader_phbkgsig (ifs_phbkgsig) ;
  histos H_phbkgsig ("phbkgsig", XS_phbkgsig) ;
  double entries_phbkgsig = fillHistos (reader_phbkgsig, H_phbkgsig, XS_phbkgsig, 500.) ;

  cout << "phantom bkg+sig events : " << entries_phbkgsig << endl ;

  //PG saving the histograms

  TFile f ("findInterference.root", "recreate") ;
  H_phbkgsig.save (f) ;
  H_phbkg.save (f) ;
  H_mg.save (f) ;
  f.Close () ;

  return 0 ;
}

/*

TCanvas c1
c1.DrawFrame (100,0.00001,2000,0.07)
h_MWW_phbkg->SetStats (0)
h_MWW_phbkgsig->SetStats (0)
h_MWW_mg->SetStats (0)
h_MWW_phbkg->SetLineColor (kOrange)
h_MWW_phbkg->SetLineWidth (2)
h_MWW_phbkg->Draw ("histsame")
h_MWW_phbkgsig->SetLineColor (kRed)
h_MWW_phbkgsig->SetLineWidth (2)
h_MWW_phbkgsig->Draw ("histsame")

TH1F * diff = (TH1F *) h_MWW_phbkgsig->Clone ("diff")
diff->SetTitle ("")
diff->Add (h_MWW_phbkg, -1) 

TH1F * ratio = (TH1F *) h_MWW_mg->Clone ("ratio") 
ratio->SetTitle ("")
ratio->Divide (diff)
cout << "scaling by " << 1. / ratio->GetBinContent (ratio->FindBin (500)) << endl ;
//h_MWW_mg->Scale (1. / ratio->GetBinContent (ratio->FindBin (500)))

h_MWW_mg->Draw ("histsame")

TCanvas c2
diff->Draw ("hist")

TCanvas c3
ratio->Draw ("hist")

TCanvas c4
diff->Draw ("hist")
h_MWW_mg->Draw ("histsame")


*/
