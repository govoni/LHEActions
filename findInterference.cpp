// c++ -o findInterference `root-config --glibs --cflags` `lhapdf-config --cppflags  --ldflags` -lm findInterference.cpp

#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>

#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"
// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

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


int fillWWhist (LHEF::Reader & reader, TH1F * h_MWW, double XS)
{
  int totalCount = 0 ;
   
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
    
      vector<TLorentzVector> v_f_Ws ;
      vector<TLorentzVector> v_f_quarks ;
      vector<TLorentzVector> v_f_leptons ;
      vector<TLorentzVector> v_f_neutrinos ;
      
      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
          TLorentzVector dummy = buildP (reader.hepeup, iPart) ;
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

      //PG the first two are the VBF jets, the following ones the W jets
      sort (v_f_quarks.rbegin (), v_f_quarks.rend (), ptsort ()) ;  
      TLorentzVector total = (v_f_leptons.at (0) + v_f_neutrinos.at (0)) + (v_f_quarks.at (2) + v_f_quarks.at (3)) ;

      h_MWW->Fill (total.M ()) ;
      ++totalCount ;
    } //PG loop over input events

  h_MWW->Scale (XS / totalCount) ;

  return totalCount ;
  
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main(int argc, char ** argv) 
{
  //PG ---- madgraph ---- signal only
  
  string filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/madgraph_500GeV_4jlv.lhe" ;
  double XS_mg = 0.009129 ; // pb
  
  std::ifstream ifs_mg (filename_mg.c_str ()) ;
  LHEF::Reader reader_mg (ifs_mg) ;
  TH1F * h_MWW_mg = new TH1F ("h_MWW_mg", "h_MWW_mg", 1000, 0, 3000) ;
  int entries_mg = fillWWhist (reader_mg, h_MWW_mg, XS_mg) ;

  //PG ---- phantom ---- background only

  string filename_phbkg = "/Users/govoni/data/lvjj_samples/interference/4jlv_h126/genh126/total.lhe" ;
  double XS_phbkg = 0.07762748 ; // 7.76274847686845293E-002 // pb
  
  std::ifstream ifs_phbkg (filename_phbkg.c_str ()) ;
  LHEF::Reader reader_phbkg (ifs_phbkg) ;
  TH1F * h_MWW_phbkg = new TH1F ("h_MWW_phbkg", "h_MWW_phbkg", 1000, 0, 3000) ;
  int entries_phbkg = fillWWhist (reader_phbkg, h_MWW_phbkg, XS_phbkg) ;

  //PG ---- phantom ---- background and signal

  string filename_phbkgsig = "/Users/govoni/data/lvjj_samples/interference/4jlv/genh500/total.lhe" ;
  double XS_phbkgsig = 0.078904216 ; // 7.890421624985394E-002 // pb
  
  std::ifstream ifs_phbkgsig (filename_phbkgsig.c_str ()) ;
  LHEF::Reader reader_phbkgsig (ifs_phbkgsig) ;
  TH1F * h_MWW_phbkgsig = new TH1F ("h_MWW_phbkgsig", "h_MWW_phbkgsig", 1000, 0, 3000) ;
  int entries_phbkgsig = fillWWhist (reader_phbkgsig, h_MWW_phbkgsig, XS_phbkgsig) ;

  //PG operations
  
  TH1F * h_subtr = (TH1F *) h_MWW_phbkgsig->Clone ("h_subtr") ;
  h_subtr->Add (h_MWW_mg, -1) ;

  TFile f ("findInterference.root", "recreate") ;
  h_MWW_mg->Write () ;
  h_MWW_phbkg->Write () ;
  h_MWW_phbkgsig->Write () ;
  h_subtr->Write () ;
  f.Close () ;

  return 0 ;
}
