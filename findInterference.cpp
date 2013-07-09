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


void fillWWhist (LHEF::Reader & reader, TH1F * h_MWW)
{
//      TLorentzVector Higgs;
//      int iPartHiggs;
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

      //PG reconstruct the final state invariant mass system
      //PG assign the two jets      
      cout << "leptons : " << v_f_leptons.size () << "\t" ;
      cout << "quarks : "  << v_f_quarks.size () << "\t" ;
      cout << "neutrinos : "  << v_f_neutrinos.size () << "\t" ;
      cout << endl ;
   
      h_MWW->Fill (total.M ()) ;
    } //PG loop over input events

}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main(int argc, char ** argv) 
{
  string filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/madgraph_500GeV_4jlv.lhe" ;
  double XS_mg = 0.009129 ; // pb
  
  std::ifstream ifs_mg (filename_mg.c_str ()) ;
  LHEF::Reader reader_mg (ifs_mg) ;
  TH1F * h_mWW_mg = new TH1F ("h_mWW_mg", "h_mWW_mg", 1000, 0, 3000) ;
  fillWWhist (reader_mg, h_mWW_mg) ;


//  TFile f ("findInterference.root", "recreate") ;
//  h_mWW.Write () ;
//  f.Close () ;

  return 0 ;
}
