// c++ -o checkMomentum_00 `root-config --glibs --cflags` -lm checkMomentum_00.cpp
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "TH1.h"
#include "TFile.h"

#include "TLorentzVector.h"

// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;
using namespace std ;

typedef LorentzVector<ROOT::Math::PxPyPzE4D<double> > lorentzVector ;


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


int main(int argc, char ** argv) 
{
  if(argc < 2)
    {
      cout << "Usage:   " << argv[0] 
           << " input.lhe " << endl ;
      return -1;
    }

  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  TH1F mass ("mass", "mass", 100, -1, 3) ;

  float k2 = 0.000510998928 * 0.000510998928 - 1.77682 * 1.77682 ; // GeV -3.15708905128
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
           // outgoing particles          
           if (reader.hepeup.ISTUP.at (iPart) == 1)
             {
               if (abs (reader.hepeup.IDUP.at (iPart)) == 15) 
                 {     
                   lorentzVector dummy_new = buildLP (reader.hepeup, iPart) ;
                   TLorentzVector dummy = buildP (reader.hepeup, iPart) ;
                   float p2 = dummy.Vect ().Mag2 () ;
                   float scale = sqrt (1 + k2 / p2) ;
                   reader.hepeup.PUP.at (iPart).at (0) *= scale ; // px
                   reader.hepeup.PUP.at (iPart).at (1) *= scale ; // px
                   reader.hepeup.PUP.at (iPart).at (2) *= scale ; // px
                   TLorentzVector dummy3 = buildP (reader.hepeup, iPart) ;
//                   cout << reader.hepeup.PUP.at (iPart).at (4)
//                        << "\t" << dummy.Mag () * dummy.Mag () / dummy_new.M2 ()
//                        << "\t" << dummy.Mag () * dummy.Mag () 
//                        << "\t" << dummy.Mag () 
//                        << "\t" << dummy_new.M2 ()
//                        << "\t" << sqrt (dummy_new.M2 ()) << "\n" ;
                   mass.Fill (dummy.Mag ()) ;
//                   << "\t(" << scale << ")\t" << dummy3.Mag () << endl ;
//                   cout << dummy.Vect ().Mag () 
//                        << "\t(" << scale << ")\t" << dummy3.Vect ().Mag () 
//                        << "\t" << (dummy3.Vect ().Mag () / dummy.Vect ().Mag ())
//                        << endl ;
//                   cout << "----\n" ;
                 }
             } // outgoing particles
        } // loop over particles in the event

    } //PG loop over input events

  TFile f ("output_checkP_00.root", "recreate") ;
  mass.Write () ;
  f.Close () ;

  return 0 ;
}
