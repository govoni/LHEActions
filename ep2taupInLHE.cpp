// c++ -o ep2taupInLHE `root-config --glibs --cflags` -lm ep2taupInLHE.cpp


#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include "TLorentzVector.h"

int main(int argc, char ** argv) {
  if(argc < 3) {
   std::cout << "Usage:   " << argv[0] << " input.lhe output.lhe" << std::endl ;
   return -1;
  }

  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  std::ofstream outputStream (argv[2]) ;
  LHEF::Writer writer (outputStream) ;

  writer.headerBlock() << reader.headerBlock ;
  writer.initComments() << reader.initComments ;
  writer.heprup = reader.heprup ;
  writer.init () ;

//  float k2 = 0.000510998928 * 0.000510998928 - 1.77682 * 1.77682 ; // GeV -3.15708905128
  //PG ele massless in phantom
  float k2 = 0. - 1.77682 * 1.77682 ; // GeV -3.15708905128
  //PG loop over input events
  while (reader.readEvent ()) {
   if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

      // loop over particles in the event
   for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) {
           // outgoing particles          
    if (reader.hepeup.ISTUP.at (iPart) == 1) {
     if ( (reader.hepeup.IDUP.at (iPart)) == -11) {    //  e+
      TLorentzVector dummy (
        reader.hepeup.PUP.at (iPart).at (0), // px
        reader.hepeup.PUP.at (iPart).at (1), // py
        reader.hepeup.PUP.at (iPart).at (2), // pz
        reader.hepeup.PUP.at (iPart).at (3) // E
      ) ;
      float p2 = dummy.Vect ().Mag2 () ;
      float scale = sqrt (1 + k2 / p2) ;
      reader.hepeup.PUP.at (iPart).at (0) *= scale ; // px
      reader.hepeup.PUP.at (iPart).at (1) *= scale ; // px
      reader.hepeup.PUP.at (iPart).at (2) *= scale ; // px

      reader.hepeup.IDUP.at (iPart) = -15 ;    //  tau+
     }
     if (reader.hepeup.IDUP.at (iPart) == 12) reader.hepeup.IDUP.at (iPart) = 16 ;
    } // outgoing particles
   } // loop over particles in the event
   writer.eventComments() << reader.eventComments;
   writer.hepeup = reader.hepeup;
   writer.writeEvent();

  } //PG loop over input events

  return 0 ;
}





