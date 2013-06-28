// c++ -o muToTauInLHE `root-config --glibs --cflags` -lm muToTauInLHE.cpp
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include "TLorentzVector.h"

using namespace std ;


int main(int argc, char ** argv) 
{
  if(argc < 3)
    {
      cout << "Usage:   " << argv[0] 
           << " input.lhe output.lhe" << endl ;
      return -1;
    }

  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  ofstream outputStream (argv[2]) ;
  LHEF::Writer writer (outputStream) ;

  writer.headerBlock() << reader.headerBlock ;
  writer.initComments() << reader.initComments ;
  writer.heprup = reader.heprup ;
  writer.init () ;

//PG mu massless in phantom
//  float k2 = 0.1056583715 * 0.1056583715 - 1.77682 * 1.77682 ; // GeV -3.14592562093
  float k2 = 0. - 1.77682 * 1.77682 ; // GeV -3.14592562093
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
               if (abs (reader.hepeup.IDUP.at (iPart)) == 13)
                 {
                   TLorentzVector dummy
                     (
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
    
                   if (reader.hepeup.IDUP.at (iPart) == 13) reader.hepeup.IDUP.at (iPart) = 15 ;
                   if (reader.hepeup.IDUP.at (iPart) == -13) reader.hepeup.IDUP.at (iPart) = -15 ;
                 }
               if (reader.hepeup.IDUP.at (iPart) == 14) reader.hepeup.IDUP.at (iPart) = 16 ;
               if (reader.hepeup.IDUP.at (iPart) == -14) reader.hepeup.IDUP.at (iPart) = -16 ;
             } // outgoing particles
        } // loop over particles in the event
      writer.eventComments() << reader.eventComments;
      writer.hepeup = reader.hepeup;
      writer.writeEvent();

    } //PG loop over input events

  return 0 ;
}
