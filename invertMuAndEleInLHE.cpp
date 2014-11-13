// c++ -o invertMuAndEleInLHE `root-config --glibs --cflags` -lm invertMuAndEleInLHE.cpp
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

//PG ele and mu massless in phantom
//  float k2 = 0.1056583715 * 0.1056583715 - 0.000510998928 * 0.000510998928 ; // GeV 0.0111634303481
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
               // mu to e
               if (reader.hepeup.IDUP.at (iPart) == 13) reader.hepeup.IDUP.at (iPart) = 11 ;
               else if (reader.hepeup.IDUP.at (iPart) == -13) reader.hepeup.IDUP.at (iPart) = -11 ;
               else if (reader.hepeup.IDUP.at (iPart) == 14) reader.hepeup.IDUP.at (iPart) = 12 ;
               else if (reader.hepeup.IDUP.at (iPart) == -14) reader.hepeup.IDUP.at (iPart) = -12 ;
               // e to mu
               else if (reader.hepeup.IDUP.at (iPart) == 11) reader.hepeup.IDUP.at (iPart) = 13 ;
               else if (reader.hepeup.IDUP.at (iPart) == -11) reader.hepeup.IDUP.at (iPart) = -13 ;
               else if (reader.hepeup.IDUP.at (iPart) == 12) reader.hepeup.IDUP.at (iPart) = 14 ;
               else if (reader.hepeup.IDUP.at (iPart) == -12) reader.hepeup.IDUP.at (iPart) = -14 ;
             } // outgoing particles
        } // loop over particles in the event
      writer.eventComments() << reader.eventComments;
      writer.hepeup = reader.hepeup;
      writer.writeEvent();

    } //PG loop over input events

  return 0 ;
}
