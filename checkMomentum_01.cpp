// c++ -o checkMomentum_01 `root-config --glibs --cflags` -lm checkMomentum_01.cpp
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <sstream>

#include "TLorentzVector.h"

using namespace std ;


// find daughters for a particle
vector<int> daughters (const LHEF::HEPEUP & event, int PID)
  {
    vector<int> result ;
    for (int iPart = 0 ; iPart < event.IDUP.size () ; ++iPart) 
      {
        pair <int, int> mothers = event.MOTHUP.at (iPart) ;
        if (mothers.first != mothers.second) continue ;
        if (mothers.first == PID) result.push_back (iPart) ;
      }
    return result ;  
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


int main(int argc, char ** argv) 
{
  if(argc < 3)
    {
      cout << "Usage:   " << argv[0] 
           << " input.lhe lepton_ID" << endl ;
      return -1;
    }

  int lepton_ID = atoi (argv[2]) ;

  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  //PG loop over input events
  while (reader.readEvent ()) 
    {
      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size () ; ++iPart) 
        {
           // outgoing particles          
           if (reader.hepeup.ISTUP.at (iPart) == 1)
             {
               if (abs (reader.hepeup.IDUP.at (iPart)) == lepton_ID) 
                 { 
                   // find mother
                   pair <int, int> mothers = reader.hepeup.MOTHUP.at (iPart) ;
                   cout << "----- " << reader.hepeup.IDUP.at (iPart) << " " << iPart <<  " \n" ; 
                   cout << "  m1 : " << mothers.first <<  ": " << reader.hepeup.IDUP.at (mothers.first) << endl ;
                   cout << "  m2 : " << mothers.second <<  ": " << reader.hepeup.IDUP.at (mothers.second) << endl ;
                   
                   if (mothers.first == mothers.second)
                     {
                       if (abs (reader.hepeup.IDUP.at (mothers.first)) != 24) continue ;
                       vector<int> daugh = daughters (reader.hepeup, mothers.first) ;
                       cout << "  d1 : " << daugh.at (0) << " " << reader.hepeup.IDUP.at (daugh.at (0)) << "\n" ;
                       cout << "  d2 : " << daugh.at (1) << " " << reader.hepeup.IDUP.at (daugh.at (1)) << "\n" ;
                       TLorentzVector moth = buildP (reader.hepeup, mothers.first) ;
                       TLorentzVector daugh1 = buildP (reader.hepeup, daugh.at (0)) ;
                       TLorentzVector daugh2 = buildP (reader.hepeup, daugh.at (1)) ;
                       TLorentzVector daughSum = daugh1 + daugh2 ;
//                       daughSum = daughSum - moth ;
                       cout << daughSum.Mag () << " " << moth.Mag () << endl ;
                     }
                   else 
                     {
                       cout << " the two mothers differ\n" ;
                     }  
                   
                   
                   
                   // find mother's daughters
                   // check momentum conservation    
//                   TLorentzVector dummy
//                     (
//                       reader.hepeup.PUP.at (iPart).at (0), // px
//                       reader.hepeup.PUP.at (iPart).at (1), // py
//                       reader.hepeup.PUP.at (iPart).at (2), // pz
//                       reader.hepeup.PUP.at (iPart).at (3) // E
//                     ) ;
//                   float p2 = dummy.Vect ().Mag2 () ;
                 }
             } // outgoing particles
        } // loop over particles in the event

    } //PG loop over input events

  return 0 ;
}
