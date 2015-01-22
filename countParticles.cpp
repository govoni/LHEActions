// c++ -o countParticles `root-config --glibs --cflags` -lm countParticles.cpp
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace std ;


int main(int argc, char ** argv) 
{
  if (argc < 2)
    {
      cout << "Usage:   " << argv[0] 
           << " input.lhe" << endl ;
      return -1;
    }


  int maxEvents = -1 ;
  if (argc >= 3)
    {
      maxEvents = atoi (argv[2]) ;
    }

  bool printList = false ;
  if (argc == 4) printList = true ;


  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  TH1F h_inp ("h_inp", "h_inp", 61, -30.5, 30.5) ;
  TH1F h_out ("h_out", "h_out", 61, -30.5, 30.5) ;
  TH1F h_int ("h_int", "h_int", 61, -30.5, 30.5) ;
   
  int counter = 0 ; 
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      if (maxEvents > 0 && counter > maxEvents) break ;
      if ( reader.outsideBlock.length ()) std::cout << reader.outsideBlock;

      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
           if (reader.hepeup.ISTUP.at (iPart) ==  1) 
             {
               h_out.Fill (reader.hepeup.IDUP.at (iPart)) ;
               if (printList) cout << reader.hepeup.IDUP.at (iPart) << "\t" ;
             }

           if (reader.hepeup.ISTUP.at (iPart) ==  2) h_int.Fill (reader.hepeup.IDUP.at (iPart)) ;
           if (reader.hepeup.ISTUP.at (iPart) == -1) h_inp.Fill (reader.hepeup.IDUP.at (iPart)) ;
             
        } // loop over particles in the event
      if (printList) cout << "\n" ;
      ++counter ;
    } //PG loop over input events

  TCanvas c1 ("c1", "c1", 100, 100, 600, 600) ;
  h_out.SetStats (0) ;
  h_out.Draw () ;
  c1.Print ("outgoing.pdf", "pdf") ;
  h_int.SetStats (0) ;
  h_int.Draw () ;
  c1.Print ("intermediate.pdf", "pdf") ;
  h_inp.SetStats (0) ;
  h_inp.Draw () ;
  c1.Print ("incoming.pdf", "pdf") ;

  TFile f ("countParticles.root", "recreate") ;
  h_inp.Write () ;
  h_out.Write () ;
  h_int.Write () ;
  f.Close () ;

  return 0 ;
}
