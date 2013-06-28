// c++ -o varyXS_01 `root-config --glibs --cflags` -lm varyXS_01.cpp
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "LHEF.h"
#include "TLorentzVector.h"
#include <algorithm>
#include <functional>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"


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


int main(int argc, char ** argv) 
{
  if(argc < 2)
    {
      cout << "Usage:   " << argv[0] 
           << " input.lhe" << endl ;
      return -1;
    }


  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  TH1F * h_higgsPt[7] ;
  TString suffix[7] = {"1_1", "h_h", "1_h", "h_1", "1_2", "2_1", "2_2"} ;
  for (int iScale = 0 ; iScale < 7 ; ++iScale)
    {
      TString name = "h_higgsPt_" + suffix[iScale] ;
      TH1F * dummy = new TH1F (name, name, 250 , 0 , 1000 ) ;
      h_higgsPt[iScale] = dummy ;
    }
    
  int counter = 0 ;  
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
      if (counter++ % 10000 == 0) 
        cerr << "reading event " << counter << endl ;

//      cout << reader.eventComments << endl ;
      vector<string> comments ;
      splitString (comments, reader.eventComments, '\n') ;
//      cout << comments.size () << endl ;
//      cout << comments.at (1) << endl ;
      vector<double> weights ;
      for (int iString = 1 ; iString < comments.size () ; ++iString)
        {
          vector<string> words ;
          splitString (words, comments.at (iString), ' ') ;
          float renFact = atof (words.at (3).c_str ()) ;
          float facFact = atof (words.at (4).c_str ()) ;
          double weight = atof (words.at (2).c_str ()) ;
//          if (renFact / facFact < 0.5 || renFact / facFact > 2) continue ;
          weights.push_back (weight) ;
        }
//      cout << weights.size () << endl ;
      if (weights.size () < 7)
        {
          cerr << "ERROR too few weigths: " << weights.size () 
               << " in event " << (counter - 1) << endl ;
          cerr << "skipping event" << endl ;
          continue ;     
        }

      TLorentzVector momentum ;
      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size () ; ++iPart) 
        {
           // outgoing particles          
           if (reader.hepeup.ISTUP.at (iPart) == 1)
             {
               if (abs (reader.hepeup.IDUP.at (iPart)) == 25) 
                 { 
                   momentum = buildP (reader.hepeup, iPart) ;
                 }
             } // outgoing particles
        } // loop over particles in the event  
      for (int iScale = 0 ; iScale < 7 ; ++iScale)
        {
          h_higgsPt[iScale]->Fill (momentum.Pt (), weights.at (iScale)) ;
//          cout << weights.at (iScale) << endl ;
        }

    } //PG loop over input events

  cout << "CROSS SECTION " << h_higgsPt[0]->Integral () << endl ;
  cout << "entries " << h_higgsPt[0]->GetEntries () << endl ;
 

  TFile out ("varyXS_01.root", "recreate") ;
  out.cd () ;
  for (int iScale = 0 ; iScale < 7 ; ++iScale)
    h_higgsPt[iScale]->Write () ;
  out.Close () ;

  return 0 ;
}
