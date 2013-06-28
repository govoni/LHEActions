// c++ -o mixLHEfiles mixLHEfiles.cpp

#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std ;

int getIndex (const vector<int> & pops, int totot) 
  {
    float pos = rand () ;
    pos /= RAND_MAX ;
    pos *= totot ; 
    int index = 0 ;
    int tot = 0 ; 
    for (int i = 0 ; i < pops.size () ; ++i) 
      {
        tot += pops.at (i) ;
        if (pos < tot) 
          {
            index = i ;
            break ;
          }
      }
    return index ;
  }
  
 
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----  


int main(int argc, char ** argv) {

  if(argc < 3)
    {
      cout << "Usage:   " << argv[0] 
           << " fileToAdd1.lhe   [fileToAdd2.lhe [...]]" << endl ;
      return -1;
    }
  
  vector<string> fileToAddNames ;
  for(int iFile = 0; iFile < argc-1; ++iFile)
    {
      fileToAddNames.push_back (argv[1+iFile]) ;
      cout << "fileToAddName = " << fileToAddNames.at (iFile) << endl ;
    }

  vector<int> totEventsNum (fileToAddNames.size (), 0) ;

  //PG get the number of events per file
  for (int iFile = 0 ; iFile < fileToAddNames.size () ; ++iFile)
    {
      cout << "opening " << fileToAddNames.at (iFile) << endl ;
      ifstream ifs (fileToAddNames.at (iFile).c_str ());

      // Create the Reader object:
      LHEF::Reader reader (ifs) ;

      //PG loop over events
      while ( reader.readEvent() ) 
        {
          if ( reader.outsideBlock.length() ) cout << reader.outsideBlock;
          ++totEventsNum.at (iFile) ;
        } //PG loop over events
    } //PG loop over input files
   
  int grandTotal = 0 ;
  for (int iFile = 0 ; iFile < fileToAddNames.size () ; ++iFile)
    {
      grandTotal += totEventsNum.at (iFile) ;
      cout << "found " << totEventsNum.at (iFile) << " in " << fileToAddNames.at (iFile) << endl ;  
    } 

  //PG re-open inputs for reading
  vector<ifstream*> ifs ;
  vector<LHEF::Reader*> readers ;
  
  for (int iFile = 0 ; iFile < fileToAddNames.size () ; ++iFile)
    {
      cout << "opening " << fileToAddNames.at (iFile) << endl ;
      ifstream * dummy = new ifstream (fileToAddNames.at (iFile).c_str ()) ;
      ifs.push_back (dummy) ;
      
      readers.push_back (new LHEF::Reader (*ifs.at (iFile))) ;
    } 
  
  //PG output file
  ofstream outputStream ("total.lhe") ;  
  LHEF::Writer writer (outputStream) ;
  writer.headerBlock () << readers.at (0)->headerBlock ;
  writer.initComments () << readers.at (0)->initComments ;
  writer.heprup = readers.at (0)->heprup ;
  writer.init () ;

  //PG fill the output
  vector<int> copiedEventsNum (fileToAddNames.size (), 0) ;
  int events = 0 ;
  int cont = 0 ;
  do {
    //PG choose randomly from which input file to read the event
    int iInput = getIndex (totEventsNum, grandTotal) ;

    if (copiedEventsNum.at (iInput) == totEventsNum.at (iInput)) continue ;

    //PG copy the event from the input to the output
    readers.at (iInput)->readEvent () ;
    if (readers.at (iInput)->outsideBlock.length ()) cout << readers.at (iInput)->outsideBlock ;
    writer.eventComments () << readers.at (iInput)->eventComments ;
    writer.hepeup = readers.at (iInput)->hepeup ;
    writer.writeEvent () ;
    ++copiedEventsNum.at (iInput) ;
    ++events ;
    if (events % 10000 == 0) cout << events << " read" << endl ; 

    //PG continue util all the reading files have been read
    cont = 0 ;
    for (int iFile = 0 ; iFile < fileToAddNames.size () ; ++iFile)
      if (copiedEventsNum.at (iFile) < totEventsNum.at (iFile)) cont += 1 ;
    
  } while (cont != 0) ;

  for (int iFile = 0 ; iFile < fileToAddNames.size () ; ++iFile)
    {
      cout << "closing " << fileToAddNames.at (iFile) << endl ;
      delete readers.at (iFile) ;
      ifs.at (iFile)->close () ;
    } 

  cout << "output file : " << "total.lhe" << " containing " << events << endl ;
  return 0 ;
}
