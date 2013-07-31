// c++ -o mixLHEfiles2 mixLHEfiles2.cpp

#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <fstream>

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


int main (int argc, char ** argv) 
{

  if(argc < 3)
    {
      cout << "Usage:   " << argv[0] 
           << " filesList.txt outputFile.lhe" << endl ;
      return -1;
    }

  vector<string> fileToAddNames ;
  vector<double> mixEventsXS ;
  double totXS = 0. ;

  string line;
  ifstream myfile (argv[1]);
  if (myfile.is_open ())
    {
      while ( myfile.good ())
        {
          getline (myfile,line) ;
          cout << line << endl ;
          stringstream ss (line) ;
          string name ;  
          ss >> name ;
          fileToAddNames.push_back (name) ;
          double xs ;    
          ss >> xs ;
          mixEventsXS.push_back (xs) ;
          totXS += xs ;
          cout << "--> " << fileToAddNames.back () << " " << mixEventsXS.back () << "\n" ;          
        }
      myfile.close () ;
    }
  else
    { 
      cout << "Unable to open file\n" ;
      exit (1) ;
    } 

  vector<int> totEventsNum (fileToAddNames.size (), 0) ;
  int grandTotal = 0 ;

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
      grandTotal += totEventsNum.at (iFile) ;
      cout << "found " << totEventsNum.at (iFile) << " in " << fileToAddNames.at (iFile) << endl ;  
    } //PG loop over input files


  //PG find the conversion factor between XS and Nevents, based on the file with the largest XS
  //PG (this should work when all files have the same stats, there's a check for having enough events though)   

  int maxXS_pos = max_element (mixEventsXS.begin (), mixEventsXS.end ()) - mixEventsXS.begin () ;
  double conv_factor = totEventsNum.at (maxXS_pos) / mixEventsXS.at (maxXS_pos) ;

  //PG determine how many events have to be picked up from each file
  vector<int> mixEventsNum (fileToAddNames.size (), 0) ;
  for (int iFile = 0 ; iFile < fileToAddNames.size () ; ++iFile)
    {
      mixEventsNum.at (iFile) = round (mixEventsXS.at (iFile) * conv_factor) ;
      if (mixEventsNum.at (iFile) > totEventsNum.at (iFile))
        {
          cout << "not enough events in file: " << fileToAddNames.at (iFile) << endl ;
        }
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
  ofstream outputStream (argv[2]) ;  
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
    int iInput = getIndex (mixEventsNum, grandTotal) ;

    if (copiedEventsNum.at (iInput) == mixEventsNum.at (iInput)) continue ;

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
      if (copiedEventsNum.at (iFile) < mixEventsNum.at (iFile)) cont += 1 ;
    
  } while (cont != 0) ;

  for (int iFile = 0 ; iFile < fileToAddNames.size () ; ++iFile)
    {
      cout << "closing " << fileToAddNames.at (iFile) << endl ;
      delete readers.at (iFile) ;
      ifs.at (iFile)->close () ;
    } 

  cout << "output file : " << argv[2] 
       << " containing " << events 
       << " for a total XS of " << totXS << endl ;
  return 0 ;
}
