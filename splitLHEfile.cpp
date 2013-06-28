// c++ -o splitLHEfile splitLHEfile.cpp 
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

using namespace std ;


int main(int argc, char ** argv) 
{
  if(argc < 4)
    {
      cout << "Usage:   " << argv[0] 
           << " fileToSplit.lhe events_per_file tag_for_splitting" << endl ;
      return -1;
    }

  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  int eventsPerFile = atoi (argv[2]) ;

  LHEF::Writer * writer ;
  ofstream outputStream ;

  int ieve = 0 ;
  int index = 0 ;
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      if (ieve == 0)
        { 
          stringstream filename ;
          filename << argv[3] << "_" << index << ".lhe" ;
          outputStream.open (filename.str ().c_str ()) ;
          cout << "opening in output : " << filename.str () << endl ;
          writer = new LHEF::Writer (outputStream) ;
          writer->headerBlock() << reader.headerBlock ;
          writer->initComments() << reader.initComments ;
          writer->heprup = reader.heprup ;
          writer->init () ;
        }

      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
      writer->eventComments() << reader.eventComments;
      writer->hepeup = reader.hepeup;
      writer->writeEvent();
      ieve++ ;

      if (ieve % eventsPerFile == 0)
        {
          ieve = 0 ;
          index++ ;
          delete writer ;
          outputStream.close () ;
        }
    } //PG loop over input events

  if (ieve % eventsPerFile != 0)  delete writer ;
  return 0 ;
}
