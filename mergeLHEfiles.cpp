// c++ -o mergeLHEfiles mergeLHEfiles.cpp 
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>

using namespace std ;

int main(int argc, char ** argv) {

  if(argc < 3)
    {
      cout << "Usage:   " << argv[0] 
           << " fileToAdd1.lhe   [fileToAdd2.lhe [...]]" << endl ;
      return -1;
    }
  
  vector<string> fileToAddNames;
  for(int fileIt = 0; fileIt < argc-1; ++fileIt)
    {
      fileToAddNames.push_back( argv[1+fileIt] );
      cout << "fileToAddName = " << fileToAddNames.at(fileIt) << endl ;
    }
  
  ofstream outputStream ("total.lhe") ;
  
  LHEF::Writer writer (outputStream) ;

  //PG loop over input files
  for (int i = 0 ; i < fileToAddNames.size () ; ++i)
    {
      cout << "opening " << fileToAddNames.at (i) << endl ;
      std::ifstream ifs(fileToAddNames.at (i).c_str ());

      // Create the Reader object:
      LHEF::Reader reader(ifs);
      if (i == 0)
        {
          writer.headerBlock() << reader.headerBlock;
          writer.initComments() << reader.initComments;
          writer.heprup = reader.heprup;
          writer.init();
        }
      //PG loop over events
      while ( reader.readEvent() ) 
        {
          if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;
          writer.eventComments() << reader.eventComments;
          writer.hepeup = reader.hepeup;
          writer.writeEvent();
        } //PG loop over events
    } //PG loop over input files

  cout << "output file : " << "total.lhe" << endl ;
  return 0 ;
}
