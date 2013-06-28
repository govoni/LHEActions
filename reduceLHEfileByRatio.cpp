// c++ -o reduceLHEfileByRatio reduceLHEfileByRatio.cpp

#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>


int main(int argc, char ** argv) {

 if(argc < 5) {
  std::cout << "Usage:   " << argv[0] << " fileToReduce.lhe reduction-factor newname" << std::endl ;
  return -1;
 }

 std::ifstream ifs (argv[1]) ;
 LHEF::Reader reader (ifs) ;

 float ratio = atof (argv[2]) ; ///--- total --> total / ratio
//  std::string doIt = "./countEvents " + std::string(argv[1]);
//  std::cout << " ratio = " << ratio << std::endl;
//  std::cout << " doIt = " << doIt << std::endl;
//  int eventsPerFile = system(doIt.c_str());
 int eventsPerFile = atoi(argv[4]);
 std::cout << " total events  = " << eventsPerFile << std::endl;
 eventsPerFile = eventsPerFile * ratio;
 std::cout << " events to be saved = " << eventsPerFile << std::endl;
 LHEF::Writer * writer ;
 std::ofstream outputStream ;

 int ieve = 0 ;
 int index = 0 ;
 //PG loop over input events
 while (reader.readEvent ()) {
  if (ieve == 0) {
   std::stringstream filename ;
   filename << argv[3] << ".lhe" ;
   outputStream.open (filename.str ().c_str ()) ;
   std::cout << "opening in output : " << filename.str () << std::endl ;
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

  if (ieve % eventsPerFile == 0) {
   ieve = 0 ;
   index++ ;
   delete writer ;
   outputStream.close () ;
   return 0 ;
  }
 } //PG loop over input events

}
