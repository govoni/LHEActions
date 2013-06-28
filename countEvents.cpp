// c++ -o countEvents -lm countEvents.cpp

#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <sstream>

using namespace std ;



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

  int counter = 0 ;
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      counter++ ;


    } //PG loop over input events

  cout << counter << endl ;
  return counter ;
}
