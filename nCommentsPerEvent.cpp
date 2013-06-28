// c++ -o nCommentsPerEvent -lm nCommentsPerEvent.cpp

#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <sstream>

using namespace std ;



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
  int commentsCounter = 0 ;
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      counter++ ;
      vector<string> comments ;
      splitString (comments, reader.eventComments, '\n') ;
      commentsCounter += comments.size () ;


    } //PG loop over input events

  cout << commentsCounter * 1. / counter << endl ;
  return counter ;
}
