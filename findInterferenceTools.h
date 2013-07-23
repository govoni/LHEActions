#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cassert>

#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"
// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "LHAPDF/LHAPDF.h"


using namespace ROOT::Math ;
using namespace std ;

typedef LorentzVector<ROOT::Math::PxPyPzE4D<double> > lorentzVector ;



struct ptsort: public std::binary_function<TLorentzVector, TLorentzVector, bool>
{
  bool operator () (const TLorentzVector & x, const TLorentzVector & y)
    {
      return  (x.Perp () < y.Perp () ) ;
    }
} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


struct etasort: public std::binary_function<TLorentzVector, TLorentzVector, bool>
{
  bool operator () (const TLorentzVector & x, const TLorentzVector & y)
    {
      return  (x.Eta () < y.Eta () ) ;
    }
} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


pair<int, int> findPairWithLargestDeta (const vector<TLorentzVector> & v_f_quarks)
{
  int minimum = min_element (v_f_quarks.begin (), v_f_quarks.end (), etasort ()) - v_f_quarks.begin () ;
  int maximum = max_element (v_f_quarks.begin (), v_f_quarks.end (), etasort ()) - v_f_quarks.begin () ;
  return pair<int, int> (minimum, maximum) ;
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


pair<int, int> findPairWithWMass (const vector<TLorentzVector> & v_f_quarks)
{
  double ref_deltaM = 100000. ;
  int one = 0 ;
  int two = 0 ;
  
  for (int iJ = 0 ; iJ < 4 ; ++iJ)
    for (int iJ2 = iJ + 1 ; iJ2 < 4 ; ++iJ2)
      {
        double deltaM = fabs ((v_f_quarks.at (iJ) + v_f_quarks.at (iJ2)).M () - 80.4) ;
        if (deltaM < ref_deltaM)
          {
            ref_deltaM = deltaM ;
            one = iJ ;
            two = iJ2 ;
          }      
      }
  return pair<int, int> (one, two) ;
  
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


lorentzVector buildLP (const LHEF::HEPEUP & event, int iPart)
{
  lorentzVector dummy ;
  dummy.SetPxPyPzE (
      event.PUP.at (iPart).at (0), // px
      event.PUP.at (iPart).at (1), // py
      event.PUP.at (iPart).at (2), // pz
      event.PUP.at (iPart).at (3) // E
    ) ;
  return dummy ;  
}


