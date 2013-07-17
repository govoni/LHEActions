// c++ -o pdfReweight `root-config --glibs --cflags` `lhapdf-config --cppflags  --ldflags` -lm pdfReweight.cpp
#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>

#include "LHAPDF/LHAPDF.h"

#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"
// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math ;
using namespace std ;
using namespace LHAPDF ;

typedef LorentzVector<ROOT::Math::PxPyPzE4D<double> > lorentzVector ;



struct ptsort: public std::binary_function<TLorentzVector, TLorentzVector, bool>
{
  bool operator () (const TLorentzVector & x, const TLorentzVector & y)
    {
      return  (x.Perp () < y.Perp () ) ;
    }
} ;


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


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main(int argc, char ** argv) 
{
  if(argc < 2)
    {
      cout << "Usage:   " << argv[0] 
           << " input.lhe " << endl ;
      return -1;
    }

  bool REW = false ;
  if (argc > 2) REW = true ;
  cout << REW << endl ;
 
  const int SUBSET = 0;
  const string NAME = "cteq6ll"; //"cteq6l1"

  LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
  const int NUMBER = LHAPDF::numberPDF();

  LHAPDF::initPDF (0) ;

  std::ifstream ifs (argv[1]) ;
  LHEF::Reader reader (ifs) ;

  TH1F h_phScale ("h_phScale", "h_phScale", 100, 0, 1000) ;
  TH1F h_scale ("h_scale", "h_scale", 100, 0, 1000) ;
  TH1F h_x ("h_x", "h_x", 100, 0, 1) ;
  TH1F h_weight ("h_weight", "h_weight", 100, 0, 2) ;
  TH1F h_phaseSp ("h_phaseSp", "h_phaseSp", 1000, 0, 10) ;

  TH1F h_ptj1 ("h_ptj1", "h_ptj1", 60, 0, 400) ;
  TH1F h_ptj2 ("h_ptj2", "h_ptj2", 50, 0, 300) ;
  TH1F h_etaj1 ("h_etaj1", "h_etaj1", 40, 0, 6) ;
  TH1F h_etaj2 ("h_etaj2", "h_etaj2", 40, 0, 6) ;
  TH1F h_mjj ("h_mjj", "h_mjj", 50, 0, 3500) ;
  TH1F h_detajj ("h_detajj", "h_detajj", 100, 30, 10) ;
  TH1F h_mll ("h_mll", "h_mll", 40, 0, 300) ;

  h_phScale.Sumw2 () ;
  h_scale.Sumw2 () ;
  h_x.Sumw2 () ;
  h_weight.Sumw2 () ;
  h_phaseSp.Sumw2 () ;
  h_ptj1.Sumw2 () ;
  h_ptj2.Sumw2 () ;
  h_etaj1.Sumw2 () ;
  h_etaj2.Sumw2 () ;
  h_mjj.Sumw2 () ;
  h_detajj.Sumw2 () ;
  h_mll.Sumw2 () ;

  int number_total = 0 ;
  int number_selec = 0 ;
  //PG loop over input events
  while (reader.readEvent ()) 
    {
      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

      ++number_total ;
      TLorentzVector Higgs;
      int iPartHiggs;
      
      std::vector<int> leptonsFlavour ;      
      std::vector<int> finalJets ;      
      std::vector<int> initialQuarks ;      
      
      vector<TLorentzVector> v_f_jets ; //PG w/ b's
      vector<TLorentzVector> v_f_quarks ; //PG w/o b's
      vector<TLorentzVector> v_f_leptons ;
      vector<TLorentzVector> v_f_neutrinos ;
      vector<TLorentzVector> v_f_intermediate ;
      
      int nele = 0 ;
      int nmu  = 0 ;
      int ntau = 0 ;

      double x[2] = {0., 0.} ;
      int flavour[2] = {0, 0} ;
    
      // loop over particles in the event
      // and fill the variables of leptons and quarks
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
           TLorentzVector dummy = buildP (reader.hepeup, iPart) ;

           // incoming particle          
           if (reader.hepeup.ISTUP.at (iPart) == -1) 
             {
               initialQuarks.push_back (iPart) ;
               x[iPart] = dummy.P () / 4000. ;
               flavour[iPart] = reader.hepeup.IDUP.at (iPart) ;

             } // incoming particle          


           if (reader.hepeup.ISTUP.at (iPart) == 2) 
             {
               if (abs (reader.hepeup.IDUP.at (iPart)) == 6)
                 {
                   v_f_intermediate.push_back (dummy) ;
                 }
             }
           
           // outgoing particles          
           if (reader.hepeup.ISTUP.at (iPart) == 1)
             {
               // quarks
               if (abs (reader.hepeup.IDUP.at (iPart)) < 7) 
                 {
                   v_f_jets.push_back (dummy) ;
                   finalJets.push_back (iPart) ;
                   if (abs (reader.hepeup.IDUP.at (iPart)) < 5)
                     {
                       v_f_quarks.push_back (dummy) ;        
                     }
                 } // quarks
               else if (abs (reader.hepeup.IDUP.at (iPart)) == 11 ||
                        abs (reader.hepeup.IDUP.at (iPart)) == 13 ||
                        abs (reader.hepeup.IDUP.at (iPart)) == 15)
                 {
                   leptonsFlavour.push_back (reader.hepeup.IDUP.at (iPart)) ;
                   v_f_leptons.push_back (dummy) ;
                   nele += (abs (reader.hepeup.IDUP.at (iPart)) == 11) ;
                   nmu  += (abs (reader.hepeup.IDUP.at (iPart)) == 13) ;
                   ntau += (abs (reader.hepeup.IDUP.at (iPart)) == 15) ;
                 }
               else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                        abs (reader.hepeup.IDUP.at (iPart)) == 14 ||
                        abs (reader.hepeup.IDUP.at (iPart)) == 16)
                 {
                   v_f_neutrinos.push_back (dummy) ;        
                 }
             } // outgoing particles
        } // loop over particles in the event

      //PG selections to equalise the two samples
      //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

      if (ntau > 0) continue ;
      
      if (v_f_quarks.size () != finalJets.size ()) continue ;    //PG no bs in the final state
      if (v_f_quarks.size () < 2) continue ;
      if (v_f_intermediate.size () > 0) continue ;               //PG no intermediate tops

      sort (v_f_quarks.rbegin (), v_f_quarks.rend (), ptsort ()) ;  
  
      //PG residual differences in the generation thresholds
      //PG selections like this work only if one expects only 2 quarks in the event
      //PG otherwise need to reduce the quarks collection first
  
      if (v_f_quarks.at (0).Pt () < 20 ||
          v_f_quarks.at (1).Pt () < 20 ||
          v_f_quarks.at (0).E () < 20 ||
          v_f_quarks.at (1).E () < 20 ||
          fabs (v_f_quarks.at (0).Eta ()) > 5 ||
          fabs (v_f_quarks.at (1).Eta ()) > 5 ||
          fabs (v_f_leptons.at (0).Eta ()) > 2.5 ||
          fabs (v_f_leptons.at (1).Eta ()) > 2.5 ||
          v_f_leptons.at (0).Pt () < 6 ||
          v_f_leptons.at (1).Pt () < 6) continue ;
      
//      if (v_f_leptons.at (0).DeltaR (v_f_leptons.at (1)) < 0.4) continue ;          
//      if (v_f_quarks.at (0).DeltaR (v_f_quarks.at (1)) < 0.4) continue ;          
//  
//      int cont = 0 ;
//      for (int iL = 0 ; iL < 2 ; ++iL)
//        for (int iJ = 0 ; iJ < 2 ; ++iJ)
//          if (v_f_quarks.at (iJ).DeltaR (v_f_leptons.at (iL)) < 0.4) cont = 1 ;
//      if (cont == 1) continue ;
      
      TLorentzVector diLepton = v_f_leptons.at (0) + v_f_leptons.at (1) ;
      if (diLepton.M () < 12) continue ;

      if (v_f_neutrinos.at (0).Pt () < 6 || v_f_neutrinos.at (1).Pt () < 6) continue ; 
//      TLorentzVector missingEnergy = v_f_neutrinos.at (0) + v_f_neutrinos.at (1) ;
//      if (missingEnergy.Pt () < 6) continue ;

      TLorentzVector diJet = v_f_quarks.at (0) + v_f_quarks.at (1) ;
      if (diJet.M () < 100) continue ;
    
      //PG only different flavour
      if (abs (leptonsFlavour[0]) == abs (leptonsFlavour[1])) continue ;
//      if (leptonsFlavour[0] == -1 * leptonsFlavour[1] &&
//          diLepton.M () < 97.5 && diLepton.M () > 83.5) continue ;

      //PG VBF cuts
      if (diJet.M () < 500 || 
          fabs (v_f_quarks.at (0).Eta () - v_f_quarks.at (1).Eta ()) < 3.5) continue ;

      //PG opposite charge leptons
      if (leptonsFlavour.at (0) * leptonsFlavour.at (1) > 0) continue ;

      sort (v_f_leptons.rbegin (), v_f_leptons.rend (), ptsort ()) ;  

      //PG analysis state cuts
      if (v_f_quarks.at (0).Pt () < 30 ||
          v_f_quarks.at (1).Pt () < 30 ||
          v_f_leptons.at (0).Pt () < 20 ||
          v_f_leptons.at (1).Pt () < 10) continue ;
     

      //PG at this point the two generations should be on the same page
      //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
      ++number_selec ;

      //PG the scale:
      float scale = reader.hepeup.SCALUP ;

      //PG determine the scale according to the phantom recipe
      
      double phantomScale = 80.385 * 80.385 + 
          (v_f_quarks.at (0).Pt () * v_f_quarks.at (0).Pt () +
           v_f_quarks.at (1).Pt () * v_f_quarks.at (1).Pt () +
           v_f_leptons.at (0).Pt () * v_f_leptons.at (0).Pt () +
           v_f_leptons.at (1).Pt () * v_f_leptons.at (1).Pt () +
           v_f_neutrinos.at (0).Pt () * v_f_neutrinos.at (0).Pt () +
           v_f_neutrinos.at (1).Pt () * v_f_neutrinos.at (1).Pt ()) / 6. ;
      phantomScale = sqrt (phantomScale) ;

      //PG calculate the weight to be applied to the event
      double weight = LHAPDF::xfx (x[0], phantomScale, flavour[0]) * LHAPDF::xfx (x[1], phantomScale, flavour[1]) /
                      (LHAPDF::xfx (x[0], scale, flavour[0]) * LHAPDF::xfx (x[1], scale, flavour[1])) ;
//      cout << weight << endl ;

      h_weight.Fill (weight) ; 
      h_phaseSp.Fill (LHAPDF::xfx (x[0], scale, flavour[0]) * LHAPDF::xfx (x[1], scale, flavour[1])) ;
      h_phScale.Fill (phantomScale) ; 

      //PG fill histos
      if (!REW) weight = 1. ;
      
      h_ptj1.Fill (v_f_quarks.at (0).Pt (), weight) ; 
      h_ptj2.Fill (v_f_quarks.at (1).Pt (), weight) ; 
      h_etaj1.Fill (v_f_quarks.at (0).Eta (), weight) ; 
      h_etaj2.Fill (v_f_quarks.at (1).Eta (), weight) ; 
      h_scale.Fill (scale, weight) ;
      h_x.Fill (x[0], weight) ;
      h_x.Fill (x[1], weight) ;
      h_mjj.Fill (diJet.M (), weight) ;
      h_detajj.Fill (fabs (v_f_quarks.at (0).Eta () - v_f_quarks.at (1).Eta ())) ;
      h_mll.Fill (diLepton.M (), weight) ;
    } //PG loop over input events

  cout << "total " << number_total << endl ;
  cout << "efficiency " << number_selec * 1. / number_total << endl ;

  TFile f ("reweight.root", "recreate") ;
  h_phaseSp.Write () ;
  h_weight.Write () ;
  h_phScale.Write () ;
  h_scale.Write () ;
  h_x.Write () ;
  h_mjj.Write () ;
  h_detajj.Write () ;
  h_mll.Write () ;
  h_ptj1.Write () ;
  h_ptj2.Write () ;
  h_etaj1.Write () ;
  h_etaj2.Write () ;
  f.Close () ;

  return 0 ;
}
