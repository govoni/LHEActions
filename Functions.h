#ifndef Functions_h
#define Functions_h

#include <iostream>

#include "TMath.h"






// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
// SIMPLE FUNCTIONS
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

/*** powerLaw ***/
double powerLaw(double* x, double* par);

/*** exponential ***/
double exponential(double* x, double* par);

/*** double exponential ***/
double doubleExponential(double* x, double* par);

/*** turn-on fermi function ***/
double antiFermi(double* x, double* par);

/*** gaussian + exponential tail smoothly joined ***/
double superGausCumCauda(double* x, double* par);

/*** gaussian + exponential tail smoothly joined ***/
double doubleSuperGausCumCauda(double* x, double* par);

/*** parabola + double exponential tail smoothly joined ***/
double superParabolaCumDoubleCauda(double* x, double* par);

/*** single crystal ball with high tail ***/
double crystalBallHigh(double* x, double* par);

/*** crystall ball with low tail ***/
double crystalBallLow(double* x, double* par);

/*** double crystall ball ***/
double crystalBallLowHigh(double* x, double* par);

/*** gaussian ***/
double gaussian(double* x, double* par);

/*** breit-wigner ***/
double breitWigner(double* x, double* par);

/*** breit-wigner convoluted with crystalBall ***/
double breitWigner_crystalBallLow(double* x, double* par);

/*** breit-wigner convoluted with a gaussian ***/
double breitWigner_gaussian(double* x, double* par);

/*** exp-polynomial of second order ***/
double expPol2order(double* x, double* par);

/*** exp-polynomial of third order ***/
double expPol3order(double* x, double* par);

/*** exp-polynomial of fourth order ***/
double expPol4order(double* x, double* par);

/*** exp-polynomial of fifth order ***/
double expPol5order(double* x, double* par);

/*** polynomial of fourth order ***/
double pol4order(double* x, double* par);

/*** polynomial of fifth order ***/
double pol5order(double* x, double* par);

/*** polynomial of fourth order in the Bernstein base ***/
double pol4orderBernstein(double* x, double* par);






// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
// COMPOSITE FUNCTIONS
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

double attenuatedPowerLaw(double* x, double* par);

double attenuatedExponential(double* x, double* par);

double attenuatedDoubleExponential(double* x, double* par);

double attenuatedDoubleExponentialCumLeadingEdge(double* x, double* par);

double attenuatedCrystalBallHigh(double* x, double* par);

double attenuatedCrystalBallLowHigh(double* x, double* par);

double attenuatedExpPol2order(double* x, double* par);
double attenuatedExpPol3order(double* x, double* par);
double attenuatedExpPol4order(double* x, double* par);
double attenuatedExpPol5order(double* x, double* par);

double antiFermiWithScale (double* x, double* par);

#endif
