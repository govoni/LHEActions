// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


/*** double crystall ball ***/
double doubleGausCrystalBallLowHigh (double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  //[5] = alpha2
  //[6] = n2

  double xx = x[0];
  double mean   = par[1] ; // mean
  double sigmaP = par[2] ; // sigma of the positive side of the gaussian
  double sigmaN = par[3] ; // sigma of the negative side of the gaussian
  double alpha  = par[4] ; // junction point on the positive side of the gaussian
  double n      = par[5] ; // power of the power law on the positive side of the gaussian
  double alpha2 = par[6] ; // junction point on the negative side of the gaussian
  double n2     = par[7] ; // power of the power law on the negative side of the gaussian

  if ((xx-mean)/sigmaP > fabs(alpha))
  {
    double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
    double B = n/fabs(alpha) - fabs(alpha);
    
    return par[0] * A * pow(B + (xx-mean)/sigmaP, -1.*n);
  }
  
  else if ((xx-mean)/sigmaN < -1.*fabs(alpha2))
  {
    double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
    double B = n2/fabs(alpha2) - fabs(alpha2);
    
    return par[0] * A * pow(B - (xx-mean)/sigmaN, -1.*n2);
  }
  
  else if ((xx-mean) > 0)
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigmaP*sigmaP) );
  }
  
  else
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigmaN*sigmaN) );
  }
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void printArray (double * array, int N)
{
  cout << "double * params[" << N << "] = {" ;
  for (int i = 0 ; i < N - 1 ; ++i)
    cout << array[i] << ", " ;
  cout << array[N-1] << "} ;\n" ;
  return ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void setParNamesdoubleGausCrystalBallLowHigh (TF1 * func)
{
  func->SetParName (1, "mean") ;
  func->SetParName (2, "sigma") ;
  func->SetParName (3, "alphaL") ;
  func->SetParName (4, "nL") ;
  func->SetParName (5, "alphaR") ;
  func->SetParName (6, "nR") ;
  return ;
}  


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t crystalBallLowHighRatio (Double_t * xx, Double_t * par)
{
  return crystalBallLowHigh (xx, par) / crystalBallLowHigh (xx, par + 7) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t withLeftAddOn (Double_t * xx, Double_t * par)
{
   double centre       = par[0] ;
   double shift        = par[1] ;
   double scale        = par[2] ;
   double slope        = par[3] ;
   double secondCentre = par[4] ;
   double secondSlope  = par[4] ;
   double x            = xx[0] - centre ;
   
   if (x > 0) 
   return scale * x * TMath::Exp (-1 * slope * x) + shift +
          (x - secondCentre) * TMath::Exp (-1 * secondSlope * (x - secondCentre)) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t doubleSlope (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double slope_left  = par[4] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return scale * x * TMath::Exp (     slope_left  * x) + shift ;
//   else       return -1 * scale * x * TMath::Exp (-1 * slope_left  * x) + shift ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t parabolicAndExp (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double power_left  = par[4] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return power_left * x * x + scale * x + shift ;
}




// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t cubicAndExp (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double pwr2_left   = par[4] ;
   double pwr3_left   = par[5] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return pwr3_left * x * x * x + pwr2_left * x * x + scale * x + shift ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t sinusAndExp (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double freq_left   = par[4] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return scale * TMath::Sin (freq_left * x) / freq_left + shift ; //pg this guarantees C1 properties
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t logAndExp (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double log_fact    = par[4] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return scale * x * TMath::Log (-1 * log_fact * x) / log_fact + shift ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int macro_findInterferece (string filename, double mass)                                                        
{        

  gSystem->Load ("Functions.cc") ;

  TFile * f = new TFile (filename.c_str ()) ;
  TH1F * h_MWW_phbkgsig = (TH1F *) f->Get ("h_MWW_phbkgsig") ;
  TH1F * h_MWW_phbkg = (TH1F *) f->Get ("h_MWW_phbkg") ;
  TH1F * h_MWW_mg = (TH1F *) f->Get ("h_MWW_mg") ;
                                                                           
  h_MWW_phbkg->SetStats (0) ;
  h_MWW_phbkgsig->SetStats (0) ;
  h_MWW_mg->SetStats (0) ;
  h_MWW_phbkg->SetLineColor (kOrange) ;
  h_MWW_phbkg->SetLineWidth (2) ;
  h_MWW_phbkgsig->SetLineColor (kRed) ;
  h_MWW_phbkgsig->SetLineWidth (2) ;

  TH1F * diff = (TH1F *) h_MWW_phbkgsig->Clone ("diff") ;
  diff->SetTitle ("") ;
  diff->Add (h_MWW_phbkg, -1) ;

  TH1F * ratio = (TH1F *) diff->Clone ("ratio") ;
  ratio->GetXaxis ()->SetTitle (ratio->GetName ()) ;
  ratio->GetYaxis ()->SetTitle ("(SBI - B) / S") ;
  ratio->SetTitle ("") ;
  ratio->Divide (h_MWW_mg) ;

  TH1F * delta = (TH1F *) h_MWW_mg->Clone ("delta") ;
  delta->SetTitle ("") ;
  delta->Add (diff, -1) ;
  delta->SetLineColor (kGreen + 2) ;

  TH1F * relDiff = delta->Clone ("relDiff") ;
  relDiff->Divide (h_MWW_mg) ;

  //PG initial spectra

  TCanvas * c1 = new TCanvas () ;
  c1->DrawFrame (200, 0.00001, 2 * mass, 0.02) ;
  h_MWW_phbkg->Draw ("histsame") ;
  h_MWW_phbkgsig->Draw ("histsame") ;
  h_MWW_mg->Draw ("histsame") ;
  c1->Print ("spectra.pdf", "pdf") ;

  //PG S only, and (SBI - B)

  TCanvas * c4 = new TCanvas () ;
  h_MWW_mg->SetTitle ("") ;
  h_MWW_mg->Draw ("hist") ;
  diff->Draw ("histsame") ;
  c4->Print ("signals.pdf", "pdf") ;

  //PG S only
  
  TCanvas * c4_mg = new TCanvas () ;
  h_MWW_mg->SetTitle ("") ;
  h_MWW_mg->Draw ("EP") ;

  TF1 * func_mg_1 = new TF1 ("func_mg_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_mg_1->SetNpx (10000) ;
  func_mg_1->SetLineWidth (1) ;
  func_mg_1->SetLineColor (kBlue + 2) ;

  setParNamesdoubleGausCrystalBallLowHigh (func_mg_1) ;

  func_mg_1->SetParameter (0, 1.) ;                  // multiplicative scale
  func_mg_1->SetParameter (1, mass) ;                // mean
  func_mg_1->SetParameter (2, h_MWW_mg->GetRMS ()) ; // gaussian sigma
  func_mg_1->SetParameter (3, 1) ;                   // right junction
  func_mg_1->SetParameter (4, 1) ;                   // right power law order
  func_mg_1->SetParameter (5, 1) ;                   // left junction
  func_mg_1->SetParameter (6, 1) ;                   // left power law order

  h_MWW_mg->Fit ("func_mg_1", "+", "", 0.5 * mass - 50, 2 * mass) ;
  c4_mg->Print ("signals_mg.pdf", "pdf") ;

  //PG (SBI - B) only

  TCanvas * c4_ph = new TCanvas () ;
  diff->SetTitle ("") ;
  diff->Draw ("EP") ;

  TF1 * gauss = new TF1 ("gauss", "gaus", 0, 2000) ;
  gauss->SetLineWidth (1) ;
  gauss->SetLineColor (kGray + 2) ;
  diff->Fit ("gauss", "+", "", 0.75 * mass, mass * 1.25) ;

  TF1 * func_ph_1 = new TF1 ("func_ph_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_ph_1->SetNpx (10000) ;
  func_ph_1->SetLineWidth (1) ;
  func_ph_1->SetLineColor (kRed + 2) ;
  
  setParNamesdoubleGausCrystalBallLowHigh (func_ph_1) ;

  func_ph_1->SetParameter (0, 1.) ;                      // multiplicative scale
  func_ph_1->SetParameter (1, gauss->GetParameter (1)) ; // mean
  func_ph_1->SetParameter (2, gauss->GetParameter (2)) ; // gaussian sigma
  func_ph_1->SetParameter (3, 1) ;                       // right junction
  func_ph_1->SetParameter (4, 2) ;                       // right power law order
  func_ph_1->SetParameter (5, 1) ;                       // left junction
  func_ph_1->SetParameter (6, 2) ;                       // left power law order

  diff->Fit ("func_ph_1", "", "", 0.5 * mass - 50, 2 * mass) ;
  func_ph_1->SetParameters (func_ph_1->GetParameters ()) ;
  diff->Fit ("func_ph_1", "+L", "", 0.5 * mass - 50, 2 * mass) ;

  c4_ph->Print ("signals_ph.pdf", "pdf") ;

  //PG S / (SBI - S)

  TF1 * f_ratio = new TF1 ("f_ratio", crystalBallLowHighRatio, 0, 2000, 14);
  f_ratio->SetLineWidth (3) ;
  f_ratio->SetLineColor (kGray + 2) ;
  Double_t f_ratio_pars [14] ;
  func_ph_1->GetParameters (f_ratio_pars) ;
  func_mg_1->GetParameters (f_ratio_pars + 7) ;
  f_ratio->SetParameters (f_ratio_pars) ;

  TCanvas * c3 = new TCanvas () ;
  c3->DrawFrame (200, -0.4, 2 * mass, 12) ;
  ratio->SetTitle ("") ;
  ratio->SetLineColor (kMagenta) ;
  ratio->GetXaxis ()->SetTitle ("mWW") ;
  ratio->GetYaxis ()->SetTitle ("(SBI - B) / S") ;
  ratio->Draw ("EPsame") ;
  f_ratio->Draw ("same") ;
  c3->Print ("ratio.pdf", "pdf") ;

  cout << "mass " << mass << "\n" ;
  cout << "fitting results:\n-----------------\n" ;
  cout << "S: \n" ;
  printArray (func_ph_1->GetParameters (), 7) ;
  cout << "SBI - B: \n" ;
  printArray (func_mg_1->GetParameters (), 7) ;

  return 0 ; //PG FIXME



  //PG (SBI - B) - S and S

  TCanvas * c5 = new TCanvas () ;
  delta->Draw ("hist") ;
  h_MWW_mg->Draw ("histsame") ;

  c5->Print ("delta.pdf", "pdf") ;


//  TF1 * func = new TF1 ("func","[2] * sin([0] * x) + [1]",0, 2000) ;
//  func->SetParameter (0, -0.005) ;
//  func->SetParameter (1, 0) ;
//  func->SetParameter (2, 0.0005) ;
//  func->Draw ("same") ;
//
//  TF1 * func1 = new TF1 ("func1","[2] * sin([0] * sqrt (x)) + [1]",0, 2000) ;
//  func1->SetParameter (0, -0.005) ;
//  func1->SetParameter (1, 0) ;
//  func1->SetParameter (2, 0.0005) ;

//  TF1 * func2 = new TF1 ("func2"," [2] * (x - [0]) * exp (-1 * [3] * abs (x - [0])) + [1]",0, 2000) ;
//
//  func2->FixParameter (0, mass)
//  func2->FixParameter (1, 0.)
//
//  func2->SetParameter (2, 0.0005) ;
//  func2->SetParameter (3, 0.05) ;
//  delta->Fit ("func2", "", "", 200, 2000) ;


  TF1 * func3 = new TF1 ("func3", doubleSlope, 0, 2000, 5) ;
  func3->SetLineWidth (3) ;
  func3->SetLineColor (kGreen + 4) ;
  func3->FixParameter (0, mass) ;
  func3->FixParameter (1, 0.) ;

  func3->SetParameter (2, 0.0005) ;
  func3->SetParameter (3, 0.05) ;
  func3->SetParameter (4, 0.05) ;

  func3->SetParName (0, "centre") ;
  func3->SetParName (1, "shift") ;
  func3->SetParName (2, "scale") ;
  func3->SetParName (3, "slope_right") ;
  func3->SetParName (4, "slope_left") ;

  delta->Fit ("func3", "+", "", 0, 650) ;

  TF1 * func4 = new TF1 ("func4", parabolicAndExp, 0, 2000, 5) ;
  func4->SetLineWidth (1) ;
  func4->SetLineColor (kRed + 2) ;
  func4->FixParameter (0, mass) ;
  func4->FixParameter (1, 0.) ;

  func4->FixParameter (2, func3->GetParameter (2)) ;
  func4->FixParameter (3, func3->GetParameter (3)) ;
  func4->SetParameter (4, 1) ;

  func4->SetParName (0, "centre") ;
  func4->SetParName (1, "shift") ;
  func4->SetParName (2, "scale") ;
  func4->SetParName (3, "slope_right") ;
  func4->SetParName (4, "power_left") ;

  delta->Fit ("func4", "+", "", 200, 2000) ;

//  TF1 * func5 = new TF1 ("func5", sinusAndExp, 0, 2000, 5) ;
//  func5->SetLineWidth (1) ;
//  func5->SetLineColor (kRed + 2) ;
//  func5->FixParameter (0, mass) ;
//  func5->FixParameter (1, 0.) ;
//
//  func5->FixParameter (2, func3->GetParameter (2)) ;
//  func5->FixParameter (3, func3->GetParameter (3)) ;
//  func5->SetParameter (4, 10) ;
//
//  func5->SetParName (0, "centre") ;
//  func5->SetParName (1, "shift") ;
//  func5->SetParName (2, "scale") ;
//  func5->SetParName (3, "slope_right") ;
//  func5->SetParName (4, "freq_left") ;
//
//  delta->Fit ("func5", "+", "", 200, 2000) ;

//  TF1 * func6 = new TF1 ("func6", withLeftAddOn, 0, 2000, 6) ;
//  func6->SetLineWidth (1) ;
//  func6->SetLineColor (kBlue + 2) ;
//  func6->FixParameter (0, mass) ;
//  func6->FixParameter (1, 0.) ;
//
//  func6->FixParameter (2, func3->GetParameter (2)) ;
//  func6->FixParameter (3, func3->GetParameter (3)) ;
//
//  func6->SetParName (0, "centre") ;
//  func6->SetParName (1, "shift") ;
//  func6->SetParName (2, "scale") ;
//  func6->SetParName (3, "slope") ;
//  func6->SetParName (4, "secondCentre") ;
//  func6->SetParName (5, "secondSlope") ;
//
//  delta->Fit ("func6", "+", "", 200, 2000) ;

  TF1 * func7 = new TF1 ("func7", cubicAndExp, 0, 2000, 6) ;
  func7->SetLineWidth (1) ;
  func7->SetLineColor (kViolet + 2) ;
  func7->FixParameter (0, mass) ;
  func7->FixParameter (1, 0.) ;

  func7->FixParameter (2, func3->GetParameter (2)) ;
  func7->FixParameter (3, func3->GetParameter (3)) ;
  func7->SetParameter (4, 1) ;

  func7->SetParName (0, "centre") ;
  func7->SetParName (1, "shift") ;
  func7->SetParName (2, "scale") ;
  func7->SetParName (3, "slope_right") ;
  func7->SetParName (4, "power_left") ;

  delta->Fit ("func7", "+", "", 200, 2000) ;


  TCanvas * c6 = new TCanvas () ;
  relDiff->Draw ("hist") ;

  c6->Print ("relDiff.pdf", "pdf") ;


}                                                                                   
                                                                                    
