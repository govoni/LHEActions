/*
r00t -q macro_findInterferece.C\(\"findInterference.350.root\",350\)
r00t -q macro_findInterferece.C\(\"findInterference.500.root\",500\)
r00t -q macro_findInterferece.C\(\"findInterference.650.root\",650\)
r00t -q macro_findInterferece.C\(\"findInterference.800.root\",800\)
r00t -q macro_findInterferece.C\(\"findInterference.1000.root\",1000\)
*/


double max (double one, double two)
{
  if (one > two) return one ;
  return two ;
}


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
  func->SetParName (3, "alphaR") ;
  func->SetParName (4, "nL") ;
  func->SetParName (5, "alphaL") ;
  func->SetParName (6, "nR") ;
  return ;
}  


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t crystalBallLowHighRatio (Double_t * xx, Double_t * par) // (SBI - B) / S
{
  return crystalBallLowHigh (xx, par) / crystalBallLowHigh (xx, par + 7) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t relativeCrystalBallLowHighRatio (Double_t * xx, Double_t * par) // [(SBI - B) - S] / S
{
  return (crystalBallLowHigh (xx, par) - crystalBallLowHigh (xx, par + 7)) / crystalBallLowHigh (xx, par + 7) ;
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


Double_t doublePeakModel (Double_t * xx, Double_t * par)
{
  double scale    = par[0] ;
  double shift    = par[1] ;
  double distance = par[2] ;
  double gamma    = par[3] ;
  double x = xx[0] - shift ;

  double norm = 1. / (shift * shift + gamma) - 1 / ((shift + 2 * distance) * (shift + 2 * distance) + gamma) ;
  return scale * (1. / norm) * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;
  return scale * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t singleSlope (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double slope_left  = par[3] ;
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


int findBin (TH1F * h, double val)
{
  int i = 1 ;
  for ( ; i <= h->GetNbinsX () ; ++i)
    {
      if (h->GetBinLowEdge (i) > val) break ;
    } 
  return i - 1 ;  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int macro_findInterferece (string filename, double mass)                                                        
{        

  gSystem->Load ("Functions.cc") ;
  gStyle->SetPadTopMargin (0.1) ;

  TFile * f = new TFile (filename.c_str ()) ;
  TH1F * h_MWW_phbkgsig = (TH1F *) f->Get ("h_MWW_phbkgsig") ;
  TH1F * h_MWW_phbkg    = (TH1F *) f->Get ("h_MWW_phbkg") ;
  TH1F * h_MWW_mg       = (TH1F *) f->Get ("h_MWW_mg") ;

  int reBin = 1 ;
  if (mass > 480) reBin = 2 ;
  if (mass > 810) reBin = 6 ;
  double rangeScale = 1.5 ;
  if (mass > 480) rangeScale = 2 ;
  
  int scaling = 1 ;
  
  h_MWW_phbkgsig->Rebin (reBin) ;
  h_MWW_phbkg   ->Rebin (reBin) ;
  h_MWW_mg      ->Rebin (reBin) ;

  double ymax= 2000. ;
  double ymin= 0. ;

  //PG graphics
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  h_MWW_phbkg->SetStats (0) ;
  h_MWW_phbkg->SetTitle ("") ;
  h_MWW_phbkg->SetLineColor (kOrange) ;
  h_MWW_phbkg->SetLineWidth (2) ;
  h_MWW_phbkg->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_phbkgsig->SetStats (0) ;
  h_MWW_phbkgsig->SetTitle ("") ;
  h_MWW_phbkgsig->SetLineColor (kRed) ;
  h_MWW_phbkgsig->SetLineWidth (2) ;
  h_MWW_phbkgsig->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_mg->SetStats (0) ;
  h_MWW_mg->SetTitle ("") ;
  h_MWW_mg->SetLineColor (kBlue) ;
  h_MWW_mg->SetLineWidth (2) ;
  h_MWW_mg->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;

  //PG histograms operations
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  //PG SBI - B
  TH1F * tempoDiff = (TH1F *) h_MWW_phbkgsig->Clone ("tempoDiff") ;
  tempoDiff->Add (h_MWW_phbkg, -1) ;
  tempoDiff->SetLineColor (kGreen + 1) ;

  if (scaling > 0)
    {
      double maxSI = tempoDiff->GetBinContent (tempoDiff->GetMaximumBin ()) ;
      double maxS = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
      h_MWW_mg->Scale (maxSI / maxS) ;
    }

  //PG SBI - B
  TH1F * diff = (TH1F *) h_MWW_phbkgsig->Clone ("diff") ;
  diff->Add (h_MWW_phbkg, -1) ;
  diff->SetLineColor (kGreen + 1) ;

  //PG (SBI - B) / S
  TH1F * ratio = (TH1F *) diff->Clone ("ratio") ;
  ratio->Divide (h_MWW_mg) ;

  //PG (SBI - B) - S
  TH1F * delta = (TH1F *) diff->Clone ("delta") ;
  delta->SetTitle ("") ;
  delta->Add (h_MWW_mg, -1) ;
  delta->SetLineColor (kViolet + 1) ;

  //PG ((SBI - B) - S) / S
  TH1F * relDiff = delta->Clone ("relDiff") ;
  relDiff->Divide (h_MWW_mg) ;



  //PG plotting
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  TString suffix = "." ; 
  suffix += mass ;
  suffix += ".pdf" ;
  
  //PG initial spectra ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TCanvas * c1 = new TCanvas ("c1", "c1") ;
  ymax = h_MWW_phbkgsig->GetBinContent (h_MWW_phbkgsig->GetMaximumBin ()) ;
  ymin = h_MWW_mg->GetBinContent (findBin (h_MWW_mg, mass + exp (mass * 0.0058461 + 0.65385))) ;
  TH1F * c1_frame = (TH1F *) c1->DrawFrame (200, 0.2 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c1_frame->SetTitle (0) ;
  c1_frame->SetStats (0) ;
  c1->SetLogy () ;
  c1_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_phbkg->Draw ("histsame") ;
  h_MWW_phbkgsig->Draw ("histsame") ;
  h_MWW_mg->Draw ("histsame") ;

  c1_leg = new TLegend (0.6,0.8,0.9,0.95) ;
  c1_leg->SetFillStyle (0) ;
  c1_leg->SetBorderSize (0) ;
  c1_leg->SetTextFont (42) ;
  c1_leg->AddEntry (h_MWW_phbkg, "B","l") ;
  c1_leg->AddEntry (h_MWW_phbkgsig, "SBI","l") ;
  c1_leg->AddEntry (h_MWW_mg, "S","l") ;
  c1_leg->Draw () ;
  
  c1->Print (TString ("spectra") + suffix, "pdf") ;
  
  //PG S only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  TCanvas * c4_mg = new TCanvas ("c4_mg", "c4_mg") ;

  TF1 * func_mg_1 = new TF1 ("func_mg_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_mg_1->SetNpx (10000) ;
  func_mg_1->SetLineWidth (1) ;
  func_mg_1->SetLineColor (kBlue + 1) ;

  setParNamesdoubleGausCrystalBallLowHigh (func_mg_1) ;
  cout << "PREFIT RMS : " << h_MWW_mg->GetRMS () << endl ;

  func_mg_1->SetParameter (0, 1.) ;                  // multiplicative scale
  func_mg_1->SetParameter (1, mass) ;                // mean
  func_mg_1->SetParameter (2, h_MWW_mg->GetRMS ()) ; // gaussian sigma
  func_mg_1->SetParLimits (2, 0.1 * h_MWW_mg->GetRMS (), 20 * h_MWW_mg->GetRMS ()) ;
  func_mg_1->SetParameter (3, 1) ;                   // right junction
  func_mg_1->SetParLimits (3, 0.1, 5) ;              // right junction
  func_mg_1->SetParameter (4, 1) ;                   // right power law order
  func_mg_1->SetParameter (5, 1) ;                   // left junction
  func_mg_1->SetParLimits (5, 0.1, 5) ;              // left junction
  func_mg_1->SetParameter (6, 1) ;                   // left power law order

  int sign = 1 ;
  if (mass < 400) sign = -2 ;
  cout << "-------------------\nFITTING THE MADGRAPH SIGNAL\n" ;
  h_MWW_mg->Fit ("func_mg_1", "+", "", 0.5 * mass + sign * 50, 2 * mass) ;
  cout << "CHI2 / NDOF = " << func_mg_1->GetChisquare () /func_mg_1->GetNDF () << endl ;
  func_mg_1->SetParameters (func_mg_1->GetParameters ()) ;
  func_mg_1->SetLineColor (kBlue + 3) ;
  cout << "-------------------\nFITTING THE MADGRAPH SIGNAL W/ LIKELIHOOD\n" ;
  h_MWW_mg->Fit ("func_mg_1", "+L", "", 0.5 * mass - 50, 2 * mass) ;
  cout << "CHI2 / NDOF = " << func_mg_1->GetChisquare () /func_mg_1->GetNDF () << endl ;

  ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  ymin = h_MWW_mg->GetBinContent (h_MWW_mg->GetMinimumBin ()) ;
  TH1F * c4_mg_frame = (TH1F *) c4_mg->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_mg_frame->SetTitle (0) ;
  c4_mg_frame->SetStats (0) ;
  c4_mg_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_mg->Draw ("EPsame") ;

  double rightTh = fabs (func_mg_1->GetParameter (3)) * func_mg_1->GetParameter (2) + func_mg_1->GetParameter (1) ;
  cout << "MG RIGHT THRESHOLD " << rightTh << endl ;
  double leftTh  = -1 * fabs (func_mg_1->GetParameter (5)) * func_mg_1->GetParameter (2) + func_mg_1->GetParameter (1) ;
  cout << "MG LEFT THRESHOLD " << leftTh << endl ;

  TLine * l_rightTh = new TLine (rightTh, 0.9 * ymin, rightTh, 1.1 * ymax) ;
  l_rightTh->SetLineColor (kRed) ;
  l_rightTh->Draw ("same") ;
  TLine * l_leftTh = new TLine (leftTh, 0.9 * ymin, leftTh, 1.1 * ymax) ;
  l_leftTh->SetLineColor (kRed) ;
  l_leftTh->Draw ("same") ;

  c4_mg->Update () ;
  c4_mg->Print (TString ("signals_mg") + suffix, "pdf") ;

  //PG (SBI - B) only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TCanvas * c4_ph = new TCanvas ("c4_ph", "c4_ph") ;

  TF1 * gauss = new TF1 ("gauss", "gaus", 0, 2000) ;
  gauss->SetLineWidth (1) ;
  gauss->SetLineColor (kGray + 2) ;
  diff->Fit ("gauss", "Q+", "", 0.75 * mass, mass * 1.25) ;

  TF1 * func_ph_1 = new TF1 ("func_ph_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_ph_1->SetNpx (10000) ;
  func_ph_1->SetLineWidth (1) ;
  func_ph_1->SetLineColor (kRed + 1) ;
  
  setParNamesdoubleGausCrystalBallLowHigh (func_ph_1) ;

  func_ph_1->SetParameter (0, 1.) ;                      // multiplicative scale
  func_ph_1->SetParameter (1, mass) ;                    // mean
  func_ph_1->SetParameter (2, gauss->GetParameter (2)) ; // gaussian sigma
//  func_ph_1->SetParLimits (2, 0.1 * gauss->GetParameter (2), 20 * gauss->GetParameter (2)) ;
  func_ph_1->SetParameter (3, 1) ;                       // right junction
//  func_ph_1->SetParLimits (3, 0.1, 5) ;                  // right junction
  func_ph_1->SetParameter (4, 1) ;                       // right power law order
  func_ph_1->SetParameter (5, 1) ;                       // left junction
//  func_ph_1->SetParLimits (5, 0.1, 5) ;                  // left junction
  func_ph_1->SetParameter (6, 1) ;                       // left power law order

  cout << "-------------------\nFITTING THE PHANTOM SIGNAL\n" ;
  diff->Fit ("func_ph_1", "", "", 0.5 * mass - 50, 2 * mass) ;
  cout << "CHI2 / NDOF = " << func_ph_1->GetChisquare () /func_ph_1->GetNDF () << endl ;
  func_ph_1->SetParameters (func_ph_1->GetParameters ()) ;
  func_ph_1->SetLineColor (kRed + 3) ;
  cout << "-------------------\nFITTING THE PHANTOM SIGNAL W/ LIKELIHOOD\n" ;
  diff->Fit ("func_ph_1", "+L", "", 0.5 * mass - 50, 2 * mass) ;
  cout << "CHI2 / NDOF = " << func_ph_1->GetChisquare () /func_ph_1->GetNDF () << endl ;

  ymax = diff->GetBinContent (diff->GetMaximumBin ()) ;
  ymax = max (ymax, func_ph_1->GetMaximum ()) ;
  ymin = diff->GetBinContent (diff->GetMinimumBin ()) ;
  ymin = max (ymin, -0.1 * ymax) ;
  TH1F * c4_ph_frame = (TH1F *) c4_ph->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_ph_frame->SetTitle (0) ;
  c4_ph_frame->SetStats (0) ;
  c4_ph_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  diff->Draw ("EPsame") ;

  double rightTh_ph = fabs (func_ph_1->GetParameter (3)) * func_ph_1->GetParameter (2) + func_ph_1->GetParameter (1) ;
  cout << "PH RIGHT THRESHOLD " << rightTh_ph << endl ;
  double leftTh_ph  = -1 * fabs (func_ph_1->GetParameter (5)) * func_ph_1->GetParameter (2) + func_ph_1->GetParameter (1) ;
  cout << "PH LEFT THRESHOLD " << leftTh_ph << endl ;

  TLine * l_rightTh_ph = new TLine (rightTh_ph, 0.9 * ymin, rightTh_ph, 1.1 * ymax) ;
  l_rightTh_ph->SetLineColor (kRed) ;
  l_rightTh_ph->Draw ("same") ;
  TLine * l_leftTh_ph = new TLine (leftTh_ph, 0.9 * ymin, leftTh_ph, 1.1 * ymax) ;
  l_leftTh_ph->SetLineColor (kRed) ;
  l_leftTh_ph->Draw ("same") ;

  c4_ph->Print (TString ("signals_ph") + suffix, "pdf") ;

  //PG (SBI - B) - S only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TCanvas * c3 = new TCanvas ("c3", "c3") ;
  double ymax = delta->GetBinContent (delta->GetMaximumBin ()) ;
  double ymin = delta->GetBinContent (delta->GetMinimumBin ()) ;
  if (ymin > 0) ymin *= 0.9 ;
  else          ymin *= 1.5 ;
  if (ymin < -2 * ymax) ymin = -2 * ymax ;
  TH1F * c3_frame = (TH1F *) c3->DrawFrame (200, ymin, rangeScale * mass, 1.1 * ymax) ;
  c3_frame->SetTitle (0) ;
  c3_frame->SetStats (0) ;
  c3_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;

  TF1 * f_doublePeakModel = new TF1 ("f_doublePeakModel", doublePeakModel, 0, 2000, 4) ;
  f_doublePeakModel->SetNpx (10000) ;
  f_doublePeakModel->SetLineWidth (1) ;
  f_doublePeakModel->SetLineColor (kRed + 2) ;

  f_doublePeakModel->SetParName (0, "scale") ;
  f_doublePeakModel->SetParName (1, "shift") ;
  f_doublePeakModel->SetParName (2, "distance") ;
  f_doublePeakModel->SetParName (3, "gamma") ; 

  f_doublePeakModel->SetParameter (0, -0.000001) ;   
  f_doublePeakModel->SetParameter (1, mass) ; 
  f_doublePeakModel->SetParameter (2, -0.01) ;  
//  f_doublePeakModel->SetParameter (2, fabs (func_ph_1->GetParameter (1) - func_mg_1->GetParameter (1))) ;  
  double aveWidth = 0.5 * sqrt (
      func_ph_1->GetParameter (2) * func_ph_1->GetParameter (2) +
      func_mg_1->GetParameter (2) * func_mg_1->GetParameter (2)  
    ) ;
  f_doublePeakModel->SetParameter (3, mass * mass * 0.25 * 0.25) ;  
//  f_doublePeakModel->SetParameter (3, 2 * aveWidth) ;  
  delta->Fit ("f_doublePeakModel", "+", "same", 0.5 * mass - 50, 2 * mass) ;

  delta->Draw ("histsame") ;
  c3_leg = new TLegend (0.5,0.8,0.9,0.95) ;
  c3_leg->SetFillStyle (0) ;
  c3_leg->SetBorderSize (0) ;
  c3_leg->SetTextFont (42) ;
  c3_leg->AddEntry (delta, "(SBI - B) - S","l") ;
  
  c3_leg->Draw () ;

  c3->Print (TString ("interf") + suffix, "pdf") ;


  //PG S only, and (SBI - B) ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TCanvas * c4 = new TCanvas ("c4", "c4") ;
  double ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  double ymax2 = diff->GetBinContent (diff->GetMaximumBin ()) ;
  if (ymax2 > ymax) ymax = ymax2 ;
  double ymin = delta->GetBinContent (delta->GetMinimumBin ()) ;
  if (ymin > 0) ymin *= 0.9 ;
  else          ymin *= 1.5 ;
  if (ymin < -1.1 * ymax) ymin = -1.1 * ymax ;
  TH1F * c4_frame = (TH1F *) c4->DrawFrame (200, ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_frame->SetTitle (0) ;
  c4_frame->SetStats (0) ;
  c4_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_mg->Draw ("histsame") ;
  diff->Draw ("histsame") ;
  delta->Draw ("histsame") ;
  c4_leg = new TLegend (0.5,0.8,0.9,0.95) ;
  c4_leg->SetFillStyle (0) ;
  c4_leg->SetBorderSize (0) ;
  c4_leg->SetTextFont (42) ;
  c4_leg->AddEntry (diff, "SBI - B","l") ;
  c4_leg->AddEntry (h_MWW_mg, "S","l") ;
  c4_leg->AddEntry (delta, "(SBI - B) - S","l") ;
  
  c4_leg->Draw () ;

  c4->Print (TString ("signals") + suffix, "pdf") ;

  //PG output of the fitting function parameters

  std::ofstream outfile;

  outfile.open ("graphs.txt", std::ios_base::app) ;
  outfile << "  \n  // ----> MASS " << mass << " ---- ---- ---- \n\n" ;
  outfile << "  // interference parametrisation:\n" ;
  outfile << "  tg_par0->SetPoint (i, " << mass << ", " << f_doublePeakModel->GetParameter (0) << ") ;\n" ;
  outfile << "  tg_par1->SetPoint (i, " << mass << ", " << f_doublePeakModel->GetParameter (1) << ") ;\n" ;
  outfile << "  tg_par2->SetPoint (i, " << mass << ", " << f_doublePeakModel->GetParameter (2) << ") ;\n" ;
  outfile << "  tg_par3->SetPoint (i, " << mass << ", " << f_doublePeakModel->GetParameter (3) << ") ;\n" ;
  outfile << "  TF1 * func_" << mass << " = new TF1 (\"func_" << mass << "\",doublePeakModel, 200, 2000, 4) ;\n" ; 
  outfile << "  double params_" << mass << "[4] = {" << f_doublePeakModel->GetParameter (0) << ", " << f_doublePeakModel->GetParameter (1) << ", " << f_doublePeakModel->GetParameter (2) << ", " << f_doublePeakModel->GetParameter (3) << " } ;\n" ;
  outfile << "  func_" << mass << "->SetParameters (params_" << mass << ") ;\n" ; 
  outfile << "  // MG signal only parametrisation:\n" ;
  outfile << "  tg_sig_par0->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (0) << ") ;\n" ;
  outfile << "  tg_sig_par1->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (1) << ") ;\n" ;
  outfile << "  tg_sig_par2->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (2) << ") ;\n" ;
  outfile << "  tg_sig_par3->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (3) << ") ;\n" ;
  outfile << "  tg_sig_par4->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (4) << ") ;\n" ;
  outfile << "  tg_sig_par5->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (5) << ") ;\n" ;
  outfile << "  tg_sig_par6->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (6) << ") ;\n" ;
  outfile << "  TF1 * func_sig_" << mass << " = new TF1 (\"func_sig_" << mass << "\",crystalBallLowHigh, 200, 2000, 7) ;\n" ; 
  outfile << "  double params_sig_" << mass << "[7] = {" << func_mg_1->GetParameter (0) << ", " << func_mg_1->GetParameter (1) << ", " << func_mg_1->GetParameter (2) << ", " << func_mg_1->GetParameter (3) << ", " << func_mg_1->GetParameter (4) << ", " << func_mg_1->GetParameter (5) << ", " << func_mg_1->GetParameter (6)  << " } ;\n" ;
  outfile << "  func_sig_" << mass << "->SetParameters (params_sig_" << mass << ") ;\n" ; 
  outfile << "  // PG SBI - B  parametrisation:\n" ;
  outfile << "  tg_sAi_par0->SetPoint (i, " << mass << ", " << func_ph_1->GetParameter (0) << ") ;\n" ;
  outfile << "  tg_sAi_par1->SetPoint (i, " << mass << ", " << func_ph_1->GetParameter (1) << ") ;\n" ;
  outfile << "  tg_sAi_par2->SetPoint (i, " << mass << ", " << func_ph_1->GetParameter (2) << ") ;\n" ;
  outfile << "  tg_sAi_par3->SetPoint (i, " << mass << ", " << func_ph_1->GetParameter (3) << ") ;\n" ;
  outfile << "  tg_sAi_par4->SetPoint (i, " << mass << ", " << func_ph_1->GetParameter (4) << ") ;\n" ;
  outfile << "  tg_sAi_par5->SetPoint (i, " << mass << ", " << func_ph_1->GetParameter (5) << ") ;\n" ;
  outfile << "  tg_sAi_par6->SetPoint (i, " << mass << ", " << func_ph_1->GetParameter (6) << ") ;\n" ;
  outfile << "  TF1 * func_sAi_" << mass << " = new TF1 (\"func_sAi_" << mass << "\",crystalBallLowHigh, 200, 2000, 7) ;\n" ; 
  outfile << "  double params_sAi_" << mass << "[7] = {" << func_ph_1->GetParameter (0) << ", " << func_ph_1->GetParameter (1) << ", " << func_ph_1->GetParameter (2) << ", " << func_ph_1->GetParameter (3) << ", " << func_ph_1->GetParameter (4) << ", " << func_ph_1->GetParameter (5) << ", " << func_ph_1->GetParameter (6)  << " } ;\n" ;
  outfile << "  func_sAi_" << mass << "->SetParameters (params_sAi_" << mass << ") ;\n" ; 
  outfile << "  i++ ;\n" ;

  outfile.close () ;

  return 0 ;


}                                                                                   
                                                                                    
