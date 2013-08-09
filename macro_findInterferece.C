w// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


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

  ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  ymin = h_MWW_mg->GetBinContent (h_MWW_mg->GetMinimumBin ()) ;
  TH1F * c4_mg_frame = (TH1F *) c4_mg->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_mg_frame->SetTitle (0) ;
  c4_mg_frame->SetStats (0) ;
  c4_mg_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_mg->Draw ("EPsame") ;

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
  func_ph_1->SetLineColor (kRed + 2) ;
  
  setParNamesdoubleGausCrystalBallLowHigh (func_ph_1) ;

  func_ph_1->SetParameter (0, 1.) ;                      // multiplicative scale
  func_ph_1->SetParameter (1, mass) ;                    // mean
  func_ph_1->SetParameter (2, gauss->GetParameter (2)) ; // gaussian sigma
  func_ph_1->SetParameter (3, 1) ;                       // right junction
  func_ph_1->SetParameter (4, 2) ;                       // right power law order
  func_ph_1->SetParameter (5, 1) ;                       // left junction
  func_ph_1->SetParameter (6, 2) ;                       // left power law order

  diff->Fit ("func_ph_1", "Q", "", 0.5 * mass - 50, 2 * mass) ;
  func_ph_1->SetParameters (func_ph_1->GetParameters ()) ;
  diff->Fit ("func_ph_1", "Q+L", "", 0.5 * mass - 50, 2 * mass) ;

  ymax = diff->GetBinContent (diff->GetMaximumBin ()) ;
  ymin = diff->GetBinContent (diff->GetMinimumBin ()) ;
  TH1F * c4_ph_frame = (TH1F *) c4_ph->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_ph_frame->SetTitle (0) ;
  c4_ph_frame->SetStats (0) ;
  c4_ph_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  diff->Draw ("EPsame") ;
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


  TF1 * f_doublePeakModel = new TF1 ("f_doublePeakModel", doublePeakModel, 0, 1000, 4) ;
  f_doublePeakModel->SetNpx (10000) ;
  f_doublePeakModel->SetLineWidth (1) ;
  f_doublePeakModel->SetLineColor (kRed + 2) ;

  f_doublePeakModel->SetParName (0, "scale") ;
  f_doublePeakModel->SetParName (1, "shift") ;
  f_doublePeakModel->SetParName (2, "distance") ;
  f_doublePeakModel->SetParName (3, "gamma") ; 

  f_doublePeakModel->SetParameter (0, 1.) ;   
  f_doublePeakModel->SetParameter (1, mass) ; 
  f_doublePeakModel->SetParameter (2, fabs (func_ph_1->GetParameter (1) - func_mg_1->GetParameter (1))) ;  
  double aveWidth = 0.5 * sqrt (
      func_ph_1->GetParameter (2) * func_ph_1->GetParameter (2) +
      func_mg_1->GetParameter (2) * func_mg_1->GetParameter (2)  
    ) ;
  f_doublePeakModel->SetParameter (3, 2 * aveWidth) ;  
  delta->Fit ("f_doublePeakModel", "Q+", "same", 0.5 * mass - 50, 2 * mass) ;

//  TF1 * func3 = new TF1 ("func3", doubleSlope, 0, 2000, 5) ;
//  func3->SetNpx (10000) ;
//  func3->SetLineWidth (1) ;
//  func3->SetLineColor (kGray + 2) ;
//
//  func3->SetParName (0, "centre") ;
//  func3->SetParName (1, "shift") ;
//  func3->SetParName (2, "scale") ;
//  func3->SetParName (3, "slope_right") ;
//  func3->SetParName (4, "slope_left") ;
//
//  //PG first fit, get the slope at the right with the mass hypothesis
//  func3->FixParameter (0, mass) ;
//  func3->FixParameter (1, 0.) ;
//  func3->SetParameter (2, 0.0005) ;
//  func3->SetParameter (3, 0.05) ;
//  func3->SetParameter (4, 0.05) ;
//  delta->Fit ("func3", "+", "same", mass, 2 * mass) ;
//
//  //PG second fit, get the centre given the slope
//  func3->SetParLimits (0, mass - sqrt (mass), mass + sqrt (mass)) ;
//  func3->FixParameter (1, 0.) ;
//  func3->SetParameter (2, 0.0005) ;
//  func3->FixParameter (3, func3->GetParameter (3)) ;
//  func3->SetParameter (4, 0.05) ;
//  delta->Fit ("func3", "+", "same", 0.5 * mass - 50, 2 * mass) ;
//
//  TF1 * func3_1 = new TF1 ("func3_1", singleSlope, 0, 2000, 5) ;
//  func3_1->SetNpx (10000) ;
//  func3_1->SetLineWidth (1) ;
//  func3_1->SetLineColor (kGreen + 2) ;
//  func3_1->SetParName (0, "centre") ;
//  func3_1->SetParName (1, "shift") ;
//  func3_1->SetParName (2, "scale") ;
//  func3_1->SetParName (3, "slope") ;
//
//  //PG third fit, fix centre and slope and get the rest
//  func3_1->SetParameter (0, func3->GetParameter (0)) ;
//  func3_1->FixParameter (1, 0.) ;
//  func3_1->SetParameter (2, 0.0005) ;
//  func3_1->FixParameter (3, func3->GetParameter (3)) ;
//
//  delta->Fit ("func3_1", "+", "same", 0.7 * mass, 2 * mass) ;

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
  outfile << "tg_par0->SetPoint (i++," << mass << ", " << f_doublePeakModel->GetParameter (0) << ") ;\n" ;
  outfile << "tg_par1->SetPoint (i++," << mass << ", " << f_doublePeakModel->GetParameter (1) << ") ;\n" ;
  outfile << "tg_par2->SetPoint (i++," << mass << ", " << f_doublePeakModel->GetParameter (2) << ") ;\n" ;
  outfile << "tg_par3->SetPoint (i++," << mass << ", " << f_doublePeakModel->GetParameter (3) << ") ;\n" ;
  outfile.close () ;

  return 0 ;










//PG (SBI - B) - S and the two signals

  TCanvas * c5_1 = new TCanvas () ;

  TF1 * f_relRatio = new TF1 ("f_relRatio", relativeCrystalBallLowHighRatio, 0, 2000, 14);
  f_relRatio->SetLineWidth (3) ;
  f_relRatio->SetLineColor (kGray + 2) ;
  Double_t f_ratio_pars [14] ;
  f_relRatio->SetParameters (f_ratio_pars) ;

  TH1F * relInterf = (TH1F *) delta->Clone ("relInterf") ;
  relInterf->Divide (h_MWW_mg) ;
  relInterf->Draw ("EP") ;
  f_relRatio->Draw ("same") ;
  
  c5_1->Print ("relRatio.pdf", "pdf") ;



  //PG (SBI - B) - S and S

  TCanvas * c5 = new TCanvas () ;
//  h_MWW_mg->Draw ("histsame") ;

  TF1 * retta = new TF1 ("retta", "pol1", 0, 2000) ;
  retta->SetLineWidth (1) ;
  retta->SetLineColor (kPink) ;
  delta->Fit ("retta", "+", "", mass - 30, mass + 30) ;

//  TF1 * expo = new TF1 ("expo", "-1 * exp ([0] * sqrt (fabs (x - [1])))", 0, 2000) ;
  TF1 * expo = new TF1 ("expo", "-1 * exp ([0] * fabs (x - [1]))", 0, 2000) ;
  expo->SetLineWidth (1) ;
  expo->SetLineColor (kPink) ;
  expo->FixParameter (1, mass) ;
//  expo->SetParameter (0, - 0.04) ;

  delta->Fit ("expo", "+", "", mass + 50, 2000) ;
  
  TF1 * func2 = new TF1 ("func2","-1 * [2] * (x - [0]) * exp (-1 * [3] * abs (x - [0])) + [1]",0, 2000) ;

  func2->SetLineWidth (1) ;
  func2->SetLineColor (kBlue) ;

  func2->FixParameter (0, mass) ;
  func2->FixParameter (1, 0.) ;

  func2->SetParameter (2, retta->GetParameter (1)) ;
  func2->SetParameter (3, expo->GetParameter (0)) ;

//  delta->Fit ("func2", "+", "", 600, 2000) ;
//  expo->Draw ("same") ;

  c5->Print ("delta.pdf", "pdf") ;
  expo->Draw ("same") ;


  //PG S / (SBI - S)

  TF1 * f_ratio = new TF1 ("f_ratio", crystalBallLowHighRatio, 0, 2000, 14);
  f_ratio->SetLineWidth (3) ;
  f_ratio->SetLineColor (kGray + 2) ;
  Double_t f_ratio_pars [14] ;
  func_ph_1->GetParameters (f_ratio_pars) ;
  func_mg_1->GetParameters (f_ratio_pars + 7) ;
  f_ratio->SetParameters (f_ratio_pars) ;

  TCanvas * c3 = new TCanvas ("c3") ;
  c3->DrawFrame (200, -0.4, rangeScale * mass, 12) ;
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
                                                                                    
