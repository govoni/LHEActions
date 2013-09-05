

/*** double crystall ball ***/
double crystalBallLowHigh (double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  //[5] = alpha2
  //[6] = n2
  
  double xx = x[0];
  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];
  double alpha2 = par[5];
  double n2 = par[6];

  if( (xx-mean)/sigma > fabs(alpha) )
  {
    double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
    double B = n/fabs(alpha) - fabs(alpha);
    
    return par[0] * A * pow(B + (xx-mean)/sigma, -1.*n);
  }
  
  else if( (xx-mean)/sigma < -1.*fabs(alpha2) )
  {
    double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
    double B = n2/fabs(alpha2) - fabs(alpha2);
    
    return par[0] * A * pow(B - (xx-mean)/sigma, -1.*n2);
  }
  
  else
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
  }
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


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


Double_t intereferenceOverSignal (Double_t * xx, Double_t * par) 
{
  if (crystalBallLowHigh (xx, par + 4) < 0.0000001) return 0 ;
  return doublePeakModel (xx, par) / crystalBallLowHigh (xx, par + 4) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int graphs () 
{
  TGraph * tg_par0 = new TGraph (5) ;
  TGraph * tg_par1 = new TGraph (5) ;
  TGraph * tg_par2 = new TGraph (5) ;
  TGraph * tg_par3 = new TGraph (5) ;
  
  TGraph * tg_sig_par0 = new TGraph (5) ;
  TGraph * tg_sig_par1 = new TGraph (5) ;
  TGraph * tg_sig_par2 = new TGraph (5) ;
  TGraph * tg_sig_par3 = new TGraph (5) ;
  TGraph * tg_sig_par4 = new TGraph (5) ;
  TGraph * tg_sig_par5 = new TGraph (5) ;
  TGraph * tg_sig_par6 = new TGraph (5) ;

  int i = 0 ; 
  
  
  
  // ----> MASS 350 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 350, -4.4362e-07) ;
  tg_par1->SetPoint (i, 350, 349.521) ;
  tg_par2->SetPoint (i, 350, -9.31673e-05) ;
  tg_par3->SetPoint (i, 350, 660.315) ;
  TF1 * func_350 = new TF1 ("func_350",doublePeakModel, 200, 2000, 4) ;
  double params_350[4] = {-4.4362e-07, 349.521, -9.31673e-05, 660.315 } ;
  func_350->SetParameters (params_350) ;
  // signal only parametrisation:
  tg_sig_par0->SetPoint (i, 350, 0.00794395) ;
  tg_sig_par1->SetPoint (i, 350, 350.343) ;
  tg_sig_par2->SetPoint (i, 350, 7.63757) ;
  tg_sig_par3->SetPoint (i, 350, 1.15603) ;
  tg_sig_par4->SetPoint (i, 350, 1.55177) ;
  tg_sig_par5->SetPoint (i, 350, 1.26848) ;
  tg_sig_par6->SetPoint (i, 350, 2.15573) ;
  TF1 * func_sig_350 = new TF1 ("func_sig_350",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_350[7] = {0.00794395, 350.343, 7.63757, 1.15603, 1.55177, 1.26848, 2.15573 } ;
  func_sig_350->SetParameters (params_sig_350) ;
  i++ ;
  
  // ----> MASS 500 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 500, -4.30542e-06) ;
  tg_par1->SetPoint (i, 500, 488.467) ;
  tg_par2->SetPoint (i, 500, -0.00537915) ;
  tg_par3->SetPoint (i, 500, 9775.1) ;
  TF1 * func_500 = new TF1 ("func_500",doublePeakModel, 200, 2000, 4) ;
  double params_500[4] = {-4.30542e-06, 488.467, -0.00537915, 9775.1 } ;
  func_500->SetParameters (params_500) ;
  // signal only parametrisation:
  tg_sig_par0->SetPoint (i, 500, 0.00142874) ;
  tg_sig_par1->SetPoint (i, 500, 503.368) ;
  tg_sig_par2->SetPoint (i, 500, 29.8165) ;
  tg_sig_par3->SetPoint (i, 500, 0.831729) ;
  tg_sig_par4->SetPoint (i, 500, 2.54529) ;
  tg_sig_par5->SetPoint (i, 500, 1.24525) ;
  tg_sig_par6->SetPoint (i, 500, 2.78472) ;
  TF1 * func_sig_500 = new TF1 ("func_sig_500",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_500[7] = {0.00142874, 503.368, 29.8165, 0.831729, 2.54529, 1.24525, 2.78472 } ;
  func_sig_500->SetParameters (params_sig_500) ;
  i++ ;
  
  // ----> MASS 650 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 650, -2.15082e-06) ;
  tg_par1->SetPoint (i, 650, 648.512) ;
  tg_par2->SetPoint (i, 650, -0.00935512) ;
  tg_par3->SetPoint (i, 650, 23441.2) ;
  TF1 * func_650 = new TF1 ("func_650",doublePeakModel, 200, 2000, 4) ;
  double params_650[4] = {-2.15082e-06, 648.512, -0.00935512, 23441.2 } ;
  func_650->SetParameters (params_650) ;
  // signal only parametrisation:
  tg_sig_par0->SetPoint (i, 650, 0.000276903) ;
  tg_sig_par1->SetPoint (i, 650, 666.593) ;
  tg_sig_par2->SetPoint (i, 650, 79.0964) ;
  tg_sig_par3->SetPoint (i, 650, 1.53166) ;
  tg_sig_par4->SetPoint (i, 650, 1.69832) ;
  tg_sig_par5->SetPoint (i, 650, 1.57689) ;
  tg_sig_par6->SetPoint (i, 650, 2.16479) ;
  TF1 * func_sig_650 = new TF1 ("func_sig_650",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_650[7] = {0.000276903, 666.593, 79.0964, 1.53166, 1.69832, 1.57689, 2.16479 } ;
  func_sig_650->SetParameters (params_sig_650) ;
  i++ ;
  
  // ----> MASS 800 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 800, -7.4816e-07) ;
  tg_par1->SetPoint (i, 800, 759.208) ;
  tg_par2->SetPoint (i, 800, 93.1791) ;
  tg_par3->SetPoint (i, 800, 27559.3) ;
  TF1 * func_800 = new TF1 ("func_800",doublePeakModel, 200, 2000, 4) ;
  double params_800[4] = {-7.4816e-07, 759.208, 93.1791, 27559.3 } ;
  func_800->SetParameters (params_800) ;
  // signal only parametrisation:
  tg_sig_par0->SetPoint (i, 800, 8.16552e-05) ;
  tg_sig_par1->SetPoint (i, 800, 803.873) ;
  tg_sig_par2->SetPoint (i, 800, 127.14) ;
  tg_sig_par3->SetPoint (i, 800, 1.63088) ;
  tg_sig_par4->SetPoint (i, 800, 1.77215) ;
  tg_sig_par5->SetPoint (i, 800, 1.58901) ;
  tg_sig_par6->SetPoint (i, 800, 3.80968) ;
  TF1 * func_sig_800 = new TF1 ("func_sig_800",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_800[7] = {8.16552e-05, 803.873, 127.14, 1.63088, 1.77215, 1.58901, 3.80968 } ;
  func_sig_800->SetParameters (params_sig_800) ;
  i++ ;
  
  // ----> MASS 1000 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 1000, -1.98075e-06) ;
  tg_par1->SetPoint (i, 1000, 797.034) ;
  tg_par2->SetPoint (i, 1000, 211.982) ;
  tg_par3->SetPoint (i, 1000, 50378) ;
  TF1 * func_1000 = new TF1 ("func_1000",doublePeakModel, 200, 2000, 4) ;
  double params_1000[4] = {-1.98075e-06, 797.034, 211.982, 50378 } ;
  func_1000->SetParameters (params_1000) ;
  // signal only parametrisation:
  tg_sig_par0->SetPoint (i, 1000, 6.33517e-05) ;
  tg_sig_par1->SetPoint (i, 1000, 937.791) ;
  tg_sig_par2->SetPoint (i, 1000, 229.464) ;
  tg_sig_par3->SetPoint (i, 1000, 1.786) ;
  tg_sig_par4->SetPoint (i, 1000, 2.05718) ;
  tg_sig_par5->SetPoint (i, 1000, 1.84084) ;
  tg_sig_par6->SetPoint (i, 1000, 28.6439) ;
  TF1 * func_sig_1000 = new TF1 ("func_sig_1000",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_1000[7] = {6.33517e-05, 937.791, 229.464, 1.786, 2.05718, 1.84084, 28.6439 } ;
  func_sig_1000->SetParameters (params_sig_1000) ;
  i++ ;


  //PG plotting
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  tg_par0->SetTitle ("") ;
  tg_par1->SetTitle ("") ;
  tg_par2->SetTitle ("") ;
  tg_par3->SetTitle ("") ;
  
  tg_par0->GetYaxis ()->SetTitle ("scale") ;
  tg_par1->GetYaxis ()->SetTitle ("shift") ;
  tg_par2->GetYaxis ()->SetTitle ("distance") ;
  tg_par3->GetYaxis ()->SetTitle ("gamma") ;
  
  tg_par0->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_par1->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_par2->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_par3->GetXaxis ()->SetTitle ("m_{WW}") ;

  TCanvas * c_par = new TCanvas () ;
  c_par->Divide (2,2) ;
  i = 0 ;
  c_par->cd (++i) ; tg_par0->Draw ("AL*") ; 
  c_par->cd (++i) ; tg_par1->Draw ("AL*") ; // tg_mass->Draw ("L") ;
  c_par->cd (++i) ; tg_par2->Draw ("AL*") ; // tg_width->Draw ("L") ;
  c_par->cd (++i) ; tg_par3->Draw ("AL*") ; // tg_width->Draw ("L") ;
  
  TCanvas * c_sig_par = new TCanvas () ;
  c_sig_par->Divide (3,2) ;
  i = 0 ;
  c_sig_par->cd (++i) ; tg_sig_par0->Draw ("AL*") ; 
  c_sig_par->cd (++i) ; tg_sig_par1->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par2->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par3->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par4->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par5->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par6->Draw ("AL*") ;

  TCanvas * c_results = new TCanvas ("c_results", "c_results", 5000, 600) ;
  c_results->Divide (5,2) ;

  i = 1 ;
  double mass = 350 ;
  double rangeScale = 1.5 ;

  // 350 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 350 ;
  rangeScale = 1.5 ;
  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_350_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_350_pars[j] = params_350[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_350_pars[j + 4] = params_sig_350[j] ;
  TF1 * f_IOS_350 = new TF1 ("f_IOS_350", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_350->SetParameters (IOS_350_pars) ;

  c_results->cd (i) ;
  gPad->DrawFrame (200, 1.1 * func_350->GetMinimum (), rangeScale * mass, 1.1 * func_sig_350->GetMaximum ()) ;
  func_sig_350->SetLineWidth (1) ;
  func_sig_350->SetNpx (10000) ;
  func_sig_350->Draw ("same") ;
  func_350->SetLineWidth (1) ;
  func_350->SetNpx (10000) ;
  func_350->Draw ("same") ;
  c_results->cd (5 + i++) ;
  gPad->DrawFrame (200, 1.1 * f_IOS_350->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_350->GetMaximum ()) ;
  f_IOS_350->SetLineWidth (1) ;
  f_IOS_350->SetNpx (10000) ;
  f_IOS_350->Draw ("same") ;
 


  // 500 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 500 ;
  rangeScale = 1.5 ;
//  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_500_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_500_pars[j] = params_500[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_500_pars[j + 4] = params_sig_500[j] ;
  TF1 * f_IOS_500 = new TF1 ("f_IOS_500", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_500->SetParameters (IOS_500_pars) ;

  c_results->cd (i) ;
  gPad->DrawFrame (200, 1.1 * func_500->GetMinimum (), rangeScale * mass, 1.1 * func_sig_500->GetMaximum ()) ;
  func_sig_500->SetLineWidth (1) ;
  func_sig_500->SetNpx (10000) ;
  func_sig_500->Draw ("same") ;
  func_500->SetLineWidth (1) ;
  func_500->SetNpx (10000) ;
  func_500->Draw ("same") ;
  c_results->cd (5 + i++) ;
  gPad->DrawFrame (200, 1.1 * f_IOS_500->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_500->GetMaximum ()) ;
  f_IOS_500->SetLineWidth (1) ;
  f_IOS_500->SetNpx (10000) ;
  f_IOS_500->Draw ("same") ;
 


  // 650 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 650 ;
  rangeScale = 1.5 ;
//  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_650_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_650_pars[j] = params_650[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_650_pars[j + 4] = params_sig_650[j] ;
  TF1 * f_IOS_650 = new TF1 ("f_IOS_650", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_650->SetParameters (IOS_650_pars) ;

  c_results->cd (i) ;
  gPad->DrawFrame (200, 1.1 * func_650->GetMinimum (), rangeScale * mass, 1.1 * func_sig_650->GetMaximum ()) ;
  func_sig_650->SetLineWidth (1) ;
  func_sig_650->SetNpx (10000) ;
  func_sig_650->Draw ("same") ;
  func_650->SetLineWidth (1) ;
  func_650->SetNpx (10000) ;
  func_650->Draw ("same") ;
  c_results->cd (5 + i++) ;
  gPad->DrawFrame (200, 1.1 * f_IOS_650->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_650->GetMaximum ()) ;
  f_IOS_650->SetLineWidth (1) ;
  f_IOS_650->SetNpx (10000) ;
  f_IOS_650->Draw ("same") ;
 


  // 800 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 800 ;
  rangeScale = 1.5 ;
//  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_800_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_800_pars[j] = params_800[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_800_pars[j + 4] = params_sig_800[j] ;
  TF1 * f_IOS_800 = new TF1 ("f_IOS_800", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_800->SetParameters (IOS_800_pars) ;

  c_results->cd (i) ;
  gPad->DrawFrame (200, 1.1 * func_800->GetMinimum (), rangeScale * mass, 1.1 * func_sig_800->GetMaximum ()) ;
  func_sig_800->SetLineWidth (1) ;
  func_sig_800->SetNpx (10000) ;
  func_sig_800->Draw ("same") ;
  func_800->SetLineWidth (1) ;
  func_800->SetNpx (10000) ;
  func_800->Draw ("same") ;
  c_results->cd (5 + i++) ;
  gPad->DrawFrame (200, 1.1 * f_IOS_800->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_800->GetMaximum ()) ;
  f_IOS_800->SetLineWidth (1) ;
  f_IOS_800->SetNpx (10000) ;
  f_IOS_800->Draw ("same") ;
 


  // 1000 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 1000 ;
  rangeScale = 1.5 ;
//  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_1000_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_1000_pars[j] = params_1000[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_1000_pars[j + 4] = params_sig_1000[j] ;
  TF1 * f_IOS_1000 = new TF1 ("f_IOS_1000", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_1000->SetParameters (IOS_1000_pars) ;

  c_results->cd (i) ;
  gPad->DrawFrame (200, 1.1 * func_1000->GetMinimum (), rangeScale * mass, 1.1 * func_sig_1000->GetMaximum ()) ;
  func_sig_1000->SetLineWidth (1) ;
  func_sig_1000->SetNpx (10000) ;
  func_sig_1000->Draw ("same") ;
  func_1000->SetLineWidth (1) ;
  func_1000->SetNpx (10000) ;
  func_1000->Draw ("same") ;
  c_results->cd (5 + i++) ;
  gPad->DrawFrame (200, 1.1 * f_IOS_1000->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_1000->GetMaximum ()) ;
  f_IOS_1000->SetLineWidth (1) ;
  f_IOS_1000->SetNpx (10000) ;
  f_IOS_1000->Draw ("same") ;
 



  TCanvas * test = new TCanvas () ;
  f_IOS_350->Draw () ;
//  func_sig_350->Draw () ;
//  func_350->Draw ("same") ;


}  
  
