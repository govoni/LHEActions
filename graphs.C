

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

  TGraph * tg_sAi_par0 = new TGraph (5) ;
  TGraph * tg_sAi_par1 = new TGraph (5) ;
  TGraph * tg_sAi_par2 = new TGraph (5) ;
  TGraph * tg_sAi_par3 = new TGraph (5) ;
  TGraph * tg_sAi_par4 = new TGraph (5) ;
  TGraph * tg_sAi_par5 = new TGraph (5) ;
  TGraph * tg_sAi_par6 = new TGraph (5) ;

  int i = 0 ; 
  
  
  
  
  
  
  // ----> MASS 350 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 350, -2e-06) ;
  tg_par1->SetPoint (i, 350, 346.879) ;
  tg_par2->SetPoint (i, 350, -0.002) ;
  tg_par3->SetPoint (i, 350, 2309.5) ;
  TF1 * func_350 = new TF1 ("func_350",doublePeakModel, 200, 2000, 4) ;
  double params_350[4] = {-2e-06, 346.879, -0.002, 2309.5 } ;
  func_350->SetParameters (params_350) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 350, 0.00794713) ;
  tg_sig_par1->SetPoint (i, 350, 350.346) ;
  tg_sig_par2->SetPoint (i, 350, 7.64082) ;
  tg_sig_par3->SetPoint (i, 350, 1.15547) ;
  tg_sig_par4->SetPoint (i, 350, 1.55098) ;
  tg_sig_par5->SetPoint (i, 350, 1.26826) ;
  tg_sig_par6->SetPoint (i, 350, 2.15497) ;
  TF1 * func_sig_350 = new TF1 ("func_sig_350",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_350[7] = {0.00794713, 350.346, 7.64082, 1.15547, 1.55098, 1.26826, 2.15497 } ;
  func_sig_350->SetParameters (params_sig_350) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 350, 0.00780341) ;
  tg_sAi_par1->SetPoint (i, 350, 349.645) ;
  tg_sAi_par2->SetPoint (i, 350, -7.26896) ;
  tg_sAi_par3->SetPoint (i, 350, 0.999194) ;
  tg_sAi_par4->SetPoint (i, 350, 3) ;
  tg_sAi_par5->SetPoint (i, 350, 1) ;
  tg_sAi_par6->SetPoint (i, 350, 3) ;
  TF1 * func_sAi_350 = new TF1 ("func_sAi_350",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_350[7] = {0.00780341, 349.645, -7.26896, 0.999194, 3, 1, 3 } ;
  func_sAi_350->SetParameters (params_sAi_350) ;
  i++ ;
  
  // ----> MASS 500 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 500, -2e-06) ;
  tg_par1->SetPoint (i, 500, 494.841) ;
  tg_par2->SetPoint (i, 500, -0.002) ;
  tg_par3->SetPoint (i, 500, 5486.16) ;
  TF1 * func_500 = new TF1 ("func_500",doublePeakModel, 200, 2000, 4) ;
  double params_500[4] = {-2e-06, 494.841, -0.002, 5486.16 } ;
  func_500->SetParameters (params_500) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 500, 0.0014435) ;
  tg_sig_par1->SetPoint (i, 500, 503.401) ;
  tg_sig_par2->SetPoint (i, 500, 29.4772) ;
  tg_sig_par3->SetPoint (i, 500, 0.850388) ;
  tg_sig_par4->SetPoint (i, 500, 2.40577) ;
  tg_sig_par5->SetPoint (i, 500, 1.13716) ;
  tg_sig_par6->SetPoint (i, 500, 3.92761) ;
  TF1 * func_sig_500 = new TF1 ("func_sig_500",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_500[7] = {0.0014435, 503.401, 29.4772, 0.850388, 2.40577, 1.13716, 3.92761 } ;
  func_sig_500->SetParameters (params_sig_500) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 500, 0.0013833) ;
  tg_sAi_par1->SetPoint (i, 500, 496.597) ;
  tg_sAi_par2->SetPoint (i, 500, -27.8762) ;
  tg_sAi_par3->SetPoint (i, 500, 1.0491) ;
  tg_sAi_par4->SetPoint (i, 500, 3) ;
  tg_sAi_par5->SetPoint (i, 500, 1) ;
  tg_sAi_par6->SetPoint (i, 500, 3) ;
  TF1 * func_sAi_500 = new TF1 ("func_sAi_500",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_500[7] = {0.0013833, 496.597, -27.8762, 1.0491, 3, 1, 3 } ;
  func_sAi_500->SetParameters (params_sAi_500) ;
  i++ ;
  
  // ----> MASS 650 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 650, -2e-06) ;
  tg_par1->SetPoint (i, 650, 650.136) ;
  tg_par2->SetPoint (i, 650, -0.002) ;
  tg_par3->SetPoint (i, 650, 22244.5) ;
  TF1 * func_650 = new TF1 ("func_650",doublePeakModel, 200, 2000, 4) ;
  double params_650[4] = {-2e-06, 650.136, -0.002, 22244.5 } ;
  func_650->SetParameters (params_650) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 650, 0.000300014) ;
  tg_sig_par1->SetPoint (i, 650, 660.681) ;
  tg_sig_par2->SetPoint (i, 650, 65.4609) ;
  tg_sig_par3->SetPoint (i, 650, 0.839612) ;
  tg_sig_par4->SetPoint (i, 650, 5.83463) ;
  tg_sig_par5->SetPoint (i, 650, 1.18495) ;
  tg_sig_par6->SetPoint (i, 650, 5.0067) ;
  TF1 * func_sig_650 = new TF1 ("func_sig_650",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_650[7] = {0.000300014, 660.681, 65.4609, 0.839612, 5.83463, 1.18495, 5.0067 } ;
  func_sig_650->SetParameters (params_sig_650) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 650, 0.000299991) ;
  tg_sAi_par1->SetPoint (i, 650, 641.439) ;
  tg_sAi_par2->SetPoint (i, 650, -62.4135) ;
  tg_sAi_par3->SetPoint (i, 650, 1.52575) ;
  tg_sAi_par4->SetPoint (i, 650, 3) ;
  tg_sAi_par5->SetPoint (i, 650, 1) ;
  tg_sAi_par6->SetPoint (i, 650, 3) ;
  TF1 * func_sAi_650 = new TF1 ("func_sAi_650",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_650[7] = {0.000299991, 641.439, -62.4135, 1.52575, 3, 1, 3 } ;
  func_sAi_650->SetParameters (params_sAi_650) ;
  i++ ;
  
  // ----> MASS 800 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 800, -2e-06) ;
  tg_par1->SetPoint (i, 800, 737.909) ;
  tg_par2->SetPoint (i, 800, -0.002) ;
  tg_par3->SetPoint (i, 800, 52295.8) ;
  TF1 * func_800 = new TF1 ("func_800",doublePeakModel, 200, 2000, 4) ;
  double params_800[4] = {-2e-06, 737.909, -0.002, 52295.8 } ;
  func_800->SetParameters (params_800) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 800, 8.56526e-05) ;
  tg_sig_par1->SetPoint (i, 800, 804.759) ;
  tg_sig_par2->SetPoint (i, 800, 121.513) ;
  tg_sig_par3->SetPoint (i, 800, 1.44252) ;
  tg_sig_par4->SetPoint (i, 800, 3.15641) ;
  tg_sig_par5->SetPoint (i, 800, 1.33041) ;
  tg_sig_par6->SetPoint (i, 800, 6.86535) ;
  TF1 * func_sig_800 = new TF1 ("func_sig_800",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_800[7] = {8.56526e-05, 804.759, 121.513, 1.44252, 3.15641, 1.33041, 6.86535 } ;
  func_sig_800->SetParameters (params_sig_800) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 800, 7.77657e-05) ;
  tg_sAi_par1->SetPoint (i, 800, 748.164) ;
  tg_sAi_par2->SetPoint (i, 800, -117.077) ;
  tg_sAi_par3->SetPoint (i, 800, 4.16525) ;
  tg_sAi_par4->SetPoint (i, 800, 3) ;
  tg_sAi_par5->SetPoint (i, 800, 1) ;
  tg_sAi_par6->SetPoint (i, 800, 3) ;
  TF1 * func_sAi_800 = new TF1 ("func_sAi_800",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_800[7] = {7.77657e-05, 748.164, -117.077, 4.16525, 3, 1, 3 } ;
  func_sAi_800->SetParameters (params_sAi_800) ;
  i++ ;
  
  // ----> MASS 1000 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 1000, -2e-06) ;
  tg_par1->SetPoint (i, 1000, 864.564) ;
  tg_par2->SetPoint (i, 1000, -0.002) ;
  tg_par3->SetPoint (i, 1000, 55615.7) ;
  TF1 * func_1000 = new TF1 ("func_1000",doublePeakModel, 200, 2000, 4) ;
  double params_1000[4] = {-2e-06, 864.564, -0.002, 55615.7 } ;
  func_1000->SetParameters (params_1000) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 1000, 6.62748e-05) ;
  tg_sig_par1->SetPoint (i, 1000, 938.553) ;
  tg_sig_par2->SetPoint (i, 1000, 217.494) ;
  tg_sig_par3->SetPoint (i, 1000, 1.56219) ;
  tg_sig_par4->SetPoint (i, 1000, 3.37535) ;
  tg_sig_par5->SetPoint (i, 1000, 1.23277) ;
  tg_sig_par6->SetPoint (i, 1000, 4.9202) ;
  TF1 * func_sig_1000 = new TF1 ("func_sig_1000",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_1000[7] = {6.62748e-05, 938.553, 217.494, 1.56219, 3.37535, 1.23277, 4.9202 } ;
  func_sig_1000->SetParameters (params_sig_1000) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 1000, 4.74452e-05) ;
  tg_sAi_par1->SetPoint (i, 1000, 769.473) ;
  tg_sAi_par2->SetPoint (i, 1000, -221.492) ;
  tg_sAi_par3->SetPoint (i, 1000, 3.61003) ;
  tg_sAi_par4->SetPoint (i, 1000, 3) ;
  tg_sAi_par5->SetPoint (i, 1000, 1) ;
  tg_sAi_par6->SetPoint (i, 1000, 3) ;
  TF1 * func_sAi_1000 = new TF1 ("func_sAi_1000",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_1000[7] = {4.74452e-05, 769.473, -221.492, 3.61003, 3, 1, 3 } ;
  func_sAi_1000->SetParameters (params_sAi_1000) ;
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

  tg_sig_par0->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par1->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par2->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par3->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par4->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par5->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par6->GetXaxis ()->SetTitle ("m_{WW}") ;

  tg_sig_par0->GetYaxis ()->SetTitle ("scale") ;
  tg_sig_par1->GetYaxis ()->SetTitle ("mean") ;
  tg_sig_par2->GetYaxis ()->SetTitle ("gaussian sigma") ;
  tg_sig_par3->GetYaxis ()->SetTitle ("right alpha") ;
  tg_sig_par4->GetYaxis ()->SetTitle ("right power law") ;
  tg_sig_par5->GetYaxis ()->SetTitle ("left alpha") ;
  tg_sig_par6->GetYaxis ()->SetTitle ("left power law") ;

  tg_sAi_par0->SetTitle ("") ;
  tg_sAi_par1->SetTitle ("") ;
  tg_sAi_par2->SetTitle ("") ;
  tg_sAi_par3->SetTitle ("") ;
  tg_sAi_par4->SetTitle ("") ;
  tg_sAi_par5->SetTitle ("") ;
  tg_sAi_par6->SetTitle ("") ;

  tg_sAi_par0->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par1->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par2->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par3->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par4->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par5->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par6->GetXaxis ()->SetTitle ("m_{WW}") ;

  tg_sAi_par0->GetYaxis ()->SetTitle ("scale") ;
  tg_sAi_par1->GetYaxis ()->SetTitle ("mean") ;
  tg_sAi_par2->GetYaxis ()->SetTitle ("gaussian sAima") ;
  tg_sAi_par3->GetYaxis ()->SetTitle ("right alpha") ;
  tg_sAi_par4->GetYaxis ()->SetTitle ("right power law") ;
  tg_sAi_par5->GetYaxis ()->SetTitle ("left alpha") ;
  tg_sAi_par6->GetYaxis ()->SetTitle ("left power law") ;

  TCanvas * c_par = new TCanvas () ;
  TF1 * sigmoid = new TF1 ("sigmoid", "[2] + [0] / (1.0 + TMath::Exp(-[1] * (x - [3])))") ;
  sigmoid->SetParameter (0, 25000) ;
  sigmoid->SetParameter (1, 0.005) ;
  sigmoid->SetParameter (2, 30000) ;
  sigmoid->SetParameter (3, 700) ;

  c_par->Divide (2,2) ;
  i = 0 ;
  double x, y ;
  cout << "---> INTERF PARAM 0\n" ;
  tg_par0->GetPoint (1, x, y) ;
  cout << y << endl ;
  c_par->cd (++i) ; tg_par0->Draw ("AL*") ; 
  cout << "---> INTERF PARAM 1\n" ;
  c_par->cd (++i) ; tg_par1->Draw ("AL*") ; tg_par1->Fit ("pol2") ;
  cout << "---> INTERF PARAM 2\n" ;
  tg_par2->GetPoint (1, x, y) ;
  cout << y << endl ;
  c_par->cd (++i) ; tg_par2->Draw ("AL*") ; 
  cout << "---> INTERF PARAM 3\n" ;
  c_par->cd (++i) ; tg_par3->Draw ("AL*") ; tg_par3->Fit (sigmoid) ;
 
  
  TCanvas * c_sig_par = new TCanvas ("c_sig_par", "c_sig_par", 4000, 600) ;
  c_sig_par->Divide (4,2) ;
  i = 0 ;
  cout << "---> SIGNAL PARAM 0\n" ;
  c_sig_par->cd (++i) ; tg_sig_par0->Draw ("AL*") ; tg_sig_par0->Fit ("expo") ;
  cout << "---> SIGNAL PARAM 1\n" ;
  c_sig_par->cd (++i) ; tg_sig_par1->Draw ("AL*") ; tg_sig_par1->Fit ("pol2") ;
  cout << "---> SIGNAL PARAM 2\n" ;
  c_sig_par->cd (++i) ; tg_sig_par2->Draw ("AL*") ; tg_sig_par2->Fit ("pol2") ;
  c_sig_par->cd (++i) ; tg_sig_par3->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par4->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par5->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par6->Draw ("AL*") ;

  c_sig_par->Print ("params_signal.pdf", "pdf") ;
  
  TCanvas * c_sAi_par = new TCanvas ("c_sAi_par", "c_sAi_par", 4000, 600) ;
  c_sAi_par->Divide (4,2) ;
  i = 0 ;
  c_sAi_par->cd (++i) ; tg_sAi_par0->SetLineColor (kRed) ; tg_sAi_par0->Draw ("AL*") ; 
  c_sAi_par->cd (++i) ; tg_sAi_par1->SetLineColor (kRed) ; tg_sAi_par1->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par2->SetLineColor (kRed) ; tg_sAi_par2->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par3->SetLineColor (kRed) ; tg_sAi_par3->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par4->SetLineColor (kRed) ; tg_sAi_par4->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par5->SetLineColor (kRed) ; tg_sAi_par5->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par6->SetLineColor (kRed) ; tg_sAi_par6->Draw ("AL*") ;

  c_sAi_par->Print ("params_phantom.pdf", "pdf") ;

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
  func_sAi_350->SetLineWidth (1) ;
  func_sAi_350->SetLineColor (kRed) ;
  func_sAi_350->SetNpx (10000) ;
  func_sAi_350->Draw ("same") ;
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
  func_sAi_500->SetLineWidth (1) ;
  func_sAi_500->SetLineColor (kRed) ;
  func_sAi_500->SetNpx (10000) ;
  func_sAi_500->Draw ("same") ;
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
  func_sAi_650->SetLineWidth (1) ;
  func_sAi_650->SetLineColor (kRed) ;
  func_sAi_650->SetNpx (10000) ;
  func_sAi_650->Draw ("same") ;
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
  func_sAi_800->SetLineWidth (1) ;
  func_sAi_800->SetLineColor (kRed) ;
  func_sAi_800->SetNpx (10000) ;
  func_sAi_800->Draw ("same") ;
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
  func_sAi_1000->SetLineWidth (1) ;
  func_sAi_1000->SetLineColor (kRed) ;
  func_sAi_1000->SetNpx (10000) ;
  func_sAi_1000->Draw ("same") ;
  func_1000->SetLineWidth (1) ;
  func_1000->SetNpx (10000) ;
  func_1000->Draw ("same") ;
  c_results->cd (5 + i++) ;
  gPad->DrawFrame (200, 1.1 * f_IOS_1000->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_1000->GetMaximum ()) ;
  f_IOS_1000->SetLineWidth (1) ;
  f_IOS_1000->SetNpx (10000) ;
  f_IOS_1000->Draw ("same") ;
 
  c_results->Print ("params_masses.pdf", "pdf") ;

}  
  
