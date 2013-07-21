

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


int macro_findInterferece (TString filename)                                                        
{        
  TFile * f = new TFile (filename) ;
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

  TH1F * ratio = (TH1F *) h_MWW_mg->Clone ("ratio") ;
  ratio->GetXaxis ()->SetTitle (ratio->GetName ()) ;
  ratio->SetTitle ("") ;
  ratio->Divide (diff) ;
  //cout << "scaling by " << 1. / ratio->GetBinContent (ratio->FindBin (500)) << endl ;
  //h_MWW_mg->Scale (1. / ratio->GetBinContent (ratio->FindBin (500))) ;

  TH1F * delta = (TH1F *) h_MWW_mg->Clone ("delta") ;
  delta->SetTitle ("") ;
  delta->Add (diff, -1) ;
  delta->SetLineColor (kGreen + 2) ;

  TH1F * relDiff = delta->Clone ("relDiff") ;
  relDiff->Divide (h_MWW_mg) ;

  //PG initial spectra

  TCanvas * c1 = new TCanvas () ;
  c1.DrawFrame (100, 0.00001, 2000, 0.02) ;
  h_MWW_phbkg->Draw ("histsame") ;
  h_MWW_phbkgsig->Draw ("histsame") ;
  h_MWW_mg->Draw ("histsame") ;
  c1->Print ("spectra.pdf", "pdf") ;

  //PG SBI - B

  TCanvas * c2 = new TCanvas () ;
  diff->SetTitle ("") ;
  diff->Draw ("hist") ;
  c2->Print ("diff.pdf", "pdf") ;

  //PG S / (SBI - S)

  TCanvas * c3 = new TCanvas () ;
  ratio->SetTitle ("") ;
  ratio->Draw ("hist") ;
  ratio->Draw ("Esame") ;
  c3->Print ("ratio.pdf", "pdf") ;

  //PG S only, and (SBI - B)

  TCanvas * c4 = new TCanvas () ;
  h_MWW_mg->SetTitle ("") ;
  h_MWW_mg->Draw ("hist") ;
  diff->Draw ("histsame") ;
  c4->Print ("signals.pdf", "pdf") ;

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
//  func2->FixParameter (0, 500.)
//  func2->FixParameter (1, 0.)
//
//  func2->SetParameter (2, 0.0005) ;
//  func2->SetParameter (3, 0.05) ;
//  delta->Fit ("func2", "", "", 200, 2000) ;


  TF1 * func3 = new TF1 ("func3", doubleSlope, 0, 2000, 5) ;
  func3->SetLineWidth (1) ;
  func3->SetLineColor (kGreen + 2) ;
  func3->FixParameter (0, 500.) ;
  func3->FixParameter (1, 0.) ;

  func3->SetParameter (2, 0.0005) ;
  func3->SetParameter (3, 0.05) ;
  func3->SetParameter (4, 0.05) ;

  func3->SetParName (0, "centre") ;
  func3->SetParName (1, "shift") ;
  func3->SetParName (2, "scale") ;
  func3->SetParName (3, "slope_right") ;
  func3->SetParName (4, "slope_left") ;

  delta->Fit ("func3", "+", "", 200, 2000) ;

  TF1 * func4 = new TF1 ("func4", parabolicAndExp, 0, 2000, 5) ;
  func4->SetLineWidth (1) ;
  func4->SetLineColor (kRed + 2) ;
  func4->FixParameter (0, 500.) ;
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

  TF1 * func5 = new TF1 ("func5", sinusAndExp, 0, 2000, 5) ;
  func5->SetLineWidth (1) ;
  func5->SetLineColor (kRed + 2) ;
  func5->FixParameter (0, 500.) ;
  func5->FixParameter (1, 0.) ;

  func5->FixParameter (2, func3->GetParameter (2)) ;
  func5->FixParameter (3, func3->GetParameter (3)) ;
  func5->SetParameter (4, 10) ;

  func5->SetParName (0, "centre") ;
  func5->SetParName (1, "shift") ;
  func5->SetParName (2, "scale") ;
  func5->SetParName (3, "slope_right") ;
  func5->SetParName (4, "freq_left") ;

  delta->Fit ("func5", "+", "", 200, 2000) ;

  TF1 * func6 = new TF1 ("func6", withLeftAddOn, 0, 2000, 6) ;
  func6->SetLineWidth (1) ;
  func6->SetLineColor (kBlue + 2) ;
  func6->FixParameter (0, 500.) ;
  func6->FixParameter (1, 0.) ;

  func6->FixParameter (2, func3->GetParameter (2)) ;
  func6->FixParameter (3, func3->GetParameter (3)) ;

  func6->SetParName (0, "centre") ;
  func6->SetParName (1, "shift") ;
  func6->SetParName (2, "scale") ;
  func6->SetParName (3, "slope") ;
  func6->SetParName (4, "secondCentre") ;
  func6->SetParName (5, "secondSlope") ;

  delta->Fit ("func6", "+", "", 200, 2000) ;

  TCanvas * c6 = new TCanvas () ;
  relDiff->Draw ("hist") ;

  c6->Print ("relDiff.pdf", "pdf") ;


}                                                                                   
                                                                                    
