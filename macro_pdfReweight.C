double mg_scale =  0.0417 / 100000 ;
double ph_scale = 0.0993 / 496500 ;

/*
void addOverflowBin (TH1F * histo)
{
  TH1F * outcome = new TH1F ("tempo", histo->GetTitle (), histo->GenNbinsX () + 1, 
  


}
*/



// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



void drawTwoHistos (TFile & f_ph, TFile & f_mg, TString name, int reverse = 0)
{
  TH1F * h_ph = (TH1F *) f_ph.Get (name) ;
  h_ph->SetName (name + "_ph") ;
  h_ph->GetXaxis ()->SetTitle (h_ph->GetTitle ()) ;
  h_ph->SetTitle ("") ;
  h_ph->SetStats (0) ;

  TH1F * h_mg = (TH1F *) f_mg.Get (name) ;
  h_mg->SetName (name + "_mg") ;
  h_mg->GetXaxis ()->SetTitle (h_mg->GetTitle ()) ;
  h_mg->SetTitle ("") ;
  h_mg->SetStats (0) ;
  h_mg->SetLineColor (kRed) ;
  
  h_mg->Scale (mg_scale) ;
  h_ph->Scale (ph_scale) ;
  leg = new TLegend (0.9, 0.8, 1, 0.9) ;
  leg->AddEntry (h_ph, "ph", "l") ;
  leg->AddEntry (h_mg, "mg", "lf") ;

  h_mg->SetFillColor (kGray) ;
  TCanvas c1 ;
  if (!reverse)
    {
      h_mg->Draw ("E2") ;
      h_ph->Draw ("Esame") ;
    }
  else
    {
      h_ph->Draw ("E") ;
      h_mg->Draw ("E2same") ;
      h_ph->Draw ("Esame") ;
    }

  leg->Draw () ;
  c1.Print (name + ".pdf", "pdf") ;
  return ;
}





int macro_pdfReweight ()
{

  TFile f_ph ("reweight_ph.root") ;
  TFile f_mg ("reweight_mg.root") ;
//  TFile f_ph ("reweight_unweight_ph.root") ;
//  TFile f_mg ("reweight_unweight_mg.root") ;
  drawTwoHistos (f_ph, f_mg, "h_scale", 1) ;
  drawTwoHistos (f_ph, f_mg, "h_x", 1) ;
  drawTwoHistos (f_ph, f_mg, "h_mjj", 0) ;
  drawTwoHistos (f_ph, f_mg, "h_mll", 1) ;
  drawTwoHistos (f_ph, f_mg, "h_ptj1", 1) ;
  drawTwoHistos (f_ph, f_mg, "h_ptj2", 0) ;
  drawTwoHistos (f_ph, f_mg, "h_etaj1", 1) ;
  drawTwoHistos (f_ph, f_mg, "h_etaj2", 0) ;
  drawTwoHistos (f_ph, f_mg, "h_weight", 1) ;
  

  TH1F * h_mll_ph = (TH1F *) f_ph.Get ("h_mll") ;
  h_mll_ph->SetName ("h_mll_ph") ;
  TH1F * h_mll_mg = (TH1F *) f_mg.Get ("h_mll") ;
  h_mll_mg->SetName ("h_mll_mg") ;
  h_mll_mg->Scale (mg_scale) ;
  h_mll_ph->Scale (ph_scale) ;

  cout << "        \t nocut    \t mll_cut\n" ;
  cout << "phantom \t" ;
  cout << h_mll_ph->Integral () << "\t" ;
  cout << h_mll_ph->Integral (h_mll_ph->FindBin (12), h_mll_ph->FindBin (50)) << endl ;
  cout << "madgraph\t" ;
  cout << h_mll_mg->Integral () << "\t" ;
  cout << h_mll_mg->Integral (h_mll_mg->FindBin (12), h_mll_mg->FindBin (50)) << endl ;
  cout << "ratio   \t" ;
  cout << h_mll_ph->Integral () / h_mll_mg->Integral () << "\t" ;
  cout << h_mll_ph->Integral (h_mll_ph->FindBin (12), h_mll_ph->FindBin (50)) / h_mll_mg->Integral (h_mll_mg->FindBin (12), h_mll_mg->FindBin (50)) << endl ;
  return 0 ;
}
