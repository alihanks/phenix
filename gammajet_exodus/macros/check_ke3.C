void check_ke3(){
  TFile *f = new TFile("cocktail_1M_pp200_run3_final.root");

  double pt, value, error;

  TCanvas *c1 = new TCanvas("c1","pp @ 200 GeV",800,800);
  c1->SetLeftMargin(0.13);
  c1->SetLogy();
  c1->SetTicks();
  c1->SetGridx();
  c1->SetGridy();
  c1->SetFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);

  TF1 *ke3 = new TF1("ke3","[0]*exp(-x/0.453)/pow(x+0.408,7.249)",0.,6.);

  TAxis *yaxis=pteKe3->GetYaxis();
  yaxis->SetTitle("Ed^{3}#sigma/d^{3}p [mb GeV^{-2}/c^{3}]");
  TAxis *xaxis=pteKe3->GetXaxis();
  xaxis->SetTitle("p_{T} [GeV/c]");
  pteKe3->SetTitleOffset(1.5,"Y");
  pteKe3->Draw();
  pteKe3->Fit("ke3","R","",0.4,6.0);

  return;
}
