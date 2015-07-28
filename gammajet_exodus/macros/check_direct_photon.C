void check_direct_photon(){
  TFile *f = new TFile("cocktail_1M_pp200_run3_final_test_direct_photon.root");

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

  TF1 *vogelsang = new TF1("vogelsang","1.385/pow(exp(0.1524*x)+x/0.2537,5.821)",0.,15.);
  vogelsang->SetLineWidth(2);
  vogelsang->SetLineColor(2);

  pteCGamma->Scale(1./42.2);
  TAxis *yaxis=pteCGamma->GetYaxis();
  yaxis->SetTitle("Ed^{3}N/d^{3}p [GeV^{-2}/c^{3}]");
  TAxis *xaxis=pteCGamma->GetXaxis();
  xaxis->SetTitle("p_{T} [GeV/c]");
  pteCGamma->SetTitleOffset(1.5,"Y");
  pteCGamma->SetLineWidth(2);
  pteCGamma->Draw("hist");
  vogelsang->Draw("same");

  TH1D *ratio;
  ratio = (TH1D*)pteCGamma->Clone();
  for (int ibin=1; ibin<=150; ibin++){
    pt = 0.05+0.1*(ibin-1);
    value = pteCGamma->GetBinContent(ibin)/vogelsang->Eval(pt);
    error = pteCGamma->GetBinError(ibin)/pteCGamma->GetBinContent(ibin);
    error = value*error;
    ratio->SetBinContent(ibin,value);
    ratio->SetBinError(ibin,error);
  }

  TCanvas *c2 = new TCanvas("c2","pp @ 200 GeV",800,800);
  c2->SetLeftMargin(0.13);
  c2->SetLogy(0);
  c2->SetTicks();
  c2->SetGridx();
  c2->SetGridy();
  c2->SetFillColor(10);
  TAxis *yaxis2=ratio->GetYaxis();
  yaxis2->SetTitle("data/fit");
  ratio->Fit("pol0","R","",.5,13.);

  return;
}
