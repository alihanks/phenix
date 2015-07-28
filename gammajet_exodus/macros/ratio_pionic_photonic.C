void ratio_pionic_photonic(){
  TFile *f = new TFile("cocktail_10M_pp200_run3_ppfinal.root");

  double pt, value, error;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas("c1","pp @ 200 GeV",800,800);
  c1->SetLeftMargin(0.13);
  c1->SetLogy(0);
  c1->SetTicks();
  c1->SetGridx();
  c1->SetGridy();
  c1->SetFillColor(10);

  TH1D *pionic;
  TH1D *photonic;
  TH1D *ratio_photonic;
  TH1D *ratio_cocktail;

  pionic = (TH1D*)ptePion->Clone();
  pionic->SetTitle("pionic");
  pionic->SetName("pionic");
  pionic->Add(pteCPion);

  photonic = (TH1D*)ptePion->Clone();
  photonic->SetTitle("photonic");
  photonic->SetName("photonic");
  photonic->Add(pteEta);
  photonic->Add(pteOmega);
  photonic->Add(pteEtaprime);
  photonic->Add(ptePhi);
  photonic->Add(pteC);

  ratio_photonic = (TH1D*)pionic->Clone();
  ratio_photonic->Divide(photonic);
  TAxis *yaxis=ratio_photonic->GetYaxis();
  yaxis->SetTitle("pionic / photonic");
  TAxis *xaxis=ratio_photonic->GetXaxis();
  xaxis->SetTitle("p_{T} [GeV/c]");
  ratio_photonic->SetTitleOffset(1.5,"Y");
  ratio_photonic->SetLineWidth(2);
  ratio_photonic->Draw("chist");

  TCanvas *c2 = new TCanvas("c2","pp @ 200 GeV",800,800);
  c2->SetLeftMargin(0.13);
  c2->SetLogy(0);
  c2->SetTicks();
  c2->SetGridx();
  c2->SetGridy();
  c2->SetFillColor(10);

  ratio_cocktail = (TH1D*)pionic->Clone();
  ratio_cocktail->Divide(pte);
  TAxis *yaxis=ratio_cocktail->GetYaxis();
  yaxis->SetTitle("pionic / cocktail");
  TAxis *xaxis=ratio_cocktail->GetXaxis();
  xaxis->SetTitle("p_{T} [GeV/c]");
  ratio_cocktail->SetTitleOffset(1.5,"Y");
  ratio_cocktail->SetLineWidth(2);
  ratio_photonic->SetLineColor(1);
  ratio_cocktail->SetLineColor(4);
  ratio_cocktail->Draw("chist");
  ratio_photonic->Draw("chistsame");

  return;
}
