void cocktail_change(){

  TH1D * cold;
  TH1D * cnew;
  TH1D * cratio;

  TFile *f = new TFile("cocktail_10M_auau200_minbias_final.root");
  TFile *g = new TFile("cocktail_10M_auau200_minbias_final_update.root");

  f->cd();
  cold = (TH1D*)pte->Clone();
  g->cd();
  cnew = (TH1D*)pte->Clone();

  cratio = (TH1D*)cnew->Clone();
  cratio->Divide(cold);

  TCanvas *c1 = new TCanvas("c1","AuAu @ 200 GeV minimum bias",800,800);

  c1->SetLeftMargin(0.13);
  c1->SetLogy(0);
  c1->SetTicks();
  c1->SetGridx();
  c1->SetGridy();
  c1->SetFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  TAxis *yaxis=cratio->GetYaxis();
  yaxis->SetTitle("cocktail ratio: update/original ");
  TAxis *xaxis=cratio->GetXaxis();
  xaxis->SetTitle("p_{T} [GeV/c]");
  cratio->SetTitleOffset(1.5,"Y");
  cratio->SetMaximum(1.2);

  cratio->Fit("pol0");

  return;

}
