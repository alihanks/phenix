{
file_biased = new TFile("cocktail_10M_pp200_run3_bias_final.root");
file_unbiased = new TFile("cocktail_10M_pp200_run3_final.root");

TCanvas *c1 = new TCanvas();

c1->SetLogy();
c1->SetTicks();
c1->SetGridx();
c1->SetGridy();
c1->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptFit(1111);
gStyle->SetOptTitle(0);

file_biased->cd();
TAxis *yaxis=pte->GetYaxis();
yaxis->SetTitle("(1/2#pi p_{T})dN/dp_{T}dy [(c/GeV)^{2}]");
TAxis *xaxis=pte->GetXaxis();
xaxis->SetTitle("p_{T} [GeV/c]");
pte->SetTitleOffset(1.5,"Y");

TH1D *biased;
biased = (TH1D*)pte->Clone();
biased->SetName("biased");
biased->SetLineColor(1);
biased->SetLineWidth(2);
biased->SetMinimum(1.0e-09);
biased->SetMaximum(5.0e+02);
biased->Draw("chist");

file_unbiased->cd();
TH1D *unbiased;
unbiased = (TH1D*)pte->Clone();
unbiased->SetName("unbiased");
unbiased->SetLineColor(2);
unbiased->SetLineWidth(2);
unbiased->Draw("chistsame");

TCanvas *c2 = new TCanvas();

c2->SetLeftMargin(0.13);
c2->SetLogy(0);
c2->SetTicks();
c2->SetGridx();
c2->SetGridy();
c2->SetFillColor(10);

TH1D *ratio = (TH1D*)biased->Clone();
ratio->SetName("ratio");
ratio->Divide(unbiased);
ratio->SetMinimum(0.6);
ratio->SetMaximum(0.8);
ratio->SetLineColor(1);
TAxis *yaxis=ratio->GetYaxis();
yaxis->SetTitle("trigger biased / unbiased cocktail");
TAxis *xaxis=ratio->GetXaxis();
xaxis->SetTitle("p_{T} [GeV/c]");
ratio->SetTitleOffset(1.5,"Y");
ratio->Draw();
TF1 *low  = new TF1("low","0.611+0.191*x-0.0613*x*x",0.,1.158);
TF1 *high = new TF1("high","0.75",1.158,6.);
low->Draw("same");
high->Draw("same");

}






