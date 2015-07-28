{
file_exodus_new = new TFile("cocktail_10M_pp200_run3_ppfinal.root");
file_exodus_old = new TFile("cocktail_10M_pp200_run3_bias_final.root");

TCanvas *c1 = new TCanvas("c1","pp @ 200 GeV",800,800);

c1->SetLeftMargin(0.13);
c1->SetLogy(0);
c1->SetTicks();
c1->SetGridx();
c1->SetGridy();
c1->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptFit(1111);
gStyle->SetOptTitle(0);

int i;

// overlay of cocktail

file_exodus_new->cd();
TH1D *ratio;
ratio = (TH1D*)pte->Clone();
ratio->Scale(0.75);
file_exodus_old->cd();
ratio->Divide(pte);

TAxis *yaxis=ratio->GetYaxis();
yaxis->SetTitle("ratio of trigger biased cocktail: new/old");
TAxis *xaxis=ratio->GetXaxis();
xaxis->SetTitle("p_{T} [GeV/c]");
ratio->SetTitleOffset(1.5,"Y");
ratio->SetMinimum(0.0);
ratio->SetMaximum(1.2);
ratio->Draw("");
ratio->Fit("pol0","R","",0.4,5.);

return;

}






