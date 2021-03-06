{
file_exodus = new TFile("cocktail_10M_auau200_minbias_final.root");

TCanvas *c1 = new TCanvas("c1","AuAu @ 200 GeV minimum bias",800,800);

c1->SetLeftMargin(0.13);
c1->SetLogy();
c1->SetTicks();
c1->SetGridx();
c1->SetGridy();
c1->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

file_exodus->cd();

TAxis *yaxis=pte->GetYaxis();
yaxis->SetTitle("(1/2#pi p_{T})dN/dp_{T}dy [(c/GeV)^{2}]");
TAxis *xaxis=pte->GetXaxis();
xaxis->SetTitle("p_{T} [GeV/c]");
pte->SetTitleOffset(1.5,"Y");

TH1D *total;
total = (TH1D*)pte->Clone();
total->SetName("total");
total->SetLineColor(1);
total->SetLineWidth(2);
total->SetMinimum(1.0e-08);
total->SetMaximum(5.0e+02);
total->Draw("chist");

file_1 = new TFile("cocktail_10M_auau200_0_10_final.root");
TH1D *total_1 = (TH1D*)pte->Clone();
total_1->SetName("total_1");
total_1->Scale(10.);
total_1->SetLineColor(2);
total_1->SetLineWidth(2);
total_1->Draw("chistsame");

file_2 = new TFile("cocktail_10M_auau200_10_20_final.root");
TH1D *total_2 = (TH1D*)pte->Clone();
total_2->SetName("total_2");
total_2->Scale(10.);
total_2->SetLineColor(3);
total_2->SetLineWidth(2);
total_2->Draw("chistsame");

file_3 = new TFile("cocktail_10M_auau200_20_40_final.root");
TH1D *total_3 = (TH1D*)pte->Clone();
total_3->SetName("total_3");
total_3->Scale(20.);
total_3->SetLineColor(4);
total_3->SetLineWidth(2);
total_3->Draw("chistsame");

file_4 = new TFile("cocktail_10M_auau200_40_60_final.root");
TH1D *total_4 = (TH1D*)pte->Clone();
total_4->SetName("total_4");
total_4->Scale(20.);
total_4->SetLineColor(5);
total_4->SetLineWidth(2);
total_4->Draw("chistsame");

file_5 = new TFile("cocktail_10M_auau200_60_92_final.root");
TH1D *total_5 = (TH1D*)pte->Clone();
total_5->SetName("total_5");
total_5->Scale(32.);
total_5->SetLineColor(6);
total_5->SetLineWidth(2);
total_5->Draw("chistsame");

TH1D *total_cent = (TH1D*)total_1->Clone();
total_cent->SetName("total_cent");
total_cent->Add(total_2);
total_cent->Add(total_3);
total_cent->Add(total_4);
total_cent->Add(total_5);
total_cent->Scale(1./92.);
total_cent->SetLineColor(7);
total_cent->SetLineWidth(2);
total_cent->Draw("chistsame");

TCanvas *c2 = new TCanvas("c2","AuAu @ 200 GeV minimum bias",800,800);

c2->SetLeftMargin(0.13);
c2->SetLogy(0);
c2->SetTicks();
c2->SetGridx();
c2->SetGridy();
c2->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TH1D *ratio = (TH1D*)total_cent->Clone();
ratio->SetName("ratio");
ratio->Divide(total);
ratio->SetMinimum(0.);
ratio->SetMaximum(1.2);
ratio->SetLineColor(1);
TAxis *yaxis=ratio->GetYaxis();
yaxis->SetTitle("centrality average / min. bias");
TAxis *xaxis=ratio->GetXaxis();
xaxis->SetTitle("p_{T} [GeV/c]");
ratio->SetTitleOffset(1.5,"Y");
ratio->Draw();


}






