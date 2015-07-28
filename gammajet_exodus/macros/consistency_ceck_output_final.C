{
TCanvas *c1 = new TCanvas("c1","AuAu @ 200 GeV minimum bias",800,800);

c1->SetLeftMargin(0.13);
c1->SetLogy();
c1->SetTicks();
c1->SetGridx();
c1->SetGridy();
c1->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TH1D *minbias = new TH1D;
TH1D *cen0 = new TH1D;

TFile *f_in = new TFile("cocktail_10M_auau200_minbias_final.root");
minbias = (TH1D*)pte->Clone();
TFile *f_in0 = new TFile("cocktail_10M_auau200_0_10_final.root");
cen0 = (TH1D*)pte->Clone();
TFile *f_in1 = new TFile("cocktail_10M_auau200_10_20_final.root");
cen1 = (TH1D*)pte->Clone();
TFile *f_in2 = new TFile("cocktail_10M_auau200_20_40_final.root");
cen2 = (TH1D*)pte->Clone();
TFile *f_in3 = new TFile("cocktail_10M_auau200_40_60_final.root");
cen3 = (TH1D*)pte->Clone();
TFile *f_in4 = new TFile("cocktail_10M_auau200_60_92_final.root");
cen4 = (TH1D*)pte->Clone();

TH1F *h1=new TH1F("h1","",320,0.,6.);
h1->SetMaximum(400.);
h1->SetMinimum(1.e-9);
h1->SetYTitle("Ed^{3}N/dp^{3} [GeV^{-2}c^{3}]");
h1->SetXTitle("p_{T} [GeV/c]");
h1->SetTitleOffset(1.5,"Y");
h1->Draw();

minbias->SetLineColor(1);
cen0->SetLineColor(2);
cen1->SetLineColor(3);
cen2->SetLineColor(4);
cen3->SetLineColor(5);
cen4->SetLineColor(6);

minbias->SetLineWidth(2);
cen0->SetLineWidth(2);
cen1->SetLineWidth(2);
cen2->SetLineWidth(2);
cen3->SetLineWidth(2);
cen4->SetLineWidth(2);

minbias->Draw("chistsame");
cen0->Draw("chistsame");
cen1->Draw("chistsame");
cen2->Draw("chistsame");
cen3->Draw("chistsame");
cen4->Draw("chistsame");

TH1D *total = cen0->Clone();
total->Add(cen1);
total->Add(cen2,2.);
total->Add(cen3,2.);
total->Add(cen4,3.2);
total->Scale(1./9.2);
total->SetLineColor(7);
total->SetLineWidth(2);
total->Draw("chistsame");

TCanvas *c2 = new TCanvas("c2","AuAu @ 200 GeV minimum bias",800,800);

c2->SetLeftMargin(0.13);
c2->SetLogy(0);
c2->SetTicks();
c2->SetGridx();
c2->SetGridy();
c2->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TH1F *h2=new TH1F("h2","",320,0.,6.);
h2->SetMinimum(0.);
h2->SetMaximum(1.2);
h2->SetYTitle("centrality average / min. bias");
h2->SetXTitle("p_{T} [GeV/c]");
h2->SetTitleOffset(1.5,"Y");
h2->Draw();

TH1D *ratio = total->Clone();
ratio->Divide(minbias);
ratio->SetLineColor(1);
ratio->SetLineWidth(2);
ratio->Draw("same");

TCanvas *c3 = new TCanvas("c3","AuAu @ 200 GeV minimum bias",800,800);

c3->SetLeftMargin(0.13);
c3->SetLogy();
c3->SetTicks();
c3->SetGridx();
c3->SetGridy();
c3->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TH1D *minbias_npart = minbias->Clone();
minbias_npart->Scale(1./109.1);
TH1D *cen0_npart = cen0->Clone();
cen0_npart->Scale(1./325.2);
TH1D *cen1_npart = cen1->Clone();
cen1_npart->Scale(1./234.6);
TH1D *cen2_npart = cen2->Clone();
cen2_npart->Scale(1./140.4);
TH1D *cen3_npart = cen3->Clone();
cen3_npart->Scale(1./60.0);
TH1D *cen4_npart = cen4->Clone();
cen4_npart->Scale(1./14.5);

TH1F *h3=new TH1F("h3","",320,0.,6.);
h3->SetMaximum(10.);
h3->SetMinimum(1.e-10);
h3->SetYTitle("Ed^{3}N/dp^{3}/<N_{part}> [GeV^{-2}c^{3}]");
h3->SetXTitle("p_{T} [GeV/c]");
h3->SetTitleOffset(1.5,"Y");
h3->Draw();

minbias_npart->SetLineColor(1);
cen0_npart->SetLineColor(2);
cen1_npart->SetLineColor(3);
cen2_npart->SetLineColor(4);
cen3_npart->SetLineColor(5);
cen4_npart->SetLineColor(6);

minbias_npart->Draw("chistsame");
cen0_npart->Draw("chistsame");
cen1_npart->Draw("chistsame");
cen2_npart->Draw("chistsame");
cen3_npart->Draw("chistsame");
cen4_npart->Draw("chistsame");

TCanvas *c4 = new TCanvas("c4","AuAu @ 200 GeV minimum bias",800,800);

c4->SetLeftMargin(0.13);
c4->SetLogy();
c4->SetTicks();
c4->SetGridx();
c4->SetGridy();
c4->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TH1D *minbias_ncoll = minbias->Clone();
minbias_ncoll->Scale(1./257.8);
TH1D *cen0_ncoll = cen0->Clone();
cen0_ncoll->Scale(1./955.4);
TH1D *cen1_ncoll = cen1->Clone();
cen1_ncoll->Scale(1./602.6);
TH1D *cen2_ncoll = cen2->Clone();
cen2_ncoll->Scale(1./296.8);
TH1D *cen3_ncoll = cen3->Clone();
cen3_ncoll->Scale(1./90.7);
TH1D *cen4_ncoll = cen4->Clone();
cen4_ncoll->Scale(1./14.5);

TH1F *h3=new TH1F("h3","",320,0.,6.);
h3->SetMaximum(10.);
h3->SetMinimum(1.e-10);
h3->SetYTitle("Ed^{3}N/dp^{3}/<N_{coll}> [GeV^{-2}c^{3}]");
h3->SetXTitle("p_{T} [GeV/c]");
h3->SetTitleOffset(1.5,"Y");
h3->Draw();

minbias_ncoll->SetLineColor(1);
cen0_ncoll->SetLineColor(2);
cen1_ncoll->SetLineColor(3);
cen2_ncoll->SetLineColor(4);
cen3_ncoll->SetLineColor(5);
cen4_ncoll->SetLineColor(6);

minbias_ncoll->Draw("chistsame");
cen0_ncoll->Draw("chistsame");
cen1_ncoll->Draw("chistsame");
cen2_ncoll->Draw("chistsame");
cen3_ncoll->Draw("chistsame");
cen4_ncoll->Draw("chistsame");


}






