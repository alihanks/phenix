{
file_sys = new TFile("cocktail_systematic_pp.root");
file_sys->cd();

TCanvas *c1 = new TCanvas("c1","electron cocktail systematics: pp @ 200 GeV",
800,800);

c1->SetLeftMargin(0.13);
c1->SetTicks();
c1->SetGridx();
c1->SetGridy();
c1->SetFillColor(10);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

rel_error_max->SetMinimum(-30.0);
rel_error_max->SetMaximum(+30.0);
rel_error_max->SetLineWidth(3.0);
rel_error_min->SetLineWidth(3.0);

TAxis *yaxis=rel_error_max->GetYaxis();
yaxis->SetTitle("total systematic uncertainty [%]");
TAxis *xaxis=rel_error_max->GetXaxis();
xaxis->SetTitle("p_{T} [GeV/c]");
rel_error_max->SetTitleOffset(1.5,"Y");

rel_error_max->Draw("chist");
rel_error_min->Scale(-1.);
rel_error_min->Draw("chistsame");

}






