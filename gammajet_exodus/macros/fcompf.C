{
TFile *f=new TFile("cocktail_2Gamma_sys_1.root","READ");
TH1D *prob6=((TH1D*)gROOT->FindObject("prob6"))->Clone("prob6");
TH1D *gammaprob6=((TH1D*)gROOT->FindObject("gammaprob6"))->Clone("gammaprob6");
prob6->Divide(gammaprob6);

TFile *g=new TFile("cocktail_2Gamma_sys_1low.root","READ");
TH1D *lowprob6=((TH1D*)gROOT->FindObject("prob6"))->Clone("lowprob6");
TH1D *lowgammaprob6=((TH1D*)gROOT->FindObject("gammaprob6"))->Clone("lowgammaprob6");
lowprob6->Divide(lowgammaprob6);

TFile *h=new TFile("cocktail_2Gamma_sys_1high.root","READ");
TH1D *highprob6=((TH1D*)gROOT->FindObject("prob6"))->Clone("highprob6");
TH1D *highgammaprob6=((TH1D*)gROOT->FindObject("gammaprob6"))->Clone("highgammaprob6");
highprob6->Divide(highgammaprob6);

prob6->Rebin(5);
prob6->Scale(0.2);
lowprob6->Rebin(5);
lowprob6->Scale(0.2);
highprob6->Rebin(5);
highprob6->Scale(0.2);

prob6->Draw();
lowprob6->SetLineColor(kBlue);
lowprob6->Draw("same");
highprob6->SetLineColor(kRed);
highprob6->Draw("same");

}
