void calculate_e_pizero(int icent=0){
  if (icent==0){
    TFile *file_exodus = new TFile("cocktail_10M_auau200_minbias_final_update.root");
    TF1 *pizero = new TF1("pizero","495.7/pow(exp(-0.5235*x-0.1475*x*x)+x/0.7392,8.253)",0.,6.);
  }
  if (icent==1){
    TFile *file_exodus = new TFile("cocktail_10M_auau200_0_10_final_update.root");
    TF1 *pizero = new TF1("pizero","1744.0/pow(exp(-0.4802*x-0.3389*x*x)+x/0.6827,8.129)",0.,6.);
  }
  if (icent==2){
    TFile *file_exodus = new TFile("cocktail_10M_auau200_10_20_final_update.root");
    TF1 *pizero = new TF1("pizero","1045.0/pow(exp(-0.5210*x-0.1885*x*x)+x/0.7328,8.247)",0.,6.);
  }
  if (icent==3){
    TFile *file_exodus = new TFile("cocktail_10M_auau200_20_40_final_update.root");
    TF1 *pizero = new TF1("pizero","600.1/pow(exp(-0.4669*x-0.08806*x*x)+x/0.7938,8.493)",0.,6.);
  }
  if (icent==4){
    TFile *file_exodus = new TFile("cocktail_10M_auau200_40_60_final_update.root");
    TF1 *pizero = new TF1("pizero","279.9/pow(exp(-0.4054*x-0.06311*x*x)+x/0.7849,8.443)",0.,6.);
  }
  if (icent==5){
    TFile *file_exodus = new TFile("cocktail_10M_auau200_60_80_final_update.root");
    TF1 *pizero = new TF1("pizero","93.57/pow(exp(-0.3500*x-0.02792*x*x)+x/0.7753,8.431)",0.,6.);
  }
  if (icent==6){
    TFile *file_exodus = new TFile("cocktail_10M_auau200_60_92_final_update.root");
    TF1 *pizero = new TF1("pizero","71.86/pow(exp(-0.06022*x-0.01968*x*x)+x/1.306,10.91)",0.,6.);
  }

  TCanvas *c1 = new TCanvas("c1","AuAu @ 200 GeV minimum bias",800,800);

  c1->SetLeftMargin(0.13);
  c1->SetLogy();
  c1->SetTicks();
  c1->SetGridx();
  c1->SetGridy();
  c1->SetFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);

  file_exodus->cd();

  TAxis *yaxis=pte->GetYaxis();
  yaxis->SetTitle("(e^{+} + e^{-})/(2#pi^{0})");
  TAxis *xaxis=pte->GetXaxis();
  xaxis->SetTitle("p_{T} [GeV/c]");

  pte->SetTitleOffset(1.5,"Y");

  for (int ibin=1; ibin<=pte->GetNbinsX(); ibin++){
    double pt_val = pte->GetBinCenter(ibin);
    double ratio  = pte->GetBinContent(ibin)/pizero->Eval(pt_val);
    pte->SetBinContent(ibin,ratio);
  }
  pte->SetMinimum(1.0e-05);
  pte->SetMaximum(1.0e-01);
  pte->SetLineWidth(2);
  pte->Draw("chist");

  for (int ibin=1; ibin<=pte->GetNbinsX(); ibin++){
    double pt_val = ptePion->GetBinCenter(ibin);
    double ratio  = ptePion->GetBinContent(ibin)/pizero->Eval(pt_val);
    //double ratio  = (ptePion->GetBinContent(ibin)+pteCPion->GetBinContent(ibin))/pizero->Eval(pt_val);
    ptePion->SetBinContent(ibin,ratio);
  }
  ptePion->SetLineWidth(2);
  ptePion->SetLineColor(2);
  ptePion->Draw("chistsame");

  for (int ibin=1; ibin<=pte->GetNbinsX(); ibin++){
    double pt_val = pteEta->GetBinCenter(ibin);
    double ratio  = pteEta->GetBinContent(ibin)/pizero->Eval(pt_val);
    //double ratio  = (pteEta->GetBinContent(ibin)+pteCEta->GetBinContent(ibin))/pizero->Eval(pt_val);
    pteEta->SetBinContent(ibin,ratio);
  }
  pteEta->SetLineWidth(2);
  pteEta->SetLineColor(4);
  pteEta->Draw("chistsame");

  for (int ibin=1; ibin<=pte->GetNbinsX(); ibin++){
    double pt_val = pteC->GetBinCenter(ibin);
    double ratio  = pteC->GetBinContent(ibin)/pizero->Eval(pt_val);
    //double ratio  = (pteC->GetBinContent(ibin)+pteCEta->GetBinContent(ibin))/pizero->Eval(pt_val);
    pteC->SetBinContent(ibin,ratio);
  }
  pteC->SetLineWidth(2);
  pteC->SetLineColor(3);
  pteC->Draw("chistsame");
  
  TLatex *minbias= new TLatex(0.,log10(0.12),"Au+Au @ #sqrt{s} = 200 GeV ; minimum bias");
  minbias->SetTextColor(4);
  minbias->Draw();

  TLatex *tsum = new TLatex(3.2,-1.7,"all electrons");
  tsum->SetTextColor(1);
  tsum->SetTextSize(0.04);
  tsum->Draw();
  TLatex *tconv = new TLatex(3.2,-1.9,"#gamma conversion");
  tconv->SetTextColor(3);
  tconv->SetTextSize(0.04);
  tconv->Draw();
  TLatex *tpi = new TLatex(3.2,-2.1,"#pi^{0} #rightarrow #gammaee");
  tpi->SetTextColor(2);
  tpi->SetTextSize(0.04);
  tpi->Draw();
  TLatex *teta = new TLatex(3.2,-2.3,"#eta #rightarrow #gammaee");
  teta->SetTextColor(4);
  teta->SetTextSize(0.04);
  teta->Draw();

  return;
}






