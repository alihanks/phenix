void systematics_auau200_minbias_hard(){
  gROOT->Reset();
    
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(1.1);
  gStyle->SetFillColor(0);

  c1 = new TCanvas();
  c1->SetGridx(1);
  c1->SetGridy(1);
  c1->SetTicks(1,1);
  c1->Draw();
  
  TH2F *hist_background = new TH2F("hist_background","Non-photonic electron invariant cross-section",10,0.,4.5,10,0.,1.2);
  hist_background->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_background->GetYaxis()->SetTitle("rel. sys. error");
  hist_background->GetYaxis()->SetTitleOffset(1.2); 
  hist_background->Draw();

  TFile *fh = new TFile("cocktail_10M_auau200_minbias_final_hard.root");
  TH1D *hard = (TH1D*) fh->Get("pteR");
  hard->SetLineColor(2);
  //  hard->Draw("chist");

  TFile *fs = new TFile("cocktail_10M_auau200_minbias_final_soft.root");
  TH1D *soft = (TH1D*) fs->Get("pteR");
  soft->SetLineColor(2);
  //  soft->Draw("chist");

  TFile *fc = new TFile("cocktail_10M_auau200_minbias_final.root");
  TH1D *best = (TH1D*) fc->Get("pteR");
  //  best->Draw("chistsame");

  TH1F *sys_pishape = (TH1F*)hard->Clone();
  // TH1F *sys_pishape = (TH1F*)soft->Clone();
  sys_pishape->SetLineColor(1);
  sys_pishape->SetLineWidth(2);
  sys_pishape->Divide(best);
  sys_pishape->Draw("chistsame");

  TH1F *sys_etapi = (TH1F*)best->Clone();
  sys_etapi->SetLineColor(2);
  sys_etapi->SetLineWidth(2);
  sys_etapi->Add(pteREta,-0.22);
  sys_etapi->Add(pteRCEta,-0.22);
  sys_etapi->Divide(best);
  sys_etapi->Draw("chistsame");

  TH1F *sys_rhopi = (TH1F*)best->Clone();
  sys_rhopi->SetLineColor(3);
  sys_rhopi->SetLineWidth(2);
  sys_rhopi->Add(pteRRho,-0.30);
  sys_rhopi->Divide(best);
  sys_rhopi->Draw("chistsame");

  TH1F *sys_omegapi = (TH1F*)best->Clone();
  sys_omegapi->SetLineColor(4);
  sys_omegapi->SetLineWidth(2);
  sys_omegapi->Add(pteROmega,-0.30);
  sys_omegapi->Divide(best);
  sys_omegapi->Draw("chistsame");

  TH1F *sys_etaprimepi = (TH1F*)best->Clone();
  sys_etaprimepi->SetLineColor(5);
  sys_etaprimepi->SetLineWidth(2);
  sys_etaprimepi->Add(pteREtaprime,-0.30);
  sys_etaprimepi->Add(pteRCEtaprime,-0.30);
  sys_etaprimepi->Divide(best);
  sys_etaprimepi->Draw("chistsame");

  TH1F *sys_phipi = (TH1F*)best->Clone();
  sys_phipi->SetLineColor(6);
  sys_phipi->SetLineWidth(2);
  sys_phipi->Add(pteRPhi,-0.30);
  sys_phipi->Divide(best);
  sys_phipi->Draw("chistsame");

  TH1F *sys_conv = (TH1F*)best->Clone();
  sys_conv->SetLineColor(7);
  sys_conv->SetLineWidth(2);
  sys_conv->Add(pteRC,-0.100);
  sys_conv->Divide(best);
  sys_conv->Draw("chistsame");

  // TH1F *sys_ke3 = (TH1F*)best->Clone();
  // sys_ke3->SetLineColor(8);
  // sys_ke3->SetLineWidth(2);
  // sys_ke3->Add(pteRKe3,-0.3);
  // sys_ke3->Divide(best);
  // sys_ke3->Draw("chistsame");

  TH1F *sys_direct = (TH1F*)best->Clone();
  sys_direct->SetLineColor(9);
  sys_direct->SetLineWidth(2);
  sys_direct->Add(pteRGamma,-0.5);
  sys_direct->Add(pteRCGamma,-0.5);
  sys_direct->Divide(best);
  sys_direct->Draw("chistsame");

  TH1F *sys_total = (TH1F*)best->Clone();
  sys_total->SetName("sys_total");
  sys_total->SetTitle("sys_total");
  double syserr, relerr;
  for (int ibin=1; ibin<=sys_pishape->GetNbinsX(); ibin++){
    syserr  = (1.-sys_pishape->GetBinContent(ibin))**2;
    syserr += (0.065)**2;
    syserr += (1.-sys_etapi->GetBinContent(ibin))**2;
    syserr += (1.-sys_rhopi->GetBinContent(ibin))**2;
    syserr += (1.-sys_omegapi->GetBinContent(ibin))**2;
    syserr += (1.-sys_etaprimepi->GetBinContent(ibin))**2;
    syserr += (1.-sys_phipi->GetBinContent(ibin))**2;
    syserr += (1.-sys_conv->GetBinContent(ibin))**2;
    // syserr += (1.-sys_ke3->GetBinContent(ibin))**2;
    syserr += (1.-sys_direct->GetBinContent(ibin))**2;
    syserr  = sqrt(syserr); 
    sys_total->SetBinContent(ibin,syserr);
    relerr = sys_pishape->GetBinError(ibin)/sys_pishape->GetBinContent(ibin);
    sys_total->SetBinError(ibin,syserr*relerr);
    cout << ibin << " " << syserr << " " << relerr << endl;
  }
  sys_total->SetLineColor(1);
  sys_total->SetLineWidth(2);
  sys_total->Draw("csamehist");
  // sys_total->Fit("pol0");


  return;

}
