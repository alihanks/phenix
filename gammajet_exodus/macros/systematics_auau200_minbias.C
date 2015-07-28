void systematics_auau200_minbias(){
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

  TFile *fs = new TFile("cocktail_10M_auau200_minbias_final_soft.root");
  TH1D *soft = (TH1D*) fs->Get("pteR");
  soft->SetLineColor(2);

  TFile *fc = new TFile("cocktail_10M_auau200_minbias_final.root");
  TH1D *best = (TH1D*) fc->Get("pteR");

  TH1F *sys_pishape = (TH1F*)soft->Clone();
  sys_pishape->SetLineColor(1);
  sys_pishape->SetLineWidth(2);
  sys_pishape->Divide(best);
  sys_pishape->Draw("chistsame");

  TH1F *sys_etapi = (TH1F*)best->Clone();
  sys_etapi->SetLineColor(2);
  sys_etapi->SetLineWidth(2);
  sys_etapi->Add(pteREta,0.22);
  sys_etapi->Add(pteRCEta,0.22);
  sys_etapi->Divide(best);
  sys_etapi->Draw("chistsame");

  TH1F *sys_rhopi = (TH1F*)best->Clone();
  sys_rhopi->SetLineColor(3);
  sys_rhopi->SetLineWidth(2);
  sys_rhopi->Add(pteRRho,0.30);
  sys_rhopi->Divide(best);
  sys_rhopi->Draw("chistsame");

  TH1F *sys_omegapi = (TH1F*)best->Clone();
  sys_omegapi->SetLineColor(4);
  sys_omegapi->SetLineWidth(2);
  sys_omegapi->Add(pteROmega,0.30);
  sys_omegapi->Divide(best);
  sys_omegapi->Draw("chistsame");

  TH1F *sys_etaprimepi = (TH1F*)best->Clone();
  sys_etaprimepi->SetLineColor(5);
  sys_etaprimepi->SetLineWidth(2);
  sys_etaprimepi->Add(pteREtaprime,0.30);
  sys_etaprimepi->Add(pteRCEtaprime,0.30);
  sys_etaprimepi->Divide(best);
  sys_etaprimepi->Draw("chistsame");

  TH1F *sys_phipi = (TH1F*)best->Clone();
  sys_phipi->SetLineColor(6);
  sys_phipi->SetLineWidth(2);
  sys_phipi->Add(pteRPhi,0.30);
  sys_phipi->Divide(best);
  sys_phipi->Draw("chistsame");

  TH1F *sys_conv = (TH1F*)best->Clone();
  sys_conv->SetLineColor(7);
  sys_conv->SetLineWidth(2);
  sys_conv->Add(pteRC,0.100);
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
  sys_direct->Add(pteRGamma,0.5);
  sys_direct->Add(pteRCGamma,0.5);
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

  c2 = new TCanvas();
  c2->SetGridx(1);
  c2->SetGridy(1);
  c2->SetTicks(1,1);
  c2->Draw();
  
  TH2F *hist_background2 = new TH2F("hist_background2","Non-photonic electron invariant cross-section",10,0.,4.5,10,0.,1.2);
  hist_background2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_background2->GetYaxis()->SetTitle("rel. sys. error");
  hist_background2->GetYaxis()->SetTitleOffset(1.2); 
  hist_background2->Draw();

  TFile *fh2 = new TFile("cocktail_10M_auau200_minbias_final_hard.root");
  TH1D *hard2 = (TH1D*) fh2->Get("pteR");
  hard2->SetLineColor(2);

  TFile *fs2 = new TFile("cocktail_10M_auau200_minbias_final_soft.root");
  TH1D *soft2 = (TH1D*) fs2->Get("pteR");
  soft2->SetLineColor(2);

  TFile *fc2 = new TFile("cocktail_10M_auau200_minbias_final.root");
  TH1D *best2 = (TH1D*) fc2->Get("pteR");

  TH1F *sys_pishape2 = (TH1F*)hard2->Clone();
  sys_pishape2->SetLineColor(1);
  sys_pishape2->SetLineWidth(2);
  sys_pishape2->Divide(best2);
  sys_pishape2->Draw("chistsame");

  TH1F *sys_etapi2 = (TH1F*)best2->Clone();
  sys_etapi2->SetLineColor(2);
  sys_etapi2->SetLineWidth(2);
  sys_etapi2->Add(pteREta,-0.22);
  sys_etapi2->Add(pteRCEta,-0.22);
  sys_etapi2->Divide(best2);
  sys_etapi2->Draw("chistsame");

  TH1F *sys_rhopi2 = (TH1F*)best2->Clone();
  sys_rhopi2->SetLineColor(3);
  sys_rhopi2->SetLineWidth(2);
  sys_rhopi2->Add(pteRRho,-0.30);
  sys_rhopi2->Divide(best2);
  sys_rhopi2->Draw("chistsame");

  TH1F *sys_omegapi2 = (TH1F*)best2->Clone();
  sys_omegapi2->SetLineColor(4);
  sys_omegapi2->SetLineWidth(2);
  sys_omegapi2->Add(pteROmega,-0.30);
  sys_omegapi2->Divide(best2);
  sys_omegapi2->Draw("chistsame");

  TH1F *sys_etaprimepi2 = (TH1F*)best2->Clone();
  sys_etaprimepi2->SetLineColor(5);
  sys_etaprimepi2->SetLineWidth(2);
  sys_etaprimepi2->Add(pteREtaprime,-0.30);
  sys_etaprimepi2->Add(pteRCEtaprime,-0.30);
  sys_etaprimepi2->Divide(best2);
  sys_etaprimepi2->Draw("chistsame");

  TH1F *sys_phipi2 = (TH1F*)best2->Clone();
  sys_phipi2->SetLineColor(6);
  sys_phipi2->SetLineWidth(2);
  sys_phipi2->Add(pteRPhi,-0.30);
  sys_phipi2->Divide(best2);
  sys_phipi2->Draw("chistsame");

  TH1F *sys_conv2 = (TH1F*)best2->Clone();
  sys_conv2->SetLineColor(7);
  sys_conv2->SetLineWidth(2);
  sys_conv2->Add(pteRC,-0.100);
  sys_conv2->Divide(best2);
  sys_conv2->Draw("chistsame");

  // TH1F *sys_ke32 = (TH1F*)best2->Clone();
  // sys_ke32->SetLineColor(8);
  // sys_ke32->SetLineWidth(2);
  // sys_ke32->Add(pteRKe3,-0.3);
  // sys_ke32->Divide(best2);
  // sys_ke32->Draw("chistsame");

  TH1F *sys_direct2 = (TH1F*)best2->Clone();
  sys_direct2->SetLineColor(9);
  sys_direct2->SetLineWidth(2);
  sys_direct2->Add(pteRGamma,-0.5);
  sys_direct2->Add(pteRCGamma,-0.5);
  sys_direct2->Divide(best2);
  sys_direct2->Draw("chistsame");

  TH1F *sys_total2 = (TH1F*)best2->Clone();
  sys_total2->SetName("sys_total2");
  sys_total2->SetTitle("sys_total2");
  double syserr2, relerr2;
  for (int ibin=1; ibin<=sys_pishape2->GetNbinsX(); ibin++){
    syserr2  = (1.-sys_pishape2->GetBinContent(ibin))**2;
    syserr2 += (0.065)**2;
    syserr2 += (1.-sys_etapi2->GetBinContent(ibin))**2;
    syserr2 += (1.-sys_rhopi2->GetBinContent(ibin))**2;
    syserr2 += (1.-sys_omegapi2->GetBinContent(ibin))**2;
    syserr2 += (1.-sys_etaprimepi2->GetBinContent(ibin))**2;
    syserr2 += (1.-sys_phipi2->GetBinContent(ibin))**2;
    syserr2 += (1.-sys_conv2->GetBinContent(ibin))**2;
    // syserr2 += (1.-sys_ke32->GetBinContent(ibin))**2;
    syserr2 += (1.-sys_direct2->GetBinContent(ibin))**2;
    syserr2  = sqrt(syserr2); 
    sys_total2->SetBinContent(ibin,syserr2);
    relerr2 = sys_pishape2->GetBinError(ibin)/sys_pishape2->GetBinContent(ibin);
    sys_total2->SetBinError(ibin,syserr2*relerr2);
    cout << ibin << " " << syserr2 << " " << relerr2 << endl;
  }
  sys_total2->SetLineColor(1);
  sys_total2->SetLineWidth(2);
  sys_total2->Draw("csamehist");

  c3 = new TCanvas();
  c3->SetGridx(1);
  c3->SetGridy(1);
  c3->SetTicks(1,1);
  c3->Draw();
  
  TH2F *hist_background3 = new TH2F("hist_background3","Non-photonic electron invariant cross-section",10,0.,4.5,10,0.,0.2);
  hist_background3->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_background3->GetYaxis()->SetTitle("rel. sys. error");
  hist_background3->GetYaxis()->SetTitleOffset(1.2); 
  hist_background3->Draw();

  TH1F *sys_final = (TH1F*)sys_total->Clone();
  sys_final->SetName("sys_final");
  sys_final->SetTitle("sys_final");
  double syserrf, relerrf;
  sys_total->SetLineColor(2);
  sys_total->Draw("csamehist");
  sys_total2->SetLineColor(4);
  sys_total2->Draw("csamehist");
  for (int ibin=1; ibin<=sys_total->GetNbinsX(); ibin++){
    syserrf = (sys_total->GetBinContent(ibin) + 
	       sys_total2->GetBinContent(ibin))/2.0;
    relerrf = (sys_total->GetBinError(ibin) + 
	       sys_total2->GetBinError(ibin))/2.0;
    sys_final->SetBinContent(ibin,syserrf);
    sys_final->SetBinError(ibin,relerrf);
  }
  sys_final->SetLineColor(1);
  sys_final->Draw("csamehist");

  return;

}
