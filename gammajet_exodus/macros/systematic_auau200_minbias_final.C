void systematic_auau200_minbias_final(){
  gROOT->Reset();
    
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(1.1);
  gStyle->SetFillColor(0);

  c1 = new TCanvas();
  c1->SetFillColor(10);
  c1->SetGridx(1);
  c1->SetGridy(1);
  c1->SetTicks(1,1);
  c1->Draw();
  
  TH2F *hist_background = new TH2F("hist_background","Non-photonic electron invariant cross-section",10,0.,5.,10,0..,20.);
  hist_background->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_background->GetYaxis()->SetTitle("sys. error [%]");
  hist_background->GetYaxis()->SetTitleOffset(1.2); 
  hist_background->Draw();

  TFile *fh = new TFile("cocktail_10M_auau200_40_60_final_update_high.root");
  TH1D *high = (TH1D*) fh->Get("pteR");
  high->SetLineColor(2);

  TFile *fs = new TFile("cocktail_10M_auau200_40_60_final_update_low.root");
  TH1D *low = (TH1D*) fs->Get("pteR");
  low->SetLineColor(2);

  TFile *fc = new TFile("cocktail_10M_auau200_40_60_final_update.root");
  TH1D *best = (TH1D*) fc->Get("pteR");

  TH1F *sys_lo_pi_input = (TH1F*)low->Clone();
  sys_lo_pi_input->SetLineColor(1);
  sys_lo_pi_input->SetLineWidth(2);
  sys_lo_pi_input->Add(best,-1.);
  sys_lo_pi_input->Divide(best);
  sys_lo_pi_input->Scale(100.);
  //sys_lo_pi_input->Draw("chistsame");

  TH1F *sys_lo_etapi = (TH1F*)best->Clone();
  sys_lo_etapi->SetLineColor(2);
  sys_lo_etapi->SetLineWidth(2);
  sys_lo_etapi->Add(pteREta,-0.22);
  sys_lo_etapi->Add(pteRCEta,-0.22);
  sys_lo_etapi->Add(best,-1.);
  sys_lo_etapi->Divide(best);
  sys_lo_etapi->Scale(100.);
  //sys_lo_etapi->Draw("chistsame");

  TH1F *sys_lo_rhopi = (TH1F*)best->Clone();
  sys_lo_rhopi->SetLineColor(3);
  sys_lo_rhopi->SetLineWidth(2);
  sys_lo_rhopi->Add(pteRRho,-0.30);
  sys_lo_rhopi->Add(best,-1.);
  sys_lo_rhopi->Divide(best);
  sys_lo_rhopi->Scale(100.);
  //sys_lo_rhopi->Draw("chistsame");

  TH1F *sys_lo_omegapi = (TH1F*)best->Clone();
  sys_lo_omegapi->SetLineColor(4);
  sys_lo_omegapi->SetLineWidth(2);
  sys_lo_omegapi->Add(pteROmega,-0.30);
  sys_lo_omegapi->Add(best,-1.);
  sys_lo_omegapi->Divide(best);
  sys_lo_omegapi->Scale(100.);
  //sys_lo_omegapi->Draw("chistsame");

  TH1F *sys_lo_etaprimepi = (TH1F*)best->Clone();
  sys_lo_etaprimepi->SetLineColor(5);
  sys_lo_etaprimepi->SetLineWidth(2);
  sys_lo_etaprimepi->Add(pteREtaprime,-0.30);
  sys_lo_etaprimepi->Add(pteRCEtaprime,-0.30);
  sys_lo_etaprimepi->Add(best,-1.);
  sys_lo_etaprimepi->Divide(best);
  sys_lo_etaprimepi->Scale(100.);
  //sys_lo_etaprimepi->Draw("chistsame");

  TH1F *sys_lo_phipi = (TH1F*)best->Clone();
  sys_lo_phipi->SetLineColor(6);
  sys_lo_phipi->SetLineWidth(2);
  sys_lo_phipi->Add(pteRPhi,-0.30);
  sys_lo_phipi->Add(best,-1.);
  sys_lo_phipi->Divide(best);
  sys_lo_phipi->Scale(100.);
  //sys_lo_phipi->Draw("chistsame");

  TH1F *sys_lo_conv = (TH1F*)best->Clone();
  sys_lo_conv->SetLineColor(7);
  sys_lo_conv->SetLineWidth(2);
  sys_lo_conv->Add(pteRC,-0.091);
  sys_lo_conv->Add(best,-1.);
  sys_lo_conv->Divide(best);
  sys_lo_conv->Scale(100.);
  //sys_lo_conv->Draw("chistsame");

  TH1F *sys_lo_ke3 = (TH1F*)best->Clone();
  sys_lo_ke3->SetLineColor(8);
  sys_lo_ke3->SetLineWidth(2);
  sys_lo_ke3->Add(pteRKe3,-0.2);
  sys_lo_ke3->Add(best,-1.);
  sys_lo_ke3->Divide(best);
  sys_lo_ke3->Scale(100.);
  //sys_lo_ke3->Draw("chistsame");

  TH1F *sys_lo_direct = (TH1F*)best->Clone();
  sys_lo_direct->SetLineColor(9);
  sys_lo_direct->SetLineWidth(2);
  sys_lo_direct->Add(pteRGamma,-0.52);
  sys_lo_direct->Add(pteRCGamma,-0.52);
  sys_lo_direct->Add(best,-1.);
  sys_lo_direct->Divide(best);
  sys_lo_direct->Scale(100.);
  //sys_lo_direct->Draw("chistsame");

  TH1F *sys_lo_total = (TH1F*)best->Clone();
  sys_lo_total->SetName("sys_lo_total");
  sys_lo_total->SetTitle("sys_lo_total");
  double syserr, relerr;
  for (int ibin=1; ibin<=sys_lo_pi_input->GetNbinsX(); ibin++){
    syserr  = 0.0;
    syserr += (sys_lo_pi_input->GetBinContent(ibin))**2;
    syserr += (sys_lo_etapi->GetBinContent(ibin))**2;
    syserr += (sys_lo_rhopi->GetBinContent(ibin))**2;
    syserr += (sys_lo_omegapi->GetBinContent(ibin))**2;
    syserr += (sys_lo_etaprimepi->GetBinContent(ibin))**2;
    syserr += (sys_lo_phipi->GetBinContent(ibin))**2;
    syserr += (sys_lo_conv->GetBinContent(ibin))**2;
    syserr += (sys_lo_ke3->GetBinContent(ibin))**2;
    syserr += (sys_lo_direct->GetBinContent(ibin))**2;
    syserr  = sqrt(syserr); 
    sys_lo_total->SetBinContent(ibin,syserr);
    relerr = sys_lo_pi_input->GetBinError(ibin)/sys_lo_pi_input->GetBinContent(ibin);
    sys_lo_total->SetBinError(ibin,syserr*relerr);
  }
  sys_lo_total->SetLineColor(4);
  sys_lo_total->SetLineWidth(2);
  //sys_lo_total->Draw("csamehist");

  TH1F *sys_hi_pi_input = (TH1F*)high->Clone();
  sys_hi_pi_input->SetLineColor(1);
  sys_hi_pi_input->SetLineWidth(2);
  sys_hi_pi_input->Add(best,-1.);
  sys_hi_pi_input->Divide(best);
  sys_hi_pi_input->Scale(100.);
  //sys_hi_pi_input->Draw("chistsame");

  TH1F *sys_hi_etapi = (TH1F*)best->Clone();
  sys_hi_etapi->SetLineColor(2);
  sys_hi_etapi->SetLineWidth(2);
  sys_hi_etapi->Add(pteREta,0.22);
  sys_hi_etapi->Add(pteRCEta,0.22);
  sys_hi_etapi->Add(best,-1.);
  sys_hi_etapi->Divide(best);
  sys_hi_etapi->Scale(100.);
  //sys_hi_etapi->Draw("chistsame");

  TH1F *sys_hi_rhopi = (TH1F*)best->Clone();
  sys_hi_rhopi->SetLineColor(3);
  sys_hi_rhopi->SetLineWidth(2);
  sys_hi_rhopi->Add(pteRRho,0.30);
  sys_hi_rhopi->Add(best,-1.);
  sys_hi_rhopi->Divide(best);
  sys_hi_rhopi->Scale(100.);
  //sys_hi_rhopi->Draw("chistsame");

  TH1F *sys_hi_omegapi = (TH1F*)best->Clone();
  sys_hi_omegapi->SetLineColor(4);
  sys_hi_omegapi->SetLineWidth(2);
  sys_hi_omegapi->Add(pteROmega,0.30);
  sys_hi_omegapi->Add(best,-1.);
  sys_hi_omegapi->Divide(best);
  sys_hi_omegapi->Scale(100.);
  //sys_hi_omegapi->Draw("chistsame");

  TH1F *sys_hi_etaprimepi = (TH1F*)best->Clone();
  sys_hi_etaprimepi->SetLineColor(5);
  sys_hi_etaprimepi->SetLineWidth(2);
  sys_hi_etaprimepi->Add(pteREtaprime,0.30);
  sys_hi_etaprimepi->Add(pteRCEtaprime,0.30);
  sys_hi_etaprimepi->Add(best,-1.);
  sys_hi_etaprimepi->Divide(best);
  sys_hi_etaprimepi->Scale(100.);
  //sys_hi_etaprimepi->Draw("chistsame");

  TH1F *sys_hi_phipi = (TH1F*)best->Clone();
  sys_hi_phipi->SetLineColor(6);
  sys_hi_phipi->SetLineWidth(2);
  sys_hi_phipi->Add(pteRPhi,0.30);
  sys_hi_phipi->Add(best,-1.);
  sys_hi_phipi->Divide(best);
  sys_hi_phipi->Scale(100.);
  //sys_hi_phipi->Draw("chistsame");

  TH1F *sys_hi_conv = (TH1F*)best->Clone();
  sys_hi_conv->SetLineColor(7);
  sys_hi_conv->SetLineWidth(2);
  sys_hi_conv->Add(pteRC,0.091);
  sys_hi_conv->Add(best,-1.);
  sys_hi_conv->Divide(best);
  sys_hi_conv->Scale(100.);
  //sys_hi_conv->Draw("chistsame");

  TH1F *sys_hi_ke3 = (TH1F*)best->Clone();
  sys_hi_ke3->SetLineColor(8);
  sys_hi_ke3->SetLineWidth(2);
  sys_hi_ke3->Add(pteRKe3,0.2);
  sys_hi_ke3->Add(best,-1.);
  sys_hi_ke3->Divide(best);
  sys_hi_ke3->Scale(100.);
  //sys_hi_ke3->Draw("chistsame");

  TH1F *sys_hi_direct = (TH1F*)best->Clone();
  sys_hi_direct->SetLineColor(9);
  sys_hi_direct->SetLineWidth(2);
  sys_hi_direct->Add(pteRGamma,0.52);
  sys_hi_direct->Add(pteRCGamma,0.52);
  sys_hi_direct->Add(best,-1.);
  sys_hi_direct->Divide(best);
  sys_hi_direct->Scale(100.);
  //sys_hi_direct->Draw("chistsame");

  TH1F *sys_hi_total = (TH1F*)best->Clone();
  sys_hi_total->SetName("sys_hi_total");
  sys_hi_total->SetTitle("sys_hi_total");
  for (int ibin=1; ibin<=sys_hi_pi_input->GetNbinsX(); ibin++){
    syserr  = 0.0;
    syserr += (sys_hi_pi_input->GetBinContent(ibin))**2;
    syserr += (sys_hi_etapi->GetBinContent(ibin))**2;
    syserr += (sys_hi_rhopi->GetBinContent(ibin))**2;
    syserr += (sys_hi_omegapi->GetBinContent(ibin))**2;
    syserr += (sys_hi_etaprimepi->GetBinContent(ibin))**2;
    syserr += (sys_hi_phipi->GetBinContent(ibin))**2;
    syserr += (sys_hi_conv->GetBinContent(ibin))**2;
    syserr += (sys_hi_ke3->GetBinContent(ibin))**2;
    syserr += (sys_hi_direct->GetBinContent(ibin))**2;
    syserr  = sqrt(syserr); 
    sys_hi_total->SetBinContent(ibin,syserr);
    relerr = sys_hi_pi_input->GetBinError(ibin)/sys_hi_pi_input->GetBinContent(ibin);
    sys_hi_total->SetBinError(ibin,syserr*relerr);
    cout << ibin << " " << syserr << " " << relerr << endl;
  }
  sys_hi_total->SetLineColor(2);
  sys_hi_total->SetLineWidth(2);
  //sys_hi_total->Draw("csamehist");

  TH1F *sys_total = (TH1F*)sys_lo_total->Clone();
  sys_total->SetName("sys_total");
  sys_total->SetTitle("sys_total");
  double average, eaverage;
  for (int ibin=1; ibin<=sys_total->GetNbinsX(); ibin++){
    average = (sys_lo_total->GetBinContent(ibin) + 
	       sys_hi_total->GetBinContent(ibin))/2.;
    eaverage = sqrt((sys_lo_total->GetBinError(ibin))**2 + 
		    (sys_hi_total->GetBinError(ibin))**2)/2.;
    sys_total->SetBinContent(ibin,average);
    sys_total->SetBinError(ibin,eaverage);
    cout << ibin << " " << average << " " << eaverage << endl;
  }
  sys_total->SetLineColor(1);
  sys_total->SetLineWidth(2);
  sys_total->Draw("same");
  TF1 *fit_p2 = new TF1("fit_p2","9.5-1.1*x+0.52*x*x",0.,5.);
  fit_p2->SetLineColor(2);
  fit_p2->Draw("same");

  return;

}
