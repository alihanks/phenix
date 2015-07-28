void Post_Exodus(char *filename)
{
  gStyle->SetOptStat(0);

  exodusfile=new TFile(filename,"READ");
  exodusfile->cd();
  
  TCanvas *c1 = new TCanvas("c1","Input Particles Spectrum");
  c1->SetLogy();

  TH1D *pi0pt=realpi0pt->Clone();
  TH1D *etapt=realetapt->Clone();
  TH1D *etaprimept=realetaprimept->Clone();
  TH1D *omegapt=realomegapt->Clone();

  pi0pt->SetMarkerStyle(20);
  pi0pt->SetMarkerColor(kBlack);
  
  etapt->SetMarkerStyle(21);
  etapt->SetMarkerColor(kRed);

  etaprimept->SetMarkerStyle(22);
  etaprimept->SetMarkerColor(kBlue);

  omegapt->SetMarkerStyle(23);
  omegapt->SetMarkerColor(kGreen);

  pi0pt->Rebin(10);
  etapt->Rebin(10);
  etaprimept->Rebin(10);
  omegapt->Rebin(10);

  etapt->Scale(0.1);
  etaprimept->Scale(0.01);
  omegapt->Scale(0.001);

  TAxis *yaxis=pi0pt->GetYaxis();
  yaxis->SetTitle("dN/dp_{T}dy");
  TAxis *xaxis=pi0pt->GetXaxis();
  xaxis->SetTitle("p_{T} [Gev/c]");

  pi0pt->Draw("e");
  etapt->Draw("esame");
  etaprimept->Draw("esame");
  omegapt->Draw("esame");

  TLatex *pi0=new TLatex(8,5000,"Min-Bias: #pi^{0},");
  pi0->SetTextColor(kBlack);

  TLatex *eta=new TLatex(16,5000,"#eta,");
  eta->SetTextColor(kRed);
  
  TLatex *etaprime=new TLatex(18,5000,"#eta^{'},");
  etaprime->SetTextColor(kBlue);

  TLatex *omega=new TLatex(20,5000,"#omega");
  omega->SetTextColor(kGreen);

  pi0->Draw();
  eta->Draw();
  etaprime->Draw();
  omega->Draw();

  TCanvas *c2 = new TCanvas("c2","eta to pi0 ratio");
  
  TH1D *etaratio = realetapt->Clone();
  
  etaratio->Rebin(10);
  etaratio->Divide(pi0pt);
  
  etaratio->SetMarkerStyle(21);
  etaratio->SetMarkerColor(kRed);

  TH1D *etaprimeratio = realetaprimept->Clone();
  
  etaprimeratio->SetMarkerStyle(22);
  etaprimeratio->SetMarkerColor(kBlue);

  etaprimeratio->Rebin(10);
  etaprimeratio->Divide(pi0pt);

  TH1D *omegaratio = realomegapt->Clone();
  omegaratio->SetMarkerStyle(23);
  omegaratio->SetMarkerColor(kGreen);

  omegaratio->Rebin(10);
  omegaratio->Divide(pi0pt);

  TAxis *yaxis=etaprimeratio->GetYaxis();
  yaxis->SetTitle("#eta(#eta^{'},#omega) / #pi^{0}");
  TAxis *xaxis=etaprimeratio->GetXaxis();
  xaxis->SetTitle("p_{T} [Gev/c]");

  TLegend *leg = new TLegend(0.5,0.2,0.6,0.4);
  
  char dummy[500];
  
  sprintf(dummy,"#eta/#pi^{0}");
  leg->AddEntry(etaratio,dummy,"p");

  sprintf(dummy,"#eta^{'}/#pi^{0}");
  leg->AddEntry(etaprimeratio,dummy,"p");
  
  sprintf(dummy,"#omega/#pi^{0}");
  leg->AddEntry(omegaratio,dummy,"p");
  
  etaprimeratio->SetTitle("#eta, #eta^{'}, #omega ratio to #pi^{0}");
  etaprimeratio->SetAxisRange(0.0,15.0);
  etaprimeratio->Draw("e");
  leg->Draw();
  omegaratio->Draw("esame");
  etaratio->Draw("esame");

  TCanvas *c3 = new TCanvas("c3","MC ratio break-down");
  c3->SetLogy();
  
  TH1D *allgamma = gammapt->Clone();

  allgamma->Rebin(10);
  pi0gamma->Rebin(10);
  etagamma->Rebin(10);
  etaprimegamma->Rebin(10);
  omegagamma->Rebin(10);

  TH1D *taggedgamma = pi0gamma->Clone();

  allgamma->Divide(pi0pt);
  pi0gamma->Divide(pi0pt);
  etagamma->Divide(pi0pt);
  etaprimegamma->Divide(pi0pt);
  omegagamma->Divide(pi0pt);

  allgamma->SetLineColor(kBlack);
  pi0gamma->SetLineColor(kRed);
  etagamma->SetLineColor(kBlue);
  etaprimegamma->SetLineColor(kGreen);
  omegagamma->SetLineColor(6);

  TAxis *yaxis=allgamma->GetYaxis();
  yaxis->SetTitle("#gamma / #pi^{0}");
  TAxis *xaxis=allgamma->GetXaxis();
  xaxis->SetTitle("p_{T} [Gev/c]");

  
  allgamma->Draw("e");
  pi0gamma->Draw("esame");
  etagamma->Draw("esame");
  etaprimegamma->Draw("esame");
  omegagamma->Draw("esame");

  TCanvas *c4 = new TCanvas("c4","Our Ratio");

  TH1D *allgamma2 = gammapt->Clone();

  allgamma2->Rebin(10);
  allgamma2->Divide(taggedgamma);

  for(int i = 1; i <= allgamma2->GetNbinsX(); i++)
    allgamma2->SetBinError(i,allgamma2->GetBinError(i)/sqrt(2));

  TAxis *yaxis=allgamma2->GetYaxis();
  yaxis->SetTitle("#gamma / #pi^{0}(#gamma)");
  TAxis *xaxis=allgamma2->GetXaxis();
  xaxis->SetTitle("p_{T} [Gev/c]");  

  TLatex *ourratio=new TLatex(6.0, 1.125, "Min-Bias (N(#gamma_{1}^{hadron})/N(#gamma_{#pi^{0}}))_{MC}");
  
  allgamma2->SetTitle("#gamma_{1}^{hadron}/#gamma_{#pi^{0}} in full space");
  allgamma2->SetAxisRange(0.0,15.0);

  allgamma2->SetMarkerStyle(20);
  allgamma2->Draw();
  ourratio->Draw();
}
