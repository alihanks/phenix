{
TCanvas *c1=new TCanvas("c1","c1",1);

  TFile *fdatagl = TFile::Open("/direct/phenix+hp/data10/mjuszkie/run7/nov172008/mergedhistos/new_calibs/pidatawidths_run7_pbgl_C1.root");
  widvspt->SetMarkerStyle(8);
  widvspt->SetMarkerColor(4);
  widvspt->SetLineColor(4);
  widvspt->Draw();

TFile *fsimboth = TFile::Open("pisimwidths_run7_both_C1.root");
  widvspt->SetMarkerStyle(24);
  widvspt->SetMarkerColor(1);
  widvspt->SetLineColor(1);
  widvspt->Draw("same");
  TFile *fsimgl = TFile::Open("pisimwidths_run7_pbgl_C1.root");
  widvspt->SetMarkerStyle(24);
  widvspt->SetMarkerColor(4);
  widvspt->SetLineColor(4);
  widvspt->Draw("same");
  TFile *fsimsc = TFile::Open("pisimwidths_run7_pbsc_C1.root");
  widvspt->SetMarkerStyle(24);
  widvspt->SetMarkerColor(2);
  widvspt->SetLineColor(2);
  widvspt->Draw("same");

TFile *fdataboth = TFile::Open("/direct/phenix+hp/data10/mjuszkie/run7/nov172008/mergedhistos/new_calibs/pidatawidths_run7_both_C1.root");
  widvspt->SetMarkerStyle(8);
  widvspt->SetMarkerColor(1);
  widvspt->SetLineColor(1);
  widvspt->Draw("same");

  TFile *fdatasc = TFile::Open("/direct/phenix+hp/data10/mjuszkie/run7/nov172008/mergedhistos/new_calibs/pidatawidths_run7_pbsc_C1.root");
  widvspt->SetMarkerStyle(8);
  widvspt->SetMarkerColor(2);
  widvspt->SetLineColor(2);
  widvspt->Draw("same");


TCanvas *c2=new TCanvas("c2","c2",1);
  fdatagl->cd();
  massvspt->SetMarkerStyle(8);
  massvspt->SetMarkerColor(4);
  massvspt->SetLineColor(4);
  massvspt->Draw();

fsimboth->cd();
  massvspt->SetMarkerStyle(24);
  massvspt->SetMarkerColor(1);
  massvspt->SetLineColor(1);
  massvspt->Draw("same");
  fsimgl->cd();
  massvspt->SetMarkerStyle(24);
  massvspt->SetMarkerColor(4);
  massvspt->SetLineColor(4);
  massvspt->Draw("same");
  fsimsc->cd();
  massvspt->SetMarkerStyle(24);
  massvspt->SetMarkerColor(2);
  massvspt->SetLineColor(2);
  massvspt->Draw("same");

fdataboth->cd();
  massvspt->SetMarkerStyle(8);
  massvspt->SetMarkerColor(1);
  massvspt->SetLineColor(1);
  massvspt->Draw("same");

  fdatasc->cd();
  massvspt->SetMarkerStyle(8);
  massvspt->SetMarkerColor(2);
  massvspt->SetLineColor(2);
  massvspt->Draw("same");

}
