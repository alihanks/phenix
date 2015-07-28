{

//TFile *fin =  new TFile("/phenix/hp/data10/manguyen/run6/fg_ert/run6pp_fg_histo_pi0_ert_taxi70_fidcut_pisacorrs_111808_merged.root");
TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4_b2w.root");
//TFile *fin =  new TFile("/phenix/hp/data06/jfrantz/merge_qm08/bf2/allhall.root");

TH2F *INVMASS = (TH2F*)fin->Get("itagC0_INVMASS_0");
/*
TH2F *htemp;

htemp = (TH2F *)fin->Get("C_1__INVMASS_1");
INVMASS->Add(htemp);
htemp = (TH2F *)fin->Get("C_1__INVMASS_2");
INVMASS->Add(htemp);
htemp = (TH2F *)fin->Get("C_1__INVMASS_3");
INVMASS->Add(htemp);
htemp = (TH2F *)fin->Get("C_1__INVMASS_4");
INVMASS->Add(htemp);
htemp = (TH2F *)fin->Get("C_1__INVMASS_5");
INVMASS->Add(htemp);
*/
//htemp = (TH2F *)fin->Get("C_0__INVMASS_6");
//INVMASS->Add(htemp);
//htemp = (TH2F *)fin->Get("C_0__INVMASS_7");
//INVMASS->Add(htemp);
/*


//TH2F *INVMASS = (TH2F*)fin->Get("piC0_INVMASS_0");
TH2F *INVMASS = (TH2F*)fin->Get("piC0_INVMASS_6");
 TH2F *htemp;
// htemp = (TH2F *)fin->Get("piC0_INVMASS_1");
// INVMASS->Add(htemp);
// htemp = (TH2F *)fin->Get("piC0_INVMASS_2");
// INVMASS->Add(htemp);
// htemp = (TH2F *)fin->Get("piC0_INVMASS_3");
// INVMASS->Add(htemp);
// htemp = (TH2F *)fin->Get("piC0_INVMASS_4");
// INVMASS->Add(htemp);
// htemp = (TH2F *)fin->Get("piC0_INVMASS_5");
// INVMASS->Add(htemp);
// htemp = (TH2F *)fin->Get("piC0_INVMASS_6");
// INVMASS->Add(htemp);
htemp = (TH2F *)fin->Get("piC0_INVMASS_7");
INVMASS->Add(htemp);
*/
TH1D * hmass = INVMASS->ProjectionY();
hmass->Draw();


TH1F * massvspt = new TH1F("massvspt","",10, 0, 10);
TH1F * widvspt = new TH1F("widvspt","",10, 0, 10);

TF1 * ff2 = new TF1("ff2","gaus(0) + pol1(3)",0,1);
ff2->SetParameter(0,1e4);
ff2->SetParameter(1,0.135);
ff2->SetParameter(2,0.010);
ff2->SetParameter(3,10);


TCanvas *c=new TCanvas("c","c",1);
c->Divide(2,3);

for (int i = 0; i < 6; i++)
{
  c->cd(i+1);
  char name[100];
  sprintf(name,"h%d",i);
  hmass = INVMASS->ProjectionY(name,40+i*10,50+i*10);
  hmass->Fit("ff2","","",0.07,0.25);
  TF1 * fgaus = hmass->GetFunction("ff2");
  float pkpos = fgaus->GetParameter(1);
  float poserr = fgaus->GetParError(1);
  float pkwid = fgaus->GetParameter(2);
  float widerr = fgaus->GetParError(2);
  cout << "pt bin " << 5+i << " pos " << pkpos << " wid " << pkwid << endl;
  massvspt->SetBinContent(5+i, pkpos);
  massvspt->SetBinError(5+i,poserr);
  widvspt->SetBinContent(5+i, pkwid);
  widvspt->SetBinError(5+i,widerr);  

}


TFile * mf2 = new TFile("pidatawidths_run4_tag0.root","RECREATE");
massvspt->Write();
widvspt->Write();



}
