{
/*
TFile *fdata = TFile::Open("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4merged/r4all.root");
piC0_CENTRALITY->GetEntries();
piC0_CENTRALITY->Integral(1,20);
piC0_CENTRALITY->Integral(1,21);
piC0_TRIGPT->Scale(1./2.15444524000000000e+08);
piC0_TRIGPT->Draw();
167*500000
TFile *fsimeff = TFile::Open("/phenix/hp/data10/mjuszkie/exodus/cent0/all_hieff.root");
pi_TRIGPT->Scale(1./83500000.);
pi_TRIGPT->Draw("same");
pi_TRIGPT->Integral(51,150);
fdata->cd();;
piC0_TRIGPT->Integral(51,150);
5.07982276890750484e-04/3.45940119760478916e-04
C0_TRIGPT->Scale(1./2.15444524000000000e+08);
C0_TRIGPT->Draw();
fsimeff->cd();
TRIGPT->Scale(1./83500000.);
TRIGPT->Draw("same");
88*1.5
TRIGPT->Integral(51,150);
fdata->cd();;
C0_TRIGPT->Integral(51,150);
5.47848690324315157e-04/3.98131736526946068e-04
fsimeff->cd();
realpi0pt->GetEntries();
7.34815770300000000e+09/83500000.
realomegapt->GetEntries();
1.01405408000000000e+08/83500000.
realetapt->GetEntries();
1.39606672400000000e+09/83500000.
1.67193619640718580e+01/8.80018886586826312e+01
1.21443602395209571/8.80018886586826312e+01
*/



//TFile *_file0 = TFile::Open("/phenix/hp/data10/mjuszkie/exodus/cent0/all_posres.root");
//TFile *_file0 = TFile::Open("/phenix/hp/data10/mjuszkie/exodus/cent0/new/all_posresmymap.root");
//TFile *_file0 = TFile::Open("/phenix/hp/data10/mjuszkie/exodus/cent0/new/set1/all_set1.root");
//TFile *_file0 = TFile::Open("/phenix/hp/data10/mjuszkie/exodus/cent0/meandmatt1.root");
//TFile *_file0 = TFile::Open("/phenix/scratch/mjuszkie/gammajet.exodus/fe_sharks/cent0_bigIII.root");
TFile *_file0 = TFile::Open("output_merged/fe/new/all_fenew.root");
float simevts= realpi0pt->GetEntries()/88.0;
TCanvas *c1 = new TCanvas("c1","c1",1);
//TRIGPT->Draw();
//pi_TRIGPT->Draw();
pi_TRIGPT->Scale(1.0/simevts);
pi_TRIGPT->SetLineColor(4);
pi_TRIGPT->Draw();
float simpi=pi_TRIGPT->Integral(51,150);
//pi_TRIGPT->Scale(1./simpi);

//TFile *fdata = TFile::Open("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4merged/r4all.root");
//TFile *fdata = TFile::Open("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4alla.root");
TFile *fdata = TFile::Open("/phenix/u/workarea/mjuszkie/run7_QM09/taxi136/rg6.root");
float datevts=piC0_CENTRALITY->GetEntries();
piC0_TRIGPT->Scale(1.0/datevts);
float datapi=piC0_TRIGPT->Integral(51,150);
float sim_dat=simpi/datapi;
piC0_TRIGPT->Draw("same");
//C0_TRIGPT->Scale(1./2.15e+08);
TCanvas *c2 = new TCanvas("c2","c2",1);
itagC0_TRIGPT->Scale(1.0/datevts);
//itagC0_TRIGPT->Scale(sim_dat);
itagC0_TRIGPT->Draw();
float datagam=itagC0_TRIGPT->Integral(51,150);


_file0->cd();
TRIGPT_0->Scale(1.0/simevts);
TRIGPT_0->SetLineColor(4);
TRIGPT_0->Draw("same");
float simgam=TRIGPT_0->Integral(51,150);

cout << "datapi/simpi = " << datapi/simpi <<endl;
cout << "datagam/simgam = " << datagam/simgam <<endl;

cout << "datapi " << datapi << " simpi " << simpi <<endl;
cout << "datagam = " << datagam << " simgam " << simgam <<endl;

}
