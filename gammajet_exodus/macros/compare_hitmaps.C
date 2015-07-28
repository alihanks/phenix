void compare_hitmaps(){
  //gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  //TFile *f1=new TFile("10m_pi0_smeared_seff_tag_momw_pbglswap.root");
  //TFile *f1=new TFile("1m_pi0_erttest_1.root");
  //TFile *f1=new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_cent0_etas.root");
  //TFile *f1=new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/cent0_etas_all.root");
  //TFile *f1 = TFile::Open("/phenix/hp/data10/mjuszkie/exodus/cent0/meandmatt1.root");
  //TFile *f1 = TFile::Open("output_merged/fe/test2.root");
  TFile *f1 = TFile::Open("/phenix/hhj/jfrantz/gammajet_exodus/fe_sharks/cent0_II_0.root");

  //  TFile *f2=new TFile("$data10/run5/fg_ert/run5pp_fg_histo_inc_allfiles_modmapaxesswapped_052407_merged.root");
  //TFile *f2=new TFile("/phenix/data10/manguyen/run5/fg_ert/run5pp_fg_histo_inc_toweriddist_053007_merged.root");
  //  TFile *f2=new TFile("/phenix/data10/manguyen/run5/fg_ert/run5pp_fg_histo_pi0_thresh4_partnerhot_eitherERT_053107_merged.root");
  //TFile *f2=new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4merged/r4_a2w.root");
  //TFile *f2=new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r7qm09_modmaps/rg3.root");
  //TFile *f2=new TFile("/phenix/hp/data10/mjuszkie/run7/ppg113/taxi230/merged/passI_all.root");
  //TFile *f2=new TFile("/phenix/hp/data10/mjuszkie/run7/ppg113/taxi208/merged/all_modmaps.root");


  TFile *f2=new TFile("/phenix/hhj/jfrantz/gammajet_exodus/bingr8/mrg1.root");
  
  TH2F *SIM[8];
  TH3F *DAT[8], *SUB[8];
  for(int i=0;i<8;i++){
    char simname[100], datname[100], subname[100];
    sprintf(simname,"hitdiag_sec%d",i);
    sprintf(datname,"MODMAP%d",i);
    sprintf(subname,"MODMAPSUB%d",i);
    //sprintf(datname,"piC0_MODMAP%d",i);
    //sprintf(subname,"piC0_MODMAPSUB%d",i);
    SIM[i]=(TH2F*)f1->Get(simname);
    DAT[i]=(TH3F*)f2->Get(datname);
    SUB[i]=(TH3F*)f2->Get(subname);
  }

  TCanvas *c1=new TCanvas("c1","sim",1);
  c1->Divide(2,4);
  for(int i=0;i<8;i++){
    c1->cd(i+1);
    SIM[i]->Draw("zcol");
  }

  return;

  TCanvas *c2=new TCanvas("c2","dat",1);
  c2->Divide(2,4);
  for(int i=0;i<8;i++){
    c2->cd(i+1);
    //c2->GetPad(i+1)->SetLogz();
    //    DAT[i]->Divide(SIM[i]);
    //DAT[i]->GetZaxis()->SetRangeUser(50,150);
    DAT[i]->Project3D("yx")->Draw("zcol");
    //DAT[i]->Draw("zcol");

    

  }

  TCanvas *c3=new TCanvas("c3","sub",1);
  c3->Divide(2,4);
  for(int i=0;i<8;i++){
    c3->cd(i+1);
    //SUB[i]->Draw("zcol");
    SUB[i]->Project3D("yx")->Draw("zcol");

  }

  
}
