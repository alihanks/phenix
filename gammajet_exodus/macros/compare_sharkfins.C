void compare_sharkfins(){
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  //TFile *f1=new TFile("sharkfin_projections_single_notmiss.root");
  //TFile *f1=new TFile("sharkfin_projections_singII_miss_sum.root");
  //TFile *f1=new TFile("sharkfin_projections_r7map1_notmiss_sum.root");
  //TFile *f1=new TFile("/phenix/u/workarea/manguyen/gammajet.exodus/histat_zbinned_pi0_sharkfinsproj_run7AA.root");
  TFile *f1=new TFile("/phenix/hp/data10/mjuszkie/sharkfins/sharkfin_projection_run7.root"); 
  //TFile *f1=new TFile("sharkfin_projections_r7mapfix_notmiss_sum.root");
  //TFile *f1=new TFile("/phenix/hp/data10/manguyen/sharkfins/1b_zbinned_pi0_sharkfinproj_rebin_run4AA.root");

  //TFile *f1=new TFile("/phenix/hp/data10/manguyen/sharkfins/1b_zbinned_pi0_sharkfinsproj_run7AA_miss.root");

  //TFile *f2=new TFile("sharkfin_projections_testsingleII_notmiss.root");
  //TFile *f2=new TFile("sharkfin_projections_single1_miss.root");
  //TFile *f2=new TFile("/phenix/hp/data10/mjuszkie/sharkfins/sharkfin_projection_run7.root");  
  TFile *f2=new TFile("/phenix/hp/data10/mjuszkie/sharkfins/orig_sharkfin_projections_run7_miss.root"); 
  //TFile *f2=new TFile("sharkfin_projections_testmattcode_miss.root");

  TH1F *SIM[33];
  TH1F *DAT[33], *SUB[33];
  //cout << "got the files" <<endl;
  for(int i=0;i<33;i++){
    char simname[100], datname[100], subname[100];
    sprintf(simname,"hshark_large_0_%d",i);
    //sprintf(simname,"hitdiag_sec%d",i);
    //sprintf(datname,"MODMAP%d",i);
    //sprintf(subname,"MODMAPSUB%d",i);
    sprintf(datname,"hshark_large_0_%d",i);
    //sprintf(datname,"piC0_MODMAP%d",i);
    //sprintf(subname,"piC0_MODMAPSUB%d",i);
    SIM[i]=(TH1F*)f1->Get(simname);
    DAT[i]=(TH1F*)f2->Get(datname);
    //SUB[i]=(TH3F*)f2->Get(subname);
    // cout << "got the " << i << "th histo" <<endl;

  }

  TCanvas *c1=new TCanvas("c1","sim",1);
  c1->Divide(4,8);
  for(int i=0;i<32;i++){
    c1->cd(i+1);
    //SIM[i]->Scale(1.0/SIM[i]->Integral(101,400));
    for(int ib=1; ib<80; ib++){ SIM[i]->SetBinContent(ib,0.0);}
    SIM[i]->Draw();
  }
  
  //TCanvas *c2=new TCanvas("c2","dat",1);
  //c2->Divide(2,4);
  for(int i=0;i<32;i++){
    c1->cd(i+1);
    //c2->GetPad(i+1)->SetLogz();

    //DAT[i]->Divide(SIM[i]);
    //DAT[i]->Scale(1.0/DAT[i]->Integral(101,400));
    //DAT[i]->Divide(SIM[i]);
 
   for(int ib=1; ib<80; ib++){ DAT[i]->SetBinContent(ib,0.0);}

   //DAT[i]->Draw();

    DAT[i]->SetLineColor(2);
    DAT[i]->SetLineStyle(2);

    DAT[i]->Draw("same");
    //DAT[i]->Draw("zcol");

    

  }

  /*
  TCanvas *c3=new TCanvas("c3","sub",1);
  c3->Divide(2,4);
  for(int i=0;i<33;i++){
    c3->cd(i+1);
    //SUB[i]->Draw("zcol");
    SUB[i]->Project3D("yx")->Draw("zcol");

  }
  */
  
}
