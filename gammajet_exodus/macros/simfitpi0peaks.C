{

//TFile *fin =  new TFile("/phenix/hp/data10/manguyen/10m_zbinned_pi0_sharkfins_rebin_run6pp_wmiss_pbsc9_pbglAN647.root");
//TFile *fin =  new TFile("/phenix/hp/data10/manguyen/sharkfins/1b_zbinned_pi0_sharkfins_rebin_run5pp_wmiss.root");
//TFile *fin =  new TFile("$data10/sharkfins/1b_zbinned_pi0_sharkfins_rebin_run4AA_wmiss_pbsc6res_pbglAN647.root");
//TFile *fin =  new TFile("$data10/10m_zbinned_pi0_sharkfins_rebin_run4AA_wmiss_pbsc6_pbglNIM.root");
//TFile *fin =  new TFile("$data10/10m_zbinned_pi0_sharkfins_rebin_run4AA_wmiss_pbscNIM_pbglNIM_weight.root");


//TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent1/cent1_posres_half.root");
//TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent0/cent0_newres_all.root");
TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/test_nopos_6_85.root");

//TH2F *INVMASS = (TH2F*)fin->Get("invmassvspt");
TH2F *INVMASS = (TH2F*)fin->Get("INVMASS_7");



TH1F * massvspt = new TH1F("massvspt","",10, 0, 10);
TH1F * widvspt = new TH1F("widvspt","",10, 0, 10);

TF1 * ff2 = new TF1("ff2","gaus(0)",0,1);
ff2->SetParameter(0,1e4);
ff2->SetParameter(1,0.135);
ff2->SetParameter(2,0.010);


//TF1 * ff2=new TF1("ff2","gaus(0)+pol0(3)",0.09,0.21);
//ff2->SetParameters(100.,0.135,0.01,1.0);




TH1D *hmass[6];
TCanvas *c=new TCanvas("c","c",1);
c->Divide(2,3);

for (int i = 0; i < 6; i++)
{

  c->cd(i+1);

  char name[100];
  sprintf(name,"projtemp%d",i);
  
  hmass[i]=(TH1D*)INVMASS->ProjectionY(name,40+i*10,50+i*10);
  /*
  hmass[i]->Reset();

  for(int iy=0;iy<INVMASS->GetNbinsY();iy++){
    float val=0;
    float err=0;
    for(int ix=40+i*10;ix<50+i*10;ix++){
      val+=INVMASS->GetBinContent(ix+1,iy+1);
      err+=INVMASS->GetBinError(ix+1,iy+1)*INVMASS->GetBinError(ix+1,iy+1);
      
    }
    err=sqrt(err);
    hmass[i]->SetBinContent(iy+1,val);
    //hmass[i]->SetBinError(iy+1,err);
  }
  */
  //hmass = INVMASS->ProjectionY("hi",40+i*10,50+i*10);
  hmass[i]->Fit("ff2","","",0.07,0.25);
  TF1 * fgaus = hmass[i]->GetFunction("ff2");
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

  TCanvas * cwid = new TCanvas("cwid","cwid",1);
  widvspt->Draw();

TFile * mf2 = new TFile("pisimwidths_nopos_6_85.root","RECREATE");
massvspt->Write();
widvspt->Write();



}
