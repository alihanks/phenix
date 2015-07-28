{
//void calculate_Rgamma_tag(){


  TH1D *tagtrue; 
  TH1D *tagged;
  TH1D *tagfalse;
  TH1D *trigpt;
  double Rtag[4]={0.0};
  int ptmin[4]={5,7,9,12};
  int ptmax[4]={7,9,12,15};
  double Ntr[4], Ntag[4], Ninc[4];

  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4_b2w.root");
  TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent0/direct/cent0_dir_all.root");


  tagtrue =(TH1D*)fin->Get("TrueTagged_0");
  tagged =(TH1D*)fin->Get("TotalTagged_0");
  tagfalse =(TH1D*)fin->Get("FalseTagged_0");
  trigpt =(TH1D*)fin->Get("TRIGPT");

  tagtrue->Draw();

  for(int i=0; i<4; i++){
    Ntr[i] = tagtrue->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Ntag[i] = tagged->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Ninc[i] = trigpt->Integral(ptmin[i]*10+1,ptmax[i]*10);
  }
  
  //tagtrue->Divide(Ntag); //eff true tag


  int centrality_bin=0;
  TFile *fTad=new TFile("/direct/phenix+u/workarea/mjuszkie/run7_QM08/postQM_taxi82/Rgamma/final_Rgamma_51808.root");
  
  char grname[100];
  sprintf(grname,"gr%d",centrality_bin);
  TGraphErrors *gRgamma=(TGraphErrors*)fTad->Get(grname);
  
  TH1F *hsysrgam;
  char sysname[100];
  sprintf(sysname,"systematic%d",centrality_bin);
  hsysrgam = (TH1F *) fTad->Get(sysname);
  
  double dum, R_fine;
  double  Rgamma;
  double  Rerrcomb;
  double  Rerrcorr;	

  for(int i=0; i<4; i++){
    
    gRgamma->GetPoint(i,dum,R_fine);
    Rgamma=  R_fine;
    Rerrcomb=gRgamma->GetErrorY(i);
    Rerrcorr=hsysrgam->GetBinContent(i+1);	
    
    //cout<<"trigbin = "<<trigpt_bin<<": Rgamma = "<<Rgamma<<", Rgamma_err "<<Rerrcomb<<endl;
       
    
    Rtag[i]=(1.0-Ntag[i]/Ninc[i])/(1.0/Rgamma+Ntr[i]*Ntag[i]/Ninc[i]);

    cout << "Rgamma: " << Rgamma << " Rgamma Tagged: " << Rtag <<endl;
  }


}
