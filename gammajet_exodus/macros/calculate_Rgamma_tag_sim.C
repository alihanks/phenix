{
//void calculate_Rgamma_tag(){


  TGraphErrors *grnew; //change to *grnew[4] for all cents
  TGraphErrors *geff;
  TGraphErrors *grsim;

  TH1D *tagtrue; 
  TH1D *tagged;
  TH1D *tagfalse;
  TH1D *trigpt;
  TH1D *dirstart;  
  TH1D *dirleft;
  TH1D *decpt;

  TH1D *realdir;
  TH1D *pigam;


  double Rsim[4]={0.0};
  double Rtag[4]={0.0};
  double eff[4]={0.0};
  double newx[4]={0.0};
  double newxerr[4]={0.0};
  double Rtagerr[4]={0.0};
  double Rsimerr[4]={0.0};
  int ptmin[4]={5,7,9,12};
  int ptmax[4]={7,9,12,15};
  double Ntr[4], Ntag[4], Ninc[4], Nfl[4], Nstart[4], Nleft[4], Ndec[4];
  double realNdec[4], realNdir[4];

  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4_b2w.root");
  //TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent0/direct/cent0_dir_all.root");
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/pwg1.root");
//TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_newrange.root");
//TFile *fin =  new TFile("/direct/phenix+u/workarea/manguyen/AllCode/awayside/wrk/megan_MC/cent0_newrange_merge.root");
//TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/cent0_etas_all.root");
//TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_hieff.root");
//TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/alltest.root");
//TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_cuts.root");
//TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/new/all_posresmymap.root");
//TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/new/all_goodset2.root");

//TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/meandmatt1.root");

TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/fe_sharks/cent0_bigIII.root");

  tagtrue =(TH1D*)fin->Get("TrueTagged_0");
  tagged =(TH1D*)fin->Get("TotalTagged_0");
  tagfalse =(TH1D*)fin->Get("FalseTagged_0");
  trigpt =(TH1D*)fin->Get("TRIGPT");
  decpt =(TH1D*)fin->Get("DECPT");

  pigam =(TH1D*)fin->Get("pi0gamma");
  realdir =(TH1D*)fin->Get("realomegapt");

  dirleft =(TH1D*)fin->Get("DIRPT_0");
  dirstart =(TH1D*)fin->Get("DIRPT");


//tagtrue->Draw();

  for(int i=0; i<4; i++){
    Ntr[i] = tagtrue->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Ntag[i] = tagged->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Ninc[i] = trigpt->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Nfl[i] = tagfalse->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Nstart[i] = dirstart->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Nleft[i] = dirleft->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Ndec[i] = decpt->Integral(ptmin[i]*10+1,ptmax[i]*10);

    realNdec[i] = pigam->Integral(ptmin[i]*10+1,ptmax[i]*10);
    realNdir[i] = realdir->Integral(ptmin[i]*10+1,ptmax[i]*10);



    cout <<"xxxxxxxx Nfl" << Nfl[i]<<endl;

  }
 
  //tagtrue->Divide(Ntag); //eff true tag


  int centrality_bin=0;
  TFile *fTad=new TFile("/direct/phenix+u/workarea/mjuszkie/run7_QM08/postQM_taxi82/Rgamma/final_Rgamma_51808.root");
  
  char grname[100];
  sprintf(grname,"gr%d",centrality_bin);
  TGraphErrors *gRgamma=(TGraphErrors*)fTad->Get(grname);
  
gRgamma->SetLineColor(1);
gRgamma->SetMarkerColor(1);
gRgamma->SetMarkerStyle(8);

  TH1F *hsysrgam;
  char sysname[100];
  sprintf(sysname,"systematic%d",centrality_bin);
  hsysrgam = (TH1F *) fTad->Get(sysname);
  
  double dum, R_fine;
  double  Rgamma;
  double  Rerrcomb;
  double  Rerrcorr;	
  double  Rerr;

  for(int i=0; i<3; i++){
    
    gRgamma->GetPoint(i,dum,R_fine);
    newx[i]=dum;
    Rgamma=  R_fine;
    Rerrcomb=gRgamma->GetErrorY(i);
    Rerrcorr=hsysrgam->GetBinContent(i+1);	
    //Rerr=sqrt(Rerrcomb*Rerrcomb+Rerrcorr*Rerrcorr);
    Rerr=Rerrcomb;

    //cout<<"trigbin = "<<trigpt_bin<<": Rgamma = "<<Rgamma<<", Rgamma_err "<<Rerrcomb<<endl;
       
    cout << "True Tagging eff= " << Ntr[i]/Ntag[i] <<endl;    
    cout << "False Tagging eff= " << Nfl[i]/Ntag[i] <<endl;    

    //eff[i]=1-Nfl[i]/Ntag[i];
    /*
    if(Nstart[i]==0)eff[i]=0.0;
    else eff[i]=1-Nleft[i]/Nstart[i];
    */

    //Rsim[i]=(1.0-Ntag[i]/Ninc[i])/(1.0/Rgamma-Ntr[i]/Ninc[i]);

    //Rtag[i]=(Ninc[i]-Ntag[i])/(Ninc[i]-Nstart[i]-(Ntag[i]-Nfl[i]));

  
    //Rsim[i]=(realNdir[i]+realNdec[i])/(realNdec[i]); //blue

    Rsim[i]=(Ninc[i])/(Ndec[i]);
    //Rtag[i]=(Ninc[i])/(Ndec[i]);
    //statistical error:
    Rsimerr[i]=sqrt((1.0/Ninc[i])+(1.0/Ndec[i]))*Rsim[i];


    //Rtag[i]=(Ninc[i])/(Ninc[i]-Nstart[i]);

    Rtag[i]=(Ninc[i]-Ntag[i])/(Ndec[i]-(Ntag[i]-Nfl[i]));
    //Rtag[i]=(Ninc[i]-Ntag[i])/(Ndec[i]-(Ntag[i]-Nfl[i]));


    //Rtag[i]=(Ninc[i]-Ntr[i])/(Ndec[i]-Ntag[i]);
    //Rtag[i]=(Ninc[i]-Ntag[i])/(Ndec[i]-Ntag[i]);


    eff[i]=Rtag[i]/Rsim[i];



    cout << "Ntag " << Ntag[i] << " Nfalse " << Nfl[i] <<endl;

    cout << "Ndec " << Ndec[i] << " Ninc-Ndir " << Ninc[i]-Nstart[i] <<endl;
    cout << "Rgamma: " << Rgamma << " Rgamma Tagged: " << Rtag[i] <<endl;

    Rtagerr[i]=Rerr*(1.0-Ntag[i]/Ninc[i])/(1.0-Ntr[i]*Rgamma/Ninc[i])/(1.0-Ntr[i]*Rgamma/Ninc[i]);

    cout << "Rgamma err: " << Rerr << " Rgamma err Tagged: " << Rtagerr[i] <<endl;

  }

 grtag = new TGraphErrors(4,newx,Rtag,newxerr,Rtagerr);

grtag->SetLineColor(2);
grtag->SetMarkerColor(2);
grtag->SetMarkerStyle(8);
 grtag->Draw("AP");
  gRgamma->Draw("P,same");

 grsim = new TGraphErrors(4,newx,Rsim,newxerr,Rsimerr);

grsim->SetLineColor(4);
grsim->SetMarkerColor(4);
grsim->SetMarkerStyle(4);

  grsim->Draw("P,same");




TCanvas *c2 = new TCanvas("c2","c2",1);
 geff = new TGraphErrors(4,newx,eff,newxerr,newxerr);

 geff = new TGraphErrors(4,newx,eff,newxerr,newxerr);
 geff->SetLineColor(2);
 geff->SetMarkerColor(2);
 geff->SetMarkerStyle(8);
geff->Draw("AP");

/*
const char* outfile = "Rgamma_tag_021309.root";

  grsim->SetName("gr0"); grsim->SetTitle("gr0");

  TFile fout(outfile,"RECREATE");
  cout << "outfile is opened" <<endl;
//for(int i = 0; i < 4; i++) {
    //grnew[i]->Draw("P");
    grsim->Write();
//  hsysnew[i]->Write();
//}
  // canvas->Write();
  cout << "wrote graph to file" <<endl;
//dummy->Write();
  fout.Close();
  cout << "closed"<<endl;
*/

}
