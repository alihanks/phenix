#include <fstream>
#include <iostream>
#include <TFile>
#include <TGraphErrors>
#include <TH1D>

void calculate_Rgamma_tag(){

  TGraphErrors *grnew; //change to *grnew[4] for all cents
  TGraphErrors *geff;
  TGraphErrors *grsim;

  TH1F *hsysrgam;
  TH1F *htagsysrgam;

  TH1D *tagtrue; 
  TH1D *tagged;
  TH1D *tagfalse;

  TH1D *datatagged;
  TH1D *trigpt;
  TH1D *inctrigpt;

  TH1D *dirstart;  
  TH1D *dirleft;
  TH1D *decpt;
  TH1D *realdir;
  TH1D *pigam;


  double Rsim[4]={0.0};
  double Rtag[4]={0.0};
  double Rtagsys[4]={0.0};
  double Rtaglo[4]={0.0};
  double Rtaghi[4]={0.0};

  double eff[4]={0.0};
  double efferr[4]={0.0};
  double newx[4]={0.0};
  double newxerr[4]={0.0};
  double Rtagerr[4]={0.0};
  int ptmin[4]={5,7,9,12};
  int ptmax[4]={7,9,12,15};
  double Ntr[4], Ntag[4], Ninc[4], Nfl[4], Nstart[4], Nleft[4], Ndec[4];
  double Nsimleft[4], Nsiminc[4];
  double realNdec[4], realNdir[4];
  
  //TFile *fdata =  new TFile("/phenix/scratch/mjuszkie/run4/fb_hpdst/merged/all_frimorn.root");
  TFile *fdata =  new TFile("/phenix/u/workarea/mjuszkie/run7_QM09/taxi136/rg6.root");
  
  datatagged=(TH1D*)fdata->Get("itagC0_TAGCOUNTER");
  trigpt =(TH1D*)fdata->Get("C0_TRIGPT");
  dirleft =(TH1D*)fdata->Get("itagC0_TRIGPT");

  //I thought that itagC0_TAGCOUNTER + itagC0_TRIGPT = C0_TRIGPT
  //But it actually does not ??? TAGCOUNTER may include 'photons' that
  //would not pass the vetocut!
  
  
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/meandmatt1.root");
  TFile *fin =  new TFile("output_merged/fe/new/all_fenew.root");  
  
  tagtrue =(TH1D*)fin->Get("TrueTagged_0");
  tagged =(TH1D*)fin->Get("TotalTagged_0");
  tagfalse =(TH1D*)fin->Get("FalseTagged_0");
  
  
    simtrigpt =(TH1D*)fin->Get("TRIGPT");
    decpt =(TH1D*)fin->Get("DECPT");
    pigam =(TH1D*)fin->Get("pi0gamma");
    realdir =(TH1D*)fin->Get("realomegapt");
    simdirleft =(TH1D*)fin->Get("DIRPT_0");
    dirstart =(TH1D*)fin->Get("DIRPT");
  
  
  //tagtrue->Draw();
  
  for(int i=0; i<4; i++){
    Ntr[i] = tagtrue->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Ntag[i] = tagged->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Nfl[i] = tagfalse->Integral(ptmin[i]*10+1,ptmax[i]*10);
    
    Nsiminc[i] = simtrigpt->Integral(ptmin[i]*10+1,ptmax[i]*10);
    Nsimleft[i] = simdirleft->Integral(ptmin[i]*10+1,ptmax[i]*10);
   
    Ninc[i] = trigpt->Integral(ptmin[i]*10+1,ptmax[i]*10);
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
  
  
  
  char sysname[100];
  sprintf(sysname,"systematic%d",centrality_bin);
  hsysrgam = (TH1F *) fTad->Get(sysname);
  

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

  for(int i=0; i<4; i++){
    
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

    //Rtag[i]=(Ninc[i])/(Ndec[i]);
    
    Rtag[i]=(Nleft[i]/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*(Ninc[i]-Nleft[i])/Ninc[i]+1.0/Rgamma);

    Rtaghi[i]=(Nleft[i]/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*(Ninc[i]-Nleft[i])/Ninc[i]+1.0/(Rgamma+Rerrcorr));
    
    Rtaglo[i]=(Nleft[i]/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*(Ninc[i]-Nleft[i])/Ninc[i]+1.0/(Rgamma-Rerrcorr));
    
    cout << "sys lo " << Rtaglo[i]-Rtag[i] << " sys hi: " << Rtaghi[i]-Rtag[i] <<endl;

    /*
    if(fabs(Rtaglo[i]-Rtag[i])>fabs(Rtaghi[i]-Rtag[i])){
      Rtagsys[i]=fabs(Rtaglo[i]-Rtag[i]);
    }else{
      Rtagsys[i]=fabs(Rtaghi[i]-Rtag[i]);
    }
    */


    // Rtag[i]=(Nleft[i]/Ninc[i])/(1.0/Rgamma+(Nfl[i]/Ntag[i]-1.0)*(Ninc[i]-Nleft[i])/Ninc[i]);
    cout << "for " << i << " Ntrue/Ninc " << (Nleft[i]/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*(Ninc[i]-Nleft[i])/Ninc[i]+1.0/Rgamma)<<endl;


    eff[i]=Rtag[i]/Rgamma;


    //Rsim[i]=(1.0-Ntag[i]/Ninc[i])/(1.0/Rgamma-Ntr[i]/Ninc[i]);

    //Rtag[i]=(Ninc[i]-Ntag[i])/(Ninc[i]-Nstart[i]-(Ntag[i]-Nfl[i]));

  
    //Rsim[i]=(realNdir[i]+realNdec[i])/(realNdec[i]); //blue

    //Rsim[i]=(Ninc[i])/(Ndec[i]);


    //Rtag[i]=(Ninc[i])/(Ninc[i]-Nstart[i]);

    Rsim[i]=(Nsiminc[i]-Ntag[i])/(Ndec[i]-(Ntag[i]-Nfl[i]));
    //Rtag[i]=(Ninc[i]-Ntag[i])/(Ndec[i]-(Ntag[i]-Nfl[i]));


    //Rtag[i]=(Ninc[i]-Ntr[i])/(Ndec[i]-Ntag[i]);
    //Rtag[i]=(Ninc[i]-Ntag[i])/(Ndec[i]-Ntag[i]);


    //eff[i]=Rtag[i]/Rsim[i];


    cout << "Ntag " << Ntag[i] << " Nfalse " << Nfl[i] <<endl;

    cout << "Ndec " << Ndec[i] << " Ninc-Ndir " << Ninc[i]-Nstart[i] <<endl;
    cout << "Rgamma: " << Rgamma << " Rgamma Tagged: " << Rtag[i] <<endl;

    //Rtagerr[i]=Rerr*(1.0-Ntag[i]/Ninc[i])/(1.0-Ntr[i]*Rgamma/Ninc[i])/(1.0-Ntr[i]*Rgamma/Ninc[i]);

    Rtagerr[i]=Rerr*(Nleft[i]/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*Rgamma*(Ninc[i]-Nleft[i])/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*Rgamma*(Ninc[i]-Nleft[i])/Ninc[i]);

    //Rtag[i]=(Nleft[i]/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*(Ninc[i]-Nleft[i])/Ninc[i]+1.0/Rgamma);


    Rtagsys[i]=Rerrcorr*(Nleft[i]/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*Rgamma*(Ninc[i]-Nleft[i])/Ninc[i])/((Nfl[i]/Ntag[i]-1.0)*Rgamma*(Ninc[i]-Nleft[i])/Ninc[i]);

    //Rtagsys[i]=Rerrcorr*(1.0-Ntag[i]/Ninc[i])/(1.0-Ntr[i]*Rgamma/Ninc[i])/(1.0-Ntr[i]*Rgamma/Ninc[i]);

    efferr[i]=eff[i]*sqrt((Rtagerr[i]/Rtag[i])*(Rtagerr[i]/Rtag[i])+(Rerrcomb/Rgamma)*(Rerrcomb/Rgamma));


    cout << "Rgamma err: " << Rerr << " Rgamma err Tagged: " << Rtagerr[i] <<endl;

    hsysrgam->SetBinContent(i,Rtagsys[i]);
  }

  grtag = new TGraphErrors(4,newx,Rtag,newxerr,Rtagerr);

 grtag->SetLineColor(2);
 grtag->SetMarkerColor(2);
 grtag->SetMarkerStyle(8);
 grtag->Draw("AP");
 gRgamma->Draw("P,same");
 
 grsim = new TGraphErrors(4,newx,Rsim,newxerr,Rtagerr);
 
 grsim->SetLineColor(4);
 grsim->SetMarkerColor(4);
 grsim->SetMarkerStyle(4);
 
 grsim->Draw("P,same");
 
 TCanvas *c2 = new TCanvas("c2","c2",1);
 //geff = new TGraphErrors(4,newx,eff,newxerr,newxerr);
 
 geff = new TGraphErrors(4,newx,eff,newxerr,efferr);
 geff->SetLineColor(2);
 geff->SetMarkerColor(2);
 geff->SetMarkerStyle(8);
 geff->Draw("AP");
 
 const char* outfile = "Rgamma_tag_030209_separate.root";
 
 //grsim->SetName("gr0"); grsim->SetTitle("gr0");
 grtag->SetName("gr0"); grtag->SetTitle("gr0");
 
 TFile fout(outfile,"RECREATE");
 cout << "outfile is opened" <<endl;
 //for(int i = 0; i < 4; i++) {
 //grnew[i]->Draw("P");
 grtag->Write();
 hsysrgam->Write();
 //htagsysrgam->Write();
 
 //  hsysnew[i]->Write();
 //}
 // canvas->Write();
 cout << "wrote graph to file" <<endl;
 //dummy->Write();
 fout.Close();
 cout << "closed"<<endl;
 
 
}
