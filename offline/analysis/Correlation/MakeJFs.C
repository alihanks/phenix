#include "MakeJFs.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

MakeJFs::MakeJFs(int type, int centbin, int trigbin, 
                 int partbin, TH1F *CFinc, 
                 double meanpart, double ntrigbg, 
                 TFile* fin, const string v2input, 
                 int nFits, int useMSMP, 
                 TH1F*& CFflowZYAM, 
                 TH1F*& CFjetZYAM)
{
  //type == 0: inclusive photon
  //type == 1: pi0
  //type == 2: dec
  if(useMSMP<2)  SetV2(v2input);
  else {
    inc_v2[centbin][trigbin] = 0;
    pi0_v2[centbin][trigbin] = 0;
    dec_v2[centbin][trigbin] = 0;
    hadron_v2[centbin][partbin] = 0;
  }
  if(type==0) trigv2 = inc_v2[centbin][trigbin];
  if(type==1) trigv2 = pi0_v2[centbin][trigbin];
  if(type==2) trigv2 = dec_v2[centbin][trigbin];
  partv2 = hadron_v2[centbin][partbin];
  c2 = trigv2 * partv2;
  //cout << "set flow constant c2 = " << c2 << endl;

  int nbins = CFinc->GetNbinsX();
  //cout<<"# of bins in CFinc: "<<nbins<<endl;

  //flowFunc = new TF1("flowFunc","[0]*(1.0+2.0*[1]*cos(2*x))",-0.5*PI,1.5*PI);
  flowFunc = new TF1("flowFunc","[0]*(1.0+2.0*[1]*cos(2*x))",0.0,PI);
  flowFunc->SetNpx(nbins*100);
  flowFunc->SetParameter(0, 1.0);
  flowFunc->SetParameter(1,c2);

  ostringstream name;
  name << "CFflowZYAM_" << useMSMP << "_" << centbin << "_" << trigbin << "_" << partbin;
  CFflowZYAM = new TH1F(name.str().c_str(),name.str().c_str(),nbins,0.0,PI);
  CFflowZYAM->Sumw2();

  for( int ibin = 1; ibin < nbins + 1; ibin++){
    //cout<<"ibin = "<<ibin<<endl;
    double binwidth = CFflowZYAM->GetBinWidth(ibin);
    //cout<<"binwidth = "<<binwidth<<endl;
    double x_min = CFflowZYAM->GetBinCenter(ibin)-0.5*binwidth;
    double x_max = CFflowZYAM->GetBinCenter(ibin)+0.5*binwidth;
    //cout<<"x_min = "<<x_min<<"; x_max = "<<x_max<<endl;

    if(CFinc->GetBinContent(ibin) !=0 && CFinc->GetBinError(ibin) !=0){
      CFflowZYAM->SetBinContent(ibin,flowFunc->Integral(x_min,x_max)/binwidth);
      CFflowZYAM->SetBinError(ibin, 0.0); 
    }
  }
  //CFflowZYAM->Sumw2();
  //cout<<"about to subtract bg..."<<endl;
  if(useMSMP==0){
    cout<<"using ZYAM method.........."<<endl;
    //********************************************
    //        Fit the CF                         *
    //********************************************
    TF1 *cfFunc = new TF1("cfFunc","[1]*(1/([2]*sqrt(2*PI)))*exp(-pow(x,2)/(2*pow([2],2)))+[1]*(1/([2]*sqrt(2*PI)))*exp(-pow(x-2*PI,2)/(2*pow([2],2)))+[3]*(1/([5]*sqrt(2*PI)))*exp(-pow((x-PI-[4]),2)/(2*pow([5],2)))+[3]*(1/([5]*sqrt(2*PI)))*exp(-pow((x-PI+[4]),2)/(2*pow([5],2)))+[3]*(1/([5]*sqrt(2*PI)))*exp(-pow((x+PI-[4]),2)/(2*pow([5],2)))+[3]*(1/([5]*sqrt(2*PI)))*exp(-pow((x+PI+[4]),2)/(2*pow([5],2)))+[0]*(1.0+2.0*[6]*cos(2.0*x))+[7]*(1/([8]*sqrt(2*PI)))*exp(-pow(x-PI,2)/(2*pow([8],2)))+[7]*(1/([8]*sqrt(2*PI)))*exp(-pow(x+PI,2)/(2*pow([8],2)))",-0.50*PI,1.50*PI);
    
    cfFunc->SetNpx(1000);
    cfFunc->SetLineWidth(1);
    
    TRandom3 *rand = new TRandom3(time(NULL));
    
    double iniPar0 = 0.0;
    double iniPar1 = 0.0;
    double iniPar2 = 0.0;
    double iniPar3 = 0.0;
    double iniPar4 = 0.0;
    double iniPar5 = 0.0;
    double iniPar6 = 0.0;
    double iniPar7 = 0.0;
    double iniPar8 = 0.0;
    
    iniPar0 = 0.5;
    iniPar1 = 0.5;
    iniPar2 = 0.20;
    iniPar3 = 1.0;
    iniPar4 = 1.1505;
    iniPar5 = 0.437188;
    iniPar6 = 0.00648251;
    iniPar7 = 0.5;
    iniPar8 = 1.20;
    
    if(CFinc->GetBinCenter(1) < 0.0)
    {
      iniPar0 = CFinc->GetBinContent(1+2*nbins/4.0);
      iniPar1 = CFinc->GetBinContent(1+1*nbins/4.0) - CFinc->GetBinContent(1+2*nbins/4.0);
      iniPar7 = CFinc->GetBinContent(1+3*nbins/4.0) - CFinc->GetBinContent(1+2*nbins/4.0);
    }
    else
    {
      iniPar1 = CFinc->GetBinContent(1) - CFinc->GetBinContent(1+2*nbins/4.0);
      iniPar7 = CFinc->GetBinContent(nbins) - CFinc->GetBinContent(1+2*nbins/4.0);
      iniPar0 = CFinc->GetBinContent(1+2*nbins/4.0);
    }

    float gauss_mean =1.0;
    float gauss_sigma=0.10;

    cfFunc->SetParameter(0, iniPar0*rand->Gaus(gauss_mean, gauss_sigma));
    cfFunc->SetParameter(1, iniPar1*rand->Gaus(gauss_mean, gauss_sigma));
    cfFunc->SetParameter(2, iniPar2*rand->Gaus(gauss_mean, gauss_sigma));
    cfFunc->SetParameter(3, iniPar3*rand->Gaus(gauss_mean, gauss_sigma));
    cfFunc->SetParameter(4, iniPar4*rand->Gaus(gauss_mean, gauss_sigma));
    cfFunc->SetParameter(5, iniPar5*rand->Gaus(gauss_mean, gauss_sigma));
    cfFunc->SetParameter(6, c2);
    cfFunc->SetParameter(7, iniPar7*rand->Gaus(gauss_mean, gauss_sigma));
    cfFunc->SetParameter(8, iniPar7*rand->Gaus(gauss_mean, gauss_sigma));

    cfFunc->SetParLimits(0, -1.0, 10);
    cfFunc->SetParLimits(1, 1e-4, 10);
    cfFunc->SetParLimits(2, 1e-4, 2.0);
    cfFunc->FixParameter(3, 0.0);
    cfFunc->FixParameter(4, 0.0);
    cfFunc->FixParameter(5, 0.0);
    cfFunc->FixParameter(6, c2);
    cfFunc->SetParLimits(7, 1e-4, 10);
    cfFunc->SetParLimits(8, 1e-4, 3.0);

    double chiSq = 9999.0;

    double fitPar[9];
    // double errFitPar[9];
    // double fitPar_ErrV2Up[9];
    // double fitPar_ErrV2Down[9];

    for ( int ipar = 0; ipar < 9; ipar++)
      fitPar[ipar] = -9999.0;

    for ( int ifit = 0; ifit < nFits; ifit++){
      cout << "Fit #" << ifit << " for CFinc" << endl;

      cfFunc->SetParameter(0, iniPar0*rand->Gaus(gauss_mean, gauss_sigma));
      cfFunc->SetParameter(1, iniPar1*rand->Gaus(gauss_mean, gauss_sigma));
      cfFunc->SetParameter(2, iniPar2*rand->Gaus(gauss_mean, gauss_sigma));
      cfFunc->SetParameter(7, iniPar7*rand->Gaus(gauss_mean, gauss_sigma));
      cfFunc->SetParameter(8, iniPar8*rand->Gaus(gauss_mean, gauss_sigma));

      CFinc->Fit(cfFunc,"I");
      if (cfFunc->GetChisquare() < chiSq){
       chiSq = cfFunc->GetChisquare();

       for ( int ipar = 0; ipar < 9; ipar++){
         fitPar[ipar] = cfFunc->GetParameter(ipar);
	       //errFitPar[ipar] = cfFunc->GetParError(ipar);
       }
      }
    }
    //*****************************************************
    //*          ZYAM Normalization                       *
    //*****************************************************
    cfFunc->SetParameter(0, fitPar[0]);
    cfFunc->SetParameter(1, fitPar[1]);
    cfFunc->SetParameter(2, fitPar[2]);
    cfFunc->FixParameter(3, 0.0);
    cfFunc->FixParameter(4, 0.0);
    cfFunc->FixParameter(5, 0.0);
    cfFunc->FixParameter(6, c2);
    cfFunc->SetParameter(7, fitPar[7]);
    cfFunc->SetParameter(8, fitPar[8]);

    //c2
    flowFunc->SetParameter(0, 1.0);
    flowFunc->SetParameter(1, c2);

    norm = GetZYAMScale(cfFunc,flowFunc); 
  }
  else if(useMSMP==1) {
    cout<<"using msmp normalization.........."<<endl;
    TFile* fcent = new TFile("CentDists_18bins_AN169.root");
    EvalXi(type,trigbin,partbin,fcent,fin,gXICORR,gXICORRLARGEBIN);

    cout<<"meanpart = "<<meanpart<<"; ntrigbg = "<<ntrigbg<<endl;
    norm = meanpart/ntrigbg/PI*CFinc->GetBinWidth(1);
    norm*=xiavglarge[centbin];
    cout<<"norm = "<<norm<<endl;
    fcent->Close();
  }
  else if (useMSMP==2) {
    CFinc->SetAxisRange(0.9,1.6,"X");
    int lbin = CFinc->GetMinimumBin()-1;//CFinc->FindBin(1.1);
    int hbin = lbin+2;
    float lphi = CFinc->GetBinCenter(lbin);
    float hphi = CFinc->GetBinCenter(hbin);
    norm = CFinc->Integral(lbin,hbin);
    norm = norm/((double)(hbin-lbin+1));
    cout << "ZYAM norm = " << CFinc->Integral(lbin,hbin) << "/(" << hphi << " - " << lphi << ") = " << norm << endl;
    CFinc->SetAxisRange(0.0,TMath::Pi(),"X");
  }
  else if (useMSMP==3) {
    CFinc->SetAxisRange(0.9,1.6,"X");
    int lbin = CFinc->GetMinimumBin()-1;//CFinc->FindBin(1.1);
    int hbin = lbin+2;
    float lphi = CFinc->GetBinCenter(lbin);
    float hphi = CFinc->GetBinCenter(hbin);
    double norm_err = 0;
    norm = CFinc->IntegralAndError(lbin,hbin,norm_err);
    norm = (norm+norm_err)/((double)(hbin-lbin+1));
    cout << "ZYAM norm+err = " << CFinc->Integral(lbin,hbin) << "/(" << hphi << " - " << lphi << ") = " << norm << endl;
    CFinc->SetAxisRange(0.0,TMath::Pi(),"X");
  }

  CFflowZYAM->Scale(norm);

  CFjetZYAM = new TH1F(*(TH1F*)CFinc);
  //CFjetZYAM->Sumw2();
  CFjetZYAM->Add(CFflowZYAM, -1.0);
  //delete CFflowZYAM;
  //cout << "testing yield: " << CFjetZYAM->GetBinContent(5) << endl;
}

void MakeJFs::InitHistos(TH1F* CFflow, string name)
{
  //CFflow = new TH1F(name.c_str(),name.c_str(),60, -PI/2, 3*PI/2);
  CFflow = new TH1F(name.c_str(),name.c_str(),30, 0.0, PI);
  CFflow->GetXaxis()->SetTitle("#Delta#phi [rad]");
}

void MakeJFs::SetV2(const string v2_inputs)
{
  v2file = new TFile(v2_inputs.c_str());
  for(int i=0; i<4; i++){//centrality
    name.str("");
    name << "gamma_inc_v2_"<<i;
    gr_inc_v2[i] = new TGraphErrors(*(TGraphErrors*)v2file->Get(name.str().c_str()));
    double *inc = gr_inc_v2[i]->GetY();
    double *inc_err = gr_inc_v2[i]->GetEY();
    for(int j=0; j<4; j++){
      inc_v2[i][j] = inc[j];
      inc_v2_err[i][j] = inc_err[j];
    }
    
    name.str("");
    name << "gamma_inc_v2sys_"<<i;
    gr_inc_v2sys[i] = new TGraphErrors(*(TGraphErrors*)v2file->Get(name.str().c_str()));
    double *inc_sys = gr_inc_v2sys[i]->GetEY();
    for(int j=0; j<4; j++) inc_v2_sys[i][j] = inc_sys[j];

      name.str("");
    name << "gamma_dec_v2_" << i;
    gr_dec_v2[i] = new TGraphErrors(*(TGraphErrors*)v2file->Get(name.str().c_str()));
    double *dec = gr_dec_v2[i]->GetY();
    double *dec_err = gr_dec_v2[i]->GetEY();
    for(int j=0; j<4; j++){
      dec_v2[i][j] = dec[j];
      dec_v2_err[i][j] = dec_err[j];
    }

    name.str("");
    name << "gamma_dec_v2sys_" << i;
    gr_dec_v2sys[i] = new TGraphErrors(*(TGraphErrors*)v2file->Get(name.str().c_str()));
    double *dec_sys = gr_dec_v2sys[i]->GetEY();
    for(int j=0; j<4; j++) dec_v2_sys[i][j] = dec_sys[j];

      name.str("");
    name << "pi0_v2_" << i;
    gr_pi0_v2[i] = new TGraphErrors(*(TGraphErrors*)v2file->Get(name.str().c_str()));
    double *pi0 = gr_pi0_v2[i]->GetY();
    double *pi0_err = gr_pi0_v2[i]->GetEY();
    for(int j=0; j<4; j++){
      pi0_v2[i][j] = pi0[j];
      pi0_v2_err[i][j] = pi0_err[j];
    }

    name.str("");
    name << "pi0_v2sys_" << i;
    gr_pi0_v2sys[i] = new TGraphErrors(*(TGraphErrors*)v2file->Get(name.str().c_str()));
    double *pi0_sys = gr_pi0_v2sys[i]->GetEY();
    for(int j=0; j<4; j++) pi0_v2_sys[i][j] = pi0_sys[j];

    name.str("");
    name << "hadron_v2_" << i;
    gr_had_v2[i] = new TGraphErrors(*(TGraphErrors*)v2file->Get(name.str().c_str()));
    double *had = gr_had_v2[i]->GetY();
    double *had_err = gr_had_v2[i]->GetEY();
    for(int j=0; j<5; j++){
      hadron_v2[i][j] = had[j];
      hadron_v2_err[i][j] = had_err[j];
    }
    
    name.str("");
    name << "hadron_v2sys_" << i;
    gr_had_v2sys[i] = new TGraphErrors(*(TGraphErrors*)v2file->Get(name.str().c_str()));
    double *had_sys = gr_had_v2sys[i]->GetEY();
    for(int j=0; j<5; j++) hadron_v2_sys[i][j] = had_sys[j];
  }
}

double MakeJFs::GetZYAMScale(TF1* topFunc, TF1* bottomFunc)
{
  norm = 0.0001 * topFunc->GetMinimumX(0.0,2.0); // Start below
  cout<<"norm = "<<norm<<endl;

  double ZYAMpoint = 1.0;

  double min = 0.0;
  double max = 2.0;
  int n = 2000;

  TF1* subFunc = new TF1("subFunc","([0]*(1.0+2.0*[1]*cos(2*x)))",-0.5*PI,1.5*PI);
  subFunc->FixParameter(0,norm);
  subFunc->FixParameter(1,bottomFunc->GetParameter(1));
  cout<<"c2_bottomFunc = "<<bottomFunc->GetParameter(1)<<endl;;
  cout<<"c2_subFunc = "<<subFunc->GetParameter(1)<<endl;;
  cout<<endl;
  double ZYAMscale = 100000;
  double d = (max-min)/(double)n;
  
  for(int i = 0; i< n; i++){
    double temp_scale = 0;
    if(subFunc->Eval(min+i*d)!=0){
      temp_scale = (topFunc->Eval(min + i * d) / subFunc->Eval(min+i*d));
    }
    if(i==0){
      ZYAMscale = temp_scale;
      cout<<"temp_scale = "<<temp_scale<<endl;   
    }
    else{
      //cout<<"i = "<<i<<", x = "<<min+ i*d<<", scale = "<<temp_scale<<endl;
      if(temp_scale < ZYAMscale){
       ZYAMscale = temp_scale;
       ZYAMpoint = min+i*d;
     }
   }
 }
 cout<<"ZYAMpoint = "<<ZYAMpoint<<", ZYAMscale = "<<ZYAMscale * norm<<endl;
 norm = topFunc->Eval(ZYAMpoint)/bottomFunc->Eval(ZYAMpoint);
 cout<<"norm = "<<norm<<", done"<<endl;
 cout<<endl;
 return norm;

}

void MakeJFs::EvalXi(int type, int trigbin, int partbin, TFile *fcentdist, TFile *finput, TGraphErrors*& XICORR, TGraphErrors*& XICORRLARGEBIN)
{
  name.str("");
  if(type == 0) name << "trigfits_inc_" << trigbin << ".txt";
  if(type == 1 || type == 2) name << "trigfits_pi0_" << trigbin << ".txt";
  cout<<"getting input file: "<<name.str().c_str()<<endl;

  ifstream fTrig(name.str().c_str());
  if(!fTrig.is_open()){
    cout << "File '" << name.str().c_str() << "' does not exist. Abort." <<endl;
    exit(0);
  }

  double tmp1, tmp2, tmp3;
  int j = 0;
  while (fTrig >> tmp1 >> tmp2 >> tmp3){
    tpars[j][0] = tmp1;
    tpars[j][1] = tmp2;
    tpars[j][2] = tmp3;
    
    cout << "read " << j << " " 
    << tpars[j][0] << " "
    << tpars[j][1] << " "
    << tpars[j][2] << " "
    << endl;
    
    j++;
  }

  name.str("");
  name << "partfits_" << partbin << ".txt";
  ifstream fPart(name.str().c_str());
  if(!fPart.is_open()){
    cout << "File '" << name.str().c_str() << "' does not exist. Abort." <<endl;
    exit(0);
  }

  j = 0;
  while (fPart >> tmp1 >> tmp2 >> tmp3){
    ppars[j][0] = tmp1;
    ppars[j][1] = tmp2;
    ppars[j][2] = tmp3;
    
    cout << "read " << j << " " 
    << ppars[j][0] << " "
    << ppars[j][1] << " "
    << ppars[j][2] << " "
    << endl;
    
    j++;
  }
  j = 0;
  cout<<"finish reading in parameters."<<endl;
  //fcentdist = new TFile(fcent.c_str());//"CentDists_18bins_AN169.root"

  ofstream xifileout;
  char txtfilename[500];
  if(type == 0) sprintf(txtfilename,"xi_inc_%d_%d.txt",trigbin,partbin);
  if(type == 1) sprintf(txtfilename,"xi_pi0_%d_%d.txt",trigbin,partbin);
  if(type == 2) sprintf(txtfilename,"xi_dec_%d_%d.txt",trigbin,partbin);
  xifileout.open(txtfilename);
  ofstream xifileoutlarge;
  char txtfilenamelarge[500];
  if(type == 0) sprintf(txtfilenamelarge,"xi_inc_large_%d_%d.txt",trigbin,partbin);
  if(type == 1) sprintf(txtfilenamelarge,"xi_pi0_large_%d_%d.txt",trigbin,partbin);
  if(type == 2) sprintf(txtfilenamelarge,"xi_dec_large_%d_%d.txt",trigbin,partbin);
  xifileoutlarge.open(txtfilenamelarge);

  int SetExpFlag;
  int NpartFlag;
  for(int icent=0; icent<18; icent++){

    for(int ifit=0; ifit<4; ifit++){
      if(ifit%2 == 0) SetExpFlag = 0;
      else SetExpFlag = 1;
      if(ifit<2) NpartFlag = 0;
      else NpartFlag = 1;
      
      name.str("");
      if(NpartFlag) name << "AuAu_Npart_" << 5*icent << "_" << 5*(icent+1);
      else name << "AuAu_Ncoll_" << 5*icent << "_" << 5*(icent+1);

      TH1F* nbdweight = new TH1F(*(TH1F*)fcentdist->Get(name.str().c_str()));
      double nbdinteg = nbdweight->Integral();
      nbdweight->Scale(1/nbdinteg);

      double seedmean = 0;
      double partmean = 0;
      double pairmean = 0;


      int ibinmax = 500;//NbinsX for AuAu_Npart_%d_%d
      if(!NpartFlag) ibinmax*=3;//NbinsX for AuAu_Ncoll_%d_%d

      for(int ibin=1;ibin<=ibinmax;ibin++){
	//double nbdw = nbdweight->GetBinContent(ibin+1);
       double nbdw = nbdweight->GetBinContent(ibin);
       if(nbdw<=0) continue;
       double nseeds = 0;
       double nparts = 0;
       double npnc = (double)ibin;

       if(SetExpFlag == 1){
         if(NpartFlag == 1){
           nseeds = 1.0-TMath::Exp(-1.0*fabs(tpars[3][1])*TMath::Power(npnc,tpars[3][2]));
           nparts = 1.0-TMath::Exp(-1.0*fabs(ppars[3][1])*TMath::Power(npnc,ppars[3][2]));
         }	
         else{
           nseeds = 1.0-TMath::Exp(-1.0*fabs(tpars[1][1])*TMath::Power(npnc,tpars[1][2]));
           nparts = 1.0-TMath::Exp(-1.0*fabs(ppars[1][1])*TMath::Power(npnc,ppars[1][2]));
         }
       }
       else{
         if(NpartFlag == 1){
           nseeds = TMath::ATan(-fabs(tpars[2][1])*TMath::Power(npnc,tpars[2][2]));
           nparts = TMath::ATan(-fabs(ppars[2][1])*TMath::Power(npnc,ppars[2][2]));
         }
         else{
           nseeds = TMath::ATan(-fabs(tpars[0][1])*TMath::Power(npnc,tpars[0][2]));
           nparts = TMath::ATan(-fabs(ppars[0][1])*TMath::Power(npnc,ppars[0][2]));
         }	    
       }
       seedmean += nbdw * nseeds;
       partmean += nbdw * nparts;
       pairmean += nbdw * nseeds * nparts;
     }
     xi[icent][ifit] = pairmean / (seedmean * partmean);
      //cout<<"xi["<<icent<<"]["<<ifit<<"]"<<" = "<<xi[icent][ifit]<<endl;
     delete nbdweight;
   }
 }

 char outfilename[100];
 if(type == 0)  sprintf(outfilename,"xi_inc_%d_%d.root",trigbin,partbin);
 if(type == 1)  sprintf(outfilename,"xi_pi0_%d_%d.root",trigbin,partbin);
 if(type == 2)  sprintf(outfilename,"xi_dec_%d_%d.root",trigbin,partbin);

 double cent[18],big[18], small[18];  

 for(int ifit=0;ifit<18;ifit++){
  cout << ifit*5 << "-" << ifit*5+5 << " ";
  cent[ifit]=2.5+ifit*5.0;
  xiavg[ifit] = (xi[ifit][0] + xi[ifit][1] + xi[ifit][2] +xi[ifit][3]) / 4.0;

    //cout << " val 1 " << xi[ifit][0] << " val 2 " << xi[ifit][1] << " val 3 " << xi[ifit][2] << " val 4 " << xi[ifit][3] << "\t\t avg: " << xiavg[ifit] << endl;

    //find error by finding max deviation between diff't fits
  int largest=0;
  int smallest=0;

  for(int jth=1;jth<4;jth++){    
    if(xi[ifit][jth]>xi[ifit][largest])largest=jth;
    if(xi[ifit][jth]<xi[ifit][smallest])smallest=jth;
  }
  xierr[ifit]=(xi[ifit][largest]-xi[ifit][smallest])/2.0;
    //cout<<" largest "<<largest<<" smallest "<<smallest<<" dev "<<xierr[ifit]<<endl;
  big[ifit]=xi[ifit][largest];
  small[ifit]=xi[ifit][smallest];  

  xifileout << xiavg[ifit]<< " "<<xierr[ifit]<<"\n";  
}

TH2D *ZVTXCENTTR;
if(type == 0) ZVTXCENTTR = new TH2D(*(TH2D*) finput->Get("h2_ptvscent_trig_inc"));
else ZVTXCENTTR = new TH2D(*(TH2D*) finput->Get("h2_ptvscent_trig_pi0"));
TH2D *ZVTXCENTPA = new TH2D(*(TH2D*) finput->Get("h2_ptvscent_part"));

TH1F *TRIGGERCENT = new TH1F(*(TH1F*) ZVTXCENTTR->ProjectionY());
TH1F *PARTNERCENT = new TH1F(*(TH1F*) ZVTXCENTPA->ProjectionY());

TRIGGERCENT->Sumw2();
PARTNERCENT->Sumw2();

double centlarge[4]={0};
for(int i=0;i<4;i++) {xiavglarge[i] = 0.0; xierrlarge[i] = 0.0;}

  for(int icl=0;icl<4;icl++){
    centlarge[icl]=20*icl+10;
    double npairstot=0;
    for(int ics=0;ics<4;ics++){
      double npairs =TRIGGERCENT->Integral(icl*20+2+5*ics,icl*20+2+5*(ics+1));
      
      npairstot+=npairs;
      xiavglarge[icl]+=npairs*xiavg[4*icl+ics];
      xierrlarge[icl]+=npairs*xierr[4*icl+ics];
    }
    xiavglarge[icl]/=npairstot;
    xierrlarge[icl]/=npairstot;
    cout<<"xiavglarge["<<icl<<"] = "<< xiavglarge[icl]<<endl;
    xifileoutlarge << xiavglarge[icl]<< " "<<xierrlarge[icl]<<"\n";
  }
  
  XICORR=new TGraphErrors(18,cent,xiavg,0,xierr);
  TGraph *XIBIG=new TGraph(18,cent,big);
  TGraph *XISMALL=new TGraph(18,cent,small);
  XICORRLARGEBIN=new TGraphErrors(4,centlarge,xiavglarge,0,xierrlarge);
  
  XICORR->SetNameTitle("XICORR","xicorr");
  XICORR->GetXaxis()->SetTitle("Centrality (%)");
  XICORR->GetYaxis()->SetTitle("#xi");
  XIBIG->SetName("XIBIG");
  XIBIG->GetXaxis()->SetTitle("Centrality (%)");
  XIBIG->GetYaxis()->SetTitle("#xi");
  XISMALL->SetName("XISMALL");
  XISMALL->GetXaxis()->SetTitle("Centrality (%)");
  XISMALL->GetYaxis()->SetTitle("#xi");
  XICORRLARGEBIN->SetName("XICORRLARGEBIN");
  XICORRLARGEBIN->GetXaxis()->SetTitle("Centrality (%)");
  XICORRLARGEBIN->GetYaxis()->SetTitle("#xi");
  // XICORR->Draw();
  // XICORRLARGEBIN->Draw();
  xifileout.close();
  cout<<" writing file "<<outfilename<<endl;
  
}
