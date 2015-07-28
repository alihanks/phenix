#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "math.h"
#include <iostream>
using std::cout;
using std::endl;


void tag_correction_pi0spect(){

  TH1D *hshark_large_sum_femiss[5];
  TH1D *N_true[5];
  TH1D *N_false[5];
  TH1D *hshark_large_sum_miss[5];
  TH1D *hshark_fe_ratio[5];
  TH1D *hshark_ratio[5];

  TH1D *hshark_fe_ratio[5];

  TH2D *PTpi_PTgam_true;
  TH2D *PTpi_PTgam_false;
  TH2D *PTpi_PTgam_miss;
  TH3D *sp_PTpi_ZEMC_PTgam_miss;
  TH2D *sp_PTpi_PTgam_miss;
  TH1D *N_sp[5];
  TH1D *N_miss[5];


  int iptlo[4]={5,7,9,12};
  int ipthi[4]={7,9,12,15};


  TFile *fshark_fe=new TFile("output_merged/fe/new/all_fenew.root");

  //PTpi_PTgam_reg=(TH3D*)fshark_fe->Get("PTpi_ZEMCpi_PTgam");

  PTpi_PTgam_true=(TH2D*)fshark_fe->Get("PTpi_PTgam_truetag");
  PTpi_PTgam_false=(TH2D*)fshark_fe->Get("PTpi_PTgam_falsetag");
  PTpi_PTgam_miss=(TH2D*)fshark_fe->Get("PTpi_PTgam_MISS_PI0");
  TH1D *realpi0pt=(TH1D*)fshark_fe->Get("realpi0pt");
  //cout << realpi0pt->GetNbinsX() << endl;


  //TFile *fshark_exodus=new TFile("/direct/phenix+scratch/mjuszkie/gammajet.exodus/merged_fins/fixed_all.root");

  //sp_PTpi_ZEMC_PTgam_miss=(TH3D*)fshark_exodus->Get("PTpi_ZEMCpi_PTgam");

  //sp_PTpi_PTgam_miss=(TH2D*)sp_PTpi_ZEMC_PTgam_miss->Project3D("xz");


  /*
  PTpi_PTgam_true->Sumw2();
  PTpi_PTgam_false->Sumw2();
  PTpi_PTgam_miss->Sumw2();
  */

  /*
  for (int iptpi=81; iptpi<400; iptpi++){
    float ntottr=0.0;
    float ntotfl=0.0;
    for(int iptgam=1; iptgam<400; iptgam++){
      ntottr+=PTpi_PTgam_true->GetBinContent(iptpi,iptgam);
      ntotfl+=PTpi_PTgam_false->GetBinContent(iptpi,iptgam);
    }
    
    for(int iptgam=1; iptgam<400; iptgam++){

      float bctr=PTpi_PTgam_true->GetBinContent(iptpi,iptgam);
      float bcfl=PTpi_PTgam_false->GetBinContent(iptpi,iptgam);
   
      if(ntottr>0) PTpi_PTgam_true->SetBinContent(iptpi,iptgam,bctr/ntottr);
      if(ntotfl>0) PTpi_PTgam_false->SetBinContent(iptpi,iptgam,bcfl/ntotfl);
    }
  }
  */

  //test norm of 'true miss'

  for(int iptpi=81; iptpi<PTpi_PTgam_miss->GetNbinsX(); iptpi++){
    float ntottr=0.0;
    float ntotfl=0.0;
    for(int iptgam=1; iptgam<PTpi_PTgam_miss->GetNbinsY(); iptgam++){
      //ntottr+=PTpi_PTgam_miss->GetBinContent(iptpi,iptgam);
      //ntotfl+=PTpi_PTgam_false->GetBinContent(iptpi,iptgam);


      //ntottr+=sp_PTpi_PTgam_miss->GetBinContent(iptpi,iptgam);
    }
    int ireal=0;
    if(iptpi%2) ireal=(iptpi-1)/2+1;
    else ireal=iptpi/2;
    float ntottr=realpi0pt->GetBinContent(ireal);    


    float ntot=ntottr+ntotfl;

    for(int iptgam=1; iptgam<PTpi_PTgam_miss->GetNbinsY(); iptgam++){

      float bctr=PTpi_PTgam_miss->GetBinContent(iptpi,iptgam);
      float bcfl=PTpi_PTgam_false->GetBinContent(iptpi,iptgam);

      //if(ntot>0) PTpi_PTgam_true->SetBinContent(iptpi,iptgam,(bctr+bcfl)/ntot);
      //else PTpi_PTgam_true->SetBinContent(iptpi,iptgam,0.0);

      if(ntot>0) PTpi_PTgam_miss->SetBinContent(iptpi,iptgam,(bctr)/ntottr);
      else PTpi_PTgam_miss->SetBinContent(iptpi,iptgam,0.0);

      //if(ntottr>0) PTpi_PTgam_true->SetBinContent(iptpi,iptgam,bctr/ntottr);
      //if(ntotfl>0) PTpi_PTgam_false->SetBinContent(iptpi,iptgam,bcfl/ntotfl);

    }

  }



  
  for(int ipt=0; ipt<3; ipt++){
    char truename[100];
    sprintf(truename,"true_%d_%d",iptlo[ipt],ipthi[ipt]);

    char falsename[100];
    sprintf(falsename,"false_%d_%d",iptlo[ipt],ipthi[ipt]);

    char missname[100];
    sprintf(missname,"miss_%d_%d",iptlo[ipt],ipthi[ipt]);

    //char spname[100];
    //sprintf(spname,"sp_miss_%d_%d",iptlo[ipt],ipthi[ipt]);
 
 
    N_true[ipt]=(TH1D*)PTpi_PTgam_true->ProjectionX(truename,iptlo[ipt]*20+1,ipthi[ipt]*20);

    //N_sp[ipt]=(TH1D*)sp_PTpi_PTgam_miss->ProjectionX(spname,iptlo[ipt]*20+1,ipthi[ipt]*20);

    
    N_false[ipt]=(TH1D*)PTpi_PTgam_false->ProjectionX(falsename,iptlo[ipt]*20+1,ipthi[ipt]*20);
    N_miss[ipt]=(TH1D*)PTpi_PTgam_miss->ProjectionX(missname,iptlo[ipt]*20+1,ipthi[ipt]*20);
    
    //N_false[ipt]->Add(N_true[ipt]);
    //N_true[ipt]->Divide(N_false[ipt]);

    //N_false[ipt]->Add(N_true[ipt]);
    //N_false[ipt]->Divide(N_true[ipt]);


    if(ipt==0) TCanvas *c1 = new TCanvas("c1","c1",1);
    if(ipt==1) TCanvas *c2 = new TCanvas("c2","c2",1);
    if(ipt==2) TCanvas *c3 = new TCanvas("c3","c3",1);
    if(ipt==3) TCanvas *c4 = new TCanvas("c4","c4",1);
 
    N_false[ipt]->Add(N_miss[ipt]);
    N_miss[ipt]->Divide(N_false[ipt]);

    N_miss[ipt]->Draw();
    
    /*
    N_miss[ipt]->Divide(N_true[ipt]);

    N_miss[ipt]->SetLineColor(4);
    // N_miss[ipt]->SetLineStyle(2);
    N_miss[ipt]->SetLineWidth(2);
    //N_true[ipt]->SetLineWidth(2);
    N_false[ipt]->SetLineColor(2);

    N_true[ipt]->Draw();
    //N_miss[ipt]->Draw("");
    //N_sp[ipt]->Draw("same");
    N_false[ipt]->Draw("same");
    */
  }
  
}


