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


void make_dec2pi0(){

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


  //TFile *fshark_fe=new TFile("output_merged/fe/new/all_fenew.root");
  TFile *fshark_fe=new TFile("/phenix/hp/data10/mjuszkie/exodus/merged/all.root");

  //PTpi_PTgam_reg=(TH3D*)fshark_fe->Get("PTpi_ZEMCpi_PTgam");

  //PTpi_PTgam_reg=(TH3D*)fshark_fe->Get("PTpi_ZEMCpi_PTgam");

  PTpi_PTgam_true=(TH2D*)fshark_fe->Get("PTpi_PTgam_truetag");
  PTpi_PTgam_falsepi=(TH2D*)fshark_fe->Get("PTpi_PTgam_falsetag");
  PTpi_PTgam_falseeta=(TH2D*)fshark_fe->Get("PTpi_PTgam_falsetag_eta");
  PTpi_PTgam_miss=(TH2D*)fshark_fe->Get("PTpi_PTgam_MISS_PI0");

  //TH1D *realpi0pt=(TH1D*)fshark_fe->Get("realpi0pt");
  //cout << realpi0pt->GetNbinsX() << endl;

  //sp_PTpi_PTgam_miss=(TH2D*)fshark_exodus->Get("PTpi_PTgam_TRUE");
  sp_PTpi_PTgam_miss=(TH2D*)fshark_fe->Get("PTpi_PTgam_TRUE");

  float bcmiss[4]={0.0};
  float bcflpi[4]={0.0};
  float bctrue[4]={0.0};
  float bcfleta[4]={0.0};
  float bctot[4]={0.0};
  
  cout << " got all the histos " <<endl;

  for(int iptgam=101; iptgam<301; iptgam++){
    int i=0; 
    for (int ibin=0; ibin<4; ibin++){
      if(iptgam>iptlo[ibin]*20 && iptgam<(ipthi[ibin]*20+1)) i=ibin;
      cout << "with in " << iptlo[ibin]*20 << " & " <<(ipthi[ibin]*20+1) <<endl;
      cout << "iptgam " << iptgam << " in bin " << i <<endl;
    }
    cout << "iptgam " << iptgam << " in bin " << i <<endl;

     cout << "i " << i << "  bctot " << bctot[i] <<endl; 
    for(int iptpi=1; iptpi<PTpi_PTgam_truetag->GetNbinsX(); iptpi++){
     bcmiss[i]+=PTpi_PTgam_miss->GetBinContent(iptpi,iptgam);
     bcflpi[i]+=PTpi_PTgam_falsepi->GetBinContent(iptpi,iptgam);
     bctrue[i]+=PTpi_PTgam_true->GetBinContent(iptpi,iptgam);
     bcfleta[i]+=PTpi_PTgam_falseeta->GetBinContent(iptpi,iptgam);
     bctot[i]+=sp_PTpi_PTgam_miss->GetBinContent(iptpi,iptgam);
    }



  }


  cout << "calculating dec2pi0..." <<endl;

  float dec2pi0[4]; //={1.23};
  float sep_eff=0;

  for(int i=0; i<4; i++){
    int trigpt_bin_mod=i;
    if(trigpt_bin_mod==0)sep_eff = 0.998;
    if(trigpt_bin_mod==1)sep_eff = 0.992;
    if(trigpt_bin_mod==2)sep_eff = 0.965;
    if(trigpt_bin_mod==3)sep_eff = 0.900;
    
    float Ntag=bctrue[i]+bcflpi[i]+bcfleta[i];
    float Npi=bctot[i];
    float Npitag=bctrue[i]+bcflpi[i];

    dec2pi0[i]=(1.23-Ntag/Npi-(1-sep_eff))/(1-Npitag/Npi-(1-sep_eff));

    cout << "ptbin " << i << "  dec2pi0: " << dec2pi0[i] << "  Ntag: " << Ntag << " Npi "<< Npi << " Npitag "<< Npitag <<endl;

  }

}
