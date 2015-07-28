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


void test_sharkfins(){
  TH2D *PTpi_PTgam[33];
  TH1D *hshark_large[5][33];
  TH1D *hshark_large_sum[5];
  TH1D *hshark_small[7][33];
  TH1D *hshark_alt[7][33];
  TH2F *RECONPI_ZEMC_PT;
  TH3F *PTpi_ZEMCpi_PTgam;
  TH3F *RECONtemp;
  TH3F *PTpi_ZEMCpi_PTgam_E;
  TH3F *PTpi_ZEMCpi_PTgam_P;


  int istrue=0;


  TFile *fshark_exodus=new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/fe_sharks/good_ak.root");

  if(istrue){
    PTpi_ZEMCpi_PTgam=(TH3F*)fshark_exodus->Get("PTpi_ZEMCpi_PTgam_TRUE");
    //RECONPI_ZEMC_PT=(TH2F*)fshark_exodus->Get("RECONPI_ZEMC_PT");
    RECONtemp=(TH3F*)fshark_exodus->Get("RECONPI_ZEMC_PTTRUE_PTRECO_real");
  }else{
    PTpi_ZEMCpi_PTgam_P=(TH3F*)fshark_exodus->Get("PTpi_ZEMCpi_PTgam_PI0");
    PTpi_ZEMCpi_PTgam_E=(TH3F*)fshark_exodus->Get("PTpi_ZEMCpi_PTgam_ETA");
    PTpi_ZEMCpi_PTgam=(TH3F*)fshark_exodus->Get("PTpi_ZEMCpi_PTgam_DIR");

    PTpi_ZEMCpi_PTgam->Add(PTpi_ZEMCpi_PTgam_P);
    PTpi_ZEMCpi_PTgam->Add(PTpi_ZEMCpi_PTgam_E);


    //PTpi_ZEMCpi_PTgam=(TH3F*)fshark_exodus->Get("PTpi_ZEMCpi_PTgam_false");
    RECONtemp=(TH3F*)fshark_exodus->Get("RECONPI_ZEMC_PTTRUE_PTRECO_false");
  }

  cout << "got histos from file" << PTpi_ZEMCpi_PTgam->GetEntries() <<endl;

  //TFile *fshark_matt=new TFile("/phenix/hp/data10/manguyen/sharkfins/1b_zbinned_pi0_sharkfins_rebin_run4AA_wmiss_pbsc6res_pbglAN647.root");
  //TFile *fshark_matt=new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/sharkfins/all_singII.root");


  //RECONPI_ZEMC_PT=(TH2F*)fshark_matt->Get("RECONPI_ZEMC_PT");


  RECONPI_ZEMC_PT=(TH2F*)RECONtemp->Project3D("xy"); 
  TCanvas *c1 = new TCanvas("c1","c1",1);
  RECONPI_ZEMC_PT->Draw("colz");

  //PTpi_ZEMCpi_PTgam->GetYaxis()->SetRange(zemcbinlo,zemcbinhi);
  PTpi_PTgam[0]=(TH2D*)PTpi_ZEMCpi_PTgam->Project3D("zx");      

  int idecl=0;

     int decbinlo=101;
      int decbinhi=140;
 

      if(idecl==1){
	decbinlo=141;
	decbinhi=180;    
      }
      if(idecl==2){
	decbinlo=181;
	decbinhi=240;    
      }
      if(idecl==3){
	decbinlo=101;
	decbinhi=200;    
      } 
      if(idecl==4){
	decbinlo=241;
	decbinhi=300;    
      } 
      char sharkname[100];
      sprintf(sharkname,"hshark_large_%d_%d",idecl,0);
    
     hshark_large[idecl][0]=(TH1D*)PTpi_PTgam[0]->ProjectionX(sharkname,decbinlo,decbinhi);

     TCanvas *c2 = new TCanvas("c2","c2",1);
  
     hshark_large[idecl][0]->Draw();


}
