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


void tagfin_correction(){

  TH1D *hshark_large_sum_femiss[5];
  TH1D *hshark_large_sum_fe[5];
  TH1D *hshark_large_sum[5];
  TH1D *hshark_large_sum_miss[5];
  TH1D *hshark_fe_ratio[5];
  TH1D *hshark_ratio[5];


  //TFile *fshark_fe=new TFile("sharkfin_projections_fe_nonmiss_singleII_reco.root");
  // TFile *fshark_fe_miss=new TFile("sharkfin_projections_fe_miss_singIIreco_sum.root");
 
  //TFile *fshark_fe=new TFile("sharkfin_projections_fe_false_goodI.root");
   TFile *fshark_fe=new TFile("sharkfin_projections_fe_true_goodI.root");

   //TFile *fshark_fe=new TFile("sharkfin_projections_fe_all_goodI.root");
   TFile *fshark_fe_miss=new TFile("sharkfin_projections_fe_all_realreco_goodI.root");

   //TFile *fshark=new TFile("sharkfin_projections_singleII_notmiss_sum.root");
   //TFile *fshark_miss=new TFile("sharkfin_projections_singII_miss_sum.root");


   int i=1;
   float norm=1.0;
   hshark_large_sum_femiss[i]=(TH1D*)fshark_fe_miss->Get("hshark_large_0_0");
   norm=hshark_large_sum_femiss[i]->Integral(81,400);
   //hshark_large_sum_femiss[i]->Scale(1.0/norm);
   hshark_large_sum_femiss[i]->SetLineColor(2);

   hshark_large_sum_fe[i]=(TH1D*)fshark_fe->Get("hshark_large_0_0");
   norm=hshark_large_sum_fe[i]->Integral(81,400);
   //hshark_large_sum_fe[i]->Scale(1.0/norm);
   hshark_large_sum_fe[i]->SetLineColor(4);

   
  TCanvas *c1=new TCanvas("c1","c1",1);
   //hshark_fe_ratio[i]->Draw();
   hshark_large_sum_femiss[i]->Draw();
   hshark_large_sum_fe[i]->Draw("same");



 
   hshark_fe_ratio[i]=(TH1D*)fshark_fe_miss->Get("hshark_large_0_0");
   norm=hshark_fe_ratio[i]->Integral(81,400);
   //hshark_fe_ratio[i]->Scale(1.0/norm);
   hshark_fe_ratio[i]->Divide(hshark_large_sum_fe[i]);
   

 
   /*
   hshark_large_sum_miss[i]=(TH1D*)fshark_miss->Get("hshark_large_0_0");
   norm=hshark_large_sum_miss[i]->Integral(81,400);
   hshark_large_sum_miss[i]->Scale(1.0/norm);

   hshark_large_sum[i]=(TH1D*)fshark->Get("hshark_large_0_0");
   norm=hshark_large_sum[i]->Integral(81,400);
   hshark_large_sum[i]->Scale(1.0/norm);

   hshark_ratio[i]=(TH1D*)fshark_miss->Get("hshark_large_0_0");
   norm=hshark_ratio[i]->Integral(81,400);
   hshark_ratio[i]->Scale(1.0/norm);
   hshark_ratio[i]->Divide(hshark_large_sum[i]);
   hshark_ratio[i]->SetLineColor(2);
   hshark_ratio[i]->Draw("same");
   */




   TCanvas *c2=new TCanvas("c2","c2",1);

   //TH1D *hshark_double_ratio = (TH1D*)hshark_fe_ratio[i]->Clone("hshark_double_ratio");
   //hshark_double_ratio->Divide(hshark_ratio[i]);
   //hshark_double_ratio->Draw();


   hshark_fe_ratio[i]->Draw();


}
