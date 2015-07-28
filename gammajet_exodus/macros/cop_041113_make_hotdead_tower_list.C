#include <iostream.h>
#include <stdlib.h>
#include <TFile.h>
#include <TText.h>


void make_hotdead_tower_list(){

  //run10

  //  TFile *fin=new TFile("run10_ppg113/combinedsimple_run10AuAu_020_merged_taxi313_total.root");
  //    TFile *fin=new TFile("bingr8/mrg1.root");
  TFile *fin=new TFile("mrg_hotr10_4060/gall.root");


  //run 7 merged
  //TFile *fin=new TFile("/direct/phenix+hp/data10/manguyen/run7AA_4QM08/combinesimpletaxi76/inchistos_22_merged.root");
  //TFile *fin=new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r7qm09_modmaps/r5.root");
  //  TFile *fin=new TFile("/phenix/hp/data10/mjuszkie/run7/ppg113/taxi230/merged/passI_all.root");

  //run4 
  //TFile *fin=new TFile("/phenix/hp/data10/manguyen/run4AA_hpdst_4QM08/all_histos_merged_0_20.root");
  //TFile *fin=new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4_a2w.root");
  //TFile *fin=new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4merged/r4all.root");

  //run 6
  //  TFile *fin=new TFile("/phenix/hp/data10/manguyen/run6/fg_ert/run6pp_fg_histo_inc_ert_taxi70_fidcut_fullhitmap_fewhotremoved_122407_merged.root");

  //run 5
  //   TFile *fin=new TFile("/phenix/hp/data10/manguyen/run5/fg_ert/run5pp_fg_histo_inc_ert_taxi50_fidcut_hitmap2edge_122307_merged.root");

  //Sim
  //TFile *fin=new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/cent0_etas_all.root");


  TH3F *MODMAP[8];
  TH2D *MODMAP_pxy[8];

  TCanvas *c=new TCanvas("c","c",1);
  c->Divide(2,4);

  for(int i=0;i<8;i++){
    
    c->cd(i+1);

    char name[100];
    //sprintf(name,"Cent_0__MODMAP%d",i);
    //sprintf(name,"MODMAP%d",i);
    //sprintf(name,"hitdiag_sec%d",i);
    sprintf(name,"C2_MODMAP%d",i);
    MODMAP[i]=(TH3F*)fin->Get(name);
    MODMAP[i]->GetZaxis()->SetRangeUser(5,7);
    MODMAP_pxy[i]=(TH2D*)MODMAP[i]->Project3D("yx");

    MODMAP_pxy[i]->Draw("zcol");


    //**********I don't think writing directly to ff1 works 
    //I just so like root -b -q 'make_hotdead_tower_list.C' > taxi230.log &
    FILE *ff1;
 
    ff1 =fopen("temp_livetowers_r10_test.dat","w");

    for(int j=0;j<MODMAP_pxy[i]->GetNbinsX();j++){
      for(int k=0;k<MODMAP_pxy[i]->GetNbinsY();k++){
	

	if(i<6){

	  //
	  if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)> 0 && MODMAP_pxy[i]->GetBinContent(j+1,k+1)< 4){ 
	    cout<<i*2592+j+72*k<<" "<<MODMAP_pxy[i]->GetBinContent(j+1,k+1)<<endl;

	    fprintf(ff1,"%d %d \n",i*2592+j+72*k,MODMAP_pxy[i]->GetBinContent(j+1,k+1));
	    }


	  //if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)>200) cout<<i*2592+j+72*k<<" "<<MODMAP_pxy[i]->GetBinContent(j+1,k+1)<<endl;
	  //if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)>100) cout<< "sector" << i << " x " << j << " y " << k<<" "<<MODMAP_pxy[i]->GetBinContent(j+1,k+1)<<endl;

	  //if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)>0) cout<<i*2592+j+72*k<<endl;
	}
	else{
	  if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)>0  && MODMAP_pxy[i]->GetBinContent(j+1,k+1)< 4){
	     cout<<15552+(i-6)*4608+j+96*k<<" "<<MODMAP_pxy[i]->GetBinContent(j+1,k+1)<<endl;
	     fprintf(ff1,"%d %d \n",15552+(i-6)*4608+j+96*k,MODMAP_pxy[i]->GetBinContent(j+1,k+1));
	     }

	  //if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)>100) cout<<15552+(i-6)*4608+j+96*k<<" "<<MODMAP_pxy[i]->GetBinContent(j+1,k+1)<<endl;
	  //if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)>0) cout<<15552+(i-6)*4608+j+96*k<<endl;
	}
      }
    }

  }

  fclose(ff1);
  gStyle->SetPalette(1);
 

}
