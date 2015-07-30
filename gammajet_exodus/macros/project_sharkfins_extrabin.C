//#include <TFile.h>
//#include <TTree.h>
//#include <TCanvas.h>

//#include "TH1.h"
//#include "TH2.h"
//#include "TH3.h"
//#include "TF1.h"
//#include "math.h"
//#include <iostream>
//using std::cout;
//using std::endl;

void project_sharkfins_extrabin(const char* input, const char* output, const char* histo_name = "PTpi_ZEMCpi_PTgam"){
  TH2D *PTpi_PTgam[33];
  TH1D *hshark_large[5][33];
  TH1D *hshark_small[7][33];
  TH1D *hshark_alt[7][33];
  TH1D* hshark_large_sum[5];

  //TFile *fshark_exodus=new TFile("/phenix/hp/data10/manguyen/sharkfins/1b_zbinned_pi0_sharkfins_rebin_run6pp_wmiss_pbsc6res_pbglAN647.root");
  //TFile *fshark_exodus=new TFile("/direct/phenix+scratch/mjuszkie/gammajet.exodus/merged_fins/fixed_all.root");
  TFile* fshark_exodus = new TFile(input);

  TH3F *PTpi_ZEMCpi_PTgam=(TH3F*)fshark_exodus->Get(histo_name);
  TH2F *RECONPI_ZEMC_PT=(TH2F*)fshark_exodus->Get("RECONPI_ZEMC_PT");

  if(!PTpi_ZEMCpi_PTgam) cout<<" shark fin not found "<<endl;

  //now take 10 cm projections out to 165 cm
  cout<<" projecting sharkfins "<<endl;

  TFile *fout=new TFile(output,"RECREATE");
  
  float ptpi_integral_0=0;
  for(int izemc=0;izemc<PTpi_ZEMCpi_PTgam->GetNbinsY();izemc++){
    ptpi_integral_0+=RECONPI_ZEMC_PT->GetBinContent(izemc+1,81);
  }
  

  //scale by the reconstruction effiency
  for(int iptpi=80;iptpi<PTpi_ZEMCpi_PTgam->GetNbinsX();iptpi++){
    cout<<" iptpi "<<iptpi<<endl;

    float ptpi_integral=0;
    for(int izemc=0;izemc<PTpi_ZEMCpi_PTgam->GetNbinsY();izemc++){
      ptpi_integral+=RECONPI_ZEMC_PT->GetBinContent(izemc+1,iptpi+1);
    }
    ptpi_integral/=ptpi_integral_0;
    for(int izemc=0;izemc<PTpi_ZEMCpi_PTgam->GetNbinsY();izemc++){
      float nrecon=RECONPI_ZEMC_PT->GetBinContent(izemc+1,iptpi+1)/ptpi_integral;
      
      for(int iptgam=0;iptgam<PTpi_ZEMCpi_PTgam->GetNbinsZ();iptgam++){
        float bc = PTpi_ZEMCpi_PTgam->GetBinContent(iptpi+1,izemc+1,iptgam+1);
        if(nrecon>0)PTpi_ZEMCpi_PTgam->SetBinContent(iptpi+1,izemc+1,iptgam+1,bc/nrecon);
        else{
          // if(bc>0) cout<<" no reconstructed pi0s here "<<" iptpi "<<iptpi<<" izemc "<<izemc<<" iptgam "<<iptgam<<endl;
        }
      }
    }
  }


  for(int izemc=0;izemc<33;izemc++){
    cout<<" izemc "<<izemc<<endl;
    char projname[100];
    sprintf(projname,"PTpi_PTgam_%d",izemc);

    //int zemcbinlo=36+10*izemc;
    //int zemcbinhi=45+10*izemc;
    int zemcbinlo=1+izemc;
    int zemcbinhi=2+izemc;


    PTpi_ZEMCpi_PTgam->GetYaxis()->SetRange(zemcbinlo,zemcbinhi);
    PTpi_PTgam[izemc]=(TH2D*)PTpi_ZEMCpi_PTgam->Project3D("zx");      
    //float shark_norm=PTpi_PTgam[izemc]->Integral();
    //if(shark_norm>0)PTpi_PTgam[izemc]->Scale(1/shark_norm);

    //now flatten shark fins in pi0 pt -- necessary?
    for(int iptpi=0;iptpi<PTpi_PTgam[izemc]->GetNbinsX();iptpi++){
      cout<<" iptpi "<<iptpi<<endl;
      float pi0pt_integral = 0.;
      for(int iptgam=0;iptgam<PTpi_PTgam[izemc]->GetNbinsY();iptgam++){
      	pi0pt_integral+=PTpi_PTgam[izemc]->GetBinContent(iptpi+1,iptgam+1);
      }
      for(int iptgam=0;iptgam<PTpi_PTgam[izemc]->GetNbinsY();iptgam++){
      	float bcorg=PTpi_PTgam[izemc]->GetBinContent(iptpi+1,iptgam+1);
      	if(pi0pt_integral>0)PTpi_PTgam[izemc]->SetBinContent(iptpi+1,iptgam+1,bcorg/pi0pt_integral);
      }
    }

    PTpi_PTgam[izemc]->SetName(projname);
    PTpi_PTgam[izemc]->Write();
    
    for(int idecl=0;idecl<5;idecl++){

      cout<<" idecl "<<idecl<<endl;
      
      // 400 bin overs 20 GeV
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
        decbinlo=241;
        decbinhi=300;
      }
      if(idecl==4){
        decbinlo=101;
        decbinhi=200;
      }
      char sharkname[100];
      sprintf(sharkname,"hshark_large_%d_%d",idecl,izemc);
      hshark_large[idecl][izemc]=(TH1D*)PTpi_PTgam[izemc]->ProjectionX(sharkname,decbinlo,decbinhi);
      hshark_large[idecl][izemc]->Write();
      if(izemc==0){
        char sumname[100];
        sprintf(sumname,"hshark_large_sum_%d",idecl);
        hshark_large_sum[idecl]=(TH1D*)PTpi_PTgam[izemc]->ProjectionX(sumname,decbinlo,decbinhi);
      }else{
        hshark_large_sum[idecl]->Add(hshark_large[idecl][izemc]);
        if(izemc==32){ 
          hshark_large_sum[idecl]->Write();
          cout << "writing sum histo " <<endl;
        }
      }
    }

    ///this loop is broken!  careful
    for(int idecs=0;idecs<7;idecs++){

      // 400 bin overs 20 GeV
      int decbinlo=101+10*idecs;
      int decbinhi=110+10*idecs;

      char sharkname[100];
      sprintf(sharkname,"hshark_small_%d_%d",idecs,izemc);
      hshark_small[idecs][izemc]=(TH1D*)PTpi_PTgam[izemc]->ProjectionX(sharkname,decbinlo,decbinhi);
      hshark_small[idecs][izemc]->Write();
    }

    for(int idecs=0;idecs<7;idecs++){

      // 400 bin overs 20 GeV
      int decbinlo=101+20*idecs;
      int decbinhi=120+20*idecs;

      if(idecs==4){
        decbinlo=181;
        decbinhi=240;
      }
      if(idecs==5){
        decbinlo=241;
        decbinhi=300;
      }
      if(idecs==6){
        decbinlo=301;
        decbinhi=400;
      }

      char sharkname[100];
      sprintf(sharkname,"hshark_alt_%d_%d",idecs,izemc);
      hshark_alt[idecs][izemc]=(TH1D*)PTpi_PTgam[izemc]->ProjectionX(sharkname,decbinlo,decbinhi);
      hshark_alt[idecs][izemc]->Write();
    }


  }

  fout->Close();

}
