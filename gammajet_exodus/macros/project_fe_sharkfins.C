#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "math.h"
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;
using namespace std;


void project_fe_sharkfins(char * filn1, char * out1){


  char line[500];
  //fin.open("pi0spectra.dat");
  char c[20];

  double A[9][30][2]={0.0}; 
  double B[4][30]={0.0}, Be[4][30]={0.0};
  double x[9][30]={0.0};
  float pt,value, err1, err2,staterr1,staterr2, syserr1,syserr2;
  float staterr,syserr,pt_uncorr_err,pt_corr_err;
  double comb_value, comb_err;
  int cent_low=0, cent_high=0;
  int cent_index = -1;
  int pt_index = 0;
  int point[9] = {0};
  float wdndpt[400]={0.0};

  int dounweight=0;

  //if(dounweight)

  TGraphErrors *gr[4];

  if(dounweight){
    ifstream fin ("pi0spectra.dat",ifstream::in);
    while(fin.getline(line,500))
      {
	if( !(line && strlen(line)) ) continue;
	if( strncmp( line,"//",2) == 0 ) continue;
	if( strncmp( line, "Cent:", 5) == 0 ) {
	  cent_index++;
	  
	  sscanf( line, "%s%i%i", &c, &cent_low, &cent_high);
	  pt_index=0;
	  cout<<"centrality bin "<<cent_index<< "  within  "<<cent_low<<" and "<<cent_high<<endl;
	  continue;
	}
	cent_index=0;
	
	if((sscanf( line, "%f%f%f%f%f%f%f%f%f", &pt, &value, &staterr,&staterr1,&pt_uncorr_err, &err1,&syserr,&pt_corr_err,&err2) !=9)) continue;
	
	//cout <<"filling arrays"<< pt << value <<endl;
	x[cent_index][pt_index] = pt;
	A[cent_index][pt_index][0] = value; 
	B[cent_index][pt_index]=value*pt;//*2*3.14;
	//B[cent_index][pt_index]=value;
	//A[cent_index][pt_index][1] = err1; //used total error
	A[cent_index][pt_index][1] = staterr; //now using stat error only
	Be[cent_index][pt_index]=staterr;
	
	//if( ispi0==0) Be[cent_index][pt_index]=err1; //total error
	//else Be[cent_index][pt_index]=staterr;
	pt_index++; 
      }
    
    int cent=0;
    
    point[cent]=23;
    //TGraphErrors *gr[4];
    double xe[5][30] = {0.0};
    gr[cent] = new TGraphErrors(point[cent],x[cent],B[cent],xe[cent],Be[cent]);
    gr[cent]->SetName("gr0"); gr[cent]->SetTitle("gr0");
    gr[cent]->SetMarkerColor(2);
    gr[cent]->SetMarkerStyle(8);
    
    //TF1 *fa = new TF1("fa","[0]/pow(x+[2],[1])",2,15); 
    TF1 *fa = new TF1("fa","([0]*x^(-1*[1])*(1-1/(1+exp((x-3.75)/0.1)))+[2]/((1+x/[3])^[4])*(1/(1+exp((x-3.75)/0.1))))*x",1,15);
    //TF1 *fa = new TF1("fa","[0]*x*exp(9.99*log(1.219/(1.219+x)))",2,15); 
    //fa->SetParameters(1,8,1);
    //fa->SetParameters(76,8,1187,2.74,15.67);
    //75.9661 8.09399 1187.04 2.74135 15.6747 5.21176 29
    //68.1934 8.15446 867.284 2.42962 14.2742 3.8447 28
    fa->SetParameters(68.1934,8.15446,867.284,2.42962,14.2742);
    
      for(int ipt=60; ipt<400; ipt++){
      wdndpt[ipt]=(fa->Integral(60./20.0,400.0/20.0))/(fa->Integral(ipt/20.0,(ipt+1)/20.0));
      cout << "for " << ipt/20.0 << " to " << (ipt+1)/20.0 << " 1/dndpt: " << wdndpt[ipt] <<endl;
      }
    
    // cout << "integral of fb" << fb->Integral(4,15) <<endl;
    

  }

  TH1D *hshark_large_sum[5];

  TH2F *PTpi_PTgam;
  TH1D *hshark_large[5][33];
  TH1D *hshark_small[7][33];
  TH1D *hshark_alt[7][33];
  TH3F *PTpi_ZEMCpi_PTgam;
  TH2F *RECONPI_ZEMC_PT;

  

  //  TFile *fshark_exodus=new TFile("/phenix/scratch/mjuszkie/test3K.root");
  TFile *fshark_exodus=new TFile(filn1);
  PTpi_ZEMCpi_PTgam=(TH3F*)fshark_exodus->Get("PTpi_ZEMCpi_PTgam");
  PTpi_PTgam=(TH2F*)fshark_exodus->Get("PTpi_PTgam_TRUE");
  //RECONPI_ZEMC_PT=(TH2F*)fshark_exodus->Get("RECONPI_ZEMC_PT");

  //PTpi_PTgam->ProjectionX("pxtest",101,140)->Draw();
  //PTpi_PTgam=(TH2F*)PTpi_ZEMCpi_PTgam->Project3D("xz");      


  if(!PTpi_ZEMCpi_PTgam) cout<<" shark fin not found "<<endl;
  
  //now take 10 cm projections out to 165 cm
  cout<<" projecting sharkfins "<<endl;
  
  //  TFile *fout=new TFile("sharkfin_projections_test3K_misspi0.root","RECREATE");
  TFile *fout=new TFile(out1,"RECREATE");
  
  /*
  float ptpi_integral_0=0;
  for(int izemc=0;izemc<PTpi_PTgam->GetNbinsY();izemc++){
    ptpi_integral_0+=RECONPI_ZEMC_PT->GetBinContent(izemc+1,101); //change 81 to 101
  }
  */

  float nrecon=0.0;
  //scale by the reconstruction effiency
  for(int iptpi=80;iptpi<PTpi_PTgam->GetNbinsX();iptpi++){
    cout<<" iptpi "<<iptpi<<endl;

    float ptpi_integral=0;
    for(int iptgam=0;iptgam<PTpi_PTgam->GetNbinsY();iptgam++){
      ptpi_integral+=PTpi_PTgam->GetBinContent(iptpi+1,iptgam+1);
    }
  
    //ptpi_integral/=ptpi_integral_0;
    for(int iptgam=0;iptgam<PTpi_PTgam->GetNbinsY();iptgam++){

      if(ptpi_integral>0) nrecon=ptpi_integral;
      else {if(iptpi<200)cout << "nreco=0 "<< iptgam <<endl; nrecon=0;}
      
      //nrecon=1.0;

      //for(int iptgam=0;iptgam<PTpi_ZEMCpi_PTgam->GetNbinsZ();iptgam++){
	float bc = PTpi_PTgam->GetBinContent(iptpi+1,iptgam+1);
	if(nrecon>0)PTpi_ZEMCpi_PTgam->SetBinContent(iptpi+1,izemc+1,iptgam+1,bc/nrecon*wdndpt[iptpi]);
	//if(nrecon>0)PTpi_PTgam->SetBinContent(iptpi+1,iptgam+1,bc/nrecon);
	else{
	  if(bc>0) cout<<" no reconstructed pi0s here "<<" iptpi "<<iptpi<<" iptgam "<<iptgam<<endl;
	}
	
    }
  }

  
  int izemc=0;
  for(izemc=0;izemc<33;izemc++){ //33
    cout<<" izemc "<<izemc<<endl;
    char projname[100];
    sprintf(projname,"PTpi_PTgam_%d",izemc);

    //int zemcbinlo=36+10*izemc;
    //int zemcbinhi=45+10*izemc;
    int zemcbinlo=1+izemc;
    int zemcbinhi=2+izemc;
    //int zemcbinlo=1;
    //int zemcbinhi=35;
    
    PTpi_ZEMCpi_PTgam->GetYaxis()->SetRange(zemcbinlo,zemcbinhi);
    PTpi_PTgam[izemc]=(TH2D*)PTpi_ZEMCpi_PTgam->Project3D("zx");      
    //float shark_norm=PTpi_PTgam[izemc]->Integral();
    //if(shark_norm>0)PTpi_PTgam[izemc]->Scale(1/shark_norm);
    
    //now flatten shark fins in pi0 pt -- necessary?
    /*
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
    */

    PTpi_PTgam[izemc]->SetName(projname);
    PTpi_PTgam[izemc]->Write();

    
    for(int idecl=0;idecl<5;idecl++){

      cout<<" idecl "<<idecl<<endl;
      //izemc=0;
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
      /*for 5-10
      if(idecl==3){
	decbinlo=101;
	decbinhi=200;    
      } 
      */ //Test full range
      if(idecl==3){
	decbinlo=1;
	decbinhi=400;    
      } 

      if(idecl==4){
	decbinlo=241;
	decbinhi=300;    
      } 
      char sharkname[100];
      sprintf(sharkname,"hshark_large_%d_%d",idecl,izemc);
      hshark_large[idecl][izemc]=(TH1D*)PTpi_PTgam->ProjectionX(sharkname,decbinlo,decbinhi);

      //if(idecl==0) hshark_large[idecl][izemc]->Draw();

      
      if(izemc==0){
	char sumname[100];
	sprintf(sumname,"hshark_large_sum_%d",idecl);
	hshark_large_sum[idecl]=(TH1D*)PTpi_PTgam->ProjectionX(sumname,decbinlo,decbinhi);
      }else{
	hshark_large_sum[idecl]->Add(hshark_large[idecl][izemc]);
	if(izemc==32){ 
	  hshark_large_sum[idecl]->Write();
	  cout << "writing sum histo " <<endl;
	}
      }
      
      hshark_large[idecl][izemc]->Write();
      
    }
    
    ///this loop is broken!  careful
    for(int idecs=0;idecs<7;idecs++){
      
      // 400 bin overs 20 GeV
      int decbinlo=101+10*idecs;
      int decbinhi=110+10*idecs;
      
      char sharkname[100];
      sprintf(sharkname,"hshark_small_%d_%d",idecs,izemc);
      hshark_small[idecs][izemc]=(TH1D*)PTpi_PTgam->ProjectionX(sharkname,decbinlo,decbinhi);
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
      hshark_alt[idecs][izemc]=(TH1D*)PTpi_PTgam->ProjectionX(sharkname,decbinlo,decbinhi);
      hshark_alt[idecs][izemc]->Write();
    }
    
    


  }
  
  fout->Close();

}
