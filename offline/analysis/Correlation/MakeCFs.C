#include "MakeCFs.h"
#include "MakeJFs.h"
#include "Correlation.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualPad.h>

#include <string>
#include <sstream>
#include <iostream>

using namespace std;


MakeCFs::MakeCFs(int type, int isfold, int dofold, int isiso, int ispertrigger, const string fin, const string fout)
{
  cout<<"start"<<endl;
  //type == 0: inclusive photon
  //type == 1: pi0
  //type == 2: decay photon
  infile = new TFile(fin.c_str());
  outfile = new TFile(fout.c_str(),"recreate");
  outfile->cd();
  SetTrigPtBinning(type);
  SetPartPtBinning();
  
  TH1D* h1_trigpt_inc[NCENTBIN];
  TH1D* h1_trigpt_pi0[NCENTBIN];
  TH1D* h1_trigpt_dec[NCENTBIN];
  TH1D* h1_trigpt_inc_mix[NCENTBIN];
  TH1D* h1_trigpt_pi0_mix[NCENTBIN];
  TH1D* h1_trigpt_dec_mix[NCENTBIN];

  // name = "h1_part_pt_tot";
  // TH1D* h1_partpt = new TH1D(*(TH1D*)infile->Get(name.c_str()));

  for(int ic=0; ic<NCENTBIN; ic++){
    bin.str("");
    bin << ic;
    name = "h1_trig_pt_inc_c" + bin.str();
    if( isiso ) name = "h1_trig_pt_inc_iso_c" + bin.str();
    h1_trigpt_inc[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
    name = "h1_trig_pt_pi0_c" + bin.str();
    if( isiso ) name = "h1_trig_pt_pi0_iso_c" + bin.str();
    h1_trigpt_pi0[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
    name = "h1_trig_pt_dec_c" + bin.str();
    h1_trigpt_dec[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
    name = "h1_trig_pt_inc_mix_c" + bin.str();
    if( isiso ) name = "h1_trig_pt_inc_iso_mix_c" + bin.str();
    h1_trigpt_inc_mix[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
    name = "h1_trig_pt_pi0_mix_c" + bin.str();
    if( isiso ) name = "h1_trig_pt_pi0_iso_mix_c" + bin.str();
    h1_trigpt_pi0_mix[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
    name = "h1_trig_pt_dec_mix_c" + bin.str();
    h1_trigpt_dec_mix[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
  }

  string dphi_title = ";#Delta#phi[rad]      ";
  if(!isfold){
    for(int ic=0; ic < NCENTBIN; ic++){
      bin.str("");
      bin << ic;
      if(type == 0){
	name = "h3_dphi_c" + bin.str();
	name_mix = "h3_dphi_mix_c" + bin.str();
      }
      if(type == 1){
	name = "h3_dphi_pi0_c" + bin.str();
	name_mix = "h3_dphi_pi0_mix_c" + bin.str();
      }
      
      if(type == 0 || type == 1){
	temp3D = new TH3D(*(TH3D*)infile->Get(name.c_str()));
	temp3D_mix = new TH3D(*(TH3D*)infile->Get(name_mix.c_str()));
	
	for(int ippt = 0; ippt < NTRIGBIN; ippt++){
	  for(int ihpt = 0; ihpt < NPARTBIN; ihpt++){
	    dphi_3d[ic][ippt][ihpt] = new TH3D(*(TH3D*)temp3D);
	    dphi_3d_mix[ic][ippt][ihpt] = new TH3D(*(TH3D*)temp3D_mix);
	    // SetPtRange(dphi_3d[ic][ippt][ihpt], trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1]);
	    // SetPtRange(dphi_3d_mix[ic][ippt][ihpt], trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1]);
	    
	    bin.str("");
	    bin << ic <<"_p"<<ippt<<"_h"<<ihpt;
	    name = "h1_dphi_c" + bin.str();
	    name_mix = "h1_dphi_mix_c" + bin.str();
	    MakeDphiProjection(dphi_3d[ic][ippt][ihpt],dphi_1d[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1],name.c_str());
	    //dphi_1d[ic][ippt][ihpt]->Sumw2();
            dphi_1d[ic][ippt][ihpt]->SetName(name.c_str());
	    SetHisto(dphi_1d[ic][ippt][ihpt],dphi_title,1);
	    
	    MakeDphiProjection(dphi_3d_mix[ic][ippt][ihpt],dphi_1d_mix[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1],name_mix.c_str());
	    //dphi_3d_mix[ic][ippt][ihpt]->Sumw2();
            dphi_1d_mix[ic][ippt][ihpt]->SetName(name_mix.c_str());
	    SetHisto(dphi_1d_mix[ic][ippt][ihpt],dphi_title,2);
	  }
	}
      }
      if(type == 2){
	name_mix = "h3_dphi_mix_c" + bin.str();
	temp3D_mix = new TH3D(*(TH3D*)infile->Get(name_mix.c_str()));
	for(int ippt=0; ippt<4; ippt++){
	  bin.str("");
	  if(ippt<3){
	    bin<<ippt<<"_c"<<ic;
	    name = "h2_dphi_dec_p" + bin.str();
	  }
	  else {
	    bin<<ippt+1<<"_c"<<ic;
	    name = "h2_dphi_dec_p" + bin.str();
	  }
	  temp2D = new TH2D(*(TH2D*)infile->Get(name.c_str()));
	  
	  for(int ihpt=0; ihpt < NPARTBIN; ihpt++){  
	    dphi_3d_mix[ic][ippt][ihpt] = new TH3D(*(TH3D*)temp3D_mix);    
	    //SetPtRange(dphi_3d_mix[ic][ippt][ihpt], trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1]);
	    dphi_2d[ic][ippt][ihpt] = new TH2D(*(TH2D*)temp2D);
	    //dphi_2d[ic][ippt][ihpt]->GetYaxis()->SetRangeUser(part_pt[ihpt],part_pt[ihpt+1]);
	    bin.str("");
	    bin << ic <<"_p"<<ippt<<"_h"<<ihpt;
	    name = "h1_dphi_c" + bin.str();
	    name_mix = "h1_dphi_mix_c" + bin.str();
	    
	    MakeDphiProjection(dphi_3d_mix[ic][ippt][ihpt],dphi_1d_mix[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1],name_mix.c_str());
	    //dphi_1d_mix[ic][ippt][ihpt]->Sumw2();
	    int ybinlo = dphi_2d[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ippt]);
	    int ybinhi = dphi_2d[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ippt+1]);
	    dphi_1d[ic][ippt][ihpt] = new TH1D (*(TH1D*)dphi_2d[ic][ippt][ihpt]->ProjectionX(name.c_str(),ybinlo,ybinhi));
	    //dphi_1d[ic][ippt][ihpt]->Sumw2();
	  }
	}
      }
    }

    for(int ic=0; ic<NCENTBIN; ic++){
      for(int ippt=0; ippt<NTRIGBIN; ippt++){
	can_dphi_name.str("");
	can_dphi_name << "can_dphi_c"<<ic<<"_p"<<ippt;
	can_dphi[ic][ippt] = new TCanvas(can_dphi_name.str().c_str(), can_dphi_name.str().c_str());
	can_dphi[ic][ippt]->Divide(3,2,0.001,0.001);
	can_corr_name.str("");
	can_corr_name << "can_corr_c"<<ic<<"_p"<<ippt;
	can_corr[ic][ippt] = new TCanvas(can_corr_name.str().c_str(), can_corr_name.str().c_str());
	can_corr[ic][ippt]->Divide(3,2,0.001,0.001);
	
	
	for(int ihpt=0; ihpt<NPARTBIN; ihpt++){
	  bin.str("");
	  bin << ic <<"_p"<<ippt<<"_h"<<ihpt;
	  name = "h1_dphi_c" + bin.str();
	  name_mix = "h1_dphi_mix_c" + bin.str();
	  
	  //*****************************************************
	  if(dofold){
	    name_fold = "dphi_fold_c" + bin.str();
	    FoldDphiDist(dphi_1d[ic][ippt][ihpt], fold[ic][ippt][ihpt],name_fold.c_str());
            fold[ic][ippt][ihpt]->SetName(name_fold.c_str());
	    SetHisto(fold[ic][ippt][ihpt],dphi_title,1);
	    name_fold_mix = "dphi_fold_mix_c" + bin.str();
	    FoldDphiDist(dphi_1d_mix[ic][ippt][ihpt], fold_mix[ic][ippt][ihpt],name_fold_mix.c_str());
            fold_mix[ic][ippt][ihpt]->SetName(name_fold_mix.c_str());
	    SetHisto(fold_mix[ic][ippt][ihpt],dphi_title,2);
	    
	    double R_fg = fold[ic][ippt][ihpt]->Integral("width")/PI;
	    fold[ic][ippt][ihpt]->Scale(1/R_fg);
	    double R_bg = fold_mix[ic][ippt][ihpt]->Integral("width")/PI;
	    fold_mix[ic][ippt][ihpt]->Scale(1/R_bg);
	    
	    fold[ic][ippt][ihpt]->Write();
	    fold_mix[ic][ippt][ihpt]->Write();
	    
	    TVirtualPad* pad = can_dphi[ic][ippt]->cd(ihpt+1);
            SetPad(pad);
	    fold[ic][ippt][ihpt]->Draw();
	    fold_mix[ic][ippt][ihpt]->Draw("same");
	    
	    legend_name.str("");
	    legend_name<<trig_pt[ippt]<<"-"<<trig_pt[ippt+1]<<" #times "<<part_pt[ihpt]<<"-"<<part_pt[ihpt+1]<<" GeV/c";
	    TLegend *l1 = new TLegend(0.5,0.7,0.8,0.9,legend_name.str().c_str(),"brNDC");
	    l1->AddEntry(fold[ic][ippt][ihpt],"real","lpf");
	    l1->AddEntry(fold_mix[ic][ippt][ihpt],"mixed","lpf");
	    l1->SetTextSize(0.05);
	    l1->Draw("same");
	    
	    //*****************************************************
	    
	    corr_name.str("");
	    corr_name << "CF_c" << ic << "_p"<<ippt <<"_h"<< ihpt; 
	    corr[ic][ippt][ihpt] = new TH1D(*(TH1D*)fold[ic][ippt][ihpt]);
	    corr[ic][ippt][ihpt]->Divide(fold_mix[ic][ippt][ihpt]);
            corr[ic][ippt][ihpt]->SetName(corr_name.str().c_str());
            SetHisto(corr[ic][ippt][ihpt],dphi_title,1);
	    //corr[ic][ippt][ihpt]->Rebin();
	    
	    double r = corr[ic][ippt][ihpt]->Integral("width")/PI;
	    corr[ic][ippt][ihpt]->Scale(1/r);
	    // if(type == 0) num_trigger = GetNTriggers(h1_trigpt_inc, trig_pt[ippt], trig_pt[ippt+1]);
	    // if(type == 1) num_trigger = GetNTriggers(h1_trigpt_pi0, trig_pt[ippt], trig_pt[ippt+1]);
	    // corr[ic][ippt][ihpt]->Scale(1/num_trigger);
	    corr[ic][ippt][ihpt]->Write();
	    
	    TLatex *la = new TLatex(0.45, 0.75, legend_name.str().c_str());
	    la->SetNDC();
            pad = can_corr[ic][ippt]->cd(ihpt+1);
            SetPad(pad);
	    corr[ic][ippt][ihpt]->Draw();
	    la->Draw("same");
	  }
	  
	  else{
	    corr_name.str("");
	    corr_name << "CF_c" << ic << "_p"<<ippt <<"_h"<< ihpt; 
	    corr[ic][ippt][ihpt] = new TH1D(*(TH1D*)dphi_1d[ic][ippt][ihpt]);
	    corr[ic][ippt][ihpt]->Divide(dphi_1d_mix[ic][ippt][ihpt]);
	    if(type == 0) num_trigger[ic][ippt] = GetNTriggers(h1_trigpt_inc[ic], trig_pt[ippt], trig_pt[ippt+1]);
	    if(type == 1) num_trigger[ic][ippt] = GetNTriggers(h1_trigpt_pi0[ic], trig_pt[ippt], trig_pt[ippt+1]);
	    corr[ic][ippt][ihpt]->Scale(1/num_trigger[ic][ippt]);
            corr[ic][ippt][ihpt]->SetName(corr_name.str().c_str());
	    SetHisto(corr[ic][ippt][ihpt],dphi_title,1);
	    
	    TLatex *la = new TLatex(0.45, 0.75, legend_name.str().c_str());
	    la->SetNDC();
	    TVirtualPad* pad = can_corr[ic][ippt]->cd(ihpt+1);
            SetPad(pad);
	    corr[ic][ippt][ihpt]->Draw();
            corr[ic][ippt][ihpt]->Write();
	    la->Draw("same");
	  }
	}
	can_dphi[ic][ippt]->Write();
	can_corr[ic][ippt]->Write();
      }
    }
  }

  else{
    cout<<"using folded histos."<<endl;
    for(int ic=0; ic < NCENTBIN; ic++){
      bin.str("");
      bin << ic;
      if(type == 0){
        name = "h3_dphi_fold_c" + bin.str();
        if( isiso ) name = "h3_dphi_iso_fold_c" + bin.str();
        name_mix = "h3_dphi_mix_fold_c" + bin.str();
        if( isiso ) name_mix = "h3_dphi_mix_iso_fold_c" + bin.str();
      }
      if(type == 1){
        name = "h3_dphi_pi0_fold_c" + bin.str();
        if( isiso ) name = "h3_dphi_pi0_iso_fold_c" + bin.str();
        name_mix = "h3_dphi_pi0_mix_fold_c" + bin.str();
        if( isiso ) name_mix = "h3_dphi_pi0_mix_iso_fold_c" + bin.str();
      }
      
      if(type == 0 || type == 1){
	temp3D = new TH3D(*(TH3D*)infile->Get(name.c_str()));
	temp3D_mix = new TH3D(*(TH3D*)infile->Get(name_mix.c_str()));
	//cout<<"get temp3D: "<<name.c_str()<<endl;
	for(int ippt = 0; ippt < NTRIGBIN; ippt++){
	  for(int ihpt = 0; ihpt < NPARTBIN; ihpt++){
	    dphi_3d[ic][ippt][ihpt] = new TH3D(*(TH3D*)temp3D);
	    dphi_3d_mix[ic][ippt][ihpt] = new TH3D(*(TH3D*)temp3D_mix);
            
            bin.str("");
            bin << ic <<"_p"<<ippt<<"_h"<<ihpt;
            name = "h1_dphi_c" + bin.str();
            name_mix = "h1_dphi_mix_c" + bin.str();
            MakeDphiProjection(dphi_3d[ic][ippt][ihpt],dphi_1d[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1],name.c_str());
            SetHisto(dphi_1d[ic][ippt][ihpt],dphi_title,1);
            dphi_1d[ic][ippt][ihpt]->SetName(name.c_str());
            
            MakeDphiProjection(dphi_3d_mix[ic][ippt][ihpt],dphi_1d_mix[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1],name_mix.c_str());
            SetHisto(dphi_1d_mix[ic][ippt][ihpt],dphi_title,2);
            dphi_1d_mix[ic][ippt][ihpt]->SetName(name_mix.c_str());
            if(ispertrigger) {
              dphi_1d_mix[ic][ippt][ihpt]->Scale(1/100.0);
              meanpart[ic][ippt][ihpt] = dphi_1d_mix[ic][ippt][ihpt]->Integral();
            }
          }
        }
      }
      if(type == 2){
        cout<<"1";
        for(int ippt=0; ippt<4; ippt++){
          bin.str("");
          if(ippt<3) bin<<ippt<<"_c"<<ic;	   
          else bin<<ippt+1<<"_c"<<ic;
          
          name = "h2_dphi_dec_fold_p" + bin.str();
          name_mix = "h2_dphi_dec_mix_fold_p" + bin.str();
          
          temp2D = new TH2D(*(TH2D*)infile->Get(name.c_str()));
	  temp2D_mix = new TH2D(*(TH2D*)infile->Get(name_mix.c_str()));
	  cout<<"2";
	  for(int ihpt=0; ihpt < NPARTBIN; ihpt++){  
	    dphi_2d[ic][ippt][ihpt] = new TH2D(*(TH2D*)temp2D);
	    dphi_2d_mix[ic][ippt][ihpt] = new TH2D(*(TH2D*)temp2D_mix);
	    bin.str("");
	    bin << ic <<"_p"<<ippt<<"_h"<<ihpt;
	    name = "h1_dphi_c" + bin.str();
	    name_mix = "h1_dphi_mix_c" + bin.str();
	    
	    cout<<"3";
	    int ymin = dphi_2d[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ihpt]);
	    int ymax = dphi_2d[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ihpt+1]);
	    dphi_1d[ic][ippt][ihpt] = new TH1D (*(TH1D*)dphi_2d[ic][ippt][ihpt]->ProjectionX(name.c_str(),ymin,ymax));
	    SetHisto(dphi_1d[ic][ippt][ihpt],dphi_title,1);
	    ymin = dphi_2d_mix[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ihpt]);
	    ymax = dphi_2d_mix[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ihpt+1]);
	    dphi_1d_mix[ic][ippt][ihpt] = new TH1D (*(TH1D*)dphi_2d_mix[ic][ippt][ihpt]->ProjectionX(name_mix.c_str(),ymin,ymax));
	    SetHisto(dphi_1d_mix[ic][ippt][ihpt],dphi_title,2);
	    if(ispertrigger) {
	      dphi_1d_mix[ic][ippt][ihpt]->Scale(1/100.0);
	      meanpart[ic][ippt][ihpt] = dphi_1d_mix[ic][ippt][ihpt]->Integral();
	    }
	    cout<<"4";
	  }
	}
      }
    }
    
    cout<<"5";
    for(int ic=0; ic<NCENTBIN; ic++){
      for(int ippt=0; ippt<NTRIGBIN; ippt++){
	can_dphi_name.str("");
	can_dphi_name << "can_dphi_c"<<ic<<"_p"<<ippt;
	can_dphi[ic][ippt] = new TCanvas(can_dphi_name.str().c_str(), can_dphi_name.str().c_str());
	can_dphi[ic][ippt]->Divide(3,2,0.001,0.001);
	can_corr_name.str("");
	can_corr_name << "can_corr_c"<<ic<<"_p"<<ippt;
	can_corr[ic][ippt] = new TCanvas(can_corr_name.str().c_str(), can_corr_name.str().c_str());
	can_corr[ic][ippt]->Divide(3,2,0.001,0.001);
	can_corr_name.str("");
	can_corr_name << "can_jet_c"<<ic<<"_p"<<ippt;
	can_jet[ic][ippt] = new TCanvas(can_corr_name.str().c_str(), can_corr_name.str().c_str());
	can_jet[ic][ippt]->Divide(3,2,0.001,0.001);
	can_corr_name.str("");
	can_corr_name << "can_abs_c"<<ic<<"_p"<<ippt;
	can_abs[ic][ippt] = new TCanvas(can_corr_name.str().c_str(), can_corr_name.str().c_str());
	can_abs[ic][ippt]->Divide(3,2,0.001,0.001);

	if(type == 0) {
	  num_trigger[ic][ippt] = GetNTriggers(h1_trigpt_inc[ic], trig_pt[ippt], trig_pt[ippt+1]);
	  num_trigger_mix[ic][ippt] = GetNTriggers(h1_trigpt_inc_mix[ic], trig_pt[ippt], trig_pt[ippt+1]);
	}
	if(type == 1) {
	  num_trigger[ic][ippt] = GetNTriggers(h1_trigpt_pi0[ic], trig_pt[ippt], trig_pt[ippt+1]);
	  num_trigger_mix[ic][ippt] = GetNTriggers(h1_trigpt_pi0_mix[ic], trig_pt[ippt], trig_pt[ippt+1]);
	}
	if(type == 2) {
	  if(ippt<3){
	    num_trigger[ic][ippt] = h1_trigpt_dec[ic]->GetBinContent(ippt+1);
	    num_trigger_mix[ic][ippt] = h1_trigpt_dec_mix[ic]->GetBinContent(ippt+1);
	  }
	  else{
	    num_trigger[ic][ippt] = h1_trigpt_dec[ic]->GetBinContent(ippt+2);
	    num_trigger_mix[ic][ippt] = h1_trigpt_dec_mix[ic]->GetBinContent(ippt+2);
	  }
	}

	cout<<"num_trigger = "<<num_trigger[ic][ippt]<<endl;
	cout<<"num_trigger_mix = "<<num_trigger_mix[ic][ippt]<<endl;
      	
	for(int ihpt=0; ihpt<NPARTBIN; ihpt++){

	  double R_bg = dphi_1d_mix[ic][ippt][ihpt]->Integral("width")/PI;
	  dphi_1d_mix[ic][ippt][ihpt]->Scale(1/R_bg);

	  double R_fg = dphi_1d[ic][ippt][ihpt]->Integral("width")/PI;
	  if(!ispertrigger) dphi_1d[ic][ippt][ihpt]->Scale(1/R_fg);
	  
	  TVirtualPad* pad = can_dphi[ic][ippt]->cd(ihpt+1);
          SetPad(pad);
	  dphi_1d[ic][ippt][ihpt]->Draw();
	  dphi_1d_mix[ic][ippt][ihpt]->Draw("same");
	  
	  legend_name.str("");
	  legend_name<<trig_pt[ippt]<<"-"<<trig_pt[ippt+1]<<" #times "<<part_pt[ihpt]<<"-"<<part_pt[ihpt+1]<<" GeV/c";
	  TLegend *l1 = new TLegend(0.5,0.7,0.8,0.9,legend_name.str().c_str(),"brNDC");
	  l1->AddEntry(dphi_1d[ic][ippt][ihpt],"real","lpf");
	  l1->AddEntry(dphi_1d_mix[ic][ippt][ihpt],"mixed","lpf");
	  l1->SetTextSize(0.05);
	  l1->Draw("same");
	  
	  //*****************************************************	  
	  corr_name.str("");
	  corr_name << "CF_c" << ic << "_p"<<ippt <<"_h"<< ihpt; 
	  corr[ic][ippt][ihpt] = new TH1D(*(TH1D*)dphi_1d[ic][ippt][ihpt]);
	  corr[ic][ippt][ihpt]->Divide(dphi_1d_mix[ic][ippt][ihpt]);
	  SetHisto(corr[ic][ippt][ihpt],dphi_title,1);
          corr[ic][ippt][ihpt]->SetName(corr_name.str().c_str());
	  
	  if(ispertrigger) corr[ic][ippt][ihpt]->Scale(1/num_trigger[ic][ippt]);
	  
	  else{
	    double r = corr[ic][ippt][ihpt]->Integral("width")/PI;
	    corr[ic][ippt][ihpt]->Scale(1/r);
	  }
	  
	  TLatex *la = new TLatex(0.45, 0.75, legend_name.str().c_str());
	  la->SetNDC();
	  pad = can_corr[ic][ippt]->cd(ihpt+1);
          SetPad(pad);
	  corr[ic][ippt][ihpt]->Draw();
          corr[ic][ippt][ihpt]->Write();
	  la->Draw("same");
	}
	can_dphi[ic][ippt]->Write();
	can_corr[ic][ippt]->Write();
      }
    }
  }
  if(ispertrigger){
    //make JFs.
    for(int ic=0; ic<4; ic++){
      for(int itrig=0; itrig<4; itrig++){
	for(int ipart=0; ipart<5; ipart++){
	  MakeJFs(type,ic,itrig,ipart,corr[ic][itrig][ipart],meanpart[ic][itrig][ipart],num_trigger_mix[ic][itrig],infile,"v2_inputs.root",5,1,flow[ic][itrig][ipart],jet[ic][itrig][ipart]);
	  bin.str("");
	  bin << ic <<"_p"<<itrig<<"_h"<<ipart;
	  name = "JF_c" + bin.str();
	  //cout<<"name: "<<name.c_str()<<endl;
	  jet[ic][itrig][ipart]->SetName(name.c_str());
	  //jet[ic][itrig][ipart]->Sumw2();
	  //if(type==0) jet[ic][itrig][ipart]->GetYaxis()->SetTitle("Arb.Unit");
	  if(type==0) jet[ic][itrig][ipart]->GetYaxis()->SetTitle("Y_{inc}");
	  if(type==1) jet[ic][itrig][ipart]->GetYaxis()->SetTitle("Y_{#pi^{0}}");
	  if(type==2) jet[ic][itrig][ipart]->GetYaxis()->SetTitle("Y_{dec}");

          double binwidth = jet[ic][itrig][ipart]->GetBinWidth(1);
	  cout<<"jet func binwidth: "<<binwidth<<endl;
	  jet[ic][itrig][ipart]->Scale(1/binwidth);
          SetHisto(jet[ic][itrig][ipart],dphi_title,1);
          TVirtualPad* pad = can_jet[ic][itrig]->cd(ipart+1);
          SetPad(pad);
          jet[ic][itrig][ipart]->Draw();
	  
	  legend_name.str("");
	  legend_name<<trig_pt[itrig]<<"-"<<trig_pt[itrig+1]<<" #times "<<part_pt[ipart]<<"-"<<part_pt[ipart+1]<<" GeV/c";
	  TLatex *la = new TLatex(0.25, 0.75, legend_name.str().c_str());
	  la->SetNDC();
	  la->Draw("same");
	  
	  name = "flow_c" + bin.str();
          pad = can_abs[ic][itrig]->cd(ipart+1);
          SetPad(pad);
	  corr[ic][itrig][ipart]->Draw();
	  flow[ic][itrig][ipart]->SetLineColor(2);
	  flow[ic][itrig][ipart]->Draw("same");
	  la->Draw("same");
	  outfile->cd();
	  //jet[ic][itrig][ipart]->Write();
	  //cout<<"jetfunc written out."<<endl;
	}
	can_abs[ic][itrig]->Write();
	can_jet[ic][itrig]->Write();
      }
    }
  }
  outfile->Close();
  
}

void MakeCFs::SetPtRange(TH3D* h3, double x_pt_min, double x_pt_max, double y_pt_min, double y_pt_max)
{
  h3->GetXaxis()->SetRangeUser(x_pt_min, x_pt_max);
  h3->GetYaxis()->SetRangeUser(y_pt_min, y_pt_max);
}

void MakeCFs::SetTrigPtBinning(int type)
{
  trig_pt[0]=5.;
  trig_pt[1]=7.;
  trig_pt[2]=9.;
  trig_pt[3]=12.;
  trig_pt[4]=15.;
}

void MakeCFs::SetPartPtBinning()
{
  part_pt[0]=0.5;
  part_pt[1]=1.;
  part_pt[2]=2.;
  part_pt[3]=3.;
  part_pt[4]=5.;
  part_pt[5]=7.;
}

void MakeCFs::MakeDphiProjection(TH3D* h3, TH1D*& h1,string hname)
{
  h1 = new TH1D(*(TH1D*)h3->Project3D("z"));
  h1->SetNameTitle(hname.c_str(),hname.c_str());
  // cout<<"h1 name: "<<hname.c_str()<<endl;
}

void MakeCFs::MakeDphiProjection(TH3D* h3, TH1D*& h1,double xmin, double xmax, double ymin, double ymax, string hname)
{
  TH1D* proj_x = (TH1D*)h3->ProjectionX("px");
  TH1D* proj_y = (TH1D*)h3->ProjectionY("py");
  int xbinlo = proj_x->FindBin(xmin);
  int xbinhi = proj_x->FindBin(xmax);
  int ybinlo = proj_y->FindBin(ymin);
  int ybinhi = proj_y->FindBin(ymax);
  // cout<<"xbinlo = "<<xbinlo<<"; xbinhi = "<<xbinhi<<endl;
  // cout<<"ybinlo = "<<ybinlo<<"; ybinhi = "<<ybinhi<<endl;
  h1 = new TH1D(*(TH1D*)h3->Project3D("z"));
  string pz = hname + "_pz";
  h1 = new TH1D(*(TH1D*)h3->ProjectionZ(pz.c_str(),xbinlo,xbinhi-1,ybinlo,ybinhi-1));
}

void MakeCFs::FoldDphiDist(TH1D* h1, TH1D*& h1_fold, string hname_fold)
{
  int Nbins = h1->GetNbinsX();
  h1_fold = new TH1D(hname_fold.c_str(),"",Nbins/2,0.0,PI);
  for(int ibin = 1; ibin <= Nbins/4; ibin++){
    double added_bin_left = h1->GetBinContent(ibin)+h1->GetBinContent(Nbins/2-ibin+1);
    double added_bin_right = h1->GetBinContent(Nbins/2+ibin)+h1->GetBinContent(Nbins-ibin+1);    
    // double lerr1 = h1->GetBinError(ibin);
    // double lerr2 = h1->GetBinError(Nbins/2-ibin+1);
    // double rerr1 = h1->GetBinError(Nbins/2+ibin);
    // double rerr2 = h1->GetBinError(Nbins-ibin+1);
    
    h1_fold->SetBinContent(16-ibin,added_bin_left);
    h1_fold->SetBinContent(ibin+15,added_bin_right);
    // h1_fold->SetBinError(16-ibin,sqrt(lerr1*lerr1+lerr2*lerr2));
    // h1_fold->SetBinError(ibin+15,sqrt(rerr1*rerr1+rerr2*rerr2));
  }
}

void MakeCFs::SetHisto(TH1D* h1, string title, int color)
{
  h1->SetName("");
  h1->SetTitle(title.c_str());
  h1->SetTitleSize(0.05,"X");
  h1->SetTitleOffset(0.95,"X");
  h1->SetLabelSize(0.05,"X");
  h1->SetMarkerColor(color);
  h1->SetLineColor(color);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.8);
}

void MakeCFs::SetPad(TVirtualPad* pad)
{
  pad->SetTopMargin(0.02);
  pad->SetRightMargin(0.02);
}

double MakeCFs::GetNTriggers(TH1D* trigpt, double trigptmin, double trigptmax)
{
  double ntrig = 0;
  for(int i = 0; i < trigpt->GetNbinsX(); i++)
  {
    double bincenter = trigpt->GetBinCenter(i+1);
    if(bincenter >=trigptmin && bincenter <= trigptmax)
    {
      double bincontent = trigpt->GetBinContent(i+1);
      ntrig += bincontent;
    }
  }
  return ntrig;
}

double MakeCFs::GetHadronEff(TH1D* hadron_pt, int ipart)
{
  TF1* fhadeff = new TF1("fhadeff","(x<=3)*(0.396-0.337*exp(-1.50*x))+(x>3&&x<=5)*(0.400+-2.13e-13*exp(5.36*x))+(x>5)*(0.396-0.337*exp(-1.50*x))",0.5,10.0);
  hadron_pt->SetAxisRange(part_pt[ipart],part_pt[ipart+1],"X");
  double meanpt = hadron_pt->GetMean();
  cout<<"mean pt for pt within ["<<part_pt[ipart]<<", "<<part_pt[ipart+1]<<"]"<<" is "<<meanpt<<endl;
  double eff = fhadeff->Eval(meanpt);
  //  cout<<"eff = "<<eff<<endl;
  return eff;
}

double MakeCFs::GetHadronEff_v2(TH1D* hadron_pt, int ipart)
{
  TF1* fhadeff = new TF1("fhadeff","(x<=3)*(0.396-0.337*exp(-1.50*x))+(x>3&&x<=5)*(0.400+-2.13e-13*exp(5.36*x))+(x>5)*(0.396-0.337*exp(-1.50*x))",0.5,10.0);
  int bin1 = hadron_pt->FindBin(part_pt[ipart]);
  int bin2 = hadron_pt->FindBin(part_pt[ipart+1]); 
  double binwidth = hadron_pt->GetBinWidth(1);
  
  cout<<"bin1 = "<<bin1<<"; bin2 = "<<bin2<<"; bin width = "<<binwidth<<endl;

  double hadeff = 0.0;
  double meanpt = 0.0;
  for(int ibin=bin1; ibin<bin2; ibin++){
    cout<<"ibin = "<<ibin<<endl;
    double pt = hadron_pt->GetBinCenter(ibin);
    cout<<"pt = "<<pt<<endl;
    double eff = fhadeff->Eval(pt);
    cout<<"eff = "<<eff<<" at pt = "<<pt<<endl;
    hadeff+=pt*eff;
    meanpt+=pt;
  }
  hadeff/=meanpt;
  //  cout<<"hadeff = "<<hadeff<<endl;
  return hadeff;
}
