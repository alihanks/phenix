#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>

using namespace std;


//("combine_pt_inc_040_taxi4759.root","combine_pt_inc_040_run10.root","combine_pt_inc_040_run7.root","eval_err.root")
void eval_err(string f11, string f10, string f7, string fout)
{
  const int NPT_TRIG = 2;
  const int NPT_PART = 5;
  const int NDPHI = 30;
  const double PI = TMath::Pi();

  double trig_pt[NPT_TRIG+1] = {5.0,9.0,15.0};
  double part_pt[NPT_PART+1] = {0.5,1.0,2.0,3.0,5.0,7.0};

  ostringstream bin;
  ostringstream latex_name;
  string name;

  TCanvas* c1 = new TCanvas("can_combined","can_combined");
  c1->Divide(2,2,0.001,0.001);
  TFile* fin11 = new TFile(f11.c_str());
  TFile* fin10 = new TFile(f10.c_str());
  TFile* fin7 = new TFile(f7.c_str());
  TFile* foutput = new TFile(fout.c_str(),"recreate");

  TH1D* h11[NPT_TRIG][NPT_PART];//[NPT_TRIG][NPT_PART]
  TH1D* h10[NPT_TRIG][NPT_PART];
  TH1D* h7[NPT_TRIG][NPT_PART];
  double yerr11[NPT_TRIG][NPT_PART][NDPHI/2];
  double yerr10[NPT_TRIG][NPT_PART][NDPHI/2];
  double yerr7[NPT_TRIG][NPT_PART][NDPHI/2];
  double sum11[NPT_TRIG][NPT_PART]={0.0};
  double sum10[NPT_TRIG][NPT_PART]={0.0};
  double sum7[NPT_TRIG][NPT_PART]={0.0};
  // double avg11[NPT_TRIG]={0.0};
  // double avg10[NPT_TRIG]={0.0};
  // double avg7[NPT_TRIG]={0.0};

  for(int itrig=0; itrig<NPT_TRIG; itrig++){
    bin.str("");
    bin << itrig;
    name = "can_p" + bin.str();
    
    for(int ipart=0; ipart<NPT_PART; ipart++){
      bin.str("");
      bin << itrig <<"_h"<<ipart;
      name = "JF_COMBTOT_c0_p" + bin.str();
      h11[itrig][ipart] = new TH1D(*(TH1D*)fin11->Get(name.c_str()));
      name = "JFrun10_COMBTOT_c0_p" + bin.str();
      h10[itrig][ipart] = new TH1D(*(TH1D*)fin10->Get(name.c_str()));
      name = "JFrun7_COMBTOT_c0_p" + bin.str();
      h7[itrig][ipart] = new TH1D(*(TH1D*)fin7->Get(name.c_str()));

      for(int idphi=0; idphi<15; idphi++){
	yerr11[itrig][ipart][idphi] = h11[itrig][ipart]->GetBinError(idphi+16);
	yerr10[itrig][ipart][idphi] = h10[itrig][ipart]->GetBinError(idphi+16);
	yerr7[itrig][ipart][idphi] = h7[itrig][ipart]->GetBinError(idphi+16);
	sum11[itrig][ipart] += pow(yerr11[itrig][ipart][idphi],2.0);
	sum10[itrig][ipart] += pow(yerr10[itrig][ipart][idphi],2.0);
	sum7[itrig][ipart] += pow(yerr7[itrig][ipart][idphi],2.0);
      }
      double interr11 = sqrt(sum11[itrig][ipart]);
      double interr10 = sqrt(sum10[itrig][ipart]);
      double interr7 = sqrt(sum7[itrig][ipart]);
      cout<<"interr11["<<itrig<<"]["<<ipart<<"]="<<interr11[itrig][ipart]<<endl;
      cout<<"interr10["<<itrig<<"]["<<ipart<<"]="<<interr10[itrig][ipart]<<endl;
      cout<<"interr7["<<itrig<<"]["<<ipart<<"]="<<interr7[itrig][ipart]<<endl;

      double ratio = sqrt(interr10*interr10+interr7*interr7)/sqrt(interr11*interr11+interr10*interr10+interr7*interr7);
      cout<<"ratio["<<itrig<<"]["<<ipart<<"] = "<<ratio<<endl;
      // avg11[itrig]+=sum11[itrig][ipart];
      // avg10[itrig]+=sum10[itrig][ipart];
      // avg7[itrig]+=sum7[itrig][ipart];

      if(ipart==3){
	// TVirtualPad* pad = 
	c1->cd(2*itrig+1);
	// pad->SetRightMargin(0.01);
	// pad->SetTopMargin(0.05);
	//h11[itrig][ipart]->GetXaxis()->SetTitle("#Delta#phi");
	//h11[itrig][ipart]->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN}{d#Delta#phi}");
	h11[itrig][ipart]->SetTitle(";#Delta#phi [rad];Arb.Unit");
	h11[itrig][ipart]->GetYaxis()->SetTitleOffset(0.5);
	h11[itrig][ipart]->GetYaxis()->SetLabelSize(0.);
	h11[itrig][ipart]->Rebin(3);
	h11[itrig][ipart]->Draw();
	latex_name.str("");
	latex_name << "0 - 40% 2011 Au + Au";
	TLatex *la = new TLatex(0.25, 0.7, latex_name.str().c_str());
	la->SetNDC();
	la->Draw("same");
	latex_name.str("");
	latex_name <<trig_pt[itrig] << " < p_{T}^{#gamma} < "<<trig_pt[itrig+1] <<" GeV/c #times "<<part_pt[ipart]<<" < p_{T}^{h} < "<<part_pt[ipart+1]<<" GeV/c";
	TLatex *la2 = new TLatex(0.25, 0.6, latex_name.str().c_str());
	la2->SetNDC();
	la2->Draw("same");
	// TVirtualPad* pad2 = 
	c1->cd(2*itrig+2);
	// pad2->SetTopMargin(0.05);
	h10[itrig][ipart]->SetMarkerColor(2);
	h10[itrig][ipart]->SetLineColor(2);
	h10[itrig][ipart]->SetTitle(";#Delta#phi [rad];Arb.Unit");
	h10[itrig][ipart]->GetYaxis()->SetLabelSize(0.);
	// h10[itrig][ipart]->GetXaxis()->SetTitle("#Delta#phi");
	// h10[itrig][ipart]->GetYaxis()->SetTitle("#frac{1}{N_{trig}}#frac{dN}{d#Delta#phi}");
	h10[itrig][ipart]->GetYaxis()->SetTitleOffset(0.5);
	h10[itrig][ipart]->Rebin(3);
	h10[itrig][ipart]->Draw();
	latex_name.str("");
	latex_name << "#color[2]{0 - 40% 2010 Au + Au}";
	TLatex *la = new TLatex(0.25, 0.7, latex_name.str().c_str());
	la->SetNDC();
	la->Draw("same");
	latex_name.str("");
	latex_name <<"#color[2]{"<<trig_pt[itrig] << " < p_{T}^{#gamma} < "<<trig_pt[itrig+1] <<" GeV/c #times "<<part_pt[ipart]<<" < p_{T}^{h} < "<<part_pt[ipart+1]<<" GeV/c}";
	TLatex *la3 = new TLatex(0.25, 0.6, latex_name.str().c_str());
	la3->SetNDC();
	la3->Draw("same");
      }
    }
    c1->Write();
    // avg11[itrig]/=5;
    // avg10[itrig]/=5;
    // avg7[itrig]/=5;
    // cout<<"avg11["<<itrig<<"]="<<avg11[itrig]<<endl;
    // cout<<"avg10["<<itrig<<"]="<<avg10[itrig]<<endl;
    // cout<<"avg7["<<itrig<<"]="<<avg7[itrig]<<endl;
  }

}
