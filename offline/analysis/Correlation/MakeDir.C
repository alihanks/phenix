#include "MakeDir.h"

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

MakeDir::MakeDir(const string Rgamma_input, const string finc, const string fdec, const string fout)
{
  fileinc = new TFile (finc.c_str());
  filedec = new TFile (fdec.c_str());
  fileout = new TFile (fout.c_str(),"recreate");
  fileout->cd();

  for(int ic=0; ic<NCENTBIN; ic++){
    SetRgamma(Rgamma_input,ic);
    for(int itrig=0; itrig<NTRIGBIN; itrig++){
      for(int ipart=0; ipart<NPARTBIN; ipart++){
      	bin.str("");
      	bin << "JF_c"<<ic<<"_p"<<itrig<<"_h"<<ipart;
      	inc_jet[ic][itrig][ipart] = new TH1D(*(TH1D*)fileinc->Get(bin.str().c_str()));
      	name = "INC_" + bin.str();
      	inc_jet[ic][itrig][ipart]->SetName(name.c_str());
      	dec_jet[ic][itrig][ipart] = new TH1D(*(TH1D*)filedec->Get(bin.str().c_str()));
      	name = "DEC_" + bin.str();
      	dec_jet[ic][itrig][ipart]->SetName(name.c_str());
      	DoSubtraction(inc_jet[ic][itrig][ipart],dec_jet[ic][itrig][ipart],rgamma[ic][itrig],dir_jet[ic][itrig][ipart]);
      	name = "DIR_" + bin.str();
      	dir_jet[ic][itrig][ipart]->SetName(name.c_str());
      	dir_jet[ic][itrig][ipart]->Write();
      }
    }
  }
}

void MakeDir::SetRgamma(string Rgamma_input, int ic)
{
  Rgammafile = new TFile(Rgamma_input.c_str());
  bin.str("");
  bin<<ic;
  name = "gr" + bin.str();
  gr[ic] = new TGraphErrors(*(TGraphErrors*)Rgammafile->Get(name.c_str()));
  double *rga = gr[ic]->GetY();
  double *rga_err = gr[ic]->GetEY();
  for(int ip=0; ip<NTRIGBIN; ip++){
    rgamma[ic][ip] = rga[ip];
    int trigpt_bin_mod = ip;
    if( trigpt_bin_mod==4 ) trigpt_bin_mod=3;
    float sep_eff=1.;
    //for n=7 (from PISA) w/ eta already exluded
    // did not put in the 5-10 bin since we don't use it
    if(trigpt_bin_mod==0) sep_eff = 0.998;
    if(trigpt_bin_mod==1) sep_eff = 0.992;
    if(trigpt_bin_mod==2) sep_eff = 0.965;
    if(trigpt_bin_mod==3) sep_eff = 0.900;
    rgamma[ic][ip] = 1/sep_eff*(rgamma[ic][ip] -1.0 +sep_eff);

    rgamma_err[ic][ip] = rga_err[ip];
  }

  bin.str("");
  bin<<ic;
  name = "stat" + bin.str();
  stat[ic] = new TGraphErrors(*(TGraphErrors*)Rgammafile->Get(name.c_str()));
  double *rga_stat = stat[ic]->GetEY();
  for(int ip=0; ip<NTRIGBIN; ip++) rgamma_stat[ic][ip] = rga_stat[ip];

  bin.str("");
  bin<<ic;
  name = "sys" + bin.str();
  sys[ic] = new TGraphErrors(*(TGraphErrors*)Rgammafile->Get(name.c_str()));
  double *rga_sys = sys[ic]->GetEY();
  for(int ip=0; ip<NTRIGBIN; ip++) rgamma_sys[ic][ip] = rga_sys[ip];
}

void MakeDir::DoSubtraction(TH1D* incjet, TH1D* decjet, double Rgamma, TH1D*& dirjet)
{
  //dir = (Rgamma*inc-dec)/(Rgamma-1)
  TH1D* incjetclone = new TH1D(*(TH1D*)incjet);
  incjetclone->Scale(Rgamma);
  incjetclone->Add(decjet, -1.0);
  incjetclone->Scale(1.0/(Rgamma-1));
  dirjet = (TH1D*)incjetclone;
}
