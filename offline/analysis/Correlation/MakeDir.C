#include <MakeDir.h>

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

MakeDir::MakeDir(const string Rgamma_input, const string finc, const string fdec, const string fout, int ic)
{
  fileinc = new TFile (finc.c_str());
  filedec = new TFile (fdec.c_str());
  fileout = new TFile (fout.c_str(),"recreate");

  ostringstream prefix, prefix_err;
  if( ic < 0 ) {
    prefix << "JF";
    prefix_err << "JFerr";
    ic = 0;
  }
  else {
    prefix << "JF_c" << ic;
    prefix_err << "JFerr_c" << ic;
  }
  SetRgamma(Rgamma_input,ic);
  fileout->cd();
  for(int itrig=0; itrig<NTRIGBIN; itrig++){
    for(int ipart=0; ipart<NPARTBIN; ipart++){
    	bin.str("");
    	bin << prefix.str()<<"_p"<<itrig<<"_h"<<ipart;
    	inc_jet[ic][itrig][ipart] = new TH1F(*(TH1F*)fileinc->Get(bin.str().c_str()));
    	name = "INC_" + bin.str();
    	inc_jet[ic][itrig][ipart]->SetName(name.c_str());
    	dec_jet[ic][itrig][ipart] = new TH1F(*(TH1F*)filedec->Get(bin.str().c_str()));
    	name = "DEC_" + bin.str();
    	dec_jet[ic][itrig][ipart]->SetName(name.c_str());
    	DoSubtraction(inc_jet[ic][itrig][ipart],dec_jet[ic][itrig][ipart],rgamma[ic][itrig],dir_jet[ic][itrig][ipart]);
    	name = "DIR_" + bin.str();
    	dir_jet[ic][itrig][ipart]->SetName(name.c_str());
    	dir_jet[ic][itrig][ipart]->Write();
      DoSubtraction(inc_jet[ic][itrig][ipart],dec_jet[ic][itrig][ipart],rgamma[ic][itrig]+rgamma_err[ic][itrig],dir_sub_err[ic][itrig][ipart]);
      name = "DIRerr_" + bin.str();
      dir_sub_err[ic][itrig][ipart]->SetName(name.c_str());
      dir_sub_err[ic][itrig][ipart]->Write();
    }
  }
  for(int itrig=0; itrig<NTRIGBIN; itrig++){
    for(int ipart=0; ipart<NPARTBIN; ipart++){
      bin.str("");
      bin << prefix_err.str()<<"_p"<<itrig<<"_h"<<ipart;
      inc_jet[ic][itrig][ipart] = new TH1F(*(TH1F*)fileinc->Get(bin.str().c_str()));
      name = "INC_" + bin.str();
      inc_jet[ic][itrig][ipart]->SetName(name.c_str());
      dec_jet[ic][itrig][ipart] = new TH1F(*(TH1F*)filedec->Get(bin.str().c_str()));
      name = "DEC_" + bin.str();
      dec_jet[ic][itrig][ipart]->SetName(name.c_str());
      DoSubtraction(inc_jet[ic][itrig][ipart],dec_jet[ic][itrig][ipart],rgamma[ic][itrig],dir_jet_err[ic][itrig][ipart]);
      name = "DIR_" + bin.str();
      dir_jet_err[ic][itrig][ipart]->SetName(name.c_str());
      dir_jet_err[ic][itrig][ipart]->Write();
    }
  }
  for( int i=0; i<2; i++) {
    for(int ipart=0; ipart<NPARTBIN; ipart++) {
      ostringstream cname;
      cname << "DIR_JF_comb_p" << i << "_h" << ipart;
      TH1F* temp = new TH1F(*dir_jet[ic][i][ipart]);
      temp->SetName(cname.str().c_str());
      CombinePtBins(dir_jet[ic][i][ipart],dir_jet[ic][i+1][ipart],temp);
      temp->Write();
      cname.str("");
      cname << "DIR_JFerr_comb_p" << i << "_h" << ipart;
      temp = new TH1F(*dir_jet_err[ic][i][ipart]);
      temp->SetName(cname.str().c_str());
      CombinePtBins(dir_jet_err[ic][i][ipart],dir_jet_err[ic][i+1][ipart],temp);
      temp->Write();
      cname.str("");
      cname << "DIRerr_JF_comb_p" << i << "_h" << ipart;
      temp = new TH1F(*dir_sub_err[ic][i][ipart]);
      temp->SetName(cname.str().c_str());
      CombinePtBins(dir_sub_err[ic][i][ipart],dir_sub_err[ic][i+1][ipart],temp);
      temp->Write();
    }
  }
}

void MakeDir::CombinePtBins(TH1F* h1, TH1F* h2, TH1F* combined)
{
  for( int ib = 1; ib <= h1->GetNbinsX(); ib++ ) {
    double err1 = h1->GetBinError(ib);
    double err2 = h2->GetBinError(ib);
    double yield1 = h1->GetBinContent(ib);
    double yield2 = h2->GetBinContent(ib);

    if( err1 > 0 ) 
      yield1 = yield1*(err1*err1);
    else {
      yield1 = 0;
    }
    if( err2 > 0 )
      yield2 = yield2*(err2*err2);
    else {
      yield2 = 0;
    }

    double err = sqrt(err1*err1 + err2*err2);
    double yield = (yield1 + yield2);
    if( err > 0 ) yield = yield/(err*err); // weighted sum to account for pT dependent statistical uncertainties
    else yield = 0;

    combined->SetBinContent(ib,yield);
    combined->SetBinError(ib,err);
  }
}

void MakeDir::SetRgamma(string Rgamma_input, int ic)
{
  Rgammafile = new TFile(Rgamma_input.c_str());
  bin.str("");
  bin << "_c" << ic;
  name = "gr" + bin.str();
  gr[ic] = new TGraphErrors(*(TGraphErrors*)Rgammafile->Get(name.c_str()));
  double *rga = gr[ic]->GetY();
  double *rga_err = gr[ic]->GetEY();
  for(int ip=0; ip<NTRIGBIN; ip++){
    cout << "Rg = " << rga[ip] << endl;
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

    rgamma_err[ic][ip] = 1/sep_eff*(rga_err[ip] -1.0 +sep_eff);
  }

  bin.str("");
  bin<< "_c" <<ic;
  name = "stat" + bin.str();
  stat[ic] = new TGraphErrors(*(TGraphErrors*)Rgammafile->Get(name.c_str()));
  double *rga_stat = stat[ic]->GetEY();
  for(int ip=0; ip<NTRIGBIN; ip++) rgamma_stat[ic][ip] = rga_stat[ip];

  bin.str("");
  bin<< "_c" <<ic;
  name = "sys" + bin.str();
  sys[ic] = new TGraphErrors(*(TGraphErrors*)Rgammafile->Get(name.c_str()));
  double *rga_sys = sys[ic]->GetEY();
  for(int ip=0; ip<NTRIGBIN; ip++) rgamma_sys[ic][ip] = rga_sys[ip];
}

void MakeDir::DoSubtraction(TH1F* incjet, TH1F* decjet, double Rgamma, TH1F*& dirjet)
{
  //dir = (Rgamma*inc-dec)/(Rgamma-1)
  TH1F* incjetclone = new TH1F(*(TH1F*)incjet);
  incjetclone->Scale(Rgamma);
  incjetclone->Add(decjet, -1.0);
  incjetclone->Scale(1.0/(Rgamma-1));
  dirjet = (TH1F*)incjetclone;
}
