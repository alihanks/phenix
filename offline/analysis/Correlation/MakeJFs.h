#ifndef __MAKEJFS_H__
#define __MAKEJFS_H__

#include <MakeCFs.h>
#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>
//#include <fstream>

class TH1F;
class TF1;
class TFile;
class TGraph;
class TGraphErrors;
class TCanvas;
class TLegend;
class TLatex;

class MakeJFs
{
public:
  MakeJFs(int type, int centbin, int trigbin, int partbin, TH1F *CFinc, double meanpart, double ntrigbg, TFile *fin,  const std::string v2input, int nFits, int useMSMP, TH1F*& CFflowZYAM, TH1F*& CFjetZYAM);
  virtual ~MakeJFs(){};

  void InitHistos(TH1F* CFflow, std::string name);
  //  double GetNTriggers(TH1F* trigpt, double trigptmin, double trigptmax);
  void SetV2(const std::string v2_inputs);
  double GetZYAMScale(TF1* topFunc, TF1* bottomFunc);
  void EvalXi(int type, int trigbin, int partbin, TFile* fcentdist, TFile *fin, TGraphErrors*& XICORR, TGraphErrors*& XICORRLARGEBIN);

private:
  std::ostringstream name;
  TFile *v2file;
  TGraphErrors* gr_inc_v2[NCENTBIN];
  TGraphErrors* gr_dec_v2[NCENTBIN];
  TGraphErrors* gr_pi0_v2[NCENTBIN];
  TGraphErrors* gr_had_v2[NCENTBIN];
  TGraphErrors* gr_inc_v2sys[NCENTBIN];
  TGraphErrors* gr_dec_v2sys[NCENTBIN];
  TGraphErrors* gr_pi0_v2sys[NCENTBIN];
  TGraphErrors* gr_had_v2sys[NCENTBIN];

  double inc_v2[NCENTBIN][NTRIGBIN]; 
  double inc_v2_err[NCENTBIN][NTRIGBIN]; 
  double inc_v2_sys[NCENTBIN][NTRIGBIN];
  double dec_v2[NCENTBIN][NTRIGBIN];
  double dec_v2_err[NCENTBIN][NTRIGBIN];
  double dec_v2_sys[NCENTBIN][NTRIGBIN];
  double pi0_v2[NCENTBIN][NTRIGBIN]; 
  double pi0_v2_err[NCENTBIN][NTRIGBIN];
  double pi0_v2_sys[NCENTBIN][NTRIGBIN];
  double hadron_v2[NCENTBIN][NPARTBIN]; 
  double hadron_v2_err[NCENTBIN][NPARTBIN]; 
  double hadron_v2_sys[NCENTBIN][NPARTBIN];

  double trigv2;
  double partv2;
  double c2;

  //TH1F* CFflowZYAM;
  //  TH1F* CFjetZYAM;
  TF1 *flowFunc;//Flow fit
  TF1* cfFunc;//CF fit
  double scale;//ZYAM scale
  double norm;//bg normalization scale

  double xi[20][4];//[centrality][fits]; fits 0-exp&npart; 1-arctan&npart; 2-exp&ncoll; 3-arctan&ncoll;
  double xiavg[18];
  double xierr[18];
  double xiavglarge[4];
  double xierrlarge[4];
  double tpars[4][3];//[fits][parameters]
  double ppars[4][3];
  TFile* fcentdist;

  TGraphErrors *gXICORR;
  // TGraph *XIBIG;
  // TGraph *XISMALL;
  TGraphErrors *gXICORRLARGEBIN;
};

#endif
