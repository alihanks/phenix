#ifndef __MAKEJFS_H__
#define __MAKEJFS_H__

#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>
//#include <fstream>

class TH1D;
class TH2D;
class TH3D;
class TF1;
class TFile;
class TGraph;
class TGraphErrors;
class TCanvas;
class TLegend;
class TLatex;

#define NCENTBIN 2
#define NTRIGBIN 4
#define NPARTBIN 5

//const double PI = TMath::Pi();
// const int NCENTBINS = 4;
// const int NTRIGBINS = 4;
// const int NPARTBINS = 5;

class MakeJFs
{
public:
  MakeJFs(int type, int centbin, int trigbin, int partbin, TH1D *CFinc, double meanpart, double ntrigbg, TFile *fin,  const std::string v2input, int nFits, int useMSMP, TH1D*& CFflowZYAM, TH1D*& CFjetZYAM);
  virtual ~MakeJFs(){};

  void InitHistos(TH1D* CFflow, std::string name);
  //  double GetNTriggers(TH1D* trigpt, double trigptmin, double trigptmax);
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

  //TH1D* CFflowZYAM;
  //  TH1D* CFjetZYAM;
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
