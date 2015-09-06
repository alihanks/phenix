#ifndef __MAKECFS_H__
#define __MAKECFS_H__

#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>

class TH1D;
class TH2D;
class TH3D;
class TFile;
class TCanvas;
class TVirtualPad;
class TLegend;
class TLatex;

const double PI = acos(-1.0);

#define NCENTBIN 4
#define NTRIGBIN 4
#define NPARTBIN 7

class MakeCFs
{
public:
  MakeCFs(const std::string fin, const std::string fout);
  virtual ~MakeCFs(){};

  void Run(int type, int ispertrigger);
  void SetTriggerName(std::string name) { trig_name = name; }
  void SetDphiNames(std::string name, std::string mix_name) { dphi_name = name; dphi_mix_name = mix_name; }

  void SetTrigPtBinning();
  void SetPartPtBinning();
  void SetPtRange(TH3D* h3, double x_pt_min, double x_pt_max, double y_pt_min, double y_pt_max);
  void MakeDphiProjection(TH3D* h3, TH1D*& h1, std::string hname);
  void MakeDphiProjection(TH3D* h3, TH1D*& h1, double xmin, double xmax, double ymin, double ymax, std::string hname);
  void FoldDphiDist(TH1D* h1, TH1D*& h1_fold, std::string hname_fold);
  void SetHisto(TH1D* h1, std::string hname, int color);
  void SetPad(TVirtualPad* pad);
  double GetNTriggers(TH1D* trigpt, double trigptmin, double trigptmax);
  double GetHadronEff(TH1D* hadron_pt, int ipart);
  double GetHadronEff_v2(TH1D* hadron_pt, int ipart);

private:
  std::string trig_name;
  std::string dphi_name;
  std::string dphi_mix_name;
  std::string name;
  std::string name_mix;
  std::ostringstream bin;
  std::ostringstream can_dphi_name;
  std::ostringstream can_corr_name;
  std::ostringstream legend_name;
  std::ostringstream corr_name;

  TFile* infile;
  TFile* outfile;
  TH3D* temp3D;
  TH3D* temp3D_mix;
  TH2D* temp2D;
  TH2D* temp2D_mix;
  double trig_pt[NTRIGBIN+1];
  double part_pt[NPARTBIN+1];

  double normint;
  double num_trigger[NCENTBIN][NTRIGBIN];
  double num_trigger_mix[NCENTBIN][NTRIGBIN];

  TH3D* dphi_3d[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH3D* dphi_3d_mix[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH2D* dphi_2d[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH2D* dphi_2d_mix[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* dphi_1d[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* dphi_1d_mix[NCENTBIN][NTRIGBIN][NPARTBIN];
  
  TH1D* fold[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* fold_mix[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* corr[NCENTBIN][NTRIGBIN][NPARTBIN];
  
  TCanvas* can_dphi[NCENTBIN][NTRIGBIN];
  TCanvas* can_corr[NCENTBIN][NTRIGBIN];
  TCanvas* can_abs[NCENTBIN][NTRIGBIN];
  TCanvas* can_jet[NCENTBIN][NTRIGBIN];

  // TGraphErrors* XICORR;
  // TGraphErrors* XICORRLARGEBIN;
  double meanpart[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* jet[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* flow[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* jet_err[NCENTBIN][NTRIGBIN][NPARTBIN];

};

#endif




