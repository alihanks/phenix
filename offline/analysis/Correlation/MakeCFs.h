#ifndef __MAKECFS_H__
#define __MAKECFS_H__

#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>

class TH1F;
class TH2F;
class TH3F;
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
  MakeCFs(const std::string fin, const std::string fout, int isxi);
  virtual ~MakeCFs(){};

  void Run(int type, int ispertrigger);
  void SetTriggerName(std::string name) { trig_name = name; }
  void SetDphiNames(std::string name, std::string mix_name) { dphi_name = name; dphi_mix_name = mix_name; }

  void SetTrigPtBinning();
  void SetPartPtBinning(int isxi);
  void SetPtRange(TH3F* h3, double x_pt_min, double x_pt_max, double y_pt_min, double y_pt_max);
  void MakeDphiProjection(TH3F* h3, TH1F*& h1, std::string hname);
  void MakeDphiProjection(TH3F* h3, TH1F*& h1, double xmin, double xmax, double ymin, double ymax, std::string hname);
  TH1F* MakeDphiProjection(TH3F* h3, float xmin, float xmax, float ymin, float ymax, std::string hname);
  void FoldDphiDist(TH1F* h1, TH1F*& h1_fold, std::string hname_fold);
  void SetHisto(TH1F* h1, std::string hname, int color);
  void SetPad(TVirtualPad* pad);
  double GetNTriggers(TH1F* trigpt, double trigptmin, double trigptmax);
  double GetHadronEff(TH1F* hadron_pt, int ipart);
  double GetHadronEff_v2(TH1F* hadron_pt, int ipart);

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
  TH3F* temp3D;
  TH3F* temp3D_mix;
  TH2F* temp2D;
  TH2F* temp2D_mix;
  double trig_pt[NTRIGBIN+1];
  double part_pt[NPARTBIN+1];
  int XiBinning;

  double normint;
  double num_trigger[NCENTBIN][NTRIGBIN];
  double num_trigger_mix[NCENTBIN][NTRIGBIN];

  TH3F* dphi_3d[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH3F* dphi_3d_mix[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH2F* dphi_2d[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH2F* dphi_2d_mix[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1F* dphi_1d[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1F* dphi_1d_mix[NCENTBIN][NTRIGBIN][NPARTBIN];
  
  TH1F* fold[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1F* fold_mix[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1F* corr[NCENTBIN][NTRIGBIN][NPARTBIN];
  
  TCanvas* can_dphi[NCENTBIN][NTRIGBIN];
  TCanvas* can_corr[NCENTBIN][NTRIGBIN];
  TCanvas* can_abs[NCENTBIN][NTRIGBIN];
  TCanvas* can_jet[NCENTBIN][NTRIGBIN];

  // TGraphErrors* XICORR;
  // TGraphErrors* XICORRLARGEBIN;
  double meanpart[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1F* jet[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1F* flow[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1F* jet_err[NCENTBIN][NTRIGBIN][NPARTBIN];

};

#endif




