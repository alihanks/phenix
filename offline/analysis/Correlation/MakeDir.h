#ifndef __MAKEDIR_H__
#define __MAKEDIR_H__

#include <MakeCFs.h>
#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>

class TH1D;
class TH2D;
class TH3D;
class TFile;
class TGraphErrors;
class TCanvas;
class TLegend;
class TLatex;

class MakeDir
{
public:
  MakeDir(const std::string Rgamma_input, const std::string finc, const std::string fdec, const std::string fout, int ic);
  virtual ~MakeDir(){};

  void SetRgamma(const std::string Rgamma_input, int ic);
  void DoSubtraction(TH1D* incjet, TH1D* decjet, double Rgamma, TH1D*& dirjet);

private:
  std::ostringstream bin;
  std::string name;
  TFile* fileinc;
  TFile* filedec;
  TFile* fileout;
  TFile* Rgammafile;
  double rgamma[NCENTBIN][NTRIGBIN];
  double rgamma_err[NCENTBIN][NTRIGBIN];//stat+sys
  double rgamma_stat[NCENTBIN][NTRIGBIN];
  double rgamma_sys[NCENTBIN][NTRIGBIN];
  TGraphErrors* gr[NCENTBIN];
  TGraphErrors* stat[NCENTBIN];
  TGraphErrors* sys[NCENTBIN];

  TH1D* inc_jet[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* dec_jet[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* dir_jet[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* dir_sub_err[NCENTBIN][NTRIGBIN][NPARTBIN];
  TH1D* dir_jet_err[NCENTBIN][NTRIGBIN][NPARTBIN];
};

#endif
