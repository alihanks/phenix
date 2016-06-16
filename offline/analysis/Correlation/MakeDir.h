#ifndef __MAKEDIR_H__
#define __MAKEDIR_H__

<<<<<<< MakeDir.h
#include <MakeWeightedJFs.h>
=======
#include <MakeCFs.h>
>>>>>>> 1.2
#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>

class TH1F;
class TH2F;
class TH3F;
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
  void DoSubtraction(TH1F* incjet, TH1F* decjet, double Rgamma, TH1F*& dirjet);
  void CombinePtBins(TH1F* h1, TH1F* h2, TH1F* combined);

private:
  std::ostringstream bin;
  std::string name;
  TFile* fileinc;
  TFile* filedec;
  TFile* fileout;
  TFile* Rgammafile;
  double rgamma[NCENT][NTRIGBIN];
  double rgamma_err[NCENT][NTRIGBIN];//stat+sys
  double rgamma_stat[NCENT][NTRIGBIN];
  double rgamma_sys[NCENT][NTRIGBIN];
  TGraphErrors* gr[NCENT];
  TGraphErrors* stat[NCENT];
  TGraphErrors* sys[NCENT];

  TH1F* inc_jet[NCENT][NTRIGBIN][NPARTBIN];
  TH1F* dec_jet[NCENT][NTRIGBIN][NPARTBIN];
  TH1F* dir_jet[NCENT][NTRIGBIN][NPARTBIN];
  TH1F* dir_sub_err[NCENT][NTRIGBIN][NPARTBIN];
  TH1F* dir_jet_err[NCENT][NTRIGBIN][NPARTBIN];
};

#endif
