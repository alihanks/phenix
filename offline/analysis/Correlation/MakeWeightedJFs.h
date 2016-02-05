#ifndef __MAKEWEIGHTEDJFS_H__
#define __MAKEWEIGHTEDJFS_H__

#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>
//#include <fstream>

class TH1F;
class TH2F;
class TH3F;
class TF1;
class TFile;
class TGraph;
class TGraphErrors;

const double PI = acos(-1.0);

#define NCENT 4
#define NTRIGBIN 4
#define NPARTBIN 7

class MakeWeightedJFs
{
public:
	MakeWeightedJFs(const std::string fin, const std::string fout);
	virtual ~MakeWeightedJFs(){};

	void Get1DOutputHistos(int type, int cbin);

private:

	void MakeDphiFrom3D(TH1F* trigpt, int cbin);
	void MakeDphiFrom2D(TH1F* trigpt, int cbin);
	void MakeDphiProjection(TH3F* h3, TH1F*& h1,double xmin, double xmax, double ymin, double ymax, string hname);
	void Make2DDphiProjection(TH2F* h3, TH1F*& h1,double ymin, double ymax, string hname);
	void MakeJetFunction(TH1F* dphi, TH1F*& correlation, TH1F* trigpt, int it, int ih, int cbin);
	void SubtractBackground(TH1F* foreground, TH1F*& signal, std::string name);
	double GetZYAMNorm(TH1F* dphi);

	double trig_pt[NTRIGBIN+1];
	double part_pt[NPARTBIN+1];
	int XiBinning;

	std::string trig_name;
	std::string dphi_name;
	std::string dphi_mix_name;

	TH1F* h1_trigpt[NCENT];
	TH1F* h1_partpt[NCENT][NTRIGBIN];
	TH1F* dphi_1d[NCENT][NTRIGBIN][NPARTBIN];
	TH1F* dphi_1d_mix[NCENT][NTRIGBIN][NPARTBIN];
	TH1F* corr[NCENT][NTRIGBIN][NPARTBIN];

}


