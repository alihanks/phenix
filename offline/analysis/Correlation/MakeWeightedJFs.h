#ifndef __MAKEWEIGHTEDJFS_H__
#define __MAKEWEIGHTEDJFS_H__

#include <TMath.h>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
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

	void GetHistos(int type);
	void SetTriggerBinning(std::vector<double> trig_bins) { trig_pt = trig_bins; }
	void SetPartnerBinning(std::vector<double> part_bins) { part_pt = part_bins; }
	void SetTriggerName(std::string name) { trig_name = name; }
	void SetDphiNames(std::string name, std::string mix_name) { dphi_name = name; dphi_mix_name = mix_name; }
	void SetTypeName(std::string type) { prefix = type; }

	int XiBinning;

private:

	void Get1DOutputHistos(int type, int cbin);
	void GetMergedHistos(int type);
	void MakeDphiFrom3D(TH1F* trigpt, int cbin);
	void MakeDphiFrom2D(TH1F* trigpt, int cbin);
	void MakeDphiProjection(TH3F* h3, TH1F*& h1, double xmin, double xmax, double ymin, double ymax, std::string hname);
	void Make2DDphiProjection(TH2F* h3, TH1F*& h1, double ymin, double ymax, std::string hname);
	void MakeJetFunction(TH1F* dphi, TH1F*& correlation, double ntrigs, int it, int ih, int cbin);
	double GetNTrigs(int type, int bin, TH1F* trigpt);
	void SubtractBackground(TH1F* foreground, TH1F*& signal, std::string name);
	double GetZYAMNorm(TH1F* dphi);

	std::vector<double> trig_pt;
	std::vector<double> part_pt;

	std::string prefix;
	std::string trig_name;
	std::string dphi_name;
	std::string dphi_mix_name;

	TFile* infile;
	TFile* outfile;

	TH1F* h1_trigpt[NCENT];
	TH1F* h1_partpt[NCENT][NTRIGBIN];
	TH1F* dphi_1d[NCENT][NTRIGBIN][NPARTBIN];
	TH1F* dphi_1d_mix[NCENT][NTRIGBIN][NPARTBIN];
	TH1F* corr[NCENT][NTRIGBIN][NPARTBIN];

};

#endif




