#include <MakeCombinedHistos.h>
#include <MakeCFs.h>
#include <MakeJFs.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualPad.h>

#include <string>
#include <sstream>
#include <iostream>

using namespace std;

MakeCombinedHistos::MakeCombinedHistos(const string fin, const string fout, const string trig_name)
{
	TFile* infile = new TFile(fin.c_str());
	TFile* outfile = new TFile(fout.c_str(),"recreate");

	cout<<"start"<<endl;
	ostringstream bin;
	string name;
	TH1D* JFhisto[NTRIGBIN][NPARTBIN];
	double pt_range[NTRIGBIN+1] = {5.0,7.0,9.0,12.0,15.0};
	double ntrig_total[NTRIGBIN] = {0};

	for(int ic=0; ic<NCENTBIN; ic++){
		bin.str("");
		bin << "_c" << ic;
		name = trig_name + bin.str();
		TH1D* trigpt = new TH1D(*(TH1D*)infile->Get(name.c_str()));

		for(int itrig=0; itrig<NTRIGBIN; itrig++){
			double ntriggers = GetNTriggers(trigpt,pt_range[itrig],pt_range[itrig+1]);
			ntrig_total[itrig] += ntriggers;
			for(int ipart=0; ipart<NPARTBIN; ipart++){
				bin.str("");
				bin << "_c" << ic <<"_p"<<itrig<<"_h"<<ipart;
				name = "JF" + bin.str();
				TH1D* temp = (TH1D*)infile->Get(name.c_str());
				temp->Scale(ntriggers);

				if(ic==0) {
					JFhisto[itrig][ipart] = new TH1D(*temp);
					bin.str("");
					bin << "_p" << itrig << "_h" << ipart;
					name = "JF" + bin.str();
					JFhisto[itrig][ipart]->SetName(name.c_str());
				}
				else {
					JFhisto[itrig][ipart]->Add(temp);
				}
				if(ic==(NCENTBIN-1)) JFhisto[itrig][ipart]->Scale(1/ntrig_total[itrig]);
			}
		}
	}
	outfile->cd();
	for(int itrig=0; itrig<NTRIGBIN; itrig++) {
		for(int ipart=0; ipart<NPARTBIN; ipart++) {
			JFhisto[itrig][ipart]->Write();
		}
	}
}

double MakeCFs::GetNTrig(TH1D* trigpt, double trigptmin, double trigptmax)
{
	double ntrig = 0;
	for(int i = 0; i < trigpt->GetNbinsX(); i++)
	{
		double bincenter = trigpt->GetBinCenter(i+1);
		if(bincenter >=trigptmin && bincenter <= trigptmax)
		{
			double bincontent = trigpt->GetBinContent(i+1);
			ntrig += bincontent;
		}
	}
	return ntrig;
}

