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

MakeCombinedHistos::MakeCombinedHistos(const string fin, const string fout, const string trig_name, int type)
{
	TFile* infile = new TFile(fin.c_str());
	TFile* outfile = new TFile(fout.c_str(),"recreate");

	cout<<"start"<<endl;
	ostringstream bin;
	string name;
	TH1D* JFhisto[NTRIGBIN][NPARTBIN];
	TH1D* JFerr[NTRIGBIN][NPARTBIN];
	TH1D* JFhisto_combined[2][NPARTBIN];
	TH1D* JFerr_combined[2][NPARTBIN];
	double pt_range[NTRIGBIN+1] = {5.0,7.0,9.0,12.0,15.0};
	double ntrig_total[NTRIGBIN] = {0};

	for(int ic=0; ic<NCENTBIN; ic++){
		bin.str("");
		bin << "_c" << ic;
		name = trig_name + bin.str();
		TH1D* trigpt = new TH1D(*(TH1D*)infile->Get(name.c_str()));
		cout << "getting trigger histo: " << name << endl;

		for(int itrig=0; itrig<NTRIGBIN; itrig++){
			double ntriggers = GetNTrig(trigpt,pt_range[itrig],pt_range[itrig+1]);
			if(type == 2) {
				if(itrig<3){
					ntriggers = trigpt->GetBinContent(itrig+1);
				}
				else{
					ntriggers = trigpt->GetBinContent(itrig+2);
				}
			}
			ntrig_total[itrig] += ntriggers;
			for(int ipart=0; ipart<NPARTBIN; ipart++){
				bin.str("");
				bin << "_c" << ic <<"_p"<<itrig<<"_h"<<ipart;
				name = "JF" + bin.str();
				TH1D* temp = (TH1D*)infile->Get(name.c_str());
				temp->Scale(ntriggers);

				bin.str("");
				bin << "_c" << ic <<"_p"<<itrig<<"_h"<<ipart;
				name = "JFerr" + bin.str();
				TH1D* tempe = (TH1D*)infile->Get(name.c_str());
				tempe->Scale(ntriggers);

				if(ic==0) {
					JFhisto[itrig][ipart] = new TH1D(*temp);
					bin.str("");
					bin << "_p" << itrig << "_h" << ipart;
					name = "JF" + bin.str();
					JFhisto[itrig][ipart]->SetName(name.c_str());

					JFerr[itrig][ipart] = new TH1D(*tempe);
					bin.str("");
					bin << "_p" << itrig << "_h" << ipart;
					name = "JFerr" + bin.str();
					JFerr[itrig][ipart]->SetName(name.c_str());
				}
				else {
					JFhisto[itrig][ipart]->Add(temp);
					JFerr[itrig][ipart]->Add(tempe);
				}
				if(ic==(NCENTBIN-1)) JFhisto[itrig][ipart]->Scale(1/ntrig_total[itrig]);
				if(ic==(NCENTBIN-1)) JFerr[itrig][ipart]->Scale(1/ntrig_total[itrig]);
			}
		}
	}
	outfile->cd();
	for( int i=0; i<2; i++) {
		for(int ipart=0; ipart<NPARTBIN; ipart++) {
			ostringstream cname;
			cname << "JF_comb_p" << i << "_h" << ipart;
			JFhisto_combined[i][ipart] = new TH1D(*JFhisto[i][ipart]);
			JFhisto_combined[i][ipart]->SetName(cname.str().c_str());
			CombinePtBins(JFhisto[i][ipart],JFhisto[i+1][ipart],JFhisto_combined[i][ipart]);
			JFhisto_combined[i][ipart]->Write();
			cname.str("");
			cname << "JFerr_comb_p" << i << "_h" << ipart;
			JFerr_combined[i][ipart] = new TH1D(*JFhisto[i][ipart]);
			JFerr_combined[i][ipart]->SetName(cname.str().c_str());
			CombinePtBins(JFerr[i][ipart],JFerr[i+1][ipart],JFerr_combined[i][ipart]);
			JFerr_combined[i][ipart]->Write();
		}
	}
	for(int itrig=0; itrig<NTRIGBIN; itrig++) {
		for(int ipart=0; ipart<NPARTBIN; ipart++) {
			JFhisto[itrig][ipart]->Write();
			JFerr[itrig][ipart]->Write();
		}
	}
}

void MakeCombinedHistos::CombinePtBins(TH1D* h1, TH1D* h2, TH1D* combined)
{
	for( int ib = 1; ib <= h1->GetNbinsX(); ib++ ) {
		double yield1 = h1->GetBinContent(ib);
		double err1 = h1->GetBinError(ib);
		double yield2 = h2->GetBinContent(ib);
		double err2 = h2->GetBinError(ib);

		double err = sqrt(err1*err1 + err2*err2);
		double yield = (yield1/(err1*err1) + yield2/(err2*err2))*err*err; // weighted sum to account for pT dependent statistical uncertainties

		combined->SetBinContent(ib,yield);
		combined->SetBinError(ib,err);
	}
}

double MakeCombinedHistos::GetNTrig(TH1D* trigpt, double trigptmin, double trigptmax)
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

