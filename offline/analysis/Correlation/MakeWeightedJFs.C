#include "MakeWeightedJFs.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TVirtualPad.h>

using namespace std;

MakeWeightedJFs::MakeWeightedJFs(const string fin, const string fout)
{
	infile = new TFile(fin.c_str());
	outfile = new TFile(fout.c_str(),"recreate");
}

void MakeWeightedJFs::GetHistos(int type)
{
	for( int ic = 0; ic < NCENT; ic++ )
	{
		cout << "Getting " << prefix << " histograms for centrality bin " << ic << endl;
		Get1DOutputHistos(type,ic);
	}
	cout << "Merging " << prefix << " centrality/pt binned histograms" << endl;
	GetMergedHistos(type);
}

void MakeWeightedJFs::Get1DOutputHistos(int type, int cbin)
{
	ostringstream bin;
	string name, mix_name;
	bin.str("");
	bin << "_c" << cbin;

	name = trig_name + bin.str();
	h1_trigpt[cbin] = (TH1F*)infile->Get(name.c_str());

	if( !type ) MakeDphiFrom3D(h1_trigpt[cbin],cbin);
	else MakeDphiFrom2D(h1_trigpt[cbin],cbin);

	h1_trigpt[cbin]->Write();
}

void MakeWeightedJFs::GetMergedHistos(int type)
{
	ostringstream bin;
	string name;
	TH1F* trigpt_combined = new TH1F(*h1_trigpt[0]);
	trigpt_combined->SetName("h1_trigpt_all");
	trigpt_combined->Reset();
	cout << "Merging " << prefix << " 1D dphi centrality histograms" << endl;
	for( int ic = 0; ic < NCENT; ic++ ) {
		trigpt_combined->Add(h1_trigpt[ic]);
		for( int it = 0; it < NTRIGBIN; it++ ) {
			for( int ih = 0; ih < NPARTBIN; ih++ ) {
				if( ic==0 ) {
					dphi_comb[it][ih] = new TH1F(*dphi_1d[ic][it][ih]);
					bin.str("");
					bin << prefix << "_p" << it << "_h" << ih;
					name = "JF_" + bin.str();
					dphi_comb[it][ih]->SetName(name.c_str());
				}
				else dphi_comb[it][ih]->Add(dphi_1d[ic][it][ih]);
			}
		}
	}
	outfile->cd();
	double ntrigs_comb[2] = {0,0};
	TH1F* dphi_pt_comb[2][NPARTBIN];
	TH1F* jf_pt_comb[2][NPARTBIN];
	for( int it = 0; it < 2; it++ ) {
		for( int ih = 0; ih < NPARTBIN; ih++ ) {
			bin.str("");
			bin << prefix << "_" << it << "_h" << ih;
			name = "JF_" + bin.str();
			dphi_pt_comb[it][ih] = new TH1F(*dphi_comb[0][ih]);
			dphi_pt_comb[it][ih]->SetName(name.c_str());
			dphi_pt_comb[it][ih]->Reset();
		}
	}

	cout << "Subtrackting background for " << prefix << " merged centrality histograms" << endl;
	for( int it = 0; it < NTRIGBIN; it++ ) {
		double ntrig_tot = GetNTrigs(type,it,trigpt_combined);
		if( it==0 || it==1 ) ntrigs_comb[0] += ntrig_tot;
		if( it==2 || it==3 ) ntrigs_comb[1] += ntrig_tot;

		for( int ih = 0; ih < NPARTBIN; ih++ ) {
			if( it==0 || it==1 ) jf_pt_comb[0][ih]->Add(dphi_comb[it][ih]);
			if( it==2 || it==3 ) jf_pt_comb[1][ih]->Add(dphi_comb[it][ih]);
			SubtractBackground(dphi_comb[it][ih], jf_comb[it][ih], dphi_comb[it][ih]->GetName());
			jf_comb[it][ih]->Scale(1/ntrig_tot);
			jf_comb[it][ih]->Write();
		}
	}

	cout << "Subtrackting background for " << prefix << " merged pt histograms" << endl;
	for( int it = 0; it < 2; it++ ){
		for( int ih = 0; ih < NPARTBIN; ih++ ) {
			SubtractBackground(dphi_pt_comb[it][ih],jf_pt_comb[it][ih], dphi_pt_comb[it][ih]->GetName());
			jf_pt_comb[it][ih]->Scale(1/ntrigs_comb[it]);
			jf_pt_comb[it][ih]->Write();
		}
	}
}

void MakeWeightedJFs::MakeDphiFrom3D(TH1F* trigpt, int cbin)
{
	ostringstream bin;
	string name, name_mix;
	bin.str("");
	bin << "_c" << cbin;

	name = dphi_name + bin.str();
	name_mix = dphi_mix_name + bin.str();
	TH3F* temp3D = (TH3F*)infile->Get(name.c_str());
	TH3F* temp3D_mix = (TH3F*)infile->Get(name_mix.c_str());
	outfile->cd();

	for(int it = 0; it < NTRIGBIN; it++){
		bin.str("");
		bin << "_p" << it << "_c" << cbin;
		name = "h1_part_pt_" + prefix + bin.str();
		TH1F* h1_trigpt = (TH1F*)temp3D->ProjectionX();
		int lbin = h1_trigpt->FindBin(trig_pt[it]);
		int hbin = h1_trigpt->FindBin(trig_pt[it+1]);
		h1_partpt[cbin][it] = (TH1F*)temp3D->ProjectionY(name.c_str(),lbin,hbin-1);
		h1_partpt[cbin][it]->Write();
		double ntrigs = GetNTrigs(0,it,trigpt);

		for(int ih = 0; ih < NPARTBIN; ih++){
			bin.str("");
			bin << "_c" << cbin << "_p"<<it<<"_h"<<ih;
			name = "h1_dphi_" + prefix + bin.str();
			name_mix = "h1_dphi_mix_" + prefix + bin.str();

			MakeDphiProjection(temp3D,dphi_1d[cbin][it][ih],trig_pt[it], trig_pt[it+1], part_pt[ih], part_pt[ih+1],name);
			MakeDphiProjection(temp3D_mix, dphi_1d_mix[cbin][it][ih], trig_pt[it], trig_pt[it+1], part_pt[ih], part_pt[ih+1], name_mix);
			dphi_1d[cbin][it][ih]->Write();
			dphi_1d_mix[cbin][it][ih]->Write();

			MakeJetFunction(dphi_1d[cbin][it][ih], corr[cbin][it][ih], ntrigs, it, ih, cbin);
			corr[cbin][it][ih]->Write();
		}
	}
}

void MakeWeightedJFs::MakeDphiFrom2D(TH1F* trigpt, int cbin)
{
	ostringstream bin;
	string name, name_mix;
	bin.str("");
	bin << "_c" << cbin;
	outfile->cd();

	for(int it = 0; it < NTRIGBIN; it++){
		bin.str("");
		if(it<3) bin <<"_p"<<it<<"_c"<<cbin;	   
		else bin<<"_p"<<it+1<<"_c"<<cbin;

		name = dphi_name + bin.str();
		name_mix = dphi_mix_name + bin.str();
		TH2F* temp2D = (TH2F*)infile->Get(name.c_str());
		TH2F* temp2D_mix = (TH2F*)infile->Get(name_mix.c_str());

		bin.str("");
		bin <<"_p"<<it<<"_c"<<cbin;
		name = "h1_part_pt_" + prefix + bin.str();
		h1_partpt[cbin][it] = (TH1F*)temp2D->ProjectionY(name.c_str());
		h1_partpt[cbin][it]->Write();
		double ntrigs = GetNTrigs(1,it,trigpt);

		for(int ih = 0; ih < NPARTBIN; ih++){
			bin.str("");
			bin << "_c" << cbin << "_p"<<it<<"_h"<<ih;
			name = "h1_dphi_" + prefix + bin.str();
			name_mix = "h1_dphi_mix_" + prefix + bin.str();

			Make2DDphiProjection(temp2D,dphi_1d[cbin][it][ih], part_pt[ih], part_pt[ih+1],name);
			Make2DDphiProjection(temp2D_mix,dphi_1d_mix[cbin][it][ih], part_pt[ih], part_pt[ih+1],name_mix);
			dphi_1d[cbin][it][ih]->Write();
			dphi_1d_mix[cbin][it][ih]->Write();

			MakeJetFunction(dphi_1d[cbin][it][ih], corr[cbin][it][ih], ntrigs, it, ih, cbin);
			corr[cbin][it][ih]->Write();
		}
	}

}

double MakeWeightedJFs::GetNTrigs(int type, int bin, TH1F* histo)
{
	if( !type ) {
		int lbin = histo->FindBin(trig_pt[bin]);
		int hbin = histo->FindBin(trig_pt[bin+1]);
		return histo->Integral(lbin,hbin);
	}
	if( type ) {
		if(type<4)
			return histo->GetBinContent(bin+1);
		else
			return histo->GetBinContent(bin+2);
	}
	return 0;
}

void MakeWeightedJFs::MakeDphiProjection(TH3F* h3, TH1F*& h1,double xmin, double xmax, double ymin, double ymax, string hname)
{
	TH1F* proj_x = (TH1F*)h3->ProjectionX("px");
	TH1F* proj_y = (TH1F*)h3->ProjectionY("py");
	int xbinlo = proj_x->FindBin(xmin);
	int xbinhi = proj_x->FindBin(xmax);
	int ybinlo = proj_y->FindBin(ymin);
	if( XiBinning ) ybinlo = proj_y->FindBin(ymax);
	int ybinhi = proj_y->FindBin(ymax);
	if( XiBinning ) ybinhi = proj_y->FindBin(ymin);
	h1 = new TH1F(*(TH1F*)h3->ProjectionZ(hname.c_str(),xbinlo,xbinhi-1,ybinlo,ybinhi-1));
	h1->SetName(hname.c_str());
}

void MakeWeightedJFs::Make2DDphiProjection(TH2F* h3, TH1F*& h1,double ymin, double ymax, string hname)
{
	TH1F* proj_y = (TH1F*)h3->ProjectionY("py");
	int ybinlo = proj_y->FindBin(ymin);
	if( XiBinning ) ybinlo = proj_y->FindBin(ymax);
	int ybinhi = proj_y->FindBin(ymax);
	if( XiBinning ) ybinhi = proj_y->FindBin(ymin);

	h1 = new TH1F(*(TH1F*)h3->ProjectionX(hname.c_str(),ybinlo,ybinhi-1));
	h1->SetName(hname.c_str());
}

void MakeWeightedJFs::MakeJetFunction(TH1F* dphi, TH1F*& correlation, double ntrigs, int it, int ih, int cbin)
{
	ostringstream name;
	name << "JF_" << prefix << "_c" << cbin << "_p" << it << "_h" << ih; 
	SubtractBackground(dphi, correlation, name.str());
	correlation->Scale(1/ntrigs);
}

void MakeWeightedJFs::SubtractBackground(TH1F* foreground, TH1F*& signal, string name)
{
	double norm = GetZYAMNorm(foreground);
	signal = new TH1F(*foreground);
	signal->SetName(name.c_str());
	TF1* bgFunc = new TF1("bgFunc","[0]",0.0,PI);
	bgFunc->SetParameter(0,norm);
	signal->Add(bgFunc,-1.0);
}

double MakeWeightedJFs::GetZYAMNorm(TH1F* dphi)
{
	dphi->SetAxisRange(1.0,1.6,"X");
	int bin = dphi->GetMinimumBin();
	dphi->SetAxisRange(0.0,TMath::Pi(),"X");
	int lbin = bin-2;//CFinc->FindBin(1.1);
	int hbin = bin+2;
	double norm = dphi->Integral(lbin,hbin);
	norm = norm/((double)(hbin-lbin+1));

	return norm;
}
