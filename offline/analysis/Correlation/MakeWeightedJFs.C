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
		name = "h1_part_pt" + bin.str();
		TH1F* h1_trigpt = (TH1F*)temp3D->ProjectionX();
		int lbin = h1_trigpt->FindBin(trig_pt[it]);
		int hbin = h1_trigpt->FindBin(trig_pt[it+1]);
		h1_partpt[cbin][it] = (TH1F*)temp3D->ProjectionY(name.c_str(),lbin,hbin-1);
		h1_partpt[cbin][it]->Write();

		for(int ih = 0; ih < NPARTBIN; ih++){
			bin.str("");
			bin << "_c" << cbin << "_p"<<it<<"_h"<<ih;
			name = "h1_dphi" + bin.str();
			name_mix = "h1_dphi_mix" + bin.str();

			MakeDphiProjection(temp3D,dphi_1d[cbin][it][ih],trig_pt[it], trig_pt[it+1], part_pt[ih], part_pt[ih+1],name);
			MakeDphiProjection(temp3D_mix, dphi_1d_mix[cbin][it][ih], trig_pt[it], trig_pt[it+1], part_pt[ih], part_pt[ih+1], name_mix);
			dphi_1d[cbin][it][ih]->Write();
			dphi_1d_mix[cbin][it][ih]->Write();

			MakeJetFunction(dphi_1d[cbin][it][ih], corr[cbin][it][ih], trigpt, it, ih, cbin);
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
		mix_name = dphi_mix_name + bin.str();
		TH2F* temp2D = (TH2F*)infile->Get(name.c_str());
		TH2F* temp2D_mix = (TH2F*)infile->Get(name_mix.c_str());

		bin.str("");
		bin <<"_p"<<it<<"_c"<<cbin;	   
		name = "h1_part_pt" + bin.str();
		h1_partpt[cbin][it] = (TH1F*)temp2D->ProjectionY(name.c_str());
		h1_partpt[cbin][it]->Write();

		for(int ih = 0; ih < NPARTBIN; ih++){
			bin.str("");
			bin << "_c" << cbin << "_p"<<it<<"_h"<<ih;
			name = "h1_dphi" + bin.str();
			name_mix = "h1_dphi_mix" + bin.str();

			Make2DDphiProjection(temp2D,dphi_1d[cbin][it][ih], part_pt[ih], part_pt[ih+1],name);
			Make2DDphiProjection(temp2D_mix,dphi_1d_mix[cbin][it][ih], part_pt[ih], part_pt[ih+1],name_mix);
			dphi_1d[cbin][it][ih]->Write();
			dphi_1d_mix[cbin][it][ih]->Write();

			MakeJetFunction(dphi_1d[cbin][it][ih], corr[cbin][it][ih], trigpt, it, ih, cbin);
		}
	}

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

void MakeWeightedJFs::MakeJetFunction(TH1F* dphi, TH1F*& correlation, TH1F* trigpt, int it, int ih, int cbin)
{
	int lbin = trigpt->FindBin(trig_pt[it]);
	int hbin = trigpt->FindBin(trig_pt[it+1]);
	double ntrigs = trigpt->Integral(lbin,hbin);
	ostringstream name;
	name << "JF_c" << cbin << "_p" << it << "_h" << ih; 
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
	float lphi = dphi->GetBinCenter(lbin);
	float hphi = dphi->GetBinCenter(hbin);
	double norm = dphi->Integral(lbin,hbin);
	norm = norm/((double)(hbin-lbin+1));
	cout << "ZYAM norm = " << dphi->Integral(lbin,hbin) << "/(" << hphi << " - " << lphi << ") = " << norm << endl;

	return norm;
}
