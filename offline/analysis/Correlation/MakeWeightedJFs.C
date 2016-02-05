#include "MakeWeightedJFs.h"

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

MakeWeightedJFs::MakeWeightedJFs(const string fin, const string fout)
{
	infile = new TFile(fin.c_str());
	outfile = new TFile(fout.c_str(),"recreate");
}

void MakeWeightedJFs::Get1DOutputHistos()
{
	TH1F* h1_trigpt = (TH1F*)infile->Get(trig_name);

	TH3F* temp3D = (TH3F*)infile->Get(inc_name);
	TH3F* temp3D_mix = (TH3F*)infile->Get(inc_mix_name);
	h1_partpt = (TH1F*)temp3D->ProjectionY("h1_part_pt");
	MakeDphiFrom3D(temp3D,temp3D_mix);

	TH2F* temp2D = (TH2F*)infile->Get(dec_name);
	TH2F* temp2D_mix = (TH2F*)infile->Get(dec_mix_name);
	cout << "Prjecting 1D FG histogram from: " << name << endl;			
	MakeDphiFrom2D(temp2D,temp2D_mix);

}

void MakeWeightedJFs::MakeDphiFrom3D(TH3F* temp3D,TH3F* temp3D_mix)
{
	ostringstream bin;
	string name, name_mix;
	for(int it = 0; it < NTRIGBIN; it++){
		for(int ih = 0; ih < NPARTBIN; ih++){

			bin.str("");
			bin << "_p"<<it<<"_h"<<ih;
			name = "h1_dphi" + bin.str();
			name_mix = "h1_dphi_mix" + bin.str();

			if( XiBinning ) MakeDphiProjection(temp3D,dphi_1d[it][ih],trig_pt[it], trig_pt[it+1], part_pt[ih+1], part_pt[ih],name.c_str());
			else MakeDphiProjection(temp3D,dphi_1d[it][ih],trig_pt[it], trig_pt[it+1], part_pt[ih], part_pt[ih+1],name.c_str());

			if( XiBinning)  MakeDphiProjection(temp3D_mix,dphi_1d_mix[it][ih], trig_pt[it], trig_pt[it+1], part_pt[ih+1], part_pt[ih],name_mix.c_str());
			else  MakeDphiProjection(temp3D_mix,trig_pt[it], dphi_1d_mix[it][ih], trig_pt[it+1], part_pt[ih], part_pt[ih+1],name_mix.c_str());

			SubtractBackground(dphi_1d[it][ih], dphi_1d_mix[it][ih], corr[it][ih]);
		}
	}
}

void MakeWeightedJFs::MakeDphiFrom2D(TH2F* temp2D, TH2F* temp2D_mix)
{
	ostringstream bin;
	string name, name_mix;

	for(int it = 0; it < NTRIGBIN; it++){
		for(int ih = 0; ih < NPARTBIN; ih++){

			bin.str("");
			bin << "_p"<<it<<"_h"<<ih;
			name = "h1_dphi" + bin.str();
			name_mix = "h1_dphi_mix" + bin.str();

			if( XiBinning ) Make2DDphiProjection(temp2D,dphi_1d[it][ih], part_pt[ih+1], part_pt[ih],name.c_str());
			else Make2DDphiProjection(temp2D,dphi_1d[it][ih], part_pt[ih], part_pt[ih+1],name.c_str());

			if( XiBinning) Make2DDphiProjection(temp2D_mix,dphi_1d_mix[it][ih], part_pt[ih+1], part_pt[ih],name_mix.c_str());
			else Make2DDphiProjection(temp2D_mix,dphi_1d_mix[it][ih], part_pt[ih], part_pt[ih+1],name_mix.c_str());
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
	int ybinhi = proj_y->FindBin(ymax);
  	//cout<<"xbinlo = "<<xbinlo<<"; xbinhi = "<<xbinhi<<endl;
  	//cout<<"ybinlo = "<<ybinlo<<"; ybinhi = "<<ybinhi<<endl;
	//h1 = new TH1F(*(TH1F*)h3->Project3D("z"));
	string pz = hname + "_pz";
	h1 = new TH1F(*(TH1F*)h3->ProjectionZ(pz.c_str(),xbinlo,xbinhi-1,ybinlo,ybinhi-1));
	h1->SetName(hname.c_str());
	//cout << "checking projection: phi=0 -> " << h1->GetBinContent(1) << endl;
}

void MakeWeightedJFs::Make2DDphiProjection(TH2F* h3, TH1F*& h1,double ymin, double ymax, string hname)
{
	TH1F* proj_y = (TH1F*)h3->ProjectionY("py");
	int ybinlo = proj_y->FindBin(ymin);
	int ybinhi = proj_y->FindBin(ymax);
	string pz = hname + "_px";
	h1 = new TH1F(*(TH1F*)h3->ProjectionX(pz.c_str(),ybinlo,ybinhi-1));
	h1->SetName(hname.c_str());
	//cout << "checking projection: phi=0 -> " << h1->GetBinContent(1) << endl;
}

void MakeWeightedJFs::SubtractBackground(TH1F* foreground, TH1F*& signal)
{
	double norm = GetZYAMNorm(foreground);
	TF1* bgFunc = new TF1("bgFunc","[0]",0.0,PI);
	bgFunc->SetParameter(0,norm);
	foreground->Add(bgFunc,-1.0);
}

double MakeWeightedJFs::GetZYAMNorm(TH1F* dphi)
{
	dphi->SetAxisRange(1.0,1.6,"X");
	int bin = dphi->GetMinimumBin();
	dphi->SetAxisRange(0.0,TMath::Pi(),"X");
	int lbin = bin-2;//CFinc->FindBin(1.1);
	int hbin = lbin+2;
	float lphi = dphi->GetBinCenter(lbin);
	float hphi = dphi->GetBinCenter(hbin);
	double norm = dphi->Integral(lbin,hbin);
	norm = norm/((double)(hbin-lbin+1));
	cout << "ZYAM norm = " << dphi->Integral(lbin,hbin) << "/(" << hphi << " - " << lphi << ") = " << norm << endl;

	return norm;
}