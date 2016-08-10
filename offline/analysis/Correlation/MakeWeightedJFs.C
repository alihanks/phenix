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
#include <TGraphErrors.h>

using namespace std;

MakeWeightedJFs::MakeWeightedJFs(const string fin, const string fout)
{
	infile = new TFile(fin.c_str());
	outfile = new TFile(fout.c_str(),"recreate");
}

void MakeWeightedJFs::GetHistos(int type, int data_type)
{
	isdAu = data_type;
	for( int ic = 0; ic < NCENT; ic++ )
	{
		cout << "Getting " << prefix << " histograms for centrality bin " << ic << endl;
		Get1DOutputHistos(type,ic);
	}
	cout << "Merging " << prefix << " centrality/pt binned histograms" << endl;
	GetMergedHistos(type);
}

void MakeWeightedJFs::Get1DOutputHistos(int type, int cbin)//0-inc,pi0; 1-dec
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
	name = "h1_trig_pt_" + prefix;
	trigpt_combined->SetName(name.c_str());
	trigpt_combined->Reset();
	cout << "Merging " << prefix << " 1D dphi centrality histograms" << endl;
	for( int ic = 0; ic < NCENT; ic++ ) {
		trigpt_combined->Add(h1_trigpt[ic]);
		for( int it = 0; it < NTRIGBIN; it++ ) {
			if( ic==0 ) {
				h1_partpt_comb[it] = new TH1F(*h1_partpt[ic][it]);
				bin.str("");
				bin << prefix << "_p" << it;
				name = "h1_part_pt_" + bin.str();
				h1_partpt_comb[it]->SetName(name.c_str());
			}
			else h1_partpt_comb[it]->Add(h1_partpt[ic][it]);
			for( int ih = 0; ih < NPARTBIN; ih++ ) {
				if( ic==0 ) {
					jf_comb[it][ih] = new TH1F(*corr[ic][it][ih]);
					bin.str("");
					bin << prefix << "_p" << it << "_h" << ih;
					name = "JF_" + bin.str();
					jf_comb[it][ih]->SetName(name.c_str());
					jf_comb[it][ih]->Scale(ntrig[it]);

					jf_comb_sys[it][ih] = new TH1F(*corr_sys[ic][it][ih]);
					bin.str("");
					bin << prefix << "_p" << it << "_h" << ih;
					name = "JF_sys_" + bin.str();
					jf_comb_sys[it][ih]->SetName(name.c_str());
					jf_comb_sys[it][ih]->Scale(ntrig[it]);
				}
				else {
					jf_comb[it][ih]->Add(corr[ic][it][ih],ntrig[it]);
					jf_comb_sys[it][ih]->Add(corr_sys[ic][it][ih],ntrig[it]);
				}
			}
		}
	}
	outfile->cd();
	double ntrigs_comb[2] = {0,0};
	TH1F* jf_pt_comb[2][NPARTBIN];
	TH1F* jf_pt_comb_sys[2][NPARTBIN];
	TH1F* h1_partpt_tot[2];
	for( int it = 0; it < 2; it++ ) {
		bin.str("");
		bin << prefix << "_" << it;
		name = "h1_part_pt_" + bin.str();
		h1_partpt_tot[it] = new TH1F(*h1_partpt_comb[0]);
		h1_partpt_tot[it]->SetName(name.c_str());
		h1_partpt_tot[it]->Reset();
		for( int ih = 0; ih < NPARTBIN; ih++ ) {
			bin.str("");
			bin << prefix << "_" << it << "_h" << ih;
			name = "JF_" + bin.str();
			jf_pt_comb[it][ih] = new TH1F(*jf_comb[0][ih]);
			jf_pt_comb[it][ih]->SetName(name.c_str());
			jf_pt_comb[it][ih]->Reset();

			bin.str("");
			bin << prefix << "_" << it << "_h" << ih;
			name = "JF_sys_" + bin.str();
			jf_pt_comb_sys[it][ih] = new TH1F(*jf_comb_sys[0][ih]);
			jf_pt_comb_sys[it][ih]->SetName(name.c_str());
			jf_pt_comb_sys[it][ih]->Reset();
		}
	}

	cout << "Subtracting background for " << prefix << " merged centrality histograms" << endl;
	for( int it = 0; it < NTRIGBIN; it++ ) {
		double ntrig_tot = GetNTrigs(type,it,trigpt_combined);
		if( it==0 || it==1 ) {
			ntrigs_comb[0] += ntrig_tot;
			h1_partpt_tot[0]->Add(h1_partpt_comb[it]);
		}
		if( it==2 || it==3 ) {
			ntrigs_comb[1] += ntrig_tot;
			h1_partpt_tot[1]->Add(h1_partpt_comb[it]);
		}

		h1_partpt_comb[it]->Write();
		for( int ih = 0; ih < NPARTBIN; ih++ ) {
			if( it==0 || it==1 ) jf_pt_comb[0][ih]->Add(jf_comb[it][ih]);
			if( it==2 || it==3 ) jf_pt_comb[1][ih]->Add(jf_comb[it][ih]);
			jf_comb[it][ih]->Scale(1/ntrig_tot);
			jf_comb[it][ih]->Write();

			jf_comb_sys[it][ih]->Scale(1/ntrig_tot);
			jf_comb_sys[it][ih]->Write();
		}
	}

	cout << "Subtracting background for " << prefix << " merged pt histograms" << endl;
	for( int it = 0; it < 2; it++ ){
		h1_partpt_tot[it]->Write();
		for( int ih = 0; ih < NPARTBIN; ih++ ) {
			jf_pt_comb[it][ih]->Scale(1/ntrigs_comb[it]);
			jf_pt_comb[it][ih]->Write();

			jf_pt_comb_sys[it][ih]->Scale(1/ntrigs_comb[it]);
			jf_pt_comb_sys[it][ih]->Write();
		}
	}
	trigpt_combined->Write();
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
	double zyam_low[NPARTBIN] = {1.1,1.3,1.1,1.3,1.3,1.3,1.3,1.3};
	double zyam_high[NPARTBIN] = {1.4,1.5,1.4,1.5,1.5,1.5,1.5,1.5};

	for(int it = 0; it < NTRIGBIN; it++){
		bin.str("");
		bin << "_p" << it << "_c" << cbin;
		name = "h1_part_pt_" + prefix + bin.str();
		TH1F* h1_trigpt = (TH1F*)temp3D->ProjectionX();
		int lbin = h1_trigpt->FindBin(trig_pt[it]);
		int hbin = h1_trigpt->FindBin(trig_pt[it+1]);
		h1_partpt[cbin][it] = (TH1F*)temp3D->ProjectionY(name.c_str(),lbin,hbin-1);
		h1_partpt[cbin][it]->Write();
		ntrig[it] = GetNTrigs(0,it,trigpt);

		for(int ih = 0; ih < NPARTBIN; ih++){
			bin.str("");
			bin << "_c" << cbin << "_p"<<it<<"_h"<<ih;
			name = "h1_dphi_" + prefix + bin.str();
			name_mix = "h1_dphi_mix_" + prefix + bin.str();

			MakeDphiProjection(temp3D,dphi_1d[cbin][it][ih],trig_pt[it], trig_pt[it+1], part_pt[ih], part_pt[ih+1],name);
			MakeDphiProjection(temp3D_mix, dphi_1d_mix[cbin][it][ih], trig_pt[it], trig_pt[it+1], part_pt[ih], part_pt[ih+1], name_mix);
			dphi_1d[cbin][it][ih]->Write();
			dphi_1d_mix[cbin][it][ih]->Scale(1.0/Nmix);
			dphi_1d_mix[cbin][it][ih]->Write();

			MakeJetFunction(isdAu, 0, dphi_1d[cbin][it][ih], dphi_1d_mix[cbin][it][ih], corr[cbin][it][ih], ntrig[it], it, ih, cbin, zyam_low[ih], zyam_high[ih]);
			MakeJetFunction(isdAu, 0, dphi_1d[cbin][it][ih], dphi_1d_mix[cbin][it][ih], corr_sys[cbin][it][ih], ntrig[it], it, ih, cbin, 1.0, 1.4);
			outfile->cd();
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
		ntrig[it] = GetNTrigs(1,it,trigpt);

		for(int ih = 0; ih < NPARTBIN; ih++){
			bin.str("");
			bin << "_c" << cbin << "_p"<<it<<"_h"<<ih;
			name = "h1_dphi_" + prefix + bin.str();
			name_mix = "h1_dphi_mix_" + prefix + bin.str();

			Make2DDphiProjection(temp2D,dphi_1d[cbin][it][ih], part_pt[ih], part_pt[ih+1],name);
			Make2DDphiProjection(temp2D_mix,dphi_1d_mix[cbin][it][ih], part_pt[ih], part_pt[ih+1],name_mix);
			dphi_1d[cbin][it][ih]->Write();
			dphi_1d_mix[cbin][it][ih]->Scale(1.0/Nmix);
			dphi_1d_mix[cbin][it][ih]->Write();

			MakeJetFunction(isdAu, 1, dphi_1d[cbin][it][ih], dphi_1d_mix[cbin][it][ih], corr[cbin][it][ih], ntrig[it], it, ih, cbin, 1.1, 1.3);
			MakeJetFunction(isdAu, 1, dphi_1d[cbin][it][ih], dphi_1d_mix[cbin][it][ih], corr_sys[cbin][it][ih], ntrig[it], it, ih, cbin, 1.0, 1.4);
			outfile->cd();
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
		if(bin<4)
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

void MakeWeightedJFs::MakeJetFunction(int isdAu, int type, TH1F* dphi, TH1F* dphi_mix, TH1F*& correlation, double ntrigs, int it, int ih, int cbin, float lphi, float hphi)
{
	ostringstream name;
	name << "JF_" << prefix << "_c" << cbin << "_p" << it << "_h" << ih; 
	if(isdAu)
	  SubtractBackground(dphi, correlation, name.str(), lphi, hphi);
	else{
	  float xi = 0.;
	  float xierr = 0.;
	  GetXi(type, it, ih, cbin, xi, xierr);
	  SubtractBackground(dphi, dphi_mix, xi, correlation, name.str());
	}
	float cutoffcorr = 1.;
	if(!type) GetCutOffCorr(it);
	correlation->Scale(1.0/cutoffcorr);
	float binwidth = correlation->GetBinWidth(1);
	correlation->Scale(1/binwidth);
	correlation->Scale(1/ntrigs);
}

void MakeWeightedJFs::SubtractBackground(TH1F* foreground, TH1F*& signal, string name, float lphi, float hphi)
{
	double norm = GetZYAMNorm(foreground, lphi, hphi);
	signal = new TH1F(*foreground);
	signal->SetName(name.c_str());
	TF1* bgFunc = new TF1("bgFunc","[0]",0.0,PI);
	bgFunc->SetParameter(0,norm);
	signal->Add(bgFunc,-1.0);
}

double MakeWeightedJFs::GetZYAMNorm(TH1F* dphi, float lphi, float hphi)
{
	int lbin = dphi->FindBin(lphi);
	int hbin = dphi->FindBin(hphi);
	double norm = 0; int count = 0;
	for( int ib = lbin; ib <= hbin; ib++ ) {
		if( dphi->GetBinContent(ib) ) norm += dphi->GetBinContent(ib);
		if( dphi->GetBinContent(ib) ) count++;
	}
	norm = norm/((double)count);

	return norm;
}

//for AuAu, use absolute normalization to determine bg level
 void MakeWeightedJFs::SubtractBackground(TH1F* foreground, TH1F* background, float norm, TH1F*& signal, string name)
{
  signal = new TH1F(*foreground);
  signal->SetName(name.c_str());
  signal->Add(background, -1.0*norm);
}

void MakeWeightedJFs::GetXi(int type, int trigptbin, int partptbin, int centbin, float & xi, float & xierr)
{
  cout<<" get xi correction "<<endl;

  char xifilename[1000];
  //testing what is used in filltime method
  if(trigptbin>3){
    if (!type) sprintf(xifilename,"/direct/phenix+u/workarea/mjuszkie/run7_QM09/my_xi/xi_corr_inc_%i_5.root",trigptbin-1);
    else sprintf(xifilename,"/direct/phenix+u/workarea/mjuszkie/run7_QM09/my_xi/xi_corr_pi0_%i_5.root",trigptbin-1);
  }
  else{
    if (!type) sprintf(xifilename,"/direct/phenix+u/workarea/mjuszkie/run7_QM09/my_xi/xi_corr_pi0_%i_5.root",trigptbin);
    else sprintf(xifilename,"/direct/phenix+u/workarea/mjuszkie/run7_QM09/my_xi/xi_corr_pi0_%i_5.root",trigptbin);
  }

  cout << "xi file: " << xifilename <<endl;
  TFile *fxicorr;
  fxicorr=new TFile(xifilename);

  TGraphErrors *gxicorr=(TGraphErrors*)fxicorr->Get("XICORRLARGEBIN");

  double xivalue, dum;
  gxicorr->GetPoint(centbin,dum,xivalue);

  xi = xivalue;
  if(xi==0) cout << "warning xi equals zero!!!!!" <<endl;
  xierr = gxicorr->GetErrorY(centbin)/xi;

  cout << "XXXXXXXXXXXXXXXXXXXXX xi corr is    ====>"  << xi  <<  "<==== " << xierr << endl;
}

float MakeWeightedJFs::GetCutOffCorr(int trig_bin)//values from view_xi_dir_filltime.C
{
  float cutoff_corr[4]={0.99964,0.997983,0.98095,0.961687};
  return cutoff_corr[trig_bin];
}

