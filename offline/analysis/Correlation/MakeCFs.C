#include "MakeCFs.h"
#include "MakeJFs.h"

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

MakeCFs::MakeCFs(const string fin, const string fout, int isxi)
{
	infile = new TFile(fin.c_str());
	outfile = new TFile(fout.c_str(),"recreate");
	XiBinning = isxi;
	SetTrigPtBinning();
	SetPartPtBinning(isxi);
}

void MakeCFs::Run(int type, int ispertrigger)
{
	cout<<"start"<<endl;
	TH1F* h1_trigpt[NCENTBIN];
	TH1F* h1_partpt[NCENTBIN];

	for(int ic=0; ic<NCENTBIN; ic++){
		bin.str("");
		bin << "_c" << ic;
		name = trig_name + bin.str();
		cout << "getting trig histo: " << name << endl;
		h1_trigpt[ic] = new TH1F(*(TH1F*)infile->Get(name.c_str()));
	}

	outfile->cd();
	string dphi_title = ";#Delta#phi[rad]      ";
	cout<<"using folded histos."<<endl;
	for(int ic=0; ic < NCENTBIN; ic++){
		bin.str("");
		bin << "_c" << ic;

		name = dphi_name + bin.str();
		name_mix = dphi_mix_name + bin.str();

		if(type == 0 || type == 1){
			temp3D = new TH3F(*(TH3F*)infile->Get(name.c_str()));
			cout << "Prjecting 1D FG histogram from: " << name << endl;			
			name = "h1_part_pt" + bin.str();
			h1_partpt[ic] = new TH1F(*(TH1F*)temp3D->ProjectionY(name.c_str()));
			temp3D_mix = new TH3F(*(TH3F*)infile->Get(name_mix.c_str()));
			cout << "Prjecting 1D BG histogram from: " << name_mix << endl;
			//cout<<"get temp3D: "<<name.c_str()<<endl;
			for(int ippt = 0; ippt < NTRIGBIN; ippt++){
				for(int ihpt = 0; ihpt < NPARTBIN; ihpt++){
					dphi_3d[ic][ippt][ihpt] = new TH3F(*(TH3F*)temp3D);
					dphi_3d_mix[ic][ippt][ihpt] = new TH3F(*(TH3F*)temp3D_mix);

					bin.str("");
					bin << ic <<"_p"<<ippt<<"_h"<<ihpt;
					name = "h1_dphi_c" + bin.str();
					name_mix = "h1_dphi_mix_c" + bin.str();
					if( XiBinning ) dphi_1d[ic][ippt][ihpt] = MakeDphiProjection(dphi_3d[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt+1], part_pt[ihpt],name.c_str());
					else dphi_1d[ic][ippt][ihpt] = MakeDphiProjection(dphi_3d[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1],name.c_str());
					SetHisto(dphi_1d[ic][ippt][ihpt],dphi_title,1);
					dphi_1d[ic][ippt][ihpt]->SetName(name.c_str());

					if( XiBinning) dphi_1d_mix[ic][ippt][ihpt] = MakeDphiProjection(dphi_3d_mix[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt+1], part_pt[ihpt],name_mix.c_str());
					else dphi_1d_mix[ic][ippt][ihpt] = MakeDphiProjection(dphi_3d_mix[ic][ippt][ihpt],trig_pt[ippt], trig_pt[ippt+1], part_pt[ihpt], part_pt[ihpt+1],name_mix.c_str());

					SetHisto(dphi_1d_mix[ic][ippt][ihpt],dphi_title,2);
					dphi_1d_mix[ic][ippt][ihpt]->SetName(name_mix.c_str());
					if(ispertrigger) {
						dphi_1d_mix[ic][ippt][ihpt]->Scale(1/500.0);
					}
				}
			}
		}
		if(type == 2){
			for(int ippt=0; ippt<NTRIGBIN; ippt++){
				bin.str("");
				if(ippt<3) bin <<"_p"<<ippt<<"_c"<<ic;	   
				else bin<<"_p"<<ippt+1<<"_c"<<ic;

				name = dphi_name + bin.str();
				name_mix = dphi_mix_name + bin.str();

				temp2D = new TH2F(*(TH2F*)infile->Get(name.c_str()));
				bin.str("");
				bin<<"_c"<<ic;
				name = "h1_part_pt" + bin.str();
				h1_partpt[ic] = new TH1F(*(TH1F*)temp2D->ProjectionY(name.c_str()));
				temp2D_mix = new TH2F(*(TH2F*)infile->Get(name_mix.c_str()));
				for(int ihpt=0; ihpt < NPARTBIN; ihpt++){  
					dphi_2d[ic][ippt][ihpt] = new TH2F(*(TH2F*)temp2D);
					dphi_2d_mix[ic][ippt][ihpt] = new TH2F(*(TH2F*)temp2D_mix);
					bin.str("");
					bin << ic <<"_p"<<ippt<<"_h"<<ihpt;
					name = "h1_dphi_c" + bin.str();
					name_mix = "h1_dphi_mix_c" + bin.str();

					int ymin = dphi_2d[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ihpt]);
					int ymax = dphi_2d[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ihpt+1]);
					dphi_1d[ic][ippt][ihpt] = new TH1F (*(TH1F*)dphi_2d[ic][ippt][ihpt]->ProjectionX(name.c_str(),ymin,ymax));
					SetHisto(dphi_1d[ic][ippt][ihpt],dphi_title,1);
					ymin = dphi_2d_mix[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ihpt]);
					ymax = dphi_2d_mix[ic][ippt][ihpt]->GetYaxis()->FindBin(part_pt[ihpt+1]);
					dphi_1d_mix[ic][ippt][ihpt] = new TH1F (*(TH1F*)dphi_2d_mix[ic][ippt][ihpt]->ProjectionX(name_mix.c_str(),ymin,ymax));
					SetHisto(dphi_1d_mix[ic][ippt][ihpt],dphi_title,2);
					if(ispertrigger) {
						dphi_1d_mix[ic][ippt][ihpt]->Scale(1/500.0);
						meanpart[ic][ippt][ihpt] = dphi_1d_mix[ic][ippt][ihpt]->Integral();
					}
				}
			}
		}
	}

	for(int ic=0; ic<NCENTBIN; ic++){
		for(int ippt=0; ippt<NTRIGBIN; ippt++){
			can_dphi_name.str("");
			can_dphi_name << "can_dphi_c"<<ic<<"_p"<<ippt;
			can_dphi[ic][ippt] = new TCanvas(can_dphi_name.str().c_str(), can_dphi_name.str().c_str());
			can_dphi[ic][ippt]->Divide(3,2,0.001,0.001);
			can_corr_name.str("");
			can_corr_name << "can_corr_c"<<ic<<"_p"<<ippt;
			can_corr[ic][ippt] = new TCanvas(can_corr_name.str().c_str(), can_corr_name.str().c_str());
			can_corr[ic][ippt]->Divide(3,2,0.001,0.001);
			can_corr_name.str("");
			can_corr_name << "can_jet_c"<<ic<<"_p"<<ippt;
			can_jet[ic][ippt] = new TCanvas(can_corr_name.str().c_str(), can_corr_name.str().c_str());
			can_jet[ic][ippt]->Divide(3,2,0.001,0.001);
			can_corr_name.str("");
			can_corr_name << "can_abs_c"<<ic<<"_p"<<ippt;
			can_abs[ic][ippt] = new TCanvas(can_corr_name.str().c_str(), can_corr_name.str().c_str());
			can_abs[ic][ippt]->Divide(3,2,0.001,0.001);

			if(type == 0 || type == 1) {
				num_trigger[ic][ippt] = GetNTriggers(h1_trigpt[ic], trig_pt[ippt], trig_pt[ippt+1]);
				num_trigger_mix[ic][ippt] = GetNTriggers(h1_trigpt[ic], trig_pt[ippt], trig_pt[ippt+1]);
			}
			if(type == 2) {
				if(ippt<3){
					num_trigger[ic][ippt] = h1_trigpt[ic]->GetBinContent(ippt+1);
					num_trigger_mix[ic][ippt] = h1_trigpt[ic]->GetBinContent(ippt+1);
				}
				else{
					num_trigger[ic][ippt] = h1_trigpt[ic]->GetBinContent(ippt+2);
					num_trigger_mix[ic][ippt] = h1_trigpt[ic]->GetBinContent(ippt+2);
				}
			}

			cout<<"num_trigger = "<<num_trigger[ic][ippt]<<endl;
			//cout<<"num_trigger_mix = "<<num_trigger_mix[ic][ippt]<<endl;

			for(int ihpt=0; ihpt<NPARTBIN; ihpt++){

				//double R_bg = dphi_1d_mix[ic][ippt][ihpt]->Integral("width")/PI;
				//if( R_bg==0 ) R_bg = 1.0;
				//dphi_1d_mix[ic][ippt][ihpt]->Scale(1/R_bg);

				//double R_fg = dphi_1d[ic][ippt][ihpt]->Integral("width")/PI;
				//if( R_fg==0 ) R_fg = 1.0;
				//if(!ispertrigger) dphi_1d[ic][ippt][ihpt]->Scale(1/R_fg);

				TVirtualPad* pad = can_dphi[ic][ippt]->cd(ihpt+1);
				SetPad(pad);
				double min = 0.75*dphi_1d_mix[ic][ippt][ihpt]->GetMinimum();
				double max = 1.1*dphi_1d[ic][ippt][ihpt]->GetMaximum();
				dphi_1d[ic][ippt][ihpt]->SetAxisRange(min,max,"Y");
				dphi_1d[ic][ippt][ihpt]->Draw();
				dphi_1d_mix[ic][ippt][ihpt]->Draw("same");

				legend_name.str("");
				legend_name<<trig_pt[ippt]<<"-"<<trig_pt[ippt+1]<<" #times "<<part_pt[ihpt]<<"-"<<part_pt[ihpt+1]<<" GeV/c";
				TLegend *l1 = new TLegend(0.45,0.7,0.75,0.9,legend_name.str().c_str(),"brNDC");
				l1->SetFillColor(0);
				l1->SetBorderSize(0);
				l1->AddEntry(dphi_1d[ic][ippt][ihpt],"real","lpf");
				l1->AddEntry(dphi_1d_mix[ic][ippt][ihpt],"mixed","lpf");
				l1->SetTextSize(0.05);
				l1->Draw("same");
				dphi_1d[ic][ippt][ihpt]->Write();
				dphi_1d_mix[ic][ippt][ihpt]->Write();

			    //*****************************************************	  
				corr_name.str("");
				corr_name << "CF_c" << ic << "_p"<<ippt <<"_h"<< ihpt; 
				corr[ic][ippt][ihpt] = new TH1F(*(TH1F*)dphi_1d[ic][ippt][ihpt]);
				corr[ic][ippt][ihpt]->Divide(dphi_1d_mix[ic][ippt][ihpt]);
				SetHisto(corr[ic][ippt][ihpt],dphi_title,1);
				corr[ic][ippt][ihpt]->SetName(corr_name.str().c_str());

				if(ispertrigger) corr[ic][ippt][ihpt]->Scale(1/num_trigger[ic][ippt]);

				else{
					double r = corr[ic][ippt][ihpt]->Integral("width")/PI;
					corr[ic][ippt][ihpt]->Scale(1/r);
				}

				TLatex *la = new TLatex(0.45, 0.75, legend_name.str().c_str());
				la->SetNDC();
				pad = can_corr[ic][ippt]->cd(ihpt+1);
				SetPad(pad);
				corr[ic][ippt][ihpt]->Draw();
				corr[ic][ippt][ihpt]->Write();
				la->Draw("same");
			}
			can_dphi[ic][ippt]->Write();
			can_corr[ic][ippt]->Write();
		}
		h1_trigpt[ic]->Write();
		h1_partpt[ic]->Write();
	}

	if(ispertrigger){
	    //make JFs.
		for(int ic=0; ic<NCENTBIN; ic++){
			for(int itrig=0; itrig<NTRIGBIN; itrig++){
				for(int ipart=0; ipart<NPARTBIN; ipart++){
					cout << "Making JFs for (c,t,p) bin (" << ic << ", " << itrig << ", " << ipart << ")" << endl;
					MakeJFs(type,ic,itrig,ipart,corr[ic][itrig][ipart],meanpart[ic][itrig][ipart],num_trigger_mix[ic][itrig],infile,"v2_inputs.root",5,ispertrigger,flow[ic][itrig][ipart],jet[ic][itrig][ipart]);
					bin.str("");
					bin << "_c" << ic <<"_p"<<itrig<<"_h"<<ipart;
					name = "JF" + bin.str();
					//cout<<"name: "<<name.c_str()<<endl;
					//jet[ic][itrig][ipart]->Sumw2();
					//if(type==0) jet[ic][itrig][ipart]->GetYaxis()->SetTitle("Arb.Unit");
					if(type==0) jet[ic][itrig][ipart]->GetYaxis()->SetTitle("Y_{inc}");
					if(type==1) jet[ic][itrig][ipart]->GetYaxis()->SetTitle("Y_{#pi^{0}}");
					if(type==2) jet[ic][itrig][ipart]->GetYaxis()->SetTitle("Y_{dec}");
					MakeJFs(type,ic,itrig,ipart,corr[ic][itrig][ipart],meanpart[ic][itrig][ipart],num_trigger_mix[ic][itrig],infile,"v2_inputs.root",5,3,flow[ic][itrig][ipart],jet_err[ic][itrig][ipart]);

					double binwidth = jet[ic][itrig][ipart]->GetBinWidth(1);
					//cout<<"jet func binwidth: "<<binwidth<<endl;
					jet[ic][itrig][ipart]->Scale(1/binwidth);
					jet_err[ic][itrig][ipart]->Scale(1/binwidth);
					SetHisto(jet[ic][itrig][ipart],dphi_title,1);
					jet[ic][itrig][ipart]->SetName(name.c_str());
					name = "JFerr" + bin.str();
					jet_err[ic][itrig][ipart]->SetName(name.c_str());
					TVirtualPad* pad = can_jet[ic][itrig]->cd(ipart+1);
					SetPad(pad);
					jet[ic][itrig][ipart]->Draw();
					jet[ic][itrig][ipart]->Write();
					jet_err[ic][itrig][ipart]->Write();

					legend_name.str("");
					if( XiBinning ) legend_name<<trig_pt[itrig]<<"-"<<trig_pt[itrig+1]<<" #times "<<part_pt[ipart+1]<<"-"<<part_pt[ipart]<<" GeV/c";
					else legend_name<<trig_pt[itrig]<<"-"<<trig_pt[itrig+1]<<" #times "<<part_pt[ipart]<<"-"<<part_pt[ipart+1]<<" GeV/c";
					TLatex *la = new TLatex(0.25, 0.75, legend_name.str().c_str());
					la->SetNDC();
					la->Draw("same");

					name = "flow" + bin.str();
					flow[ic][itrig][ipart]->SetName(name.c_str());
					pad = can_abs[ic][itrig]->cd(ipart+1);
					SetPad(pad);
					corr[ic][itrig][ipart]->Draw();
					flow[ic][itrig][ipart]->SetLineColor(2);
					flow[ic][itrig][ipart]->Draw("same");
					la->Draw("same");
					outfile->cd();
					//jet[ic][itrig][ipart]->Write();
					//cout<<"jetfunc written out."<<endl;
				}
				can_abs[ic][itrig]->Write();
				can_jet[ic][itrig]->Write();
			}
		}
	}
	outfile->Close();

}

void MakeCFs::SetPtRange(TH3F* h3, double x_pt_min, double x_pt_max, double y_pt_min, double y_pt_max)
{
	h3->GetXaxis()->SetRangeUser(x_pt_min, x_pt_max);
	h3->GetYaxis()->SetRangeUser(y_pt_min, y_pt_max);
}

void MakeCFs::SetTrigPtBinning()
{
	trig_pt[0]=5.;
	trig_pt[1]=7.;
	trig_pt[2]=9.;
	trig_pt[3]=12.;
	trig_pt[4]=15.;
}

void MakeCFs::SetPartPtBinning(int isxi)
{
	if( isxi ) {
		part_pt[7]=0.;
		part_pt[6]=0.4;
		part_pt[5]=0.8;
		part_pt[4]=1.2;
		part_pt[3]=1.6;
		part_pt[2]=2.0;
		part_pt[1]=2.4;
		part_pt[0]=2.8;
	}
	else {
		part_pt[0] = 0.0;
		part_pt[1] = 0.5;
		part_pt[2] = 1.0;
		part_pt[3] = 2.0;
		part_pt[4] = 3.0;
		part_pt[5] = 5.0;
		part_pt[6] = 7.0;
	}
}

void MakeCFs::MakeDphiProjection(TH3F* h3, TH1F*& h1,string hname)
{
	h1 = new TH1F(*(TH1F*)h3->Project3D("z"));
	h1->SetNameTitle(hname.c_str(),hname.c_str());
  // cout<<"h1 name: "<<hname.c_str()<<endl;
}

void MakeCFs::MakeDphiProjection(TH3F* h3, TH1F*& h1,double xmin, double xmax, double ymin, double ymax, string hname)
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
	//cout << "checking projection: phi=0 -> " << h1->GetBinContent(1) << endl;
}

TH1F* MakeCFs::MakeDphiProjection(TH3F* h3, float xmin, float xmax, float ymin, float ymax, string hname)
{
  TH1F* proj_x = (TH1F*)h3->ProjectionX("px");
  TH1F* proj_y = (TH1F*)h3->ProjectionY("py");
  int xbinlo = proj_x->FindBin(xmin);
  int xbinhi = proj_x->FindBin(xmax);
  int ybinlo = proj_y->FindBin(ymin);
  int ybinhi = proj_y->FindBin(ymax);
  string pz = hname + "_pz";
  cout << "Projecting from xbin " << xbinlo << " - " << xbinhi << " and ybin " << ybinlo << " - " << ybinhi << endl;
  TH1F* proj_hist = (TH1F*)h3->ProjectionZ(pz.c_str(),xbinlo,xbinhi-1,ybinlo,ybinhi-1);
  return (TH1F*)proj_hist;
}

void MakeCFs::FoldDphiDist(TH1F* h1, TH1F*& h1_fold, string hname_fold)
{
	int Nbins = h1->GetNbinsX();
	h1_fold = new TH1F(hname_fold.c_str(),"",Nbins/2,0.0,PI);
	for(int ibin = 1; ibin <= Nbins/4; ibin++){
		double added_bin_left = h1->GetBinContent(ibin)+h1->GetBinContent(Nbins/2-ibin+1);
		double added_bin_right = h1->GetBinContent(Nbins/2+ibin)+h1->GetBinContent(Nbins-ibin+1);    
    // double lerr1 = h1->GetBinError(ibin);
    // double lerr2 = h1->GetBinError(Nbins/2-ibin+1);
    // double rerr1 = h1->GetBinError(Nbins/2+ibin);
    // double rerr2 = h1->GetBinError(Nbins-ibin+1);

		h1_fold->SetBinContent(16-ibin,added_bin_left);
		h1_fold->SetBinContent(ibin+15,added_bin_right);
    // h1_fold->SetBinError(16-ibin,sqrt(lerr1*lerr1+lerr2*lerr2));
    // h1_fold->SetBinError(ibin+15,sqrt(rerr1*rerr1+rerr2*rerr2));
	}
}

void MakeCFs::SetHisto(TH1F* h1, string title, int color)
{
	h1->SetName("");
	h1->SetTitle(title.c_str());
	h1->SetTitleSize(0.05,"X");
	h1->SetTitleOffset(0.95,"X");
	h1->SetLabelSize(0.05,"X");
	h1->SetMarkerColor(color);
	h1->SetLineColor(color);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(0.8);
}

void MakeCFs::SetPad(TVirtualPad* pad)
{
	pad->SetTopMargin(0.02);
	pad->SetRightMargin(0.02);
}

double MakeCFs::GetNTriggers(TH1F* trigpt, double trigptmin, double trigptmax)
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

double MakeCFs::GetHadronEff(TH1F* hadron_pt, int ipart)
{
	TF1* fhadeff = new TF1("fhadeff","(x<=3)*(0.396-0.337*exp(-1.50*x))+(x>3&&x<=5)*(0.400+-2.13e-13*exp(5.36*x))+(x>5)*(0.396-0.337*exp(-1.50*x))",0.5,10.0);
	hadron_pt->SetAxisRange(part_pt[ipart],part_pt[ipart+1],"X");
	double meanpt = hadron_pt->GetMean();
	cout<<"mean pt for pt within ["<<part_pt[ipart]<<", "<<part_pt[ipart+1]<<"]"<<" is "<<meanpt<<endl;
	double eff = fhadeff->Eval(meanpt);
  //  cout<<"eff = "<<eff<<endl;
	return eff;
}

double MakeCFs::GetHadronEff_v2(TH1F* hadron_pt, int ipart)
{
	TF1* fhadeff = new TF1("fhadeff","(x<=3)*(0.396-0.337*exp(-1.50*x))+(x>3&&x<=5)*(0.400+-2.13e-13*exp(5.36*x))+(x>5)*(0.396-0.337*exp(-1.50*x))",0.5,10.0);
	int bin1 = hadron_pt->FindBin(part_pt[ipart]);
	int bin2 = hadron_pt->FindBin(part_pt[ipart+1]); 
	double binwidth = hadron_pt->GetBinWidth(1);

	cout<<"bin1 = "<<bin1<<"; bin2 = "<<bin2<<"; bin width = "<<binwidth<<endl;

	double hadeff = 0.0;
	double meanpt = 0.0;
	for(int ibin=bin1; ibin<bin2; ibin++){
		cout<<"ibin = "<<ibin<<endl;
		double pt = hadron_pt->GetBinCenter(ibin);
		cout<<"pt = "<<pt<<endl;
		double eff = fhadeff->Eval(pt);
		cout<<"eff = "<<eff<<" at pt = "<<pt<<endl;
		hadeff+=pt*eff;
		meanpt+=pt;
	}
	hadeff/=meanpt;
  //  cout<<"hadeff = "<<hadeff<<endl;
	return hadeff;
}
