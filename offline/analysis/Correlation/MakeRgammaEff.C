
double GetNTriggers(TH1D* trigpt, double trigptmin, double trigptmax)
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

void MakeRgammaEff(const char* Rgamma_input, const char* trigger_input)
{
	// Rg_eff = Rg*alpha_miss/(1-e_niso,dec)(1-e_tag,dec)
	// alpha_miss = N_ntag,iso/N_inc
	// e_niso,dec = 1 - W_ntag,iso/W_ntag
	// e_tag,dec = Rg*(1-N_ntag/N_inc)

	TFile* infile = new TFile(trigger_input);
	TFile* Rgfile = new TFile(Rgamma_input,"recreate");

	const int NTRIGBIN = 4;
	double Rgamma[NTRIGBIN];
    Rgamma[0]=1.18352;
    Rgamma[1]=1.33204;
    Rgamma[2]=1.5254;
    Rgamma[3]=1.20248;
    double Rgamma_eff[NTRIGBIN];

	double trigpt_bins[NTRIGBIN+1] = {5.0,7.0,9.0,12.0,15.0};
	double trigpt[NTRIGBIN] = {5.0,7.0,9.0,12.0};

	const int NCENTBIN = 1;
	TH1D* h1_trigpt_iso[NCENTBIN];
	TH1D* h1_trigpt_all[NCENTBIN];
	TH1D* h1_trigpt[NCENTBIN];
	TH1D* h1_trigpt_dec[NCENTBIN];
	TH1D* h1_trigpt_dec_iso[NCENTBIN];

	TGraphErrors* gr[NCENTBIN];
	TGraphErrors* stat[NCENTBIN];
	TGraphErrors* sys[NCENTBIN];

	std::ostringstream bin;
	for(int ic=0; ic<NCENTBIN; ic++){
		bin.str("");
		bin << "_c" << ic;

		name = "h1_trig_pt_all" + bin.str();
		h1_trigpt_all[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
		name = "h1_trig_pt_inc" + bin.str();
		h1_trigpt[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
		name = "h1_trig_pt_inc_iso" + bin.str();
		h1_trigpt_iso[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
		name = "h1_trig_pt_dec" + bin.str();
		h1_trigpt_dec[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));
		name = "h1_trig_pt_dec_iso" + bin.str();
		h1_trigpt_dec_iso[ic] = new TH1D(*(TH1D*)infile->Get(name.c_str()));

		for(int ip=0; ip<NTRIGBIN; ip++){
			double ntrig_all = GetNTriggers(h1_trigpt_all[ic],trigpt_bins[ip],trigpt_bins[ip+1]);
			double ntrig_ntag = GetNTriggers(h1_trigpt[ic],trigpt_bins[ip],trigpt_bins[ip+1]);
			double ntrig_iso = GetNTriggers(h1_trigpt_iso[ic],trigpt_bins[ip],trigpt_bins[ip+1]);
			double alpha = ntrig_iso/ntrig_all;
			cout << "alpha = " << alpha << endl;
			double e_tag = Rgamma[ip]*(1-ntrig_ntag/ntrig_all);
			cout << "e-tag = " << e_tag << endl;
			double e_niso = 1 - h1_trigpt_dec_iso[ic]->GetBinContent(ip+1)/h1_trigpt_dec[ic]->GetBinContent(ip+1);
			cout << "e-niso = " << e_niso << endl;
			if(ip==0) e_tag = e_tag*1.2;
			Rgamma_eff[ip] = Rgamma[ip]*alpha/((1.0-e_niso)*(1.0-e_tag));
			cout << "Rg = " << Rgamma_eff[ip] << endl;
		}

		Rgfile->cd();
	    name = "gr" + bin.str();
	    gr[ic] = new TGraphErrors(NTRIGBIN,trigpt,0,Rgamma_eff,0);
	    gr[ic]->Write();
	    name = "stat" + bin.str();
	    stat[ic] = new TGraphErrors(NTRIGBIN,trigpt,0,Rgamma_eff,0);
	    stat[ic]->Write();
	    name = "sys" + bin.str();
	    sys[ic] = new TGraphErrors(NTRIGBIN,trigpt,0,Rgamma_eff,0);
	    sys[ic]->Write();
	}

}


