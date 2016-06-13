// type = 0(inc), 1(pi0), 2(dec)
  //int ispertrigger = 2; // 2 - pp/dAu ZYAM method
void make_cfs(int type = 0, int isiso = 0, int atype = 0, const char* fin = "dAu_merged.root", const char* fout = "test_corr.root", int ispertrigger = 0)
{
  gSystem->Load("libCorrelationPlots.so");

  int isfold = 2;
  ostringstream trig_name, dphi_name, dphi_mix_name;
  if( type==0 ) {
  	trig_name << "h1_trig_pt_inc";
  	if( atype == 0 ) dphi_name << "h3_dphi";
        else dphi_name << "h3_ptxidphi";
  	dphi_mix_name << dphi_name.str() << "_mix";
  }
  if( type==1 ) {
  	trig_name << "h1_trig_pt_pi0";
        if( atype == 0 ) dphi_name << "h3_dphi_pi0";
        else dphi_name << "h3_ptxidphi_pi0";
  	dphi_mix_name << dphi_name.str() << "_mix";
  }
  if( type==2 ) {
  	trig_name << "h1_trig_pt_dec";
        if( atype == 0 ) dphi_name << "h2_dphi_dec";
        else dphi_name << "h2_dphixi_dec";
  	dphi_mix_name << dphi_name.str() << "_mix";
  }

  if( isiso==1 ) {
  	trig_name << "_iso";
  	//dphi_name << "_iso";
  	//dphi_mix_name << "_iso";
  }

  if( isfold ) {
  	dphi_name << "_fold";
  	dphi_mix_name << "_fold";
  }

  MakeCFs* cfs = new MakeCFs(fin,fout,atype);
  cfs->SetTriggerName(trig_name.str());
  cfs->SetDphiNames(dphi_name.str(),dphi_mix_name.str());
  cfs->Run(type,ispertrigger);
  cout<<"finish making CFs for type "<<type<<endl;
}

