// type = 0 for pt bins, 1 for xi, 2 for zt
// data_type = 0 for AuAu and = 1 for dAu
void make_jet_function(int type, int data_type, int do_iso, const char* fin = "dAu_merged.root", const char* fout = "test_corr.root")
{
  gSystem->Load("libCorrelationPlots.so");

  ostringstream trig_name, dphi_name, dphi_mix_name;
  ostringstream suffix; ostringstream iso;
  suffix << "_fold";
  iso.str("");
  if( do_iso ) iso << "_iso";

  vector<double> trig_bins;
  trig_bins.push_back(5.0); trig_bins.push_back(7.0); trig_bins.push_back(9.0); trig_bins.push_back(12.0); trig_bins.push_back(15.0);
  vector<double> part_bins;
  if( type == 0 ) {
    part_bins.push_back(0.61); part_bins.push_back(1.01); part_bins.push_back(2.01); part_bins.push_back(3.01);
    part_bins.push_back(5.01); part_bins.push_back(7.01); part_bins.push_back(0.0); part_bins.push_back(0.0);
  }
  if( type == 1 ) {
    part_bins.push_back(2.42); part_bins.push_back(2.01); part_bins.push_back(1.61);
    part_bins.push_back(1.22); part_bins.push_back(0.81); part_bins.push_back(0.41); part_bins.push_back(0.0); part_bins.push_back(0.0);
  }

  MakeWeightedJFs* jet_functions = new MakeWeightedJFs(fin,fout);
  jet_functions->SetTriggerBinning(trig_bins);
  jet_functions->SetPartnerBinning(part_bins);
  jet_functions->XiBinning = type%2;
  jet_functions->Nmix = 500;

  trig_name << "h1_trig_pt_inc" << iso.str();
  if( type == 0 ) dphi_name << "h3_dphi" << iso.str();
  if( type == 1 ) dphi_name << "h3_ptxidphi" << iso.str();
  dphi_mix_name << dphi_name.str() << "_mix";
  dphi_name << suffix.str();  dphi_mix_name << suffix.str();

  jet_functions->SetTriggerName(trig_name.str());
  jet_functions->SetDphiNames(dphi_name.str(),dphi_mix_name.str());
  jet_functions->SetTypeName("inc");
  jet_functions->GetHistos(0,data_type);

  trig_name.str(""); dphi_name.str(""); dphi_mix_name.str("");
	trig_name << "h1_trig_pt_pi0" << iso.str();
  if( type == 0 ) dphi_name << "h3_dphi_pi0" << iso.str();
  if( type == 1 ) dphi_name << "h3_ptxidphi_pi0" << iso.str();
	dphi_mix_name << dphi_name.str() << "_mix";
  dphi_name << suffix.str();  dphi_mix_name << suffix.str();

  jet_functions->SetTriggerName(trig_name.str());
  jet_functions->SetDphiNames(dphi_name.str(),dphi_mix_name.str());
  jet_functions->SetTypeName("pi0");
  jet_functions->GetHistos(0,data_type);

  trig_name.str(""); dphi_name.str(""); dphi_mix_name.str("");
	trig_name << "h1_trig_pt_dec" << iso.str();
  if( type == 0 ) dphi_name << "h2_dphi_dec" << iso.str();
  if( type == 1 ) dphi_name << "h2_dphixi_dec" << iso.str();
	dphi_mix_name << dphi_name.str() << "_mix";
  dphi_name << suffix.str();  dphi_mix_name << suffix.str();

  jet_functions->SetTriggerName(trig_name.str());
  jet_functions->SetDphiNames(dphi_name.str(),dphi_mix_name.str());
  jet_functions->SetTypeName("dec");
  jet_functions->GetHistos(1,data_type);

  cout<<"finish making CFs for type "<<type<<endl;
}

