// type = 0 for pt bins, 1 for xi, 2 for zt
void make_jet_function(int type, const char* fin = "dAu_merged.root", const char* fout = "test_corr.root")
{
  gSystem->Load("libCorrelationPlots.so");

  ostringstream trig_name, dphi_name, dphi_mix_name;
  ostringstream suffix;
  suffix << "_fold";

  vector trig_bins;
  trig_bins.push_back(5.0); trig_bins.push_back(7.0); trig_bins.push_back(9.0); trig_bins.push_back(12.0); trig_bins.push_back(15.0);
  vector part_bins;
  if( type == 0 ) {
    part_bins.push_back(0.5); part_bins.push_back(1.0); part_bins.push_back(2.0); part_bins.push_back(3.0);
    part_bins.push_back(5.0); part_bins.push_back(7.0); part_bins.push_back(0.0); part_bins.push_back(0.0);
  }
  if( type == 1 ) {
    part_bins.push_back(2.8); part_bins.push_back(2.4); part_bins.push_back(2.0); part_bins.push_back(1.6);
    part_bins.push_back(1.2); part_bins.push_back(0.8); part_bins.push_back(0.4); part_bins.push_back(0.0);
  }

  MakeWeightedJFs* jet_functions = new MakeWeightedJFs(fin,fout);
  jet_functions->SetTriggerBinning(trig_bins);
  jet_functions->SetPartnerBinning(part_bins);
  jet_functions->XiBinning = type%2;

  trig_name << "h1_trig_pt_inc_iso";
  if( type == 0 ) dphi_name << "h3_dphi";
  if( type == 1 ) dphi_name << "h3_ptxidphi";
  dphi_mix_name << dphi_name.str() << "_mix";
  dphi_name << suffix.str();  dphi_mix_name << suffix.str();

  jet_functions->SetTriggerName(trig_name.str());
  jet_functions->SetDphiNames(dphi_name.str(),dphi_mix_name.str());
  jet_functions->GetCentralityHistos(0);

  trig_name.str(""); dphi_name.str(""); dphi_mix_name.str("");
	trig_name << "h1_trig_pt_pi0_iso";
  if( type == 0 ) dphi_name << "h3_dphi_pi0";
  if( type == 1 ) dphi_name << "h3_ptxidphi_pi0";
	dphi_mix_name << dphi_name.str() << "_mix";
  dphi_name << suffix.str();  dphi_mix_name << suffix.str();

  jet_functions->SetTriggerName(trig_name.str());
  jet_functions->SetDphiNames(dphi_name.str(),dphi_mix_name.str());
  jet_functions->GetCentralityHistos(0);

  trig_name.str(""); dphi_name.str(""); dphi_mix_name.str("");
	trig_name << "h1_trig_pt_dec_iso";
  if( type == 0 ) dphi_name << "h2_dphi_dec";
  if( type == 1 ) dphi_name << "h2_dphixi_dec";
	dphi_mix_name << dphi_name.str() << "_mix";
  dphi_name << suffix.str();  dphi_mix_name << suffix.str();

  jet_functions->SetTriggerName(trig_name.str());
  jet_functions->SetDphiNames(dphi_name.str(),dphi_mix_name.str());
  jet_functions->GetCentralityHistos(1);

  cout<<"finish making CFs for type "<<type<<endl;
}

