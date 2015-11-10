void Run_d_Au_Correlation(const char* outFile = "test.root")
{
  std::cout << "Run_d_Au_Correlation: loading libCorrelation.so" << std::endl;
  int status = gSystem->Load("libCorrelation.so");
  std::cout << "Run_d_Au_Correlation: Loading libCorrelation.so returned status = " << status << std::endl;

  Correlation *co = new Correlation(outFile);
  //co->Verbosity(1);
  co->SetDiagFlag(0);
  co->SetNCentBins(10,4);
  co->SetNMixEvents(500);
  co->SetHadronEffFileName("inputs/run8/chhadron_eff_dAu_C0.root");
  co->SetPi0EffFileName("inputs/run8/dA_run8_trig_eff_0_v0.root","inputs/run8/dA_run8_trig_eff_1_v0.root","inputs/run8/dA_run8_trig_eff_2_v0.root","inputs/run8/dA_run8_trig_eff_3_v0.root");
  co->SetSharkFinFileName("inputs/run8/sharkfin_miss.root");
  co->SetWeightFileNames("inputs/run8/filltimeweighting.root","inputs/run11/v2_inputs.root");

  std::vector<std::string> warnmaps;
  std::vector<double> map_cuts;
  //warnmaps.push_back("Warnmap/d_Au_warnmap.txt"); map_cuts.push_back(999.0);
  warnmaps.push_back("Warnmap/d_Au_warnmap_0to5GeV_cut.txt"); map_cuts.push_back(5.0);
  warnmaps.push_back("Warnmap/d_Au_warnmap_above5GeV_cut.txt"); map_cuts.push_back(999.0);
  co->SetWarnmaps(warnmaps,map_cuts);
  co->SetZVertexCut(30.0);
  co->Rcut = 0.4;

  Fun4AllServer *se = Fun4AllServer::instance();

  Fun4AllSyncManager *insync1 = new Fun4AllSyncManager("INSYNC1");
  se->registerSyncManager(insync1);
  insync1->Repeat();
  Fun4AllInputManager *in = new Fun4AllDstInputManager("MIXCNT","DST","MIXCNT");
  insync1->registerInputManager(in);
  if (!gSystem->AccessPathName("cnt_mb.list"))
  {
    in->AddListFile("cnt_mb.list");
  }
  else
  {
    cout << "cnt_mb.list does not exist" << endl;
    gSystem->Exit(1);
  }
  
  in = new Fun4AllDstInputManager("MIXPWG","DST","MIXPWG");
  insync1->registerInputManager(in);
  if (!gSystem->AccessPathName("pwg_mb.list"))
  {
    in->AddListFile("pwg_mb.list");
  }
  else
  {
    cout << "pwg_mb.list does not exist" << endl;
    gSystem->Exit(1);
  }
  
  se->registerSubsystem(co);

}

void InputData(std::vector<std::string>& inData)
{
  inData.push_back("CNT");
  inData.push_back("PWG");
  inData.push_back("EWG");
  return;
}

