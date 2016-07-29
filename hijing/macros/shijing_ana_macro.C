void shijing_ana_macro(char* input="hijing.dat", char* output="test.root", int events = 0, int cone_size = 3)
{
  gSystem->Load("libfun4allfuncs");
  gSystem->Load("libsimreco");
  gSystem->Load("libphhepmc.so");
  gSystem->Load("libCorrelation");
  cout << "Loading libhijing_analysis from " << gSystem->DynamicPathName("libhijing_analysis") << endl;
  gSystem->Load("libhijing_analysis");
  
  Fun4AllServer *se = Fun4AllServer::instance();
  //se->Verbosity(3);

  hijing_analysis *hijing = new hijing_analysis(output);
  //hijing->Verbosity(1);
  hijing->SetSharkFin("/phenix/hhj/ahanks/gammajet_exodus/sharkfin_run8_miss.root");
  hijing->Rcut = cone_size*0.1;
  hijing->_MinTrigPt = 3.0;
  //cout << "Set cone cut to " << cone_size*0.1 << " -> " << hijing->Rcut << endl;
  hijing->Centrality = -1;
  se->registerSubsystem(hijing);

  Fun4AllInputManager *in1 = new Fun4AllHepMCInputManager("DSTIN");
  //in1->Verbosity(3);
  se->registerInputManager(in1);
  se->fileopen(in1->Name(), input);
  cout << input << endl;

  cout << "running" << endl;
  se->run(events);
  cout << "done running" << endl;
  se->End();
}
