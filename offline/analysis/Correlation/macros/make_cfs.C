// type = 0(inc), 1(pi0), 2(dec)
void make_cfs(int type = 0, int isiso = 0, const char* fin = "dAu_merged.root", const char* fout = "test_corr.root", int ispertrigger = 0)
{
  gSystem->Load("libCorrelationPlots.so");

  int isfold = 2;
  int dofold = 0;
  //int ispertrigger = 2; // 2 - pp/dAu ZYAM method
  //int dosub = 0;
  string rgamma_file = "/phenix/u/workarea/ahanks/devel/offline/AnalysisTrain/combinesimple/wrk/run10/Rgamma/Rgamma_final_test.root";

  MakeCFs(type,isfold,dofold,isiso,ispertrigger,fin,fout);
  cout<<"finish making CFs for type "<<type<<endl;
}

