#include <iostream>

using namespace std;

void make_cfs()
{
  gSystem->Load("libCorrelation.so");

  int type = 1;
  int isfold = 2;
  int dofold = 0;
  int ispertrigger = 0;
  //int dosub = 0;
  string rgamma_file = "/phenix/u/workarea/ahanks/devel/offline/AnalysisTrain/combinesimple/wrk/run10/Rgamma/Rgamma_final_test.root";
  // string fin = "/direct/phenix+hhj/hge/taxi/Run11AuAu200MinBias/3645/data/Merged_3645.root";
  // string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4448/data/correlation_taxi4448.root";
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4509/data/correlation_taxi4509.root";
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4528/data/correlation_taxi4528.root";
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4560/data/correlation_taxi4560.root";
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4569/data/correlation_taxi4569.root";
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4594/data/correlation_taxi4594.root";
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4587/data/correlation_taxi4587.root";//2sigma
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4615/data/correlation_taxi4615.root";//2sigma
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4627/data/correlation_taxi4627.root";//2sigma, master recalibrator, missing higher pt hadrons!
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4651/data/correlation_taxi4651.root";//2sigma, master recalibrator, with extra warm tower cleaned up, run10 hadeff applied in fg
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4662/data/correlation_taxi4662.root";//2sigma, master recalibrator, with extra warm tower cleaned up, run10 had eff applied fg+bg
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4680/data/correlation_taxi4680.root";//pc3||emc 3sigma, run10 had eff taken out. 
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4686/data/correlation_taxi4686.root";//pc3||emc 2sigma
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4694/data/correlation_taxi4694.root";//pc3||emc 2sigma for pt 0.5-3, 1.5sigma for pt 3-5, 1sigma for pt 5-7.
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4702/data/correlation_taxi4702.root";//pc3||emc 2sigma for pt 0.5-3, 1.5sigma for pt 3-5, 1sigma for pt 5-7.add decay trigger counters, some update on recalibrator?
  // string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4734/data/correlation_taxi4734.root";//pc3 2sigma cut || emc 2sigma cut where there is no pc3 acceptance
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4755/data/correlation_taxi4755.root";//add pi0 asymmetry cut
  //string fin = "/phenix/hhj/hge/taxi/Run11AuAu200MinBias/4759/data/correlation_taxi4759.root";//fix decay counter
  //string fin = "/direct/phenix+hhj2/hge/taxi/Run11AuAu200MinBias/5288/data/correlation_taxi5288.root";//fix a bug in event selection, |zvtx|<30cm, centrality only 0-40%
  string fin = "/direct/phenix+hhj2/hge/taxi/Run11AuAu200MinBias/5286/data/correlation_taxi5286.root";
  string fout;
  if(!ispertrigger){
    if(isfold){
      if(type == 0) fout = "makecfs_inc_fold_taxi5286.root";
      if(type == 1) fout = "makecfs_pi0_fold_taxi5286.root";
      if(type == 2) fout = "makecfs_dec_fold_taxi5286.root";
    }
    else{
      if(type == 0) fout = "makecfs_inc.root";
      if(type == 1) fout = "makecfs_pi0.root";
      if(type == 2) fout = "makecfs_dec.root";
    }
  }
  else{
    if(isfold){
      if(type == 0) fout = "makejfs_inc_fold_taxi5286.root";
      if(type == 1) fout = "makejfs_pi0_fold_taxi5286.root";
      if(type == 2) fout = "makejfs_dec_fold_taxi5286.root";
    }
    else{
      if(type == 0) fout = "makejfs_inc.root";
      if(type == 1) fout = "makejfs_pi0.root";
      if(type == 2) fout = "makejfs_dec.root";
    }
  }

  MakeCFs(type,isfold,dofold,ispertrigger,fin.c_str(),fout.c_str());
  cout<<"finish making CFs for type "<<type<<endl;
}

