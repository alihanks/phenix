void make_dirs()
{
  gSystem->Load("libCorrelation.so");

  string rgamma_file = "/phenix/u/workarea/ahanks/devel/offline/AnalysisTrain/combinesimple/wrk/run10/Rgamma/Rgamma_final_test.root";
  string incfile = "makejfs_inc_fold_taxi4759.root";
  string decfile = "makejfs_dec_fold_taxi4759.root";
  string outfile = "makedirs_taxi4759.root";
  MakeDir(rgamma_file.c_str(),incfile.c_str(),decfile.c_str(),outfile.c_str());
}
