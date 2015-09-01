void make_cfs(int type = 0, int isiso = 0, const char* fin = "dAu_merged.root", const char* fout = "test_jf.root")
{
  gSystem->Load("libCorrelationPlots.so");

  ostringstream trig_name;
  if( type==0 ) {
  	trig_name << "h1_trig_pt_inc";
  }
  if( type==1 ) {
  	trig_name << "h1_trig_pt_pi0";
  }
  if( type==2 ) {
  	trig_name << "h1_trig_pt_dec";
  }

  if( isiso ) {
  	trig_name << "_iso";
  }

  MakeCombinedHistos* combine = new MakeCombinedHistos(fin,fout,trig_name.str());
  cout<<"finished combining JFs for type "<<type<<endl;
}