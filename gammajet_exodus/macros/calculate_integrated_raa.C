void calculate_integrated_raa(int cent = 0, double pt_min = 2.5){
  if (pt_min<2.4) cout << "R_AA integrated from 2 GeV/c" << endl;
  if (pt_min<2.6 && pt_min>2.4) 
    cout << "R_AA integrated from 2.5 GeV/c" << endl;
  if (pt_min>2.5) cout << "R_AA integrated from 3 GeV/c" << endl;
  if (cent == 0){
    // minimum bias
    cout << "minimum bias" << endl;
    TFile *f_noshift = TFile::Open("ralf/1_MIN_graphs.root");  
    TFile *f_shift = TFile::Open("ralf/1_MIN_bin_shifted_non_photonic.root");  
    double embed_rel_error = 0.064;
    double t_ab = 6.14, t_ab_rel = 0.073;
  }
  if (cent == 1){
    cout << "0-10 %" << endl;
    // 0-10 %
    TFile *f_noshift = TFile::Open("ralf/1_C0110_graphs.root");  
    TFile *f_shift = TFile::Open("ralf/1_C0110_bin_shifted_non_photonic.root");  
    double embed_rel_error = 0.067;
    double t_ab = 22.75, t_ab_rel = 0.069;
  }
  if (cent == 2){
    cout << "10-20 %" << endl;
    // 10-20 %
    TFile *f_noshift = TFile::Open("ralf/1_C1120_graphs.root");  
    TFile *f_shift = TFile::Open("ralf/1_C1120_bin_shifted_non_photonic.root");  
    double embed_rel_error = 0.081;
    double t_ab = 14.35, t_ab_rel = 0.070;
  }
  if (cent == 3){
    cout << "20-40 %" << endl;
    // 20-40 %
    TFile *f_noshift = TFile::Open("ralf/1_C2140_graphs.root");  
    TFile *f_shift = TFile::Open("ralf/1_C2140_bin_shifted_non_photonic.root");  
    double embed_rel_error = 0.060;
    double t_ab = 7.07, t_ab_rel = 0.082;
  }
  if (cent == 4){
    cout << "40-60 %" << endl;
    // 40-60 %
    TFile *f_noshift = TFile::Open("ralf/1_C4160_graphs.root");  
    TFile *f_shift = TFile::Open("ralf/1_C4160_bin_shifted_non_photonic.root");  
    double embed_rel_error = 0.037;
    double t_ab = 2.16, t_ab_rel = 0.118;
  }
  if (cent == 5){
    cout << "60-80 %" << endl;
    // 60-80 %
    TFile *f_noshift = TFile::Open("ralf/1_C6180_graphs.root");  
    TFile *f_shift = TFile::Open("ralf/1_C6180_bin_shifted_non_photonic.root");  
    double embed_rel_error = 0.013;
    double t_ab = 0.49, t_ab_rel = 0.29;
  }

  double pi = acos(-1);
  
  TGraphAsymmErrors *np_stat;
  TGraphAsymmErrors *np_sys;

  TGraphAsymmErrors *inc_stat;
  TGraphAsymmErrors *inc_sys;
  TGraphAsymmErrors *sub_stat;
  TGraphAsymmErrors *sub_sys;
  TGraphAsymmErrors *coc_stat;
  TGraphAsymmErrors *coc_sys;
  TGraphAsymmErrors *pp_stat;
  TGraphAsymmErrors *pp_sys;

  np_stat = (TGraphAsymmErrors*)f_shift->Get("non_photonic_stat")->Clone();
  np_sys  = (TGraphAsymmErrors*)f_shift->Get("non_photonic_sys")->Clone();
  f_shift->Close();
  inc_stat = (TGraphAsymmErrors*)f_noshift->Get("inclusive_stat")->Clone();
  inc_sys  = (TGraphAsymmErrors*)f_noshift->Get("inclusive_sys")->Clone();
  coc_stat = (TGraphAsymmErrors*)f_noshift->Get("cocktail_stat")->Clone();
  coc_sys  = (TGraphAsymmErrors*)f_noshift->Get("cocktail_sys")->Clone();
  sub_stat = (TGraphAsymmErrors*)f_noshift->Get("subtracted_stat")->Clone();
  sub_sys  = (TGraphAsymmErrors*)f_noshift->Get("subtracted_sys")->Clone();
  pp_stat = (TGraphAsymmErrors*)f_noshift->Get("sum_plot")->Clone();
  pp_sys  = (TGraphAsymmErrors*)f_noshift->Get("sys_plot")->Clone();
  f_noshift->Close();

  int    nval;
  double *pt;
  nval = np_stat->GetN();
  pt = np_stat->GetX();

  double *np, *npstat, *npsys;
  np = np_stat->GetY();
  npstat = np_stat->GetEYhigh();
  npsys  = np_sys->GetEYhigh();
  double *inc, *incstat, *incsys;
  inc = inc_stat->GetY();
  incstat = inc_stat->GetEYhigh();
  incsys  = inc_sys->GetEYhigh();
  double *coc, *cocstat, *cocsys;
  coc = coc_stat->GetY();
  cocstat = coc_stat->GetEYhigh();
  cocsys  = coc_sys->GetEYhigh();
  double *sub, *substat, *subsys;
  sub = sub_stat->GetY();
  substat = sub_stat->GetEYhigh();
  subsys  = sub_sys->GetEYhigh();
  double *pp, *ppstat, *ppsys;
  pp = pp_stat->GetY();
  ppstat = pp_stat->GetEYhigh();
  ppsys  = pp_sys->GetEYhigh();

  if ( pt_min>2.6 ){
    np[9]  = 0.; npstat[9]  = 0.; npsys[9]  = 0.;
    inc[9] = 0.; incstat[9] = 0.; incsys[9] = 0.;
    coc[9] = 0.; cocstat[9] = 0.; cocsys[9] = 0.;
    sub[9] = 0.; substat[9] = 0.; subsys[9] = 0.;
    pp[10] = 0.; ppstat[10] = 0.; ppsys[10] = 0.;
  }

  if ( pt_min>2.4 ){
    np[8]  = 0.; npstat[8]  = 0.; npsys[8]  = 0.;
    inc[8] = 0.; incstat[8] = 0.; incsys[8] = 0.;
    coc[8] = 0.; cocstat[8] = 0.; cocsys[8] = 0.;
    sub[8] = 0.; substat[8] = 0.; subsys[8] = 0.;
    pp[9]  = 0.; ppstat[9]  = 0.; ppsys[9]  = 0.;
  }

  double npsysCent[12],  npsysNoCent[12];
  double incsysCent[12], incsysNoCent[12];
  double cocsysCent[12], cocsysNoCent[12];
  double subsysCent[12], subsysNoCent[12];

  double int_inc, int_inc_stat, int_inc_sys_cent, int_inc_sys_nocent;
  double int_coc, int_coc_stat, int_coc_sys_cent, int_coc_sys_nocent;
  double int_sub, int_sub_stat, int_sub_sys_cent, int_sub_sys_nocent;
  double int_np, int_np_stat, int_np_sys_cent, int_np_sys_nocent;
  double int_pp, int_pp_stat, int_pp_sys;
  double int_np_shift, int_np_stat_shift, int_np_sys_shift;

  double pt_width[12] = {0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 
			 0.4, 0.5, 0.5, 1., 1.}; 

  int_inc=0.; int_inc_stat=0.; int_inc_sys_cent=0.; int_inc_sys_nocent=0.;
  int_coc=0.; int_coc_stat=0.; int_coc_sys_cent=0.; int_coc_sys_nocent=0.;
  int_sub=0.; int_sub_stat=0.; int_sub_sys_cent=0.; int_sub_sys_nocent=0.;
  int_pp=0.;  int_pp_stat=0.;  int_pp_sys=0.;
  int_np_shift=0.; int_np_stat_shift=0.;  int_np_sys_shift=0.;
  for (int ibin=8; ibin<=11; ibin++){
    npsysCent[ibin] = 0.; npsysNoCent[ibin] = 0.;
    incsysCent[ibin]   = embed_rel_error * inc[ibin];
    incsysNoCent[ibin] = 0.1555 * inc[ibin];
    cocsysCent[ibin]   = 0.14 * coc[ibin];
    cocsysNoCent[ibin] = 0.05 * coc[ibin];
    int_inc += 2.0*pi*pt[ibin]*pt_width[ibin]*inc[ibin];
    int_coc += 2.0*pi*pt[ibin]*pt_width[ibin]*coc[ibin];
    int_sub += 2.0*pi*pt[ibin]*pt_width[ibin]*sub[ibin];
    int_inc_stat += (2.0*pi*pt[ibin]*pt_width[ibin]*incstat[ibin])**2;
    int_coc_stat += (2.0*pi*pt[ibin]*pt_width[ibin]*cocstat[ibin])**2;
    int_sub_stat += (2.0*pi*pt[ibin]*pt_width[ibin]*substat[ibin])**2;
    int_sub_sys_cent += 2.0*pi*pt[ibin]*pt_width[ibin]*subsys[ibin];
    int_pp += 2.0*pi*pt[ibin]*pt_width[ibin]*pp[ibin+1];
    int_pp_stat += (2.0*pi*pt[ibin]*pt_width[ibin]*ppstat[ibin+1])**2;
    int_pp_sys  += 2.0*pi*pt[ibin]*pt_width[ibin]*ppsys[ibin+1];
    int_np_shift += 2.0*pi*pt[ibin]*pt_width[ibin]*np[ibin];
    int_np_stat_shift += (2.0*pi*pt[ibin]*pt_width[ibin]*npstat[ibin])**2;
    int_np_sys_shift  += 2.0*pi*pt[ibin]*pt_width[ibin]*npsys[ibin];
  }
  int_inc_stat = sqrt(int_inc_stat);
  int_coc_stat = sqrt(int_coc_stat);
  int_sub_stat = sqrt(int_sub_stat);
  int_inc_sys_cent   = embed_rel_error * int_inc;
  int_inc_sys_nocent = 0.1555 * int_inc;
  int_coc_sys_cent   = 0.14 * int_coc;
  int_coc_sys_nocent = 0.05 * int_coc;
  int_pp_stat = sqrt(int_pp_stat);
  int_np_stat_shift = sqrt(int_np_stat_shift);

  int_np = int_inc - int_coc; 
  int_np_stat = sqrt(int_inc_stat**2 + int_coc_stat**2);
  int_np_sys_cent = sqrt(int_inc_sys_cent**2 + int_coc_sys_cent**2);
  int_np_sys_nocent = sqrt(int_inc_sys_nocent**2 + int_coc_sys_nocent**2);

  double int_np_sys_tot = sqrt(int_np_sys_cent**2 + int_np_sys_nocent**2);
  double raa;
  double raa_stat_rel, raa_sys_cent_rel, raa_sys_nocent_rel;
  double raa_stat_abs, raa_sys_cent_abs, raa_sys_nocent_abs;

  raa                = int_np_shift/(t_ab*int_pp);
  raa_stat_rel       = int_np_stat_shift/int_np_shift;

  raa_sys_cent_rel   = sqrt((int_np_sys_cent/int_np_sys_tot * 
			     int_np_sys_shift/int_np_shift)**2 + 
			    t_ab_rel**2);

  raa_sys_nocent_rel   = sqrt((int_np_sys_nocent/int_np_sys_tot * 
			       int_np_sys_shift/int_np_shift)**2 + 
			      (int_pp_stat/int_pp)**2 +
			      (int_pp_sys/int_pp)**2);

  cout << "R_AA = " << raa << endl;
  cout << "rel.stat.error " << raa_stat_rel 
       << " abs.stat.error " << raa_stat_rel*raa << endl;
  cout << "rel.sys.error (cent.dep.) " << raa_sys_cent_rel 
       << " abs.sys.error (cent.dep.) " << raa_sys_cent_rel*raa << endl; 
  cout << "rel.sys.error (cent.indep.) " << raa_sys_nocent_rel 
       << " abs.sys.error (cent.indep.) " << raa_sys_nocent_rel*raa << endl;

  return;
}
