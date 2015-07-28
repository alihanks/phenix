{
  float pt, mt;
  float weight = 1.0;
  int nbins = 100000;
  float ptmin = 0.0;
  float ptmax = 10.0;
  float binwidth = (ptmax-ptmin)/(double)nbins;
  TH1F * pthist = new TH1F("pthist","pthist",nbins,ptmin,ptmax);

  double pi    = acos(-1.0);
  double mass  = 0.135;
  double n_pow = 386.;
  double power = 9.99;
  double p0    = 1.219;

  double integral = 0.0;

  for (int ibin=1; ibin<=nbins; ibin++ )
  {
    pt = ptmin+(double)(ibin-1)*binwidth+binwidth/2.0;
    mt = sqrt(pt*pt+0.135*0.135);

    weight = 2.0*pi*pt*(n_pow/exp(power*log(1.0+pt/p0)));

    pthist->AddBinContent(ibin,weight);
    integral = integral + weight*binwidth;
  }

  cout << "Integral (dN/dy): " << integral << endl;
 
}

