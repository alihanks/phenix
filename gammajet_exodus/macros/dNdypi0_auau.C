void dNdypi0_auau(double c, double p0, double a, double b, double n)
{
  float pt, mt;
  float weight = 1.0;
  int nbins = 10000;
  float ptmin = 0.0;
  float ptmax = 10.;
  float binwidth = (ptmax-ptmin)/(double)nbins;
  TH1F * pthist = new TH1F("pthist","pthist",nbins,ptmin,ptmax);

  double pi    = acos(-1.0);

  double integral = 0.0;

  for (int ibin=1; ibin<=nbins; ibin++ )
  {
    pt = ptmin+(double)(ibin-1)*binwidth+binwidth/2.0;

    // minimum bias
    weight = 2.0*pi*pt*c/pow(exp(-a*pt-b*pt*pt)+pt/p0,n);

    pthist->AddBinContent(ibin,weight);
    integral = integral + weight*binwidth;
  }

  cout << "Integral (dN/dy): " << integral << endl;
 
}

