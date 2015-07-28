{
  float pt, mt;
  float weight = 1.0;
  int nbins = 1000;
  float ptmin = 0.0;
  float ptmax = 10.;
  float binwidth = (ptmax-ptmin)/(double)nbins;
  TH1F * pthist = new TH1F("pthist","pthist",nbins,ptmin,ptmax);

  double pi    = acos(-1.0);

  double integral = 0.0;

  for (int ibin=1; ibin<=nbins; ibin++ )
  {
    pt = ptmin+(double)(ibin-1)*binwidth+binwidth/2.0;

    TF1 *trans = new TF1("trans","(1./(1.+exp(5.*x-17.5)))",0.,15.);
    // minimum bias:
    //    TF1 *hagedorn = new TF1("hagedorn","3.53e+07/pow(2.28+x,13.7)",0.,15.);
    //    TF1 *power = new TF1("power","33.9/pow(x,8.16)",0.,15.);
    // 0 - 10 % central
    // TF1 *hagedorn = new TF1("hagedorn","1.41e+13/pow(3.48+x,18.7)",0.,15.);
    // TF1 *power = new TF1("power","137.3/pow(x,8.44)",0.,15.);
    // 10 - 20 % central
    // TF1 *hagedorn = new TF1("hagedorn","2.176e+09/pow(2.651+x,15.170)",0.,15.);
    // TF1 *power = new TF1("power","92.46/pow(x,8.326)",0.,15.);
    // 20 - 40 % central
    // TF1 *hagedorn = new TF1("hagedorn","2.99e+06/pow(2.00+x,12.45)",0.,15.);
    // TF1 *power = new TF1("power","51.28/pow(x,8.21)",0.,15.);
    // 40 - 60 % central
    // TF1 *hagedorn = new TF1("hagedorn","3.00e+04/pow(1.56+x,10.76)",0.,15.);
    // TF1 *power = new TF1("power","15.75/pow(x,7.97)",0.,15.);
    // 60 - 92 % central
    // TF1 *hagedorn = new TF1("hagedorn","8.73e+02/pow(1.29+x,9.94)",0.,15.);
    // TF1 *power = new TF1("power","5.28/pow(x,8.28)",0.,15.);
    // if ( pt<1.5 ) weight = 2.0 * pi * pt * hagedorn->Eval(pt);
    // if ( pt>=1.5 && pt<=5.5 ) 
    //   weight = 2.0 * pi * pt * (trans->Eval(pt)*hagedorn->Eval(pt) 
    //				+ (1.-trans->Eval(pt))*power->Eval(pt));
    // if ( pt>5.5 ) weight = 2.0*pi*pt*power->Eval(pt);

    // pp (Xinhua's Run-2 parameterization)
    // weight = 2.0*pi*pt * (485.0/42.) / pow(exp(-0.5187*pt) + pt/0.6631,8.163);
    // pp (Run-2 hard parameterization)
    weight = 2.0*pi*pt * (450.8/42.) / pow(exp(-0.5553*pt) + pt/0.6467,8.046);
    // pp (average Run-2/3 parameterization)
    // weight = 2.0*pi*pt * (546.2/42.2) / pow(exp(-0.2664*pt) + pt/0.7703,8.718);
    // pp (Run-3 parameterization)
    // weight = 2.0*pi*pt * (413.5/42.2) / pow(exp(-0.2614*pt) + pt/0.8144,8.807);
    // weight = weight*(1+0.6362/exp(7.736*pt));

    pthist->AddBinContent(ibin,weight);
    integral = integral + weight*binwidth;
  }

  cout << "Integral (dN/dy): " << integral << endl;
 
}

