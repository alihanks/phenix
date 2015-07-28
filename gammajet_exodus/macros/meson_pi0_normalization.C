void meson_pi0_normalization(){
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);

  c1->SetFillColor(10);
  c1->SetGrid();
  c1->SetLogy();
  c1->SetTicks();

  float pt, mtpion, mtother;
  float weight = 1.0;
  int   nbins  = 5000;
  float ptmin  = 0.0;
  float ptmax  = 5.0;
  float binwidth = (ptmax-ptmin)/(double)nbins;
  TH1F * ptpion = new TH1F("ptpion","ptpion",nbins,ptmin,ptmax);
  TH1F * ptother = new TH1F("ptother","ptother",nbins,ptmin,ptmax);
  TH1F * ptothermt = new TH1F("ptothermt","ptothermt",nbins,ptmin,ptmax);
  TH1F * ratio = new TH1F("ratio","ratio",nbins,ptmin,ptmax);
  TH1F * ratiomt = new TH1F("ratiomt","ratiomt",nbins,ptmin,ptmax);

  TMath math1;

  double pi     = acos(-1.0);
  double pimass = 0.1349766;
  double omass  = 0.54730;
  double expected = 0.45;
  // double omass  = 0.7711;
  // double expected = 1.0;
  // double omass  = 0.78257;
  // double expected = 1.0;
  // double omass  = 0.95778;
  // double expected = 0.25;
  // double omass  = 1.019456;
  // double expected = 0.4;

  double t_fo = 0.122;
  double beta = 0.70;
  double norm = 1.0;
 
  double intpion    = 0.0;
  double intother   = 0.0;
  double intothermt = 0.0;

  for (int ibin=1; ibin<=nbins; ibin++ )
  {
    pt = ptmin+(double)(ibin-1)*binwidth+binwidth/2.0;
    mtpion  = sqrt(pt*pt+pimass*pimass);
    mtother = sqrt(pt*pt+omass*omass);

    // minimum bias:
    // double trans      = 1./(1.+exp(5.*pt-17.5));
    // double pihagedorn = 3.53e+07/pow(2.28+pt,13.7);
    // double pipower    = 33.9/pow(pt,8.16);
    // double otherhagedorn = 3.53e+07 / 
    //   pow(2.28+sqrt(mtother*mtother-pimass*pimass),13.7);
    // double otherpower = 33.9/pow(sqrt(mtother*mtother-pimass*pimass),8.16);
    // 0 - 10 % central
    // double trans      = 1./(1.+exp(5.*pt-17.5));
    // double pihagedorn = 1.41e+13/pow(3.48+pt,18.7);
    // double pipower    = 137.3/pow(pt,8.44);
    // double otherhagedorn = 1.41e+13 / 
    //   pow(3.48+sqrt(mtother*mtother-pimass*pimass),18.7);
    // double otherpower = 137.3/pow(sqrt(mtother*mtother-pimass*pimass),8.44);
    // 10 - 20 % central
    // double trans      = 1./(1.+exp(5.*pt-17.5));
    // double pihagedorn = 2.176e+09/pow(2.651+pt,15.170);
    // double pipower    = 92.46/pow(pt,8.326);
    // double otherhagedorn = 2.176e+09 / 
    //   pow(2.651+sqrt(mtother*mtother-pimass*pimass),15.170);
    // double otherpower = 92.46/pow(sqrt(mtother*mtother-pimass*pimass),8.326);    // 20 - 40 % central
    // double trans      = 1./(1.+exp(5.*pt-17.5));
    // double pihagedorn = 2.99e+06/pow(2.00+pt,12.45);
    // double pipower    = 51.28/pow(pt,8.21);
    // double otherhagedorn = 2.99e+06 / 
    //   pow(2.00+sqrt(mtother*mtother-pimass*pimass),12.45);
    // double otherpower = 51.28/pow(sqrt(mtother*mtother-pimass*pimass),8.21);
    // 40 - 60 % central
    // double trans      = 1./(1.+exp(5.*pt-17.5));
    // double pihagedorn = 3.00e+04/pow(1.56+pt,10.76);
    // double pipower    = 15.75/pow(pt,7.97);
    // double otherhagedorn = 3.00e+04 / 
    //   pow(1.56+sqrt(mtother*mtother-pimass*pimass),10.76);
    // double otherpower = 15.75/pow(sqrt(mtother*mtother-pimass*pimass),7.97);
    // 60 - 92 % central
    // double trans      = 1./(1.+exp(5.*pt-17.5));
    // double pihagedorn = 8.73e+02/pow(1.29+pt,9.94);
    // double pipower    = 5.28/pow(pt,8.28);
    // double otherhagedorn = 8.73e+02 / 
    //   pow(1.29+sqrt(mtother*mtother-pimass*pimass),9.94);
    // double otherpower = 5.28/pow(sqrt(mtother*mtother-pimass*pimass),8.28);
    // if ( pt<1.5 ) weight = 2.0 * pi * pt * pihagedorn;
    // if ( pt>=1.5 && pt<=5.5 ) 
    //   weight = 2.0 * pi * pt * (trans*pihagedorn + (1.-trans)*pipower);
    // if ( pt>5.5 ) weight = 2.0*pi*pt*pipower;

    // pp (Xinhua's Run-2 parameterization)
    // weight = 2.0*pi*pt * (485.0/42.) / pow(exp(-0.5187*pt) + pt/0.6631,8.163);
    // pp (Run-2 hard parameterization)
    weight = 2.0*pi*pt * (450.8/42.) / pow(exp(-0.5553*pt) + pt/0.6467,8.046);
    // pp (Run-3 parameterization)
    // weight = 2.0*pi*pt * (413.5/42.2) / pow(exp(-0.2614*pt) + pt/0.8144,8.807);
    // pp (average Run-2/3 parameterization)
    // weight = 2.0*pi*pt * (546.2/42.2) / pow(exp(-0.2664*pt) + pt/0.7703,8.718);
    //    weight = weight*(1+0.6362/exp(7.736*pt));

    ptpion->AddBinContent(ibin,weight);
    intpion = intpion + weight*binwidth;

    // if ( pt<1.5 ) weight = 2.0 * pi * pt * otherhagedorn;
    // if ( pt>=1.5 && pt<=5.5 ) 
    //   weight = 2.0 * pi * pt * (trans*otherhagedorn + (1.-trans)*otherpower);
    // if ( pt>5.5 ) weight = 2.0*pi*pt*otherpower;

    weight = 2.0*pi*pt * (450.8/42.2) / 
      pow(exp(-0.5553*sqrt(mtother*mtother-pimass*pimass)) + 
	  sqrt(mtother*mtother-pimass*pimass)/0.6467,8.046);

    ptothermt->AddBinContent(ibin,weight);
    intothermt = intothermt + weight*binwidth;

    weight = norm*InitializePtHydro(omass,t_fo,beta,pt);
    ptother->AddBinContent(ibin,weight);
    intother = intother + weight*binwidth;
  }

  cout << "Integral pion    (dN/dy): " << intpion << endl;
  cout << "Integral other   (dN/dy): " << intother << endl;
  cout << "Integral othermt (dN/dy): " << intothermt << endl;

  norm = intpion/intothermt;
  ptothermt->Scale(norm);
 
  for (int ibin=1; ibin<=nbins; ibin++ )
  {
    weight = ptothermt->GetBinContent(ibin)/ptpion->GetBinContent(ibin);
    ratiomt->AddBinContent(ibin,weight);
  }

  norm = expected/ratiomt->GetMaximum();
  ratiomt->Scale(norm);
  ptothermt->Scale(norm);
  cout << "meson/pion ratio should be set to: " << norm << endl;

  norm = expected*intothermt/intother;
  ptother->Scale(norm);
 
  for (int ibin=1; ibin<=nbins; ibin++ )
  {
    weight = ptother->GetBinContent(ibin)/ptpion->GetBinContent(ibin);
    ratio->AddBinContent(ibin,weight);
  }

  ptpion->SetXTitle("p_{t} [GeV/c]");
  ptpion->SetYTitle("(1/2#pi p_{t})dN/dp_{t}dy [(c/GeV)^{2}]");
  ptpion->SetLineWidth(5);
  ptpion->SetLineColor(1);
  ptother->SetLineWidth(5);
  ptother->SetLineColor(2);
  ptothermt->SetLineWidth(5);
  ptothermt->SetLineColor(4);
  ptpion->Draw("");
  ptothermt->Draw("same");
  ptother->Draw("same");

  TLatex *text_pion0 = new TLatex(3.2,1.5,"#pi^{0}");
  text_pion0->SetTextColor(1);
  text_pion0->Draw();
  TLatex *text_pionm = new TLatex(3.2,1.1,"#eta m_{t} scaled");
  text_pionm->SetTextColor(4);
  text_pionm->Draw();
  TLatex *text_pionp = new TLatex(3.2,0.7,"#eta hydro");
  text_pionp->SetTextColor(2);
  text_pionp->Draw();

   c2 = new TCanvas("c2","A Simple Graph with error bars",200,10,700,500);

   c2->SetFillColor(10);
   c2->SetGrid();
   c2->SetLogy(0);
   c2->SetTicks();
 
   c2->cd();

  ratiomt->SetXTitle("p_{t} [GeV/c]");
  ratiomt->SetYTitle("#eta / #pi^{0}");
  ratio->SetLineWidth(5);
  ratio->SetLineColor(2);
  ratiomt->SetLineWidth(5);
  ratiomt->SetLineColor(4);
  ratiomt->Draw();
  ratio->Draw("same");

  TLatex *text_pionm = new TLatex(3.2,0.42,"#eta m_{t} scaled");
  text_pionm->SetTextColor(4);
  text_pionm->Draw();
  TLatex *text_pionp = new TLatex(3.2,0.35,"#eta hydro");
  text_pionp->SetTextColor(2);
  text_pionp->Draw();

  return;
}

float InitializePtHydro(double mass, double t_fo, double beta, float pt) 
{
  float sum    = 0.0;
  float weight = 1.0;
  double profx, valprofx, mt;
  double besi, besk, besiarg, beskarg, arg;

  TMath math;

  for (int profbin=1; profbin<101; profbin++)
  {
    profx    = ((double)(profbin-1) + 0.5) / 100.0;
    mt       = sqrt(mass*mass+pt*pt);
    arg      = math.ATanH(beta*profx);
    besiarg  = pt*sinh(arg)/t_fo;;
    beskarg  = mt*cosh(arg)/t_fo;;
    besi     = math.BesselI0(besiarg); 
    besk     = math.BesselK1(beskarg); 
    valprofx = mt*profx*besi*besk;

    sum = sum + valprofx;
  }

  weight = pt*sum*0.01; 

  return weight;

}

