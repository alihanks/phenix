void meson_pi0_normalization_auau(double c, double p0, double a, double b, double n, int particle_flag){
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
  int   nbins  = 1000;
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
  double omass, expected;
  if ( particle_flag==1 ) {
    omass  = 0.54730;
    expected = 0.45;
  }
  if ( particle_flag==2 ) {
    omass  = 0.7711;
    expected = 1.0;
  }
  if ( particle_flag==3 ) {
    omass  = 0.78257;
    expected = 1.0;
  }
  if ( particle_flag==4 ) {
    omass  = 0.95778;
    expected = 0.25;
  }
  if ( particle_flag==5 ) {
    omass  = 1.019456;
    expected = 0.4;
  }

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

    // minimum bias
    weight = 2.0*pi*pt*c/pow(exp(-a*pt-b*pt*pt)+pt/p0,n);

    ptpion->AddBinContent(ibin,weight);
    intpion = intpion + weight*binwidth;

    double mt = sqrt(mtother*mtother-pimass*pimass);

    weight = 2.0*pi*pt*c/pow(exp(-a*mt-b*mt*mt)+mt/p0,n);

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
  //ptother->Draw("same");

  TLatex *text_pion0 = new TLatex(3.2,1.5,"#pi^{0}");
  text_pion0->SetTextColor(1);
  text_pion0->Draw();
  TLatex *text_pionm = new TLatex(3.2,1.1,"#eta m_{t} scaled");
  text_pionm->SetTextColor(4);
  text_pionm->Draw();
  TLatex *text_pionp = new TLatex(3.2,0.7,"#eta hydro");
  text_pionp->SetTextColor(2);
  //text_pionp->Draw();

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
  ratiomt->SetLineWidth(2);
  ratiomt->SetLineColor(4);
  ratiomt->Draw();
  //ratio->Draw("same");

  TLatex *text_pionm = new TLatex(3.2,0.42,"#eta m_{t} scaled");
  text_pionm->SetTextColor(4);
  // text_pionm->Draw();
  TLatex *text_pionp = new TLatex(3.2,0.35,"#eta hydro");
  text_pionp->SetTextColor(2);
  // text_pionp->Draw();

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

