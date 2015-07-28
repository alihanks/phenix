void test_random(double lambda=1.){

  double costheta;
  TH1D *hist = new TH1D("hist","hist",20,-1.,1.);
  TF1 *fitfunc = new TF1("fitfunc","[0]*(1.+[1]*x*x)",-1.,1.);
  fitfunc->SetParameter(0,1.);
  fitfunc->SetParameter(1,0.);

  for (int ievent=1; ievent<=1000000; ievent++){    
    do{
      costheta = (2.0*gRandom->Rndm())-1.;
    }
    while ( (1.0+lambda*costheta*costheta)<(2.0*gRandom->Rndm()) );
    hist->Fill(costheta);
  }

  hist->Fit("fitfunc");

  return;
}
