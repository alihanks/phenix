{
TFile *out = new TFile("cocktail_systematic_pp.root","RECREATE");

TH1D *pion_total     = new TH1D("pion_total","pion_total",60,0.,6.);
TH1D *pion_max       = new TH1D("pion_max","pion_max",60,0.,6.);
TH1D *pion_min       = new TH1D("pion_min","pion_min",60,0.,6.);
pion_total->Sumw2();
pion_max->Sumw2();
pion_min->Sumw2();

TH1D *eta_total      = new TH1D("eta_total","eta_total",60,0.,6.);
TH1D *eta_max        = new TH1D("eta_max","eta_max",60,0.,6.);
TH1D *eta_min        = new TH1D("eta_min","eta_min",60,0.,6.);
eta_total->Sumw2();
eta_max->Sumw2();
eta_min->Sumw2();

TH1D *etaprime_total = new TH1D("etaprime_total","etaprime_total",60,0.,6.);
TH1D *etaprime_max   = new TH1D("etaprime_max","etaprime_max",60,0.,6.);
TH1D *etaprime_min   = new TH1D("etaprime_min","etaprime_min",60,0.,6.);
etaprime_total->Sumw2();
etaprime_max->Sumw2();
etaprime_min->Sumw2();

TH1D *rho_total      = new TH1D("rho_total","rho_total",60,0.,6.);
TH1D *rho_max        = new TH1D("rho_max","rho_max",60,0.,6.);
TH1D *rho_min        = new TH1D("rho_min","rho_min",60,0.,6.);
rho_total->Sumw2();
rho_max->Sumw2();
rho_min->Sumw2();

TH1D *omega_total    = new TH1D("omega_total","omega_total",60,0.,6.);
TH1D *omega_max      = new TH1D("omega_max","omega_max",60,0.,6.);
TH1D *omega_min      = new TH1D("omega_min","omega_min",60,0.,6.);
omega_total->Sumw2();
omega_max->Sumw2();
omega_min->Sumw2();

TH1D *phi_total      = new TH1D("phi_total","phi_total",60,0.,6.);
TH1D *phi_max        = new TH1D("phi_max","phi_max",60,0.,6.);
TH1D *phi_min        = new TH1D("phi_min","phi_min",60,0.,6.);
phi_total->Sumw2();
phi_max->Sumw2();
phi_min->Sumw2();

TH1D *total_opt      = new TH1D("total_opt","total_opt",60,0.,6.);
TH1D *total_max      = new TH1D("total_max","total_max",60,0.,6.);
TH1D *total_min      = new TH1D("total_min","total_min",60,0.,6.);
total_opt->Sumw2();
total_max->Sumw2();
total_min->Sumw2();

TH1D *rel_error_max = new TH1D("rel_error_max","rel_error_max",60,0.,6.);
TH1D *rel_error_min = new TH1D("rel_error_min","rel_error_min",60,0.,6.);

double pt = 0.;
double error = 0.;
double average, max, min, rel_err, rel_err_max, rel_err_min;
double f = 0.;
double rel_err_f = 0.13/0.77;
double rel_err_meson = 0.50/3.;
double dsigmady_pi0 = 1.195*42.;

TFile *fbest = new TFile("cocktail_100M_final.root");
fbest->cd();

cout << "pion" << endl;
pion_total->Add(ptePion,pteCPion,dsigmady_pi0,dsigmady_pi0);
for (int i=1; i<=60; i++) {
  pt = (i-1)*0.1+0.05;
  f = 0.77*(1.0-0.0718*exp(-0.760*log(pt)));
  average = ptePion->GetBinContent(i);
  max = ptePionMax->GetBinContent(i);
  min = ptePionMin->GetBinContent(i);
  if ( average != 0 ) {
    rel_err = (max-average)/average;
    rel_err_max = (rel_err+f*sqrt(rel_err_f*rel_err_f+rel_err*rel_err))/(1+f);
    pion_max->SetBinContent(i,pion_total->GetBinContent(i)*(1+rel_err_max));
    rel_err = (average-min)/average;
    rel_err_min = (rel_err+f*sqrt(rel_err_f*rel_err_f+rel_err*rel_err))/(1+f);
    pion_min->SetBinContent(i,pion_total->GetBinContent(i)*(1-rel_err_min));
    cout << pt << " " << 100.0*rel_err_max << " " 
	 << 100.0*rel_err_min << endl;
  }
}

cout << endl;
cout << "eta" << endl;
eta_total->Add(pteEta,pteCEta,dsigmady_pi0,dsigmady_pi0);
for (int i=1; i<=60; i++) {
  pt = (i-1)*0.1+0.05;
  f = 0.77*(1.0-0.0718*exp(-0.760*log(pt)));
  average = pteEta->GetBinContent(i);
  max = pteEtaMax->GetBinContent(i);
  min = pteEtaMin->GetBinContent(i);
  if ( average != 0 ) {
    rel_err = (max-average)/average;
    rel_err_max = (rel_err+f*sqrt(rel_err_f*rel_err_f+rel_err*rel_err))/(1+f);
    rel_err_max = sqrt(rel_err_max*rel_err_max + rel_err_meson*rel_err_meson);
    eta_max->SetBinContent(i,eta_total->GetBinContent(i)*(1+rel_err_max));
    rel_err = (average-min)/average;
    rel_err_min = (rel_err+f*sqrt(rel_err_f*rel_err_f+rel_err*rel_err))/(1+f);
    rel_err_min = sqrt(rel_err_min*rel_err_min + rel_err_meson*rel_err_meson);
    eta_min->SetBinContent(i,eta_total->GetBinContent(i)*(1-rel_err_min));
    cout << pt << " " << 100.0*rel_err_max << " " 
	 << 100.0*rel_err_min << endl;
  }
}

cout << endl;
cout << "etaprime" << endl;
etaprime_total->Add(pteEtaprime,pteCEtaprime,dsigmady_pi0,dsigmady_pi0);
for (int i=1; i<=60; i++) {
  pt = (i-1)*0.1+0.05;
  f = 0.77*(1.0-0.0718*exp(-0.760*log(pt)));
  average = pteEtaprime->GetBinContent(i);
  max = pteEtaprimeMax->GetBinContent(i);
  min = pteEtaprimeMin->GetBinContent(i);
  if ( average != 0 ) {
    rel_err = (max-average)/average;
    rel_err_max = (rel_err+f*sqrt(rel_err_f*rel_err_f+rel_err*rel_err))/(1+f);
    rel_err_max = sqrt(rel_err_max*rel_err_max + rel_err_meson*rel_err_meson);
    etaprime_max->SetBinContent(i,etaprime_total->GetBinContent(i)*(1+rel_err_max));
    rel_err = (average-min)/average;
    rel_err_min = (rel_err+f*sqrt(rel_err_f*rel_err_f+rel_err*rel_err))/(1+f);
    rel_err_min = sqrt(rel_err_min*rel_err_min + rel_err_meson*rel_err_meson);
    etaprime_min->SetBinContent(i,etaprime_total->GetBinContent(i)*(1-rel_err_min));
    cout << pt << " " << 100.0*rel_err_max << " " 
	 << 100.0*rel_err_min << endl;
  }
}

cout << endl;
cout << "rho" << endl;
rho_total->Add(pteRho,pteRho,dsigmady_pi0,0.0);
for (int i=1; i<=60; i++) {
  pt = (i-1)*0.1+0.05;
  average = pteRho->GetBinContent(i);
  max = pteRhoMax->GetBinContent(i);
  min = pteRhoMin->GetBinContent(i);
  if ( average != 0 ) {
    rel_err = (max-average)/average;
    rel_err_max = sqrt(rel_err*rel_err + rel_err_meson*rel_err_meson);
    rho_max->SetBinContent(i,rho_total->GetBinContent(i)*(1+rel_err_max));
    rel_err = (average-min)/average;
    rel_err_min = sqrt(rel_err*rel_err + rel_err_meson*rel_err_meson);
    rho_min->SetBinContent(i,rho_total->GetBinContent(i)*(1-rel_err_min));
    cout << pt << " " << 100.0*rel_err_max << " " 
	 << 100.0*rel_err_min << endl;
  }
}

cout << endl;
cout << "omega" << endl;
omega_total->Add(pteOmega,pteCOmega,dsigmady_pi0,dsigmady_pi0);
for (int i=1; i<=60; i++) {
  pt = (i-1)*0.1+0.05;
  f = 0.77*(1.0-0.0718*exp(-0.760*log(pt)));
  average = pteOmega->GetBinContent(i);
  max = pteOmegaMax->GetBinContent(i);
  min = pteOmegaMin->GetBinContent(i);
  if ( average != 0 ) {
    rel_err = (max-average)/average;
    rel_err_max = (rel_err+f*sqrt(rel_err_f*rel_err_f+rel_err*rel_err))/(1+f);
    rel_err_max = sqrt(rel_err_max*rel_err_max + rel_err_meson*rel_err_meson);
    omega_max->SetBinContent(i,omega_total->GetBinContent(i)*(1+rel_err_max));
    rel_err = (average-min)/average;
    rel_err_min = (rel_err+f*sqrt(rel_err_f*rel_err_f+rel_err*rel_err))/(1+f);
    rel_err_min = sqrt(rel_err_min*rel_err_min + rel_err_meson*rel_err_meson);
    omega_min->SetBinContent(i,omega_total->GetBinContent(i)*(1-rel_err_min));
    cout << pt << " " << 100.0*rel_err_max << " " 
	 << 100.0*rel_err_min << endl;
  }
}

cout << endl;
cout << "phi" << endl;
phi_total->Add(ptePhi,ptePhi,dsigmady_pi0,0.0);
for (int i=1; i<=60; i++) {
  pt = (i-1)*0.1+0.05;
  average = ptePhi->GetBinContent(i);
  max = ptePhiMax->GetBinContent(i);
  min = ptePhiMin->GetBinContent(i);
  if ( average != 0 ) {
    rel_err = (max-average)/average;
    rel_err_max = sqrt(rel_err*rel_err + rel_err_meson*rel_err_meson);
    phi_max->SetBinContent(i,phi_total->GetBinContent(i)*(1+rel_err_max));
    rel_err = (average-min)/average;
    rel_err_min = sqrt(rel_err*rel_err + rel_err_meson*rel_err_meson);
    phi_min->SetBinContent(i,phi_total->GetBinContent(i)*(1-rel_err_min));
    cout << pt << " " << 100.0*rel_err_max << " " 
	 << 100.0*rel_err_min << endl;
  }
}

total_opt->Add(pion_total,eta_total,1.,1.);
total_opt->Add(total_opt,etaprime_total,1.,1.);
total_opt->Add(total_opt,rho_total,1.,1.);
total_opt->Add(total_opt,omega_total,1.,1.);
total_opt->Add(total_opt,phi_total,1.,1.);

cout << endl;
cout << "total systematic errors" << endl;
for (int i=1; i<=60; i++) {
  error = 0.0;
  error += (pion_max->GetBinContent(i)-pion_total->GetBinContent(i)); 
  error += (eta_max->GetBinContent(i)-eta_total->GetBinContent(i)); 
  error += (etaprime_max->GetBinContent(i)-etaprime_total->GetBinContent(i)); 
  error += (rho_max->GetBinContent(i)-rho_total->GetBinContent(i)); 
  error += (omega_max->GetBinContent(i)-omega_total->GetBinContent(i)); 
  error += (phi_max->GetBinContent(i)-phi_total->GetBinContent(i));
  total_max->SetBinContent(i,total_opt->GetBinContent(i)+error);
  error = 0.0;
  error += (pion_total->GetBinContent(i)-pion_min->GetBinContent(i)); 
  error += (eta_total->GetBinContent(i)-eta_min->GetBinContent(i)); 
  error += (etaprime_total->GetBinContent(i)-etaprime_min->GetBinContent(i)); 
  error += (rho_total->GetBinContent(i)-rho_min->GetBinContent(i)); 
  error += (omega_total->GetBinContent(i)-omega_min->GetBinContent(i)); 
  error += (phi_total->GetBinContent(i)-phi_min->GetBinContent(i)); 
  total_min->SetBinContent(i,total_opt->GetBinContent(i)-error);
  rel_err_max = 100.0*(total_max->GetBinContent(i)-total_opt->GetBinContent(i))
                /total_opt->GetBinContent(i);
  rel_err_min = 100.0*(total_opt->GetBinContent(i)-total_min->GetBinContent(i))
                /total_opt->GetBinContent(i);
  cout << (i-1)*0.1+0.05 << " " << rel_err_max << " " << rel_err_min << endl;
  rel_error_max->SetBinContent(i,rel_err_max);
  rel_error_min->SetBinContent(i,rel_err_min);
}

out->Write();
out->Close();

}




