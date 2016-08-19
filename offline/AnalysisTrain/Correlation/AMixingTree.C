#include "AMixingTree.h"
#include <TTree.h>
#include <iostream>

AMixingTree::AMixingTree()
{
  _ttrig = new TTree("trigger_tree", "trigger tree");
  _ttrig->SetDirectory(0);
  _tpart = new TTree("partner_tree", "partner tree");
  _tpart->SetDirectory(0);
  _evt = 0;
  _zvtx = 0.;
  _cent = 0.;
  _mix_evt = 0;
  _mix_zvtx = 0.;
  _mix_cent = 0.;
  _ntrig_photons = 0;
  _ntrig_pi0s = 0;
  _ntrig = 0;
  _npart = 0;
  _nallphotons = 0;

  memset(_pt_trig, 0., DEPTH);
  memset(_phi_trig, 0., DEPTH);
  memset(_eta_trig, 0., DEPTH);
  memset(_e_trig, 0., DEPTH);
  memset(_x, 0., DEPTH);
  memset(_y, 0., DEPTH);
  memset(_z, 0., DEPTH);
  memset(_iso, 0, DEPTH);

  memset(_pt_part, 0., DEPTH);
  memset(_phi_part, 0., DEPTH);
  memset(_eta_part, 0., DEPTH);
  memset(_e_part, 0., DEPTH);
  memset(_pemcx, 0., DEPTH);
  memset(_pemcy, 0., DEPTH);
  memset(_pemcz, 0., DEPTH);
  memset(_pt_clus, 0., DEPTH);
  memset(_phi_clus, 0., DEPTH);
  memset(_eta_clus, 0., DEPTH);
  memset(_e_clus, 0., DEPTH);
}

void AMixingTree::SetEventData(int evt, float zvtx, float cent, int nphotons, int npi0s, int ntracks)
{
  _evt = evt;
  _zvtx = zvtx;
  _cent = cent;
  _mix_evt = evt;
  _mix_zvtx = zvtx;
  _mix_cent = cent;
  _ntrig_photons = nphotons;
  _ntrig_pi0s = npi0s;
  _ntrig = nphotons+npi0s;
  _npart = ntracks;
}

void AMixingTree::SetMBEventData(int evt, float zvtx, float cent, int nphotons, int ntracks)
{
  _mix_evt = evt;
  _mix_zvtx = zvtx;
  _mix_cent = cent;
  _nallphotons = nphotons;
  _npart = ntracks;
}

void AMixingTree::SetTriggerData(float pt, float phi, float eta, float e, float x, float y, float z, int iso, int index)
{
  _pt_trig[index] = pt;
  _phi_trig[index] = phi;
  _eta_trig[index] = eta;
  _e_trig[index] = e;
  _x[index] = x;
  _y[index] = y;
  _z[index] = z;
  _iso[index] = iso;
}

void AMixingTree::SetPartnerData(float pt, float phi, float eta, float e, float pemcx, float pemcy, float pemcz, int index)
{
  _pt_part[index] = pt;
  std::cout << "AMixingTree::SetPartnerData: pt = " << std::endl;
  _phi_part[index] = phi;
  _eta_part[index] = eta;
  _e_part[index] = e;
  _pemcx[index] = pemcx;
  _pemcy[index] = pemcy;
  _pemcz[index] = pemcz;
}

void AMixingTree::SetPhotonData(float pt, float phi, float eta, float e, int index)
{
  _pt_clus[index] = pt;
  _phi_clus[index] = phi;
  _eta_clus[index] = eta;
  _e_clus[index] = e;
}

void AMixingTree::SetTriggerBranches()
{
  _ttrig->Branch("evt", &_evt, "evt/I");
  _ttrig->Branch("zvtx", &_zvtx, "zvtx/F");
  _ttrig->Branch("cent", &_cent, "cent/F");
  _ttrig->Branch("nphotons", &_ntrig_photons, "nphotons/I");
  _ttrig->Branch("npi0s", &_ntrig_pi0s, "npi0s/I");
  _ttrig->Branch("ntrig", &_ntrig, "ntrig/I");
  _ttrig->Branch("pt", _pt_trig, "pt[ntrig]/F");
  _ttrig->Branch("phi", _phi_trig, "phi[ntrig]/F");
  _ttrig->Branch("eta", _eta_trig, "eta[ntrig]/F");
  _ttrig->Branch("e", _e_trig, "e[ntrig]/F");
  _ttrig->Branch("x", _x, "x[ntrig]/F");
  _ttrig->Branch("y", _y, "y[ntrig]/F");
  _ttrig->Branch("z", _z, "z[ntrig]/F");
  _ttrig->Branch("iso", _iso, "iso[ntrig]/I");
}

void AMixingTree::SetPartnerBranches()
{
  _tpart->Branch("evt", &_mix_evt, "evt/I");
  _tpart->Branch("zvtx", &_mix_zvtx, "zvtx/F");
  _tpart->Branch("cent", &_mix_cent, "cent/F");
  _tpart->Branch("ntracks", &_npart, "ntracks/I");
  _tpart->Branch("nphotons", &_nallphotons, "nphotons/I");
  _tpart->Branch("pt", _pt_part, "pt[ntracks]/F");
  _tpart->Branch("phi", _phi_part, "phi[ntracks]/F");
  _tpart->Branch("eta", _eta_part, "eta[ntracks]/F");
  _tpart->Branch("e", _e_part, "e[ntracks]/F");
  _tpart->Branch("pemcx", _pemcx, "pemcx[ntracks]/F");
  _tpart->Branch("pemcy", _pemcy, "pemcy[ntracks]/F");
  _tpart->Branch("pemcz", _pemcz, "pemcz[ntracks]/F");
  _tpart->Branch("pt_clus", _pt_clus, "pt_clus[nphotons]/F");
  _tpart->Branch("phi_clus", _phi_clus, "phi_clus[nphotons]/F");
  _tpart->Branch("eta_clus", _eta_clus, "eta_clus[nphotons]/F");
  _tpart->Branch("e_clus", _e_clus, "e_clus[nphotons]/F");
}


