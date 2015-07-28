#ifndef __AMIXINGTREE_H__
#define __AMIXINGTREE_H__

#include <TTree.h>
#include <iostream>

#define DEPTH 500

class AMixingTree {
public:
  AMixingTree();
  virtual ~AMixingTree(){ delete _ttrig; delete _tpart; }

  TTree* _ttrig;
  TTree* _tpart;

  void SetEventData(int evt, float zvtx, float cent, int nphotons, int npi0s, int ntracks);
  void SetMBEventData(int evt, float zvtx, float cent, int nphotons, int ntracks);
  void SetTriggerData(float pt, float phi, float eta, float e, float x, float y, float z, int iso, int index);
  void SetPartnerData(float pt, float phi, float eta, float e, float pemcx, float pemcy, float pemcz, int index);
  void SetPhotonData(float pt, float phi, float eta, float e, int index);

  void SetTriggerBranches();
  void SetPartnerBranches();

private:
  
  //event info
  int _evt;
  float _zvtx;
  float _cent;
  int _mix_evt;
  float _mix_zvtx;
  float _mix_cent;

  //trig info
  float _pt_trig[DEPTH];
  float _phi_trig[DEPTH];
  float _eta_trig[DEPTH];
  float _e_trig[DEPTH];
  float _x[DEPTH];
  float _y[DEPTH];
  float _z[DEPTH];
  int _iso[DEPTH];
  int _ntrig_photons;
  int _ntrig_pi0s;
  int _ntrig;

  //partner info
  float _pt_part[DEPTH];
  float _phi_part[DEPTH];
  float _eta_part[DEPTH];
  float _e_part[DEPTH];
  float _pemcx[DEPTH];
  float _pemcy[DEPTH];
  float _pemcz[DEPTH];
  float _pt_clus[DEPTH];
  float _phi_clus[DEPTH];
  float _eta_clus[DEPTH];
  float _e_clus[DEPTH];
  int _npart;
  int _nallphotons;
};

#endif /* __AMIXINGTREE_H__ */
