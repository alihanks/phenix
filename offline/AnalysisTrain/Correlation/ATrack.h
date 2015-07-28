#ifndef __ATRACK_H__
#define __ATRACK_H__

#include "AParticle.h"
#include <PHCentralTrack.h>

class ATrack : public AParticle {

public:
  ATrack();
  virtual ~ATrack() {};

  void SetPHCentralTrack(PHCentralTrack* phcentraltrack) { _phcentraltrack = phcentraltrack; }
  void SetPxPyPz(float px, float py, float pz) { _px = px; _py = py; _pz = pz; }
  void SetPc3Match(float sdphi, float sdz) { _pc3sdphi = sdphi; _pc3sdz = sdz; }
  void SetEmcMatch(float sdphi, float sdz) { _emcsdphi = sdphi; _emcsdz = sdz; }
  void SetPc3dz(float dz) { _pc3dz = dz;}
  void SetPc3dphi(float dphi) { _pc3dphi = dphi; }
  void SetCharge(float charge) { _charge = charge; }
  void SetN0(int n0) { _n0 = n0; }
  void SetQuality(float qual) { _qual = qual; }
  void SetPhiD(float phid) { _phid = phid; }
  void SetPhi(float phi) { _phi = phi; }
  void SetTheta(float theta) { _theta = theta; }
  void SetAlpha(float alpha) { _alpha = alpha; }
  void SetZed(float zed) { _zed = zed; }
  void SetPpc3(float ppc3x, float ppc3y, float ppc3z) { _ppc3x = ppc3x; _ppc3y = ppc3y; _ppc3z = ppc3z; }
  void SetPpc1(float ppc1x, float ppc1y, float ppc1z) { _ppc1x = ppc1x; _ppc1y = ppc1y; _ppc1z = ppc1z; }
  void SetPemc(float pemcx, float pemcy, float pemcz) { _pemcx = pemcx; _pemcy = pemcy; _pemcz = pemcz; }
  void SetDCArm(float dcarm) { _dcarm = dcarm; }
  void SetBoard(float phi, int arm);
  void SetEcore(float ecore) { _ecore = ecore; }
  bool IsIso() const { return 0; }
  bool IsTagged() const { return 0; }
  bool IsFiducial() const { return true; }

  PHCentralTrack* GetPHCentralTrack() const { return _phcentraltrack; }
  float GetPpc3x() const { return _ppc3x; }
  float GetPpc3y() const { return _ppc3y; }
  float GetPpc3z() const { return _ppc3z; }
  float GetPpc1x() const { return _ppc1x; }
  float GetPpc1y() const { return _ppc1y; }
  float GetPpc1z() const { return _ppc1z; }
  float GetPemcx() const { return _pemcx; }
  float GetPemcy() const { return _pemcy; }
  float GetPemcz() const { return _pemcz; }

  float GetPc3dphi() const { return _pc3dphi; }
  float GetPc3dz() const { return _pc3dz; }
  float GetPc3sdphi() const { return _pc3sdphi; }
  float GetPc3sdz() const { return _pc3sdz; }
  float GetEmcsdphi() const { return _emcsdphi; }
  float GetEmcsdz() const { return _emcsdz; }

  
  int   GetN0() const { return _n0; }
  float GetCharge() const { return _charge; }
  float GetQuality() const { return _qual; }
  float GetPhiD() const { return _phid; }
  float GetPhi() const { return _phi; }
  float GetTheta() const { return _theta; }
  float GetAlpha() const { return _alpha; }
  float GetZed() const { return _zed; }
  int   GetDCArm() const { return _dcarm; }
  float GetBoard() const { return _board; }
  float GetEcore() const { return _ecore; }

  ATrack* clone() const;

private:
  PHCentralTrack *_phcentraltrack;
  float _px;
  float _py;
  float _pz;
  float _pc3dphi;
  float _pc3dz;
  float _pc3sdphi;
  float _pc3sdz;
  float _emcsdphi;
  float _emcsdz;
  int   _n0;
  float _charge;
  float _qual;
  float _alpha;
  float _phid;
  float _phi;
  float _theta;
  float _zed;
  int   _dcarm;
  float _board;
  float _ppc3x, _ppc3y, _ppc3z;
  float _ppc1x, _ppc1y, _ppc1z;
  float _pemcx, _pemcy, _pemcz;
  float _ecore;


  ClassDef(ATrack,1);
  
};

#endif /* __ATrack_H__ */
