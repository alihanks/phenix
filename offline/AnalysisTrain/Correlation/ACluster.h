#ifndef __ACLUSTER_H__
#define __ACLUSTER_H__

#include "AParticle.h"
#include <emcClusterContent.h>

class ACluster : public AParticle{

public:
  ACluster();
  virtual ~ACluster() {};

  void SetemcClusterContent(emcClusterContent* clus) { _emcClusterContent = clus; }
  void SetEcore(const float ecore) { _ecore = ecore; }
  void SetProb(const float prob) { _prob = prob; }
  void SetArm(const int arm) { _arm = arm; }
  void SetSec(const int sec) { _sec = sec; }
  void SetSm(const int sm) { _sm = sm; }
  void SetSm(const int iy, const int iz, const int sm);
  void SetIypos(const int iy) { _ypos = iy; }
  void SetIzpos(const int iz) { _zpos = iz; }
  void SetEmctrkdz(const float emctrkdz) { _emctrkdz = emctrkdz; }
  void SetEmctrkdphi(const float emctrkdphi) { _emctrkdphi = emctrkdphi; }
  void SetEmcpc3dz(const float emcpc3dz) { _emcpc3dz = emcpc3dz; }
  void SetEmcpc3dphi(const float emcpc3dphi) { _emcpc3dphi = emcpc3dphi; }
  void SetWarnDead(const unsigned int warnmap, const unsigned int deadmap) { _warnmap = warnmap; _deadmap = deadmap; }
  void SetTag(const bool tag) { _is_tagged = tag; }
  void SetTrigger(const bool trigger ) { _is_ert = trigger; }
  void SetFiducial(const bool fiducial ) { _is_fiducial = fiducial; }
  
  emcClusterContent* GetemcClusterContent() const { return _emcClusterContent; }
  float GetEcore() const { return _ecore; }
  float GetProb() const { return _prob; }
  int GetArm() const { return _arm; }
  int GetSec() const { return _sec;}
  int GetArmSect() const { return _arm * 4 + _sec;}
  int GetIypos() const { return _ypos; }
  int GetIzpos() const { return _zpos; }
  float GetEmctrkdz() const { return _emctrkdz; }
  float GetEmctrkdphi() const { return _emctrkdphi; }
  float GetEmcpc3dz() const { return _emcpc3dz; }
  float GetEmcpc3dphi() const { return _emcpc3dphi; }

  bool IsTrigger() const { return _is_ert; }
  bool IsTagged() const { return _is_tagged; }
  bool IsFiducial() const { return _is_fiducial; }
  unsigned int GetWarnStat() const { return _warnmap; }
  unsigned int GetDeadStat() const { return _deadmap; }

  ACluster* clone() const;

private:
  emcClusterContent *_emcClusterContent;
  float _ecore;
  float _prob;
  int _ypos;
  int _zpos;
  float _emctrkdz;
  float _emctrkdphi;
  float _emcpc3dz;
  float _emcpc3dphi;
  unsigned int _warnmap;
  unsigned int _deadmap;
  bool _is_tagged;
  bool _is_ert;
  bool _is_fiducial;

  ClassDef(ACluster,1);
};

#endif /* __ACLUSTER_H__ */
