#ifndef __APARTICLE_H__
#define __APARTICLE_H__

#include <TLorentzVector.h>

class AParticle : public TLorentzVector {

public:
  AParticle(int pid);
  AParticle(AParticle* daughter1, AParticle* daughter2, int pid);
  virtual ~AParticle() {}
  
  int GetPid() const { return _pid; }
  AParticle* Daughter1() const { return _daughter1; }
  AParticle* Daughter2() const { return _daughter2; }

  void SetEmcX(const float x) { _x = x; }
  void SetEmcY(const float y) { _y = y; }
  void SetEmcZ(const float z) { _z = z; }

  void SetEvent(int event) { _event = event; }
  void SetZvtx(float vtx) { _zvtx = vtx; }
  void SetCent(float cent) { _cent = cent; }

  float GetX() const { return _x; }
  float GetY() const { return _y; }
  float GetZ() const { return _z; }

  int   GetEvent() const { return _event; }
  float GetZvtx() const { return _zvtx; }
  float GetCent() const { return _cent; }

  unsigned int GetArm() const { return _arm; }
  unsigned int GetSec() const { return _sec; }
  unsigned int GetSm() const { return _sm; }
  
  void SetIso(const bool iso) { _is_iso = iso; }
  bool IsIso() const { return _is_iso; }
  virtual bool IsTagged() const =0;
  virtual bool IsFiducial() const=0;

  virtual AParticle* clone() const =0;

protected:
  AParticle* _daughter1;
  AParticle* _daughter2;
  float _x;
  float _y;
  float _z;

  int   _event;
  float _zvtx;
  float _cent;

  int   _pid;

  unsigned int _arm;
  unsigned int _sec;
  unsigned int _sm;

  bool _is_iso;

  ClassDef(AParticle,1);
};

#endif /* __APARTICLE_H__ */
