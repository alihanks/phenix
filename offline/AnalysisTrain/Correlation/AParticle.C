#include "AParticle.h"

ClassImp(AParticle);

AParticle::AParticle(int pid) : TLorentzVector(0,0,0,0) {
  _is_iso = false;
  _daughter1 = NULL;
  _daughter2 = NULL;
  _pid = pid;
  _sec = 9999;
  _sm = 9999;
  _arm = 9999;
  _x = -9999.;
  _y = -9999.;
  _z = -9999.;
  _event = 0;
  _zvtx = 0.;
  _cent = 0.;
}

AParticle::AParticle(AParticle* daughter1, AParticle* daughter2, int pid) : TLorentzVector(*daughter1+*daughter2) {

  _pid = pid;

  if( daughter1 ) {
    _daughter1 = daughter1->clone();
  }
  if( daughter2 ){
    _daughter2 = daughter2->clone();
  }

  if( _daughter1 ) {
    _sec = _daughter1->GetSec();
    _sm = _daughter1->GetSm();
    _arm = _daughter1->GetArm();
  }
}
