#include "ACluster.h"
#include <iostream>

ClassImp(ACluster);

ACluster::ACluster() : AParticle(22)
{
  _emcClusterContent = NULL;
  _is_ert = true;
  _is_iso = false;
  _is_tagged = false;
  _is_fiducial = true;
  _prob = -9999;
  _ypos = -9999;
  _zpos = -9999;
  _x = -9999;
  _y = -9999;
  _z = -9999;
  _warnmap = 0;
  _deadmap = 0;
  _ecore = -9999;
  _arm = -1;
  _sec = -1;
  _emctrkdz = -9999;
  _emctrkdphi = -9999;
  _emcpc3dz = -9999;
  _emcpc3dphi = -9999;
}

void ACluster::SetSm(int iy, int iz, int sm)
{
  if( sm < 32 && sm > 0 ) _sm = sm;
  else _sm = int((floor(iy/12)+1)*floor(iz/12));
}

ACluster* ACluster::clone() const
{
  return new ACluster(*this);
}


