#include "ATrack.h"

ClassImp(ATrack);

ATrack::ATrack() : AParticle(0){
  _phcentraltrack = NULL;
  _pc3sdphi = -9999.;  _pc3sdz = -9999.;
  _emcsdphi = -9999.;  _emcsdz = -9999.;
  _pc3dphi = -9999.;   _pc3dz = -9999.;
  _px = -9999.; _py = -9999.; _pz = -9999.;
  _n0 = -9999;
  _charge = -9999.;
  _qual = -9999;
  _alpha = -9999.;
  _phid = -9999.;
  _phi = -9999.;
  _theta = -9999.;
  _zed = -9999.;
  _dcarm = -1;
  _board = -9999.;
  _ecore = -9999.;
  _ppc3x = -9999.;  _ppc3y = -9999.; _ppc3z = -9999.;
  _ppc1x = -9999.;  _ppc1y = -9999.; _ppc1z = -9999.;
  _pemcx = -9999.;  _pemcy = -9999.;  _pemcz = -9999.;
}

void ATrack::SetBoard(float phi, int arm)
{
  if (arm == 0) _board = (3.72402 - phi + 0.008047 * cos(phi + 0.87851)) / 0.01963496;
  if (arm == 1) _board = (0.573231 + phi - 0.0046 * cos(phi + 0.05721)) / 0.01963496;
  
}

ATrack* ATrack::clone() const
{
  return new ATrack(*this);
}
