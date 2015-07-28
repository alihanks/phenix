//-----------------------------------------------------------------------------
//
//  Check whether a particle hits PHENIX and return the sector ID
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>

#include "Momentum.h"
#include "Particle.h"

#define PI acos(-1.)

double thetaof(Mom3);
double phiof(Mom3);

bool PHENIXPhotonFilter(Particle *PParticle, int fidcut = 0)
{
  Mom4 mom;
  bool accept;
  double theta, eta, phi;

  accept = true;

  mom    = PParticle->Get4mom();
  theta  = thetaof(mom.Getp());
  eta    = -1.0*log(tan(theta/2.));
  if ( fabs(eta)>0.35 ) accept = false;

  // very simple minded!
  phi    = 180.0*phiof(mom.Getp())/PI - 90.0;
  //  if ( !((phi>-34. && phi<56.) || (phi>124. && phi<214.)) )  //whole range
  if ( !((phi>-34. && phi<11.) || (phi>124. && phi<214.)) )  //no PbGl
    accept = false;
  
  if (fidcut)
    {
      if ( fabs(eta)>0.20 ) accept = false;
      if ( !((phi>-25. && phi<46.) || (phi>134. && phi<204.)) )  //no PbGl
	accept = false;
    }

  return accept;
}

