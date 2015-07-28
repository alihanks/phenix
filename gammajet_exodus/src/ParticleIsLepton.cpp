//-----------------------------------------------------------------------------
//
//  Check whether a perticle is a lepton or not
//
//-----------------------------------------------------------------------------

#include "Momentum.h"
#include "Particle.h"
#include <math.h>
bool ParticleIsLepton(Particle *PParticle)
{
  int pid = PParticle->GetID();
  
  if ( fabs(pid)==11 || fabs(pid)==13 )
  {
    return(true);
  }
  else
  {
    return(false);
  }

}
