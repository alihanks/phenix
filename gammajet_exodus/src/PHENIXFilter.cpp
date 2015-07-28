//-----------------------------------------------------------------------------
//
//  Check whether a particle hits PHENIX and return the sector ID
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticlePropertyList.h"
#include "InPHENIXAcceptance.h"

#define PI  3.141592653589793238

double thetaof(Mom3);
double phiof(Mom3);

int PHENIXFilter(Particle *PParticle, ParticlePropertyList *PPList)
{
  Mom4 mom;
  int pid, charge, sector;
  double p, theta, phi, mass;

  mom    = PParticle->Get4mom();
  pid    = PParticle->GetID();
  p      = sqrt(mom.GetE()*mom.GetE()-mom*mom);
  theta  = 180.0*thetaof(mom.Getp())/PI;
  phi    = 180.0*phiof(mom.Getp())/PI;
  mass   = PPList->GetByID(pid)->GetMass();
  charge = PPList->GetByID(pid)->GetCharge();
  sector = InPHENIXAcceptance(p,theta,phi,charge,mass);

  return sector;
}

