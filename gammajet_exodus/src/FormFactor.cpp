//-----------------------------------------------------------------------------
//
// Calculate the form factor for Dalitz decays according to the
//
//               Lepton-G parametrization
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>

#include "Momentum.h"
#include "Particle.h"

#define SQR(a) ((a)*(a))

double FormFactor(double q2, ParticleProperty* PParent)
{
  double formfactor;
  int ParentID;

  ParentID = PParent->GetID();

  switch( ParentID )
  {
    case 111:
    {
      formfactor = 1.0/(1.0-5.5*q2);
      formfactor = SQR(formfactor);
    } break;
    case 221:
    {
      formfactor = 1.0/(1.0-1.9*q2);
      formfactor = SQR(formfactor);
    } break;
    case 331:
    {
      formfactor = (SQR(SQR(0.764))) / 
                   (SQR(SQR(0.764)-q2)+SQR(0.1020)*SQR(0.764));
    } break;
    case 223:
    {
      formfactor = (SQR(SQR(0.6519))) / 
                   (SQR(SQR(0.6519)-q2)+SQR(0.04198)*SQR(0.6519));
    } break;

    default: formfactor = 1.0;
  }

  return formfactor;

}
