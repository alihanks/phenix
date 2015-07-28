//-----------------------------------------------------------------------------
//
//  Apply CERES resolution to a four momentum
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>
#include <TH1.h>
#include <TRandom.h>
#include "Momentum.h"

#define PI  3.141592653589793238

double thetaof(Mom3);
double phiof(Mom3);

Mom4 ApplyCERESResolution(Mom4 mom4)
{
  double theta, phi, p, px, py, pz, E, mass, sigma;
  Mom4   result;

  mass  = sqrt(mom4*mom4);
  E     = mom4.GetE();
  p     = sqrt(E*E-mass*mass);
  theta = thetaof(mom4.Getp());
  phi   = phiof(mom4.Getp());

  sigma = p*sqrt(0.053*p*0.053*p+0.041*0.041);
  p     = gRandom->Gaus(p,sigma);
  sigma = 0.6e-03;
  theta = gRandom->Gaus(theta,sigma);
  sigma = 3.0e-03;
  phi   = gRandom->Gaus(phi,sigma);

  E     = sqrt(p*p+mass*mass);
  px    = p*sin(theta)*cos(phi);
  py    = p*sin(theta)*sin(phi);
  pz    = p*cos(theta);

  result.SetE(E);
  result.Setp(px,py,pz);

  return (result);

}















