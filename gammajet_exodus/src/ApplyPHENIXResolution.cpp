//-----------------------------------------------------------------------------
//
//  Apply PHENIX resolution to a four momentum
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>
#include <TH1.h>
#include <TRandom.h>
#include "Momentum.h"

#define PI  3.141592653589793238

double thetaof(Mom3);
double phiof(Mom3);

Mom4 ApplyPHENIXResolution(Mom4 mom4)
{
  double theta, phi, p, px, py, pz, E, mass, sigma;
  Mom4   result;

  px = mom4.Getp().Getpx();
  py = mom4.Getp().Getpy();
  pz = mom4.Getp().Getpz();
  p  = sqrt(px*px+py*py+pz*pz);

  mass = mom4*mom4;
  if ( mass<0.0001 )
  {
    mass = 0.51099906e-3;
  }
  else
  {
    mass = sqrt(mass);
  }

  theta = thetaof(mom4.Getp());
  phi   = phiof(mom4.Getp());

  sigma = p*sqrt(p*p*(0.01)*(0.01)+0.007*0.007);
  //  sigma = p*sqrt(p*p*(3./84.)*(3./84.)+0.006*0.006);
  p     = gRandom->Gaus(p,sigma);

  E     = sqrt(p*p+mass*mass);
  px    = p*sin(theta)*cos(phi);
  py    = p*sin(theta)*sin(phi);
  pz    = p*cos(theta);

  result.SetE(E);
  result.Setp(px,py,pz);

  return (result);

}















