//-----------------------------------------------------------------------------
//
// Generate transverse-momentum distributions according to hydrodynamics
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>
#include <TH1.h>
#include <TMath.h>

float InitializePtHydro(double mass, double t_fo, double beta, float pt) 
{
  float sum    = 0.0;
  float weight = 1.0;
  double profx, valprofx, mt;
  double besi, besk, besiarg, beskarg, arg;

  for (int profbin=1; profbin<101; profbin++)
  {
    profx    = ((double)(profbin-1) + 0.5) / 100.0;
    mt       = sqrt(mass*mass+pt*pt);
    arg      = atanh(beta*profx);
    besiarg  = pt*sinh(arg)/t_fo;;
    beskarg  = mt*cosh(arg)/t_fo;;
    besi     = TMath::BesselI0(besiarg); 
    besk     = TMath::BesselK1(beskarg); 
    valprofx = mt*profx*besi*besk;

    sum = sum + valprofx;
  }

  weight = pt*sum*0.01; 

  return weight;

}





