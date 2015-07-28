//-----------------------------------------------------------------------------
//
//  Implementation of some useful tools acting on momenta and particles
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <math.h>

#include "Momentum.h" 
#include "Particle.h"
#include "ParticlePropertyList.h"

#define SQR(a) ((a)*(a))
#define PI  3.141592653589793238

Mom3 invert(Mom3);
double abs(Mom4);
Mom4 boost(Mom4, Mom4);
double thetaof(Mom3);
double phiof(Mom3);
Mom3 Rotate(Mom3, double, double);
double GetMass(int, ParticlePropertyList *);
double GetWidth(int, ParticlePropertyList *);

Mom3 invert(Mom3 input)
{
  Mom3 output;
  output.Setpx(-input.Getpx());
  output.Setpy(-input.Getpy());
  output.Setpz(-input.Getpz());
  return output;
}

double abs(Mom4 input)
{
  return sqrt(input*input);
}

Mom4 boost(Mom4 pinit, Mom4 pframe)
{
  Mom4 pfinal;
  Mom3 eta, pl;
  double gamma;
  double massofparent = abs(pframe);
  eta = pframe.Getp()*(-1./massofparent);
  gamma = pframe.GetE()/massofparent;
  if ( eta*eta==0 )
    pl.Set(0.,0.,0.);
  else
    pl = eta*((pinit.Getp()*eta)/(eta*eta));
  pfinal.SetE(gamma*pinit.GetE()-eta*pinit.Getp());
  pfinal.Setp(pinit.Getp()+pl*(gamma-1)-eta*pinit.GetE());
  return pfinal;
}

double thetaof(Mom3 p)
{
  double absp = sqrt(p*p);
  if ( absp>0.0 ) 
    return acos(p.Getpz()/absp);
  else 
    return 0.0;
}

double phiof (Mom3 p)
{
  double phi;
  if ( p.Getpy()>=0 )
    phi = atan2(p.Getpy(),p.Getpx());
  else 
    phi = atan2(p.Getpy(),p.Getpx())+2.0*PI;
  return phi;
}

Mom3 Rotate(Mom3 pold, double costheta, double sintheta,
                       double cosphi,   double sinphi)
{
  Mom3 pnew;
  pnew.Setpx(pold.Getpx()*costheta*cosphi
            -pold.Getpy()*sinphi
            +pold.Getpz()*sintheta*cosphi);
  pnew.Setpy(pold.Getpx()*costheta*sinphi
            +pold.Getpy()*cosphi
            +pold.Getpz()*sintheta*sinphi);
  pnew.Setpz(-pold.Getpx()*sintheta
             +pold.Getpz()*costheta);
  return pnew;
}

double GetMass(int id, ParticlePropertyList *Index)
{
  return Index->GetByID(id)->GetMass();
}

double GetWidth(int id, ParticlePropertyList *Index)
{
  return Index->GetByID(id)->GetWidth();
}





