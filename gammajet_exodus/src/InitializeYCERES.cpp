//-----------------------------------------------------------------------------
//
// Generate rapidity distributions for the CERES generator
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>

#include <TH1.h>
#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticlePropertyList.h"

#define PI 3.141592653589793238

TH1F * InitializeRandomHist(TH1F*);
double y_shift(double, double, double, double);

TH1F * InitializeYCERES(int ParticleID, ParticlePropertyList *PPList)
{
  int    nbins, ibin;
  float  binwidth;
  float  y, ymin, ymax, weight;
  double mass;
  double sqrts, y0;
  double gamma, sigma_landau, sig, sigma;
  double emax_pi0, pmax_pi0, ymax_pi0;
  double emax_part, pmax_part, ymax_part, ymin_part;
  double ymean;
  
  double mnuc  = 0.938;
  double ebeam = 450.0;
  double Aproj = 1.0;
  double Atarg = 9.0;
  double bpar  = 0.0;
  double m_pi0 = 0.1349743;

  sqrts        = sqrt(2.0*mnuc*ebeam+2.0*mnuc*mnuc);
  y0           = log(2.0*ebeam/mnuc)/2.0;
  gamma        = sqrts/(2.0*mnuc);
  sigma_landau = sqrt(log(gamma));
  emax_pi0     = (sqrts*sqrts-4.0*mnuc*mnuc+m_pi0*m_pi0)/(2.0*sqrts);
  pmax_pi0     = sqrt(emax_pi0*emax_pi0-m_pi0*m_pi0);
  ymax_pi0     = log(sqrts/m_pi0);

  ParticleProperty * PParticle = PPList->GetByID(ParticleID);
  mass      = PParticle->GetMass();
  emax_part = (sqrts*sqrts-4.0*mnuc*mnuc+mass*mass)/(2.0*sqrts);
  pmax_part = sqrt(emax_pi0*emax_pi0-mass*mass);
  ymax_part = log(sqrts/mass);
  ymin_part = -1.0*ymax_part;
  sig       = sigma_landau*(ymax_part/ymax_pi0);

  ymax_part = ymax_part*y0;
  ymin_part = ymin_part*y0;

  ymean = y0-y_shift(ebeam,Aproj,Atarg,bpar);

  nbins = 1000;
  ymin  = -2.0;
  ymax  = 8.0;
  binwidth = (ymax-ymin)/(double)nbins;
  TH1F * yhistogram = new TH1F("y","y",nbins,ymin,ymax);

  for ( ibin=1; ibin<=nbins; ibin++ )
  {
    y  = ymin+(double)(ibin-1)*binwidth+binwidth/2.0;
    if ( y<=ymean )
    {
      sigma = sig*(1.0-((y0-ymean)*(y0-ymean)/(ymax_part-ymin_part)));
    }
    else
    {
      sigma = sig*(1.0+((y0-ymean)*(y0-ymean)/(ymax_part-ymin_part)));
    }
    weight = exp(0.75*(y-ymean)/sigma)+exp(-0.75*(y-ymean)/sigma);
    weight = 1.0/(weight*weight);
    yhistogram->AddBinContent(ibin,weight);
  }
  
  yhistogram = InitializeRandomHist(yhistogram);

  return yhistogram;

}

double y_shift(double ebeam, double Aproj, double Atarg, double bpar)
{
  double density = 0.168;
  double rmin, rmax;
  double hmax, rcyl, vmax, hmic, rapb, vmin, rapc;
  double t_part, p_part, shift;

  if ( Aproj<2.0 )
  {
    if ( Atarg<12.0 ) return 0.0;
    return 1.0;
  }
  else
  {
    if ( sqrt((Atarg/Aproj-1.0)*(Atarg/Aproj-1.0))<0.1 ) return 0.0;
  }

  rmin = exp((1.0/3.0)*log(0.75*Aproj/PI/density));
  rmax = exp((1.0/3.0)*log(0.75*Atarg/PI/density));
  if ( bpar>(rmin+rmax) ) return 0.0;

  hmax = 2.0*sqrt(rmax*rmax-rmin*rmin);

  if ( bpar<(rmax-rmin) )
  {
    rcyl   = rmin;
    vmax   = 2.0*PI*rcyl*rcyl*hmax;
    hmic   = rmax-sqrt(rmax*rmax-rmin*rmin);
    rapb   = (rmax-hmic)/rmax;
    vmin   = 2.0/3.0*PI*rmax*rmax*rmax*(1.0+0.5*rapb*rapb*rapb-1.5*rapb);
    t_part = density*(2.0*vmin+vmax);
    p_part = Aproj;
  }
  else
  {
    rcyl   = 0.5*(rmax+rmin-bpar);
    vmax   = 2.0*PI*rcyl*rcyl*hmax;
    hmic   = rmax-sqrt(rmax*rmax-rmin*rmin);
    rapb   = (rmax-hmic)/rmax;
    vmin   = 2.0/3.0*PI*rmax*rmax*rmax*(1.0+0.5*rapb*rapb*rapb-1.5*rapb);
    rapc   = (rmin-2.0*rcyl)/rmin;
    t_part = density*(2.0/3.0*PI*rcyl*rcyl*hmic+vmax);
    p_part = density*2.0/3.0*PI*rmin*rmin*rmin*
             (1.0+0.5*rapc*rapc*rapc-1.5*rapc);
  }

  shift = 0.0055*(t_part-p_part);
  return shift;

}









