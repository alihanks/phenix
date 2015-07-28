//-----------------------------------------------------------------------------
//
// Generate rapidity distributions for the ISR generator
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

TH1F * InitializeYISR(int ParticleID, ParticlePropertyList *PPList)
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

  ymean = y0;

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










