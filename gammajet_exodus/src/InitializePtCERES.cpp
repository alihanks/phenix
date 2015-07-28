//-----------------------------------------------------------------------------
//
// Generate transverse-momentum distributions for the CERES generator
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

TH1F * InitializeRandomHist(TH1F*);

TH1F * InitializePtCERES(int ParticleID, ParticlePropertyList *PPList)
{
  int    nbins, ibin;
  float  binwidth;
  float  pt, ptmin, ptmax, mt, weight;
  double mass;
  double a1, a2, a3, T1, T2, T3;
  double T;

  ParticleProperty * PParticle = PPList->GetByID(ParticleID);
  mass = PParticle->GetMass();

  a1 = 1.0;
  a2 = 0.139;
  a3 = 0.107;
  T1 = 0.1;
  T2 = 0.23;
  T3 = 0.102;

  T  = 0.175+0.115*mass;

  nbins = 1000;
  ptmin = 0.0;
  ptmax = 4.0;
  binwidth = (ptmax-ptmin)/(double)nbins;
  TH1F * pthistogram = new TH1F("pt","pt",nbins,ptmin,ptmax);

  for ( ibin=1; ibin<=nbins; ibin++ )
  {
    pt = ptmin+(double)(ibin-1)*binwidth+binwidth/2.0;
    mt = sqrt(pt*pt+mass*mass);
    if ( ParticleID==111 )
    {
      weight = pt*(a1*exp(-1.0*mt/T1)+a2*exp(-1.0*mt/T2)+a3*exp(-1.0*mt/T3));
    }
    else
    {
      weight = pt*exp(-1.0*mt/T);
    }
    pthistogram->AddBinContent(ibin,weight);
  }
  
  pthistogram = InitializeRandomHist(pthistogram);

  return pthistogram;

}





