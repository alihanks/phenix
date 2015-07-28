//-----------------------------------------------------------------------------
//
// Generate mass distributions
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
double GounarisSakurai(float, double, double, double); 
double Lorentz(float, double, double); 
TH1F * InitializeM(int ParticleID, float mass_min, float mass_max, 
		   ParticlePropertyList *PPList)
{
  int    nbins, ibin;
  float  binwidth;
  float  mass_bin, weight;
  double vmass, vwidth;

  double pimass = 0.13956995;
  double emass  = 0.00051099906;

  ParticleProperty * PParticle = PPList->GetByID(ParticleID);
  vmass  = PParticle->GetMass();
  vwidth = PParticle->GetWidth();

  nbins = 10000;
  if ( mass_min == 0.0 && mass_max == 0.0 )
  {
    mass_min  = 2.0*pimass;
    mass_max  = 5.0;
  }
  binwidth  = (mass_max-mass_min)/(double)nbins;
  TH1F * mhistogram = new TH1F("mass","mass",nbins,mass_min,mass_max);

  for ( ibin=1; ibin<=nbins; ibin++ )
  {
    mass_bin = mass_min+(double)(ibin-1)*binwidth+binwidth/2.0;
    if ( ParticleID==113 || ParticleID==223 || ParticleID==333 )
    {
      weight   = (float)GounarisSakurai(mass_bin,vmass,vwidth,emass);
    }
    else	
    {
      weight = (float)Lorentz(mass_bin,vmass,vwidth);
    }	
    mhistogram->AddBinContent(ibin,weight);
  }

  mhistogram = InitializeRandomHist(mhistogram);

  return mhistogram;

}

double GounarisSakurai(float mass, double vmass, double vwidth, double emass)
{
  double corr = 0.0;
  double epsilon = 0.0;	
  double weight = 0.0;

  double pimass = 0.13956995;

  // corr = vwidth*(vmass/mass) * 
  //   exp(1.5*log((mass*mass/4.0-emass*emass)/(vmass*vmass/4.0-emass*emass)));
  corr = vwidth*(vmass/mass) * 
    exp(1.5*log((mass*mass/4.0-pimass*pimass)/
		(vmass*vmass/4.0-pimass*pimass)));
  epsilon = (emass/mass)*(emass/mass);
  
  if ( 1.0-4.0*epsilon>=0.0 ) 
  {
    weight = sqrt(1.0-4.0*epsilon)*(1.0+2.0*epsilon)/
             ((vmass*vmass-mass*mass)*(vmass*vmass-mass*mass)+
	      (vmass*corr)*(vmass*corr));
  }

  return weight;

}

double Lorentz(float mass, double vmass, double vwidth)
{	
  double weight;
  
  weight = (vwidth*vwidth/4.0)/(vwidth*vwidth/4.0+(vmass-mass)*(vmass-mass));

  return weight;

}








