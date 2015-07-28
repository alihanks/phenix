//-----------------------------------------------------------------------------
//
// Generate rapidity distributions for the PHENIX generator
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

TH1F * InitializeYPHENIX(int ParticleID, ParticlePropertyList *PPList)
{
  int    nbins, ibin;
  float  binwidth;
  float  y, ymin, ymax, weight;

  nbins = 1000;
  ymin  = -0.5;
  ymax  = 0.5;
  binwidth = (ymax-ymin)/(double)nbins;
  TH1F * yhistogram = new TH1F("y","y",nbins,ymin,ymax);

  for ( ibin=1; ibin<=nbins; ibin++ )
  {
    y  = ymin+(double)(ibin-1)*binwidth+binwidth/2.0;
    weight = 1.0;
    yhistogram->AddBinContent(ibin,weight);
  }
  
  //  yhistogram = InitializeRandomHist(yhistogram);

  return yhistogram;

}










