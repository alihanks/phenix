//-----------------------------------------------------------------------------
//
// Generate white rapidity distributions 
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>

#include <TH1.h>
#include <TRandom.h>

TH1F * InitializeRandomHist(TH1F*);

TH1F * InitializeWhiteY(double YMin, double YMax)
{
  int    nbins, ibin;
  float  binwidth;
  float  weight;
  
  nbins = 1000;
  binwidth = (YMax-YMin)/(double)nbins;
  TH1F * yhistogram = new TH1F("y","y",nbins,YMin,YMax);

  for ( ibin=1; ibin<=nbins; ibin++ )
  {
    weight = 1.0;
    yhistogram->AddBinContent(ibin,weight);
  }
  
  yhistogram = InitializeRandomHist(yhistogram);

  return yhistogram;

}










