//-----------------------------------------------------------------------------
//
// Generate white transverse-momentum distribution
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>

#include <TH1.h>
#include <TRandom.h>

TH1F * InitializeRandomHist(TH1F*);

TH1F * InitializeWhitePt(double PtMin, double PtMax)
{
  int    nbins, ibin;
  float  binwidth;
  float  weight;

  nbins = 1000;
  binwidth = (PtMax-PtMin)/(double)nbins;
  TH1F * pthistogram = new TH1F("pt","pt",nbins,PtMin,PtMax);

  for ( ibin=1; ibin<=nbins; ibin++ )
  {
    weight = 1.0;
    pthistogram->AddBinContent(ibin,weight);
  }
  
  pthistogram = InitializeRandomHist(pthistogram);

  return pthistogram;

}





