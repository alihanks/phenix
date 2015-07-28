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

TH1F * InitializeExpPt(double PtMin, double PtMax, 
		       double Mass, double InvSlope)
{
  int    nbins, ibin;
  float  binwidth;
  float  weight;
  double Pt, Mt;

  nbins = 1000;
  binwidth = (PtMax-PtMin)/(double)nbins;
  TH1F * pthistogram = new TH1F("pt","pt",nbins,PtMin,PtMax);

  for ( ibin=1; ibin<=nbins; ibin++ )
  {
    Pt = PtMin + (ibin-1)*binwidth + binwidth/2.0;
    Mt = sqrt(Mass*Mass+Pt*Pt);
    weight = Pt*exp(-1.0*Mt/InvSlope);
    pthistogram->AddBinContent(ibin,weight);
  }
  
  pthistogram = InitializeRandomHist(pthistogram);

  return pthistogram;

}





