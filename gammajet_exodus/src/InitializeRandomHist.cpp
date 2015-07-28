#include <stdlib.h>
#include <math.h>

#include <TH1.h>
#include <TRandom.h>

TH1F * InitializeRandomHist(TH1F *histogram)
{
  for ( int i=1; i<=(int)(100000.0*gRandom->Rndm()+1.0) ; i++ )
  {
    histogram->GetRandom();
  }
  
  return histogram;
}
