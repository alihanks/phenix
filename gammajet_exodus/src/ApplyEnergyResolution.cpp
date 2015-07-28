 //-----------------------------------------------------------------------------
//
//  Apply PHENIX resolution to a four momentum
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <TH1.h>
#include <TRandom.h>
#include "Momentum.h"

#define PI  3.141592653589793238

double thetaof(Mom3);
double phiof(Mom3);

using namespace std;

Mom4 ApplyEnergyResolution(Mom4 mom4, int pbsc_pbgl)
{
  double p, px, py, pz, E, theta, phi;
  double sigma=0;
  Mom4   result;

  px = mom4.Getp().Getpx();
  py = mom4.Getp().Getpy();
  pz = mom4.Getp().Getpz();
  p  = sqrt(px*px+py*py+pz*pz);
  E  = mom4.GetE();

  theta = thetaof(mom4.Getp());
  phi   = phiof(mom4.Getp());

  /*
  if(pbsc_pbgl==0){ 
    cout << "applying_pbsc_resolution" <<endl;
  }else{
    cout << "applying_pbgl_resolution" <<endl;
  }
  */

  if(pbsc_pbgl==0)  
    //    sigma = E*sqrt((0.035)*(0.035)+(0.081/sqrt(E))*(0.081/sqrt(E))); //NIM values w/ test beam
  //  sigma = E*sqrt((0.021)*(0.021)+(0.081/sqrt(E))*(0.081/sqrt(E))); //NIM values w/ test beam

  //if(pbsc_pbgl==0)
  sigma = E*sqrt((0.06)*(0.06)+(0.081/sqrt(E))*(0.081/sqrt(E))); //more realistic
  //if(pbsc_pbgl==0)sigma = E*sqrt((0.085)*(0.085)+(0.081/sqrt(E))*(0.081/sqrt(E))); 

  //****For the pbgl switched to newest resolution used by Matt on 12/09/08
  //if(pbsc_pbgl==1)sigma = E*sqrt((0.008)*(0.008)+(0.059/sqrt(E))*(0.059/sqrt(E))); //latest and greatest? -man
 //AN647
  if(pbsc_pbgl==1)
    sigma = E*sqrt((0.05)*(0.05)+(0.06/sqrt(E))*(0.06/sqrt(E))); 
  
  //  if(pbsc_pbgl==1)sigma = E*sqrt((0.06)*(0.06)+(0.09/sqrt(E))*(0.09/sqrt(E))); 


  // cout<<" pbsc_pbgl "<<pbsc_pbgl<<" sigma "<<sigma<<endl;
  E     = gRandom->Gaus(E,sigma);

  //E = 1.02*E; //shift energy up by 2%
  //cout << "Applied Energy Shift" <<endl;

  px    = E*sin(theta)*cos(phi);
  py    = E*sin(theta)*sin(phi);
  pz    = E*cos(theta);

  //  px    = px*E/p;;
  //  py    = py*E/p;
  //  pz    = pz*E/p;

  result.SetE(E);
  result.Setp(px,py,pz);

  return (result);

}

