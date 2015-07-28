//-----------------------------------------------------------------------------
//
//  Check if lepton pair is in the CERES acceptance
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>

#include "Momentum.h"

#define PI  3.141592653589793238

double thetaof(Mom3);
double phiof(Mom3);

bool InCERESAcceptance(Mom4 mom1, Mom4 mom2)
{
  double m1, E1, p1, px1, py1, pz1, pt1, y1;
  double m2, E2, p2, px2, py2, pz2, pt2, y2;
  double angle_12;
  bool   result;

  result = true;

  m1  = sqrt(mom1*mom1);
  E1  = mom1.GetE();
  p1  = sqrt(E1*E1-m1*m1);
  px1 = mom1.Getp().Getpx();
  py1 = mom1.Getp().Getpy();
  pz1 = mom1.Getp().Getpz();
  pt1 = sqrt(px1*px1+py1*py1);
  y1  = log((E1+pz1)/(E1-pz1))/2.0;

  m2  = sqrt(mom2*mom2);
  E2  = mom2.GetE();
  p2  = sqrt(E2*E2-m2*m2);
  px2 = mom2.Getp().Getpx();
  py2 = mom2.Getp().Getpy();
  pz2 = mom2.Getp().Getpz();
  pt2 = sqrt(px2*px2+py2*py2);
  y2  = log((E2+pz2)/(E2-pz2))/2.0;

  angle_12 = acos((mom1.Getp()*mom2.Getp())/(p1*p2));

  if ( pt1<0.05 )           result=false;
  if ( pt2<0.05 )           result=false;
  if ( y1<2.10 || y1>2.65 ) result=false;
  if ( y2<2.10 || y2>2.65 ) result=false;
  if ( angle_12<0.035 )     result=false;

  return (result);

}
