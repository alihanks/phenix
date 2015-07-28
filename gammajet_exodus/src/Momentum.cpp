//-----------------------------------------------------------------------------
//
//  Implementation of the classes Mom3 and Mom4
//
//-----------------------------------------------------------------------------

#include <iostream>
#include "Momentum.h"

Mom3::Mom3(){}

Mom3::Mom3(double initpx, double initpy, double initpz)
{
  itsPx = initpx; 
  itsPy = initpy; 
  itsPz = initpz;
}

Mom3::~Mom3(){}

void Mom3::Set(double px, double py, double pz)
{
  itsPx = px;
  itsPy = py;
  itsPz = pz;
}

Mom3 Mom3::operator+ (const Mom3 & rhs)
{
  return Mom3( itsPx+rhs.Getpx(), itsPy+rhs.Getpy(), itsPz+rhs.Getpz() );
}

Mom3 Mom3::operator- (const Mom3 & rhs)
{
  return Mom3( itsPx-rhs.Getpx(), itsPy-rhs.Getpy(), itsPz-rhs.Getpz() );
}

double Mom3::operator* (const Mom3 & rhs)
{
  return itsPx*rhs.Getpx() + itsPy*rhs.Getpy() + itsPz*rhs.Getpz();
}

Mom3 Mom3::operator* (const double & rhs)
{
  return Mom3( itsPx*rhs, itsPy*rhs, itsPz*rhs);
}

Mom4::Mom4(){}

Mom4::Mom4(double initE, Mom3 initp)
{
  itsE = initE;
  itsP = initp;
}

Mom4::Mom4(double initE, double initpx, double initpy, double initpz)
{
  itsE = initE;
  itsP.Setpx(initpx);
  itsP.Setpy(initpy);
  itsP.Setpz(initpz);
}

Mom4::~Mom4(){}

Mom4 Mom4::operator+ (const Mom4 & rhs)
{
  return Mom4( itsE+rhs.GetE(), itsP+rhs.Getp());
}

double Mom4::operator* (const Mom4&rhs)
{
  return itsE*rhs.GetE() - itsP*rhs.Getp();
}















