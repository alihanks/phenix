//-----------------------------------------------------------------------------
//
//  Implementation of the classes Particle and ParticleProperty
//
//-----------------------------------------------------------------------------

#include <iostream>
#include "Momentum.h"
#include "Particle.h"
#include "TMath.h"

Particle::Particle(){}

Particle::Particle(int inittype, Mom4 initmom4, double initdecaysum, 
                   int initGeneration)
{
  itsType       = inittype;
  itsMom4       = initmom4;
  itsDecaysum   = initdecaysum;
  itsGeneration = initGeneration;
}

Particle::Particle(int inittype, double initE, 
                   double initpx, double initpy, double initpz,
                   double initdecaysum, int initGeneration)
{
  itsType       = inittype;
  itsDecaysum   = initdecaysum;
  itsGeneration = initGeneration;
  itsMom4.SetE(initE);
  itsMom4.Setp(initpx,initpy,initpz);
}

Particle::~Particle(){}

void Particle::Show()
{
  cout << "| PID: " << itsType << " |";
  cout.width(10);
  cout<<itsMom4.GetE()<<"|";
  cout.width(10);
  cout<<itsMom4.Getp().Getpx()<<"|";
  cout.width(10);
  cout<<itsMom4.Getp().Getpy()<<"|";
  cout.width(10);
  cout<<itsMom4.Getp().Getpz()<<"|";
  cout.width(4);
  cout<<itsDecaysum<<"|";
  cout.width(10);
  cout<<itsGeneration<<"|"<<endl;
}

long double Particle::GetMass()
{
  if (itsMom4*itsMom4<0)
    {
      cout << "Error in relativistic kinematics: mass < 0" << endl;
      return 0.;
    }
  else
    return TMath::Sqrt(itsMom4*itsMom4);
    //return sqrt(itsMom4*itsMom4);
}

void ParticleProperty::Set(int ID, double mass, double width, 
                           int spin, int charge)
{
  itsID     = ID;
  itsMass   = mass;
  itsWidth  = width;
  itsSpin   = spin;
  itsCharge = charge;
}










