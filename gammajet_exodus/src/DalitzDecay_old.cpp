//-----------------------------------------------------------------------------
//
//  Generate Dalitz decay
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <stdlib.h>
#include <math.h>

#include <TH1.h>
#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticlePropertyList.h"
#include "DecayList.h"

#define SQR(a) ((a)*(a))
#define PI  3.141592653589793238

double Cos2(double *, double *);
Mom3 invert(Mom3);
double thetaof(Mom3);
double abs(Mom4);
Mom4 boost(Mom4, Mom4);
Mom3 Rotate(Mom3, double, double, double, double);
double GetMass(int, ParticlePropertyList *);

void DalitzDecay(Particle *PParent, Particle *PLepton1, Particle *PLepton2,
		 Particle *POther, ParticlePropertyList *PPList, Decay *PDecay)
{
  double pmass, lmass, omass, lpmass;
  double E1, p1, E3, p3;
  double costheta, sintheta, cosphi, sinphi, phi;

  int new_generation;

  pmass=GetMass(PDecay->GetParentID(),PPList);
  lmass=0;
  omass=0;
  int ibody=0;
  for( ibody=1; ibody<=3; ibody++)
  {
    int ID;
    ID = PDecay->GetChildID(ibody);
    if ( abs(ID)==11 || abs(ID)==13 )
      lmass = GetMass(PDecay->GetChildID(ibody), PPList);
    else
      omass = GetMass(PDecay->GetChildID(ibody), PPList);
  }
 
  TH1F * LeptonPairMass = PDecay->GetHistogram();

  for ( ibody=1; ibody<=3; ibody++ )
  { 
    switch( PDecay->GetChildID(ibody) )
    {
       case -11:   PLepton1->SetID(-11); break;
       case  11:   PLepton2->SetID(11); break;
       case -13:   PLepton1->SetID(-13); break;
       case  13:   PLepton2->SetID(13); break;
       default :   POther->SetID(PDecay->GetChildID(ibody));
    }
  }

  for( ;; ) 
  {
    lpmass = LeptonPairMass->GetRandom();
    if ( pmass-omass>lpmass && lpmass/2.>lmass ) break;
  }

  E1 = lpmass/2.;
  p1 = sqrt((E1+lmass)*(E1-lmass));
  if ( omass<0.01 )
  {
    do
    {
      costheta = (2.0*gRandom->Rndm())-1.;
    }
    while ( (1.0+costheta*costheta)<(2.0*gRandom->Rndm()) );
  }
  else
  {
    costheta = (2.0*gRandom->Rndm())-1.;
  }
  sintheta = sqrt((1.+costheta)*(1.-costheta));
  phi      = 2.0*PI*gRandom->Rndm();
  sinphi   = sin(phi);
  cosphi   = cos(phi);

  new_generation = PParent->GetGeneration()+1;
  PLepton1->SetGeneration(new_generation);
  PLepton2->SetGeneration(new_generation);
  POther->SetGeneration(new_generation);

  PLepton1->SetDecaysum(0.0);
  PLepton2->SetDecaysum(0.0);
  POther->SetDecaysum(0.0);

  PLepton1->Set4mom(E1,
		    p1*sintheta*cosphi,
		    p1*sintheta*sinphi,
		    p1*costheta);
  PLepton2->Set4mom(E1,invert(PLepton1->Get4mom().Getp()));
  
  E3       = (SQR(pmass)+SQR(omass)-SQR(lpmass))/(2.*pmass);
  p3       = sqrt((E3+omass)*(E3-omass));
  costheta = (2.0*gRandom->Rndm())-1.;
  sintheta = sqrt((1.+costheta)*(1.-costheta));
  phi      = 2.0*PI*gRandom->Rndm();
  sinphi   = sin(phi);
  cosphi   = cos(phi);	    

  POther->Set4mom(E3,
		  p3*sintheta*cosphi,
		  p3*sintheta*sinphi,
		  p3*costheta);

  PLepton1->Set4mom(E1,
    Rotate(PLepton1->Get4mom().Getp(),costheta,-sintheta,-cosphi,-sinphi));
  PLepton2->Set4mom(E1,
    Rotate(PLepton2->Get4mom().Getp(),costheta,-sintheta,-cosphi,-sinphi));
		    
  Mom4 lp_boost(sqrt(SQR(p3)+SQR(lpmass)),invert(POther->Get4mom().Getp()));
  PLepton1->Set4mom(boost(PLepton1->Get4mom(),lp_boost));
  PLepton2->Set4mom(boost(PLepton2->Get4mom(),lp_boost));
  
  PLepton1->Set4mom(boost(PLepton1->Get4mom(),PParent->Get4mom()));
  PLepton2->Set4mom(boost(PLepton2->Get4mom(),PParent->Get4mom()));
  POther->Set4mom(boost(POther->Get4mom(),PParent->Get4mom()));
  
  return;

}









