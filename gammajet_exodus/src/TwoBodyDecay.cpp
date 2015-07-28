//-----------------------------------------------------------------------------
//
//  Generate two-body decay
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <stdlib.h>
#include <math.h>

#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticlePropertyList.h"
#include "DecayList.h"

#define SQR(a) ((a)*(a))
#define PI  3.141592653589793238

Mom3 invert(Mom3);
Mom4 boost(Mom4, Mom4);

double GetMass(int,ParticlePropertyList *);
double GetWidth(int,ParticlePropertyList *);
double thetaof(Mom3);
double phiof(Mom3);
Mom3   Rotate(Mom3, double,double,double,double);

void TwoBodyDecay(Particle *PParent, Particle *PChild1, Particle *PChild2,
			    ParticlePropertyList * PPList, Decay *PDecay)
{
  Mom4   mom4;
  double wp, mp, md1, md2, Ed1, Ed2, pd1, pd2;
  double costheta, sintheta, theta, cosphi,sinphi, phi;

  int new_generation;

  wp  = GetWidth(PDecay->GetParentID(),PPList);
  if ( wp==0.0 )
    {
      mp = GetMass(PDecay->GetParentID(),PPList);
    }
  else
    {
      mp = sqrt(SQR(PParent->Get4mom()));
    }
  md1 = GetMass(PDecay->GetChildID(1),PPList);
  md2 = GetMass(PDecay->GetChildID(2),PPList);

  if ( mp<md1+md2 )
  {
    //cout << "Decay kinematically impossible for particle " 
    //	 << PDecay->GetParentID() << "!" << endl;
    return;
  }

  Ed1 = (SQR(mp)+SQR(md1)-SQR(md2))/(2.*mp);
  Ed2 = mp-Ed1;
  pd1 = sqrt((Ed1+md1)*(Ed1-md1));
  pd2 = sqrt((SQR(mp)-SQR(md1+md2))*(SQR(mp)-SQR(md1-md2)))/(2.*mp);
    
  costheta = (2.0*gRandom->Rndm())-1.0;
  sintheta = sqrt((1.+costheta)*(1.-costheta));
  phi      = 2.0*PI*gRandom->Rndm();

  PChild1->SetID(PDecay->GetChildID(1));
  PChild2->SetID(PDecay->GetChildID(2));

  new_generation = PParent->GetGeneration()+1;
  PChild1->SetGeneration(new_generation);
  PChild2->SetGeneration(new_generation);

  PChild1->SetDecaysum(0.0);
  PChild2->SetDecaysum(0.0);

  PChild1->Set4mom(Ed1,
		   pd1*sintheta*cos(phi), 
		   pd1*sintheta*sin(phi), 
		   pd1*costheta);
  PChild2->Set4mom(Ed2,invert(PChild1->Get4mom().Getp()));

  theta    = thetaof(PParent->Get4mom().Getp());
  phi      = phiof(PParent->Get4mom().Getp());
  costheta = cos(theta); 
  sintheta = sin(theta);
  cosphi   = cos(phi); 
  sinphi   = sin(phi);
  
  PChild1->Set4mom(Ed1,
    Rotate(PChild1->Get4mom().Getp(),costheta,sintheta,cosphi,sinphi));
  PChild2->Set4mom(Ed2,
    Rotate(PChild2->Get4mom().Getp(),costheta,sintheta,cosphi,sinphi));
  
  PChild1->Set4mom(boost(PChild1->Get4mom(),PParent->Get4mom()));
  PChild2->Set4mom(boost(PChild2->Get4mom(),PParent->Get4mom()));
  
  return;

}
  



    



