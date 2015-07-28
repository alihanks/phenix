//-----------------------------------------------------------------------------
//
//  Generate three-body decay (except Dalitz decay)
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>

#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "DecayList.h"
#include "ParticlePropertyList.h"

#define SQR(a) ((a)*(a))
#define PI  3.141592653589793238

Mom3 invert(Mom3);
double abs(Mom4);
Mom4 boost(Mom4, Mom4);
ParticleList * ThreeBodyDecay(Particle *,ParticlePropertyList *,Decay *);
double GetMass(int,ParticlePropertyList *);
double thetaof(Mom3);
double phiof(Mom3);
Mom3 Rotate(Mom3, double,double);

void ThreeBodyDecay(Particle *PParent, Particle *PChild1, Particle *PChild2,
		    Particle *PChild3, ParticlePropertyList *PPList,
		    Decay *PDecay)
{
  int i;
  double mp, md[3], Ed1, Ed2, Ed3, pd1, pd3, E12;
  double m12, m13, t1, e1s, e3s;
  double m12lo, m12hi, m13lo, m13hi, m13max, m13min;
  double costheta, sintheta, phi, sinphi, cosphi;

  int new_generation;
 
  mp=GetMass(PParent->GetID(),PPList); 
  
  for ( i=1; i<=3; i++)
  {
    md[i-1] = GetMass(PDecay->GetChildID(i),PPList);
  }

  m12lo = SQR(md[0]+md[1]);
  m12hi = SQR(mp-md[2]);
  m13lo = SQR(md[0]+md[2]);
  m13hi = SQR(mp-md[1]);
  
  for ( ; ; )
  {
    m12 = sqrt(m12lo+(m12hi-m12lo)*gRandom->Rndm());
    m13 = sqrt(m13lo+(m13hi-m13lo)*gRandom->Rndm());
    e1s = (m12*m12+md[0]*md[0]-md[1]*md[1])/(2.0*m12);
    e3s = (mp*mp-m12*m12-md[2]*md[2])/(2.0*m12);
      
    t1 = e1s*e1s-md[0]*md[0];
    if ( t1<0 ) continue;
      
    m13max = sqrt(SQR(e1s+e3s)-SQR(sqrt(t1)-sqrt(SQR(e3s)-SQR(md[2]))));
    m13min = sqrt(SQR(e1s+e3s)-SQR(sqrt(t1)+sqrt(SQR(e3s)-SQR(md[2]))));
      
    if( m13<=m13max && m13>=m13min ) break;
  }
  
  Ed3 = (SQR(mp)+SQR(md[2])-SQR(m12))/(2.0*mp);
  pd3 = sqrt((Ed3+md[2])*(Ed3-md[2]));
  
  costheta = 2.0*gRandom->Rndm()-1.;
  sintheta = sqrt((1.+costheta)*(1.-costheta));
  phi      = 2.0*PI*gRandom->Rndm();
  sinphi   = sin(phi);
  cosphi   = cos(phi);

  PChild1->SetID(PDecay->GetChildID(1));
  PChild2->SetID(PDecay->GetChildID(2));
  PChild3->SetID(PDecay->GetChildID(3));
  
  new_generation = PParent->GetGeneration()+1;
  PChild1->SetGeneration(new_generation);
  PChild2->SetGeneration(new_generation);
  PChild3->SetGeneration(new_generation);

  PChild1->SetDecaysum(0.0);
  PChild2->SetDecaysum(0.0);
  PChild3->SetDecaysum(0.0);

  PChild3->Set4mom(Ed3,pd3*sintheta*cosphi,pd3*sintheta*sinphi,pd3*costheta);
   
  Ed1      = (m12*m12+md[0]*md[0]-md[1]*md[1])/(2.*m12);
  Ed2      = m12-Ed1;
  pd1      = sqrt((Ed1+md[0])*(Ed1-md[0]));
  costheta = 2.0*gRandom->Rndm()-1.0;
  sintheta = sqrt((1.-costheta)*(1.+costheta));
  phi      = 2.0*PI*gRandom->Rndm();
  
  PChild1->Set4mom(Ed1,
		   pd1*sintheta*cos(phi), 
		   pd1*sintheta*sin(phi), 
		   pd1*costheta);
  PChild2->Set4mom(sqrt(SQR(PChild1->Get4mom().Getp())+SQR(md[1])),
		   invert(PChild1->Get4mom().Getp()));
    
  E12 = sqrt(SQR(pd3)+SQR(m12));
  Mom4 p4boost(E12,invert(PChild3->Get4mom().Getp()));

  PChild1->Set4mom(boost(PChild1->Get4mom(),p4boost));
  PChild2->Set4mom(boost(PChild2->Get4mom(),p4boost));

  PChild1->Set4mom(boost(PChild1->Get4mom(),PParent->Get4mom()));
  PChild2->Set4mom(boost(PChild2->Get4mom(),PParent->Get4mom()));
  PChild3->Set4mom(boost(PChild3->Get4mom(),PParent->Get4mom()));

  return;

}
















