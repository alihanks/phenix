#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticlePropertyList.h"
#include "ParticleGeneratorList.h"

#define PI  3.141592653589793238

void GenerateFullEvent(ParticleList *PList, int dnch_dy, 
		       ParticleGeneratorList *PGList,
		       ParticlePropertyList *PPList)
{
  double pt, y, mass, phi;
  double mt, E,px,py, pz, weight, decaysum;
  double sample_pt = 0.0;
  double sample_y = 0.0;
  int    generation = 1;
  int    multiplicity = 1;
  int    ievent, iGenerator, itotal;
  int    GeneratorID;
  double GeneratorWeight;
  TF2*   GeneratorPtYFunc;
  TH1F*  GeneratorPtHist;
  TH1F*  GeneratorYHist;
  TH1F*  GeneratorMHist;

  Particle         * PParticle = 0;
  ParticleProperty * PProperty = 0;
  PLNode           * Current   = PList->GetHeadNode();
  
  ParticleGenerator * Generator = 0;
  itotal     = 0;
  iGenerator = 1;
  Generator  = PGList->Get(iGenerator);

  while ( Generator->GetID()!=0 )
  {
    GeneratorID      = Generator->GetID();
    PProperty        = PPList->GetByID(GeneratorID);
    mass             = PProperty->GetMass();
    GeneratorWeight  = Generator->GetWeight();
    GeneratorPtYFunc = Generator->GetPtYFunction();
    GeneratorPtHist  = Generator->GetPtHistogram();
    GeneratorYHist   = Generator->GetYHistogram();
    GeneratorMHist   = Generator->GetMHistogram();

    //To use takao's dn/dy for pi0s we need to remove this factor of 2
    //multiplicity = gRandom->Poisson((double)(2*dnch_dy)*GeneratorWeight);
    multiplicity = gRandom->Poisson((double)(dnch_dy)*GeneratorWeight);

   if ( GeneratorID==113 )
    {
      multiplicity = 1;
      GeneratorWeight = 0.1;
    }
    if ( GeneratorID==223 )
    {
      multiplicity = 1;
      GeneratorWeight = 0.08;
    }
    if ( GeneratorID==333 )
    {
      multiplicity = 1;
      GeneratorWeight = 0.02;
    }
    if ( GeneratorID==443 )
    {
      multiplicity = 1;
      GeneratorWeight = 2.18e-05;
    }
    if ( GeneratorID==553 )
    {
      multiplicity = 1;
      GeneratorWeight = 1.267e-07;
    }
    if ( abs(GeneratorID)==211 ||
	 abs(GeneratorID)==321 ||
	 abs(GeneratorID)==2212 ||
	 GeneratorID==111 ||
	 GeneratorID==221 ||
	 GeneratorID==331 )
    {
      GeneratorWeight = 1.0; //was 1.0 changed to 0.5 on 01/12/09 MEC
      //changed back to 1.0 on 01/27/09
    }


    //Run 5 p+p parameterization:
    //double gmean= 1.2;
    //double gwidth = 34.32;
	
    //Run 7 A+A parameterization:
    double gmean= -0.13;
    double gwidth = 23.8;

    //Added zVTX smearing on 01/16/09 MEC
    double zVTX  = gRandom->Gaus(gmean,gwidth);

    while (fabs(zVTX)>30.0){
      zVTX  = gRandom->Gaus(gmean, gwidth);
    }


    for ( ievent=1; ievent<=multiplicity; ievent++ )
    {
      cout << "on event " << ievent << endl;
      itotal++;
      if ( !GeneratorPtYFunc )
      {

        pt = GeneratorPtHist->GetRandom();
	//pt = gRandom->Rndm()*160.0 + 4.0;
        y  = GeneratorYHist->GetRandom();
      }
      else
      {

	GeneratorPtYFunc->GetRandom2(sample_pt,sample_y);
	pt = sample_pt;
	y  = sample_y;
      }
      if ( GeneratorMHist!=0 ) mass = GeneratorMHist->GetRandom();
      phi      = 2.0*PI*gRandom->Rndm();
      mt       = sqrt(pt*pt+mass*mass);
      E        = mt*cosh(y);
      px       = pt*cos(phi);
      py       = pt*sin(phi);
      pz       = mt*sinh(y);
      weight   = Generator->GetWeight();
      decaysum = 0.0;
      PParticle = new Particle;
      PParticle->SetID(GeneratorID);
      PParticle->Set4mom(E,px,py,pz);
      PParticle->SetDecaysum(decaysum);
      PParticle->SetGeneration(generation);
      PParticle->SetWeight(GeneratorWeight);
      PParticle->SetVertex(0.0,0.0,zVTX);
      PList->InsertAfter(Current, PParticle);
      Current = Current->GetNextNode();
    }
    iGenerator++;
    Generator = PGList->Get(iGenerator);
  }
  delete Generator;
  Generator = 0;

  return;
}






















