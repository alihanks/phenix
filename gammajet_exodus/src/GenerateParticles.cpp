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

void GenerateParticles(ParticleList *PList, int events, 
		       ParticleGeneratorList *PGList,
		       ParticlePropertyList *PPList)
{

  double pt, y, mass, phi;
  double mt, E,px,py, pz, weight, decaysum;
  double  sample_pt = 0.0;
  double  sample_y = 0.0;
  int    generation = 1;
  int    ievent, iGenerator, itotal;
  int    GeneratorID;
  double GeneratorWeight;
  int    ptbin;
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
    for ( ievent=1; ievent<=events; ievent++ )
    {
      itotal++;

      //cout<<" genid "<<GeneratorID<<" Genptyfun "<<GeneratorPtYFunc<<endl;
      if ( fmod((double)itotal,100000.0)==0.0 )
	cout << itotal << " particles generated" << endl;

      //Run 5 p+p parameterization:
      double gmean= 1.2;
      double gwidth = 34.32;
	
      //Run 7 A+A parameterization:
      //double gmean= -0.13;
      //double gwidth = 23.8;


      double zVTX  = gRandom->Gaus(gmean,gwidth);

      while (fabs(zVTX)>30.0){
	zVTX  = gRandom->Gaus(gmean, gwidth);
      }

      if ( !GeneratorPtYFunc )
      {
        if ( GeneratorID==-111 ) {
	  //cout << "there was a direct photon in the generators" << endl;
	  pt = 30.0*gRandom->Rndm();
	}
	else {
	  //pt = 4.0 + 7.5*gRandom->Rndm();
	  pt = 30.0*gRandom->Rndm();
	}
        y  = GeneratorYHist->GetRandom();
      }
      else
      {
	//cout << "JF JF JF JF using gen functions" << endl;
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
      ptbin = (int)(1000.0*pt) + 1;
      weight   = Generator->GetWeight()*GeneratorPtHist->GetBinContent(ptbin);
      decaysum = 0.0;
      PParticle = new Particle;
      PParticle->SetID(GeneratorID);
      PParticle->Set4mom(E,px,py,pz);
      PParticle->SetDecaysum(decaysum);
      PParticle->SetGeneration(generation);
      PParticle->SetVertex(0.0,0.0,zVTX);
      PParticle->SetWeight(weight);
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






















