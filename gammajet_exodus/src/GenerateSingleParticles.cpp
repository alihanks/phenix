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

#define PI  3.141592653589793238

TH1F * InitializeRandomHist(TH1F *);
int    PHENIXFilter(Particle*, ParticlePropertyList*);

void GenerateSingleParticles(ParticleList *PList, int events, 
		             ParticlePropertyList *PPList)
{
  int    pID, setup;
  int    sector;
  double mass, phi, zVTXmax, zVTX;
  double ptmax, pt, ymin, ymax, y;
  double phimin, phimax;
  double mt, E,px,py, pz;
  int    nbins;
  double min, max, binwidth, weight;
  int    ievent, itotal;

  Particle         * PParticle = 0;
  ParticleProperty * PProperty = 0;
  ParticleProperty * Property = 0;
  PLNode           * Current   = PList->GetHeadNode();
  
  do
  {
    cout << "Choose the shape of the pt distribution:" << endl;
    cout << "----------------------------------------" << endl;
    cout << endl;
    cout << "1) flat" << endl;
    cout << "2) low pt enhancement" << endl;
    cout << endl;
    cout << "Your choice (1-2): "; 
    cin  >> setup; 
    cout << endl;
  } while ( setup<1 || setup>2 );

  cout << "Select particle by ID from the following list:" << endl << endl;

  int iProperty = 1;
  Property  = PPList->Get(iProperty);
  while ( Property->GetID()!=0 )
    {
      cout << "ID: " << Property->GetID() << " is a " << Property->GetName() << endl;
      iProperty++;
      Property = PPList->Get(iProperty);
    }
  cout << endl;
  delete Property;
  Property = 0;

  cout << "ID of particle to generate: ";
  cin  >> pID;
  cout << endl;
  cout<<" over here in generate singles "<<endl;
  PProperty        = PPList->GetByID(pID);

  mass             = PProperty->GetMass();
  ptmax            = 10.0;
  zVTXmax          = 20.0;
  ymin             = -0.6;
  ymax             = 0.6;
  // TOF: 1.92<phi<4.36
  //  phimin           = 1.92;
  //  phimax           = 4.36;
  // TOF (new):
  phimin           = PI-(80./180.)*PI;
  phimax           = PI+(80./180.)*PI;
  // EMC: -pi/2<phi<pi/2
  //  phimin           = -PI/2.;
  //  phimax           =  PI/2.;
  itotal           = 0;

  nbins = 10000;
  min = 0.0;
  max = 8.0;
  binwidth = (max-min)/(double)nbins;
  TH1F * pthistogram = new TH1F("pt","pt",nbins,min,max);

  for ( int ibin=1; ibin<=nbins; ibin++ )
  {
    pt = min+(double)(ibin-1)*binwidth+binwidth/2.0;
    weight = 1.0;
    if ( setup==2 )
    {
      weight = 1.0 + 10.0*exp(-2.0*pt);
    }
    pthistogram->AddBinContent(ibin,weight);
  }
  pthistogram = InitializeRandomHist(pthistogram);

  for ( ievent=1; ievent<=events; ievent++ )
  {
    itotal++;
    if ( fmod((double)itotal,1000.0)==0.0 )
      cout << itotal << " particles generated" << endl;
    //    phi      = 2.0*PI*gRandom->Rndm();
    phi      = phimin+(phimax-phimin)*gRandom->Rndm();
    if ( phi<0. ) phi = 2.0*PI+phi;

    //    pt       = ptmax*gRandom->Rndm();
    pt       = pthistogram->GetRandom(); 

    zVTX     = 2.0*zVTXmax*(gRandom->Rndm()-0.5);
    y        = ymin+(ymax-ymin)*gRandom->Rndm();

    zVTX = 0.0;
    y    = 0.0;
    pt   = 10.0;
    phi  = 2.0*PI*gRandom->Rndm();

    mt       = sqrt(pt*pt+mass*mass);
    E        = mt*cosh(y);
    px       = pt*cos(phi);
    py       = pt*sin(phi);
    pz       = mt*sinh(y);

    PParticle = new Particle;
    PParticle->SetID(pID);
    PParticle->Set4mom(E,px,py,pz);
    PParticle->SetDecaysum(0.0);
    PParticle->SetWeight(1.0);
    PParticle->SetGeneration(1);
    PParticle->SetVertex(0.0,0.0,zVTX);
    sector = PHENIXFilter(PParticle,PPList);
    PParticle->SetAccept(sector);
    PList->InsertAfter(Current, PParticle);
    Current = Current->GetNextNode();
  }

  return;
}
