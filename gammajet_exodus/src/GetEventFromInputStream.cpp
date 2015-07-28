//-----------------------------------------------------------------------------
//
//  Read a stream of particles from an ASCII input file and insert
//  them into a ParticleList
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <math.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"

ParticleList * GetEventFromInputStream(ifstream * input_file)
{
  ParticleList * ReadParticleList = 0; 
  PLNode       * CurrentNode      = 0;

  char   character;
  int    pid; 
  double E, px, py, pz, weight;

  ReadParticleList = new ParticleList;
  CurrentNode      = ReadParticleList->GetHeadNode();

  if ( input_file->get(character) ) 
  {
    pid    = 0;
    E      = 0.0;
    px     = 0.0;
    py     = 0.0;
    pz     = 0.0;
    weight = 0.0;
    *input_file >> pid;
    *input_file >> E;
    *input_file >> px;
    *input_file >> py;
    *input_file >> pz;
    *input_file >> weight;
    Particle * pParticle = new Particle;
    pParticle->Set4mom(E,px,py,pz);
    pParticle->SetID(pid);
    pParticle->SetWeight(weight);
    pParticle->SetDecaysum(0.0);
    pParticle->SetGeneration(1);
    if ( E>0.0 ) 
    { 
      ReadParticleList->InsertAfter(CurrentNode, pParticle);
      CurrentNode = CurrentNode->GetNextNode();
     }
  }

  if ( E==0.0 ) ReadParticleList=0;  
  return ReadParticleList;

}





