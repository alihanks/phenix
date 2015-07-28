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

ParticleList * ReadParticleListFromFile(char * file)
{
  ParticleList * ReadParticleList = 0; 
  PLNode       * CurrentNode      = 0;

  char   character;
  int    pid, itotal; 
  double E, px, py, pz, weight;

  itotal = 0;
  ReadParticleList = new ParticleList;
  CurrentNode      = ReadParticleList->GetHeadNode();

  cout << "Reading particle list from file: " 
       << file << endl;
  cout << endl;

  ifstream fin(file);
  while (fin.get(character)) 
  {
    fin >> pid;
    fin >> E;
    fin >> px;
    fin >> py;
    fin >> pz;
    fin >> weight;
    Particle * pParticle = new Particle;
    pParticle->Set4mom(E,px,py,pz);
    pParticle->SetID(pid);
    pParticle->SetWeight(weight);
    pParticle->SetDecaysum(0.0);
    pParticle->SetGeneration(1);
    if ( E>0.0 ) 
    { 
      itotal++;
      if ( fmod((double)itotal,10000.0)==0.0 )
	cout << itotal << " particles read from file" << endl;
      ReadParticleList->InsertAfter(CurrentNode, pParticle);
      CurrentNode = CurrentNode->GetNextNode();
     }
    E = 0.0;
  }

  fin.close();

  cout << endl;

  return ReadParticleList;
}





