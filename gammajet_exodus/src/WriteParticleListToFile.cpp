//-----------------------------------------------------------------------------
//
//  Write a stream of particles in an ASCII file 
//
//-----------------------------------------------------------------------------

#include <fstream>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"

void WriteParticleListToFile(char * file, ParticleList * PList)
{
  PLNode   * CurrentNode     = PList->GetHeadNode();
  Particle * CurrentParticle = 0;
  int    pid; 
  Mom4   four_momentum;
  double E, px, py, pz, weight;

  cout << "Writing particle list to file: " 
       << file << endl;

  ofstream fout(file);

  fout.precision(6);
  fout << endl;
  for ( int i=1; i<=PList->GetLength(); i++) 
  {
    CurrentNode     = CurrentNode->GetNextNode();
    CurrentParticle = CurrentNode->Get(0);
    pid             = CurrentParticle->GetID();
    four_momentum   = CurrentParticle->Get4mom();
    E               = four_momentum.GetE();
    px              = four_momentum.Getp().Getpx();
    py              = four_momentum.Getp().Getpy();
    pz              = four_momentum.Getp().Getpz();
    weight          = CurrentParticle->GetWeight();
    fout << pid    << " ";
    fout.precision(12);
    fout << E      << " ";
    fout << px     << " ";
    fout << py     << " ";
    fout << pz     << " ";
    fout.precision(6);
    fout << weight << " ";
    fout << endl;
  }

  fout.close();

  return;

}





