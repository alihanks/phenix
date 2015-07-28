//-----------------------------------------------------------------------------
//
//  Write a stream of particles in an ASCII file 
//
//-----------------------------------------------------------------------------

#include <fstream>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"

void WriteEventToOutputStream(ofstream * output_file, ParticleList * PList)
{
  PLNode   * CurrentNode     = PList->GetHeadNode();
  Particle * CurrentParticle = 0;
  int    pid; 
  Mom4   four_momentum;
  double E, px, py, pz, weight;

  output_file->precision(6);

  //  *output_file << endl;
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
    //    *output_file << pid    << " ";
    output_file->precision(10);
    *output_file << E      << " ";
    *output_file << px     << " ";
    *output_file << py     << " ";
    *output_file << pz     << " ";
    //    output_file->precision(6);
    //    *output_file << weight << " ";
    //    *output_file << endl;
  }
    *output_file << endl;

  return;

}





