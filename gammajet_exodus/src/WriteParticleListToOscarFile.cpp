//-----------------------------------------------------------------------------
//
//  Write a stream of particles in an OSCAR compliant ASCII file 
//
//-----------------------------------------------------------------------------

#include <fstream>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticlePropertyList.h"

void WriteParticleListToOscarFile(char * file, ParticleList * PList, 
				  ParticlePropertyList * PPList)
{
  PLNode   * CurrentNode     = PList->GetHeadNode();
  Particle * CurrentParticle = 0;
  ParticleProperty * Species = 0;
  int    pid; 
  Mom4   four_momentum;
  double E, px, py, pz, mass;
  double xvtx, yvtx, zvtx;

  cout << "Writing particle list to OSCAR compliant file: " 
       << file << endl;

  ofstream fout(file);

  fout << "# OSC1999A" << endl;
  fout << "# final_id_p_x" << endl;
  fout << "# EXODUS event generator in single particle mode" << endl;
  fout << "#" << endl;

  fout.precision(6);
  for ( int i=1; i<=PList->GetLength(); i++) 
  {
    CurrentNode     = CurrentNode->GetNextNode();
    CurrentParticle = CurrentNode->Get(0);
    pid             = CurrentParticle->GetID();
    Species         = PPList->GetByID(pid);
    mass            = Species->GetMass();
    four_momentum   = CurrentParticle->Get4mom();
    E               = four_momentum.GetE();
    px              = four_momentum.Getp().Getpx();
    py              = four_momentum.Getp().Getpy();
    pz              = four_momentum.Getp().Getpz();
    xvtx            = CurrentParticle->GetxVertex();
    yvtx            = CurrentParticle->GetyVertex();
    zvtx            = CurrentParticle->GetzVertex();

    if ( CurrentParticle->GetAccept()>=0 )
    {
      fout << "1 0" << endl;
      fout << "0 " << pid << " 0 ";
      fout << px   << " ";
      fout << py   << " ";
      fout << pz   << " ";
      fout << E    << " ";
      fout << mass << " ";
      fout << xvtx << " ";
      fout << yvtx << " ";
      fout << zvtx << " ";
      fout << "0"  << endl;
      fout << "0 0" << endl;
    }

  }

  fout.close();

  return;

}





