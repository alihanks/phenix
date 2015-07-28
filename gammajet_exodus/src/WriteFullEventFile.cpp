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

void WriteFullEventFile(ofstream * fout, int eventID, ParticleList * PList, 
			ParticlePropertyList * PPList)
{
  ParticleProperty * Species = 0;
  int    pid, parent; 
  int    final_state_particles = 0;
  int    particle_in_file = 0;
  Mom4   four_momentum;
  double E, px, py, pz, mass, weight;
  double xvtx, yvtx, zvtx;

  PLNode   * CurrentNode     = PList->GetHeadNode();
  PLNode   * NextNode        = 0;
  Particle * CurrentParticle = 0;
  Particle * NextParticle    = 0;
  Particle * CurrentPrimary  = 0;

  CurrentNode     = CurrentNode->GetNextNode();
  CurrentParticle = CurrentNode->Get(0);
  if ( PList->GetLength()>1 )
  {
    NextNode        = CurrentNode->GetNextNode();
    NextParticle    = NextNode->Get(0);
  }
  CurrentPrimary  = CurrentParticle;
  for ( int i=1; i<=PList->GetLength(); i++) 
  {
    bool final_state = false;
    if ( i != PList->GetLength() ) {
      if ( CurrentParticle->GetGeneration() >= NextParticle->GetGeneration() ) 
	final_state = true;
    }
    else {
      final_state = true;
    }
    if ( final_state ) final_state_particles++;
    CurrentNode     = NextNode;
    CurrentParticle = NextParticle;
    if ( i<PList->GetLength()-1 )
    {
      NextNode = CurrentNode->GetNextNode();
      NextParticle    = NextNode->Get(0);
    }  
  }

  CurrentNode     = PList->GetHeadNode();
  NextNode        = 0;
  CurrentParticle = 0;
  NextParticle    = 0;
  CurrentPrimary  = 0;

  CurrentNode     = CurrentNode->GetNextNode();
  CurrentParticle = CurrentNode->Get(0);
  if ( PList->GetLength()>1 )
  {
    NextNode        = CurrentNode->GetNextNode();
    NextParticle    = NextNode->Get(0);
  }
  CurrentPrimary  = CurrentParticle;

  *fout << "2 " << final_state_particles << endl;
  for ( int i=1; i<=PList->GetLength(); i++) 
  {
    if ( CurrentParticle->GetGeneration()==1 )
    {
      CurrentPrimary = CurrentParticle;
      pid            = CurrentParticle->GetID();
    }
    bool final_state = false;
    if ( i != PList->GetLength() ) {
      if ( CurrentParticle->GetGeneration() >= NextParticle->GetGeneration() ) 
	final_state = true;
    }
    else {
      final_state = true;
    }
    parent = 0;
    if ( final_state ) parent = CurrentPrimary->GetID();
    if ( final_state )
    {
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
      weight          = CurrentParticle->GetWeight();
      particle_in_file++;
      *fout << particle_in_file << " "; 
      *fout << pid              << " ";
      *fout << eventID          << " ";
      *fout << px               << " ";
      *fout << py               << " ";
      *fout << pz               << " ";
      *fout << E                << " ";
      *fout << mass             << " ";
      *fout << xvtx             << " ";
      *fout << yvtx             << " ";
      *fout << zvtx             << " ";
      *fout << "0"              << " ";
      *fout << weight           << " ";
      *fout << parent           << endl;
    }
    CurrentNode     = NextNode;
    CurrentParticle = NextParticle;
    if ( i<PList->GetLength()-1 )
    {
      NextNode = CurrentNode->GetNextNode();
      NextParticle    = NextNode->Get(0);
    }  
  }
  *fout << "0 0 " << endl;

  return;

}





