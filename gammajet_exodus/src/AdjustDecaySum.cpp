//-----------------------------------------------------------------------------
//
//  Adjust decay sum of particles in a particle list according to the
//  decays in a decay list
//
//-----------------------------------------------------------------------------

#include <stdlib.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "DecayList.h"

ParticleList * AdjustDecaySum(ParticleList * PList, DecayList * DList)
{
  int i, pid;
  PLNode * CurrentNode = PList->GetHeadNode();
  Particle * CurrentParticle = 0;

  for ( i=1; i<=PList->GetLength(); i++ )
  {
    CurrentNode = CurrentNode->GetNextNode();
    CurrentParticle = CurrentNode->Get(0);
    pid = CurrentParticle->GetID();
    for ( int j=1; j<=DList->GetSize(); j++ )
    {
      Decay * CurrentDecay = DList->Get(j);
      if ( pid==CurrentDecay->GetParentID() )
      {
	CurrentParticle->SetDecaysum(CurrentDecay->GetBRSum());
	break;
      }
    }
  }

  return PList;

}

Particle * AdjustDecaySum(Particle * PParticle, DecayList * DList)
{
  int pid;

  pid = PParticle->GetID();

  for ( int j=1; j<=DList->GetSize(); j++ )
  {
    Decay * CurrentDecay = DList->Get(j);
    if ( pid==CurrentDecay->GetParentID() )
    {
      PParticle->SetDecaysum(CurrentDecay->GetBRSum());
      break;
    }
  }

  return PParticle;

}


