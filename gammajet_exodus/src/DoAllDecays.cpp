//-----------------------------------------------------------------------------
//
//  Main steering routine of the decay machine
//
//-----------------------------------------------------------------------------

#include <math.h>

#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticlePropertyList.h"
#include "DecayList.h"

void TwoBodyDecay(Particle *, Particle *, Particle *, 
		  ParticlePropertyList *, Decay *);
void ThreeBodyDecay(Particle *, Particle *, Particle *, Particle *,
		    ParticlePropertyList *, Decay *);
void DalitzDecay(Particle *, Particle *, Particle *, Particle *,
		 ParticlePropertyList *, Decay *);
Particle     * AdjustDecaySum(Particle *, DecayList *);
ParticleList * AdjustDecaySum(ParticleList *, DecayList *);

ParticleList * DoAllDecays(ParticleList * PList, 
                           ParticlePropertyList * PPList, DecayList * DList)
{
  Particle     * PParent   = 0;
  Particle     * PChild1   = 0;
  Particle     * PChild2   = 0;
  Particle     * PChild3   = 0;
  Decay        * PDecay    = 0;

  int itotal, particle_counter, decay_counter, selected_decay, nb, nlep;
  float rndm, new_weight;
  float decay_matrix[5][2];

  PLNode * CurrentNode = PList->GetHeadNode();
  itotal           = 0;
  particle_counter = 0;
  for (;;)
  {
    particle_counter++;
    if ( particle_counter>PList->GetLength() ) break;

    CurrentNode = CurrentNode->GetNextNode();
    PParent = CurrentNode->Get(0);

    if ( PParent->GetGeneration()==1 )
    {
      itotal++;
      if ( fmod((double)itotal,10000.0)==0.0 )
	cout << itotal << " primary particles done" << endl;
    }

    if ( PParent->GetDecaysum()!=0.0 )
    {
      int possible_decays = 0;
      selected_decay = 0;
      for ( decay_counter=1; decay_counter<=DList->GetSize(); decay_counter++ )
      {
	PDecay = DList->Get(decay_counter);
	if ( (PParent->GetID())!=(PDecay->GetParentID()) )
	  continue;
       	decay_matrix[possible_decays][0] = (float)decay_counter;
	if ( possible_decays == 0 )
       	  decay_matrix[possible_decays][1] 
            = PDecay->GetBranchingRatio()/PDecay->GetBRSum();
        else
       	  decay_matrix[possible_decays][1] 
            = decay_matrix[possible_decays-1][1]
              + PDecay->GetBranchingRatio()/PDecay->GetBRSum();
	possible_decays++;    
      }

      if ( possible_decays )
      {
	rndm = gRandom->Rndm();
	for ( decay_counter=1; decay_counter<=possible_decays; 
              decay_counter++ )
	{       
	  if ( !selected_decay && rndm<=decay_matrix[decay_counter-1][1] )
	    selected_decay = (int)decay_matrix[decay_counter-1][0];
	}
	PDecay = DList->Get(selected_decay);

	if ( PDecay->GetNBody()==2 )
	{
	  PChild1 = new Particle;
	  PChild2 = new Particle;
	  TwoBodyDecay(PParent,PChild1,PChild2,PPList,PDecay); 
	  new_weight = PParent->GetWeight()*PParent->GetDecaysum();
	  PParent->SetDecaysum(0.0);
	  if ( !PDecay->GetChildrenStable() )
	  {
	    PChild1 = AdjustDecaySum(PChild1, DList);
	    PChild2 = AdjustDecaySum(PChild2, DList);
	  }
	  PChild1->SetWeight(new_weight);
	  PChild2->SetWeight(new_weight);
	  PList->InsertAfter(CurrentNode,PChild1);
	  PList->InsertAfter(CurrentNode,PChild2);
	}
	else if ( PDecay->GetNBody()==3 )
	{
	  PChild1 = new Particle;
	  PChild2 = new Particle;
	  PChild3 = new Particle;
	  nlep=0;
	  for( nb=1;nb<=3;nb++ )
	  {
	    if ( fabs(PDecay->GetChildID(nb))==11 ) 
	      nlep++;
	  }
	  if ( nlep==2 )
	    DalitzDecay(PParent,PChild1,PChild2,PChild3,
				  PPList,PDecay);
	  else
	    ThreeBodyDecay(PParent,PChild1,PChild2,PChild3,
				     PPList,PDecay);
	  new_weight = PParent->GetWeight()*PParent->GetDecaysum();
	  PParent->SetDecaysum(0.0);
	  if ( !PDecay->GetChildrenStable() )
	  {
	    PChild1 = AdjustDecaySum(PChild1, DList);
	    PChild2 = AdjustDecaySum(PChild2, DList);
	    PChild3 = AdjustDecaySum(PChild3, DList);
	  }
	  PChild1->SetWeight(new_weight);
	  PChild2->SetWeight(new_weight);
	  PChild3->SetWeight(new_weight);
	  PList->InsertAfter(CurrentNode,PChild1);
	  PList->InsertAfter(CurrentNode,PChild2);
	  PList->InsertAfter(CurrentNode,PChild3);
       	}
	else
	  cout << "Decays with more than 3 decay products are currently "
               << "not implemented" << endl;
      }
    }
  }

  return PList;

}









