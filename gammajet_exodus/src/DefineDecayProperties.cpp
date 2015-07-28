//-----------------------------------------------------------------------------
//
//  Read decay properties from file 
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <stdlib.h>
#include <string>

#include <TH1.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticlePropertyList.h"
#include "DecayList.h"

TH1F * CreateDalitzLeptonPairMass(Decay*, ParticlePropertyList*);

DecayList * DefineDecayProperties(ParticlePropertyList *PPList)
{
  double branchingratio;
  int id, i;
  int nbody, parentid, childid[4];
  char ch;
  bool enabled, children_stable, dalitz;
  string label, name;

  DecayList * DefinedDecays = new DecayList;

  cout << "Reading particle definitions from file: defined_decays.txt"
       << endl;

  ifstream fin("defined_decays.txt");
  while ( fin.get(ch) ) 
    {
      i=0;
      fin >> label;
      if ( label=="@{" )
      {
        fin >> name;
        for(;;)
	{
	  fin >> label;
	  if ( label=="DecayId=" )        fin >> id; 
	  if ( label=="NBody=" )          fin >> nbody;
	  if ( label=="BranchingRatio=" ) fin >> branchingratio;
	  if ( label=="ParentID=" )       fin >> parentid;
	  if ( label=="ChildID=" )
	  {
	    i++;
	    fin >> childid[i];
	  }
	  if ( label=="Enabled=" )        fin >> enabled;
	  if ( label=="ChildrenStable=" ) fin >> children_stable;
	  if ( label=="}" )               break;
	}
      }
      else continue;
      if ( nbody!=2 && nbody!= 3 )
      {
	cout << "Only 2 and 3 body decays are implemented at the moment."
             << endl
             << "Check the file 'defined_decays.txt' for invalid decays."
	     << endl;
	return 0;
      }
      if ( !enabled ) continue;
      
      Decay *pDec = new Decay;

      pDec->SetID(id);
      pDec->SetNBody(nbody);
      pDec->SetBranchingRatio(branchingratio);
      pDec->SetParentID(parentid);
      for ( i=1; i<=nbody; i++)
        pDec->SetChildID(i,childid[i]);
      pDec->SetChildrenStable(children_stable);
      DefinedDecays->Insert(pDec);

      dalitz = false;
      if ( nbody==3 )
      {
	bool lepton_n = false;
        bool lepton_p = false;
	for ( i=1; i<=3; i++ )
	{
	  if ( childid[i]==-11 || childid[i]==-13 ) lepton_p = true;
	  if ( childid[i]==11  || childid[i]==13 )  lepton_n = true;
	}
	if ( lepton_p && lepton_n ) dalitz = true;
      }

      if ( dalitz )
      {
        TH1F * hdalitz = 0;
	hdalitz = CreateDalitzLeptonPairMass(pDec, PPList);
	pDec->SetHistogram(hdalitz);
      }
      else
      {
	pDec->SetHistogram(0);
      }

      cout << "Defined: " << name << endl;
    }
  fin.close();
    
  cout << endl;

  double pid;
  double br;
  int MAX_PARTICLE_TYPES = 20;
  double brsum_array[MAX_PARTICLE_TYPES][2];
  for ( i=1; i<=MAX_PARTICLE_TYPES; i++ )
  {
    brsum_array[i-1][0] = 0.0;
    brsum_array[i-1][1] = 0.0;
  } 

  for ( i=1; i<=DefinedDecays->GetSize(); i++ )
  {
    Decay * CurrentDecay = DefinedDecays->Get(i);
    pid = (double)CurrentDecay->GetParentID();
    br = CurrentDecay->GetBranchingRatio();
    for ( int j=1; j<=MAX_PARTICLE_TYPES; j++ )
    {
      if ( brsum_array[j-1][0]==pid )
      {
	brsum_array[j-1][1]+=br;
	break;
      }
      if ( brsum_array[j-1][0]==0.0 ) 
      {
	brsum_array[j-1][0] = pid;
	brsum_array[j-1][1] = br;
	break;
      }
    }
  }

  for ( i=1; i<=DefinedDecays->GetSize(); i++ )
  {
    Decay * CurrentDecay = DefinedDecays->Get(i);
    pid = (double)CurrentDecay->GetParentID();
    for ( int j=1; j<=MAX_PARTICLE_TYPES; j++ )
    {
      if ( pid==brsum_array[j-1][0] )
      {
	CurrentDecay->SetBRSum(brsum_array[j-1][1]);
	break;
      }
    }
  }    

  return DefinedDecays;

}


