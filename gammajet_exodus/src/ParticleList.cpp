//-----------------------------------------------------------------------------
//
//  Implementation of the class ParticleList
//
//-----------------------------------------------------------------------------

#include <iostream>
#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"


PLInternalNode::PLInternalNode(Particle *theParticle, PLNode *next)
{
  ParticleInList = theParticle;
  Next           = next;
}

PLNode * PLInternalNode::Insert(Particle *theParticle)
{
   Next = Next->Insert(theParticle);
   return this;
}

void PLInternalNode::InsertAfter(Particle *theParticle)
{
   PLInternalNode * dataNode = new PLInternalNode(theParticle,Next);
   Next                      = dataNode;
}

void PLInternalNode::ShowOne(int n)
{
  if ( n==0 )
    ParticleInList->Show();
  else 
    Next->ShowOne(n-1);
}

Particle * PLInternalNode::Get(int n)
{
  if ( n==0 )
    return ParticleInList;
  else 
    return Next->Get(n-1);
}

void PLInternalNode::Put(Particle *newpart, int n)
{
  if ( n==0 )
    ParticleInList = newpart;
  else
   Next->Put(newpart,n-1);
}


PLNode * PLTailNode::Insert(Particle *theParticle)
{
  PLInternalNode * dataNode = new PLInternalNode(theParticle,this);
  return dataNode;
}

void PLTailNode::InsertAfter(Particle *theParticle)
{
  cout << "I don't know to insert something after the tailnode of the list!";
  cout << endl;
}

PLHeadNode::PLHeadNode()
{
  Next = new PLTailNode;
}

PLNode * PLHeadNode::Insert(Particle *theParticle)
{
  Next = Next->Insert(theParticle);
  return this;
}

void PLHeadNode::InsertAfter(Particle *theParticle)
{
  PLInternalNode * dataNode = new PLInternalNode(theParticle,Next);
  Next                      = dataNode;
}

ParticleList::ParticleList()
{
  Head   = new PLHeadNode; 
  Length = 0;
}

void ParticleList::Insert(Particle *pParticle)
{
  Head->Insert(pParticle); 
  Length++;
}

void ParticleList::InsertAfter(PLNode *theCurrent, Particle *theParticle)
{
  theCurrent->InsertAfter(theParticle); 
  Length++;
}

void ParticleList::ShowAll()
{
  cout << endl;
  cout << "____________________________________________________________________________" << endl;
  cout << "|___Type_______|__Energy__|___px_____|___py_____|___pz_____|Decaysum|Generation|" << endl;
  Head->Show();
  cout << "|______________|__________|__________|__________|__________|____|__________|" << endl;
  cout << "Total Particles:";
  cout.width(9); 
  cout << Length; 
  cout << endl;
}


