//-----------------------------------------------------------------------------
//
//  Implementation of the class ParticleGeneratorList
//
//-----------------------------------------------------------------------------

#include <iostream>
#include "Momentum.h"
#include "Particle.h"
#include "ParticleGeneratorList.h"

PGLInternalNode::PGLInternalNode(ParticleGenerator *theGen, PGLNode *theNext)
{
  Next      = theNext;
  Generator = theGen;
}

PGLNode * PGLInternalNode::Insert(ParticleGenerator *theGen)
{
 Next = Next->Insert(theGen);
 return this;
}
  
ParticleGenerator * PGLInternalNode::Get(int n)
{
  if ( n==0 )
    return Generator;
  else
    return Next->Get(n-1);
}

ParticleGenerator * PGLInternalNode::GetByID(int n)
{
  if (n == Generator->GetID())
    return Generator;
  else
    return Next->GetByID(n);
}

PGLNode * PGLTailNode::Insert(ParticleGenerator *gen)
{
  PGLInternalNode * gNode = new PGLInternalNode(gen,this);
  return gNode;
}

ParticleGenerator * PGLTailNode::Get(int n)
{
  ParticleGenerator * zero = new ParticleGenerator;
  TH1F * Pzero = 0;
  zero->SetID(0);
  zero->SetWeight(0.0);
  zero->SetPtHistogram(Pzero);
  zero->SetYHistogram(Pzero);
  zero->SetMHistogram(Pzero);
  return zero;
}

ParticleGenerator * PGLTailNode::GetByID(int n)
{
  ParticleGenerator * zero = new ParticleGenerator;
  TH1F * Pzero = 0;
  zero->SetID(0);
  zero->SetWeight(0.0);
  zero->SetPtHistogram(Pzero);
  zero->SetYHistogram(Pzero);
  zero->SetMHistogram(Pzero);
  return zero;
}

PGLHeadNode::PGLHeadNode()
{
  Next = new PGLTailNode; 
}

PGLNode * PGLHeadNode::Insert(ParticleGenerator *theGen)
{
  Next = Next->Insert(theGen);
  return this;
}

ParticleGeneratorList::ParticleGeneratorList()
{
  Head = new PGLHeadNode;
}

ParticleGeneratorList::~ParticleGeneratorList()
{
  delete Head;
  Head = 0;
}

void ParticleGeneratorList::Insert(ParticleGenerator *pGen)
{
  Head->Insert(pGen);
}
  
















