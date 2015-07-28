//-----------------------------------------------------------------------------
//
//  Implementation of the class DecayList
//
//-----------------------------------------------------------------------------

#include <iostream>
#include "Momentum.h"
#include "DecayList.h"

DLInternalNode::DLInternalNode(Decay *thePart, DLNode *theNext)
{
  myNext = theNext;
  myPart = thePart;
}

DLNode * DLInternalNode::Insert(Decay *thePart)
{
 myNext = myNext->Insert(thePart);
 return this;
}
  
Decay *DLInternalNode::Get(int n)
{
  if ( n==0 )
    return myPart;
  else
    return myNext->Get(n-1);
}

Decay *DLInternalNode::GetByID(int n)
{
  if ( n==myPart->GetID() )
    return myPart;
  else
    return myNext->GetByID(n);
}

DLNode * DLTailNode::Insert(Decay *part)
{
  DLInternalNode * dNode = new DLInternalNode(part,this);
  return dNode;
}

Decay *DLTailNode::Get(int n)
{
  Decay * zero = new Decay;
  zero->SetEnabled(0);
  return zero;
}
Decay *DLTailNode::GetByID(int n)
{
  Decay * zero = new Decay;
  zero->SetEnabled(0);
  return zero;
}

DLHeadNode::DLHeadNode()
{
  myNext = new DLTailNode; 
}

DLNode * DLHeadNode::Insert(Decay *thePart)
{
  myNext = myNext->Insert(thePart);
  return this;
}

DecayList::DecayList()
{
  myHead=new DLHeadNode;
}

DecayList::~DecayList()
{
  delete myHead;
  myHead = 0;
}

void DecayList::Insert(Decay *pPart)
{
  myHead->Insert(pPart);
  mySize++;
}
  










