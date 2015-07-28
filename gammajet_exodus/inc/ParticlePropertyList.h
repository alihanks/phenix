//-----------------------------------------------------------------------------
//
//  Declaration of the class ParticlePropertyList
//
//-----------------------------------------------------------------------------

#ifndef PARTICLEPROPERTYLIST_H
#define PARTICLEPROPERTLIST_H

#include <iostream>

class PPLNode;
class PPLHeadNode;
class PPLTailNode;
class PPLInternalNode;
 
class PPLNode
{
 public:
  PPLNode(){}
  virtual ~PPLNode(){}
  virtual PPLNode * Insert(ParticleProperty *pDec)=0;
  virtual ParticleProperty * Get(int)=0;
  virtual ParticleProperty * GetByID(int)=0;
 private:
};

class PPLInternalNode: public PPLNode
{
 public:
  PPLInternalNode(ParticleProperty *thePart, PPLNode *theNext);
  ~PPLInternalNode(){delete Next;Next=0;}
  virtual PPLNode * Insert(ParticleProperty *pDec);
  ParticleProperty * Get(int n);
  ParticleProperty * GetByID(int n);
 private:
  PPLNode          * Next;
  ParticleProperty * Particle;
};

class PPLTailNode: public PPLNode
{
 public:
  PPLTailNode(){}
  ~PPLTailNode(){}
 virtual PPLNode * Insert(ParticleProperty *part);
 virtual ParticleProperty * Get(int n);
 virtual ParticleProperty * GetByID(int n);
 private:
};

class PPLHeadNode: public PPLNode
{
 public:
  PPLHeadNode();
  ~PPLHeadNode(){delete Next; Next=0;}
  virtual PPLNode * Insert(ParticleProperty *pDec);
  ParticleProperty * Get(int n){return Next->Get(n-1);}
  ParticleProperty * GetByID(int n){return Next->GetByID(n);}
 private:
  PPLNode * Next;
};

class ParticlePropertyList
{
 public:
  ParticlePropertyList();
  ~ParticlePropertyList();
  void Insert(ParticleProperty *pDec);
  ParticleProperty * Get(int n){return Head->Get(n);}
  ParticleProperty * GetByID(int n){return Head->GetByID(n);}
 private:
  PPLHeadNode * Head;
};

#endif /* PARTICLEPROPERTYLIST_H */

