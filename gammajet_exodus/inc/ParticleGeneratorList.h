//-----------------------------------------------------------------------------
//
//  Declaration of the class ParticleGeneratorList
//
//-----------------------------------------------------------------------------

#ifndef PARTICLEGENERATORLIST_H
#define PARTICLEGENERATORLIST_H

#include <iostream>

class PGLNode;
class PGLHeadNode;
class PGLTailNode;
class PGLInternalNode;
 
class PGLNode
{
 public:
  PGLNode(){}
  virtual ~PGLNode(){}
  virtual PGLNode * Insert(ParticleGenerator *pGen)=0;
  virtual ParticleGenerator * Get(int)=0;
  virtual ParticleGenerator * GetByID(int)=0;
 private:
};

class PGLInternalNode: public PGLNode
{
 public:
  PGLInternalNode(ParticleGenerator *theGen, PGLNode *theNext);
  ~PGLInternalNode(){delete Next;Next=0;}
  virtual PGLNode * Insert(ParticleGenerator *pGen);
  ParticleGenerator * Get(int n);
  ParticleGenerator * GetByID(int n);
 private:
  PGLNode           * Next;
  ParticleGenerator * Generator;
};

class PGLTailNode: public PGLNode
{
 public:
  PGLTailNode(){}
  ~PGLTailNode(){}
 virtual PGLNode * Insert(ParticleGenerator *gen);
 virtual ParticleGenerator * Get(int n);
 virtual ParticleGenerator * GetByID(int n);
 private:
};

class PGLHeadNode: public PGLNode
{
 public:
  PGLHeadNode();
  ~PGLHeadNode(){delete Next; Next=0;}
  virtual PGLNode * Insert(ParticleGenerator *pGen);
  ParticleGenerator * Get(int n){return Next->Get(n-1);}
  ParticleGenerator * GetByID(int n){return Next->GetByID(n);}
 private:
  PGLNode * Next;
};

class ParticleGeneratorList
{
 public:
  ParticleGeneratorList();
  ~ParticleGeneratorList();
  void Insert(ParticleGenerator *pGen);
  ParticleGenerator * Get(int n){return Head->Get(n);}
  ParticleGenerator * GetByID(int n){return Head->GetByID(n);}
 private:
  PGLHeadNode * Head;
};

#endif /* PARTICLEGENERATORLIST_H */



