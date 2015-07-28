//-----------------------------------------------------------------------------
//
//  Declaration of the class ParticleList
//
//-----------------------------------------------------------------------------

#ifndef PARTICLELIST_H
#define PARTICLELIST_H 

#include <iostream>

class PLNode;
class PLHeadNode;
class PLTailNode;
class PLInternalNode;

class PLNode
{
public:
  PLNode(){}
  virtual ~PLNode(){}
  virtual PLNode * Insert(Particle *theParticle)=0;
  virtual void InsertAfter(Particle *theParticle)=0;
  virtual void Show()=0;
  virtual void ShowOne(int)=0;
  virtual Particle * Get(int)=0;
  virtual void Put(Particle*, int)=0;
  virtual PLNode * GetNextNode()=0;
private:
};

class PLInternalNode : public PLNode
{
public:
  PLInternalNode(Particle *theParticle, PLNode *next);
  ~PLInternalNode(){delete Next; Next=0; 
                    delete ParticleInList; ParticleInList=0;}
  virtual PLNode * Insert(Particle *theParticle);
  virtual void InsertAfter(Particle *theParticle);
  virtual void Show(){ParticleInList->Show(); Next->Show();}
  virtual void ShowOne(int n);
  PLNode * GetNextNode(){return Next;}
  Particle * Get(int n); 
  void Put(Particle *newpart, int n);
private:
  Particle * ParticleInList;
  PLNode   * Next;
};

class PLTailNode : public PLNode
{
public:
  PLTailNode(){}
  ~PLTailNode(){}
  virtual PLNode * Insert(Particle *theParticle);
  virtual void InsertAfter(Particle *theParticle);
  virtual void Show(){};
  virtual void ShowOne(int n){cout << "End of particle list reached" << endl;}
  PLNode * GetNextNode(){return this;}
  virtual Particle * Get(int n){cout << "End of particle list reached" << endl;
                                Particle *p = new Particle; return p;}
  virtual void Put(Particle *newpart, int n) {}
private:
};

class PLHeadNode : public PLNode
{
public:
  PLHeadNode();
  ~PLHeadNode(){delete Next;}
  virtual PLNode * Insert(Particle *theParticle);
  virtual void InsertAfter(Particle *theParticle);
  virtual void Show(){Next->Show();}
  virtual void ShowOne(int n){Next->ShowOne(n-1);}
  PLNode * GetNextNode(){return Next;}
  virtual Particle * Get(int n){return Next->Get(n-1);}
  virtual void Put(Particle *newpart, int n){Next->Put(newpart,n-1);}
private:
  PLNode * Next;
};

class ParticleList
{
public:
  ParticleList();
  ~ParticleList(){delete Head; Head=0;}
  void Insert(Particle *theParticle);
  void InsertAfter(PLNode *theCurrent, Particle *theParticle);
  void ShowAll();
  void ShowOne(int n){Head->ShowOne(n);}
  Particle * Get(int n){return Head->Get(n);}
  PLNode * GetHeadNode(){return Head;}
  void Put(Particle *newpart, int n){Head->Put(newpart,n);}
  int GetLength(){return Length;}
private:
  PLHeadNode * Head;
  int          Length;
};

#endif /* PARTICLELIST_H */





