//-----------------------------------------------------------------------------
//
//  Declaration of the classes Decay and DecayList
//
//-----------------------------------------------------------------------------

#ifndef DECAYLIST_H
#define DECAYLIST_H

#include <iostream>
#include <TH1.h>

class Decay
{
 public:
  Decay(){}
  ~Decay(){}
  void SetID(int id){itsID=id;}
  void SetNBody(int nb){itsNBody=nb;}
  void SetBranchingRatio(double br){itsBranchingRatio=br;}
  void SetBRSum(double brsum){itsBRSum=brsum;}
  void SetParentID(int pID){itsParentID=pID;}
  void SetChildID(int i, int cID){itsChildID[i]=cID;}
  void SetEnabled(bool en){itsEnabled=en;}
  void SetChildrenStable(bool stable){itsChildrenStable=stable;}
  void SetHistogram(TH1F *histogram){itsHistogram=histogram;}
  int GetID(){return itsID;}
  int GetNBody() const {return itsNBody;}
  double GetBranchingRatio(){return itsBranchingRatio;}
  double GetBRSum(){return itsBRSum;}
  int GetParentID(){return itsParentID;}
  int GetChildID(int i){return itsChildID[i];}
  bool GetEnabled(){return itsEnabled;}
  bool GetChildrenStable(){return itsChildrenStable;}
  TH1F* GetHistogram(){return itsHistogram;}
 private:
  int    itsID;
  int    itsNBody;
  double itsBranchingRatio;
  double itsBRSum;
  int    itsParentID;
  int    itsChildID[10];
  bool   itsEnabled;
  bool   itsChildrenStable;
  TH1F*  itsHistogram;
};  

class DLNode;
class DLHeadNode;
class DLTailNode;
class DLInternalNode;
 
class DLNode
{
 public:
  DLNode(){}
  virtual ~DLNode(){}
  virtual DLNode * Insert(Decay *pDec)=0;
  virtual Decay * Get(int)=0;
  virtual Decay * GetByID(int)=0;
 private:
};

class DLInternalNode: public DLNode
{
 public:
  DLInternalNode(Decay *thePart, DLNode *theNext);
  ~DLInternalNode(){delete myNext; myNext=0;}
  virtual DLNode * Insert(Decay * pDec);
  Decay * Get(int n);
  Decay * GetByID(int n);
 private:
  DLNode * myNext;
  Decay  * myPart;
};

class DLTailNode: public DLNode
{
 public:
  DLTailNode(){}
  ~DLTailNode(){}
 virtual DLNode * Insert(Decay *part);
 virtual Decay * Get(int n);
 virtual Decay * GetByID(int n);
 private:

};

class DLHeadNode: public DLNode
{
 public:
  DLHeadNode();
  ~DLHeadNode(){delete myNext; myNext=0;}
  virtual DLNode * Insert(Decay *pDec);
  Decay * Get(int n){return myNext->Get(n-1);}
  Decay * GetByID(int n){return myNext->GetByID(n);}
 private:
  DLNode * myNext;
};

class DecayList
{
 public:
  DecayList();
  ~DecayList();
  void Insert(Decay *pDec);
  Decay * Get(int n){return myHead->Get(n);}
  Decay * GetByID(int n){return myHead->GetByID(n);}
  DLNode * GetHeadNode(){return myHead;}
  int GetSize() const {return mySize;} 
 private:
  DLHeadNode * myHead;
  int          mySize;
};

#endif /* DECAYLIST_H */
