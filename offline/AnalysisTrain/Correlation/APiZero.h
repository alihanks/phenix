#ifndef __APIZERO_H__
#define __APIZERO_H__

#include <TLorentzVector.h>
#include <vector>
#include <AParticle.h>

class APiZero : public AParticle {

public:
  APiZero();
  APiZero(AParticle* g1, AParticle* g2);
  APiZero(const APiZero&);
  virtual ~APiZero() { if(_daughter1) delete _daughter1; if(_daughter2) delete _daughter2; };

  std::vector<float> GetDecayWeights() const { return _mwweights; }
  void SetDecayWeights(std::vector<float> weights);

  bool IsTagged() const { return _daughter1->IsTagged(); }
  bool IsFiducial() const { return _daughter1->IsFiducial(); }
 
  APiZero* clone() const;
  //APiZero* clone_vector() const;

private:
  std::vector<float> _mwweights;

  ClassDef(APiZero,1);
};

#endif /* __APIZERO_H__ */

