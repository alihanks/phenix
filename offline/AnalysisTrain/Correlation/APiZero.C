#include "APiZero.h"
#include "ACluster.h"

ClassImp(APiZero);

APiZero::APiZero() : AParticle(111)
{
  _daughter1 = new ACluster();
  _daughter2 = new ACluster();
  _mwweights.clear();
}

APiZero::APiZero(AParticle* g1, AParticle* g2) : AParticle(g1,g2,111)
{
  _mwweights.clear();
}

APiZero::APiZero(const APiZero& pizero) : AParticle(pizero.Daughter1(),pizero.Daughter2(),111)
{
  SetDecayWeights(pizero.GetDecayWeights());
  _is_iso = pizero.IsIso();
  _pid = pizero.GetPid();
  _event = pizero.GetEvent();
  _zvtx = pizero.GetZvtx();
  _cent = pizero.GetCent();
  _arm = pizero.GetArm();
  _sec = pizero.GetSec();
  _sm = pizero.GetSm();
}

APiZero* APiZero::clone() const
{
  return new APiZero(*this);
}

// APiZero* APiZero::clone_vector() const
// {
//   APiZero* api0 = new APiZero();
//   api0->SetPtEtaPhiE(Perp(),PseudoRapidity(),fP.Phi(),T());
//   api0->SetDecayWeights(_mwweights);
//   return api0;
// }

void APiZero::SetDecayWeights(std::vector<float> weights)
{ 
  for( unsigned int i = 0; i < weights.size(); i++ )
  {
    _mwweights.push_back(weights[i]);
  }
}
