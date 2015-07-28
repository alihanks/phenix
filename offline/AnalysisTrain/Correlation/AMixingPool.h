#ifndef __AMIXINGPOOL_H__
#define __AMIXINGPOOL_H__

#include <deque>
#include <vector>

class AParticle;
class ACluster;
class ATrack;
class AEvent;

class AMixingPool {
public:
  AMixingPool(int nzvtx, int ncent);
  virtual ~AMixingPool();

  void AddEvent(AEvent* event);
  bool IsFullPool(int ivtx, int icent);

  void SetPoolSize(const unsigned int poolsize) { _poolsize = poolsize; }
  unsigned int GetNEvents(int ivtx, int icent) const { return _events[ivtx][icent].size(); }

  AEvent* GetEvent(int ivtx, int icent, unsigned int ievent){ return _events[ivtx][icent][ievent]; }

private:
  double _centbinwidth;
  double _vtxbinwidth;
  unsigned int _poolsize;
  int NZVTXBINS;
  int NCENTBINS;

  std::vector<std::vector<std::deque<AEvent*> > > _events;

  template<class T> void ClearParticleList(std::deque<T*>& particles){
    unsigned int npart = particles.size();
    for(unsigned int ipart=0; ipart<npart; ipart++){
      delete particles[ipart];
    }
    particles.clear();
  }

};

#endif /* __AMIXINGPOOL_H__ */
