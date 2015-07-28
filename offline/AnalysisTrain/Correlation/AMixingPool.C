#include "AMixingPool.h"
#include "AParticle.h"
#include "ACluster.h"
#include "ATrack.h"
#include "AEvent.h"

#include <iostream>

AMixingPool::AMixingPool(int nvtx, int ncent)
{
  _centbinwidth = 0.;
  _vtxbinwidth = 0.;
  _poolsize = 0;
  NZVTXBINS = nvtx;
  NCENTBINS = ncent;
  
  //initialize vector of vectors of deques
  //[NZVTXBINS][NCENTBINS]
  for( int i=0; i<NZVTXBINS; i++ ){
    _events.push_back(std::vector<std::deque<AEvent*> >(NCENTBINS) );
  }
  
}

AMixingPool::~AMixingPool()
{
  for(int i = 0; i < NZVTXBINS; i++){
    for(int j = 0; j < NCENTBINS; j++){
      ClearParticleList<AEvent>(_events[i][j]);
    }
  }
}

void AMixingPool::AddEvent(AEvent* event)
{
  int ivtx = event->GetVertexBin();
  int icent = event->GetCentralityBin();
  _events[ivtx][icent].push_back(event->clone());
  if(_events[ivtx][icent].size() > _poolsize){
    delete _events[ivtx][icent].front();
    _events[ivtx][icent].pop_front();
  }
}

bool AMixingPool::IsFullPool(int ivtx, int icent)
{
  bool full = true;
  if(_events[ivtx][icent].size() < _poolsize) full = false;
  
  return full;
}

