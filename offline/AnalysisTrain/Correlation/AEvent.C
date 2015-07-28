#include "AEvent.h"
#include "ACluster.h"
#include "ATrack.h"
#include "PHCompositeNode.h"

using namespace std;

ClassImp(AEvent);

AEvent::AEvent() : _clusters(), _tracks(), _lessqualitytracks()
{
  _zvtx = -9999;
  _cent = -9999;
  _vtxbin = -1;
  _centbin = -1;
}

AEvent::AEvent(const AEvent& evt)
{
  _zvtx = evt._zvtx;
  _cent = evt._cent;
  _vtxbin = evt._vtxbin;
  _centbin = evt._centbin;

  SetTracks(evt._tracks);
  SetClusters(evt._clusters);
  SetLessQualityTracks(evt._lessqualitytracks);
}

AEvent::~AEvent(){

  for(unsigned int i = 0; i < _clusters.size(); i++){
    delete _clusters[i];
  }
  _clusters.clear();
  for(unsigned int i = 0; i < _tracks.size(); i++){
    delete _tracks[i];
  }
  _tracks.clear();
  for(unsigned int i = 0; i < _lessqualitytracks.size(); i++){
    delete _lessqualitytracks[i];
  }
  _lessqualitytracks.clear();
}

void AEvent::SetCentrality(int ncent, double bin, float cent)
{
  if(cent < 0 || cent > 100) _cent = -1;
  else _cent = cent;
  float CentBinWidth = bin / ncent;
  int centbin = int( ( cent - 1.0 ) / CentBinWidth ); // -1.0 in order to be consisitent with the centrality binning as combinesimple -042915
  if(centbin < 0 || centbin > (ncent-1)) _centbin = -1;
  else _centbin = centbin;
  //cout<<"AEvent2: cent = "<<cent<<"; centbin = "<<centbin<<"; _centbin = "<<_centbin<<endl;
}

void AEvent::SetVertex(int nvtx, double bin, float vtx)
{
  _zvtx = vtx;
  float VtxBinWidth = bin / nvtx;
  int vtxbin = int( (vtx + bin/2) / VtxBinWidth );//potentially buggy
  if(vtx < -30.0 || vtx > 30.0) _vtxbin = -1;
  else _vtxbin = vtxbin;
  // cout<<"AEvent: vtx = "<<vtx<<"; vtxbin = "<<vtxbin<<endl;
}

void AEvent::SetClusters(const vector<ACluster*>& clusters) {
  for(unsigned int i = 0; i < clusters.size(); i++)
    _clusters.push_back(clusters[i]->clone());
}

void AEvent::SetTracks(const vector<ATrack*>& tracks) {
  for(unsigned int i = 0; i < tracks.size(); i++)
    _tracks.push_back(tracks[i]->clone());
}

void AEvent::SetLessQualityTracks(const vector<ATrack*>& lessqualtracks) {
  _lessqualitytracks.clear();
  for(unsigned int i = 0; i < lessqualtracks.size(); i++)
    _lessqualitytracks.push_back(lessqualtracks[i]->clone());
}


AEvent* AEvent::clone() const
{
  return new AEvent(*this);
}

void AEvent::ClearEvent()
{
  for (unsigned int i = 0; i < _clusters.size(); i++){
    delete _clusters[i];
  }
  _clusters.clear();
  for (unsigned int i = 0; i < _tracks.size(); i++){
    delete _tracks[i];
  }
  _tracks.clear();
  for (unsigned int i = 0; i < _lessqualitytracks.size(); i++){
    delete _lessqualitytracks[i];
  }
  _lessqualitytracks.clear();
}
