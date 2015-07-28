#ifndef __AEVENT_H__
#define __AEVENT_H__

#include <TObject.h>
#include <vector>

class ACluster;
class ATrack;
class PHCompositeNode;

class AEvent : public TObject {

public:
  AEvent();
  AEvent(const AEvent& evt);
  virtual ~AEvent();

  void SetVertex(int nvtx, double bin, float vtx);
  void SetCentrality(int ncent, double bin, float cent);
  void SetClusters(const std::vector<ACluster*>& clusters);
  void SetTracks(const std::vector<ATrack*>& tracks);
  void SetLessQualityTracks(const std::vector<ATrack*>& lessqualtracks);//for tracks with quality>7
  void AddTrack(ATrack* track) { _tracks.push_back(track); }

  float GetVertex() const { return _zvtx; }
  float GetCentrality() const { return _cent; }
  int GetCentralityBin() const { return _centbin; }
  int GetVertexBin() const { return _vtxbin; }
  unsigned int GetNClusters() { return _clusters.size(); }
  unsigned int GetNTracks() { return _tracks.size(); }
  unsigned int GetNLessQualityTracks() { return _lessqualitytracks.size(); }
  ACluster* GetCluster(int index){ return _clusters[index]; }
  ATrack* GetTrack(int index){ return _tracks[index]; }
  ATrack* GetLessQualityTrack(int index){ return _lessqualitytracks[index]; }//for tracks with quality>7
  std::vector<ACluster*> GetClusters() { return _clusters; }
  std::vector<ATrack*> GetTracks() { return _tracks; }
  std::vector<ATrack*> GetLessQualityTracks() { return _lessqualitytracks; }
  void ClearEvent();

  AEvent* clone() const;


private:
  float _zvtx;
  float _cent;
  int _vtxbin;
  int _centbin;

  std::vector<ACluster*> _clusters;
  std::vector<ATrack*> _tracks;
  std::vector<ATrack*> _lessqualitytracks;

  ClassDef(AEvent,1);
};

#endif /* __AEVENT_H__ */
