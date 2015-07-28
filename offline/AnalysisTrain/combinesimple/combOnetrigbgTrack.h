#ifndef __COMBONETRIGBGTRACK_H__
#define __COMBONETRIGBGTRACK_H__
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <TString.h>
#include <CombinedSimple.h>

#include <stdio.h>
#include <stdlib.h>

class PHCompositeNode;
class PHGlobal;
class PHCentralTrack;
class TTree;
class TH1;
class TH2;
class TH3;
class TF1;
class TFile;
class TFormula;
class TNtuple;

using namespace std;


class combOnetrigbgTrack: public CombinedSimple
{
 public:
  combOnetrigbgTrack(int, int, int, int, int, int, int);
  virtual ~combOnetrigbgTrack() {}
  int Init(PHCompositeNode *);
  virtual int process_event(PHCompositeNode *topNode);
  virtual int End(PHCompositeNode *topNode); 
  
  int process_event( int);
  void Print(const char *what) const {return;}
  int End(char *);
  void ResetMultCounter(){used_events_trig=0;}

  void SetTrigProps(float, float, float, float, float, float, float, float, float, float, float, float, float, float, float);
  int IsFinished();
  void ResetFinished();

  void SetMixingChainFile(char *);
  void DeleteMixingChainFile();
  void Mix(TNtuple * trig, TTree * partners = 0);
  void set_segmentFile(int segfile) { _segFile = segfile; } 
  void set_mixFilepath(char * newpath) { _mixFilepath = newpath; }
  void set_mixTrigNode(CombinedSimple * cst) { _mixTrig = cst; }
  void set_mixSingNode(CombinedSimple * cst) { _mixSing = cst; }
  
 protected:
  
  CombinedSimple * _mixTrig, *_mixSing; // when registering these hold mix ntuples

  int _segFile;
  TString _mixFilepath;
  TTree * partners;

  int event_counter;
  int minbias_counter;

  int tagflag, minbias_ert, iscompressed, species;
  
  vector<int> trigger_list;
 
  //for now just hard code this:
  int total_minbias;
  
  int isfinished;

  float lorun, hirun, multiplier;
  int total_events, used_events_run, used_events_trig;

  TriggerHolder trigger;
  float trig_photon_vert, trig_photon_runno, trig_photon_cent;

};

#endif /* __COMBONETRIGBGTRACK_H__ */
