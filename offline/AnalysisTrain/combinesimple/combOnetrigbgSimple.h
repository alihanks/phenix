#ifndef __combONETRIGBGSIMPLE_H__
#define __combONETRIGBGSIMPLE_H__
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <CombinedSimple.h>

#include <stdio.h>
#include <stdlib.h>

class PHCompositeNode     ;
class PHGlobal      ;
class PHCentralTrack      ;
class emcClusterContainer      ;
class ReactionPlaneObject;
class TNtuple;
class TH1;
class TH2;
class TH3;
class TF1;
class TFile;
class TFormula;

using namespace std;

class combOnetrigbgSimple: public CombinedSimple
{
 public:

  combOnetrigbgSimple(int, int, int, int, float, float, int);
  virtual ~combOnetrigbgSimple() {}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode) {return 0;}
  void SetInputName(char *triggerfile) { filename=triggerfile;}
  int process_event(PHCompositeNode *topNode);
  int Reset(PHCompositeNode *topNode) {return 0;}
  void Print(const char *what) const {return;}
  int ResetEvent(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void ResetMultCounter(){used_events_trig=0;}  

  void SetTrigProps(float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float);
  int IsFinished();
  void ResetFinished();


  void SetupProximityMixing();

  TString * _cntnodename;
  TString * _emcnodename;
  void  set_cntnodename(char * fn) 
    { if (_cntnodename) delete _cntnodename; _cntnodename = new TString(fn);}
  void  set_emcnodename(char * fn)
    { if (_emcnodename) delete _emcnodename; _emcnodename = new TString(fn);}


 protected:

  PHGlobal      *global;     
  PHCentralTrack      *particle;     
  ReactionPlaneObject *reacplane;
  emcClusterContainer *emccluster;

  int event_counter;
  int minbias_counter;

  int tagflag, minbias_ert, species, mixrxn;
  
  float centlo, centhi, locent, hicent;

  vector<int> trigger_list;
 
  //for now just hard code this:
  int total_minbias;
  
  char *filename;
  
  int isfinished;

  float lorun, hirun, multiplier;
  int total_events, used_events_run, used_events_trig;

  TriggerHolder trigger;
  float trig_photon_vert, trig_photon_cent, trig_photon_rp, trig_photon_runno, trig_photon_seqno;

  TNtuple * trigger_ntuple;


  TH2F *h3pc[25][25];  
  TH1F *hpctest[25][25];
  TH1F *htrtest[25];

  float * buf_dphis[17][1000];
  int buf_dphiscounter[1000];
  int last_buf;
  int cur_max_buf;
  
  int seqno_bypass;
  
    
};

#endif /* __combONETRIGBGSIMPLE_H__ */
