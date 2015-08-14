#ifndef hijing_analysis_h
#define hijing_analysis_h

#include <SubsysReco.h>
#include <string>
#include <vector>
#include <CorrelationFunctions.h>

class PHCompositeNode;
class Fun4AllHistoManager;
class TH1F;
class TH2F;
class TH3F;
class TTree;
class THmulf;
class ACluster;
class ATrack;
class APiZero;
class AMixingTree;
class TLorentzVector;

namespace HepMC
{
  class GenEvent;
  class GenVertex;
  class GenParticle;
};

class hijing_analysis: public SubsysReco
{

public:
  hijing_analysis(const char* output="test.root", const char* name = "NONAMECUTTER");
  virtual ~hijing_analysis();
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  void Verbosity(int v) { verbosity = v; }
  int End(PHCompositeNode *topNode);
  void SetSharkFin(const char* filename);
  const char* outfile;
  
  Fun4AllHistoManager *manager;
  float Rcut;
  int Centrality;
  double _MaxEta;
  double _MinAssocPt;
  double _MinTrigPt;

protected:

  bool MakeCluster(HepMC::GenParticle* p, ACluster* clus);
  bool MakePiZero(HepMC::GenParticle* p, APiZero* clus);
  bool MakeTrack(HepMC::GenParticle* p, ATrack* clus);
  int GetCentrality(int mult);
  void DoMixing(TTree* trig, TTree* assoc, int size);
  void MakeClusterObject(ACluster* aclus, float pt, float phi, float eta, float e);
  void MakePi0Object(APiZero* api0, float pt, float phi, float eta, float e);
  void MakeTrackObject(ATrack* atrk, float pt, float phi, float eta, float e);
  int CheckPool(int nenpart, int j, int pooldepth, int size, int& nloop);
  void EvalDecWeights(APiZero* pi0trigger, std::vector<float>& mwweight);
  void ApplyEnergyResolution(TLorentzVector* mom4, int pbsc_pbgl);

  void Init1DHisto(TH1F*& h1, std::string name, std::string xtitle, int nxbin, double xmin, double xmax);
  void Init2DHisto(TH2F*& h2, std::string name, std::string xtitle, int nxbin, double xmin, double xmax, std::string ytitle, int nybin, double ymin, double ymax);
  void Init3DHisto(TH3F*& h3, std::string name, std::string xtitle, int nxbin, double xmin, double xmax, std::string ytitle, int nybin, double ymin, double ymax, std::string ztitle, int nzbin, double zmin, double zmax);
  bool OutsideAcceptance(double phi);
  
  std::string ThisName;
  int verbosity;
  int nevents;
  AMixingTree* atree;
  int NMIX;

  TH1D* hshark[5];
  
  TH1F* hmult;
  TH1F* h1_mass;
  TH1F* h1_trigger_pt;
  TH1F* h1_trigger_pi0_pt;
  TH1F* h1_trigger_dir_pt;
  TH1F* h1_trigger_iso_pt;
  TH1F* h1_trigger_iso_pi0_pt;
  TH1F* h1_trigger_iso_dir_pt;
  TH1F* h1_trigger_dec_pt;
  TH1F* h1_trigger_dec_iso_pt;
  TH3F* h3_dphi;
  TH3F* h3_dphi_iso;
  TH3F* h3_dphi_dir;
  TH3F* h3_dphi_dir_iso;
  TH3F* h3_dphi_pi0;
  TH3F* h3_dphi_pi0_iso;
  TH3F* h3_dphi_mix;
  TH3F* h3_dphi_mix_iso;
  TH3F* h3_dphi_pi0_mix;
  TH3F* h3_dphi_pi0_mix_iso;

  std::vector<TH2F*> h2_dphi_dec;
  std::vector<TH2F*> h2_dphi_dec_iso;
  std::vector<TH2F*> h2_dphi_dec_mix;
  std::vector<TH2F*> h2_dphi_dec_mix_iso;

  TH2F* h2_cluster_wdR;
  TH3F* h3_cluster_dR;
  TH3F* h3_cluster_etot;
  TH2F* h2_cluster_etot;
  TH2F* h2_cluster_pi0_wdR;
  TH3F* h3_cluster_pi0_dR;
  TH3F* h3_cluster_pi0_etot;
  TH2F* h2_cluster_pi0_etot;
  TH2F* h2_cluster_dir_wdR;
  TH3F* h3_cluster_dir_dR;
  TH3F* h3_cluster_dir_etot;
  TH2F* h2_cluster_dir_etot;

};

#endif //hijing_analysis_h
