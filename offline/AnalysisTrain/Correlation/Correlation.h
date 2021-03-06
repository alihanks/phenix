#ifndef __CORRELATION_H__
#define __CORRELATION_H__

#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH3.h>
#include <TH2.h>
#include <TF1.h>
#include <PHAngle.h>
#include <CorrelationFunctions.h>

#include <APiZero.h>
#include <ACluster.h>
#include <ATrack.h>

#include <SubsysReco.h>
#include <Fun4AllHistoManager.h>

class Fun4allServer;
class PHCompositeNode;
class PHGlobal;
class emcClusterContainer;
class emcClusterContent;
class PHCentralTrack;
class TOAD;

class TFile;
class TH1;
class THmulf;
class TGraphErrors;
class TTree;

class Warnmap;

class AParticle;
class ATrack;
class AEvent;
class AMixingPool;
class AMixingTree;

const int N_ARMSECT = 8;
const int N_YPOS_PBGL = 48;
const int N_YPOS_PBSC = 36;
const int N_ZPOS_PBGL = 96;
const int N_ZPOS_PBSC = 72;
const int NTRIGBINS = 4;
const int NPARTBINS = 5;
const int NXIBINS = 6;

const int DIM = 500;

class Correlation: public SubsysReco 
{
public:
  Correlation(const char* outfile);
  virtual ~Correlation();

  virtual int Init(PHCompositeNode* topNode);
  virtual int process_event(PHCompositeNode* topNode);
  virtual int ResetEvent(PHCompositeNode* topNode);
  virtual int End(PHCompositeNode* topNode);

  void SetHadronEfficiency(const char* filename);
  void SetTriggerEfficiency(const std::string* filenames);
  void SetV2(const char* v2_inputs);
  void SetSharkFin(const char* filename);
  void SetWarnmaps(std::vector<std::string> filenames, std::vector<double> ptcuts);
  void SetZVertexCut(float cut) { zVertexCut = cut; }
  void SetNVtxBins(int nvtx) { NZVTX = nvtx; }
  void SetNCentBins(int ncent, int ncbin) { NCENT = ncent; NCBINS = ncbin; }
  void SetNMixEvents(int nmix) { NMIX = nmix; }

  void SetHadronEffFileName(std::string fn) { _hadeffFilename = fn; }
  void SetPi0EffFileName(std::string fn0, std::string fn1, std::string fn2, std::string fn3) { _pi0effFilename_0 = fn0; _pi0effFilename_1 = fn1; _pi0effFilename_2 = fn2; _pi0effFilename_3 = fn3;}
  void SetDoFillTime(int doft) { dofilltime = doft; }
  void SetWeightFileNames(std::string fn, std::string fflow) { _accfilename = fn; _flowfilename = fflow; }
  void SetSharkFinFileName(std::string fn) { _sharkfinname = fn;}
  void SetDiagFlag(int flag) { DiagFlag = flag; }
  void GetAcceptanceWeights(std::string filename);

  void Clear();

  enum DataSet { Run8pp, Run8dAu, Run10AuAu, Run11AuAu, Run14AuAu, INVALID };
  enum PairType { REAL, MIX, REALPI, MIXPI, DEC, MIXDEC, DIAGNOSTIC };
  DataSet GetDataSet(int RunNumber);
  float Rcut;
  int cbin;
  float econe_min[4];
  
  template<class T, class A> void MakePairs(std::vector<T*> triggers,
                                            std::vector<A*> associated,
                                            PairType type, DataSet dataset,
                                            int use_iso = 1,
                                            TH3F* h3dphi = NULL,
                                            TH3F* h3dphi_fold = NULL,
                                            TH3F* h3ptxidphi = NULL,
                                            TH3F* h3ptxidphi_fold = NULL,
                                            TH3F* h3ptztdphi = NULL,
                                            TH3F* h3ptztdphi_fold = NULL,
                                            std::vector<TH2F*> h2dphi_dec = std::vector<TH2F*>(),
                                            std::vector<TH2F*> h2dphi_dec_fold = std::vector<TH2F*>(),
                                            std::vector<TH2F*> h2dphixi_dec = std::vector<TH2F*>(),
                                            std::vector<TH2F*> h2dphixi_dec_fold = std::vector<TH2F*>(),
                                            std::vector<TH2F*> h2dphizt_dec = std::vector<TH2F*>(),
                                            std::vector<TH2F*> h2dphizt_dec_fold = std::vector<TH2F*>(),
                                            TH2F* h2paircut_bf = NULL,
                                            TH2F* h2paircut_af = NULL,
                                            TH2F* h2dphi = NULL, 
                                            TH2F* h2dphi_xi = NULL,
                                            TH2F* h2dphiaccw = NULL, 
                                            TH2F* h2dphiaccw_xi = NULL,
                                            TH3F* h3dphiaccw = NULL, 
                                            TH3F* h3dphiaccw_xi = NULL,
                                            TH2F* h2partptxi = NULL
                                            )
  {
    for(unsigned int it = 0; it < triggers.size(); it++){
      if( (dataset==Run8dAu||dataset==Run8pp) && !triggers[it]->IsIso() && use_iso ) continue;
      float trig_pt = triggers[it]->Pt();
      float trig_phi = PHAngle(triggers[it]->Phi());
      int itrigpt = GetPtBin(trig_pt,1);
      float mindist = 99.;
      float ptofvetotrack = 0.;
      if( verbosity > 1 )
        std::cout << "Correlation::MakePairs Event " << evt << " - found a trigger of type " << type << std::endl;

      for(unsigned int ia = 0; ia < associated.size(); ia++){
        //if(type == MIXPI) std::cout<<"ia = "<<ia<<std::endl;
        float assoc_pt = associated[ia]->Pt();
        float assoc_phi = associated[ia]->Phi();
        
        if(associated[ia]->Theta()<-99) continue;
        
        //adopted from loopthrutrack in old code -032215 this might not have any effect as the pair cut was implemented when selecting the triggers already
        float dist = 99.;
        if(assoc_pt>vetoPtCut){
          if(type==REAL) dist=FindTrackDistance((ACluster*)triggers[it], associated[ia]);
          if(type==REALPI) dist=FindTrackDistance((ACluster*)triggers[it]->Daughter1(), associated[ia]);
        }
        if(dist<8) continue;
        
        int ipartpt = GetPtBin(assoc_pt,0);

        // Apply pair cut for mixed pairs
        if( type==MIX||type==MIXPI ) {
          if( h2paircut_bf ) h2paircut_bf->Fill(itrigpt,ipartpt);
          if(type==MIX) dist=FindTrackDistance((ACluster*)triggers[it], associated[ia]);
          if(type==MIXPI) dist=FindTrackDistance((ACluster*)triggers[it]->Daughter1(), associated[ia]);
          if(dist < mindist && assoc_pt > vetoPtCut){
              mindist = dist;
              ptofvetotrack = assoc_pt;
          }
          
          if( dist < 8.0 && assoc_pt > vetoPtCut ) continue;
          if( h2paircut_af ) h2paircut_af->Fill(itrigpt,ipartpt);
        }
        
        float deltaphi = CalculateDphi(trig_phi,assoc_phi);
        float dphifold = CalculateFoldedDphi(trig_phi,assoc_phi);
        if (verbosity > 3) std::cout<<"Correlation::MakePairs - " << type << " - deltaphi = "<< deltaphi << ", deltaphiFold = " << dphifold << std::endl;
        if(dphifold<0||dphifold>PI) std::cout<<" dphifold out of bounds "<<std::endl;
        
        float zt = assoc_pt/trig_pt;
        float xi = log(1.0/zt);

        // Fill acceptance histograms with efficiency correction weighting
        float seffcorr = GetHadronEfficiencyCorr(assoc_pt);
        if( h3dphi ) h3dphi->Fill(trig_pt, assoc_pt, deltaphi, seffcorr);
        if( h3ptxidphi ) h3ptxidphi->Fill(trig_pt, xi, deltaphi, seffcorr);
        if( h3ptztdphi ) h3ptztdphi->Fill(trig_pt, zt, deltaphi, seffcorr);
        
        float filltimeflow = 1.;
        float filltimeflowxi = 1.;
        
        if( dofilltime ){
          int tbin = GetPtBin(trig_pt, 1);
          int pbin = GetPtBin(assoc_pt, 0);
          int xbin = GetXiBin(xi);
          filltimeflow = GetFilltimeWeight(type,deltaphi,assoc_pt,pbin,tbin,0,use_iso);
          filltimeflowxi = GetFilltimeWeight(type,deltaphi,assoc_pt,xbin,tbin,1,use_iso);
        }
        if( h3dphi_fold ) h3dphi_fold->Fill(trig_pt, assoc_pt, dphifold, filltimeflow);
        
        if( h3ptxidphi_fold ) {
          h3ptxidphi_fold->Fill(trig_pt, xi, dphifold, filltimeflowxi);
          if( verbosity > 1 ) std::cout<<"h3ptxidphi_fold filled."<<std::endl;
        }
        if( h3ptztdphi_fold ) {
          h3ptztdphi_fold->Fill(trig_pt, zt, dphifold, filltimeflow);
          if( verbosity > 1 ) std::cout<<"h3ptztdphi_fold filled."<<std::endl;
        }
        //filltime debugging histos
        if((trig_pt >= 5.) && (trig_pt<7.0) && DiagFlag ){
          if( h2_dphi ) h2dphi->Fill(assoc_pt,dphifold);
          if( h2dphiaccw ) h2dphiaccw->Fill(assoc_pt,dphifold,filltimeflow);
          if( h3dphiaccw ) h3dphiaccw->Fill(assoc_pt,dphifold,filltimeflow);
          if( h2dphi_xi ) h2dphi_xi->Fill(xi,dphifold);
          if( h2dphiaccw_xi ) h2dphiaccw_xi->Fill(xi,dphifold,filltimeflow);
          if( h3dphiaccw_xi ) h3dphiaccw_xi->Fill(xi,dphifold,filltimeflow);
          if( h2partptxi ) h2partptxi->Fill(assoc_pt,xi);
        }
        if(type==REAL&&DiagFlag) h3_EoverP[cbin]->Fill(assoc_pt,associated[ia]->GetEcore()/assoc_pt,dphifold);
        //******************************************
        //*  Make decay photon-h pairs             *
        //******************************************
        if(type==REALPI ) {
          if (verbosity > 1) std::cout<<"Correlation::MakePairs - making real decay pairs" << std::endl;
                  MakeDecays(DEC,use_iso,deltaphi,dphifold,assoc_pt,trig_pt,(APiZero*)triggers[it],((APiZero*)triggers[it])->GetDecayWeights(),h2dphi_dec,h2dphi_dec_fold,h2dphixi_dec,h2dphixi_dec_fold,h2dphizt_dec,h2dphizt_dec_fold);
        }
        
        if(type==MIXPI ) {
          if (verbosity > 1) std::cout<<"Correlation::MakePairs - making mixed decay pairs" << std::endl;
                  MakeDecays(MIXDEC,use_iso,deltaphi,dphifold,assoc_pt,trig_pt,(APiZero*)triggers[it],((APiZero*)triggers[it])->GetDecayWeights(),h2dphi_dec,h2dphi_dec_fold,h2dphixi_dec,h2dphixi_dec_fold,h2dphizt_dec,h2dphizt_dec_fold);
        }
      }
      if( type==MIX && DiagFlag ) {
        h3_mintrackdist_bg_allcent->Fill(mindist,trig_pt, ptofvetotrack);
        h3_mintrackdist_bg[cbin]->Fill(mindist,trig_pt, ptofvetotrack);
      }
    }
  }
  
private:
  float GetHadronEfficiencyCorr(float pt);
  float GetFilltimeWeight(PairType type, float dphi, float partpt, int pbin, int tbin, int isxi, int isiso);
  double GetAcc(TH1F* hist, float dphi);
  double GetAcceptance(PairType type, int cbin, int tbin, int pbin, float dphi, int isxi, int isiso);
  float GetFlowWeights(PairType type, int tbin, int pbin,  float partpt, float dphifold);
  //void GetXi(int decayflag, int trigptbin, int partptbin, int centbin, float& xi, float& xierr);
  TH1F* MakeDphiProjection(TH3F* h3, float xmin, float xmax, float ymin, float ymax, float norm);
  TH1F* MakeDphiProjection2D(TH2F* h2, float ymin, float ymax, float norm);
  void MakeAccHistos(TH1F* h1in, TH1F* h1out);
  int IsPbGl(int armsect){ return ((armsect==4||armsect==5)? 1 : 0); }
  int VetoTracks(ACluster* aclus, std::vector<ATrack*> lessqualtrk_vector);
  int VetoTracks(ACluster* aclus, std::vector<ATrack*> lessqualtrk_vector, float& mindist, float& vetopt);
  bool IsGoodTower(ACluster* aclus);
  bool IsGoodCluster(ACluster* aclus, DataSet dataset);
  bool IsGoodTrack(ATrack* trk, DataSet dataset);
  bool PassMatchingCut(ATrack* atrk, float pc3nsig);
  bool PassMatchingCut(ATrack* atrk, float pc3nsig, float emcnsig);
  bool Chi2Cut(emcClusterContent* sngl_emc, float zvertex, TH1F* h1_chi2 = NULL, TH3F* h3_disp = NULL);
  bool CheckPhiFiducial(float phi);
  void CalcPbglDisp(const float& m1z, const float& m1y, const float& m2z, const float& m2y, float* returnVals_z_y);
  double IncidentAngle(const int arm, const int sect, const double& fX, const double& fY, const double& fZ, const double& fBbcZEvt);
  
  void MakeEventObject(PHGlobal* glob, AEvent* aevt);
  void MakeClusterObject(emcClusterContent* clus, ACluster* aclus);
  void MakeClusterObject(ACluster* aclus, float pt, float phi, float eta, float e, float x, float y, float z, float zvtx);
  void MakeTrackObject(PHCentralTrack* trk, int itrk, ATrack* atrk);
  void MakeTrackObject(ATrack* atrk, float pt, float phi, float eta, float e, float pemcx, float pemcy, float pemcz, float zvtx);
  void MakePi0s(std::vector<ACluster*> all_clusters, int cbin, float zvertex, DataSet dataset);
  void MakePi0Object(APiZero* api0, float pt, float phi, float eta, float e, float x, float y, float z, float zvtx);
  float GetNTriggers(TH1F* trigpt, float trigptmin, float trigptmax);
  int GetPtBin(float pt, int istrig);
  int GetCentBin(int centbin);
  int GetXiBin(float xi);
  int GetPi0ZEMCBin(APiZero* pizero);
  float FindTrackDistance(ACluster* clus, ATrack* trk);
  void EvalDecWeights(APiZero* pi0trigger, float zvertex, int cbin, std::vector<float>& mwweight);
  
  void MakeDecays(PairType type, int isiso, float dphi, float dphifold, float partpt, float trigpt, APiZero* pi0, std::vector<float> weight, std::vector<TH2F*> hdphi, std::vector<TH2F*> hdphi_fold, std::vector<TH2F*> hdphixi, std::vector<TH2F*> hdphixi_fold, std::vector<TH2F*> hdphizt, std::vector<TH2F*> hdphizt_fold);
  float GetDecayXiWeights(float decxibw, int itdec, float xidec, int ipi0zemc, float trigpt, float partpt);

  void AddMBEvent(DataSet data_set);
  
  void InitHistos();
  void Init1DHisto(TH1F*& h1, std::string name, std::string xtitle, int nxbin, double xmin, double xmax);
  void Init2DHisto(TH2F*& h2, std::string name, std::string xtitle, int nxbin, double xmin, double xmax, std::string ytitle, int nybin, double ymin, double ymax);
  void Init3DHisto(TH3F*& h3, std::string name, std::string xtitle, int nxbin, double xmin, double xmax, std::string ytitle, int nybin, double ymin, double ymax, std::string ztitle, int nzbin, double zmin, double zmax);
  void Init3DHisto(TH3D*& h3, std::string name, std::string xtitle, int nxbins, double* xbins, std::string ytitle, int nybins, double* ybins, std::string ztitle, int nzbins, double* zbins);
  
  void InitPhotonCutChecker(std::string name, THmulf*& histo);
  void InitHadronCutChecker(std::string name, THmulf*& histo);
  
  void FillClusterQAHistos(int isafter, float pt, float emctrkdphi, float emctrkdz, float emcpc3dphi, float emcpc3dz);
  void FillTrackQAHistos(int isafter, float ppc3x, float ppc3y, float ppc3z, float ppc1z, float ppc1y, float zed, float phid, float pt, float quality, int n0, float pc3sdphi, float pc3sdz);
  void DoMixing(TTree* trig, TTree* assoc, int size);
  int CheckPool(int nenpart, int j, int pooldepth, int size, int& nloop);
  
  Fun4AllHistoManager* manager;
  PHGlobal* global;
  emcClusterContainer* emcclustercontainer;
  PHCentralTrack* particle;
  TOAD* toad_loader;
  
  Warnmap* warnmap;
  
  AMixingTree* atree;
  
  std::string output;
  std::vector<std::string> warnmap_filenames;
  std::vector<double> warnmap_cuts;
  
  TH1F* IncAcc[4][NTRIGBINS][NPARTBINS];//[cent][ntrig][npart]
  TH1F* Pi0Acc[4][NTRIGBINS][NPARTBINS];
  TH1F* DecAcc[4][NTRIGBINS][NPARTBINS];
  TH1F* IncAccIso[4][NTRIGBINS][NPARTBINS];//[cent][ntrig][npart]
  TH1F* Pi0AccIso[4][NTRIGBINS][NPARTBINS];
  TH1F* DecAccIso[4][NTRIGBINS][NPARTBINS];

  TH1F* IncAccXi[4][NTRIGBINS][NXIBINS];//[cent][ntrig][nxi]
  TH1F* Pi0AccXi[4][NTRIGBINS][NXIBINS];
  TH1F* DecAccXi[4][NTRIGBINS][NXIBINS];
  TH1F* IncAccXiIso[4][NTRIGBINS][NXIBINS];//[cent][ntrig][nxi]
  TH1F* Pi0AccXiIso[4][NTRIGBINS][NXIBINS];
  TH1F* DecAccXiIso[4][NTRIGBINS][NXIBINS];
  
  DataSet data_set;
  unsigned int Cut3x3Map;
  int NZVTX;
  int NCENT;
  int NCBINS;
  int NMIX;
  int evt;
  int event;
  int DiagFlag;
  int RecalFlag;
  float event_z;
  float event_c;
  float PC3_NSIGMA;
  float EMC_NSIGMA;
  float vetoPtCut;
  float minAsym;
  float zVertexCut;
  float fieldPolarity;
  int dofilltime;  
  float photon_pt_min;
  float photon_pt_max;
  float hadron_pt_min;
  float hadron_pt_max;
  float photon_ecore_min;
  float pi0_pt_min;
  float pi0_pt_max;

  TGraphErrors* gr_had_v2[4];//save hadron v2 graph for interpolating points in weighting method
  int v2sysflag;//0-no sys err; 1-sys up; 2-sys dn

  double trig_v2[3][4][NTRIGBINS]; //[type][cent][trig]
  double trig_v2_sys[3][4][NTRIGBINS];//stat and sys summed in quadrature
  double part_v2[4][NPARTBINS]; //[cent][trig]
  double part_v2_sys[4][NPARTBINS];//stat and sys summed in quadrature
  
  float num_bgtrig[3][4][NTRIGBINS];
  double meanpart[3][4][NTRIGBINS][NPARTBINS];
  
  std::vector<ACluster*> clus_everything;//to store every cluster for pi0 looping (same as in old code) -032215
  //std::vector<ACluster*> clus_vector_novetotracks;
  std::vector<ACluster*> clus_vector;
  std::vector<ACluster*> all_clus_vector;
  std::vector<ACluster*> bgclus_vector;
  std::vector<ATrack*> trk_vector;//my actual hadrons
  //std::vector<ATrack*> my_trk_vector;//matching cut at 2 sigma for hadron pt < 3 GeV, 1.5 sigma for 3<pt<5, 1.0 sigma for 5<pt<7
  std::vector<ATrack*> trk_vector_05sig;//matching cut at 0.5 sigma
  std::vector<ATrack*> trk_vector_1sig;//matching cut at 1 sigma
  std::vector<ATrack*> trk_vector_15sig;//matching cut at 1.5 sigma
  std::vector<ATrack*> trk_vector_2sig;//matching cut at 2 sigma
  std::vector<ATrack*> trk_vector_3sig;//matching cut at 2 sigma
  std::vector<ATrack*> lessqualtrk_vector;//quality>7 tracks
  std::vector<APiZero*> pi0_vector;
  std::vector<APiZero*> bgpi0_vector;
  
  TH1F* h1_n0;//n0 distribution
  TH1F* h1_chi2[N_ARMSECT];//pbsc
  TH3F* h3_disp[N_ARMSECT];//pbgl
  TH1F* h1_trigger[N_ARMSECT];//ert trigger
  
  std::vector<TH3F*> h3_dphi;
  std::vector<TH3F*> h3_dphi_iso;
  std::vector<TH3F*> h3_dphi_mix;
  std::vector<TH3F*> h3_dphi_iso_mix;
  std::vector<TH3F*> h3_dphi_fold;
  std::vector<TH3F*> h3_dphi_iso_fold;
  std::vector<TH3F*> h3_dphi_mix_fold;
  std::vector<TH3F*> h3_dphi_mix_fold_v2up;
  std::vector<TH3F*> h3_dphi_mix_fold_v2dn;
  std::vector<TH3F*> h3_dphi_iso_mix_fold;
  std::vector<TH3F*> h3_dphi_pi0;
  std::vector<TH3F*> h3_dphi_pi0_iso;
  std::vector<TH3F*> h3_dphi_pi0_mix;
  std::vector<TH3F*> h3_dphi_pi0_iso_mix;
  std::vector<TH3F*> h3_dphi_pi0_fold;
  std::vector<TH3F*> h3_dphi_pi0_iso_fold;
  std::vector<TH3F*> h3_dphi_pi0_mix_fold;
  std::vector<TH3F*> h3_dphi_pi0_mix_fold_v2up;
  std::vector<TH3F*> h3_dphi_pi0_mix_fold_v2dn;
  std::vector<TH3F*> h3_dphi_pi0_iso_mix_fold;
  
  std::vector<TH3F*> h3_ptxidphi;
  std::vector<TH3F*> h3_ptxidphi_iso;
  std::vector<TH3F*> h3_ptxidphi_fold;
  std::vector<TH3F*> h3_ptxidphi_iso_fold;
  std::vector<TH3F*> h3_ptxidphi_mix;
  std::vector<TH3F*> h3_ptxidphi_iso_mix;
  std::vector<TH3F*> h3_ptxidphi_mix_fold;
  std::vector<TH3F*> h3_ptxidphi_mix_fold_v2up;
  std::vector<TH3F*> h3_ptxidphi_mix_fold_v2dn;
  std::vector<TH3F*> h3_ptxidphi_iso_mix_fold;
  std::vector<TH3F*> h3_ptxidphi_pi0;
  std::vector<TH3F*> h3_ptxidphi_pi0_iso;
  std::vector<TH3F*> h3_ptxidphi_pi0_fold;
  std::vector<TH3F*> h3_ptxidphi_pi0_iso_fold;
  std::vector<TH3F*> h3_ptxidphi_pi0_mix;
  std::vector<TH3F*> h3_ptxidphi_pi0_iso_mix;
  std::vector<TH3F*> h3_ptxidphi_pi0_mix_fold;
  std::vector<TH3F*> h3_ptxidphi_pi0_mix_fold_v2up;
  std::vector<TH3F*> h3_ptxidphi_pi0_mix_fold_v2dn;
  std::vector<TH3F*> h3_ptxidphi_pi0_iso_mix_fold;

  std::vector<TH3F*> h3_ptztdphi;
  std::vector<TH3F*> h3_ptztdphi_fold;
  std::vector<TH3F*> h3_ptztdphi_mix;
  std::vector<TH3F*> h3_ptztdphi_mix_fold;
  std::vector<TH3F*> h3_ptztdphi_pi0;
  std::vector<TH3F*> h3_ptztdphi_pi0_fold;
  std::vector<TH3F*> h3_ptztdphi_pi0_mix;
  std::vector<TH3F*> h3_ptztdphi_pi0_mix_fold;

  std::vector<TH3F*> h3_dphi_pi0_bg;
  std::vector<TH3F*> h3_dphi_pi0_bg_mix;
  std::vector<TH3F*> h3_dphi_1sig;
  std::vector<TH3F*> h3_dphi_2sig;
  std::vector<TH3F*> h3_dphi_fold_1sig;
  std::vector<TH3F*> h3_dphi_fold_2sig;
  
  TH3D* h3_nhit[N_ARMSECT];
  TH3D* h3_nhit_mywarn[N_ARMSECT];
  
  std::vector<TH2F*> h2_pi0mass;//cent
  TH2F* h2_pi0mass_as[N_ARMSECT];
  TH2F* h2_pi0mass_PbGl;
  TH2F* h2_pi0mass_PbSc;

  TH2F* h2_pi0_bg;
  TH2F* h2_pi0_bg_as[N_ARMSECT];
  TH2F* h2_pi0_bg_PbGl;
  TH2F* h2_pi0_bg_PbSc;

  TH1F* h1_trig_pt_inc_tot;//to count # of triggers(partners)
  TH1F* h1_trig_pt_pi0_tot;
  TH1F* h1_trig_pt_dec_tot;
  TH1F* h1_part_pt_tot;
  TH1F* h1_part_pt_05sig;
  TH1F* h1_part_pt_1sig;
  TH1F* h1_part_pt_15sig;
  TH1F* h1_part_pt_2sig;
  TH1F* h1_part_pt_3sig;

  TH1F* h1_trig_pt_inc_mix_tot;
  TH1F* h1_trig_pt_pi0_mix_tot;
  TH1F* h1_trig_pt_dec_mix_tot;
  std::vector<TH1F*> h1_trig_pt_all;
  std::vector<TH1F*> h1_trig_pt_all_iso;
  std::vector<TH1F*> h1_trig_pt_inc;//[centbin]
  std::vector<TH1F*> h1_trig_pt_inc_iso;
  std::vector<TH1F*> h1_trig_pt_pi0;
  std::vector<TH1F*> h1_trig_pt_pi0_iso;
  std::vector<TH1F*> h1_trig_pt_dec;
  std::vector<TH1F*> h1_trig_pt_dec_iso;
  std::vector<TH1F*> h1_trig_pt_dec_niso;
                  
  std::vector<TH1F*> h1_part_pt;
  std::vector<TH1F*> h1_trig_pt_inc_mix;
  std::vector<TH1F*> h1_trig_pt_pi0_mix;
  std::vector<TH1F*> h1_trig_pt_dec_mix;
  std::vector<TH1F*> h1_trig_pt_inc_iso_mix;
  std::vector<TH1F*> h1_trig_pt_pi0_iso_mix;
  std::vector<TH1F*> h1_trig_pt_dec_iso_mix;
  TH2F* h2_ptvscent_trig_inc;
  TH2F* h2_ptvscent_trig_pi0;
  TH2F* h2_ptvscent_part;

  TH2F* h2_trig_pt_zvtx_inc;
  TH2F* h2_part_pt_zvtx;
  TH1F* h1_phi;
  TH1F* h1_phi_mix;
  TH2F* h2_pool_counter;
  TH2F* h2_event_counter;

  //Some QA plots
  //cluster
  
  TH2F* h2_emctrkdphi_bf;
  TH2F* h2_emctrkdz_bf;
  TH2F* h2_emcpc3dphi_bf;
  TH2F* h2_emcpc3dz_bf;

  TH2F* h2_emctrkdphi_aft;
  TH2F* h2_emctrkdz_aft;
  TH2F* h2_emcpc3dphi_aft;
  TH2F* h2_emcpc3dz_aft;
  
  std::vector<TH2F*> h2_cluster_wdR;
  std::vector<TH3F*> h3_cluster_dR;
  std::vector<TH3F*> h3_cluster_etot;
  std::vector<TH2F*> h2_cluster_etot;
  std::vector<TH2F*> h2_cluster_pi0_wdR;
  std::vector<TH3F*> h3_cluster_pi0_dR;
  std::vector<TH3F*> h3_cluster_pi0_etot;
  std::vector<TH2F*> h2_cluster_pi0_etot;
  std::vector<TH2F*> h2_cluster_mix_wdR;
  std::vector<TH3F*> h3_cluster_mix_dR;
  std::vector<TH3F*> h3_cluster_mix_etot;
  std::vector<TH2F*> h2_cluster_mix_etot;
  std::vector<TH3F*> h3_iso_acc;
  std::vector<TH3F*> h3_iso_mix_acc;
  std::vector<TH3F*> h3_iso_pi0_acc;

  //track
  TH2F* h2_pc3sdphi_bf;
  TH2F* h2_pc3sdz_bf;
  TH2F* h2_pc3sdphi_aft;
  TH2F* h2_pc3sdz_aft;
     
  TH2F* h2_ppc3_west_bf;
  TH2F* h2_ppc3_east_bf;
  TH2F* h2_ppc1_west_bf;
  TH2F* h2_ppc1_east_bf;
  TH2F* h2_ppc3_west_aft;
  TH2F* h2_ppc3_east_aft;
  TH2F* h2_ppc1_west_aft;
  TH2F* h2_ppc1_east_aft;
     
  TH2F* h2_phi_zed_bf;
  TH2F* h2_phi_zed_aft;
  TH2F* h2_phi_zed_mix;
  TH2F* h2_bfpaircut_inc;
  TH2F* h2_aftpaircut_inc;
  TH2F* h2_bfpaircut_pi0;
  TH2F* h2_aftpaircut_pi0;
  TH2F* h2_EoverPvspt;
  std::vector<TH3F*> h3_EoverP;//[centbin]

  THmulf* photon_cut_check;
  THmulf* hadron_cut_check;
  
  TH1F* h1_zvertex;
  std::vector<TH1F*> h1_centrality;

  TH3F* h3_mintrackdist_fg_allcent;
  TH3F* h3_mintrackdist_bg_allcent;
  std::vector<TH3F*> h3_mintrackdist_fg;
  std::vector<TH3F*> h3_mintrackdist_bg;
  TH1D* hshark_large[5][33];
  TH2D* ptpivsptgam[33];
 
  std::vector<std::vector<TH2F*> > h2_dphi_dec;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_iso;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_mix;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_iso_mix;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_fold;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_iso_fold;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_mix_fold;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_mix_fold_v2up;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_mix_fold_v2dn;
  std::vector<std::vector<TH2F*> > h2_dphi_dec_iso_mix_fold;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_iso;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_fold;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_iso_fold;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_mix;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_iso_mix;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_mix_fold;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_mix_fold_v2up;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_mix_fold_v2dn;
  std::vector<std::vector<TH2F*> > h2_dphixi_dec_iso_mix_fold;

  std::vector<std::vector<TH2F*> > h2_dphizt_dec;
  std::vector<std::vector<TH2F*> > h2_dphizt_dec_fold;
  std::vector<std::vector<TH2F*> > h2_dphizt_dec_mix;
  std::vector<std::vector<TH2F*> > h2_dphizt_dec_mix_fold;

  std::string _accfilename;
  std::string _flowfilename;
  std::string _hadeffFilename;
  std::string _pi0effFilename_0;
  std::string _pi0effFilename_1;
  std::string _pi0effFilename_2;
  std::string _pi0effFilename_3;
  std::string _sharkfinname;
  
  TF1* fexemb;
  TF1* fhadroneff;
  TH1F* ratio;
  
  TFile* fpi0eff[4];
  TGraphErrors* grpi0eff[4];

  //debug tools
  //for type REAL
  TH2F* h2_dphi;//x-part pt; y-dphi
  TH2F* h2_dphi_xi;//x-xi; y-dphi
  TH2F* h2_dphi_accw;//fill with weights
  TH2F* h2_dphi_accw_xi;
  TH3F* h3_dphi_accw;//x-part pt; y-dphi; z-weights
  TH3F* h3_dphi_accw_xi;//x-xi; y-dphi; z-weights
  TH2F* h2_partpt_xi;//x-assoc pt; y-xi

  //for type MIX
  TH2F* h2_dphi_mix;
  TH2F* h2_dphi_xi_mix;
  TH2F* h2_dphi_accw_mix;
  TH2F* h2_dphi_accw_xi_mix;
  TH3F* h3_dphi_accw_mix;
  TH3F* h3_dphi_accw_xi_mix;
  TH2F* h2_partpt_xi_mix;

  //for type REALPI
  TH2F* h2_dphi_pi0;
  TH2F* h2_dphi_xi_pi0;
  TH2F* h2_dphi_accw_pi0;
  TH2F* h2_dphi_accw_xi_pi0;
  TH3F* h3_dphi_accw_pi0;
  TH3F* h3_dphi_accw_xi_pi0;
  TH2F* h2_partpt_xi_pi0;

  //for type MIXPI
  TH2F* h2_dphi_pi0_mix;
  TH2F* h2_dphi_xi_pi0_mix;
  TH2F* h2_dphi_accw_pi0_mix;
  TH2F* h2_dphi_accw_xi_pi0_mix;
  TH3F* h3_dphi_accw_pi0_mix;
  TH3F* h3_dphi_accw_xi_pi0_mix;
  TH2F* h2_partpt_xi_pi0_mix;

  TH3F* h3_pt_phi_eta_clus;
  TH3F* h3_pt_phi_eta_trk;
};

#endif /* __CORRELATION_H__ */
