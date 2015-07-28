#ifndef __COMBINEDSIMPLE_H__
#define __COMBINEDSIMPLE_H__
#include <string>
#include <vector>
#include <SubsysReco.h>
#include <TString.h>
#include <TFile.h>

class THmulf;
class PHCompositeNode;
class PHGlobal;
class PHCentralTrack;
class PHSnglCentralTrack;
class TNtuple;
class TH1;
class TH1F;
class TH2;
class TH3;
class TH2F;
class PhCglList;
class emcClusterContainer;
class emcClusterContent;
class PhPhotonList; 
class ReactionPlaneObject;
class RunHeader;
class TF1;
class EventHeader;
class TGraphErrors;
class Lvl2OutArray;
class SpinDataEventOut;
class TGraphErrors;
class Lvl2OutArray;
class ErtOut;
class TTree;
class TFile;

// Justin's vars
const double _MAXTRIGPT = 14.0; //14.0 for pp & 12.5 for au
const int _TRBINS = 24; // bins between Maxtrig- 2.0


using namespace std;

class CombinedSimple: public SubsysReco
{
 public:  

  CombinedSimple(char *, float, float, int, int, int, int);
  CombinedSimple() {}
  virtual ~CombinedSimple();

  virtual  int Init(PHCompositeNode *topNode);
  virtual int InitRun(PHCompositeNode *topNode) {return 0;}
  virtual int process_event(PHCompositeNode *topNode);
  virtual int Reset(PHCompositeNode *topNode) {return 0;}
  virtual void Print(const char *what) const {return;}
  virtual int ResetEvent(PHCompositeNode *topNode);
  virtual int End(PHCompositeNode *topNode); 

  TString* _cntnodename;
  TString* _emcnodename;
  TString* _emcBadMapFilename;
  TString* _pi0effFilename;
  TString* _filltimeFilename;
  TString* _filltimeFilenameup;
  TString* _filltimeFilenamedown;
  TString* _sharkfinname;

  TString* _histPrefix;
  TString  _assocFilename; // I think this is a better way

  TFile* fpi0eff;
  TFile* ffilltimecorrs;
  TFile* ffilltimesysup;
  TFile* ffilltimesysdown;

  void  set_cntnodename(char * fn) 
    { if (_cntnodename) delete _cntnodename; _cntnodename = new TString(fn);}
  void  set_emcnodename(char * fn)
    { if (_emcnodename) delete _emcnodename; _emcnodename = new TString(fn);}
  void  set_emcBadMapFilename(char * fn) 
    { if (_emcBadMapFilename) delete _emcBadMapFilename; _emcBadMapFilename = new TString(fn);}

  void  set_pi0effFilename(char * fn)
    { if (_pi0effFilename) delete _pi0effFilename; _pi0effFilename = new TString(fn);}
  void  set_filltimeFilename(char * fn) 
    { if (_filltimeFilename) delete _filltimeFilename; _filltimeFilename = new TString(fn);}
  void  set_filltimeFilenameup(char * fn) 
    { if (_filltimeFilenameup) delete _filltimeFilenameup; _filltimeFilenameup = new TString(fn);}
  void  set_filltimeFilenamedown(char * fn) 
    { if (_filltimeFilenamedown) delete _filltimeFilenamedown; _filltimeFilenamedown = new TString(fn);}

  void  set_histPrefix(char * fn) 
    { if (_histPrefix) delete _histPrefix; _histPrefix = new TString(fn);}
  void  set_assocFilename(char * fn) {_assocFilename = fn;}

  void  set_sharkfinname(char * fn) 
    { if (_sharkfinname) delete _sharkfinname; _sharkfinname = new TString(fn);}



  void set_savAssoc(short le)
    {_savAssoc = le;}

  void set_anticuts(short le)
    {_anticuts = le;}


  TTree * get_assocNtuple()
    {return _assocNtuple;}

  TNtuple * get_triggaz()

    {return triggaz;}

  void set_useLessHistos(int inflag)
    {_useLessHistos = inflag;}

  void set_vetoPtCut(float inpt)
    {_vetoPtCut = inpt;} // default is 1.0 (GeV)

  void set_removetags(int inflag)
    {removetags = inflag;}

  void set_dofilltime(int inflag)
    {dofilltime = inflag;}

  void set_pouthistos(int inflag)
    {pouthistos = inflag;}

  void set_pouthistos_more(int inflag)
    {pouthistos_more = inflag;}

  void set_makepi0cone(int inflag)
    {makepi0cone = inflag;}


  int twr2id(int armsect,int ytwr,int ztwr);
  void SetTagFlag(int);
  void SetSpeciesFlag(int);
  float SumEcorePtInCone(float, float, float, int, int, float, float);
  //float SumConeE(float, float, float, int, int, float, float);

 protected:

  void SetTriggerEfficiency(const char* filename);
  void SetSharkFin(const char* filename);
  void SetFilltimeCorrs(const char* filltime_name, const char* filltimeup_name, const char* filltimedown_name);

  struct PartnerHolder{
    //keep them float so we can put them in ntuple
    //these are properties of the photon OR pi0
    //int size;
    float vthe0;
    float vmom;
    float vpc3dphi;
    float vpc3dz;
    float vemcdphi;
    float vemcdz;
    int vn0;
    float vm2emc;
    float vm2tof;
    float vphi0;
    int vdcarm;
    int vcharge;
    int vquality;
    float vppc3x;
    float vppc3y;
    float vppc3z;
    float vpc3sdz;
    float vpc3sdphi;
    float vemcsdz;
    float vemcsdphi;

    PartnerHolder(){
      //     E=0;	       px=0;	       py=0;	       pz=0;	       the=0;	       eta=0;	       phi=0;	       z=0;	       arm=0;	       pt=0;	       x=0;	          y=0;	         z=0;	       sector=0;    pc3dr = 0;   emctof=0;

     vthe0=0;
     vmom=0;
     vpc3dphi=0;
     vpc3dz=0;
     vemcdphi=0;
     vemcdz=0;
     vn0=0;
     vm2emc=0;
     vm2tof=0;
     vphi0=0;
     vdcarm=0;
     vcharge=0;
     vquality=0;
     vppc3x=0;
     vppc3y=0;
     vppc3z=0;
     vpc3sdz=0;
     vpc3sdphi=0;
     vemcsdz=0;
     vemcsdphi=0;


    }
  };


  int species; 
  PHGlobal      *global;     
  PHCentralTrack      *particle;     
  emcClusterContainer *emccluster;
  emcClusterContent *photon, *photon2;
  ReactionPlaneObject *reacplane;
  RunHeader *d_runhdr;
  EventHeader *evtsync;


  bool Chi2Cut(emcClusterContent* sngl_emc, float zvertex);

  struct TriggerHolder{
    //keep them float so we can put them in ntuple
    //these are properties of the photon OR pi0
    float E;
    float px;
    float py;
    float pz;
    float the;
    float eta;
    float phi;
    float arm;
    float pt;
    //these are properties of the leading gamma
    float x;
    float y;
    float z;
    float sector;
    float pc3dr;
    float emctof;
    float zvertex;
    
    TriggerHolder(){
      E=0;	       px=0;	       py=0;	       pz=0;	       the=0;	       eta=0;	       phi=0; 	       arm=0;	       pt=0;	       x=0;	          y=0;	         z=0;	       sector=0;    pc3dr = 0;   emctof=0; zvertex=0;
    }
  };
  
  float LoopThruTracks(TriggerHolder, PHCentralTrack *, emcClusterContainer *emc = 0, int thres=0);
  int LoopThruTrack(TriggerHolder, PartnerHolder, int tag_thres=0);
  int VetoTracks(TriggerHolder, PHCentralTrack *, emcClusterContainer *emc = 0, int thres=0);

  short _anticuts;   
  short _savAssoc; 
   
  TTree * _assocNtuple;
  short _assocNum;

  static const int asNt=80;
  //  PartnerHolderPhoton _assocHolder[asNt];
  //PartnerHolder _assocHolderTrk[asNt];

  void initAssocNtuple();
  
  void makeMNSingles(PHCentralTrack *, int doNtup = 1); 
  //  int evtno, int seqno, int runno, 
  
  Int_t as_runno;  
  Int_t as_segno;  
  Long64_t as_evtno;  
  Float_t as_zvertex;
  Int_t as_centrality;
  Int_t as_N;

  Float_t as_px[asNt];     
  Float_t as_py[asNt];     
  Float_t as_pz[asNt];  
  Float_t as_pt[asNt]; // not to be saved in ntup
  Float_t as_ppc3x[asNt];  
  Float_t as_ppc3y[asNt];  
  Float_t as_ppc3z[asNt];  
  Float_t as_ecore[asNt];  
  Int_t as_charge[asNt]; 
  Float_t as_m2emc[asNt];  
  Float_t as_m2tof[asNt];  
  //  Int_t as_quality[asNt];
  Float_t as_n0[asNt];     

  int useiso;  // set for isolation cut study
  int removetags;  // set for tagging study
  int dofilltime;  // set for using fill time weights

  float FindTrackDistance(TriggerHolder, float xtrk, float ytrk, float ztrk, float& dtheta, float& dphi);
  float FindTrackDistance(TriggerHolder,PartnerHolder, float& dtheta, float& dphi);

  void closestApproachPPP(float zvertex, float p1x, float p1y, float p1z,float p2x, float p2y, float p2z,float& p3x, float& p3y, float& p3z);
 
  
  float foldphi(float inphi); // justins way

  int CheckNewBadTowers(int, int, int);  
  void InitializeBadRunList();
  int CheckBadRunList(int);

   
  void CalcPbglDisp(const float& m1z, const float& m1y,
		    const float& m2z, const float& m2y,
		    float* returnVals_z_y
		    );
  
  double IncidentAngle(const int arm,
		       const int sect,
		       const double& fX, 
		       const double& fY,
		       const double& fZ,
		       const double& fBbcZEvt);

  //Methods for pp
  int CheckLvl1Fired(emcClusterContent * emc, ErtOut * Ert, int * ertTrigs);
  //  int CheckLvl2Fired(emcClusterContent * phot, Lvl2OutArray * lvl2outarray);
  void ReadHotMap();
  void ReadMattHotMap();
  int GetTowerEnergyCorrections();
  float GetEcore(emcClusterContent *);
  
  void InitializeHotTowerMap(int species = 0, int addTowers = 0);
  int CheckHotTowerMap(int sect, int iz, int iy, int do3x3);
  int isEdgeTower(int sect, int iz, int iy);

  //run8 hot tower
  bool isHotTower(int towerid, int species);
   

  vector <int> HotTowers;

  
  int event_counter;
  int badevent_counter;
  
  char *ntuple_filename;
  TFile *ntuple_file;
  
  int nevent;
  int skip_fg_histos;
  int _useLessHistos;

  int runno,lastrun;
    
  float locent, hicent;
  int tagflag, run_5_6;

  float _vetoPtCut;

  float poutlothresh, pouthithresh;
  int pouthistos;
  int pouthistos_more;
  int makepi0cone;


  TH1 * tilespectra;
  TH1 * lvl1failed;
  

  TGraphErrors *greffman;
  TGraphErrors *greff;

  TGraphErrors *grpi0eff;
  TGraphErrors *grpi0efflo;
  TGraphErrors *grpi0effhi;

  TH2 *PI0ZEMC;

  TH2 *PC3SDPHI;
  TH2 *PC3SDZ;
  TH3 *PT_PC3SDPHI_PC3SDZ;
  
  TH2 * AWPHIPT[5];

  TH2 * MWPHIPT[5];
  TH2 * MWPHIPTPLUS[5];
  TH2 * MWPHIPTMINUS[5];
  TH2 * MWPHIPT_MES_P_TOF[5];
  TH2 * MWPHIPT_MES_P_EMC[5];
  TH2 * MWPHIPT_MES_M_TOF[5];
  TH2 * MWPHIPT_MES_M_EMC[5];
  TH2 * MWPHIPT_BAR_P_TOF[5];
  TH2 * MWPHIPT_BAR_P_EMC[5];
  TH2 * MWPHIPT_BAR_M_TOF[5];
  TH2 * MWPHIPT_BAR_M_EMC[5];

  TH2 * MWPHIFOLDPT[5];
  TH2 * MWPHIFOLDPTQ2[5];
  TH2 * MWPHIFOLDPTPLUS[5];
  TH2 * MWPHIFOLDPTMINUS[5];
  TH2 * MWPHIFOLDPT_MES_P_TOF[5];
  TH2 * MWPHIFOLDPT_MES_P_EMC[5];
  TH2 * MWPHIFOLDPT_MES_M_TOF[5];
  TH2 * MWPHIFOLDPT_MES_M_EMC[5];
  TH2 * MWPHIFOLDPT_BAR_P_TOF[5];
  TH2 * MWPHIFOLDPT_BAR_P_EMC[5];
  TH2 * MWPHIFOLDPT_BAR_M_TOF[5];
  TH2 * MWPHIFOLDPT_BAR_M_EMC[5];

  TH2 * TRIG_LOCATION;
  TH2 * ASSOC_LOCATION;
  TH2 * PI0_LOCATION;
  TH2 * INVMASS[8];
  TH2 * DECINVMASS;
  TH1 * PI0_PT;
  TH1 * TRIG_COUNTER;
  TH1 * RPLANE;
  TH1 * CENTRALITY;

  //add by hge
  TH2 *PT_QUAL;
  TH2 *PT_N0;
  TH1 *PT_PC3SDPHISDZ;
  TH1 *PT_AFTN0;
  TH1 *PT_AFTN0QUAL;
  TH1 *PTAFTCUTS;
  THmulf *HADRONCUTCHECK;
  TH3 *DCH_N0;
  TH3 *DCH_N0_QUAL;
  TH3 *DCH_N0_QUAL_PC3;
  TH1 *PT_PC3SDPHI_UNDEFINE;
  TH1 *PT_PC3SDZ_UNDEFINE;
  TH1 *PT_PC3DPHI_UNDEFINE;
  TH1 *PT_PC3DZ_UNDEFINE;
 
  TH3 *MINTRACKDIST;
  TH3 *MINRADDIST;
  TH3 *MINPHIDIST;
  TH3 *MINTHETADIST;
  TH3 *DCH;



  TH2 *DPHIPT_NORM[4];
  TH2 *DPHIPT_DECNORM[4];
  TH2 *DPHIPT_FLOW[4];
  TH2 *DPHIPT_DECFLOW[4];
  TH2 *DPHIPT_NORM_PLUS[4];
  TH2 *DPHIPT_DECNORM_PLUS[4];
  TH2 *DPHIPT_NORM_MINUS[4];
  TH2 *DPHIPT_DECNORM_MINUS[4];
  TH2 *DPHIPT_FLOW_PLUS[4];
  TH2 *DPHIPT_DECFLOW_PLUS[4];
  TH2 *DPHIPT_FLOW_MINUS[4];
  TH2 *DPHIPT_DECFLOW_MINUS[4];

  double pout2sum[4][7];
  double pout2weight[4][7];
  double decpout2sum[4][7];
  double decpout2weight[4][7];
  double pout2sum_large[4][4];
  double pout2weight_large[4][4];
  double decpout2sum_large[4][4];
  double decpout2weight_large[4][4];

  TH2 * XEDIST[4][7];
  TH2 * ZTDIST[4][7];
  TH2 * XIDIST[4][7];
  TH2 * POUTDIST[4][7];
  TH2 * POUT2DIST[4][7];
  TH1 * POUT2SUM[4][7];
  TH2 * XEVSPOUT2[4][7];
  TH2 * ZTVSPOUT2[4][7];
  TH2 * DECXEDIST[4][7];
  TH2 * DECZTDIST[4][7];
  TH2 * DECXIDIST[4][7];
  TH2 * DECPOUTDIST[4][7];
  TH2 * DECPOUT2DIST[4][7];
  TH1 * DECPOUT2SUM[4][7];
  TH2 * DECXEVSPOUT2[4][7];
  TH2 * DECZTVSPOUT2[4][7];


  TH2 * XEDIST_LARGE[4][4];
  TH2 * XIDIST_LARGE[4][4];
  TH2 * ZTDIST_LARGE[4][4];

  TH2 * POUTDIST_LARGE[4][4];
  TH2 * POUT2DIST_LARGE[4][4];
  TH1 * POUT2SUM_LARGE[4][4];
  TH2 * DECXEDIST_LARGE[4][4];
  TH2 * DECXIDIST_LARGE[4][4];
  TH2 * DECZTDIST_LARGE[4][4];
  TH2 * DECPOUTDIST_LARGE[4][4];
  TH2 * DECPOUT2DIST_LARGE[4][4];
  TH1 * DECPOUT2SUM_LARGE[4][4];


  TH2 * XEDIST_STATERR[4][4];
  TH2 * XIDIST_STATERR[4][4];
  TH2 * ZTDIST_STATERR[4][4];

  TH2 * POUTDIST_STATERR[4][4];
  TH2 * POUT2DIST_STATERR[4][4];
  TH1 * POUT2SUM_STATERR[4][4];
  TH2 * DECXEDIST_STATERR[4][4];
  TH2 * DECXIDIST_STATERR[4][4];
  TH2 * DECZTDIST_STATERR[4][4];
  TH2 * DECPOUTDIST_STATERR[4][4];
  TH2 * DECPOUT2DIST_STATERR[4][4];
  TH1 * DECPOUT2SUM_STATERR[4][4];


  TH2 * XEDIST_LARGE_NOW[2][4];
  TH2 * XIDIST_LARGE_NOW[2][4];
  TH2 * ZTDIST_LARGE_NOW[2][4];

  TH2 * DECXEDIST_LARGE_NOW[2][4];
  TH2 * DECXIDIST_LARGE_NOW[2][4];
  TH2 * DECZTDIST_LARGE_NOW[2][4];

  TH2 * DECPT2ZT[4];
  TH2 * DECPT2XI[4];



  TH2 * XEDIST_LARGE_PLUS[4][4];
  TH2 * XIDIST_LARGE_PLUS[4][4];
  TH2 * ZTDIST_LARGE_PLUS[4][4];
  TH2 * ZTDIST_PLUS[4][4];
  TH2 * POUTDIST_LARGE_PLUS[4][4];
  TH2 * DECXEDIST_LARGE_PLUS[4][4];
  TH2 * DECXIDIST_LARGE_PLUS[4][4];
  TH2 * DECZTDIST_LARGE_PLUS[4][4];
  TH2 * DECZTDIST_PLUS[4][4];
  TH2 * DECPOUTDIST_LARGE_PLUS[4][4];

  TH2 * XEDIST_LARGE_MINUS[4][4];
  TH2 * XIDIST_LARGE_MINUS[4][4];
  TH2 * ZTDIST_LARGE_MINUS[4][4];
  TH2 * ZTDIST_MINUS[4][4];
  TH2 * POUTDIST_LARGE_MINUS[4][4];
  TH2 * DECXEDIST_LARGE_MINUS[4][4];
  TH2 * DECXIDIST_LARGE_MINUS[4][4];
  TH2 * DECZTDIST_LARGE_MINUS[4][4];
  TH2 * DECZTDIST_MINUS[4][4];
  TH2 * DECPOUTDIST_LARGE_MINUS[4][4];


  TH2 * PC3S_DPHI_PT_PLUS;
  TH2 * PC3S_DZ_PT_PLUS;
  TH2 * PC3S_DPHI_PT_MINUS;
  TH2 * PC3S_DZ_PT_MINUS;


  TH3 *PT1PT2DPHIFOLDCF;
  TH3 *PT1PT2DPHIFLOW;

  TH2 *MWPHIFOLDPTFLOW[5];
  TH2 *MWPHIFOLDPTCF[5];


  TH3 *PT1PT2ZT;
  TH3 *PT1PT2XI;


  TH2 *M2TOF;
  TH2 *M2EMC;

  TH2 *ptvscent_trig;
  TH2 *ptvscent_part;
  TH2 *ptvscent_dec;

  TH1 *VETO_COUNTER;
  TH1 *VETO_COUNTER_0;

  TH1 *TRIGPT;
  TH1 *TRIGPTEFFW;
  TH1 *TRIGPT1;
  TH1 *TRIGPT2;
  TH1 *TRIGPT3;
  TH1 *TRIGPT4;
  TH1 *TRIGPTSC;
  TH1 *TRIGPTGL;

  TH1 *AWTRIGPT;
  TH1 *MWTRIGPT;
  TH1 *MWFINETRIGPT;

  TH3 *PT1PT2DPHI;
  TH3 *PT1PT2DPHI1;
  TH3 *PT1PT2DPHI2;
  TH3 *PT1PT2DPHI3;
  TH3 *PT1PT2DPHI4;

  TH3 *PiCone;
  TH2 *trigptvse;
  TH1 *EvsR1;
  TH1 *EvsR2;
  TH1 *EvsR3;
  TH1 *EvsR4;

  TH2 *ISOMAP1;
  TH2 *ISOMAP2;
  TH2 *ISOMAP3;
  TH2 *ISOMAP4;

  TH2 *ISOMAPSH1;
  TH2 *ISOMAPSH2;
  TH2 *ISOMAPSH3;
  TH2 *ISOMAPSH4;


  TH3 *PT1PT2DPHIFOLD;
  TH3 *PT1PT2DPHIFOLDGL;
  TH3 *PT1PT2DPHIFOLDSC;

  TH3 *PT1PT2DPHIPLUS;
  TH3 *PT1PT2DPHIMINUS;

  TH3 *PT1PT2DPHIFOLDPLUS;
  TH3 *PT1PT2DPHIFOLDMINUS;

  TH3 *PT1PT2DPHI_FINE;
  TH3 *PT1PT2DPHI_E;
  TH3 *PT1PT2DPHI_W;
  TH3 *PT1PT2DPHI_BAR_P_TOF;
  TH3 *PT1PT2DPHI_BAR_P_EMC;
  TH3 *PT1PT2DPHI_BAR_M_TOF;
  TH3 *PT1PT2DPHI_BAR_M_EMC;
  TH3 *PT1PT2DPHI_MES_P_TOF;
  TH3 *PT1PT2DPHI_MES_P_EMC;
  TH3 *PT1PT2DPHI_MES_M_TOF;
  TH3 *PT1PT2DPHI_MES_M_EMC;

  TH3 *PT1PT2DPHIFOLD_BAR_P_TOF;
  TH3 *PT1PT2DPHIFOLD_BAR_P_EMC;
  TH3 *PT1PT2DPHIFOLD_BAR_M_TOF;
  TH3 *PT1PT2DPHIFOLD_BAR_M_EMC;
  TH3 *PT1PT2DPHIFOLD_MES_P_TOF;
  TH3 *PT1PT2DPHIFOLD_MES_P_EMC;
  TH3 *PT1PT2DPHIFOLD_MES_M_TOF;
  TH3 *PT1PT2DPHIFOLD_MES_M_EMC;

  TH3 *PT1PT2DPHI_BG_NORM;
  TH1 *heff;
  TF1 *feff;
  TF1 *fexemb;


  TH3 *MODMAP[8];
  TH3 *MODMAPSUB[8];

  TH1 *ASSYM57;
  TH1 *ASSYM79;
  TH1 *ASSYM912;
  TH1 *ASSYM1220;

  TH1F *TWRID;
  TH1F *TWRIDMYMAP;

  TH1F *PHILEAD;
  TH1F *PHISUB;

  
  TH1 *ZVERTEX;


  TH1 *TAGCOUNTER; 
  TH1 *TAGCOUNTER1; 
  TH1 *TAGCOUNTER2; 
  
  TH2 *pivsgampi;

  TH2 *ptpivsptgam[33];

  TH1 *hshark_large[5][33];
  TH1 *hshark_small[7][33];
  TH1 *hshark_alt[7][33];



  TNtuple * triggaz;

  int HotTower[82];
  int HotTowerChris[6175];
  int HotTowerMatt[61];
  float tower_corr[8][48][96];
  

};

// declare Wei's lvl2 ert class as a private (sub) class
class smTileModule
{
 private:
  bool debug;
 public:
  smTileModule() {debug = false;}
  ~smTileModule() {;}
  bool   get_smID(int towerkey, int &arm, int &sector, int &smID);
  int  get_towerkey(int arm, int sector, int iy, int iz);
  int  GetModuleFromTower(int towerkey);
  int  GetEmcttFromModule(int moduleid);
  void GetSMidFromEMCTrgTile(int &arm, int &sector,int &smID, int emctt);
  
  void set_debug(bool input) { debug = input;}
  void GetSMidFromCrkTrgTile(int &arm, int &sector, int &smID, int crktt);
  int GetCrkttFromPMT(int pmt);
  void get_CRKsmID(int pmt, int &armRICH, int &sectRICH,int &smRICH);
  

};

#endif /* __COMBINEDSIMPLE_H__ */
