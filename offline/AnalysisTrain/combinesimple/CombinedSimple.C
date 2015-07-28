
/*
Looks at the away-side from a high energy hadron or photon.
This code is designed to run on (justin's) filtered hard pDSTs
-man 6/2/05
tagflag --
0 - inclusive
1 - pi0
2 - 'direct'
10/31/05

species flag
0 pp
1 run 7 AuAu
2 run 4 AuAu
3 run 8 dAu MB
4 run 8 dAu ERT
5 run 10 AuAu
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <PHGlobal.h>
#include <PhCglList.h>
#include <PHCentralTrack.h>
#include <PHCentralTrackv23.h>
#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <Fun4AllServer.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TFile.h>
#include <CombinedSimple.h>
#include <PHObject.h>
#include <PHAngle.h>
#include <PHSnglCentralTrack.h>
#include <emcClusterContainer.h>
#include <emcClusterContent.h>
#include <utiCentrality.h>
#include <TLorentzVector.h>
#include <ReactionPlaneObject.h>
#include <getClass.h>
#include <TF1.h>
#include <THmulf.h>
#include <RunHeader.h>
#include <EventHeader.h>
#include <Lvl2OutArray.h>
#include <SpinDataEventOut.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <stdio.h>
#include <stdlib.h>
#include <TriggerHelper.h>
#include <ErtOut.h>
#include <TOAD.h>

typedef PHIODataNode <PHObject> PHObjectNode_t;
typedef PHIODataNode <PHGlobal> PHGlobalNode_t;
typedef PHIODataNode <PHCentralTrack> PHPCentralNode_t;
typedef PHIODataNode <PhCglList> PHCglNode_t;
typedef PHIODataNode <emcClusterContainer> PHEmcNode_t;


using namespace std;
using namespace findNode;

const double DEG_PER_RAD = 180.0 / M_PI;
const double pi = M_PI;
const double JFMYPI = M_PI;


//some constants for the pbgl
static const float vnx[2][4] = { { 0.92619,     0.999986,     0.921968,     0.702254 }, {-0.92627,     -0.999996,    -0.921187,    -0.702855 } };
static const float vny[2][4] = { {-0.377058,    0.00528539,   0.387266,     0.711927 }, {-0.37686,      0.00268898,   0.389118,     0.711332 } };
static const float vnz[2][4] = { { 4.61055e-05, 7.46704e-05, -0.000252431, -0.000332109}, {-0.000400881, -0.000738336, -0.000926837, -0.000577316} };

vector <int> sbadruns; //for Hpdst

CombinedSimple::CombinedSimple(char *name_ntuple_file, float centlo, float centhi, int flagtag, int setspecies, int isMixEv, int set_run_5_6)
{
  //many of these flags are not currently in use but may be re-implemented again in the future
  ThisName = "CombinedSimple";

  fexemb = NULL;
  ntuple_filename = NULL;
  if (name_ntuple_file)
  {
    ntuple_filename = new char[1000];
    sprintf(ntuple_filename, "%s", name_ntuple_file);
  }

  locent=centlo;
  hicent=centhi;
  tagflag=flagtag;

  skip_fg_histos=isMixEv;
  _useLessHistos = 0;
  _vetoPtCut = 1.0;
  removetags = 0;
  dofilltime = 0;
  pouthistos=0;
  pouthistos_more=0;

  run_5_6 = set_run_5_6;

  //defaults 
  _cntnodename = new TString("PHCentralTrack");
  _emcnodename = new TString("emcClusterContainer");
  _emcBadMapFilename = 0;
  _pi0effFilename = 0;
  _filltimeFilename = 0;
  _filltimeFilenameup = 0;
  _filltimeFilenamedown = 0;
  _sharkfinname = 0;
  _histPrefix = 0;  
  fpi0eff=NULL;
  ffilltimecorrs=NULL;
  ffilltimesysup=NULL;
  ffilltimesysdown=NULL;
  _savAssoc = 0;
  _assocNtuple = 0;
  _assocNum = -1;
  _anticuts = 0;

  species = setspecies; //0 for pp & 1 for run7 AuAu & 2 for run4 AuAu & 3 for run8dAu MB & 4 for run8dAu ERT#CHC
  cout <<PHWHERE<< "species= " << species <<endl;

  //For ztdist stuff
  poutlothresh=2.0;
  //poutlothresh=0.5;
  pouthithresh=7.0;
  
}

CombinedSimple::~CombinedSimple()
{
  if (_cntnodename) delete _cntnodename;
  if (_sharkfinname) delete _sharkfinname;
  if (_emcnodename) delete _emcnodename;
  if (_emcBadMapFilename) delete _emcBadMapFilename;
  if (_pi0effFilename) delete _pi0effFilename;
  if (_filltimeFilename) delete _filltimeFilename;
  if (_filltimeFilenameup) delete _filltimeFilenameup;
  if (_filltimeFilenamedown) delete _filltimeFilenamedown;
  if (_histPrefix) delete _histPrefix;
  if( ntuple_filename) delete[] ntuple_filename;
  if(fexemb) delete fexemb;
  if( fpi0eff ){
    fpi0eff->Close();
    delete fpi0eff;
  }
  if( grpi0eff ) delete grpi0eff;
  if( ffilltimecorrs ){
    ffilltimecorrs->Close();
    delete ffilltimecorrs;
  }
  if( ffilltimesysdown ){
    ffilltimesysdown->Close();
    delete ffilltimesysdown;
  }
  if( ffilltimesysup ){
    ffilltimesysup->Close();
    delete ffilltimesysup;
  }
}

int CombinedSimple::Init(PHCompositeNode *topNode)
{
  // Initialize input file loader
  TOAD *toad_loader = new TOAD("combinesimple");

  useiso=0;  

  for(int iside=0;iside<4;iside++){
    for(int ipt=0;ipt<7;ipt++){
      pout2sum[iside][ipt]=0.;
      pout2weight[iside][ipt]=0.;
      decpout2sum[iside][ipt]=0.;
      decpout2weight[iside][ipt]=0.;
    }
  }
  for(int iside=0;iside<4;iside++){
    for(int ipt=0;ipt<4;ipt++){
      pout2sum_large[iside][ipt]=0.;
      pout2weight_large[iside][ipt]=0.;
      decpout2sum_large[iside][ipt]=0.;
      decpout2weight_large[iside][ipt]=0.;
    }
  }


  if(_anticuts) cout <<"!!!!!!!!!!WARNING!!!!!!!!! USING ANTI CUTS !!!!!!!!!!! " << _anticuts <<endl;
  //Load hot tower map

  if(species==1||species==2) InitializeHotTowerMap(2,1);
  if(species==5) InitializeHotTowerMap(2,0);

  Fun4AllServer * se = Fun4AllServer::instance();
  
  cout<<PHWHERE<<" running CombinedSimple for: "<<" locent "<<locent<<" hicent "<<hicent<<" tagflag "<<tagflag<<" species "<<species<<endl;
  

  vector <string> str_cent;
  for(int i=0;i<4;i++){ 
    ostringstream centnamer;
    centnamer<<i; 
    str_cent.push_back(centnamer.str());    
  }


  int centbin=0;
  if (locent==0) centbin=0;
  if (locent==20) centbin=1;
  if (locent==40) centbin=2;
  if (locent==60) centbin=3;

  string centname="C"+str_cent[centbin];
  
  if (_histPrefix) centname = _histPrefix->Data() + centname;
     
     
  event_counter=0;  

  //book ntuples, histos
  
  cout<<PHWHERE<< " booking histos "<<endl;

  if(!skip_fg_histos){
    triggaz=new TNtuple("triggaz","triggaz","runno:seqno:evtno:E:Esub:eta:the:phi:vertex:ngoodpi0s:px:py:pz:x:y:z:sector:pc3dr:emctof:percent:thetaDEG:towerid1:towerid2");


    string triglocation = centname+"_TRIG_LOCATION";
    se->registerHisto(triglocation.c_str(), TRIG_LOCATION = new TH2F(triglocation.c_str(),triglocation.c_str(),100,-0.5,0.5,100, -1.0, 5.2));
 
    /*   
    string assoclocation = centname+"_ASSOC_LOCATION";
    se->registerHisto(assoclocation.c_str(),ASSOC_LOCATION = new TH2F(assoclocation.c_str(),"assoc location",100,-0.5,0.5,100,-1.0, 5.2));
    */
    //  se->registerHisto("DETADPHI",DETADPHI = new TH2F("DETADPHI","assoc location wrt trigger axis",30,0,1.5,170,-0.1 ,1.6));
    
    if(tagflag>0){
      string pilocation = centname+"_PI0_LOCATION";
      se->registerHisto(pilocation.c_str(),PI0_LOCATION = new TH2F(pilocation.c_str(),"pi0 location",1000,-0.5,0.5,1000, -1.0, 5.2));
      string pipt = centname+"_PI0_PT";
      se->registerHisto(pipt.c_str(), PI0_PT = new TH1F(pipt.c_str(),"pi0 pt",600,4,10));
    }
    string trigcount = centname+"_TRIG_COUNTER";
    se->registerHisto(trigcount.c_str(),TRIG_COUNTER = new TH1F(trigcount.c_str(),"pi0 counter",7,-0.5,6.5));
  }
  
  
  if(species>0){
    string centralityname = centname+"_CENTRALITY";
    se->registerHisto(centralityname.c_str(),CENTRALITY = new TH1F(centralityname.c_str(),"centrality",101,-0.5,100.5));
    string rplane = centname+"_RPLANE";
    se->registerHisto(rplane.c_str(),RPLANE = new TH1F(rplane.c_str(),"reaction plane",361,-180.5,180.5));  
  }
  
  //originally only for pp --this is why we change above argument to vertex
  string zvertexname = centname+"_ZVERTEX";
  se->registerHisto(zvertexname.c_str(),ZVERTEX = new TH1F(zvertexname.c_str(),"zvertex dist.",60,-30,30));
  
  
  
  string alldphi = centname+"_PT1PT2DPHI";
  se->registerHisto(alldphi.c_str(),PT1PT2DPHI = new TH3F(alldphi.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
  
  string onedphi = centname+"_PT1PT2DPHI1";
  se->registerHisto(onedphi.c_str(),PT1PT2DPHI1 = new TH3F(onedphi.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));

  if(removetags)
    {
      string twodphi = centname+"_PT1PT2DPHI2";
      se->registerHisto(twodphi.c_str(),PT1PT2DPHI2 = new TH3F(twodphi.c_str(),"trig pt vs. part pt vs dphi",200,0,20,100,0,10,60,-pi/2.,3*pi/2.));
      
      string threedphi = centname+"_PT1PT2DPHI3";
      se->registerHisto(threedphi.c_str(),PT1PT2DPHI3 = new TH3F(threedphi.c_str(),"trig pt vs. part pt vs dphi",200,0,20,100,0,10,60,-pi/2.,3*pi/2.));
      string fourdphi = centname+"_PT1PT2DPHI4";
      se->registerHisto(fourdphi.c_str(),PT1PT2DPHI4 = new TH3F(fourdphi.c_str(),"trig pt vs. part pt vs dphi",200,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    }
  
  if (!_useLessHistos)
    {
      string dphiplus = centname+"_PT1PT2DPHIPLUS";
      se->registerHisto(dphiplus.c_str(),PT1PT2DPHIPLUS = new TH3F(dphiplus.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
      string dphiminus = centname+"_PT1PT2DPHIMINUS";
      se->registerHisto(dphiminus.c_str(),PT1PT2DPHIMINUS = new TH3F(dphiminus.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
      
    }
  
  string fold = centname+"_PT1PT2DPHIFOLD";
  se->registerHisto(fold.c_str(),PT1PT2DPHIFOLD = new TH3F(fold.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
  
  string foldgl = centname+"_PT1PT2DPHIFOLDGL";
  se->registerHisto(foldgl.c_str(),PT1PT2DPHIFOLDGL = new TH3F(foldgl.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
  string foldsc = centname+"_PT1PT2DPHIFOLDSC";
  se->registerHisto(foldsc.c_str(),PT1PT2DPHIFOLDSC = new TH3F(foldsc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
  
  if (!_useLessHistos)
    {
      string pfold = centname+"_PT1PT2DPHIFOLDPLUS";
      se->registerHisto(pfold.c_str(),PT1PT2DPHIFOLDPLUS = new TH3F(pfold.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));    
      string mfold = centname+"_PT1PT2DPHIFOLDMINUS";
      se->registerHisto(mfold.c_str(),PT1PT2DPHIFOLDMINUS = new TH3F(mfold.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    }
  
  string trigpt = centname+"_TRIGPT";
  se->registerHisto(trigpt.c_str(),TRIGPT = new TH1F(trigpt.c_str(),"trig pt ",200,0,20));
  string trigptsc = centname+"_TRIGPTSC";
  se->registerHisto(trigptsc.c_str(),TRIGPTSC = new TH1F(trigptsc.c_str(),"trig pt ",200,0,20));
  string trigptgl = centname+"_TRIGPTGL";
  se->registerHisto(trigptgl.c_str(),TRIGPTGL = new TH1F(trigptgl.c_str(),"trig pt ",200,0,20));

  //add by hge********************************************************************************
  string ptqual = centname+"_PT_QUAL";
  se->registerHisto(ptqual.c_str(), PT_QUAL = new TH2F (ptqual.c_str(),"pt vs quality",35,0,7,3,0,3));
  string ptn0 = centname+"_PT_N0";
  se->registerHisto(ptn0.c_str(), PT_N0 = new TH2F (ptn0.c_str(),"pt vs n0",35,0,7,10,0,10));
  string ptpc3only = centname+"_PT_PC3SDPHISDZ";
  se->registerHisto(ptpc3only.c_str(), PT_PC3SDPHISDZ = new TH1F(ptpc3only.c_str(),"hadron yield after pc3 cut",35,0,7));
  string ptaftn0 = centname+"_PT_AFTN0";
  se->registerHisto(ptaftn0.c_str(), PT_AFTN0 = new TH1F (ptaftn0.c_str(),"hadron yield after n0 cut",35,0,7));
  string ptaftn0qual = centname+"_PT_AFTN0QUAL";
  se->registerHisto(ptaftn0qual.c_str(), PT_AFTN0QUAL = new TH1F (ptaftn0qual.c_str(),"hadron yield after n0+quality cut",35,0,7));
  string ptaftcuts = centname+"_PTAFTCUTS";
  se->registerHisto(ptaftcuts.c_str(), PTAFTCUTS = new TH1F (ptaftcuts.c_str(),"hadron yield after all cuts",35,0,7));
  string dchn0 = centname+"_DCH_N0";
  se->registerHisto(dchn0.c_str(), DCH_N0 = new TH3F (dchn0.c_str(),"dch acceptance after n0 cut",9,0.5,5.,500,-pi/2.,3*pi/2.,80,-80,80));
  string dchn0qual = centname+"_DCH_N0_QUAL";
  se->registerHisto(dchn0qual.c_str(), DCH_N0_QUAL = new TH3F (dchn0qual.c_str(),"dch acceptance after n0+quality cuts",9,0.5,5.,500,-pi/2.,3*pi/2.,80,-80,80));
  string dchn0qualpc3 = centname+"_DCH_N0_QUAL_PC3";
  se->registerHisto(dchn0qualpc3.c_str(), DCH_N0_QUAL_PC3 = new TH3F (dchn0qualpc3.c_str(),"dch acceptance after all cuts",9,0.5,5.,500,-pi/2.,3*pi/2.,80,-80,80));

  string pc3sdphi_undef = centname+"_PT_PC3SDPHI_UNDEF";
  se->registerHisto(pc3sdphi_undef.c_str(), PT_PC3SDPHI_UNDEFINE = new TH1F (pc3sdphi_undef.c_str(),"tracks with undefined pc3sdphi",35,0,7));
  string pc3sdz_undef = centname+"_PT_PC3SDZ_UNDEF";
  se->registerHisto(pc3sdz_undef.c_str(), PT_PC3SDZ_UNDEFINE = new TH1F (pc3sdz_undef.c_str(),"tracks with undefined pc3sdz",35,0,7));
  string pc3dphi_undef = centname+"_PT_PC3DPHI_UNDEF";
  se->registerHisto(pc3dphi_undef.c_str(), PT_PC3DPHI_UNDEFINE = new TH1F (pc3dphi_undef.c_str(),"tracks with undefined pc3dphi",35,0,7));
  string pc3dz_undef = centname+"_PT_PC3DZ_UNDEF";
  se->registerHisto(pc3dz_undef.c_str(), PT_PC3DZ_UNDEFINE = new TH1F (pc3dz_undef.c_str(),"tracks with undefined pc3dz",35,0,7));

  string hadroncheck = centname+"_HADRONCUTCHECK";		
  se->registerHisto(hadroncheck.c_str(), HADRONCUTCHECK = new THmulf(hadroncheck.c_str(),"hadron cuts check"));		
  HADRONCUTCHECK->AddAxis("pt", "p_{T}", 10, 0.5, 7.0);		
  HADRONCUTCHECK->AddAxis("qual", "qual", 2, 0, 2);		
  HADRONCUTCHECK->AddAxis("n0", "n0", 2, 0, 2);		
  HADRONCUTCHECK->AddAxis("pc3sdphi", "pc3sdphi", 40, -10., 10.);		
  HADRONCUTCHECK->AddAxis("pc3sdz", "pc3sdz", 40, -10., 10.);		
  // HADRONCUTCHECK->AddAxis("pc3sdphi", "pc3sdphi", 2, 0, 2);		
  // HADRONCUTCHECK->AddAxis("pc3sdz", "pc3sdz", 2, 0, 2);	
  HADRONCUTCHECK->Sumw2();
  //**********************************************************************************************
  string ptone = centname+"_TRIGPT1";
  se->registerHisto(ptone.c_str(),TRIGPT1 = new TH1F(ptone.c_str(),"trig pt ",200,0,20));

  if(removetags)
    {      
      string pttwo = centname+"_TRIGPT2";
      se->registerHisto(pttwo.c_str(),TRIGPT2 = new TH1F(pttwo.c_str(),"trig pt ",200,0,20));
      string ptthree = centname+"_TRIGPT3";
      se->registerHisto(ptthree.c_str(),TRIGPT3 = new TH1F(ptthree.c_str(),"trig pt ",200,0,20));
      string ptfour = centname+"_TRIGPT4";
      se->registerHisto(ptfour.c_str(),TRIGPT4 = new TH1F(ptfour.c_str(),"trig pt ",200,0,20));
    }  

    string fpbaremc = centname+"_PT1PT2DPHIFOLD_BAR_P_EMC";
    se->registerHisto(fpbaremc.c_str(),PT1PT2DPHIFOLD_BAR_P_EMC = new TH3F(fpbaremc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
    string fmbaremc = centname+"_PT1PT2DPHIFOLD_BAR_M_EMC";
    se->registerHisto(fmbaremc.c_str(),PT1PT2DPHIFOLD_BAR_M_EMC = new TH3F(fmbaremc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));


  if(species<1){
    string pbartof = centname+"_PT1PT2DPHI_BAR_P_TOF";
    se->registerHisto(pbartof.c_str(),PT1PT2DPHI_BAR_P_TOF = new TH3F(pbartof.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    string pbaremc = centname+"_PT1PT2DPHI_BAR_P_EMC";
    se->registerHisto(pbaremc.c_str(),PT1PT2DPHI_BAR_P_EMC = new TH3F(pbaremc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    string mbartof = centname+"_PT1PT2DPHI_BAR_M_TOF";
    se->registerHisto(mbartof.c_str(),PT1PT2DPHI_BAR_M_TOF = new TH3F(mbartof.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    string mbaremc = centname+"_PT1PT2DPHI_BAR_M_EMC";
    se->registerHisto(mbaremc.c_str(),PT1PT2DPHI_BAR_M_EMC = new TH3F(mbaremc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    string pmestof = centname+"_PT1PT2DPHI_MES_P_TOF";
    se->registerHisto(pmestof.c_str(),PT1PT2DPHI_MES_P_TOF = new TH3F(pmestof.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    string pmesemc = centname+"_PT1PT2DPHI_MES_P_EMC";
    se->registerHisto(pmesemc.c_str(),PT1PT2DPHI_MES_P_EMC = new TH3F(pmesemc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    string mmestof = centname+"_PT1PT2DPHI_MES_M_TOF";
    se->registerHisto(mmestof.c_str(),PT1PT2DPHI_MES_M_TOF = new TH3F(mmestof.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    string mmesemc = centname+"_PT1PT2DPHI_MES_M_EMC";
    se->registerHisto(mmesemc.c_str(),PT1PT2DPHI_MES_M_EMC = new TH3F(mmesemc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,60,-pi/2.,3*pi/2.));
    
    
    string fpbartof = centname+"_PT1PT2DPHIFOLD_BAR_P_TOF";
    se->registerHisto(fpbartof.c_str(),PT1PT2DPHIFOLD_BAR_P_TOF = new TH3F(fpbartof.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
    string fmbartof = centname+"_PT1PT2DPHIFOLD_BAR_M_TOF";
    se->registerHisto(fmbartof.c_str(),PT1PT2DPHIFOLD_BAR_M_TOF = new TH3F(fmbartof.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
    string fpmestof = centname+"_PT1PT2DPHIFOLD_MES_P_TOF";
    se->registerHisto(fpmestof.c_str(),PT1PT2DPHIFOLD_MES_P_TOF = new TH3F(fpmestof.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
    string fpmesemc = centname+"_PT1PT2DPHIFOLD_MES_P_EMC";
    se->registerHisto(fpmesemc.c_str(),PT1PT2DPHIFOLD_MES_P_EMC = new TH3F(fpmesemc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
    string fmmestof = centname+"_PT1PT2DPHIFOLD_MES_M_TOF";
    se->registerHisto(fmmestof.c_str(),PT1PT2DPHIFOLD_MES_M_TOF = new TH3F(fmmestof.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
    string fmmesemc = centname+"_PT1PT2DPHIFOLD_MES_M_EMC";
    se->registerHisto(fmmesemc.c_str(),PT1PT2DPHIFOLD_MES_M_EMC = new TH3F(fmmesemc.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
  }


  ///HISTOS TO TEST POSSIBLE ISO CUTS --meg 093009
  if(useiso){
    string emapo = centname+"_EvsR1";
    se->registerHisto(emapo.c_str(),EvsR1 = new TH1F(emapo.c_str(),"radius vs percent E",60,0.0,2*pi));
    string emapt = centname+"_EvsR2";
    se->registerHisto(emapt.c_str(),EvsR2 = new TH1F(emapt.c_str(),"radius vs percent E",60,0.0,2*pi));
    string emapth = centname+"_EvsR3";
    se->registerHisto(emapth.c_str(),EvsR3 = new TH1F(emapth.c_str(),"radius vs percent E",60,0.0,2*pi));
    string emapf = centname+"_EvsR4";
    se->registerHisto(emapf.c_str(),EvsR4 = new TH1F(emapf.c_str(),"radius vs percent E",60,0.0,2*pi));
    
    
    string oneiso = centname+"_ISOMAP1";
    se->registerHisto(oneiso.c_str(),ISOMAP1 = new TH2F(oneiso.c_str(),"radius vs percent E",20,0.05,1.05,500,0,500));
    
    string twoiso = centname+"_ISOMAP2";
    se->registerHisto(twoiso.c_str(),ISOMAP2 = new TH2F(twoiso.c_str(),"radius vs percent E",20,0.05,1.05,500,0,500));
    string threeiso = centname+"_ISOMAP3";
    se->registerHisto(threeiso.c_str(),ISOMAP3 = new TH2F(threeiso.c_str(),"radius vs percent E",20,0.05,1.05,500,0,500));
    string fouriso = centname+"_ISOMAP4";
    se->registerHisto(fouriso.c_str(),ISOMAP4 = new TH2F(fouriso.c_str(),"radius vs percent E",20,0.05,1.05,500,0,500));
    
    
    string oneisosh = centname+"_ISOMAPSH1";
    se->registerHisto(oneisosh.c_str(),ISOMAPSH1 = new TH2F(oneisosh.c_str(),"radius vs percent E",20,0.05,1.05,500,0,500)); 
    string twoisosh = centname+"_ISOMAPSH2";
    se->registerHisto(twoisosh.c_str(),ISOMAPSH2 = new TH2F(twoisosh.c_str(),"radius vs percent E",20,0.05,1.05,500,0,500));
    string threeisosh = centname+"_ISOMAPSH3";
    se->registerHisto(threeisosh.c_str(),ISOMAPSH3 = new TH2F(threeisosh.c_str(),"radius vs percent E",20,0.05,1.05,500,0,500));
    string fourisosh = centname+"_ISOMAPSH4";
    se->registerHisto(fourisosh.c_str(),ISOMAPSH4 = new TH2F(fourisosh.c_str(),"radius vs percent E",20,0.05,1.05,500,0,500));
  } 
  

  PT1PT2DPHI->Sumw2();
  if(removetags)
    {
      PT1PT2DPHI1->Sumw2();
      PT1PT2DPHI2->Sumw2();
      PT1PT2DPHI3->Sumw2();
      PT1PT2DPHI4->Sumw2();
    }
  PT1PT2DPHIFOLD->Sumw2();
  PT1PT2DPHIFOLDGL->Sumw2();
  PT1PT2DPHIFOLDSC->Sumw2();
  
  PT1PT2DPHIFOLD_BAR_P_EMC->Sumw2();
  PT1PT2DPHIFOLD_BAR_M_EMC->Sumw2();
  
  if(species<1){
    PT1PT2DPHI_BAR_P_TOF->Sumw2();
    PT1PT2DPHI_BAR_P_EMC->Sumw2();
    PT1PT2DPHI_BAR_M_TOF->Sumw2();
    PT1PT2DPHI_BAR_M_EMC->Sumw2();
    PT1PT2DPHI_MES_P_TOF->Sumw2();
    PT1PT2DPHI_MES_P_EMC->Sumw2();
    PT1PT2DPHI_MES_M_TOF->Sumw2();
    PT1PT2DPHI_MES_M_EMC->Sumw2();
    
    PT1PT2DPHIFOLD_BAR_P_TOF->Sumw2();
    PT1PT2DPHIFOLD_BAR_M_TOF->Sumw2();
    PT1PT2DPHIFOLD_MES_P_TOF->Sumw2();
    PT1PT2DPHIFOLD_MES_P_EMC->Sumw2();
    PT1PT2DPHIFOLD_MES_M_TOF->Sumw2();
    PT1PT2DPHIFOLD_MES_M_EMC->Sumw2();
  }

  if(useiso){
    EvsR1->Sumw2();
    EvsR2->Sumw2();
    EvsR3->Sumw2();
    EvsR4->Sumw2();
    
    ISOMAP1->Sumw2();
    ISOMAP2->Sumw2();
    ISOMAP3->Sumw2();
    ISOMAP4->Sumw2();
    
    ISOMAPSH1->Sumw2();
    ISOMAPSH2->Sumw2();
    ISOMAPSH3->Sumw2();
    ISOMAPSH4->Sumw2();
  }

  string teffw = centname+"_TRIGTEFFW";
  se->registerHisto(teffw.c_str(),TRIGPTEFFW = new TH1F(teffw.c_str(),"trig pt ",200,0,20));
  
  string vetocount = centname+"_VETO_COUNTER";
  se->registerHisto(vetocount.c_str(),VETO_COUNTER = new TH1F(vetocount.c_str(),"veto counter ",200,0.0,20.));

  string vetocount0 = centname+"_VETO_COUNTER_0";
  se->registerHisto(vetocount0.c_str(),VETO_COUNTER_0 = new TH1F(vetocount0.c_str(),"veto counter test",200,0.0,20.));


  string trignorm = centname+"_ptvscent_trig";
  se->registerHisto(trignorm.c_str(),ptvscent_trig=new TH2F(trignorm.c_str(),"ptT vs cent ",200,0.0,20.0,101,0,100));

  string decnorm = centname+"_ptvscent_dec";
  se->registerHisto(decnorm.c_str(),ptvscent_dec=new TH2F(decnorm.c_str(),"ptT vs cent ",200,0.0,20.0,101,0,100));

  string partnorm = centname+"_ptvscent_part";
  se->registerHisto(partnorm.c_str(),ptvscent_part=new TH2F(partnorm.c_str(),"ptp vs cent ",100,0.0,10.0,101,0,100));


  string mtof = centname+"_M2TOF";
  se->registerHisto(mtof.c_str(),M2TOF=new TH2F(mtof.c_str(),"m2 tof ",200,-0.5,1.5,100,0,10));
  string memc = centname+"_M2EMC";
  se->registerHisto(memc.c_str(),M2EMC=new TH2F(memc.c_str(),"m2 emc ",200,-0.5,1.5,100,0,10));

  if( (species==5 || species==1) )
  {
    string pcsdphi = centname+"_PT_PC3SDPHI_PC3SDZ";
    se->registerHisto(pcsdphi.c_str(),PT_PC3SDPHI_PC3SDZ=new TH3F(pcsdphi.c_str(),"PC3SDPHI_PC3SDZ",20,0,10,100,-10,10,100,-10,10));
    PT_PC3SDPHI_PC3SDZ->Sumw2();
  }

  //Include invmass histo again MEC 9/16/08 
  //for tagging 01/07/09 and reduced bins
    //if(!skip_fg_histos&&tagflag>0){
  if((!skip_fg_histos&&tagflag>0 ) || removetags){
    string masszero = centname+"_INVMASS_0";
    se->registerHisto(masszero.c_str(),INVMASS[0] = new TH2F(masszero.c_str(),"INVMASS_0",200,0,20,200,0,1));
    string massone = centname+"_INVMASS_1";
    se->registerHisto(massone.c_str(),INVMASS[1] = new TH2F(massone.c_str(),"INVMASS_1",200,0,20,200,0,1));
    string masstwo = centname+"_INVMASS_2";
    se->registerHisto(masstwo.c_str(),INVMASS[2] = new TH2F(masstwo.c_str(),"INVMASS_2",200,0,20,200,0,1));
    string massthree = centname+"_INVMASS_3";
    se->registerHisto(massthree.c_str(),INVMASS[3] = new TH2F(massthree.c_str(),"INVMASS_3",200,0,20,200,0,1));
    string massfour = centname+"_INVMASS_4";
    se->registerHisto(massfour.c_str(),INVMASS[4] = new TH2F(massfour.c_str(),"INVMASS_4",200,0,20,200,0,1));
    string massfive = centname+"_INVMASS_5";
    se->registerHisto(massfive.c_str(),INVMASS[5] = new TH2F(massfive.c_str(),"INVMASS_5",200,0,20,200,0,1));
    string masssix = centname+"_INVMASS_6";
    se->registerHisto(masssix.c_str(),INVMASS[6] = new TH2F(masssix.c_str(),"INVMASS_6",200,0,20,200,0,1));
    string masssev = centname+"_INVMASS_7";
    se->registerHisto(masssev.c_str(),INVMASS[7] = new TH2F(masssev.c_str(),"INVMASS_7",200,0,20,200,0,1));
    
    string decmassname = centname+"_DECINVMASS";
    se->registerHisto(decmassname.c_str(),DECINVMASS = new TH2F(decmassname.c_str(),"DECINVMASS",7,0,7,200,0,1));

    }


  if(dofilltime){
    string ztmap = centname+"_PT1PT2ZT";
    se->registerHisto(ztmap.c_str(),PT1PT2ZT = new TH3F(ztmap.c_str(),"trig pt vs. part pt vs zt",20,0,20,100,0,10,40,0.0,2.0));
    
    string ximap = centname+"_PT1PT2XI";
    se->registerHisto(ximap.c_str(),PT1PT2XI = new TH3F(ximap.c_str(),"trig pt vs. part pt vs xi",20,0,20,100,0,10,25,-1.0,4.0));
    
    string wfold = centname+"_PT1PT2DPHIFOLDCF";
    se->registerHisto(wfold.c_str(),PT1PT2DPHIFOLDCF = new TH3F(wfold.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
    
    string wflow = centname+"_PT1PT2DPHIFLOW";
    se->registerHisto(wflow.c_str(),PT1PT2DPHIFLOW = new TH3F(wflow.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
    
    
    PT1PT2ZT->Sumw2();
    PT1PT2XI->Sumw2();
    
    
    PT1PT2DPHIFOLDCF->Sumw2();
    PT1PT2DPHIFLOW->Sumw2();
    
    
      
      vector <string> str_side;
    str_side.push_back("NEAR");
    str_side.push_back("AWAY");
    str_side.push_back("NEARBG");
    str_side.push_back("AWAYBG");
    
    cout<<" vectors "<<endl;
    
    //the binning for these histograms is {5,6,7,8,9,10,11,12} however other binnings may be projected out of the 3D histogram
    vector <string> str_pt;
    for(int i=0;i<7;i++){ //used to be 8?
      ostringstream ptnamer;
      ptnamer<<i+5;
      str_pt.push_back(ptnamer.str());    
    }
    
    
    cout<< " book more histos "<<endl;
    
    for(int l=0;l<4;l++){
      
      if(tagflag>0){
	string decztmap = centname+"_DECPT2ZT_"+str_pt[l];
	se->registerHisto(decztmap.c_str(),DECPT2ZT[l] = new TH2F(decztmap.c_str(),"part pt vs zt",100,0,10,40,0.0,2.0));
	
	string decximap = centname+"_DECPT2XI_"+str_pt[l];
	se->registerHisto(decximap.c_str(),DECPT2XI[l] = new TH2F(decximap.c_str(),"part pt vs xi",100,0,10,25,-1.0,4.0));
	
	/*
	  string decwfold = centname+"_DECPT1PT2DPHIFOLDCF";
	  se->registerHisto(decwfold.c_str(),DECPT1PT2DPHIFOLDCF = new TH3F(decwfold.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
	  
	  string decwflow = centname+"_DECPT1PT2DPHIFLOW";
	  se->registerHisto(decwflow.c_str(),DECPT1PT2DPHIFLOW = new TH3F(decwflow.c_str(),"trig pt vs. part pt vs dphi",20,0,20,100,0,10,30,0,pi));
	*/
	
	DECPT2ZT[l]->Sumw2();
	DECPT2XI[l]->Sumw2();
	
	
	

	//DECPT1PT2DPHIFOLDCF->Sumw2();
	//DECPT1PT2DPHIFLOW->Sumw2();
	
      }
    }



    
    //i==0-> trigger above 5 GeV, i==3 -> above 8 GeV 
    //j==0 -> west, j==1-> east
    
    

    for(int j=0;j<(int)str_side.size();j++){
      
      
      //energy threshold
      for(int l=0;l<7;l++){
	
	
	//error histograms not necessary since bg error negligible
	
	string xedistname=centname+"XEDIST_"+str_side[j]+"_"+str_pt[l];
	string decxedistname=centname+"DECXEDIST_"+str_side[j]+"_"+str_pt[l];

	string xidistname=centname+"XIDIST_"+str_side[j]+"_"+str_pt[l];
	string decxidistname=centname+"DECXIDIST_"+str_side[j]+"_"+str_pt[l];
	
	string ztdistname=centname+"ZTDIST_"+str_side[j]+"_"+str_pt[l];
	string decztdistname=centname+"DECZTDIST_"+str_side[j]+"_"+str_pt[l];


	string poutdistname=centname+"POUTDIST_"+str_side[j]+"_"+str_pt[l];
	string decpoutdistname=centname+"DECPOUTDIST_"+str_side[j]+"_"+str_pt[l];

	string pout2distname=centname+"POUT2DIST_"+str_side[j]+"_"+str_pt[l];
	string decpout2distname=centname+"DECPOUT2DIST_"+str_side[j]+"_"+str_pt[l];

	string pout2sumname=centname+"POUT2SUM_"+str_side[j]+"_"+str_pt[l];
	string decpout2sumname=centname+"DECPOUT2SUM_"+str_side[j]+"_"+str_pt[l];

	string xevspout2name=centname+"XEVSPOUT2_"+str_side[j]+"_"+str_pt[l];
	string decxevspout2name=centname+"DECXEVSPOUT2_"+str_side[j]+"_"+str_pt[l];
	string ztvspout2name=centname+"ZTVSPOUT2_"+str_side[j]+"_"+str_pt[l];
	string decztvspout2name=centname+"DECZTVSPOUT2_"+str_side[j]+"_"+str_pt[l];

	
     
	se->registerHisto(xedistname.c_str(),XEDIST[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),40,0.,2.,30,0,pi));      
	//se->registerHisto(xedistname.c_str(),XEDIST[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),25,0.,2.,30,0,pi));      
	se->registerHisto(xidistname.c_str(),XIDIST[j][l]=new TH2F(xidistname.c_str(),xidistname.c_str(),25,-1.0,4.,30,0,pi));      	
	se->registerHisto(ztdistname.c_str(),ZTDIST[j][l]=new TH2F(ztdistname.c_str(),ztdistname.c_str(),40,0.,2.,30,0,pi));      


	XEDIST[j][l]->Sumw2();
	XIDIST[j][l]->Sumw2();
	ZTDIST[j][l]->Sumw2();



	if(pouthistos){

	  //se->registerHisto(poutdistname.c_str(),POUTDIST[j][l]=new TH2F(poutdistname.c_str(),poutdistname.c_str(),200,0.,5.,30,0,pi));      
	  se->registerHisto(poutdistname.c_str(),POUTDIST[j][l]=new TH2F(poutdistname.c_str(),poutdistname.c_str(),20,0.,10.,30,0,pi));      
	  se->registerHisto(pout2distname.c_str(),POUT2DIST[j][l]=new TH2F(pout2distname.c_str(),pout2distname.c_str(),1000,0.,100.,30,0,pi));      
	  se->registerHisto(pout2sumname.c_str(),POUT2SUM[j][l]=new TH1F(pout2sumname.c_str(),pout2sumname.c_str(),2,0.,2.));      
	  //se->registerHisto(xevspout2name.c_str(),XEVSPOUT2[j][l]=new TH2F(xevspout2name.c_str(),xevspout2name.c_str(),40,0.,2.,1000,0.,10.));      
	  //se->registerHisto(ztvspout2name.c_str(),ZTVSPOUT2[j][l]=new TH2F(ztvspout2name.c_str(),ztvspout2name.c_str(),40,0.,2.,1000,0.,10.));      
	  
	  POUTDIST[j][l]->Sumw2();
	  POUT2DIST[j][l]->Sumw2();
	  POUT2SUM[j][l]->Sumw2();
	  //XEVSPOUT2[j][l]->Sumw2();
	  //ZTVSPOUT2[j][l]->Sumw2();
	}	


	if(tagflag>0){ 
	  se->registerHisto(decxedistname.c_str(),DECXEDIST[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),40,0.,2.,30,0,pi));      	  
	  //se->registerHisto(decxedistname.c_str(),DECXEDIST[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),25,0.,2.,30,0,pi));      	  
	  //se->registerHisto(decztdistname.c_str(),DECZTDIST[j][l]=new TH2F(decztdistname.c_str(),decztdistname.c_str(),40,0.,2.,30,0,pi));      
	  se->registerHisto(decztdistname.c_str(),DECZTDIST[j][l]=new TH2F(decztdistname.c_str(),decztdistname.c_str(),25,0.,2.,30,0,pi));      
	  se->registerHisto(decxidistname.c_str(),DECXIDIST[j][l]=new TH2F(decxidistname.c_str(),decxidistname.c_str(),25,-1.0,4.,30,0,pi));      

	  DECXEDIST[j][l]->Sumw2();
	  DECZTDIST[j][l]->Sumw2();
	  DECXIDIST[j][l]->Sumw2();

	  if(pouthistos_more){
	  //se->registerHisto(decpoutdistname.c_str(),DECPOUTDIST[j][l]=new TH2F(decpoutdistname.c_str(),decpoutdistname.c_str(),200,0.,5.,30,0,pi));      
	  se->registerHisto(decpoutdistname.c_str(),DECPOUTDIST[j][l]=new TH2F(decpoutdistname.c_str(),decpoutdistname.c_str(),20,0.,10.,30,0,pi));      
	  se->registerHisto(decpout2distname.c_str(),DECPOUT2DIST[j][l]=new TH2F(decpout2distname.c_str(),decpout2distname.c_str(),1000,0.,100.,30,0,pi));      
	  se->registerHisto(decpout2sumname.c_str(),DECPOUT2SUM[j][l]=new TH1F(decpout2sumname.c_str(),decpout2sumname.c_str(),2,0.,2.));
	  //se->registerHisto(decxevspout2name.c_str(),DECXEVSPOUT2[j][l]=new TH2F(decxevspout2name.c_str(),decxevspout2name.c_str(),40,0.,2.,1000,0.,10.));      	  
	  //se->registerHisto(decztvspout2name.c_str(),DECZTVSPOUT2[j][l]=new TH2F(decztvspout2name.c_str(),decztvspout2name.c_str(),40,0.,2.,1000,0.,10.));      

	    DECPOUTDIST[j][l]->Sumw2();
	    DECPOUT2DIST[j][l]->Sumw2();
	    DECPOUT2SUM[j][l]->Sumw2();
	    //DECXEVSPOUT2[j][l]->Sumw2();
	    //DECZTVSPOUT2[j][l]->Sumw2();
	  }
	}				
      }
      
      
          
    
      
      //energy threshold
      for(int l=0;l<4;l++){
	                            	
	
	string xedistname=centname+"XEDIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string decxedistname=centname+"DECXEDIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string xidistname=centname+"XIDIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string decxidistname=centname+"DECXIDIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string ztdistname=centname+"ZTDIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string decztdistname=centname+"DECZTDIST_LARGE_"+str_side[j]+"_"+str_pt[l];

	
	string nwxedistname=centname+"XEDIST_LARGE_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwdecxedistname=centname+"DECXEDIST_LARGE_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwxidistname=centname+"XIDIST_LARGE_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwdecxidistname=centname+"DECXIDIST_LARGE_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwztdistname=centname+"ZTDIST_LARGE_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwdecztdistname=centname+"DECZTDIST_LARGE_NOW_"+str_side[j]+"_"+str_pt[l];
	

	
	string poutdistname=centname+"POUTDIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string decpoutdistname=centname+"DECPOUTDIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string pout2distname=centname+"POUT2DIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string decpout2distname=centname+"DECPOUT2DIST_LARGE_"+str_side[j]+"_"+str_pt[l];
	string pout2sumname=centname+"POUT2SUM_LARGE_"+str_side[j]+"_"+str_pt[l];
	string decpout2sumname=centname+"DECPOUT2SUM_LARGE_"+str_side[j]+"_"+str_pt[l];
	
	
	se->registerHisto(xedistname.c_str(),XEDIST_LARGE[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),40,0.,2.,30,0,pi));      
	//se->registerHisto(xedistname.c_str(),XEDIST_LARGE[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),25,0.,2.,30,0,pi));      
	se->registerHisto(xidistname.c_str(),XIDIST_LARGE[j][l]=new TH2F(xidistname.c_str(),xidistname.c_str(),25,-1.0,4.,30,0,pi));      
	//se->registerHisto(ztdistname.c_str(),ZTDIST_LARGE[j][l]=new TH2F(ztdistname.c_str(),ztdistname.c_str(),40,0.,2.,30,0,pi));      
	se->registerHisto(ztdistname.c_str(),ZTDIST_LARGE[j][l]=new TH2F(ztdistname.c_str(),ztdistname.c_str(),25,0.,2.,30,0,pi));      
	if(pouthistos){
	  //se->registerHisto(poutdistname.c_str(),POUTDIST_LARGE[j][l]=new TH2F(poutdistname.c_str(),poutdistname.c_str(),200,0.,5.,30,0,pi));      
	  se->registerHisto(poutdistname.c_str(),POUTDIST_LARGE[j][l]=new TH2F(poutdistname.c_str(),poutdistname.c_str(),20,0.,10.,30,0,pi));      
	  se->registerHisto(pout2distname.c_str(),POUT2DIST_LARGE[j][l]=new TH2F(pout2distname.c_str(),pout2distname.c_str(),1000,0.,100.,30,0,pi));      
	  se->registerHisto(pout2sumname.c_str(),POUT2SUM_LARGE[j][l]=new TH1F(pout2sumname.c_str(),pout2sumname.c_str(),2,0.,2.));
	}


	XEDIST_LARGE[j][l]->Sumw2();	
	XIDIST_LARGE[j][l]->Sumw2();	
	ZTDIST_LARGE[j][l]->Sumw2();	


	if(j<2){
	  se->registerHisto(nwxedistname.c_str(),XEDIST_LARGE_NOW[j][l]=new TH2F(nwxedistname.c_str(),nwxedistname.c_str(),40,0.,2.,30,0,pi));      
	  //se->registerHisto(xedistname.c_str(),XEDIST[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),25,0.,2.,30,0,pi));      
	  se->registerHisto(nwxidistname.c_str(),XIDIST_LARGE_NOW[j][l]=new TH2F(nwxidistname.c_str(),nwxidistname.c_str(),25,-1.0,4.,30,0,pi));      	
	  se->registerHisto(nwztdistname.c_str(),ZTDIST_LARGE_NOW[j][l]=new TH2F(nwztdistname.c_str(),nwztdistname.c_str(),40,0.,2.,30,0,pi));      
	  
	  
	  XEDIST_LARGE_NOW[j][l]->Sumw2();
	  XIDIST_LARGE_NOW[j][l]->Sumw2();
	  ZTDIST_LARGE_NOW[j][l]->Sumw2();
	}


	if(pouthistos){
	  POUTDIST_LARGE[j][l]->Sumw2();
	  POUT2DIST_LARGE[j][l]->Sumw2();
	  POUT2SUM_LARGE[j][l]->Sumw2();
	}

	if(tagflag>0){
	  se->registerHisto(decxedistname.c_str(),DECXEDIST_LARGE[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),40,0.,2.,30,0,pi));      	  
	  //se->registerHisto(decxedistname.c_str(),DECXEDIST_LARGE[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),25,0.,2.,30,0,pi));      	  
	  se->registerHisto(decxidistname.c_str(),DECXIDIST_LARGE[j][l]=new TH2F(decxidistname.c_str(),decxidistname.c_str(),25,-1.0,4.,30,0,pi));      	  
	  //se->registerHisto(decztdistname.c_str(),DECZTDIST_LARGE[j][l]=new TH2F(decztdistname.c_str(),decztdistname.c_str(),40,0.,2.,30,0,pi));      
	  se->registerHisto(decztdistname.c_str(),DECZTDIST_LARGE[j][l]=new TH2F(decztdistname.c_str(),decztdistname.c_str(),25,0.,2.,30,0,pi));      


	  if(j<2){
	    se->registerHisto(nwdecxedistname.c_str(),DECXEDIST_LARGE_NOW[j][l]=new TH2F(nwdecxedistname.c_str(),nwdecxedistname.c_str(),40,0.,2.,30,0,pi));      	  
	    se->registerHisto(nwdecxidistname.c_str(),DECXIDIST_LARGE_NOW[j][l]=new TH2F(nwdecxidistname.c_str(),nwdecxidistname.c_str(),25,-1.0,4.,30,0,pi));      	  
	    se->registerHisto(nwdecztdistname.c_str(),DECZTDIST_LARGE_NOW[j][l]=new TH2F(nwdecztdistname.c_str(),nwdecztdistname.c_str(),25,0.,2.,30,0,pi));      
	  }



	



	  if(pouthistos){
	    //se->registerHisto(decpoutdistname.c_str(),DECPOUTDIST_LARGE[j][l]=new TH2F(decpoutdistname.c_str(),decpoutdistname.c_str(),200,0.,5.,30,0,pi));    
	    se->registerHisto(decpoutdistname.c_str(),DECPOUTDIST_LARGE[j][l]=new TH2F(decpoutdistname.c_str(),decpoutdistname.c_str(),20,0.,10.,30,0,pi));      
	    se->registerHisto(decpout2distname.c_str(),DECPOUT2DIST_LARGE[j][l]=new TH2F(decpout2distname.c_str(),decpout2distname.c_str(),1000,0.,100.,30,0,pi));
	    se->registerHisto(decpout2sumname.c_str(),DECPOUT2SUM_LARGE[j][l]=new TH1F(decpout2sumname.c_str(),decpout2sumname.c_str(),2,0.,2.));
	  }

	  DECXEDIST_LARGE[j][l]->Sumw2();
	  DECXIDIST_LARGE[j][l]->Sumw2();
	  DECZTDIST_LARGE[j][l]->Sumw2();
	  if(pouthistos){
	    DECPOUTDIST_LARGE[j][l]->Sumw2();
	    DECPOUT2DIST_LARGE[j][l]->Sumw2();
	    DECPOUT2SUM_LARGE[j][l]->Sumw2();
	  }
	}
      	
      }



      //same but for stat errors
      for(int l=0;l<4;l++){
	                            	
	
	string xedistname=centname+"XEDIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string decxedistname=centname+"DECXEDIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string xidistname=centname+"XIDIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string decxidistname=centname+"DECXIDIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string ztdistname=centname+"ZTDIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string decztdistname=centname+"DECZTDIST_STATERR_"+str_side[j]+"_"+str_pt[l];

	
	string nwxedistname=centname+"XEDIST_STATERR_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwdecxedistname=centname+"DECXEDIST_STATERR_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwxidistname=centname+"XIDIST_STATERR_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwdecxidistname=centname+"DECXIDIST_STATERR_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwztdistname=centname+"ZTDIST_STATERR_NOW_"+str_side[j]+"_"+str_pt[l];
	string nwdecztdistname=centname+"DECZTDIST_STATERR_NOW_"+str_side[j]+"_"+str_pt[l];
	

	
	string poutdistname=centname+"POUTDIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string decpoutdistname=centname+"DECPOUTDIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string pout2distname=centname+"POUT2DIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string decpout2distname=centname+"DECPOUT2DIST_STATERR_"+str_side[j]+"_"+str_pt[l];
	string pout2sumname=centname+"POUT2SUM_STATERR_"+str_side[j]+"_"+str_pt[l];
	string decpout2sumname=centname+"DECPOUT2SUM_STATERR_"+str_side[j]+"_"+str_pt[l];
	
	
	se->registerHisto(xedistname.c_str(),XEDIST_STATERR[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),40,0.,2.,30,0,pi));      
	//se->registerHisto(xedistname.c_str(),XEDIST_STATERR[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),25,0.,2.,30,0,pi));      
	se->registerHisto(xidistname.c_str(),XIDIST_STATERR[j][l]=new TH2F(xidistname.c_str(),xidistname.c_str(),25,-1.0,4.,30,0,pi));      
	//se->registerHisto(ztdistname.c_str(),ZTDIST_STATERR[j][l]=new TH2F(ztdistname.c_str(),ztdistname.c_str(),40,0.,2.,30,0,pi));      
	se->registerHisto(ztdistname.c_str(),ZTDIST_STATERR[j][l]=new TH2F(ztdistname.c_str(),ztdistname.c_str(),25,0.,2.,30,0,pi));      
	if(pouthistos){
	  //se->registerHisto(poutdistname.c_str(),POUTDIST_STATERR[j][l]=new TH2F(poutdistname.c_str(),poutdistname.c_str(),200,0.,5.,30,0,pi));      
	  se->registerHisto(poutdistname.c_str(),POUTDIST_STATERR[j][l]=new TH2F(poutdistname.c_str(),poutdistname.c_str(),20,0.,10.,30,0,pi));      
	  se->registerHisto(pout2distname.c_str(),POUT2DIST_STATERR[j][l]=new TH2F(pout2distname.c_str(),pout2distname.c_str(),1000,0.,100.,30,0,pi));      
	  se->registerHisto(pout2sumname.c_str(),POUT2SUM_STATERR[j][l]=new TH1F(pout2sumname.c_str(),pout2sumname.c_str(),2,0.,2.));
	}


	XEDIST_STATERR[j][l]->Sumw2();	
	XIDIST_STATERR[j][l]->Sumw2();	
	ZTDIST_STATERR[j][l]->Sumw2();	

	/*
	if(j<2){
	  se->registerHisto(nwxedistname.c_str(),XEDIST_STATERR_NOW[j][l]=new TH2F(nwxedistname.c_str(),nwxedistname.c_str(),40,0.,2.,30,0,pi));      
	  //se->registerHisto(xedistname.c_str(),XEDIST[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),25,0.,2.,30,0,pi));      
	  se->registerHisto(nwxidistname.c_str(),XIDIST_STATERR_NOW[j][l]=new TH2F(nwxidistname.c_str(),nwxidistname.c_str(),25,-1.0,4.,30,0,pi));      	
	  se->registerHisto(nwztdistname.c_str(),ZTDIST_STATERR_NOW[j][l]=new TH2F(nwztdistname.c_str(),nwztdistname.c_str(),40,0.,2.,30,0,pi));      
	  
	  
	  XEDIST_STATERR_NOW[j][l]->Sumw2();
	  XIDIST_STATERR_NOW[j][l]->Sumw2();
	  ZTDIST_STATERR_NOW[j][l]->Sumw2();
	}
	*/

	if(pouthistos){
	  POUTDIST_STATERR[j][l]->Sumw2();
	  POUT2DIST_STATERR[j][l]->Sumw2();
	  POUT2SUM_STATERR[j][l]->Sumw2();
	}

	if(tagflag>0){
	  se->registerHisto(decxedistname.c_str(),DECXEDIST_STATERR[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),40,0.,2.,30,0,pi));      	  
	  //se->registerHisto(decxedistname.c_str(),DECXEDIST_STATERR[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),25,0.,2.,30,0,pi));      	  
	  se->registerHisto(decxidistname.c_str(),DECXIDIST_STATERR[j][l]=new TH2F(decxidistname.c_str(),decxidistname.c_str(),25,-1.0,4.,30,0,pi));      	  
	  //se->registerHisto(decztdistname.c_str(),DECZTDIST_STATERR[j][l]=new TH2F(decztdistname.c_str(),decztdistname.c_str(),40,0.,2.,30,0,pi));      
	  se->registerHisto(decztdistname.c_str(),DECZTDIST_STATERR[j][l]=new TH2F(decztdistname.c_str(),decztdistname.c_str(),25,0.,2.,30,0,pi));      

	  /*
	  if(j<2){
	    se->registerHisto(nwdecxedistname.c_str(),DECXEDIST_STATERR_NOW[j][l]=new TH2F(nwdecxedistname.c_str(),nwdecxedistname.c_str(),40,0.,2.,30,0,pi));      	  
	    se->registerHisto(nwdecxidistname.c_str(),DECXIDIST_STATERR_NOW[j][l]=new TH2F(nwdecxidistname.c_str(),nwdecxidistname.c_str(),25,-1.0,4.,30,0,pi));      	  
	    se->registerHisto(nwdecztdistname.c_str(),DECZTDIST_STATERR_NOW[j][l]=new TH2F(nwdecztdistname.c_str(),nwdecztdistname.c_str(),25,0.,2.,30,0,pi));      
	  }
	  */


	  if(pouthistos){
	    //se->registerHisto(decpoutdistname.c_str(),DECPOUTDIST_STATERR[j][l]=new TH2F(decpoutdistname.c_str(),decpoutdistname.c_str(),200,0.,5.,30,0,pi));    
	    se->registerHisto(decpoutdistname.c_str(),DECPOUTDIST_STATERR[j][l]=new TH2F(decpoutdistname.c_str(),decpoutdistname.c_str(),20,0.,10.,30,0,pi));      
	    se->registerHisto(decpout2distname.c_str(),DECPOUT2DIST_STATERR[j][l]=new TH2F(decpout2distname.c_str(),decpout2distname.c_str(),1000,0.,100.,30,0,pi));
	    se->registerHisto(decpout2sumname.c_str(),DECPOUT2SUM_STATERR[j][l]=new TH1F(decpout2sumname.c_str(),decpout2sumname.c_str(),2,0.,2.));
	  }

	  DECXEDIST_STATERR[j][l]->Sumw2();
	  DECXIDIST_STATERR[j][l]->Sumw2();
	  DECZTDIST_STATERR[j][l]->Sumw2();
	  if(pouthistos){
	    DECPOUTDIST_STATERR[j][l]->Sumw2();
	    DECPOUT2DIST_STATERR[j][l]->Sumw2();
	    DECPOUT2SUM_STATERR[j][l]->Sumw2();
	  }
	}
      	
      }



      
     
      //same but for plus charge...which is sys up for AuAu 
      for(int l=0;l<4;l++){
	                            	
	
	string xedistname=centname+"XEDIST_LARGE_PLUS_"+str_side[j]+"_"+str_pt[l];
	string xidistname=centname+"XIDIST_LARGE_PLUS_"+str_side[j]+"_"+str_pt[l];
	string decxedistname=centname+"DECXEDIST_LARGE_PLUS_"+str_side[j]+"_"+str_pt[l];
	string decxidistname=centname+"DECXIDIST_LARGE_PLUS_"+str_side[j]+"_"+str_pt[l];
	string ztdistname=centname+"ZTDIST_LARGE_PLUS_"+str_side[j]+"_"+str_pt[l];
	string decztdistname=centname+"DECZTDIST_LARGE_PLUS_"+str_side[j]+"_"+str_pt[l];
	string ztfinedistname=centname+"ZTDIST_PLUS_"+str_side[j]+"_"+str_pt[l];
	string decztfinedistname=centname+"DECZTDIST_PLUS_"+str_side[j]+"_"+str_pt[l];
	string poutdistname=centname+"POUTDIST_LARGE_PLUS_"+str_side[j]+"_"+str_pt[l];
	string decpoutdistname=centname+"DECPOUTDIST_LARGE_PLUS_"+str_side[j]+"_"+str_pt[l];
	

	se->registerHisto(xidistname.c_str(),XIDIST_LARGE_PLUS[j][l]=new TH2F(xidistname.c_str(),xidistname.c_str(),25,-1.0,4.,30,0,pi));      	
	
	se->registerHisto(xedistname.c_str(),XEDIST_LARGE_PLUS[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),40,0.,2.,30,0,pi));      
	//se->registerHisto(xedistname.c_str(),XEDIST_LARGE_PLUS[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),25,0.,2.,30,0,pi));      
	se->registerHisto(ztfinedistname.c_str(),ZTDIST_PLUS[j][l]=new TH2F(ztfinedistname.c_str(),ztfinedistname.c_str(),40,0.,2.,30,0,pi));      
	se->registerHisto(ztdistname.c_str(),ZTDIST_LARGE_PLUS[j][l]=new TH2F(ztdistname.c_str(),ztdistname.c_str(),25,0.,2.,30,0,pi));      
	if(pouthistos) se->registerHisto(poutdistname.c_str(),POUTDIST_LARGE_PLUS[j][l]=new TH2F(poutdistname.c_str(),poutdistname.c_str(),20,0.,10.,30,0,pi));      

	XEDIST_LARGE_PLUS[j][l]->Sumw2();	
	XIDIST_LARGE_PLUS[j][l]->Sumw2();	
	ZTDIST_LARGE_PLUS[j][l]->Sumw2();	
	ZTDIST_PLUS[j][l]->Sumw2();	

	if(pouthistos) POUTDIST_LARGE_PLUS[j][l]->Sumw2();	
	
	if(tagflag>0){
	  se->registerHisto(decxidistname.c_str(),DECXIDIST_LARGE_PLUS[j][l]=new TH2F(decxidistname.c_str(),decxidistname.c_str(),25,-1.0,4.,30,0,pi));      	  
	  se->registerHisto(decxedistname.c_str(),DECXEDIST_LARGE_PLUS[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),40,0.,2.,30,0,pi)); 
     	  //se->registerHisto(decxedistname.c_str(),DECXEDIST_LARGE_PLUS[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),25,0.,2.,30,0,pi));      	  
	  se->registerHisto(decztfinedistname.c_str(),DECZTDIST_PLUS[j][l]=new TH2F(decztfinedistname.c_str(),decztfinedistname.c_str(),40,0.,2.,30,0,pi));      
	  se->registerHisto(decztdistname.c_str(),DECZTDIST_LARGE_PLUS[j][l]=new TH2F(decztdistname.c_str(),decztdistname.c_str(),25,0.,2.,30,0,pi));      
	  if(pouthistos) se->registerHisto(decpoutdistname.c_str(),DECPOUTDIST_LARGE_PLUS[j][l]=new TH2F(decpoutdistname.c_str(),decpoutdistname.c_str(),20,0.,10.,30,0,pi));      
	  DECXEDIST_LARGE_PLUS[j][l]->Sumw2();
	  DECXIDIST_LARGE_PLUS[j][l]->Sumw2();
	  DECZTDIST_LARGE_PLUS[j][l]->Sumw2();
	  DECZTDIST_PLUS[j][l]->Sumw2();
	  if(pouthistos) DECPOUTDIST_LARGE_PLUS[j][l]->Sumw2();
	  
	  
	}
      	
      }


      //same but for minus charge 
      for(int l=0;l<4;l++){
	                            	
	string decxidistname=centname+"DECXIDIST_LARGE_MINUS_"+str_side[j]+"_"+str_pt[l];
	string xidistname=centname+"XIDIST_LARGE_MINUS_"+str_side[j]+"_"+str_pt[l];
	string xedistname=centname+"XEDIST_LARGE_MINUS_"+str_side[j]+"_"+str_pt[l];
	string decxedistname=centname+"DECXEDIST_LARGE_MINUS_"+str_side[j]+"_"+str_pt[l];
	string ztdistname=centname+"ZTDIST_LARGE_MINUS_"+str_side[j]+"_"+str_pt[l];
	string decztdistname=centname+"DECZTDIST_LARGE_MINUS_"+str_side[j]+"_"+str_pt[l];
	string ztfinedistname=centname+"ZTDIST_MINUS_"+str_side[j]+"_"+str_pt[l];
	string decztfinedistname=centname+"DECZTDIST_MINUS_"+str_side[j]+"_"+str_pt[l];
	string poutdistname=centname+"POUTDIST_LARGE_MINUS_"+str_side[j]+"_"+str_pt[l];
	string decpoutdistname=centname+"DECPOUTDIST_LARGE_MINUS_"+str_side[j]+"_"+str_pt[l];
	


	se->registerHisto(xidistname.c_str(),XIDIST_LARGE_MINUS[j][l]=new TH2F(xidistname.c_str(),xidistname.c_str(),25,-1.0,4.,30,0,pi));      	

	se->registerHisto(xedistname.c_str(),XEDIST_LARGE_MINUS[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),40,0.,2.,30,0,pi));      
	se->registerHisto(ztfinedistname.c_str(),ZTDIST_MINUS[j][l]=new TH2F(ztfinedistname.c_str(),ztfinedistname.c_str(),40,0.,2.,30,0,pi));      
	//se->registerHisto(xedistname.c_str(),XEDIST_LARGE_MINUS[j][l]=new TH2F(xedistname.c_str(),xedistname.c_str(),25,0.,2.,30,0,pi));      
	se->registerHisto(ztdistname.c_str(),ZTDIST_LARGE_MINUS[j][l]=new TH2F(ztdistname.c_str(),ztdistname.c_str(),25,0.,2.,30,0,pi));      
	if(pouthistos) se->registerHisto(poutdistname.c_str(),POUTDIST_LARGE_MINUS[j][l]=new TH2F(poutdistname.c_str(),poutdistname.c_str(),20,0.,10.,30,0,pi));      

	XIDIST_LARGE_MINUS[j][l]->Sumw2();	
	XEDIST_LARGE_MINUS[j][l]->Sumw2();	
	ZTDIST_LARGE_MINUS[j][l]->Sumw2();	
	ZTDIST_MINUS[j][l]->Sumw2();	
	if(pouthistos) POUTDIST_LARGE_MINUS[j][l]->Sumw2();	
	
	if(tagflag>0){
	  se->registerHisto(decxidistname.c_str(),DECXIDIST_LARGE_MINUS[j][l]=new TH2F(decxidistname.c_str(),decxidistname.c_str(),25,-1.0,4.,30,0,pi));      	  
	  se->registerHisto(decxedistname.c_str(),DECXEDIST_LARGE_MINUS[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),40,0.,2.,30,0,pi));  
	  se->registerHisto(decztfinedistname.c_str(),DECZTDIST_MINUS[j][l]=new TH2F(decztfinedistname.c_str(),decztfinedistname.c_str(),40,0.,2.,30,0,pi));      
	  //se->registerHisto(decxedistname.c_str(),DECXEDIST_LARGE_MINUS[j][l]=new TH2F(decxedistname.c_str(),decxedistname.c_str(),25,0.,2.,30,0,pi)); 
	  se->registerHisto(decztdistname.c_str(),DECZTDIST_LARGE_MINUS[j][l]=new TH2F(decztdistname.c_str(),decztdistname.c_str(),25,0.,2.,30,0,pi));      
 	  if(pouthistos) se->registerHisto(decpoutdistname.c_str(),DECPOUTDIST_LARGE_MINUS[j][l]=new TH2F(decpoutdistname.c_str(),decpoutdistname.c_str(),20,0.,10.,30,0,pi));      
	  DECXEDIST_LARGE_MINUS[j][l]->Sumw2();
	  DECXIDIST_LARGE_MINUS[j][l]->Sumw2();
	  DECZTDIST_LARGE_MINUS[j][l]->Sumw2();
	  DECZTDIST_MINUS[j][l]->Sumw2();
	  if(pouthistos)DECPOUTDIST_LARGE_MINUS[j][l]->Sumw2();
	  
	}
      	
      }


    }

  }
 
  
  if(!skip_fg_histos){   
    string modmapz = centname+"_MODMAP0";
    se->registerHisto(modmapz.c_str(),MODMAP[0]=new TH3F(modmapz.c_str(),"MODMAP0",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string modmapone = centname+"_MODMAP1";
    se->registerHisto(modmapone.c_str(),MODMAP[1]=new TH3F(modmapone.c_str(),"MODMAP1",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string modmaptwo = centname+"_MODMAP2";
    se->registerHisto(modmaptwo.c_str(),MODMAP[2]=new TH3F(modmaptwo.c_str(),"MODMAP2",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string modmapthree = centname+"_MODMAP3";
    se->registerHisto(modmapthree.c_str(),MODMAP[3]=new TH3F(modmapthree.c_str(),"MODMAP3",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string modmapfour = centname+"_MODMAP4";
    se->registerHisto(modmapfour.c_str(),MODMAP[4]=new TH3F(modmapfour.c_str(),"MODMAP4",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string modmapfive = centname+"_MODMAP5";
    se->registerHisto(modmapfive.c_str(),MODMAP[5]=new TH3F(modmapfive.c_str(),"MODMAP5",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string modmapsix = centname+"_MODMAP6";
    se->registerHisto(modmapsix.c_str(),MODMAP[6]=new TH3F(modmapsix.c_str(),"MODMAP6",96,-0.5,95.5,48,-0.5,47.5,7,-0.5,6.5));
    string modmapsev = centname+"_MODMAP7";
    se->registerHisto(modmapsev.c_str(),MODMAP[7]=new TH3F(modmapsev.c_str(),"MODMAP7",96,-0.5,95.5,48,-0.5,47.5,7,-0.5,6.5));
    
    string smodmapz = centname+"_MODMAPSUB0";
    se->registerHisto(smodmapz.c_str(),MODMAPSUB[0]=new TH3F(smodmapz.c_str(),"MODMAPSUB0",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string smodmapone = centname+"_MODMAPSUB1";
    se->registerHisto(smodmapone.c_str(),MODMAPSUB[1]=new TH3F(smodmapone.c_str(),"MODMAPSUB1",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string smodmaptwo = centname+"_MODMAPSUB2";
    se->registerHisto(smodmaptwo.c_str(),MODMAPSUB[2]=new TH3F(smodmaptwo.c_str(),"MODMAPSUB2",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string smodmapthree = centname+"_MODMAPSUB3";
    se->registerHisto(smodmapthree.c_str(),MODMAPSUB[3]=new TH3F(smodmapthree.c_str(),"MODMAPSUB3",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string smodmapfour = centname+"_MODMAPSUB4";
    se->registerHisto(smodmapfour.c_str(),MODMAPSUB[4]=new TH3F(smodmapfour.c_str(),"MODMAPSUB4",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string smodmapfive = centname+"_MODMAPSUB5";
    se->registerHisto(smodmapfive.c_str(),MODMAPSUB[5]=new TH3F(smodmapfive.c_str(),"MODMAPSUB5",72,-0.5,71.5,36,-0.5,35.5,7,-0.5,6.5));
    string smodmapsix = centname+"_MODMAPSUB6";
    se->registerHisto(smodmapsix.c_str(),MODMAPSUB[6]=new TH3F(smodmapsix.c_str(),"MODMAPSUB6",96,-0.5,95.5,48,-0.5,47.5,7,-0.5,6.5));
    string smodmapsev = centname+"_MODMAPSUB7";
    se->registerHisto(smodmapsev.c_str(),MODMAPSUB[7]=new TH3F(smodmapsev.c_str(),"MODMAPSUB7",96,-0.5,95.5,48,-0.5,47.5,7,-0.5,6.5));
    
    string asfs = centname+"_ASSYM57";
    se->registerHisto(asfs.c_str(),ASSYM57=new TH1F(asfs.c_str(),"ASSYM57",100,0,1));
    string assn = centname+"_ASSYM79";
    se->registerHisto(assn.c_str(),ASSYM79=new TH1F(assn.c_str(),"ASSYM79",100,0,1));
    string asnt = centname+"_ASSYM912";
    se->registerHisto(asnt.c_str(),ASSYM912=new TH1F(asnt.c_str(),"ASSYM912",100,0,1));
    string astt = centname+"_ASSYM1220";
    se->registerHisto(astt.c_str(),ASSYM1220=new TH1F(astt.c_str(),"ASSYM1220",100,0,1));
    
    
    string twid = centname+"_TWRID";
    se->registerHisto(twid.c_str(),TWRID=new TH1F(twid.c_str(),"TWRID",24768,-0.5,24767.5));
    string twidmm = centname+"_TWRIDMYMAP";
    se->registerHisto(twidmm.c_str(),TWRIDMYMAP=new TH1F(twidmm.c_str(),"TWRIDMYMAP",24768,-0.5,24767.5));
  }
    
    string lead = centname+"_PHILEAD";
    se->registerHisto(lead.c_str(),PHILEAD=new TH1F(lead.c_str(),"PHILEAD",500,-pi/2.,3*pi/2.));
    string sub = centname+"_PHISUB";
    se->registerHisto(sub.c_str(),PHISUB=new TH1F(sub.c_str(),"PHISUB",500,-pi/2.,3*pi/2.));
    
    string dch = centname+"_DCH";
  se->registerHisto(dch.c_str(),DCH = new TH3F(dch.c_str(),"DCH",9,0.5,5.,500,-pi/2.,3*pi/2.,80,-80,80));

  if(!skip_fg_histos){
  
    string mind = centname+"_MINTRACKDIST";
    se->registerHisto(mind.c_str(),MINTRACKDIST=new TH3F(mind.c_str(),"MINTRACKDIST",100,0.,100.,20,0.,20.,40,0.,20.));
    string minthetadistname  = centname+"_MINTHETADIST";
    se->registerHisto(minthetadistname.c_str(),MINTHETADIST=new TH3F(minthetadistname.c_str(),minthetadistname.c_str(),100,0.,0.5,20,0.,20.,40,0.,20.));
    string minphidistname  = centname+"_MINPHIDIST";
    se->registerHisto(minphidistname.c_str(),MINPHIDIST=new TH3F(minphidistname.c_str(),minphidistname.c_str(),100,0.,0.5,20,0.,20.,40,0.,20.));
    string minraddistname  = centname+"_MINRADDIST";
    se->registerHisto(minraddistname.c_str(),MINRADDIST=new TH3F(minraddistname.c_str(),minraddistname.c_str(),100,0.,5.0,20,0.,20.,40,0.,20.));
    
  }


  for(int it=0;it<5;it++){
	
    char  awphiptname[100],  mwphiptname[100], mesptofname[100],mespemcname[100],mesmtofname[100],mesmemcname[100],barptofname[100],barpemcname[100],barmtofname[100],barmemcname[100];
    
    char   mwphifoldptname[100], mesptoffoldname[100],mespemcfoldname[100],mesmtoffoldname[100],mesmemcfoldname[100],barptoffoldname[100],barpemcfoldname[100],barmtoffoldname[100],barmemcfoldname[100];
    
    if(dofilltime){
      sprintf(mwphifoldptname,"_MWPHIFOLDPTCF_%d",it);
      string tempmwphifoldptnamef = centname+mwphifoldptname;
      se->registerHisto(tempmwphifoldptnamef.c_str(),MWPHIFOLDPTCF[it] = new TH2F(tempmwphifoldptnamef.c_str(),tempmwphifoldptnamef.c_str(),30,0,pi,100,0,10.));
      MWPHIFOLDPTCF[it]->Sumw2();
      
      sprintf(mwphifoldptname,"_MWPHIFOLDPTFLOW_%d",it);
      string tempmwphifoldptname = centname+mwphifoldptname;
      se->registerHisto(tempmwphifoldptname.c_str(),MWPHIFOLDPTFLOW[it] = new TH2F(tempmwphifoldptname.c_str(),tempmwphifoldptname.c_str(),30,0,pi,100,0,10.));
      MWPHIFOLDPTFLOW[it]->Sumw2();
    }


    sprintf(awphiptname,"_AWPHIPT%d",it);
    string cawphiptname = centname+awphiptname;
  se->registerHisto(cawphiptname.c_str(),AWPHIPT[it] = new TH2F(cawphiptname.c_str(),cawphiptname.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
    AWPHIPT[it]->Sumw2();
    
    sprintf(mwphiptname,"_MWPHIPT%d",it);
    string cmwphiptname = centname+mwphiptname;
  se->registerHisto(cmwphiptname.c_str(),MWPHIPT[it] = new TH2F(cmwphiptname.c_str(),mwphiptname,60,-pi/2.,3*pi/2.,100,0,10.));
    MWPHIPT[it]->Sumw2();

    sprintf(barpemcfoldname,"_MWPHIFOLDPT_BAR_P_EMC%d",it);
    string cbarpemcfoldname = centname+barpemcfoldname;
    se->registerHisto(cbarpemcfoldname.c_str(),MWPHIFOLDPT_BAR_P_EMC[it] = new TH2F(cbarpemcfoldname.c_str(),cbarpemcfoldname.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPT_BAR_P_EMC[it]->Sumw2();

    sprintf(barmemcfoldname,"_MWPHIFOLDPT_BAR_M_EMC%d",it);
    string cbarmemcfoldnamem = centname+barmemcfoldname;
    se->registerHisto(cbarmemcfoldnamem.c_str(),MWPHIFOLDPT_BAR_M_EMC[it] = new TH2F(cbarmemcfoldnamem.c_str(),cbarmemcfoldnamem.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPT_BAR_M_EMC[it]->Sumw2();

    if (!_useLessHistos)
      {
	sprintf(mwphiptname,"_MWPHIPTPLUS%d",it);
	string cmwphiptnamep = centname+mwphiptname;
	se->registerHisto(cmwphiptnamep.c_str(),MWPHIPTPLUS[it] = new TH2F(cmwphiptnamep.c_str(),cmwphiptnamep.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
	MWPHIPTPLUS[it]->Sumw2();
	
	sprintf(mwphiptname,"_MWPHIPTMINUS%d",it);
	string cmwphiptnamem = centname+mwphiptname;
	se->registerHisto(cmwphiptnamem.c_str(),MWPHIPTMINUS[it] = new TH2F(cmwphiptnamem.c_str(),cmwphiptnamem.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
	MWPHIPTMINUS[it]->Sumw2();
      }
    
    if(species<1){
      sprintf(mesptofname,"_MWPHIPT_MES_P_TOF%d",it);
      string cmesptofname = centname+mesptofname;
      se->registerHisto(cmesptofname.c_str(),MWPHIPT_MES_P_TOF[it] = new TH2F(cmesptofname.c_str(),cmesptofname.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
      MWPHIPT_MES_P_TOF[it]->Sumw2();
      
      sprintf(mespemcname,"_MWPHIPT_MES_P_EMC%d",it);
      string cmespemcname = centname+mespemcname;
      se->registerHisto(cmespemcname.c_str(),MWPHIPT_MES_P_EMC[it] = new TH2F(cmespemcname.c_str(),cmespemcname.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
      MWPHIPT_MES_P_EMC[it]->Sumw2();
      
      sprintf(mesmtofname,"_MWPHIPT_MES_M_TOF%d",it);
      string cmesmtofnamem = centname+mesmtofname;
      se->registerHisto(cmesmtofnamem.c_str(),MWPHIPT_MES_M_TOF[it] = new TH2F(cmesmtofnamem.c_str(),cmesmtofnamem.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
      MWPHIPT_MES_M_TOF[it]->Sumw2();
      
      sprintf(mesmemcname,"_MWPHIPT_MES_M_EMC%d",it);
      string cmesmemcnamem = centname+mesmemcname;
      se->registerHisto(cmesmemcnamem.c_str(),MWPHIPT_MES_M_EMC[it] = new TH2F(cmesmemcnamem.c_str(),cmesmemcnamem.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
      MWPHIPT_MES_M_EMC[it]->Sumw2();
      
      
      sprintf(barptofname,"_MWPHIPT_BAR_P_TOF%d",it);
      string cbarptofname = centname+barptofname;
      se->registerHisto(cbarptofname.c_str(),MWPHIPT_BAR_P_TOF[it] = new TH2F(cbarptofname.c_str(),cbarptofname.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
      MWPHIPT_BAR_P_TOF[it]->Sumw2();
      
      sprintf(barpemcname,"_MWPHIPT_BAR_P_EMC%d",it);
      string cbarpemcname = centname+barpemcname;
      se->registerHisto(cbarpemcname.c_str(),MWPHIPT_BAR_P_EMC[it] = new TH2F(cbarpemcname.c_str(),cbarpemcname.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
      MWPHIPT_BAR_P_EMC[it]->Sumw2();
      
      sprintf(barmtofname,"_MWPHIPT_BAR_M_TOF%d",it);
      string cbarmtofnamem = centname+barmtofname;
      se->registerHisto(cbarmtofnamem.c_str(),MWPHIPT_BAR_M_TOF[it] = new TH2F(cbarmtofnamem.c_str(),cbarmtofnamem.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
      MWPHIPT_BAR_M_TOF[it]->Sumw2();
      
      sprintf(barmemcname,"_MWPHIPT_BAR_M_EMC%d",it);
      string cbarmemcnamem = centname+ barmemcname;
      se->registerHisto(cbarmemcnamem.c_str(),MWPHIPT_BAR_M_EMC[it] = new TH2F(cbarmemcnamem.c_str(),cbarmemcnamem.c_str(),60,-pi/2.,3*pi/2.,100,0,10.));
      MWPHIPT_BAR_M_EMC[it]->Sumw2();
    }
    
   	
    sprintf(mwphifoldptname,"_MWPHIFOLDPT%d",it);
      string cmwphifoldptnamef = centname+mwphifoldptname;
  se->registerHisto(cmwphifoldptnamef.c_str(),MWPHIFOLDPT[it] = new TH2F(cmwphifoldptnamef.c_str(),cmwphifoldptnamef.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPT[it]->Sumw2();

    sprintf(mwphifoldptname,"_MWPHIFOLDPTQ2%d",it);
    string cmwphifoldptnamefq = centname+mwphifoldptname;
  se->registerHisto(cmwphifoldptnamefq.c_str(),MWPHIFOLDPTQ2[it] = new TH2F(cmwphifoldptnamefq.c_str(),cmwphifoldptnamef.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPTQ2[it]->Sumw2();


    if (!_useLessHistos){
      
    sprintf(mwphifoldptname,"_MWPHIFOLDPTPLUS%d",it);
    string cmwphifoldptnamep = centname+mwphifoldptname;
  se->registerHisto(cmwphifoldptnamep.c_str(),MWPHIFOLDPTPLUS[it] = new TH2F(cmwphifoldptnamep.c_str(),mwphifoldptname,30,0,pi,100,0,10.));
    MWPHIFOLDPTPLUS[it]->Sumw2();
    
    sprintf(mwphifoldptname,"_MWPHIFOLDPTMINUS%d",it);
      string cmwphifoldptnamefm = centname+mwphifoldptname;
  se->registerHisto(cmwphifoldptnamefm.c_str(),MWPHIFOLDPTMINUS[it] = new TH2F(cmwphifoldptnamefm.c_str(),cmwphifoldptnamefm.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPTMINUS[it]->Sumw2();

    }

    if(species<1){
      sprintf(mesptoffoldname,"_MWPHIFOLDPT_MES_P_TOF%d",it);
      string cmesptoffoldname = centname+mesptoffoldname;
      se->registerHisto(cmesptoffoldname.c_str(),MWPHIFOLDPT_MES_P_TOF[it] = new TH2F(cmesptoffoldname.c_str(),mesptoffoldname,30,0,pi,100,0,10.));
      MWPHIFOLDPT_MES_P_TOF[it]->Sumw2();
      
      sprintf(mespemcfoldname,"_MWPHIFOLDPT_MES_P_EMC%d",it);
      string cmespemcfoldname = centname+mespemcfoldname;
      se->registerHisto(cmespemcfoldname.c_str(),MWPHIFOLDPT_MES_P_EMC[it] = new TH2F(cmespemcfoldname.c_str(),cmespemcfoldname.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPT_MES_P_EMC[it]->Sumw2();
    
    sprintf(mesmtoffoldname,"_MWPHIFOLDPT_MES_M_TOF%d",it);
    string cmesmtoffoldnamem = centname+mesmtoffoldname;
    se->registerHisto(cmesmtoffoldnamem.c_str(),MWPHIFOLDPT_MES_M_TOF[it] = new TH2F(cmesmtoffoldnamem.c_str(),cmesmtoffoldnamem.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPT_MES_M_TOF[it]->Sumw2();
    
    sprintf(mesmemcfoldname,"_MWPHIFOLDPT_MES_M_EMC%d",it);
    string cmesmemcfoldnamem = centname+mesmemcfoldname;
    se->registerHisto(cmesmemcfoldnamem.c_str(),MWPHIFOLDPT_MES_M_EMC[it] = new TH2F(cmesmemcfoldnamem.c_str(),cmesmemcfoldnamem.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPT_MES_M_EMC[it]->Sumw2();
    
    
    sprintf(barptoffoldname,"_MWPHIFOLDPT_BAR_P_TOF%d",it);
    string cbarptoffoldname = centname+barptoffoldname;
    se->registerHisto(cbarptoffoldname.c_str(),MWPHIFOLDPT_BAR_P_TOF[it] = new TH2F(cbarptoffoldname.c_str(),cbarptoffoldname.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPT_BAR_P_TOF[it]->Sumw2();
    
    
    sprintf(barmtoffoldname,"_MWPHIFOLDPT_BAR_M_TOF%d",it);
    string cbarmtoffoldnamefm = centname+barmtoffoldname;
    se->registerHisto(cbarmtoffoldnamefm.c_str(),MWPHIFOLDPT_BAR_M_TOF[it] = new TH2F(cbarmtoffoldnamefm.c_str(),cbarmtoffoldnamefm.c_str(),30,0,pi,100,0,10.));
    MWPHIFOLDPT_BAR_M_TOF[it]->Sumw2();
    
    }
  }

  
  string awtrig = centname+"_AWTRIGPT";
  se->registerHisto(awtrig.c_str(),AWTRIGPT=new TH1F(awtrig.c_str(),"_AWTRIGPT",5,-0.5,4.5));
  string mwtrig = centname+"_MWTRIGPT";
  se->registerHisto(mwtrig.c_str(),MWTRIGPT=new TH1F(mwtrig.c_str(),"MWTRIGPT",5,-0.5,4.5));
  string mwfine = centname+"_MWFINETRIGPT";
  se->registerHisto(mwfine.c_str(),MWFINETRIGPT=new TH1F(mwfine.c_str(),"MWFINETRIGPT",7,-0.5,6.5));
  
  string pi0zemcname = centname+"_PI0ZEMC";
  se->registerHisto(pi0zemcname.c_str(),PI0ZEMC=new TH2F(pi0zemcname.c_str(),pi0zemcname.c_str(),400,-200,200,20,4,20));
  
  string tagcntname = centname + "_TAGCOUNTER";
  se->registerHisto(tagcntname.c_str() ,TAGCOUNTER=new TH1F(tagcntname.c_str(),tagcntname.c_str() ,200,0,20));
  string tag1cntname = centname + "_TAGCOUNTER1";
  se->registerHisto(tag1cntname.c_str() ,TAGCOUNTER1=new TH1F(tag1cntname.c_str(),tagcntname.c_str() ,200,0,20));
  string tag2cntname = centname + "_TAGCOUNTER2";
  se->registerHisto(tag2cntname.c_str() ,TAGCOUNTER2=new TH1F(tag2cntname.c_str(),tagcntname.c_str() ,200,0,20));

  // determine if this is a hard pDST for AuAu
  bool isHpdst = false;
  if ( (species==1||species==5) && _emcnodename->Contains("Hard"))
    isHpdst = true;

  string feff_location = toad_loader->location(_pi0effFilename->Data());
  SetTriggerEfficiency(feff_location.c_str());
  cout << "pi0 trigger efficiency loaded" << endl;

  string sharkfin_file = toad_loader->location(_sharkfinname->Data());
  SetSharkFin(sharkfin_file.c_str());
  cout << "sharkfin input loaded" << endl;

  //fill time corrs
  if(dofilltime)
  {
    string filltime_file = toad_loader->location(_filltimeFilename->Data());
    string filltimeup_file = toad_loader->location(_filltimeFilenameup->Data());
    string filltimedown_file = toad_loader->location(_filltimeFilenamedown->Data());

    SetFilltimeCorrs(filltime_file.c_str(),filltimeup_file.c_str(),filltimedown_file.c_str());
  }//filltimecorrs

  delete toad_loader;
  cout<< " finish init "<<endl;
  return 0;
}


int CombinedSimple::process_event(PHCompositeNode *topNode)
{

  // Get Data nodes
  // Global
  global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");

  if (!global) cout << PHWHERE << "CombinedSimple:: PHGlobal not in Node Tree" << endl;

  // TString member varibles
  particle = findNode::getClass<PHCentralTrack>(topNode, _cntnodename->Data());
  if (!particle) 
  {
    cout << PHWHERE << "CombinedSimple:: PHCentralTrack not in Node Tree" << endl;
    return -2;
  }

  emccluster = findNode::getClass<emcClusterContainer>(topNode,_emcnodename->Data() );      
  if(!emccluster)
  {
    cout << PHWHERE << "CombinedSimple::emcClusterContainer not in Node Tree"<<endl;
    return -2;
  }

  if(species > 0)
  {
    reacplane=NULL;
    //reacplane = getClass<ReactionPlaneObject>(topNode,"ReactionPlaneObject");
    //if(!reacplane){
    //cout << PHWHERE << "CombinedSimple::ReactionPlaneObject not in Node Tree"<<endl;
    //	return -2;
    //}
  }

  d_runhdr = getClass<RunHeader> (topNode, "RunHeader");
  if(!d_runhdr) {
    cout << PHWHERE << "CombinedSimple:: Run Header not in Node Tree" << endl;
    return 0;
  }
  evtsync = getClass<EventHeader> (topNode, "EventHeader");
  if(!evtsync){
    cout << PHWHERE << "CombinedSimple:: EventHeader not in Node Tree" << endl;
    return 0;
  }

  ErtOut * ertOut = NULL; 
  if(species==0 || species == 4){
    ertOut = getClass<ErtOut> (topNode, "ErtOut");
    if(!ertOut) {
      cout << PHWHERE << "CombinedSimple:: Run ErtOut not in Node Tree" << endl;
      return 0;        
    }
  }

  //Can't use trigger helper, assume 4x4a and 4x4b
  int ertTrigs[3] = {1,1,0};
  int isERT = 1;

  runno=d_runhdr->get_RunNumber();

  float zvertex = 0.0;

  if(runno >= 106935 && runno <= 122223)
    zvertex = global->getZVertex();
  else if(runno >= 168705)
    zvertex = global->getBbcZVertex();

  if (fabs(zvertex)>30) return 0;
  ZVERTEX->Fill(zvertex);

  as_zvertex = zvertex;
  as_runno=runno;

  int seqno = -1;
  float evtno = -1; 
  evtno = (float) evtsync->get_EvtSequence();
  as_evtno = evtsync->get_EvtSequence();
  seqno = as_evtno/100000; //SyncObject->SegmentNumber() depricated
  as_segno = seqno;
  int ncluster = emccluster->size();

  //Changed to standard warnmap for both
  //const static unsigned int sccut3x3Map = 0xffe1ce70;    
  //const static unsigned int glcut3x3Map = 0x1ce70;     
  //const static unsigned int sccut3x3Map = 0x1ce70;  
  //Changed to 0x7de1ce70 to cut on cold and warm too
  const static unsigned int glcut3x3Map = 0x7de1ce70;     
  const static unsigned int sccut3x3Map = 0x7de1ce70;

  unsigned int cut3x3Map;
  
  float ngoodpi0s=0;
  
  float mass_window[2]={0.12,0.16};
  if(tagflag==2){
    mass_window[0]=0.53;
    mass_window[1]=0.58;
  }
  if(tagflag==3){ // loband
    mass_window[0]=0.065;
    mass_window[1]=0.115;
  }
  if(tagflag==4){ // hiband
    mass_window[0]=0.165;
    mass_window[1]=0.200;
  }


  float bbcqs, bbcqn, bbcq;
  float zdcen, zdces;
  float percent = 0;
  float thetaBBC=0;
  float thetaDEG=0;
  //int icbin_pw=0;

  //include run4 AuAu as species == 2, run 8 dAu == 3 (MB) and 4(ERT)
  if(species==1 || species==2 || species==3 || species==4 || species==5)
  {
    bbcqs = global->getBbcChargeS();
    bbcqn = global->getBbcChargeN();
    bbcq = 0.5 * (bbcqs + bbcqn);
    zdcen = global->getZdcEnergyN();
    zdces = global->getZdcEnergyS();


    if(species==2)percent = PhUtilities::getCentralityByClockRun4(bbcqn,bbcqs,zdcen,zdces,runno);
    else
      percent = global->getCentrality();

    if (percent<=0)
      percent = 0;
    if (percent>100)
      percent = -1;

    if(!(percent > locent && percent <= hicent))
      return 0;

    if(reacplane)
      thetaBBC =reacplane->getBBCrp12();

    float thetaDEG = thetaBBC * DEG_PER_RAD;

    // todo:  uncomment this when we start using RxP binning in Au+Au again
    if(thetaDEG<-90||thetaDEG>90) //return 0;
      cout << "reaction plane angle thetaDEG not in [-90,90] degrees" << endl;

    RPLANE->Fill(thetaDEG);
    CENTRALITY->Fill(percent);
  }

  as_centrality = (int) percent;
  event_counter++;
  _assocNum = -1;

  // even if _savAssoc == 0, we want makeMN called, because
  // this creates the arrays now needed for loop thru tracks:  however
  // it's ok to not call it here if _savAssoc == 0, because it will be 
  // called there again if it hasn't been  here  
  if (_savAssoc > 0 && (event_counter % _savAssoc == 0 || event_counter <3)) 
    makeMNSingles(particle);


  for( int j=0; j<ncluster; j++) 
  {
    int ispi0=0;
    photon = emccluster->getCluster(j);

    float ecore1=photon->ecore();
    if(ecore1<1.0) continue;

    float pc3dz = photon->emcpc3dz();
    float pc3dphi = photon->emcpc3dphi();

    //circular cut of 10 dz, 0.02 dphi
    float pc3dr=sqrt(pc3dphi*pc3dphi/0.0004+pc3dz*pc3dz/100);

    int sector1 = photon->sector();
    //arm = 0 -> west
    int arm1 = photon->arm();
    if (arm1==1) sector1 = 7 - sector1;
    if (sector1 < 6) cut3x3Map = sccut3x3Map;
    else cut3x3Map = glcut3x3Map;

    if (!Chi2Cut(photon,zvertex))
      continue;

    int iypos1 = photon->iypos();
    int izpos1 = photon->izpos();
    int towerid1 = photon->towerid(0);

    float x1 = photon->x();
    float y1 = photon->y();
    float z1 = photon->z();
    float radius1 = sqrt( x1*x1 + y1*y1 + (z1-zvertex)*(z1-zvertex) );

    //find z at reference radius
    float zemcsub = z1*510.0/sqrt(x1*x1+y1*y1);

    if(isEdgeTower(sector1,izpos1,iypos1))
      continue;

    if(species==1||species==2||species==5)
      if(CheckHotTowerMap(sector1, izpos1, iypos1, 0))
	continue;

    //check hot tower map
    if(species == 3 || species ==4)
      if(isHotTower(towerid1, species))
	continue;

    if(ecore1>2)
      TWRIDMYMAP->Fill(towerid1); 
    if((photon->warnmap() & cut3x3Map) ==0 && (photon->deadmap() & cut3x3Map) ==0)
    {
      //dead and cold towers: deadmap() &  0x7de1ce70
      //hot and warm towers:  warnmap() &  0x7de1ce70

      float px1 = ecore1 * x1 / radius1;
      float py1 = ecore1 * y1 / radius1;
      float pz1 = ecore1 * (z1-zvertex) / radius1;

      float emctof=photon->tofcorr();

      int iecore1 = 0;
      if(ecore1>1&&ecore1<3) iecore1=1;
      if(ecore1>3&&ecore1<5) iecore1=2;
      if(ecore1>5&&ecore1<7) iecore1=3;
      if(ecore1>7&&ecore1<9) iecore1=4;
      if(ecore1>9&&ecore1<12) iecore1=5;
      if(ecore1>12&&ecore1<20) iecore1=6;

      if(tagflag==0)
      {
	float pt1=sqrt(px1*px1+py1*py1);
	if(pt1>2)TWRID->Fill(towerid1);

	MODMAP[sector1]->Fill(izpos1,iypos1,iecore1);

	//fiducial cut on inclusive/leading photon 
	if(fabs(zemcsub)>155.0) continue;

	float phitemp=atan2(py1,px1);
	phitemp=PHAngle(phitemp);

	if(pt1>5&&pt1<20)
	{
	  //if (species==0&&run_5_6==1)
	  if (species==0 || species == 4)//pp and dAu ERT
	  {
	    if (isERT && CheckLvl1Fired(photon, ertOut, ertTrigs) == 0)
	    {
	      //cout << "trig photon failed lvl1 " << endl;
	      //lvl1failed->Fill(sqrt(px1*px1 + py1*py1));//remove for run8, CHC
	      continue;
	    }
	  }

	  PHILEAD->Fill(phitemp);

	  TriggerHolder trig_gamma;
	  trig_gamma.pt=pt1;
	  trig_gamma.E=ecore1;
	  trig_gamma.px=px1;
	  trig_gamma.py=py1;
	  trig_gamma.pz=pz1;
	  trig_gamma.the=acos(pz1/sqrt(px1*px1+py1*py1+pz1*pz1));;
	  trig_gamma.eta=-log(tan(trig_gamma.the/2.));
	  trig_gamma.phi=phitemp;
	  trig_gamma.arm=arm1;
	  trig_gamma.sector=sector1;
	  trig_gamma.x=x1;
	  trig_gamma.y=y1;
	  trig_gamma.z=z1;
	  trig_gamma.pc3dr=pc3dr;
	  trig_gamma.emctof=emctof;
	  trig_gamma.zvertex=zvertex;

	  int passvetocut = VetoTracks(trig_gamma, NULL, emccluster,0);

	  if(useiso)
	  {
	    float sumincone =SumEcorePtInCone(trig_gamma.the, trig_gamma.phi, zvertex, j, -1, trig_gamma.pt, trig_gamma.E);
	    sumincone=SumEcorePtInCone(trig_gamma.the, trig_gamma.phi+0.4, zvertex, j, -1, trig_gamma.pt+20.0, trig_gamma.E);

	    if(useiso && sumincone>1.1*trig_gamma.E)
	      continue;
	  }

	  //remove tagged photons
	  if(removetags && passvetocut)
	  {
	    for(int k=0; k<ncluster; k++)
	    {
	      if(k==j) continue;
	      photon2 = emccluster->getCluster(k);

	      if(photon2->prob_photon()<0.02) continue;

	      int sector2 = photon2->sector();
	      //arm = 0 -> west
	      int arm2 = photon2->arm();
	      if(arm1==arm2)
	      {
		if (arm2==1) sector2 = 7 - sector2;
		if (sector2 < 6) cut3x3Map = sccut3x3Map;
		else cut3x3Map = glcut3x3Map;

		float ecore2 = photon2->ecore();
		if(ecore2<1.0) continue;

		int iypos2 = photon2->iypos();
		int izpos2 = photon2->izpos();
		int towerid2=photon2->towerid(0);

		if(species==1||species==2||species==5) if(CheckHotTowerMap(sector2, izpos2, iypos2, 0)) continue;

		//check hot tower map
		if(species == 3 || species ==4) if(isHotTower(towerid2, species)) continue;

		if((photon2->warnmap() & cut3x3Map) ==0 && (photon2->deadmap() & cut3x3Map) ==0 )
		{
		  float x2 = photon2->x();
		  float y2 = photon2->y();
		  float z2 = photon2->z();
		  float radius2 = sqrt( x2*x2 + y2*y2 + (z2-zvertex)*(z2-zvertex) );
		  float px2 = ecore2 * x2 / radius2;
		  float py2 = ecore2 * y2 / radius2;
		  float pz2 = ecore2 * (z2-zvertex) / radius2;

		  float mass2= 2.0*ecore1*ecore2 - 2.0*(px1*px2+py1*py2+pz1*pz2);
		  float invMass=0;
		  invMass=sqrt(mass2);

		  TLorentzVector pho1, pho2, pi0;
		  pho1.SetPxPyPzE(px1,py1,pz1,ecore1);
		  pho2.SetPxPyPzE(px2,py2,pz2,ecore2);
		  pi0=pho1+pho2;
		  // cout<<" pi0pt "<<pi0.Pt()<<" pt "<<sqrt((px1+px2)*(px1+px2)+(py1+py2)*(py1+py2))<<endl;
		  //INVMASS[arm1]->Fill(pi0.Pt(),invMass);

		  //INVMASS[0]->Fill(pt1,invMass);

		  if(ecore2>1.0) INVMASS[0]->Fill(pt1,invMass);
		  if(ecore2>1.5) INVMASS[1]->Fill(pt1,invMass);
		  if(ecore2>2.0) INVMASS[2]->Fill(pt1,invMass);
		  if(ecore2>2.5) INVMASS[3]->Fill(pt1,invMass);
		  if(ecore2>3.0) INVMASS[4]->Fill(pt1,invMass);

		  //make wide invariant mass for isolation cut
		  if((invMass>0.118&&invMass<0.162))
		  {  // JEF switched to just pi0 for AuAu for now

		    if(ispi0<1 && ecore2>1.0) ispi0=1; 
		    if(ispi0<2 && ecore2>1.5) ispi0=2;
		    if(ispi0<3 && ecore2>2.0) ispi0=3;
		    if(ispi0<4 && ecore2>2.5) ispi0=4;
		    if(ispi0<5 && ecore2>3.0) ispi0=5;

		    //if(ispi0==5) break;
		  }
		}
	      }
	    }
	  }

	  float sumptaway = LoopThruTracks(trig_gamma, NULL, emccluster,ispi0);

	  if(sumptaway>=0 && !ispi0)
	  {
	    TRIG_COUNTER->Fill(1);
	    float ntuplefill[]={runno,seqno,evtno,ecore1,-1,trig_gamma.eta,trig_gamma.the,trig_gamma.phi,zvertex,ngoodpi0s,trig_gamma.px,trig_gamma.py,trig_gamma.pz,x1,y1,z1,sector1,pc3dr,emctof,percent,thetaDEG,towerid1,-1};

	    TRIG_LOCATION->Fill(trig_gamma.eta, trig_gamma.phi);
	    triggaz->Fill(ntuplefill);
	  }
	}
      }
      else if(tagflag>0){
	// pi0 loop
	for(int k=0; k<ncluster; k++){
	  if(j==k) continue;
	  photon2 = emccluster->getCluster(k);
	  float ecore2 = photon2->ecore();
	  if(ecore2<1.0||ecore1<ecore2|| ecore1 == ecore2) continue;
	  int sector2 = photon2->sector();

	  float sume12 = ecore1+ecore2;
	  if (sume12 < 4.0) continue;

	  float asym12 =  (ecore1 - ecore2)/ sume12; // fabs not neccesary because upper > check
	  const float mina = 0.15;  //this is the smallest asym value used for the cut
	  //for the most central, which should increase to 1 at sume = 5.25 (this effectively means
	  // that at some point below 5.25, this cut is doing nothing due to the e > 1.0 cut above

	  // energy dep asym cut in au + au central events
	  if (percent <= 40 && 
	      sume12 < 5.25 &&
	      species > 0 &&
	      asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25))
	      )
	    continue;

	  //arm = 0 -> west
	  int arm2 = photon2->arm();
	  if(arm1==arm2){
	    if (arm2==1) sector2 = 7 - sector2;
	    if (sector2 < 6) cut3x3Map = sccut3x3Map;
	    else cut3x3Map = glcut3x3Map;

	    int iypos2 = photon2->iypos();
	    int izpos2 = photon2->izpos();

	    int towerid2=photon2->towerid(0);

	    if(isEdgeTower(sector2,izpos2,iypos2)) continue;

	    if(species==1||species==2||species==5) if(CheckHotTowerMap(sector2, izpos2, iypos2, 0)) continue;

	    //check hot tower map
	    if(species == 3 || species ==4) if(isHotTower(towerid2, species)) continue;

	    if( (photon2->warnmap() & cut3x3Map) ==0 && (photon2->deadmap() & cut3x3Map) ==0 ){

	      float x2 = photon2->x();
	      float y2 = photon2->y();
	      float z2 = photon2->z();
	      float radius2 = sqrt( x2*x2 + y2*y2 + (z2-zvertex)*(z2-zvertex) );

	      float px2 = ecore2 * x2 / radius2;
	      float py2 = ecore2 * y2 / radius2;
	      float pz2 = ecore2 * (z2-zvertex) / radius2;

	      float mass2= 2.0*ecore1*ecore2 - 2.0*(px1*px2+py1*py2+pz1*pz2);
	      float invMass=0;
	      invMass=sqrt(mass2);

	      TLorentzVector pho1, pho2, pi0;
	      pho1.SetPxPyPzE(px1,py1,pz1,ecore1);
	      pho2.SetPxPyPzE(px2,py2,pz2,ecore2);
	      pi0=pho1+pho2;
	      float pi0pt=pi0.Pt();
	      //float pi0p=pi0.P();


	      //new for veto tracks for invmass histos
	      TriggerHolder trig_pi0_test;
	      trig_pi0_test.E=pi0.E();
	      trig_pi0_test.px=pi0.Px();
	      trig_pi0_test.py=pi0.Py();
	      trig_pi0_test.pz=pi0.Pz();
	      trig_pi0_test.the=acos(pi0.Pz()/sqrt(pi0.Px()*pi0.Px()+pi0.Py()*pi0.Py()+pi0.Pz()*pi0.Pz()));;
	      trig_pi0_test.eta=-log(tan(trig_pi0_test.the/2.));
	      float phitemp=atan2(pi0.Py(),pi0.Px());
	      trig_pi0_test.phi=PHAngle(phitemp);
	      trig_pi0_test.pt=sqrt(pi0.Px()*pi0.Px()+pi0.Py()*pi0.Py());

	      //just hold one of the photons stuff here
	      trig_pi0_test.x=x1;
	      trig_pi0_test.y=y1;
	      trig_pi0_test.z=z1;
	      trig_pi0_test.arm=arm1;
	      trig_pi0_test.sector=sector1;
	      trig_pi0_test.pc3dr=pc3dr;
	      trig_pi0_test.emctof=emctof;
	      trig_pi0_test.zvertex=zvertex;

	      int passvetocut = VetoTracks(trig_pi0_test, NULL, emccluster,0);
	      if(passvetocut){
		if(removetags)
		{
		  float ptg1=sqrt(px1*px1+py1*py1);
		  float ptg2=sqrt(px2*px2+py2*py2);
		  INVMASS[0]->Fill(ptg1,invMass);
		  INVMASS[1]->Fill(ptg2,invMass);
		  INVMASS[2]->Fill(ptg1,invMass);
		  INVMASS[2]->Fill(ptg2,invMass);
		}else{
		  INVMASS[sector1]->Fill(pi0pt,invMass);
		}
	      }

	      if(pi0pt>4&&pi0pt<17)
	      {
		float pi0trigeff=grpi0eff->Eval(pi0pt);
		float pi0zemc  = 510.0*pi0.Pz()/pi0pt+zvertex;   
		int ipi0zemc = (int)TMath::Floor((pi0zemc+165.0)/10.);
		if(ipi0zemc<0)ipi0zemc=0;
		if(ipi0zemc>32)ipi0zemc=32;

		for(int idecs=0;idecs<7;idecs++){
		  //cout<<" pi0pt "<<pi0pt<<" pi0zemc "<<pi0zemc<<" idecs "<<idecs<<" ipi0zemc "<<ipi0zemc<<endl;
		  int trigptbin = hshark_alt[0][0]->FindBin(pi0pt);
		  float mattshark=hshark_alt[idecs][ipi0zemc]->GetBinContent(trigptbin); 

		  float mwweightfine=mattshark*pi0trigeff;
		  if(mwweightfine>0){
		    DECINVMASS->Fill(idecs,invMass,mwweightfine);
		  }
		}

		//if(species==0&&run_5_6==1){
		if(species==0 ||species ==4){
		  emcClusterContent * phocheck = photon;
		  if (pho1.E() < pho2.E()) phocheck = photon2;

		  if (isERT && CheckLvl1Fired(phocheck, ertOut, ertTrigs) == 0)
		  {
		    //cout << "trig pi0 failed lvl1 " << endl;
		    //lvl1failed->Fill(pi0pt); remove for run8, CHC
		    continue;
		  }
		}

		if(invMass>mass_window[0]&&invMass<mass_window[1]){

		  int iecore2 = 0;
		  if(ecore2>1&&ecore2<3) iecore2=1;
		  if(ecore2>3&&ecore2<5) iecore2=2;
		  if(ecore2>5&&ecore2<7) iecore2=3;
		  if(ecore2>7&&ecore2<9) iecore2=4;
		  if(ecore2>9&&ecore2<12) iecore2=5;
		  if(ecore2>12&&ecore2<20) iecore2=6;

		  MODMAP[sector1]->Fill(izpos1,iypos1,iecore1);
		  MODMAPSUB[sector2]->Fill(izpos2,iypos2,iecore2);

		  ngoodpi0s++;
		  TriggerHolder trig_pi0;
		  trig_pi0.E=pi0.E();
		  trig_pi0.px=pi0.Px();
		  trig_pi0.py=pi0.Py();
		  trig_pi0.pz=pi0.Pz();
		  trig_pi0.the=acos(pi0.Pz()/sqrt(pi0.Px()*pi0.Px()+pi0.Py()*pi0.Py()+pi0.Pz()*pi0.Pz()));;
		  trig_pi0.eta=-log(tan(trig_pi0.the/2.));
		  float phitemp=atan2(pi0.Py(),pi0.Px());
		  trig_pi0.phi=PHAngle(phitemp);
		  trig_pi0.pt=sqrt(pi0.Px()*pi0.Px()+pi0.Py()*pi0.Py());

		  //just hold one of the photons stuff here
		  trig_pi0.x=x1;
		  trig_pi0.y=y1;
		  trig_pi0.z=z1;
		  trig_pi0.arm=arm1;
		  trig_pi0.sector=sector1;
		  trig_pi0.pc3dr=pc3dr;
		  trig_pi0.emctof=emctof;
		  trig_pi0.zvertex=zvertex;

		  float sumptaway = LoopThruTracks(trig_pi0, NULL, emccluster,0);

		  if(sumptaway>=0){

		    TRIG_COUNTER->Fill(1);

		    float ntuplefill[23]={runno,seqno,evtno,trig_pi0.E,ecore2,trig_pi0.eta,trig_pi0.the,trig_pi0.phi,zvertex,ngoodpi0s,trig_pi0.px,trig_pi0.py,trig_pi0.pz,x1,y1,z1,trig_pi0.sector,pc3dr,emctof,percent,thetaDEG,towerid1,towerid2};

		    TRIG_LOCATION->Fill(trig_pi0.eta, trig_pi0.phi);
		    triggaz->Fill(ntuplefill);

		    if(trig_pi0.pt>5&&trig_pi0.pt<7)ASSYM57->Fill(fabs(ecore1-ecore2)/(ecore1+ecore2));
		    if(trig_pi0.pt>7&&trig_pi0.pt<9)ASSYM79->Fill(fabs(ecore1-ecore2)/(ecore1+ecore2));
		    if(trig_pi0.pt>9&&trig_pi0.pt<12)ASSYM912->Fill(fabs(ecore1-ecore2)/(ecore1+ecore2));
		    if(trig_pi0.pt>12&&trig_pi0.pt<20)ASSYM1220->Fill(fabs(ecore1-ecore2)/(ecore1+ecore2));
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  event_counter++;
  return 0;
}


int CombinedSimple::VetoTracks(TriggerHolder trigger, PHCentralTrack *bgparticle, emcClusterContainer * emccluster, int thres){

  //int tag_thres = thres;
  
  int npart=0;

  if(!bgparticle)
    {
      //cout << "using particle and running fg" <<endl;
      npart = particle->get_npart();
    }else{
      //cout << "must be running bg code" <<endl;
      particle = bgparticle;
      npart = particle->get_npart();
    }


  //int ngoodtracks=0;
  //float sumptaway=0.;
  
  float phi_r =  trigger.phi;
  phi_r = PHAngle(phi_r);
  //float eta_r =  trigger.eta;  
  
  //99 is the edge of the histogram, I want to see count how many land outside
  float mindist=99;
  float ptofvetotrack=0.;
  float minthetadist=99.;
  float minphidist=99.;

  for(int i=0;i<npart;i++)
    {
      //PHSnglCentralTrack* trk= particle->get_track(i);
 
      int quality=particle->get_quality(i);
      float the = particle->get_the0(i);
      if(the<-99) continue;
      float mom = particle->get_mom(i);
      float pt = mom*sin(the);      
      if(quality>7){
	float dist=99;
	//stay away from recal edge
	//float phidch = particle->get_phi(i);
	float zed = particle->get_zed(i);
	int side = 0;
	if(zed>0) side=1;
	//DCH->Fill(pt,phidch,zed); //DCH->Fill(pt,phidch,side);

	//try using the pc3 projection
	float xtrk = particle->get_ppc3x(i);      
	float ytrk = particle->get_ppc3y(i);
	float ztrk = particle->get_ppc3z(i);

	if(pt>_vetoPtCut){
	  float thetadist, phidist;
	  dist=FindTrackDistance(trigger, xtrk, ytrk, ztrk, thetadist, phidist);
	  //add a histogram to check dist distribution

	  if(dist<mindist) {
	    mindist=dist;
	    ptofvetotrack = pt;
	    minthetadist=thetadist;
	    minphidist=phidist;
	  }
	}
      }
    }

  float minovalcut =sqrt(minphidist*minphidist/(0.014*0.014)+minthetadist*minthetadist/(0.021*0.021));
  
  //does it necessary for dAu?
  if(species==0){
    if(minovalcut<1.1){
      VETO_COUNTER_0->Fill(trigger.pt);
      return 0; 
    }
  }
  else{
    if(mindist<8){      
      VETO_COUNTER_0->Fill(trigger.pt);
      return 0; 
    }
  }

  //passed veto cut:
  return 1;
  
}

float CombinedSimple::LoopThruTracks(TriggerHolder trigger, PHCentralTrack *bgparticle, emcClusterContainer * emccluster, int thres){
  
  int tag_thres = thres;
  float radpt[20]={0.0};
  float radcut[20]={0.0};

  
  int ispbgl=0;
  int sector =  (int)trigger.sector;
  //int garm=trigger.arm;
  //if(garm==1) sector=7-sector;
  if(sector>5) ispbgl=1;
  
  int npart=0;

  if(!bgparticle)
    {
      //cout << "using particle and running fg" <<endl;
      npart = particle->get_npart();
    }else{
      //cout << "must be running bg code" <<endl;
      particle = bgparticle;
      npart = particle->get_npart();
    
    }


  int ngoodtracks=0;
  float sumptaway=0.;
  
 
  
  float phi_r =  trigger.phi;
  phi_r = PHAngle(phi_r);
  //float eta_r =  trigger.eta;  
  
  
  
  //99 is the edge of the histogram, I want to see count how many land outside
  float mindist=99;
  float ptofvetotrack=0.;
  float minthetadist=99.;
  float minphidist=99.;


  for(int i=0;i<npart;i++)
    {
      //PHSnglCentralTrack* trk= particle->get_track(i);
 
      //cout << i << " of " << npart << " quality " << particle->get_quality(i) << " " << particle->get_mom(i) << " " << particle->get_the0(i) <<endl;

      int quality=particle->get_quality(i);
      float the = particle->get_the0(i);
      if(the<-99) continue;
      float mom = particle->get_mom(i);
      float pt = mom*sin(the);      

      //cout << "quality " << quality << " " << mom << " " << the <<endl;

      if(quality>7){
	float dist=99;
	//stay away from recal edge
	float phidch = particle->get_phi(i);
	float zed = particle->get_zed(i);
	int side = 0;
	if(zed>0) side=1;
	DCH->Fill(pt,phidch,zed); //DCH->Fill(pt,phidch,side);

	//cout << "pt: " << pt << " veto: " << _vetoPtCut <<endl;

	//try using the pc3 projection
	float xtrk = particle->get_ppc3x(i);      
	float ytrk = particle->get_ppc3y(i);
	float ztrk = particle->get_ppc3z(i);


	//if(pt>_vetoPtCut){
	  float thetadist, phidist;
	  dist=FindTrackDistance(trigger, xtrk, ytrk, ztrk, thetadist, phidist);

	  //float raddist =sqrt(phidist*phidist/(0.014*0.014)+thetadist*thetadist/(0.021*0.021));
	  float raddist =sqrt(phidist*phidist+thetadist*thetadist);

	  //cout << "raddist " <<raddist << " phidist " << phidist << "thetadist" << thetadist <<endl;

	  for(int b=0; b<20; b++){
	    radcut[b]=1.0/20.0*(b+1);
	    if(raddist<radcut[b]){
	      radpt[b]=radpt[b]+pt;
	      //cout << raddist << " < " << radcut[b] << " pt: " << pt << "->" << radpt[b] <<endl;
	    }
	  }
	
	  if(dist<mindist && pt>_vetoPtCut) {
	    mindist=dist;
	    ptofvetotrack = pt;
	    minthetadist=thetadist;
	    minphidist=phidist;
	  }

	  //}//pt cut
      }

    }
  

  if(!skip_fg_histos){
    MINPHIDIST->Fill(minphidist,trigger.pt,ptofvetotrack);
    MINTHETADIST->Fill(minthetadist,trigger.pt,ptofvetotrack);
    MINTRACKDIST->Fill(mindist,trigger.pt,ptofvetotrack);
  }

  float minovalcut =sqrt(minphidist*minphidist/(0.014*0.014)+minthetadist*minthetadist/(0.021*0.021));
  if(!skip_fg_histos) MINRADDIST->Fill(minovalcut,trigger.pt,ptofvetotrack);


  if(species==0){
    if(minovalcut<1.1){
      VETO_COUNTER->Fill(trigger.pt);
      return -1; //used to be return 0;
    }
  }
  else{
    if(mindist<8){      
      VETO_COUNTER->Fill(trigger.pt);
      return -1; //used to be return 0;
    }
  }
  
  if(removetags){
    if(tag_thres<1) TRIGPT->Fill(trigger.pt);
    if(ispbgl && tag_thres<1) TRIGPTGL->Fill(trigger.pt);
    if(!ispbgl && tag_thres<1) TRIGPTSC->Fill(trigger.pt);

    if(tag_thres<2) TRIGPT1->Fill(trigger.pt);
    if(tag_thres<3) TRIGPT2->Fill(trigger.pt);
    if(tag_thres<4) TRIGPT3->Fill(trigger.pt);
    if(tag_thres<5) TRIGPT4->Fill(trigger.pt);

    if(tag_thres>0) TAGCOUNTER->Fill(trigger.pt);
    if(tag_thres>1) TAGCOUNTER1->Fill(trigger.pt);
    if(tag_thres>2) TAGCOUNTER2->Fill(trigger.pt);
  }else{
    TRIGPT->Fill(trigger.pt); 
    if(ispbgl) TRIGPTGL->Fill(trigger.pt);
    if(!ispbgl) TRIGPTSC->Fill(trigger.pt);
    ptvscent_trig->Fill(trigger.pt,as_centrality);
  }    

  /*
  if(trigger.pc3dr>1)TRIGPT1->Fill(trigger.pt);
  if(trigger.pc3dr>1.25)TRIGPT2->Fill(trigger.pt);
  if(trigger.emctof<6.5)TRIGPT3->Fill(trigger.pt);
  if(trigger.pc3dr>1.25&&trigger.emctof<6.5)TRIGPT4->Fill(trigger.pt);	
  */  

  float awweight[5]={0}; 
  float mwweight[5]={0};//, mwetaweight[4];
  
  float mwweightfine[7];

  float pi0trigeff=1.; 
  //float pi0trigefflo=grpi0efflo->Eval(trigger.pt);
  //float pi0trigeffhi=grpi0effhi->Eval(trigger.pt);  

    
  if(tagflag>0)pi0trigeff=grpi0eff->Eval(trigger.pt);

  
  if(tagflag>0)TRIGPTEFFW->Fill(trigger.pt, pi0trigeff);


  int trigptbin = hshark_large[0][0]->FindBin(trigger.pt);
  if(trigptbin>400) trigptbin=400;

  float pi0zemc  = 510.0*trigger.pz/trigger.pt+trigger.zvertex;
  int ipi0zemc = (int)TMath::Floor((pi0zemc+165.0)/10.);

  PI0ZEMC->Fill(pi0zemc,trigger.pt);

  if(ipi0zemc<0||ipi0zemc>32){
    //cout<<" zemc out of range "<<pi0zemc<<" bin = "<<ipi0zemc<<endl;
  }

  if(ipi0zemc<0)ipi0zemc=0;
  if(ipi0zemc>32)ipi0zemc=32;


  if(tagflag>0){
    
    for(int idecl=0;idecl<5;idecl++){
      
      float ptblo=5;
      float ptbhi=7;
      if(idecl==1){
	ptblo=7;
	ptbhi=9;	    
      }
      if(idecl==2){
	ptblo=9;
	ptbhi=12;	    
      }
      if(idecl==3){
	ptblo=5;
	ptbhi=10;	    
      }
      if(idecl==4){
	ptblo=12;
	ptbhi=15;	    
      }
      
      float anashark=0;
      
      if(trigger.pt>ptblo&&trigger.pt<ptbhi){
	anashark=1-ptblo/trigger.pt;
      }
      else if(trigger.pt>ptbhi){
	anashark=(ptbhi-ptblo)/trigger.pt;
      }
 
      
      float mattshark=hshark_large[idecl][ipi0zemc]->GetBinContent(trigptbin);
    
      
      if(mattshark>0) {
	
	awweight[idecl]=anashark*pi0trigeff;      
	mwweight[idecl]=mattshark*pi0trigeff;
	
	if(awweight[idecl]>0)AWTRIGPT->Fill(idecl,awweight[idecl]);
	if(mwweight[idecl]>0)MWTRIGPT->Fill(idecl,mwweight[idecl]);


	//ptvscent_dec->Fill();   
      }
    }
    for(int idecs=0;idecs<7;idecs++){      
      float mattshark=hshark_small[idecs][ipi0zemc]->GetBinContent(trigptbin); 
      //float mattshark=0; 
      mwweightfine[idecs]=mattshark*pi0trigeff;
      if(mwweightfine[idecs]>0)MWFINETRIGPT->Fill(idecs,mwweightfine[idecs]);
    }
    
  }
  
  if (_assocNum == -1) // then makeMN wasn't called
    {  
	makeMNSingles(particle,0);
    }	 	    

  

  PartnerHolder partnertrack;

  for(int i=0;i<_assocNum;i++){

        
    float mom = sqrt(as_px[i]*as_px[i]+as_py[i]*as_py[i]+as_pz[i]*as_pz[i]);
    
    float the0 = acos(as_pz[i]/mom);
    float phi0 = atan2(as_py[i],as_px[i]);
    int dcarm = 0;
    if(phi0<1.5) dcarm=1;
  
 
    partnertrack.vm2emc=as_m2emc[i];
    partnertrack.vm2tof=as_m2tof[i];
    partnertrack.vcharge=as_charge[i] / abs(as_charge[i]);
    partnertrack.vppc3x=as_ppc3x[i];      
    partnertrack.vppc3y=as_ppc3y[i];
    partnertrack.vppc3z=as_ppc3z[i];
    partnertrack.vthe0=the0;
    partnertrack.vmom=mom;
    partnertrack.vphi0=phi0;
    partnertrack.vdcarm=dcarm;
    partnertrack.vn0=(int)as_n0[i];
    partnertrack.vquality=abs(as_charge[i]);

    partnertrack.vpc3dphi=0.;
    partnertrack.vpc3dz=0.;
    partnertrack.vemcdphi=0.;
    partnertrack.vemcdz=0.;
    //      partnertrack.vn0=0;
    //	partnertrack.vquality=63;
    partnertrack.vpc3sdz=0.;
    partnertrack.vpc3sdphi=0.;
    partnertrack.vemcsdz=0.;
    partnertrack.vemcsdphi=0.;
    


    int trackused=LoopThruTrack(trigger, partnertrack,tag_thres);
    
    if(trackused) ngoodtracks++;
    
  }

  ngoodtracks = _assocNum;
  
  
  return sumptaway;
}


  

int CombinedSimple::LoopThruTrack(TriggerHolder trigger,  PartnerHolder partner,int tag_thres){
  

  //This is for ztdist filltime stuff
  int trigbinbounds_small[8]={5,6,7,8,9,12,15,20};  
  int trigbinbounds_large[6]={5,7,9,12,15,20};  

  int its=-1;
  int itl=-1;
  

  for(int it=0;it<7;it++) {
    if(trigger.pt>trigbinbounds_small[it]&&trigger.pt<trigbinbounds_small[it+1]) {
      its=it;	
      break;
    }
  }


  for(int it=0;it<5;it++) {
    if(trigger.pt>trigbinbounds_large[it]&&trigger.pt<trigbinbounds_large[it+1]) {
      itl=it;	
      break;
    }
  }



  //Start of old stuff
  float phi_r =  trigger.phi;
  phi_r = PHAngle(phi_r);
  //float eta_r =  trigger.eta;  
  
  int ispbgl=0;
  int sector = (int)trigger.sector;
  //int garm=trigger.arm;
  //if(garm==1) sector=7-sector;
  if(sector>5) ispbgl=1;
  
  float mwweight[5]={0};
  float awweight[5]={0};
  
  float mwweightfine[7]={0};

  float pi0trigeff=0.;
  //cout<<" tagflag "<<tagflag<<endl;
  if(tagflag>0)pi0trigeff=grpi0eff->Eval(trigger.pt);


  int trigptbin = hshark_large[0][0]->FindBin(trigger.pt);
  if(trigptbin>400) trigptbin=400;

  float pi0zemc  = 510.0*trigger.pz/trigger.pt+trigger.zvertex;
  int ipi0zemc = (int)TMath::Floor((pi0zemc+165.0)/10.);
  if(ipi0zemc<0||ipi0zemc>32){
    //cout<<" zemc out of range "<<pi0zemc<<" bin = "<<ipi0zemc<<endl;
  }
  if(ipi0zemc<0)ipi0zemc=0;
  if(ipi0zemc>32)ipi0zemc=32;      


  if(tagflag>0){
    
    for(int idecl=0;idecl<5;idecl++){
      
      float ptblo=5;
      float ptbhi=7;
      if(idecl==1){
	ptblo=7;
	ptbhi=9;	    
      }
      if(idecl==2){
	ptblo=9;
	ptbhi=12;	    
      }
      if(idecl==3){
	ptblo=5;
	ptbhi=10;	    
      }
      if(idecl==4){
	ptblo=12;
	ptbhi=15;	    
      }
      
      float anashark=0;
      
      if(trigger.pt>ptblo&&trigger.pt<ptbhi){
	anashark=1-ptblo/trigger.pt;
      }
      else if(trigger.pt>ptbhi){
	anashark=(ptbhi-ptblo)/trigger.pt;
      }
      
      float mattshark=hshark_large[idecl][ipi0zemc]->GetBinContent(trigptbin);
      
      if(mattshark>0) {
	
	awweight[idecl]=anashark*pi0trigeff;      
	mwweight[idecl]=mattshark*pi0trigeff;
	
// 	if(awweight[idecl]>0)AWTRIGPT->Fill(idecl,awweight[idecl]);
// 	if(mwweight[idecl]>0)MWTRIGPT->Fill(idecl,mwweight[idecl]);   
      }
    }
    for(int idecs=0;idecs<7;idecs++){      
      float mattshark=hshark_small[idecs][ipi0zemc]->GetBinContent(trigptbin); 
      mwweightfine[idecs]=mattshark*pi0trigeff;
//       if(mwweightfine[idecs]>0)MWFINETRIGPT->Fill(idecs,mwweightfine[idecs]);
    }


    
  }
  
  
  //loop thru single track only

  int quality=partner.vquality;
  float the = partner.vthe0;
  if(the<-99) return 1;
  float mom = partner.vmom;

  float pt = mom*sin(the);
  //partner cuts no longer applied here
  //float pc3sdphi=partner.vpc3sdphi;
  //float pc3sdz=partner.vpc3sdz;
  //float emcsdphi=partner.vemcsdphi;
  //float emcsdz=partner.vemcsdz;


  if (quality > 7) {
    // fabs(pc3sdphi)<3.0&&fabs(pc3sdz)<3.0&&fabs(emcsdphi)<3.0&&fabs(emcsdz)<3.0&&pt>0.5&&pt<10.0 
    // matching and basic cuts are in makeMNSingles now
    
    //test effect of n0 cut in run4 // it's ok to leave this
    // cut in for run7 even if n0's aren't in the mnsingles for r7 yet
    //  because if no n0, the default init value of n0 should pass the cut

    float n0 = partner.vn0;

    if(n0>=0 && pt<5.0) return 1; 
    //testing removing the n0 cut like Andrew
    //n0 cut added back June 7, 2010
    
    float m2emc=partner.vm2emc;
    float m2tof=partner.vm2tof;
    
    M2TOF->Fill(m2tof,pt);
    M2EMC->Fill(m2emc,pt);
  
    float dist=99;
    float thetadist, phidist;
    if(pt>_vetoPtCut){
      //	  dist=FindTrackDistance(trigger,trk);
      dist=FindTrackDistance(trigger,partner, thetadist, phidist);
    }

    float phi = partner.vphi0;
    phi = PHAngle(phi);
    //float eta = -log(tan(the/2.0));
    //float arm = partner.vdcarm;
    //int iarm = (int)arm;
    
    //unfolded version
    float dphi=phi-phi_r;
    dphi=PHAngle(dphi);

    // test for Justin  --nonpair cut pt1pt2:  todo: make a real histo to contain non-pair cut histos?  -- done, Matt
    if(!removetags)PT1PT2DPHI1->Fill(trigger.pt,pt,dphi);
  

    //todo:  Justin 3/16/07 shouldn't we make this cut on the photons, not the pi0 for tagflag > 0?
    // for now leave it alone.
    //Actually, the veto cut is applied to the leading photon since x,y,z are filled with leading photon values instead of pi0.  zvertex is calculated slightly wrong I think but I doubt it's important - Matt 03/19/08
    //Justin:  cool.  I guess we can erase this comment.

    if(species==0){
      if(sqrt(phidist*phidist/(0.014*0.014)+thetadist*thetadist/(0.021*0.021))<1.1) return 1;
    }
    else{
      //cout << "track distance is: "<< dist <<endl;
      if(dist<8){
	VETO_COUNTER->Fill(pt); 
	return 1;       
      }
    }
    
   //folded version	    
    float dphi_r=dphi;	
    if(dphi_r<0)dphi_r+=2*pi;
    else if(dphi_r>2*pi) dphi_r-=2*pi;	
    if(dphi_r>pi) dphi_r=2*pi-dphi_r;
    

     
    int iside=0;
    if(dphi_r>pi/2.) iside=1;

    int charge=(int)partner.vcharge;



     
    if(removetags)
      {
	if(tag_thres<1) PT1PT2DPHI->Fill(trigger.pt,pt,dphi);
	if(tag_thres<2) PT1PT2DPHI1->Fill(trigger.pt,pt,dphi);   
	if(tag_thres<3) PT1PT2DPHI2->Fill(trigger.pt,pt,dphi);    
	if(tag_thres<4) PT1PT2DPHI3->Fill(trigger.pt,pt,dphi);    
	if(tag_thres<5) PT1PT2DPHI4->Fill(trigger.pt,pt,dphi);
	//if(tag_thres<1 && !ispbgl) PT1PT2DPHI3->Fill(trigger.pt,pt,dphi);    
	//if(tag_thres<1 && ispbgl) PT1PT2DPHI4->Fill(trigger.pt,pt,dphi);

	if(charge==1 && tag_thres<1) PT1PT2DPHIPLUS->Fill(trigger.pt,pt,dphi);
	if(charge==-1 && tag_thres<1)PT1PT2DPHIMINUS->Fill(trigger.pt,pt,dphi);

      }else{
	PT1PT2DPHI->Fill(trigger.pt,pt,dphi);
	if(charge==1)PT1PT2DPHIPLUS->Fill(trigger.pt,pt,dphi);
	if(charge==-1)PT1PT2DPHIMINUS->Fill(trigger.pt,pt,dphi);

      }

    
    //if (!_useLessHistos) {
    //if(charge==1)PT1PT2DPHIPLUS->Fill(trigger.pt,pt,dphi);
    //if(charge==-1)PT1PT2DPHIMINUS->Fill(trigger.pt,pt,dphi);
    //}
  
    float dphifold=dphi;
    if(dphifold<0) dphifold+=2*pi;
    if(dphifold>2*pi) dphifold-=2*pi;
    if(dphifold>pi) dphifold=2*pi-dphifold;
    
    if(dphifold<0||dphifold>pi) cout<<" dphifold out of bounds "<<endl;

    if(removetags)
      {
	if(tag_thres<1) PT1PT2DPHIFOLD->Fill(trigger.pt,pt,dphifold);
	if(tag_thres<1 && !ispbgl) PT1PT2DPHIFOLDSC->Fill(trigger.pt,pt,dphifold);    
	if(tag_thres<1 && ispbgl) PT1PT2DPHIFOLDGL->Fill(trigger.pt,pt,dphifold);

	if(charge==1 && tag_thres<1) PT1PT2DPHIFOLDPLUS->Fill(trigger.pt,pt,dphifold);
	if(charge==-1 && tag_thres<1) PT1PT2DPHIFOLDMINUS->Fill(trigger.pt,pt,dphifold);
      }else{
	PT1PT2DPHIFOLD->Fill(trigger.pt,pt,dphifold);
	if(charge==1)PT1PT2DPHIFOLDPLUS->Fill(trigger.pt,pt,dphifold);
	if(charge==-1)PT1PT2DPHIFOLDMINUS->Fill(trigger.pt,pt,dphifold);

	if(!ispbgl) PT1PT2DPHIFOLDSC->Fill(trigger.pt,pt,dphifold);    
	if(ispbgl) PT1PT2DPHIFOLDGL->Fill(trigger.pt,pt,dphifold);

      }



    //if (!_useLessHistos) {
    //if(charge==1)PT1PT2DPHIFOLDPLUS->Fill(trigger.pt,pt,dphifold);
    //if(charge==-1)PT1PT2DPHIFOLDMINUS->Fill(trigger.pt,pt,dphifold);
    //}


    if(species<1){
      if(m2tof>-0.05&&m2tof<0.35){
	if(charge==1)PT1PT2DPHI_MES_P_TOF->Fill(trigger.pt,pt,dphi);
	if(charge==-1)PT1PT2DPHI_MES_M_TOF->Fill(trigger.pt,pt,dphi);
      }
      else if(m2emc>-0.05&&m2emc<0.35){
	if(charge==1)PT1PT2DPHI_MES_P_EMC->Fill(trigger.pt,pt,dphi);
	if(charge==-1)PT1PT2DPHI_MES_M_EMC->Fill(trigger.pt,pt,dphi);
      }
      if(m2tof>0.7&&m2tof<1.05){
	if(charge==1)PT1PT2DPHI_BAR_P_TOF->Fill(trigger.pt,pt,dphi);
	if(charge==-1)PT1PT2DPHI_BAR_M_TOF->Fill(trigger.pt,pt,dphi);
      }
      else if(m2emc>0.7&&m2emc<1.05){
	if(charge==1)PT1PT2DPHI_BAR_P_EMC->Fill(trigger.pt,pt,dphi);
	if(charge==-1)PT1PT2DPHI_BAR_M_EMC->Fill(trigger.pt,pt,dphi);
      }
      
      if(m2tof>-0.05&&m2tof<0.35){
	if(charge==1)PT1PT2DPHIFOLD_MES_P_TOF->Fill(trigger.pt,pt,dphifold);
	if(charge==-1)PT1PT2DPHIFOLD_MES_M_TOF->Fill(trigger.pt,pt,dphifold);
      }
      else if(m2emc>-0.05&&m2emc<0.35){
	if(charge==1)PT1PT2DPHIFOLD_MES_P_EMC->Fill(trigger.pt,pt,dphifold);
	if(charge==-1)PT1PT2DPHIFOLD_MES_M_EMC->Fill(trigger.pt,pt,dphifold);
      }
      if(m2tof>0.7&&m2tof<1.05){
	if(charge==1)PT1PT2DPHIFOLD_BAR_P_TOF->Fill(trigger.pt,pt,dphifold);
	if(charge==-1)PT1PT2DPHIFOLD_BAR_M_TOF->Fill(trigger.pt,pt,dphifold);
      }
      else if(m2emc>0.7&&m2emc<1.05){
	if(charge==1)PT1PT2DPHIFOLD_BAR_P_EMC->Fill(trigger.pt,pt,dphifold);
	if(charge==-1)PT1PT2DPHIFOLD_BAR_M_EMC->Fill(trigger.pt,pt,dphifold);
      }
    }



    if(tagflag>0){
      for(int ipw=0;ipw<5;ipw++){
        
        //cout<<" ipw "<<ipw<<" mww "<<mwweight[ipw]<<endl;
        if(mwweight[ipw]>0) {
          
          MWPHIPT[ipw]->Fill(dphi,pt,mwweight[ipw]);
          
          if (!_useLessHistos) {
            if(charge==1)MWPHIPTPLUS[ipw]->Fill(dphi,pt,mwweight[ipw]);
            if(charge==-1)MWPHIPTMINUS[ipw]->Fill(dphi,pt,mwweight[ipw]);
          }
          
          MWPHIFOLDPT[ipw]->Fill(dphifold,pt,mwweight[ipw]);
          
          if (!_useLessHistos) {
            if(charge==1)MWPHIFOLDPTPLUS[ipw]->Fill(dphifold,pt,mwweight[ipw]);
            if(charge==-1)MWPHIFOLDPTMINUS[ipw]->Fill(dphifold,pt,mwweight[ipw]);
          }
          
          if(species<1){
            if(m2tof>-0.05&&m2tof<0.35){
              if(charge==1)MWPHIPT_MES_P_TOF[ipw]->Fill(dphi,pt,mwweight[ipw]);
              if(charge==-1)MWPHIPT_MES_M_TOF[ipw]->Fill(dphi,pt,mwweight[ipw]);
            }
            else if(m2emc>-0.05&&m2emc<0.35){
              if(charge==1)MWPHIPT_MES_P_EMC[ipw]->Fill(dphi,pt,mwweight[ipw]);
              if(charge==-1)MWPHIPT_MES_M_EMC[ipw]->Fill(dphi,pt,mwweight[ipw]);
            }
            if(m2tof>0.7&&m2tof<1.05){
              if(charge==1)MWPHIPT_BAR_P_TOF[ipw]->Fill(dphi,pt,mwweight[ipw]);
              if(charge==-1)MWPHIPT_BAR_M_TOF[ipw]->Fill(dphi,pt,mwweight[ipw]);
            }
            else if(m2emc>0.7&&m2emc<1.05){
              if(charge==1)MWPHIPT_BAR_P_EMC[ipw]->Fill(dphi,pt,mwweight[ipw]);
              if(charge==-1)MWPHIPT_BAR_M_EMC[ipw]->Fill(dphi,pt,mwweight[ipw]);
            }
            
            if(m2tof>-0.05&&m2tof<0.35){
              if(charge==1)MWPHIFOLDPT_MES_P_TOF[ipw]->Fill(dphifold,pt,mwweight[ipw]);
              if(charge==-1)MWPHIFOLDPT_MES_M_TOF[ipw]->Fill(dphifold,pt,mwweight[ipw]);
            }
            else if(m2emc>-0.05&&m2emc<0.35){
              if(charge==1)MWPHIFOLDPT_MES_P_EMC[ipw]->Fill(dphifold,pt,mwweight[ipw]);
              if(charge==-1)MWPHIFOLDPT_MES_M_EMC[ipw]->Fill(dphifold,pt,mwweight[ipw]);
            }
            if(m2tof>0.7&&m2tof<1.05){
              if(charge==1)MWPHIFOLDPT_BAR_P_TOF[ipw]->Fill(dphifold,pt,mwweight[ipw]);
              if(charge==-1)MWPHIFOLDPT_BAR_M_TOF[ipw]->Fill(dphifold,pt,mwweight[ipw]);
            }
            else if(m2emc>0.7&&m2emc<1.05){
              if(charge==1)MWPHIFOLDPT_BAR_P_EMC[ipw]->Fill(dphifold,pt,mwweight[ipw]);
              if(charge==-1)MWPHIFOLDPT_BAR_M_EMC[ipw]->Fill(dphifold,pt,mwweight[ipw]);
            }
            
          }//species<1
        }//mwweight>0
      }//ipw
      
    }//tagflag>0
    
  
    //In matt's code this is in loopthrutracks but I think it's okay here
    

    // filltime stuff
    if(dofilltime){
      int centbin=0;

      //if(!(percent > locent && percent <= hicent))
      if(as_centrality > 20 && as_centrality <= 40) centbin=1;
      else if(as_centrality > 40 && as_centrality <= 60) centbin=2;
      else if(as_centrality > 60 && as_centrality <= 93) centbin=3;


      float xe=fabs(pt*cos(dphi)/trigger.pt);
      
      float pout = fabs(pt*sin(dphi));
      
      float zt=pt/trigger.pt;
      
      //changed definition of xi
      //float xi=log(1.0/xe);
      float xi=log(1.0/zt);
      
      //cout<<" xe "<<xe<<endl;
      
      int itbin = 0;
      
      if(trigger.pt>5&&trigger.pt<7) itbin=0;
      if(trigger.pt>7&&trigger.pt<9) itbin=1;
      if(trigger.pt>9&&trigger.pt<12) itbin=2;
      if(trigger.pt>12&&trigger.pt<15) itbin=3;
      
      int bin3d = DPHIPT_NORM[itbin]->FindBin(dphi,pt);
      
      float richcorr[4]={0.680,0.835,0.925,0.975};
      //float richcorr=0.680;
      
      //From Andrew's thesis for 0-20% without Rich embedding
      fexemb->SetParameters(0.761,1.640,-4.734);
      
      float embcorr[4]={(0.779+0.851)/2.,(0.887+.938)/2.,(0.968+0.988)/2.,(.986+.993+.998)/3.};

      if(pt>5.0 && centbin==0)embcorr[0]=fexemb->Eval(pt);

      
      float seffcorr = 1.0;
      //seffcorr = heff->GetBinContent(bin1d);
      seffcorr = 2.0/feff->Eval(pt);
      if(pt>5.0) seffcorr = seffcorr/embcorr[centbin];
      else seffcorr = seffcorr/richcorr[centbin];
      
      float filltimeweight=1.;
      float filltimeweighterr=0.;
      float filltimeflowerr=0.;
      float filltimeup=1.;
      float filltimedown=1.;
      
      float accw = DPHIPT_NORM[itbin]->GetBinContent(bin3d);
      if(accw>0)filltimeweight=seffcorr/accw;
      accw= DPHIPT_FLOW[itbin]->GetBinContent(bin3d);
      float filltimeflow=seffcorr*accw;

      float weighterr = 0.;
      if(accw>0) weighterr = seffcorr*DPHIPT_NORM[itbin]->GetBinError(bin3d)/accw/accw;
      float flowerr= seffcorr*DPHIPT_FLOW[itbin]->GetBinError(bin3d);


      float accwup = DPHIPT_NORM_PLUS[itbin]->GetBinContent(bin3d);
      if(accwup>0)filltimeup=seffcorr/accwup;
      accwup= DPHIPT_FLOW_PLUS[itbin]->GetBinContent(bin3d);
      float filltimeflowup=seffcorr*accwup;

      float accwdown = DPHIPT_NORM_MINUS[itbin]->GetBinContent(bin3d);
      if(accwdown>0)filltimedown=seffcorr/accwdown;
      accwdown= DPHIPT_FLOW_MINUS[itbin]->GetBinContent(bin3d);
      float filltimeflowdown=seffcorr*accwdown;
      

      PT1PT2ZT->Fill(trigger.pt,pt,zt);
      PT1PT2XI->Fill(trigger.pt,pt,xi);
      
      
      if(its>-1&&its<7){	    

	XEDIST[iside][its]->Fill(xe,dphifold,filltimeweight);
	ZTDIST[iside][its]->Fill(zt,dphifold,filltimeweight);
	XIDIST[iside][its]->Fill(xi,dphifold,filltimeweight);
	
	XEDIST[iside+2][its]->Fill(xe,dphifold,filltimeflow);
	ZTDIST[iside+2][its]->Fill(zt,dphifold,filltimeflow);
	XIDIST[iside+2][its]->Fill(xi,dphifold,filltimeflow);
	
	if(pouthistos_more){
	  if(pt>poutlothresh)POUTDIST[iside][its]->Fill(pout,dphifold,filltimeweight);
	  if(pt>poutlothresh)POUTDIST[iside+2][its]->Fill(pout,dphifold,filltimeflow);
	  
	  if(pt>poutlothresh&&pt<pouthithresh){
	    POUT2DIST[iside][its]->Fill(pout*pout,dphifold,filltimeweight);
	    pout2sum[iside][its]+=pout*pout*filltimeweight;
	    pout2weight[iside][its]+=filltimeweight;
	    XEVSPOUT2[iside][its]->Fill(xe,pout*pout,filltimeweight);
	    ZTVSPOUT2[iside][its]->Fill(zt,pout*pout,filltimeweight);	    
	    
	  }
	  
	  if(pt>poutlothresh&&pt<pouthithresh){
	    POUT2DIST[iside+2][its]->Fill(pout*pout,dphifold,filltimeflow);
	    pout2sum[iside+2][its]+=pout*pout*filltimeflow;
	    pout2weight[iside+2][its]+=filltimeflow;
	    XEVSPOUT2[iside+2][its]->Fill(xe,pout*pout,filltimeflow);
	    ZTVSPOUT2[iside+2][its]->Fill(zt,pout*pout,filltimeflow);	    
	    
	  }	    
	}
	
      }
      
      if(itl>-1&&itl<4){	    
	PT1PT2DPHIFOLDCF->Fill(trigger.pt,pt,dphifold,filltimeweight);
	PT1PT2DPHIFLOW->Fill(trigger.pt,pt,dphifold,filltimeflow);
	
	XEDIST_LARGE[iside][itl]->Fill(xe,dphifold,filltimeweight);
	XIDIST_LARGE[iside][itl]->Fill(xi,dphifold,filltimeweight);
	ZTDIST_LARGE[iside][itl]->Fill(zt,dphifold,filltimeweight);
	
	XEDIST_LARGE[iside+2][itl]->Fill(xe,dphifold,filltimeflow);
	XIDIST_LARGE[iside+2][itl]->Fill(xi,dphifold,filltimeflow);
	ZTDIST_LARGE[iside+2][itl]->Fill(zt,dphifold,filltimeflow); 


	XEDIST_STATERR[iside][itl]->Fill(xe,dphifold,weighterr);
	XIDIST_STATERR[iside][itl]->Fill(xi,dphifold,weighterr);
	ZTDIST_STATERR[iside][itl]->Fill(zt,dphifold,weighterr);
	
	XEDIST_STATERR[iside+2][itl]->Fill(xe,dphifold,flowerr);
	XIDIST_STATERR[iside+2][itl]->Fill(xi,dphifold,flowerr);
	ZTDIST_STATERR[iside+2][itl]->Fill(zt,dphifold,flowerr); 
	
	
	XEDIST_LARGE_PLUS[iside][itl]->Fill(xe,dphifold,filltimeup);
	XIDIST_LARGE_PLUS[iside][itl]->Fill(xi,dphifold,filltimeup);
	ZTDIST_LARGE_PLUS[iside][itl]->Fill(zt,dphifold,filltimeup);
	ZTDIST_PLUS[iside][itl]->Fill(zt,dphifold,filltimeup);
	
	XEDIST_LARGE_PLUS[iside+2][itl]->Fill(xe,dphifold,filltimeflowup);
	XIDIST_LARGE_PLUS[iside+2][itl]->Fill(xi,dphifold,filltimeflowup);
	ZTDIST_LARGE_PLUS[iside+2][itl]->Fill(zt,dphifold,filltimeflowup); 
	ZTDIST_PLUS[iside+2][itl]->Fill(zt,dphifold,filltimeflowup); 
	
	XEDIST_LARGE_MINUS[iside][itl]->Fill(xe,dphifold,filltimedown);
	XIDIST_LARGE_MINUS[iside][itl]->Fill(xi,dphifold,filltimedown);
	ZTDIST_LARGE_MINUS[iside][itl]->Fill(zt,dphifold,filltimedown);
	ZTDIST_MINUS[iside][itl]->Fill(zt,dphifold,filltimedown);
	
	XEDIST_LARGE_MINUS[iside+2][itl]->Fill(xe,dphifold,filltimeflowdown);
	XIDIST_LARGE_MINUS[iside+2][itl]->Fill(xi,dphifold,filltimeflowdown);
	ZTDIST_LARGE_MINUS[iside+2][itl]->Fill(zt,dphifold,filltimeflowdown); 
	ZTDIST_MINUS[iside+2][itl]->Fill(zt,dphifold,filltimeflowdown); 
	
	XEDIST_LARGE_NOW[iside][itl]->Fill(xe,dphifold);
	XIDIST_LARGE_NOW[iside][itl]->Fill(xi,dphifold);
	ZTDIST_LARGE_NOW[iside][itl]->Fill(zt,dphifold);
	
	
	if(pouthistos){
	  if(pt>poutlothresh){
	    POUTDIST_LARGE[iside][itl]->Fill(pout,dphifold,filltimeweight);
	    POUTDIST_LARGE[iside+2][itl]->Fill(pout,dphifold,filltimeflow);
	    
	    POUTDIST_LARGE_PLUS[iside][itl]->Fill(pout,dphifold,filltimeup);
	    POUTDIST_LARGE_MINUS[iside][itl]->Fill(pout,dphifold,filltimedown);
	    
	    POUTDIST_LARGE_PLUS[iside+2][itl]->Fill(pout,dphifold,filltimeflowup);
	    POUTDIST_LARGE_MINUS[iside+2][itl]->Fill(pout,dphifold,filltimeflowdown);
	    
	  }
	  if(pt>poutlothresh&&pt<pouthithresh){	      	      	    
	    POUT2DIST_LARGE[iside][itl]->Fill(pout*pout,dphifold,filltimeweight);
	    pout2sum_large[iside][itl]+=pout*pout*filltimeweight;
	    pout2weight_large[iside][itl]+=filltimeweight;	      
	  }
	  
	  
	  if(pt>poutlothresh&&pt<pouthithresh){	      	      	    
	    POUT2DIST_LARGE[iside+2][itl]->Fill(pout*pout,dphifold,filltimeflow);
	    pout2sum_large[iside+2][itl]+=pout*pout*filltimeflow;
	    pout2weight_large[iside+2][itl]+=filltimeflow;	      
	  }
	  
	}
      }
    
    if(tagflag>0){	
      //need to loop over xe bins 2x.  
      //first time we count up all the decay weights
      //normalize decay weights to # of triggers entered in MWFINETRIGPT
      //loop over again and fill
      
      float decxebw = DECXEDIST[0][0]->GetBinWidth(1);
      float decztbw = DECZTDIST[0][0]->GetBinWidth(1);
      float decxibw = DECXIDIST[0][0]->GetBinWidth(1);	    
      
      
      for(int itdec=0;itdec<7;itdec++){
	
	float trigweightnorm = mwweightfine[itdec];	    
	float sumtrigweight =0.;
	
	int itdecl = 0;
	if(itdec==0||itdec==1) itdecl=0;
	if(itdec==2||itdec==3) itdecl=1;
	if(itdec==4) itdecl=2;
	if(itdec==5) itdecl=3;
	
	
	int bin3d = DPHIPT_DECNORM[itdecl]->FindBin(dphi,pt);	     
	


	filltimeweight=1.;
	
	
	accw = DPHIPT_DECNORM[itdecl]->GetBinContent(bin3d);
	if(accw>0)filltimeweighterr= seffcorr*DPHIPT_DECNORM[itdecl]->GetBinError(bin3d)/accw/accw;
	if(accw>0)filltimeweight=seffcorr/accw;
	
	accw = DPHIPT_DECFLOW[itdecl]->GetBinContent(bin3d);
	filltimeflowerr = seffcorr*DPHIPT_DECFLOW[itdecl]->GetBinError(bin3d);
	filltimeflow=accw*seffcorr;
	
	accwup = DPHIPT_DECNORM_PLUS[itbin]->GetBinContent(bin3d);
	if(accwup>0)filltimeup=seffcorr/accwup;
	accwup= DPHIPT_DECFLOW_PLUS[itbin]->GetBinContent(bin3d);
	filltimeflowup=seffcorr*accwup;
	
	accwdown = DPHIPT_DECNORM_MINUS[itbin]->GetBinContent(bin3d);
	if(accwdown>0)filltimedown=seffcorr/accwdown;
	accwdown= DPHIPT_DECFLOW_MINUS[itbin]->GetBinContent(bin3d);
	filltimeflowdown=seffcorr*accwdown;


	sumtrigweight =0.;	      

	MWPHIFOLDPTFLOW[itdecl]->Fill(pt,dphifold,filltimeflow);
	MWPHIFOLDPTCF[itdecl]->Fill(pt,dphifold,filltimeweight);


	
	for(int ixedecbin=0;ixedecbin<DECXEDIST[iside][itdec]->GetNbinsX();ixedecbin++){
	  float xedec=DECXEDIST[iside][itdec]->GetBinCenter(ixedecbin+1);
	  float xedeclo=xedec-decxebw/2.;
	  float xedechi=xedec+decxebw/2.;
	  
	  //float ptdec=fabs(pt / xedec * cos(dphi));
	  float ptdeclo=fabs(pt / xedechi * cos(dphi));
	  float ptdechi=fabs(pt / xedeclo * cos(dphi));
	  
	  if(ptdechi<trigbinbounds_small[itdec]||ptdeclo>trigbinbounds_small[itdec+1]) continue;
	  if(ptdeclo<trigbinbounds_small[itdec]) ptdeclo=trigbinbounds_small[itdec];
	  if(ptdechi>trigbinbounds_small[itdec+1]) ptdechi=trigbinbounds_small[itdec+1];
	  //int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
	  int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
	  int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
	  float fineweightave=0.;
	  for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
	    fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
	  }
	  
	  fineweightave*=(ptdechi-ptdeclo);
	  sumtrigweight+=fineweightave;
	}
      
	for(int ixedecbin=0;ixedecbin<DECXEDIST[iside][0]->GetNbinsX();ixedecbin++){
	  float xedec=DECXEDIST[iside][itdec]->GetBinCenter(ixedecbin+1);
	  float xedeclo=xedec-decxebw/2.;
	  float xedechi=xedec+decxebw/2.;
	  
	  //float ptdec=fabs(pt / xedec * cos(dphi));
	  float ptdeclo=fabs(pt / xedechi * cos(dphi));
	  float ptdechi=fabs(pt / xedeclo * cos(dphi));
	  
	  
	  if(ptdechi<trigbinbounds_small[itdec]||ptdeclo>trigbinbounds_small[itdec+1]) continue;
	  if(ptdeclo<trigbinbounds_small[itdec]) ptdeclo=trigbinbounds_small[itdec];
	  if(ptdechi>trigbinbounds_small[itdec+1]) ptdechi=trigbinbounds_small[itdec+1];
	  //int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
	  int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
	  int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
	  float fineweightave=0.;
	  for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
	    fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
	  }
	  
	  
	  fineweightave*=(ptdechi-ptdeclo);
	  if(sumtrigweight>0)fineweightave*=trigweightnorm/sumtrigweight;
	  
	  
	  DECXEDIST[iside][itdec]->Fill(xedec,dphifold,filltimeweight*fineweightave);
	  if(pouthistos_more) DECXEVSPOUT2[iside][itdec]->Fill(xedec,pout*pout,filltimeweight*fineweightave);
	  
	  DECXEDIST[iside+2][itdec]->Fill(xedec,dphifold,filltimeflow*fineweightave);
	  if(pouthistos_more) DECXEVSPOUT2[iside+2][itdec]->Fill(xedec,pout*pout,filltimeflow*fineweightave);
	  
	}
	
	
	sumtrigweight =0.;
	
	for(int iztdecbin=0;iztdecbin<DECZTDIST[iside][itdec]->GetNbinsX();iztdecbin++){
	  float ztdec=DECZTDIST[iside][itdec]->GetBinCenter(iztdecbin+1);
	  float ztdeclo=ztdec-decztbw/2.;
	  float ztdechi=ztdec+decztbw/2.;
	  
	  //need to include the effect of decay kinematics here:	    
	  //float ptdec=fabs(pt / ztdec);
	  float ptdeclo=fabs(pt / ztdechi);
	  float ptdechi=fabs(pt / ztdeclo);
	  if(ptdechi<trigbinbounds_small[itdec]||ptdeclo>trigbinbounds_small[itdec+1]) continue;
	  
	  
	  
	  if(ptdeclo<trigbinbounds_small[itdec]) ptdeclo=trigbinbounds_small[itdec];
	  if(ptdechi>trigbinbounds_small[itdec+1]) ptdechi=trigbinbounds_small[itdec+1];
	  //int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
	  int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
	  int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
	  float fineweightave=0.;
	  for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
	    fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
	  }
	  
	  //take into account the phase space for decay
	  fineweightave*=(ptdechi-ptdeclo);
	  sumtrigweight+=fineweightave;
	}
	
	
	for(int iztdecbin=0;iztdecbin<DECZTDIST[iside][0]->GetNbinsX();iztdecbin++){
	  float ztdec=DECZTDIST[iside][itdec]->GetBinCenter(iztdecbin+1);
	  float ztdeclo=ztdec-decztbw/2.;
	  float ztdechi=ztdec+decztbw/2.;
	  
	  //need to include the effect of decay kinematics here:	    
	  //float ptdec=fabs(pt / ztdec);
	  float ptdeclo=fabs(pt / ztdechi);
	  float ptdechi=fabs(pt / ztdeclo);
	  
	  if(ptdechi<trigbinbounds_small[itdec]||ptdeclo>trigbinbounds_small[itdec+1]) continue;
	  if(ptdeclo<trigbinbounds_small[itdec]) ptdeclo=trigbinbounds_small[itdec];
	  if(ptdechi>trigbinbounds_small[itdec+1]) ptdechi=trigbinbounds_small[itdec+1];
	  
	  //int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
	  int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
	  int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
	  float fineweightave=0.;
	  for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
	    fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
	  }
	  
	  fineweightave*=(ptdechi-ptdeclo);
	  if(sumtrigweight>0)fineweightave*=trigweightnorm/sumtrigweight;
	  
	  DECZTDIST[iside][itdec]->Fill(ztdec,dphifold,filltimeweight*fineweightave);
	  if(pouthistos_more) DECZTVSPOUT2[iside][itdec]->Fill(ztdec,pout*pout,filltimeweight*fineweightave);		
	  
	  
	  DECZTDIST[iside+2][itdec]->Fill(ztdec,dphifold,filltimeflow*fineweightave);
	  if(pouthistos_more) DECZTVSPOUT2[iside+2][itdec]->Fill(ztdec,pout*pout,filltimeflow*fineweightave);		
	  
	}
	
	
	sumtrigweight =0.;
	
	for(int ixidecbin=0;ixidecbin<DECXIDIST[iside][itdec]->GetNbinsX();ixidecbin++){
	  float xidec=DECXIDIST[iside][itdec]->GetBinCenter(ixidecbin+1);
	  float xideclo=xidec-decxibw/2.;
	  float xidechi=xidec+decxibw/2.;
	  
	  //need to include the effect of decay kinematics here:	    
	  //float ptdec=pt*exp(xidec);
	  //float ptdeclo=fabs(pt*cos(dphifold)*exp(xideclo));
	  //float ptdechi=fabs(pt*cos(dphifold)*exp(xidechi));

	  float ptdeclo=fabs(pt*exp(xideclo));
	  float ptdechi=fabs(pt*exp(xidechi));

	  
	  if(ptdechi<trigbinbounds_small[itdec]||ptdeclo>trigbinbounds_small[itdec+1]) continue;
	  
	  if(ptdeclo<trigbinbounds_small[itdec]) ptdeclo=trigbinbounds_small[itdec];
	  if(ptdechi>trigbinbounds_small[itdec+1]) ptdechi=trigbinbounds_small[itdec+1];
	  //int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
	  int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
	  int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
	  float fineweightave=0.;
	  for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
	    fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
	  }
	  
	  //take into account the phase space for decay
	  fineweightave*=(ptdechi-ptdeclo);
	  sumtrigweight+=fineweightave;
	}
	
	
	for(int ixidecbin=0;ixidecbin<DECXIDIST[iside][0]->GetNbinsX();ixidecbin++){
	  float xidec=DECXIDIST[iside][itdec]->GetBinCenter(ixidecbin+1);
	  float xideclo=xidec-decxibw/2.;
	  float xidechi=xidec+decxibw/2.;
	  
	  //float xi=log(1/zt);
	  //need to include the effect of decay kinematics here:	    
	  //float ptdec=pt*exp(xidec);
	  //float ptdeclo=fabs(pt*cos(dphifold)*exp(xideclo));
	  //float ptdechi=fabs(pt*cos(dphifold)*exp(xidechi));
	  
	  float ptdeclo=fabs(pt*exp(xideclo));
	  float ptdechi=fabs(pt*exp(xidechi));


	  if(ptdechi<trigbinbounds_small[itdec]||ptdeclo>trigbinbounds_small[itdec+1]) continue;
	  if(ptdeclo<trigbinbounds_small[itdec]) ptdeclo=trigbinbounds_small[itdec];
	  if(ptdechi>trigbinbounds_small[itdec+1]) ptdechi=trigbinbounds_small[itdec+1];
	  
	  //int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
	  int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
	  int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
	  float fineweightave=0.;
	  for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
	    fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
	  }
	  
	  fineweightave*=(ptdechi-ptdeclo);
	  if(sumtrigweight>0)fineweightave*=trigweightnorm/sumtrigweight;
	  
	  DECXIDIST[iside][itdec]->Fill(xidec,dphifold,filltimeweight*fineweightave);
	  
	  DECXIDIST[iside+2][itdec]->Fill(xidec,dphifold,filltimeflow*fineweightave);
	  
	}
	
	if(pouthistos_more){
	  if(pt>poutlothresh)DECPOUTDIST[iside][itdec]->Fill(pout,dphifold,filltimeweight*mwweightfine[itdec]);
	  if(pt>poutlothresh)DECPOUTDIST[iside+2][itdec]->Fill(pout,dphifold,filltimeflow*mwweightfine[itdec]);
	  // change limits?
	  if(pt>poutlothresh&&pt<pouthithresh){
	    
	    DECPOUT2DIST[iside][itdec]->Fill(pout*pout,dphifold,filltimeweight*mwweightfine[itdec]);
	    
	    decpout2sum[iside][itdec]+=pout*pout*filltimeweight*mwweightfine[itdec];
	    decpout2weight[iside][itdec]+=filltimeweight*mwweightfine[itdec];
	    
	    
	    DECPOUT2DIST[iside+2][itdec]->Fill(pout*pout,dphifold,filltimeflow*mwweightfine[itdec]);
		  
	    decpout2sum[iside+2][itdec]+=pout*pout*filltimeflow*mwweightfine[itdec];
	    decpout2weight[iside+2][itdec]+=filltimeflow*mwweightfine[itdec];
	    
	    
	    //if(charge==1)DECPOUTDIST_PLUS[iside][itdec]->Fill(pout,dphifold,filltimeweight*mwweightfine[itdec]);
	    //if(charge==-1)DECPOUTDIST_MINUS[iside][itdec]->Fill(pout,dphifold,filltimeweight*mwweightfine[itdec]);
	  }  
	}	
	
      }
    }
    
    // same as above but now for large bins
    
    if(tagflag>0){
      
      float decxebw = DECXEDIST_LARGE[0][0]->GetBinWidth(1);
      float decxibw = DECXIDIST_LARGE[0][0]->GetBinWidth(1);
      float decztbw = DECZTDIST_LARGE[0][0]->GetBinWidth(1);
      
      
      for(int itdec=0;itdec<4;itdec++){
	
	//skip over 5-10 bin
	int itdecl = itdec;
	if(itdecl==3) itdecl = 4;
	
	float sumtrigweight =0.;  
	
	float trigweightnorm = mwweight[itdecl];	    
	
	
	int bin3d = DPHIPT_DECNORM[itdec]->FindBin(dphi,pt);	     
	filltimeweight=1.;
	

	accw = DPHIPT_DECNORM[itdec]->GetBinContent(bin3d);
	filltimeweighterr = seffcorr*DPHIPT_DECNORM[itdec]->GetBinError(bin3d)/accw/accw;
	if(accw>0)filltimeweight=seffcorr/accw;
	
	accw = DPHIPT_DECFLOW[itbin]->GetBinContent(bin3d);
	filltimeflowerr = seffcorr*DPHIPT_DECFLOW[itbin]->GetBinError(bin3d);
	if(accw>0)filltimeflow=seffcorr*accw;
	
	accwup = DPHIPT_DECNORM_PLUS[itbin]->GetBinContent(bin3d);
	if(accwup>0)filltimeup=seffcorr/accwup;
	accwup= DPHIPT_DECFLOW_PLUS[itbin]->GetBinContent(bin3d);
	filltimeflowup=seffcorr*accwup;
	
	accwdown = DPHIPT_DECNORM_MINUS[itbin]->GetBinContent(bin3d);
	if(accwdown>0)filltimedown=seffcorr/accwdown;
	accwdown= DPHIPT_DECFLOW_MINUS[itbin]->GetBinContent(bin3d);
	filltimeflowdown=seffcorr*accwdown;


	
	sumtrigweight =0.;	      
	
	for(int ixedecbin=0;ixedecbin<DECXEDIST_LARGE[iside][itdec]->GetNbinsX();ixedecbin++){
	  
	  
	  float xedec=DECXEDIST_LARGE[iside][itdec]->GetBinCenter(ixedecbin+1);
	  float xedeclo=xedec-decxebw/2.;
	  float xedechi=xedec+decxebw/2.;
	  
	  //float ptdec=fabs(pt / xedec * cos(dphi));
	  float ptdeclo=fabs(pt / xedechi * cos(dphi));
	  float ptdechi=fabs(pt / xedeclo * cos(dphi));
	  
	  if(ptdechi<trigbinbounds_large[itdec]||ptdeclo>trigbinbounds_large[itdec+1]) continue;
	  if(ptdeclo<trigbinbounds_large[itdec]) ptdeclo=trigbinbounds_large[itdec];
	  if(ptdechi>trigbinbounds_large[itdec+1]) ptdechi=trigbinbounds_large[itdec+1];
	  //int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
	  int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
	  int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
	  float fineweightave=0.;
	  for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
	    fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
	  }
	  
	  fineweightave*=(ptdechi-ptdeclo);
	  sumtrigweight+=fineweightave;
	}
	      

	      for(int ixedecbin=0;ixedecbin<DECXEDIST_LARGE[iside][0]->GetNbinsX();ixedecbin++){
		float xedec=DECXEDIST_LARGE[iside][itdec]->GetBinCenter(ixedecbin+1);
		float xedeclo=xedec-decxebw/2.;
		float xedechi=xedec+decxebw/2.;
		
		//float ptdec=fabs(pt / xedec * cos(dphi));
		float ptdeclo=fabs(pt / xedechi * cos(dphi));
		float ptdechi=fabs(pt / xedeclo * cos(dphi));
		
		
		if(ptdechi<trigbinbounds_large[itdec]||ptdeclo>trigbinbounds_large[itdec+1]) continue;
		if(ptdeclo<trigbinbounds_large[itdec]) ptdeclo=trigbinbounds_large[itdec];
		if(ptdechi>trigbinbounds_large[itdec+1]) ptdechi=trigbinbounds_large[itdec+1];
		//int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
		int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
		int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
		float fineweightave=0.;
		for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
		  fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
		}
		 
		fineweightave*=(ptdechi-ptdeclo);
		if(sumtrigweight>0)fineweightave*=trigweightnorm/sumtrigweight;
		
			
		DECXEDIST_LARGE[iside][itdec]->Fill(xedec,dphifold,filltimeweight*fineweightave);

		DECXEDIST_LARGE[iside+2][itdec]->Fill(xedec,dphifold,filltimeflow*fineweightave);

		DECXEDIST_STATERR[iside][itdec]->Fill(xedec,dphifold,filltimeweighterr*fineweightave);

		DECXEDIST_STATERR[iside+2][itdec]->Fill(xedec,dphifold,filltimeflowerr*fineweightave);

		DECXEDIST_LARGE_NOW[iside][itdec]->Fill(xedec,dphifold,fineweightave);
	      


		DECXEDIST_LARGE_PLUS[iside][itdec]->Fill(xedec,dphifold,filltimeup*fineweightave);
		DECXEDIST_LARGE_MINUS[iside][itdec]->Fill(xedec,dphifold,filltimedown*fineweightave);
		DECXEDIST_LARGE_PLUS[iside+2][itdec]->Fill(xedec,dphifold,filltimeflowup*fineweightave);
		DECXEDIST_LARGE_MINUS[iside+2][itdec]->Fill(xedec,dphifold,filltimeflowdown*fineweightave);

	      }

	      
	      sumtrigweight =0.;	      
	      
	      for(int ixidecbin=0;ixidecbin<DECXIDIST_LARGE[iside][itdec]->GetNbinsX();ixidecbin++){
		float xidec=DECXIDIST_LARGE[iside][itdec]->GetBinCenter(ixidecbin+1);
		float xideclo=xidec-decxibw/2.;
		float xidechi=xidec+decxibw/2.;
		
		//float ptdeclo=fabs(pt*cos(dphifold)*exp(xideclo));
		//float ptdechi=fabs(pt*cos(dphifold)*exp(xidechi));

		float ptdeclo=fabs(pt*exp(xideclo));
		float ptdechi=fabs(pt*exp(xidechi));

		if(ptdechi<trigbinbounds_large[itdec]||ptdeclo>trigbinbounds_large[itdec+1]) continue;
		if(ptdeclo<trigbinbounds_large[itdec]) ptdeclo=trigbinbounds_large[itdec];
		if(ptdechi>trigbinbounds_large[itdec+1]) ptdechi=trigbinbounds_large[itdec+1];
		//int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
		int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
		int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
		float fineweightave=0.;
		for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
		  fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
		}
		
		fineweightave*=(ptdechi-ptdeclo);
		sumtrigweight+=fineweightave;
	      }
	      
	      for(int ixidecbin=0;ixidecbin<DECXIDIST_LARGE[iside][0]->GetNbinsX();ixidecbin++){
		float xidec=DECXIDIST_LARGE[iside][itdec]->GetBinCenter(ixidecbin+1);
		float xideclo=xidec-decxibw/2.;
		float xidechi=xidec+decxibw/2.;
		
		//float ptdeclo=fabs(pt*cos(dphifold)*exp(xideclo));
		//float ptdechi=fabs(pt*cos(dphifold)*exp(xidechi));
		float ptdeclo=fabs(pt*exp(xideclo));
		float ptdechi=fabs(pt*exp(xidechi));

		
		if(ptdechi<trigbinbounds_large[itdec]||ptdeclo>trigbinbounds_large[itdec+1]) continue;
		if(ptdeclo<trigbinbounds_large[itdec]) ptdeclo=trigbinbounds_large[itdec];
		if(ptdechi>trigbinbounds_large[itdec+1]) ptdechi=trigbinbounds_large[itdec+1];
		//int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
		int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
		int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
		float fineweightave=0.;
		for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
		  fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
		}
		
		fineweightave*=(ptdechi-ptdeclo);
		if(sumtrigweight>0)fineweightave*=trigweightnorm/sumtrigweight;
		
			
		DECXIDIST_LARGE[iside][itdec]->Fill(xidec,dphifold,filltimeweight*fineweightave);
		DECXIDIST_LARGE[iside+2][itdec]->Fill(xidec,dphifold,filltimeflow*fineweightave);

		DECXIDIST_STATERR[iside][itdec]->Fill(xidec,dphifold,filltimeweighterr*fineweightave);
		DECXIDIST_STATERR[iside+2][itdec]->Fill(xidec,dphifold,filltimeflowerr*fineweightave);

		DECXIDIST_LARGE_NOW[iside][itdec]->Fill(xidec,dphifold,fineweightave);

		DECPT2XI[itdec]->Fill(pt,xidec,fineweightave);


		DECXIDIST_LARGE_PLUS[iside][itdec]->Fill(xidec,dphifold,filltimeup*fineweightave);
		DECXIDIST_LARGE_MINUS[iside][itdec]->Fill(xidec,dphifold,filltimedown*fineweightave);
		DECXIDIST_LARGE_PLUS[iside+2][itdec]->Fill(xidec,dphifold,filltimeflowup*fineweightave);
		DECXIDIST_LARGE_MINUS[iside+2][itdec]->Fill(xidec,dphifold,filltimeflowdown*fineweightave);


	      }

    	  	    
	      sumtrigweight =0.;
	      
	      for(int iztdecbin=0;iztdecbin<DECZTDIST_LARGE[iside][itdec]->GetNbinsX();iztdecbin++){
		float ztdec=DECZTDIST_LARGE[iside][itdec]->GetBinCenter(iztdecbin+1);
		float ztdeclo=ztdec-decztbw/2.;
		float ztdechi=ztdec+decztbw/2.;
		
		//float ptdec=fabs(pt / ztdec);
		float ptdeclo=fabs(pt / ztdechi);
		float ptdechi=fabs(pt / ztdeclo);
		if(ptdechi<trigbinbounds_large[itdec]||ptdeclo>trigbinbounds_large[itdec+1]) continue;
		
		
		
		if(ptdeclo<trigbinbounds_large[itdec]) ptdeclo=trigbinbounds_large[itdec];
		if(ptdechi>trigbinbounds_large[itdec+1]) ptdechi=trigbinbounds_large[itdec+1];
		//int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
		int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
		int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
		float fineweightave=0.;
		for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
		  fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
		}
		
		//take into account the phase space for decay
		fineweightave*=(ptdechi-ptdeclo);
		sumtrigweight+=fineweightave;
	      }
	      
	      
	      for(int iztdecbin=0;iztdecbin<DECZTDIST_LARGE[iside][0]->GetNbinsX();iztdecbin++){
		float ztdec=DECZTDIST_LARGE[iside][itdec]->GetBinCenter(iztdecbin+1);
		float ztdeclo=ztdec-decztbw/2.;
		float ztdechi=ztdec+decztbw/2.;
		
		//need to include the effect of decay kinematics here:	    
		//float ptdec=fabs(pt / ztdec);
		float ptdeclo=fabs(pt / ztdechi);
		float ptdechi=fabs(pt / ztdeclo);
		
		if(ptdechi<trigbinbounds_large[itdec]||ptdeclo>trigbinbounds_large[itdec+1]) continue;
		if(ptdeclo<trigbinbounds_large[itdec]) ptdeclo=trigbinbounds_large[itdec];
		if(ptdechi>trigbinbounds_large[itdec+1]) ptdechi=trigbinbounds_large[itdec+1];
		
		//int ipdecbin = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdec);	    
		int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);	    
		int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);	    
		float fineweightave=0.;
		for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++){
		  fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin,idecbin);
	      }
		
		fineweightave*=(ptdechi-ptdeclo);
		if(sumtrigweight>0)fineweightave*=trigweightnorm/sumtrigweight;
		
		DECZTDIST_LARGE[iside][itdec]->Fill(ztdec,dphifold,filltimeweight*fineweightave);
		DECZTDIST_LARGE[iside+2][itdec]->Fill(ztdec,dphifold,filltimeflow*fineweightave);

		DECZTDIST_STATERR[iside][itdec]->Fill(ztdec,dphifold,filltimeweighterr*fineweightave);
		DECZTDIST_STATERR[iside+2][itdec]->Fill(ztdec,dphifold,filltimeflowerr*fineweightave);

		DECZTDIST_LARGE_NOW[iside][itdec]->Fill(ztdec,dphifold,fineweightave);

		DECPT2ZT[itdec]->Fill(pt,ztdec,fineweightave);

		
		DECZTDIST_LARGE_PLUS[iside][itdec]->Fill(ztdec,dphifold,filltimeup*fineweightave);
		DECZTDIST_LARGE_MINUS[iside][itdec]->Fill(ztdec,dphifold,filltimedown*fineweightave);
		DECZTDIST_LARGE_PLUS[iside+2][itdec]->Fill(ztdec,dphifold,filltimeflowup*fineweightave);
		DECZTDIST_LARGE_MINUS[iside+2][itdec]->Fill(ztdec,dphifold,filltimeflowdown*fineweightave);

		
		DECZTDIST_PLUS[iside][itdec]->Fill(ztdec,dphifold,filltimeup*fineweightave);
		DECZTDIST_MINUS[iside][itdec]->Fill(ztdec,dphifold,filltimedown*fineweightave);
		DECZTDIST_PLUS[iside+2][itdec]->Fill(ztdec,dphifold,filltimeflowup*fineweightave);
		DECZTDIST_MINUS[iside+2][itdec]->Fill(ztdec,dphifold,filltimeflowdown*fineweightave);


	      }
	      
	      	      
	      if(pouthistos){
		if(pt>poutlothresh){
		  DECPOUTDIST_LARGE[iside][itdec]->Fill(pout,dphifold,filltimeweight*mwweight[itdecl]);
		  DECPOUTDIST_LARGE[iside+2][itdec]->Fill(pout,dphifold,filltimeflow*mwweight[itdecl]);

		  DECPOUTDIST_STATERR[iside][itdec]->Fill(pout,dphifold,filltimeweighterr*mwweight[itdecl]);
		  DECPOUTDIST_STATERR[iside+2][itdec]->Fill(pout,dphifold,filltimeflowerr*mwweight[itdecl]);
		  
		  DECPOUTDIST_LARGE_PLUS[iside][itdec]->Fill(pout,dphifold,filltimeup*mwweight[itdecl]);
		  DECPOUTDIST_LARGE_MINUS[iside][itdec]->Fill(pout,dphifold,filltimedown*mwweight[itdecl]);
		  DECPOUTDIST_LARGE_PLUS[iside+2][itdec]->Fill(pout,dphifold,filltimeflowup*mwweight[itdecl]);
		  DECPOUTDIST_LARGE_MINUS[iside+2][itdec]->Fill(pout,dphifold,filltimeflowdown*mwweight[itdecl]);

		}
		
		if(pt>poutlothresh&&pt<pouthithresh){
		  
		  DECPOUT2DIST_LARGE[iside][itdec]->Fill(pout*pout,dphifold,filltimeweight*mwweight[itdecl]);
		  decpout2sum_large[iside][itdec]+=pout*pout*filltimeweight*mwweight[itdecl];
		  decpout2weight_large[iside][itdec]+=filltimeweight*mwweight[itdecl];
		  
		  DECPOUT2DIST_LARGE[iside+2][itdec]->Fill(pout*pout,dphifold,filltimeflow*mwweight[itdecl]);
		  
		  
		  decpout2sum_large[iside+2][itdec]+=pout*pout*filltimeflow*mwweight[itdecl];
		  
		  decpout2weight_large[iside+2][itdec]+=filltimeflow*mwweight[itdecl];
		  
		  
		}
	      }
	    }     	      

	  }
    }//dofilltime flag
    




    return 1;
  }      // if track good (quailty == 33 || ...		

  return 0;
}


float CombinedSimple::FindTrackDistance(TriggerHolder gamma, PartnerHolder partner, float& dtheta, float& dphi)
{
  //right now I look for matches amongst the same tracks that pass into my analysis  -- that's only tracks above 1 GeV for hpdsts
  //for pi0s mixing pi0
 
 
  float xyz[3]={gamma.x,gamma.y,gamma.z};
 
  //use x/r = px/E
  float hitdist = xyz[0]*gamma.E/gamma.px;
  //now use (z-zvertex)/r = pz/E to get zvertex
  float zVtx=xyz[2]-hitdist*gamma.pz/gamma.E;

  //try using the pc3 projection
   float xtrk = partner.vppc3x;      
   float ytrk = partner.vppc3y;
   float ztrk = partner.vppc3z;

 //  float xtrk =0;
//   float ytrk =0;
//   float ztrk =0;
  
  float epx, epy, epz;
  dphi = fabs(atan2(ytrk,xtrk) - atan2(xyz[1],xyz[0]));
  if(dphi>1.57)
    {
      epx=xtrk;
      epy=ytrk;
      epz=ztrk;
      
    }else{
      closestApproachPPP(zVtx, xtrk,ytrk,ztrk,  xyz[0], xyz[1], xyz[2],epx,epy, epz);
    }

  float dx = epx - xyz[0];
  float dy = epy - xyz[1];
  float dz = epz - xyz[2];
 
  //cout<<" xtrk "<<xtrk<<" ytrk "<<ytrk<<" ztrk "<<ztrk<<endl;
  //cout << "zVtx "<<zVtx <<endl;
  //cout<<" epx "<<epx<<" epy "<<epy<<" epz "<<epz<<endl;
  //cout<<" xyz "<<xyz[0]<<"  "<<xyz[1]<<" "<<xyz[2]<<endl;
  //cout<<" dx "<<dx<<" dy "<<dy<<" dz "<<dz<<endl;

  dtheta = fabs(atan2((double)510.0,(double)epz) - atan2((double)510.0,(double)xyz[2]));
  dphi = fabs(atan2(epy,epx) - atan2(xyz[1],xyz[0]));
 

  float dist = (510. / hitdist) * sqrt( dx*dx + dy*dy + dz*dz );  
  return dist;
}




  int CombinedSimple::ResetEvent(PHCompositeNode *topNode)
{
  global = 0;     
  particle = 0;     
  emccluster =0;
  reacplane = 0;
  d_runhdr =0;
  evtsync =0;
  return 0;
}

int CombinedSimple::End(PHCompositeNode *topNode){

  // this could need changed in the future
  cout<<"event_counter "<<event_counter<<endl;
  return 0;
}

bool CombinedSimple::Chi2Cut(emcClusterContent* sngl_emc, float zvertex)
{
  int arm = sngl_emc->arm();
  int sector = sngl_emc->sector();
  int iypos = sngl_emc->iypos();
  int izpos = sngl_emc->izpos();

  int index = arm * 100000 + sector * 10000 + iypos * 100 + izpos;

  if (index < 100000 || index >=120000) // pbsc
    {
      if(sngl_emc->chi2()<3.0) return true;
    }      
  else // pbgl...dispersion cut
    {
      
      float yz_cg[2], disp[2], pbgl_disp_z_y[2];
      
      yz_cg[0] = sngl_emc->zcg();
      yz_cg[1] = sngl_emc->ycg();
      disp[0] = sngl_emc->dispz(); 
      disp[1] = sngl_emc->dispy(); 
      
      CalcPbglDisp(yz_cg[0], // /cg_z
		   yz_cg[1], //cg_y,
		   disp[0],
		   disp[1],
		   pbgl_disp_z_y);
      
      float fDispYCor = pbgl_disp_z_y[1];
      float fDispZCor = pbgl_disp_z_y[0];
      
      double fX = sngl_emc->x();
      double fY = sngl_emc->y();
      double fZ = sngl_emc->z();
      double fBbcZEvt = zvertex;
      
      // taken straight from offline/ana/PbGl/gampi0 
      // ---> Flag 3 : dispersion cut
      // angle dependent dispersion cuts
      // the cut values are taken from PbGl and are not tuned for PbSc
      static const double p1 =  0.270;
      static const double p2 = -0.0145;
      static const double p3 =  0.00218;
      
      int arm = sngl_emc->arm();
      
      // Incident Angle in degrees (0 = perpendicular to sector surface)
      double theta = IncidentAngle(arm, (index % 100000 / 10000), fX,fY,fZ,fBbcZEvt);
      
      // angle dependent dispersion cut
      double dispCut = p1 + p2*theta + p3*theta*theta;
      
      // maximal dispersion
      double dispMax = fDispZCor;
      if (dispMax < fDispYCor) dispMax = fDispYCor;
      
      // apply dispersion cut
      if (dispMax < dispCut) 
	return true;
      // cout<< "dispMax "<<dispMax<< " dispCut "<<dispCut<<endl;
    } //if PbGl
     
  return false;  
}



float CombinedSimple::GetEcore(emcClusterContent *photon)
{  
  float E_tot=0;

  int multi = photon->multiplicity();
  float ecore = photon->ecore();
  //      cout<<" ecore "<<ecore<<endl;
  float e = photon->e();  
  double E[multi]; // To be filled with corrected E in each tower.
  int tsec, ty, tz;
  //  First do for tower 0
  if (photon->towerid(0)<15552){
    //  First 15552 towers are PbSc in section 0,1,2,3,6,7
    tsec=photon->towerid(0)/2592;  // get sector
    if (tsec > 3) {tsec = tsec+2;}  // correct for sector shift 
    ty=(photon->towerid(0)%2592)/72;  //get iy
    tz=(photon->towerid(0)%2592)%72;  //get iz
  }else{
    //  Rest are from PbGl sectors 4,5
    tsec=((photon->towerid(0)-15552)/4608) + 4;  // get sector (+4 to get PbGl)
    ty=((photon->towerid(0)-15552)%4608)/96;  // get iy
    tz=((photon->towerid(0)-15552)%4608)%96;  // get iz
  }
  E[0] = photon->partesum(0)*tower_corr[tsec][ty][tz]*(ecore/e);
  E_tot = E[0];
  //  Second do for other towers
  for (int i=1; i<multi; i++){
    if (photon->towerid(i)<15552){  
      tsec=photon->towerid(i)/2592;  // get sector
      if (tsec > 3) {tsec = tsec+2;}  // correct for sector shift 
      ty=(photon->towerid(i)%2592)/72;  //get iy
      tz=(photon->towerid(i)%2592)%72;  //get iz
    }else{
      tsec=((photon->towerid(i)-15552)/4608) + 4;  // get sector
      ty=((photon->towerid(i)-15552)%4608)/96;  //get iy
      tz=((photon->towerid(i)-15552)%4608)%96;  //get iz
    }      
    E[i] = 
      (photon->partesum(i)-photon->partesum(i-1))*tower_corr[tsec][ty][tz]*(ecore/e);

    E_tot +=E[i];
  }
  
  //  if(ecore>5)cout<<" ecore "<<ecore<<" E_tot "<<E_tot<<endl;
  return E_tot;
}


void CombinedSimple::CalcPbglDisp(const float& m1z, const float& m1y,
				      const float& m2z, const float& m2y,
				            float* returnVals_z_y
				      ) 
{

  double zModSize, yModSize;

  // modules sizes extracted from CVS
  // offline/packages/emc/mEmcGeometryModule.C (revision 2.8)
  // -> mEmcGeometryModule::BuildGeometry()

  //   if ((iSecId == idSecE0) || (iSecId == idSecE1)) {
  zModSize = 4.09168;
  yModSize = 4.10615;
  //   }
  //   else {
  //     zModSize = 5.58003;
  //     yModSize = 5.57091;
  //   }

  //=====> dispersion stuff
  // correct simple periodic structure of the dispersion
 float fPosZ  = m1z / zModSize;
 float fPosY  = m1y / yModSize;
 
 float fDispZ = m2z / (zModSize * zModSize);
 float fDispY = m2y / (yModSize * yModSize);
 
  double zPosMod = fPosZ - floor(fPosZ);
  double yPosMod = fPosY - floor(fPosY);
  
  float fDispZCor = fDispZ - (zPosMod - zPosMod*zPosMod);
  float fDispYCor = fDispY - (yPosMod - yPosMod*yPosMod);

  returnVals_z_y[0] = fDispZCor;
  returnVals_z_y[1] = fDispYCor;

}

double CombinedSimple::IncidentAngle(const int arm,
					 const int sect,
					 const double& fX, 
					 const double& fY,
					 const double& fZ,
					 const double& fBbcZEvt)
{

  double CosPhiNum =   vnx[arm][sect] * fX 
                     + vny[arm][sect] * fY 
    + vnz[arm][sect] * (fZ-fBbcZEvt); 

  double CosPhiDen = sqrt(fX*fX + fY*fY + (fZ-fBbcZEvt)*(fZ-fBbcZEvt));
 
  double Phi;
  if (CosPhiDen != 0) Phi = 180./3.14159 * acos(CosPhiNum / CosPhiDen);
  else                Phi = -999.;

  return Phi;  // in degrees 

}




void CombinedSimple::InitializeBadRunList()
{  

  //jiamin's bad emcal run list  
  sbadruns.push_back(110544);
  sbadruns.push_back(110545);
  sbadruns.push_back(110568);
  sbadruns.push_back(110601);
  sbadruns.push_back(110602);
  sbadruns.push_back(110604);
  sbadruns.push_back(110609);
  sbadruns.push_back(110628);
  sbadruns.push_back(110634);
  sbadruns.push_back(110651);
  sbadruns.push_back(110652);
  sbadruns.push_back(110655);
  sbadruns.push_back(110663);
  sbadruns.push_back(110665);
  sbadruns.push_back(110667);
  sbadruns.push_back(110669);
  sbadruns.push_back(110671);
  sbadruns.push_back(110683);
  sbadruns.push_back(110687);
  sbadruns.push_back(110689);
  sbadruns.push_back(110691);
  sbadruns.push_back(110698);
  sbadruns.push_back(110700);
  sbadruns.push_back(110704);
  sbadruns.push_back(110706);
  sbadruns.push_back(110708);
  sbadruns.push_back(111824);
  sbadruns.push_back(111830);
  sbadruns.push_back(111831);
  sbadruns.push_back(120408);
  sbadruns.push_back(120410);
  sbadruns.push_back(120411);
  sbadruns.push_back(120416);
  sbadruns.push_back(120419);
  sbadruns.push_back(120420);
  sbadruns.push_back(120422);
  sbadruns.push_back(120426);
  sbadruns.push_back(120427);
  sbadruns.push_back(120428);
  sbadruns.push_back(120478);
  sbadruns.push_back(120479);
  sbadruns.push_back(121510);
  sbadruns.push_back(121954);
  sbadruns.push_back(121956);
  sbadruns.push_back(121959);
  sbadruns.push_back(121961);
  sbadruns.push_back(121967);
  sbadruns.push_back(121968);
  sbadruns.push_back(122041);
  sbadruns.push_back(122212);
  sbadruns.push_back(122213);
  sbadruns.push_back(122214);
  sbadruns.push_back(122215);
  sbadruns.push_back(122220);
  sbadruns.push_back(122221);
  sbadruns.push_back(122223);
  

  return;
}



int CombinedSimple::CheckBadRunList(int runnumber)
{
  
  int size = sbadruns.size();  
  for(int i=0;i<size;i++){
    if(runnumber==sbadruns[i]) return 1;   
  }
  
  return 0;
}



float CombinedSimple::FindTrackDistance(TriggerHolder gamma, float xtrk, float ytrk, float ztrk, float& dtheta, float& dphi)
{
  //right now I look for matches amongst the same tracks that pass into my analysis  -- that's only tracks above 1 GeV for hpdsts
 
 
 
  float xyz[3]={gamma.x,gamma.y,gamma.z};
 
  //use x/r = px/E
  float hitdist = xyz[0]*gamma.E/gamma.px;
  //now use (z-zvertex)/r = pz/E to get zvertex
  float zVtx=xyz[2]-hitdist*gamma.pz/gamma.E;

  //try using the pc3 projection
  //float xtrk = track->get_ppc3x();      
  //float ytrk = track->get_ppc3y();
  //float ztrk = track->get_ppc3z();


 //  float xtrk =0;
//   float ytrk =0;
//   float ztrk =0;
  
  float epx, epy, epz;
  dphi = fabs(atan2(ytrk,xtrk) - atan2(xyz[1],xyz[0]));
  if(dphi>1.57)
    {
      epx=xtrk;
      epy=ytrk;
      epz=ztrk;
      
    }else{
      closestApproachPPP(zVtx, xtrk,ytrk,ztrk,  xyz[0], xyz[1], xyz[2],epx,epy, epz);
    }
  float dx = epx - xyz[0];
  float dy = epy - xyz[1];
  float dz = epz - xyz[2];
  
  
  //cout<<" epx "<<epx<<" epy "<<epy<<" epz "<<epz<<endl;
  //cout<<" xyz "<<xyz[0]<<"  "<<xyz[1]<<" "<<xyz[2]<<endl;
  //cout<<" dx "<<dx<<" dy "<<dy<<" dz "<<dz<<endl;
  
  dtheta = fabs(atan2((double)510.0,(double)epz) - atan2((double)510.0,(double)xyz[2]));
  dphi = fabs(atan2(epy,epx) - atan2(xyz[1],xyz[0]));
  
  float dist = (510. / hitdist) * sqrt( dx*dx + dy*dy + dz*dz );  
  return dist;
}


void CombinedSimple::closestApproachPPP(float vtxZ, float p1x, float p1y, float p1z,
						float p2x, float p2y, float p2z,
						float& p3x, float& p3y, float& p3z)
  
{
  //finds the point on the line formed by p1 and vertex closest to p2: 
  // outputs the result to p3
  p1z = p1z - vtxZ; 
  p2z = p2z - vtxZ;
  static float p2costhet_overp1;
  
  p2costhet_overp1 = (p1x*p2x + p1y*p2y + p1z*p2z)/ (p1x*p1x +p1y *p1y + p1z*p1z);
  
  p3x = p2costhet_overp1 * p1x;
  p3y = p2costhet_overp1 * p1y;
  p3z = vtxZ + p2costhet_overp1 * p1z;
  
}



float CombinedSimple::foldphi(float inphi)
{ 
  // this only being used for ph-ph diffs,thus we don't need to worry about
  // -2*pi -2pi = -4pi cases 
  //  inphi += TMath::Pi()/2.0;
  //  if (inphi > 2.0 *TMath::Pi() || inphi < -2.0*TMath::Pi()) cout << "warning!!!" << endl;
  if (inphi < 0) inphi += 2.0*TMath::Pi();
  if (inphi>TMath::Pi()) inphi = 2.0*TMath::Pi() - inphi;
  return inphi;
} 



int CombinedSimple::CheckLvl1Fired(emcClusterContent * emc, ErtOut * Ert, int * ertTrigs)
{
  static smTileModule* smtilemodule = new smTileModule;
  
  int towerkey = 100000*emc->arm() +  10000*emc->sector() +  100*emc->iypos() +  emc->izpos();
  int arm= -1;
  int sector= -1;
  int smid = -1;
  
  if(smtilemodule->get_smID(towerkey, arm, sector, smid)) {
    
    for (int i = 0; i < 3; i++)
      {
	if (ertTrigs[i] && Ert->get_ERTbit(0, arm, sector, smid))
	  return 1;
      }
    
    // 	  Albit4x4a1[i] = Ert->get_ERTbit(0, arm, sector, smid);
    // 	  Albit4x4b1[i] = Ert->get_ERTbit(1, arm, sector, smid);
    // 	  Albit4x4c1[i] = Ert->get_ERTbit(2, arm, sector, smid);
    
  } // get smID worked
  

  return 0;
  
  
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//below is the impl's of  wei's lvl1 util class which we've included as a private (sub) class of ppAnaAwaysideSimple
//
//
//



//
int smTileModule::get_towerkey(int arm, int sector, int iy, int iz) 
{
  return (100000 * arm + 10000 * sector + 100 * iy + iz);
}
//
//... get smID ...
//
bool smTileModule::get_smID(int towerkey, int &arm, int &sector, int &smID)
{
  int moduleID = GetModuleFromTower(towerkey);
  if(moduleID==-1) return false;

  int emctt = GetEmcttFromModule(moduleID);
  if(emctt==-1) return false;

  GetSMidFromEMCTrgTile(arm, sector,smID, emctt);
  
  if(debug) {
    cout<<" moduleID = "<<moduleID<<" emctt = "<<emctt<<endl;
    cout<<" arm = "<<arm<<" sector = "<<sector<<" smID = "<<smID<<endl;
  }

  return true;

}
//
//... steal functions from ert trigger simulation
//
int smTileModule::GetModuleFromTower(int towerkey)
{//.. get module ID from Tower Key

  int iz=towerkey%100;
  towerkey/=100;
  int iy=towerkey%100;
  towerkey/=100;
  int sector=towerkey%10;
  towerkey/=10;
  int arm=towerkey;

  if (0==arm && sector<4 && iy<36 && iz<72)
            return (int) iz/2+(iy/2)*36+sector*648;

  if (1==arm && sector<2 && iy<48 && iz<96) 
            return (int) iz/2+(iy/2)*48+sector*1152+2592;

  if (1==arm && sector<4 && iy<36 && iz<72)
            return (int) iz/2+(iy/2)*36+(sector-2)*648+4896;

  return -1;
}
//........................
int smTileModule::GetEmcttFromModule(int moduleid)
{//.. get trigger tile ID from module id 

  if (moduleid>=0 && moduleid<2592) 
    return (int) (35-moduleid%36)/6+((moduleid/36)/6)*6;

  if (moduleid>=2592 && moduleid<4896) { 
    moduleid-=2592;
    return (int) (moduleid%48)/6-((moduleid/48)/6)*8+164;
  }
      
  if (moduleid>=4896 && moduleid<6192) {
    moduleid-=4896;
    return (int) (moduleid%36)/6-((moduleid/36)/6)*6+102;
  }

  return -1;
}
//....................
void smTileModule::GetSMidFromEMCTrgTile(int &arm, int &sector, int &smID, int emctt)
{ //  convert ERT EMCal simulated trigger tile to
  //  the hardware supermodule position

   if(emctt<72) {//...west arm

      if(emctt>=0  && emctt<18)        sector = 0; 
      else if (emctt>=18 && emctt<36)  sector = 1;
      else if (emctt>=36 && emctt<54)  sector = 2;
      else if (emctt>=54 && emctt<72)  sector = 3;

      arm = 0; //.. west arm 

      smID = ((emctt - 18*sector)/6 + 1)*6 - (emctt - 18*sector)%6 - 1;   

   } else if(emctt<172) { //.. east arm

     if(emctt>=72 && emctt<90)          sector = 3;
     else if(emctt>=90 && emctt<108)    sector = 2;
     else if(emctt>=108 && emctt<140)   sector = 1;
     else if(emctt>=140 && emctt<172)   sector = 0;

     arm = 1; //.. east arm 

     if(sector<2) {
        
         smID = (emctt - 140 + sector*32) + 24 - (emctt - 140 + sector*32)/8*16;

     } else if(sector <4) {

         smID = (emctt - 126 + sector*18) + 12 - (emctt - 126 + sector*18)/6*12;
        
     } else {
        cout<<"+++  sector should not >= 4 , something is wrong +++ "<<endl;
     }

   } else {
      cout<<"++++ emctt = "<<emctt<<" goes beyond the range"<<endl; 
   }
}
//..............................................
void smTileModule::get_CRKsmID(int pmt, int &armRICH, int &sectRICH,int &smRICH)
{
  int crktt=  GetCrkttFromPMT(pmt);
  GetSMidFromCrkTrgTile(armRICH, sectRICH, smRICH, crktt);
}
//..........................................
int smTileModule::GetCrkttFromPMT(int pmt)
{
  int crktt=0;
  int crk_sector = pmt/1280;
  int npmt = pmt%1280;
  
  if (0==crk_sector) {
    crktt=(15-npmt%16)/4+((npmt/16)/5)*8;
  }
  if (1==crk_sector) {
    crktt=(npmt%16)/4+((npmt/16)/5)*8+4;
  }
  if (2==crk_sector) {
    crktt=(15-npmt%16)/4+((79-npmt/16)/5)*8+128;
  }
  if (3==crk_sector) {
    crktt=(npmt%16)/4+((79-npmt/16)/5)*8+132;
  }
  return(crktt);
}
//......................................
void smTileModule::GetSMidFromCrkTrgTile(int &arm, int &sector, int &smID, int crktt)
{
  if(crktt<128) { //.. west arm
    arm = 0;
    sector = crktt/32;
    smID = 8*((2*(crktt/8)+1)-4*sector)-1 - crktt;
  } else if (crktt<256) { //.. east arm
    arm = 1;
    sector = 7 - crktt/32;
    int tmp = crktt - 104 - (3-sector)*32;
    smID = tmp - (tmp/8 - 3)*16;
  } else {
    cout<<"error:  crktt = "<<crktt<<" goes beyond range "<<endl;
  }
}


// todo remove
// int isGoodRun(int runnum)

//todo remove
//void LoadNONInterpCorrs() 

int CombinedSimple::CheckNewBadTowers(int lsec, int lIypos, int lIzpos){
  // > > and I get lsec, lIypos and lIzpos by the following:
  // > >    lIypos    = clus->iypos();
  // > >    lIzpos    = clus->izpos();
  // > >    lsec      = clus->sector();
  //> >    if(!clus->arm()) lsec+=4;

  //me:      0, 1, 2 ,3, 4, 5, 6, 7
  //jiamin:  4, 5, 6, 7, 3, 2, 1, 0

  if( lsec==7&&((lIzpos==6&&lIypos==4)||(lIzpos==7&&lIypos==4)||(lIzpos==54&&lIypos==5)||(lIzpos==55&&lIypos==5)) ) return 1;  //jiamin -- lsec=0
  if( lsec==6&&((lIzpos==86&&lIypos==2)||(lIzpos=87&&lIypos==2)||(lIzpos==54&&lIypos==7)||(lIzpos=55&&lIypos==7)||(lIzpos==4&&lIypos==12)||(lIzpos=5&&lIypos==12)||(lIzpos==10&&lIypos==35)||(lIzpos=11&&lIypos==35) ) ) return 1; //jiamin -- lsec=1
  if( lsec==5&&((lIzpos==70&&lIypos==23)||(lIzpos=71&&lIypos==23)||(lIzpos==54&&lIypos==4)||(lIzpos=55&&lIypos==4)||(lIzpos==0&&lIypos==28)||(lIzpos=1&&lIypos==28) ) ) return 1; //jiamin -- lsec=2
  if( lsec==2&&((lIzpos==0&&lIypos==27)||(lIzpos=1&&lIypos==27)))return 1; //jiamin -- lsec=6
  if( lsec==3&&((lIzpos==8&&lIypos==19)||(lIzpos=9&&lIypos==19)||(lIzpos==18&&lIypos==21)||(lIzpos=19&&lIypos==21)||(lIzpos==48&&lIypos==21)||(lIzpos=49&&lIypos==21) ) ) return 1; //jiamin -- lsec=7
  return 0;
}

int CombinedSimple::twr2id(int armsect, int ytwr,int ztwr)
{
  int itw,id;
  if(armsect<6){ // PbSc
    if(ytwr<0||ytwr>35) return -1;
    if(ztwr<0||ztwr>71) return -1;
    itw=ytwr*72+ztwr;
    id=armsect*2592+itw;
  }else{
    if(ytwr<0||ytwr>47) return -1;
    if(ztwr<0||ztwr>95) return -1;
    itw=ytwr*96+ztwr;
    id=15552+(armsect-6)*4608+itw;
  }
  

  return id;
} 



void CombinedSimple::SetTagFlag(int flag){
  tagflag=flag;
  return;
}


void CombinedSimple::SetSpeciesFlag(int flag){

  species=flag;
  return;
}



void CombinedSimple::InitializeHotTowerMap(int species_ht, int addTowers){

  ifstream fhot;
  TOAD *toad_hottower = new TOAD("combinesimple");

  if (!_emcBadMapFilename)
    cout << "CombinedSimple:: init: No hot tower map loaded" << endl;
  else{
    string hot_file = toad_hottower->location(_emcBadMapFilename->Data());
    fhot.open(hot_file.c_str());
    cout << PHWHERE << "loading hottower map " << hot_file << endl;
    if (!fhot.is_open())
      cout<<"No Hot File!"<<endl;

    int t1;
    while(!fhot.eof()){
      fhot>>t1;
      HotTowers.push_back(t1);
    }
  }
  cout << PHWHERE << " masking " << HotTowers.size() << " hot towers!" << endl;
  delete toad_hottower;

  if(addTowers) {
    //Added some additional hot towers from Matt 2/26
    HotTowers.push_back(10088);
    HotTowers.push_back(10160);
    HotTowers.push_back(10017);
    HotTowers.push_back(10089);
    HotTowers.push_back(10161);
    HotTowers.push_back(10092);
    HotTowers.push_back(10164);
    HotTowers.push_back(10093);
    HotTowers.push_back(10165);
    HotTowers.push_back(9591 );
    HotTowers.push_back(10023);
    HotTowers.push_back(11722);
    HotTowers.push_back(14724);
    HotTowers.push_back(14872);
    HotTowers.push_back(14944);
    HotTowers.push_back(15016);
    HotTowers.push_back(15160);
    HotTowers.push_back(15376);
    HotTowers.push_back(15018);
    HotTowers.push_back(15090);
    HotTowers.push_back(15091);
    HotTowers.push_back(14732);
    HotTowers.push_back(15021);
    HotTowers.push_back(15381);
    HotTowers.push_back(14735);
    HotTowers.push_back(17024);
    HotTowers.push_back(20764);
    HotTowers.push_back(14783);
    HotTowers.push_back(15019);
  }

  return;
}

int CombinedSimple::CheckHotTowerMap(int sect, int iz, int iy, int do3x3){
  
  int ishot=0;
  
  int imin=0;
  int imax=1;
  int jmin=0;
  int jmax=1;

  if(do3x3){
    imin=-1;
    imax=2;
    jmin=-1;
    jmax=2;
  }

  for(int i=imin;i<imax;i++)
  {
    for(int j=jmin;j<jmax;j++)
    {
      int towernum;
      if(sect<6) towernum=sect*2592+(iz+i)+72*(iy+j);
      else  towernum=15552+(sect-6)*4608+(iz+i)+96*(iy+j);

      int nhottowers=0;
      nhottowers=HotTowers.size();

      for(int k=0;k<nhottowers;k++)
      {
	if(towernum==HotTowers[k])
	{
	  //cout<<" towernum "<<towernum<<" sect "<<sect<<" iz "<<iz<<" iy "<<iy<<" i "<<i<<" j "<<j<<endl;
	  ishot=1;
	  break;
	}
      }
    }
  }
  return ishot;
}

void CombinedSimple::SetTriggerEfficiency(const char* filename)
{
  fpi0eff=new TFile(filename);
  cout << "using pi0 efficiency from file " << fpi0eff->GetName() << endl;

  if(species == 0 || species == 3 || species == 4)//use pp for dAu for now
    fpi0eff->GetObject("ratio_graph_pp",grpi0eff);

  if(species == 1 ||species ==2 || species == 5 )
  {
    int icentbin=0;
    if(locent==20) icentbin=1;
    if(locent==40) icentbin=2;
    if(locent==60) icentbin=3;
    if(locent==80) icentbin=3;

    char grname[100];
    sprintf(grname,"ratio_graph_AA_%d",icentbin);
    fpi0eff->GetObject(grname,grpi0eff);
  }

  return;
}

void CombinedSimple::SetSharkFin(const char* filename)
{
  TFile *fshark_exodus=NULL;
  fshark_exodus = new TFile(filename);

  //now take 10 cm projections out to 165 cm
  //right now use bins out to 170 since we have two cm bins instead of 1.
  cout<<PHWHERE<<" loading sharkfins "<< fshark_exodus->GetName() <<endl;

  for(int izemc=0;izemc<33;izemc++){
    char sharkname[100];
    sprintf(sharkname,"PTpi_PTgam_%d",izemc);
    ptpivsptgam[izemc]=(TH2D*)fshark_exodus->Get(sharkname);
    
    for(int idecl=0;idecl<5;idecl++){
      sprintf(sharkname,"hshark_large_%d_%d",idecl,izemc);
      hshark_large[idecl][izemc]=(TH1D*)fshark_exodus->Get(sharkname);
    }
    for(int idecs=0;idecs<7;idecs++){
      sprintf(sharkname,"hshark_small_%d_%d",idecs,izemc);
      hshark_small[idecs][izemc]=(TH1D*)fshark_exodus->Get(sharkname);
    }
    for(int idecs=0;idecs<7;idecs++){
      sprintf(sharkname,"hshark_alt_%d_%d",idecs,izemc);
      hshark_alt[idecs][izemc]=(TH1D*)fshark_exodus->Get(sharkname);
    }
  }

  if(!removetags)
  {
    cout << "applying pisacorr " << endl;
    TF1 * pisacorr_large[5];
    pisacorr_large[0]=new TF1("pisacorr_large0","-7.85e-2*x+1.85",4,20);
    pisacorr_large[1]=new TF1("pisacorr_large1","-8.7e-2*x+2.0",4,20);
    pisacorr_large[2]=new TF1("pisacorr_large2","-9.65e-2*x+2.2",4,20);
    pisacorr_large[3]=new TF1("pisacorr_large3","-12.5e-2*x+2.8",4,20);
    pisacorr_large[4]=new TF1("pisacorr_large4","-12.5e-2*x+2.8",4,20);
    
    for (int ipzemc = 0; ipzemc < 33; ipzemc++){
      for (int ipdecs = 0; ipdecs < 5; ipdecs++){
	for (int ipbin = 1; ipbin < hshark_large[ipdecs][ipzemc]->GetNbinsX()+1; ipbin++)
	{
	  float mattshark = hshark_large[ipdecs][ipzemc]->GetBinContent(ipbin);
	  float pi0pt =  hshark_large[ipdecs][ipzemc]->GetBinCenter(ipbin);
	  float pisacorr=pisacorr_large[ipdecs]->Eval(pi0pt);
	  if(pisacorr>1.0) pisacorr=1.0;
	  float mwweightfine=mattshark*pisacorr;
	  hshark_large[ipdecs][ipzemc]->SetBinContent(ipbin,mwweightfine);
	}
      }
    }
    for(int ij=0; ij<5; ij++)
      delete pisacorr_large[ij];
  }

  return;
}

void CombinedSimple::SetFilltimeCorrs(const char* filltime_name, const char* filltimeup_name, const char* filltimedown_name)
{

  ffilltimecorrs = new TFile(filltime_name);
  ffilltimesysup = new TFile(filltimeup_name);
  ffilltimesysdown = new TFile(filltimedown_name);

  cout<<PHWHERE<<" loading filltime weights "<< ffilltimecorrs->GetName() <<endl;
  cout<<PHWHERE<<" loading filltime sysup "<< ffilltimesysup->GetName() <<endl;
  cout<<PHWHERE<<" loading filltime sysdown "<< ffilltimesysdown->GetName() <<endl;

  fexemb = new TF1("fexemb","[0]+[1]*exp([2]*x)",5.0,10.0);
  fexemb->SetParameters(0.761,1.640,-4.734);

  feff = (TF1*)ffilltimecorrs->Get("feff");

  char normname[100], decnormname[100];
  char flowname[100], decflowname[100];
  
  for(int it=0;it<4;it++){
    if(tagflag==0){
      if(useiso)
	sprintf(normname,"PTDPHI_NORM_INC_ISO_%d",it);
      else
	sprintf(normname,"PTDPHI_NORM_INC_%d",it);

      sprintf(flowname,"PTDPHI_FLOW_INC_%d",it);
    }
    else if(tagflag==1||tagflag==5){
      if(useiso)
	sprintf(normname,"PTDPHI_NORM_PI0_ISO_%d",it);
      else
	sprintf(normname,"PTDPHI_NORM_PI0_%d",it);
      sprintf(flowname,"PTDPHI_FLOW_PI0_%d",it);

      if(useiso)
	sprintf(decnormname,"PTDPHI_NORM_DECPI0_ISO_%d",it);
      else
	sprintf(decnormname,"PTDPHI_NORM_DEC_%d",it);
      sprintf(decflowname,"PTDPHI_FLOW_DEC_%d",it);
    }
    else if(tagflag==2){
      if(useiso)
	sprintf(normname,"PTDPHI_NORM_ETA_ISO_%d",it);
      else
	sprintf(normname,"PTDPHI_NORM_ETA_%d",it);
      if(useiso)
	sprintf(decnormname,"PTDPHI_NORM_DECETA_ISO_%d",it);
      else
	sprintf(decnormname,"PTDPHI_NORM_DECETA_%d",it);
    }
    else if(tagflag==3){
      if(useiso)
	sprintf(normname,"PTDPHI_NORM_SID_ISO_%d",it);
      else
	sprintf(normname,"PTDPHI_NORM_SID_%d",it);
      if(useiso)
	sprintf(decnormname,"PTDPHI_NORM_DECSID_ISO_%d",it);
      else
	sprintf(decnormname,"PTDPHI_NORM_DECSID_%d",it);
    }

    DPHIPT_NORM[it]=(TH2F*)ffilltimecorrs->Get(normname);
    if(tagflag>0){
      DPHIPT_DECNORM[it]=(TH2F *)ffilltimecorrs->Get(decnormname);
      cout << "getting " << decnormname <<endl;
    }
    DPHIPT_FLOW[it]=(TH2F*)ffilltimecorrs->Get(flowname);
    if(tagflag>0)
      DPHIPT_DECFLOW[it]=(TH2F *)ffilltimecorrs->Get(decflowname);

    DPHIPT_NORM_PLUS[it]=(TH2F*)ffilltimesysup->Get(normname);
    if(tagflag>0)
      DPHIPT_DECNORM_PLUS[it]=(TH2F *)ffilltimesysup->Get(decnormname);

    DPHIPT_FLOW_PLUS[it]=(TH2F*)ffilltimesysup->Get(flowname);
    if(tagflag>0)
      DPHIPT_DECFLOW_PLUS[it]=(TH2F *)ffilltimesysup->Get(decflowname);

    //sysdown
    DPHIPT_NORM_MINUS[it]=(TH2F*)ffilltimesysdown->Get(normname);
    if(tagflag>0)
      DPHIPT_DECNORM_MINUS[it]=(TH2F *)ffilltimesysdown->Get(decnormname);

    DPHIPT_FLOW_MINUS[it]=(TH2F*)ffilltimesysdown->Get(flowname);
    if(tagflag>0)
      DPHIPT_DECFLOW_MINUS[it]=(TH2F *)ffilltimesysdown->Get(decflowname);
 }

  return;
}

int CombinedSimple::isEdgeTower(int sector, int iz, int iy){
  if(sector<6){
    if(iz>=70||iy>=34||iy<2||iz<2) return 1;
  }
  else{
    if(iz>=94||iy>=46||iy<2||iz<2) return 1;
  }

  return 0;
}


void CombinedSimple::makeMNSingles(PHCentralTrack * particle, int doNtup)
{
  
  if (_assocNum >= 0) return;  // if < 0, then never called, see below; 
                               //_assocNum is more of a flag to say 	    
                               // if this is first time makeMN was called 
  as_N = 0;
  for(int Nt=0;Nt<asNt;Nt++){    
    as_px[Nt]=0;     
    as_py[Nt]=0;     
    as_pz[Nt]=0;     
    as_ppc3x[Nt]=0;  
    as_ppc3y[Nt]=0;  
    as_ppc3z[Nt]=0;     
    as_charge[Nt]=0; 
    as_m2emc[Nt]=0;  
    as_m2tof[Nt]=0; 
    // as_quality[Nt]=0;
    as_n0[Nt]=0;
  }

  //remove && doNtup 
  if (_savAssoc && doNtup)
  {
    if (!_assocNtuple)
    {    
      _assocNtuple = new TTree("singles", "single particle info");
      _assocNtuple->SetDirectory(0);

      // event branches...
      _assocNtuple->Branch("runno",    &as_runno,     "runno/I");      
      _assocNtuple->Branch("segno",    &as_segno,     "segno/I");
      _assocNtuple->Branch("evtno",    &as_evtno,     "evtno/L");
      _assocNtuple->Branch("zvertex",  &as_zvertex,   "zvertex/F");          
      _assocNtuple->Branch("N",        &as_N,         "N/I");          

      if(species>0) _assocNtuple->Branch("centrality",   &as_centrality,    "centrality/I");

      //  track branches...arrays of length N
      _assocNtuple->Branch("px",       as_px,         "px[N]/F");
      _assocNtuple->Branch("py",       as_py,         "py[N]/F");
      _assocNtuple->Branch("pz",       as_pz,         "pz[N]/F");
      _assocNtuple->Branch("charge",   as_charge,     "charge[N]/I");
      _assocNtuple->Branch("m2emc",    as_m2emc,      "m2emc[N]/F");
      _assocNtuple->Branch("m2tof",    as_m2tof,      "m2tof[N]/F");
      _assocNtuple->Branch("ppc3x",    as_ppc3x,      "ppc3x[N]/F");
      _assocNtuple->Branch("ppc3y",    as_ppc3y,      "ppc3y[N]/F");
      _assocNtuple->Branch("ppc3z",    as_ppc3z,      "ppc3z[N]/F");
      _assocNtuple->Branch("n0",       as_n0,         "n0[N]/I");  

    }
  }


  float npart = particle->get_npart();

  for(int i=0;i<npart;i++)
  {      
    //PHSnglCentralTrack* trk = particle->get_track(i); 

    if (particle->get_phi0(i) <= -99) continue;

    Float_t px_temp = particle->get_px(i);
    Float_t py_temp = particle->get_py(i);

    as_pt[as_N] = sqrt(px_temp*px_temp+py_temp*py_temp); 
    //changed lower pT to 0.5 instead of 1GeV 
    //changed upper pT to 7.0 for taxi 218     
    PT_N0->Fill(as_pt[as_N],particle->get_n0(i));

    Int_t quality_temp=particle->get_quality(i);
    int pass_quality = 0;
    if(quality_temp > 7) pass_quality = 1;
    if(quality_temp==31 || quality_temp==63) pass_quality = 2;
    PT_QUAL->Fill(as_pt[as_N],pass_quality);

    if(as_pt[as_N] < 0.5 || as_pt[as_N] > 7.0) continue;
    HADRONCUTCHECK->Fill(1, as_pt[as_N], quality_temp==63||quality_temp==31, as_pt[as_N]<5.0&&particle->get_n0(i)>=0, particle->get_pc3sdphi(i), particle->get_pc3sdz(i));//add by hge

    float pc3sdphi=particle->get_pc3sdphi(i);
    float pc3sdz=particle->get_pc3sdz(i);
    float pc3dphi=particle->get_pc3dphi(i);
    float pc3dz=particle->get_pc3dz(i);
    if ( pc3sdphi>-9999 && (sqrt(pc3sdz*pc3sdz+pc3sdphi*pc3sdphi)<2.0) ) PT_PC3SDPHISDZ->Fill(as_pt[as_N]);

    if(as_pt[as_N]<5.0&&particle->get_n0(i)>=0) continue;
    DCH_N0->Fill(as_pt[as_N],particle->get_phi(i),particle->get_zed(i));
    PT_AFTN0->Fill(as_pt[as_N]);

    if(!(quality_temp == 31 || quality_temp == 63)) continue;
    DCH_N0_QUAL->Fill(as_pt[as_N],particle->get_phi(i),particle->get_zed(i));
    PT_AFTN0QUAL->Fill(as_pt[as_N]);
 
    
    //if(quality_temp<=7)continue;
    // float pc3sdphi=particle->get_pc3sdphi(i);
    // float pc3sdz=particle->get_pc3sdz(i);
    // float pc3dphi=particle->get_pc3dphi(i);
    // float pc3dz=particle->get_pc3dz(i);
    // float emcsdphi=particle->get_emcsdphi(i);
    // float emcsdz=particle->get_emcsdz(i);
    //float emcdphi=particle->get_emcdphi(i);
    //float emcdz=particle->get_emcdz(i);
    //float charge=particle->get_charge(i);

    if(species==5||species==1)
    {
      PT_PC3SDPHI_PC3SDZ->Fill(as_pt[as_N],pc3sdphi,pc3sdz);
    }
    if(pc3sdphi == -9999){
      PT_PC3SDPHI_UNDEFINE->Fill(as_pt[as_N]);
      if(pc3dphi == -9999){
	PT_PC3DPHI_UNDEFINE->Fill(as_pt[as_N]);
      }
    }
    if(pc3sdz == -9999){
      PT_PC3SDZ_UNDEFINE->Fill(as_pt[as_N]);
      if(pc3dz == -9999){
	PT_PC3DZ_UNDEFINE->Fill(as_pt[as_N]);
      }
    }
    //Changed to tighter pc3 sigma cut
    //if (!(fabs(pc3sdphi)<3.0&&fabs(pc3sdz)<3.0&&fabs(emcsdphi)<3.0&&fabs(emcsdz)<3.0))
    //if (!_anticuts&&!(fabs(pc3sdphi)<2.0&&fabs(pc3sdz)<2.0&&fabs(emcsdphi)<2.0&&fabs(emcsdz)<2.0))	continue;
    //Circular cut in pc3 and no emc cut
    //cout << "about to make matching cut" <<endl;
    if ( /*pc3sdphi>-9999 && !_anticuts &&*/ !(sqrt(pc3sdz*pc3sdz+pc3sdphi*pc3sdphi)<2.0) ) continue;
    //else if( (species==5||species==1) && !_anticuts && !(fabs(pc3dz+1)<2 && fabs(pc3dphi)<0.006) ) continue;

    if( verbosity ) {
      std::cout<<"pc3sdphi = "<<pc3sdphi<<" pc3dphi = "<<pc3dphi<<std::endl;
      std::cout<<"pc3sdz = "<<pc3sdz<<" pc3dz = "<<pc3dz<<std::endl;
    }
    
    //if (_anticuts==1&&!(fabs(pc3sdphi)>3.0&&fabs(pc3sdz)>3.0&&fabs(emcsdphi)>3.0&&fabs(emcsdz)>3.0&&fabs(pc3sdphi)<30.0&&fabs(pc3sdz)<30.0&&fabs(emcsdphi)<30.0&&fabs(emcsdz)<30.0)) continue;

    //if (_anticuts==2&&!(fabs(pc3sdphi)>1.0&&fabs(pc3sdz)>1.0&&fabs(emcsdphi)>1.0&&fabs(emcsdz)>1.0&&fabs(pc3sdphi)<4.0&&fabs(pc3sdz)<4.0&&fabs(emcsdphi)<4.0&&fabs(emcsdz)<4.0))
      //continue;
      //if (_anticuts==2&&!(sqrt(pc3sdz*pc3sdz+pc3sdphi*pc3sdphi)>2.0 && sqrt(pc3sdz*pc3sdz+pc3sdphi*pc3sdphi)<3.0))  continue;

    //cout << "passed matching cut" <<endl;
 
    PTAFTCUTS->Fill(as_pt[as_N]);
    DCH_N0_QUAL_PC3->Fill(as_pt[as_N],particle->get_phi(i),particle->get_zed(i));
    ptvscent_part->Fill(as_pt[as_N],as_centrality);

    as_px[as_N] = px_temp;
    as_py[as_N] = py_temp;
    as_pz[as_N] = particle->get_pz(i);
    as_charge[as_N] = particle->get_charge(i) * quality_temp;
    as_m2emc[as_N] = particle->get_m2emc(i);
    as_m2tof[as_N] = particle->get_m2tof(i);
    as_ppc3x[as_N] = particle->get_ppc3x(i);
    as_ppc3y[as_N] = particle->get_ppc3y(i);
    as_ppc3z[as_N] = particle->get_ppc3z(i);
    as_n0[as_N] = particle->get_n0(i);
    as_N++;
  }

  //Should require &&doNtup if that is what is required above 
  if (_savAssoc)
  {
    //       as_runno=runno;  
    //       as_segno=seqno;  
    //       as_evtno=evtno;  
    //       as_zvertex=zvertex;   
    //      as_centrality= (int) percent;
    //all filled in process event
    float fdoNtup = (float) doNtup;      
    float oneOrNeg = (0.5-fdoNtup)/-0.5;

    as_centrality *= (int)oneOrNeg;

    _assocNtuple->Fill();

    as_centrality *= (int)oneOrNeg; // oneOrNeg^2 returns val back to normal
    // as_centrality will serve as a "trigger bit" : positive cent values 
    // will be for "minbias events", negative for trigger events (this will allow
    // us to run the fg and the bg from ntuples only)
  }

  _assocNum = as_N;
  
}


float CombinedSimple::SumEcorePtInCone(float thetrig,float phitrig, float zvertex, int i1, int i2, float trigpt, float trigE){

  const static unsigned int sccut3x3Map = 0xffe1ce70;    
  const static unsigned int glcut3x3Map = 0x1ce70;     
  
  unsigned int cut3x3Map;

  double sum=0;
  double sump[20]={0.};
  double sume[20]={0.};
  double cone[20]={0.};
  int npart = particle->get_npart();
  int ncluster = emccluster->size();  

  for(int i=0;i<npart;i++)
  {

    PHSnglCentralTrack* trk = particle->get_track(i); 
    int quality=trk->get_quality();
    float mom = trk->get_mom();
    float the = trk->get_the0();
    float phi = trk->get_phi0();
    phi=PHAngle(phi);

    for(int b=0; b<20; b++)
      cone[b]=1.0/20.0*(b+1);

    float Rdist=sqrt((thetrig-the)*(thetrig-the)+(phitrig-phi)*(phitrig-phi));

    if(quality>3&&mom>0.2&&mom<15)
    {
      if(trigpt<7) EvsR1->Fill(Rdist,mom/trigE);
      if(trigpt<9) EvsR2->Fill(Rdist,mom/trigE);
      if(trigpt<12) EvsR3->Fill(Rdist,mom/trigE);
      if(trigpt<15) EvsR4->Fill(Rdist,mom/trigE);
    }

    for(int j=0; j<20; j++)
    {
      if(Rdist<cone[j])
	if(quality>3&&mom>0.2&&mom<15) sump[j]+=mom;
    }
  }


  for( int j=0; j<ncluster; j++) 
  {            
    if(j==i1||j==i2) continue;
    photon = emccluster->getCluster(j);             

    float ecore=photon->ecore();
    if(ecore<0.5) continue;

    int sector = photon->sector();
    //arm = 0 -> west
    int arm = photon->arm();	    
    if (arm==1) sector = 7 - sector;    
    if (sector < 6) cut3x3Map = sccut3x3Map;
    else cut3x3Map = glcut3x3Map;         

    if(photon->prob_photon()<0.02) continue;
    /*
    int iypos = photon->iypos();
    int izpos = photon->izpos();
    if(species==0){
      if(run_5_6==0) if( hotdead[twr2id(sector,iypos,izpos)]!=1) continue;

      if(run_5_6==1)
	if(CheckHotTowerMap(sector, izpos, iypos, 1)) continue;
      }
    */

    if((photon->warnmap() & cut3x3Map) ==0 && (photon->deadmap() & cut3x3Map) ==0)
    {
      float x = photon->x();
      float y = photon->y();
      float z = photon->z();
      float radius = sqrt( x*x + y*y + (z-zvertex)*(z-zvertex) );	 
      float px = ecore * x / radius;
      float py = ecore * y / radius;
      float pz = ecore * (z-zvertex) / radius;

      float phi = atan2(py,px);
      phi = PHAngle(phi);
      float the=acos(pz/sqrt(px*px+py*py+pz*pz));

      //charge veto --  This is different than what Kensuke and Rob do -- need to look into this
      //int isveto= 0;

      /*
      for(int i=0;i<npart;i++)
      {
	  PHSnglCentralTrack* trk = particle->get_track(i); 
	  int quality=trk->get_quality();

	  if(quality<=3) continue;

	  float trkpx = trk->get_px();
	  float trkpy = trk->get_py();
	  float trkpt = sqrt(trkpx*trkpx+trkpy*trkpy);      

	  float dist = 99;	     

	  if(quality>3){
	    //stay away from recal edge
	    if(trkpt>1.0)
	    {
	      TriggerHolder sumtrig; 			 
	      sumtrig.E=ecore;
	      sumtrig.px=px;
	      sumtrig.py=py;
	      sumtrig.pz=pz;		    
	      sumtrig.x=x;		    
	      sumtrig.y=y;		    
	      sumtrig.z=z;		    
	      float phidist=1.0;
	      float thetadist=1.0;

	      dist=FindTrackDistance(sumtrig,trk,thetadist,phidist);	     
	      if(fabs(dist)<8)
	      {
		//cout<<" sumtrig.x "<<sumtrig.x<<" sumtrig.px "<<sumtrig.px<<endl;
		//cout<<" vetoed photon "<<j<<" with dist "<<dist<<endl;
		isveto=1;
		break;
	      }
	    }
	  }
	}
      */
      //if(isveto) continue;

      float Rdist=sqrt((thetrig-the)*(thetrig-the)+(phitrig-phi)*(phitrig-phi));
      //EvsR->Fill(Rdist,ecore/trigE);
      if(trigpt<7) EvsR1->Fill(Rdist,ecore/trigE);
      if(trigpt<9) EvsR2->Fill(Rdist,ecore/trigE);
      if(trigpt<12) EvsR3->Fill(Rdist,ecore/trigE);
      if(trigpt<15) EvsR4->Fill(Rdist,ecore/trigE);

      for(int jj=0; jj<20; jj++)
      {
	if(sqrt((thetrig-the)*(thetrig-the)+(phitrig-phi)*(phitrig-phi))<cone[jj])
	  sume[jj]+=ecore;
      }
    }
  }

  for(int j=0; j<20; j++)
  {
    //float percentradpt = (sump[j]+sume[j]-trigE)/trigE*100.0;
    float percentradpt = (sump[j]+sume[j])/trigE*100.0; //don't subtract trigE

    if(percentradpt>499.99) percentradpt=499.99;
    if(trigpt<7) ISOMAP1->Fill(cone[j]-0.025,percentradpt);
    else if(trigpt<9) ISOMAP2->Fill(cone[j]-0.025,percentradpt);
    else if(trigpt<12) ISOMAP3->Fill(cone[j]-0.025,percentradpt);
    else if(trigpt<15) ISOMAP4->Fill(cone[j]-0.025,percentradpt);
    if(trigpt>25)
    {
      percentradpt = (sump[j]+sume[j])/trigE*100.0; //don't subtract trigE
      if(trigpt<27) ISOMAPSH1->Fill(cone[j]-0.025,percentradpt);
      else if(trigpt<29) ISOMAPSH2->Fill(cone[j]-0.025,percentradpt);
      else if(trigpt<32) ISOMAPSH3->Fill(cone[j]-0.025,percentradpt);
      else if(trigpt<35) ISOMAPSH4->Fill(cone[j]-0.025,percentradpt);
    }
  }

  sum=sume[0]+sump[0];

  return sum;
}

bool CombinedSimple::isHotTower(int towerid, int species)
{
  vector <int> hottower;

  if(species == 3)//run8 dAu MB
  {
    //5 sigma
    hottower.push_back(1830);
    hottower.push_back(1851);
    hottower.push_back(1922);
    hottower.push_back(1994); 
    hottower.push_back(2138);
    hottower.push_back(2139);
    hottower.push_back(2211);
    hottower.push_back(2282);
    hottower.push_back(2355);
    hottower.push_back(2426);
    hottower.push_back(2427);
    hottower.push_back(3737);
    hottower.push_back(7566);
    hottower.push_back(10017);
    hottower.push_back(10088);
    hottower.push_back(10089);
    hottower.push_back(10092);
    hottower.push_back(10093);
    hottower.push_back(10160);
    hottower.push_back(10161);
    hottower.push_back(10164);
    hottower.push_back(10165);
    hottower.push_back(10232);
    hottower.push_back(10233);
    hottower.push_back(11388);
    hottower.push_back(12427);
    hottower.push_back(12498);
    hottower.push_back(14314);
    hottower.push_back(16963);
    hottower.push_back(19419);
    hottower.push_back(19784);
    hottower.push_back(21155);
    hottower.push_back(21258);
    hottower.push_back(21768);
    hottower.push_back(22607);
    hottower.push_back(24286);

    hottower.push_back(1759);
    hottower.push_back(1995);
    hottower.push_back(2190);
    hottower.push_back(2210);
    hottower.push_back(2354);
    hottower.push_back(2498);
    hottower.push_back(7207);
    hottower.push_back(7278);
    hottower.push_back(7423);
    hottower.push_back(8916);
    hottower.push_back(9204);
    hottower.push_back(9205);
    hottower.push_back(9285);
    hottower.push_back(9356);
    hottower.push_back(12132);
    hottower.push_back(12281);
    hottower.push_back(12785);
    hottower.push_back(15381);
    hottower.push_back(19489);
    hottower.push_back(24382);

    hottower.push_back(1923);
    hottower.push_back(5690);
    hottower.push_back(6918);
    hottower.push_back(7062);
    hottower.push_back(7495);
    hottower.push_back(8917);
    hottower.push_back(9336);
    hottower.push_back(9767);
    hottower.push_back(10016);
    hottower.push_back(10148);
    hottower.push_back(10149);
    hottower.push_back(10198);
    hottower.push_back(10199);
    hottower.push_back(10270);
    hottower.push_back(10271);
    hottower.push_back(11853);
    hottower.push_back(12714);
    hottower.push_back(16087);
    hottower.push_back(17709);
    hottower.push_back(21156);
    hottower.push_back(21284);
    hottower.push_back(21730);
    hottower.push_back(21794);
    hottower.push_back(21860);
    hottower.push_back(22485);
    hottower.push_back(24112);
    hottower.push_back(24163);
    hottower.push_back(24627);

    hottower.push_back(1758);
    hottower.push_back(1974);
    hottower.push_back(6919);
    hottower.push_back(7135);
    hottower.push_back(7206);
    hottower.push_back(7279);
    hottower.push_back(7350);
    hottower.push_back(9932);
    hottower.push_back(10018);
    hottower.push_back(10021);
    hottower.push_back(10090);
    hottower.push_back(10094);
    hottower.push_back(10162);
    hottower.push_back(10166);
    hottower.push_back(10234);
    hottower.push_back(10292);

    hottower.push_back(6991);
    hottower.push_back(7063);
    hottower.push_back(7134);
    hottower.push_back(7494);
    hottower.push_back(7567);
    hottower.push_back(7983);
    hottower.push_back(9945);
    hottower.push_back(10076);
    hottower.push_back(10237);

    hottower.push_back(7422);

    hottower.push_back(1778);
    hottower.push_back(3832);
    hottower.push_back(6675);
    hottower.push_back(7450);
    hottower.push_back(9164);
    hottower.push_back(15445);
    hottower.push_back(19595);

    hottower.push_back(1975);
    hottower.push_back(4381);
    hottower.push_back(8493);
    hottower.push_back(9856);
    hottower.push_back(23428);

    /*
    //10sigma
    hottower.push_back(1830);
    hottower.push_back(1851);
    hottower.push_back(2138);
    hottower.push_back(2282);
    hottower.push_back(2426);
    hottower.push_back(7566);
    hottower.push_back(10088);
    hottower.push_back(10089);
    hottower.push_back(10092);
    hottower.push_back(10093);
    hottower.push_back(10160);
    hottower.push_back(10161);
    hottower.push_back(10164);
    hottower.push_back(10165);
    hottower.push_back(10232);
    hottower.push_back(10233);
    hottower.push_back(11388);
    hottower.push_back(12427);
    hottower.push_back(12498);
    hottower.push_back(14314);
    hottower.push_back(16963);
    hottower.push_back(19419);
    hottower.push_back(19784);
    hottower.push_back(21258);
    hottower.push_back(21768);
    hottower.push_back(22607);
    hottower.push_back(24286);

    hottower.push_back(1830);
    hottower.push_back(1851);
    hottower.push_back(2138);
    hottower.push_back(2282);
    hottower.push_back(2426);
    hottower.push_back(7566);
    hottower.push_back(10088);
    hottower.push_back(10089);
    hottower.push_back(10092);
    hottower.push_back(10093);
    hottower.push_back(10160);
    hottower.push_back(10161);
    hottower.push_back(10164);
    hottower.push_back(10165);
    hottower.push_back(10232);
    hottower.push_back(10233);
    hottower.push_back(11388);
    hottower.push_back(12427);
    hottower.push_back(12498);
    hottower.push_back(14314);
    hottower.push_back(16963);
    hottower.push_back(19419);
    hottower.push_back(19784);
    hottower.push_back(21258);
    hottower.push_back(21768);
    hottower.push_back(22607);
    hottower.push_back(24286);

    hottower.push_back(1922);
    hottower.push_back(1994);
    hottower.push_back(2139);
    hottower.push_back(2211);
    hottower.push_back(2355);
    hottower.push_back(2427);
    hottower.push_back(7207);
    hottower.push_back(7278);
    hottower.push_back(7423);
    hottower.push_back(10017);
    hottower.push_back(12132);
    hottower.push_back(12281);
    hottower.push_back(12785);
    hottower.push_back(19489);
    hottower.push_back(21155);
    hottower.push_back(24382);

    hottower.push_back(5690);
    hottower.push_back(7062);
    hottower.push_back(7495);
    hottower.push_back(8916);
    hottower.push_back(9204);
    hottower.push_back(9205);
    hottower.push_back(9285);
    hottower.push_back(12714);

    hottower.push_back(6918);
    hottower.push_back(6919);
    hottower.push_back(8917);
    hottower.push_back(9356);
    hottower.push_back(10199);

    hottower.push_back(7135);
    hottower.push_back(7206);
    hottower.push_back(7279);
    hottower.push_back(7350);
    hottower.push_back(9767);
    hottower.push_back(10148);
    hottower.push_back(10149);
    hottower.push_back(10198);
    hottower.push_back(10270);
    hottower.push_back(10271);

    hottower.push_back(6991);
    hottower.push_back(7134);
    hottower.push_back(7494);
    hottower.push_back(9336);
    hottower.push_back(10016);
    hottower.push_back(10094);

    hottower.push_back(10021);
    hottower.push_back(10166);
    hottower.push_back(10234);
    hottower.push_back(10292);

    hottower.push_back(9932);
    hottower.push_back(10162);

    hottower.push_back(10018);
    hottower.push_back(10090);

    hottower.push_back(7567);

    hottower.push_back(1778);
    */
  }
  else if(species == 4)//run8 dAu ERT
  {
    //5 sigma
    hottower.push_back(2292);
    hottower.push_back(3737);
    hottower.push_back(4393);
    hottower.push_back(5690);
    hottower.push_back(7566);
    hottower.push_back(8916);
    hottower.push_back(9204);
    hottower.push_back(9205);
    hottower.push_back(9285);
    hottower.push_back(9336);
    hottower.push_back(9356);
    hottower.push_back(11388);
    hottower.push_back(13079);
    hottower.push_back(14314);
    hottower.push_back(15309);
    hottower.push_back(15381);
    hottower.push_back(19419);
    hottower.push_back(21155);
    hottower.push_back(21258);
    hottower.push_back(21768);
    hottower.push_back(22607);

    hottower.push_back(7423);
    hottower.push_back(8917);
    hottower.push_back(10017);
    hottower.push_back(10088);
    hottower.push_back(10089);
    hottower.push_back(10160);
    hottower.push_back(10232);
    hottower.push_back(10233);
    hottower.push_back(10270);
    hottower.push_back(10945);
    hottower.push_back(11377);
    hottower.push_back(11449);
    hottower.push_back(12302);
    hottower.push_back(12511);
    hottower.push_back(12564);
    hottower.push_back(13321);
    hottower.push_back(13667);
    hottower.push_back(16963);
    hottower.push_back(19489);
    hottower.push_back(19784);
    hottower.push_back(20698);
    hottower.push_back(21060);
    hottower.push_back(21156);
    hottower.push_back(21553);
    hottower.push_back(21665);
    hottower.push_back(21680);
    hottower.push_back(21730);
    hottower.push_back(21794);
    hottower.push_back(21860);
    hottower.push_back(22485);
    hottower.push_back(23460);
    hottower.push_back(24163);
    hottower.push_back(24286);

    hottower.push_back(10092);
    hottower.push_back(10093);
    hottower.push_back(10161);
    hottower.push_back(10164);
    hottower.push_back(10165);
    hottower.push_back(15912);
    hottower.push_back(16033);
    hottower.push_back(17285);
    hottower.push_back(17709);
    hottower.push_back(17799);
    hottower.push_back(19494);
    hottower.push_back(19688);
    hottower.push_back(19777);
    hottower.push_back(19873);
    hottower.push_back(19980);
    hottower.push_back(20775);
    hottower.push_back(24100);

    hottower.push_back(19070);
    hottower.push_back(20955);
    hottower.push_back(22714);

    hottower.push_back(1390);
    hottower.push_back(5465);
    hottower.push_back(5742);
    hottower.push_back(5743);
    hottower.push_back(7368);
    hottower.push_back(9767);
    hottower.push_back(10198);
    hottower.push_back(10585);
    hottower.push_back(13102);
    hottower.push_back(13177);
    hottower.push_back(21026);

    hottower.push_back(7278);
    hottower.push_back(10199);
    hottower.push_back(10271);
    hottower.push_back(14974);
    hottower.push_back(18095);
    hottower.push_back(21594);
    hottower.push_back(23412);
    hottower.push_back(24382);

    hottower.push_back(875);
    hottower.push_back(6775);
    hottower.push_back(7567);
    hottower.push_back(10016);
    hottower.push_back(10021);
    hottower.push_back(10148);

    /*
    //10sigma
    hottower.push_back(5690);
    hottower.push_back(7566);
    hottower.push_back(8916);
    hottower.push_back(9204);
    hottower.push_back(9205);
    hottower.push_back(9285);
    hottower.push_back(9356);
    hottower.push_back(11388);
    hottower.push_back(14314);
    hottower.push_back(19419);
    hottower.push_back(21155);
    hottower.push_back(21258);
    hottower.push_back(21768);
    hottower.push_back(22607);

    hottower.push_back(8917);
    hottower.push_back(9336);
    hottower.push_back(10088);
    hottower.push_back(10232);
    hottower.push_back(10270);
    hottower.push_back(12302);
    hottower.push_back(16963);
    hottower.push_back(19489);
    hottower.push_back(19784);
    hottower.push_back(21156);
    hottower.push_back(21730);
    hottower.push_back(21860);
    hottower.push_back(24286);

    hottower.push_back(10017);
    hottower.push_back(10089);
    hottower.push_back(10160);

    hottower.push_back(7423);
    hottower.push_back(21680);

    hottower.push_back(7278);
    hottower.push_back(10092);
    hottower.push_back(10164);
    hottower.push_back(10165);
    hottower.push_back(10198);
    hottower.push_back(10233);
    hottower.push_back(24163);

    hottower.push_back(9767);
    hottower.push_back(10093);
    hottower.push_back(10161);
    hottower.push_back(10271);

    hottower.push_back(10016);
    hottower.push_back(10199);
    hottower.push_back(24382);

    hottower.push_back(10021);
    hottower.push_back(23460);
    */
  }
  else
  {
    return false;
  }

  for(unsigned int i = 0; i < hottower.size(); i++)
  {
    //if it is hot tower
    if(towerid == hottower[i])
    {
      /*
      if(verbosity)
      {
	cout<<"tower id = "<<towerid<<" reject!!! ***************************"<<endl;
      }
      */
      //cout<<"tower id = "<<towerid<<" reject!!! ***************************"<<endl;

      return true;
      break;
    }
  }
  return false;
}

