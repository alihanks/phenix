/*
Construct a background to subtract from the output of AnaAwaySide
-man 4/22/05

code is adjusted to do in-run mixing. 

Runs on CNT/PWG but replicates cuts used in the foreground
Uses private version of preco/recal that is hacked to replicate the cuts
in the hard pDSTs.

*/



#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <PHGlobal.h>
#include <PHCentralTrack.h>
#include <emcClusterContainer.h>
#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <Fun4AllServer.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TFormula.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <combOnetrigbgSimple.h>
#include <PHObject.h>
#include <PHAngle.h>
#include <PHSnglCentralTrack.h>
#include <emcClusterContent.h>
#include <utiCentrality.h>
#include <ReactionPlaneObject.h>
#include <getClass.h>
#include <RunHeader.h>
#include <EventHeader.h>


typedef PHIODataNode <PHObject> PHObjectNode_t;
typedef PHIODataNode <PHGlobal> PHGlobalNode_t;
typedef PHIODataNode <PHCentralTrack> PHPNode_t;


using namespace std;
using namespace findNode;

const double DEG_PER_RAD = 180.0 / M_PI;
const double pi = M_PI;

combOnetrigbgSimple::combOnetrigbgSimple( int species, int flag, int mult, int minbiasert, float centlo, float centhi, int rxnmix) : CombinedSimple(NULL, centlo, centhi, flag, 0, species, 1)
{
  minbias_ert=minbiasert;
  tagflag=flag;

  ThisName = "combOnetrigbgSimple";
  multiplier=mult;  
  
  locent=centlo;
  hicent=centhi;
  mixrxn=rxnmix;

  for (int i = 0; i < 25; i++) {
    htrtest[i] = 0;    
    for (int j = 0; j < 25; j++) {
      h3pc[i][j] = 0;
      hpctest[i][j] = 0;
    }
  }

  //defaults 
  _cntnodename = new TString("PHCentralTrack");
  _emcnodename = new TString("emcClusterContainer");  

}



int combOnetrigbgSimple::Init(PHCompositeNode *topNode)
{
  //  Fun4AllServer * se = Fun4AllServer::instance();
  CombinedSimple::Init(topNode);  
  CombinedSimple::SetTagFlag(tagflag);

  //counters
  total_events=0;
  used_events_run=0; 
  used_events_trig=0; 

  if(mixrxn) cout<<PHWHERE<<": rxn mixing turned ON "<<endl;
  else cout<<PHWHERE<<": rxn mixing turned OFF "<<endl;

  isfinished=0;
  return 0;

}



int combOnetrigbgSimple::process_event(PHCompositeNode *topNode)
{

  RunHeader* d_runhdr = findNode::getClass<RunHeader>(topNode, "RunHeader");
  EventHeader* sync = findNode::getClass<EventHeader>(topNode, "EventHeader");
  int runno=d_runhdr->get_RunNumber();
  int evtno = (float)sync->get_EvtSequence();
  int seqno = evtno/100000;



  if(seqno_bypass){
    if(seqno<trig_photon_seqno) return 0;
  }
  
  seqno_bypass=1;

  //finish, if we've run over enough minbias events  -- set to 999 for train
  if(used_events_trig==multiplier){
    isfinished=1;
    return -2;  
  }
  total_events++;
  
  global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if (!global) cout <<PHWHERE<< ": no global!!" <<endl;

  particle = findNode::getClass<PHCentralTrack>(topNode, _cntnodename->Data());
  if (!particle) 
  {
    cout << PHWHERE << ": PHCentralTrack not in Node Tree" << endl;
    return -2;
  }

  if (species>0){
    PHTypedNodeIterator<ReactionPlaneObject> riter(topNode);
    reacplane = findNode::getClass<ReactionPlaneObject>(topNode, "ReactionPlaneObject");
    if( !reacplane ) {
      cout << PHWHERE << ": Could not find " << "ReactionPlaneObject" << " Node" << endl;
      total_events++;
      return 0;
    }
    
    emccluster = findNode::getClass<emcClusterContainer>(topNode,_emcnodename->Data() );
    
    if(!emccluster){
      cout << PHWHERE << ": emcClusterContainer not in Node Tree"<<endl;
      return -2;
    }
  }
  
  
  if(total_events%1000==0)cout<<PHWHERE<<": total events = "<<total_events<<" used_events_trig "<<used_events_trig<<" seqno "<<seqno<<endl; 
  
  float zvertex;

  if (species==2)
    {
      zvertex = global->getZVertex();
    }else{
      zvertex = global->getBbcZVertex();
    }


  if (fabs(zvertex)>30||zvertex==0) return 0;
  
  float bbcqs=0;
  float bbcqn=0;
  float bbcq=0;
  float zdcen=0; 
  float zdces=0;
  float thetaBBC =0;
  float thetaDEG = 0;
  int percent = 0;
  
  
  if(species>0)
    {

      bbcqs = global->getBbcChargeS();
      bbcqn = global->getBbcChargeN();
      bbcq = 0.5 * (bbcqs + bbcqn);
      zdcen = global->getZdcEnergyN();
      zdces = global->getZdcEnergyS();
      
      if (species == 1)
        percent = (int)global->getCentrality();
      else 
        percent = PhUtilities::getCentralityByClockRun4(bbcqn,bbcqs,zdcen,zdces,runno);
      
      if (percent<=0) percent = 0;
      if (percent>100) percent = -1;
      
      
      thetaBBC =reacplane->getBBCrp12(); 
      thetaDEG = thetaBBC * DEG_PER_RAD;
            
      if(thetaDEG<-90||thetaDEG>90) return 0;  
      
    }
  

  float ncent_fg=TMath::Floor((percent-1)/5);
  float ncent_bg=TMath::Floor((trig_photon_cent-1)/5);
  
  //find out whether this minbias file matches the trigger
  float nvert_fg=TMath::Floor(zvertex/5.0);
  float nvert_bg=TMath::Floor(trig_photon_vert/5.0);
  float nrp_fg=TMath::Floor(thetaDEG/30.0);
  float nrp_bg=TMath::Floor(trig_photon_rp/30.0);
  
  if (species>0) 
    {
      if(ncent_fg != ncent_bg) 
	{
	  cout<<PHWHERE<<": ncent_fg "<<nvert_fg<<" ncent_bg "<<nvert_fg<<endl;
	  return 0; 
	}
    }
  if (nvert_fg != nvert_bg)return 0;
  
  if(mixrxn)
    {
      if(nrp_fg != nrp_bg)return 0; 
    } 
  
  
  
  if(fabs(zvertex-trig_photon_vert)<0.001)
    {
      cout<<PHWHERE<<": same fg & bg "<<endl;
      return 0;
    }
  
  
  
  if(fabs(trig_photon_vert-zvertex)>5)
    cout<<PHWHERE<<": fg vert "<<trig_photon_vert<<" bg vert "<<zvertex<<endl;
  
    
  used_events_trig++;
  used_events_run++;
 
  RPLANE->Fill(thetaDEG);
  CENTRALITY->Fill(percent); 
 
  if(total_events%10000==0) cout<<PHWHERE<<": total  " <<total_events<<endl; 
  if(used_events_run%10000==0) cout<<PHWHERE<<": used  " <<used_events_run<<endl; 
  
  LoopThruTracks(trigger, particle);
  
  return 0;
}

int combOnetrigbgSimple::ResetEvent(PHCompositeNode *topNode)
{
  
  CombinedSimple::ResetEvent(topNode);

  return 0;
}

int combOnetrigbgSimple::End(PHCompositeNode *topNode){
  cout<<PHWHERE<<"total_counter "<<total_events<<endl;
  return 0;
}

void combOnetrigbgSimple::SetTrigProps(float E, float px, float py, float pz, float the, float eta, float phi, float x, float y, float z, float sector, float pc3dr, float emctof,  float vertex, float centrality, float reactionplane, float runno, float seqno){
 
  if(species==0)
    {
      centrality=0;
      reactionplane=0;
    }
  
  trig_photon_vert=vertex;
  trig_photon_cent=centrality;
  trig_photon_rp=reactionplane;
  
  trig_photon_runno=runno;
  trig_photon_seqno=seqno;

  trigger.E=E;
  trigger.px=px;
  trigger.py=py;
  trigger.pt=sqrt(px*px+py*py);
  trigger.pz=pz;
  trigger.the=the;
  trigger.eta=eta;
  trigger.phi=phi;
  trigger.x=x;
  trigger.y=y;
  trigger.z=z;
  trigger.sector=sector;
  trigger.pc3dr=pc3dr;
  trigger.emctof=emctof;

  cur_max_buf = 0;
  last_buf = -1;

   return;
}


int combOnetrigbgSimple::IsFinished(){
  return isfinished;

}


void combOnetrigbgSimple::ResetFinished(){
  used_events_trig=0;
  isfinished=0;
  return;
}

  
void combOnetrigbgSimple::SetupProximityMixing(){

  seqno_bypass=1;
  return;
}
