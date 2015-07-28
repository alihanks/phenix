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
#include <fstream>
#include <string>
#include <vector>
#include <PHGlobal.h>
#include <PHCentralTrack.h>
#include <PHCentralTrackv14.h>
#include <PHCompositeNode.h>
#include <PHIODataNode.h>
#include <Fun4AllServer.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <combOnetrigbgTrack.h>
#include <PHObject.h>
#include <PHAngle.h>
#include <PHSnglCentralTrackv14.h>
#include <emcClusterContainer.h>
#include <emcClusterContent.h>
#include <TLorentzVector.h>
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

typedef PHIODataNode <PHObject> PHObjectNode_t;
typedef PHIODataNode <PHGlobal> PHGlobalNode_t;
typedef PHIODataNode <PHCentralTrack> PHPNode_t;


using namespace std;
using namespace findNode;

const double DEG_PER_RAD = 180.0 / M_PI;
const double pi = M_PI;

combOnetrigbgTrack::combOnetrigbgTrack( int flag, int mult, int minbiasert, int set_compressed, int set_species, int centlo, int centhi) : CombinedSimple(NULL, centlo, centhi, flag, set_species, 1,0)
{
  minbias_ert=minbiasert;
  tagflag=flag;
  cout << "note: iscompressed == 0 mode is no longer supported in combOnetrigTrack, setting to 1" << endl; 
  iscompressed = 1;

  species = set_species;
  
  ThisName = "combOnetrigbgTrack";
  // mixFilepath defaults to run7
  _mixFilepath = "XXXXXXXXXXXXXXXX";
  multiplier=mult;  
  _segFile = -1;
  _mixTrig = 0;
  _mixSing = 0;
  
}




int combOnetrigbgTrack::Init(PHCompositeNode *topNode)
{
  cout<<" Initializing combOnetrigbgTrack "<< endl;

  CombinedSimple::Init(topNode);
  //delete topNode;

  CombinedSimple::SetTagFlag(tagflag);
  CombinedSimple::SetSpeciesFlag(species);

  total_events=0;
  used_events_run=0; 
  used_events_trig=0; 
  isfinished=0;

  trig_photon_vert = 0;
  trig_photon_runno = 0;
  trig_photon_cent = 0;
  return 0;
}


int combOnetrigbgTrack::process_event(int trackindex)
{
 
  //since this is a subsysreco but I'm not actually reading a dst
  //I overloaded the return code to the number of tracks
  //  -2 still means abort run
  //Important -- return codes don't really do anything right now except the -2


  //finish, if we've run over enough minbias events  -- set to 999 for train
  if(used_events_trig==multiplier){
    isfinished=1;
    return -2;  
  }
  total_events++;


  //global properties
  partners->SetBranchAddress("zvertex",  &as_zvertex);
  partners->SetBranchAddress("N",        &as_N);
  if(species>0) partners->SetBranchAddress("centrality",   &as_centrality);
  if(species>0) partners->SetBranchAddress("n0", &as_n0);

  //  track branches...arrays of length N
  partners->SetBranchAddress("px",       &as_px);
  partners->SetBranchAddress("py",       &as_py);
  partners->SetBranchAddress("pz",       &as_pz);
  partners->SetBranchAddress("charge",   &as_charge);
  partners->SetBranchAddress("m2emc",    &as_m2emc);
  partners->SetBranchAddress("m2tof",    &as_m2tof);
  partners->SetBranchAddress("ppc3x",    &as_ppc3x);
  partners->SetBranchAddress("ppc3y",    &as_ppc3y);
  partners->SetBranchAddress("ppc3z",    &as_ppc3z);

  partners->GetEntry(trackindex);

  if (species>0 && as_centrality < 0) return 0;
  
  //find out whether this minbias event matches the trigger
  float nvert_fg=TMath::Floor(as_zvertex/5.0);
  float nvert_bg=TMath::Floor(trig_photon_vert/5.0);


  if( nvert_fg != nvert_bg){
    return 0;
  }

  if(species>0){
    float ncent_fg=TMath::Floor((as_centrality-1)/5);
    float ncent_bg=TMath::Floor((trig_photon_cent-1)/5);
    
    if(ncent_fg != ncent_bg) 
    {
      return 0; 
    }
  }

  if(species==0){
    if(fabs(as_zvertex-trig_photon_vert)<0.001)
    {
      if(verbosity) cout<<PHWHERE<<": same fg & bg "<<endl;
      return 0;
    }
  }
  if(species>0){
    float cent_compare=(float)(as_centrality-trig_photon_cent);
    if(fabs(as_zvertex-trig_photon_vert)<0.001&&fabs(cent_compare)<1.0)
    {
      if(verbosity) cout<<PHWHERE<<": same fg & bg "<<endl;
      return 0;
    }
  }

  TRIGPT->Fill(trigger.pt);

  float mwweight[5]={0};
  float awweight[5]={0};
  
  float mwweightfine[7]={0};
  
  float pi0trigeff=0.;
  if(tagflag>0)pi0trigeff=grpi0eff->Eval(trigger.pt);


  int trigptbin = hshark_large[0][0]->FindBin(trigger.pt);
  if(trigptbin>400) trigptbin=400;

  float pi0zemc  = 510.0*trigger.pz/trigger.pt+trigger.zvertex;
  int ipi0zemc = (int)TMath::Floor((pi0zemc+165.0)/10.);
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
      }
    }
    for(int idecs=0;idecs<7;idecs++){
      float mattshark=hshark_small[idecs][ipi0zemc]->GetBinContent(trigptbin);
      mwweightfine[idecs]=mattshark*pi0trigeff;
      if(mwweightfine[idecs]>0)MWFINETRIGPT->Fill(idecs,mwweightfine[idecs]);
    }
  }


  PartnerHolder partnertrack;

  int ngoodtracks=0;

  for(int i=0;i<as_N;i++){
    
    float mom = sqrt(as_px[i]*as_px[i]+as_py[i]*as_py[i]+as_pz[i]*as_pz[i]);
    
    float the0 = acos(as_pz[i]/mom);
    float phi0 = atan2(as_py[i],as_px[i]);
    int dcarm = 0;
    if(phi0<1.5) dcarm=1;

    partnertrack.vm2emc=as_m2emc[i];
    partnertrack.vm2tof=as_m2tof[i];
    partnertrack.vcharge=as_charge[i]/abs(as_charge[i]);
    partnertrack.vppc3x=as_ppc3x[i];
    partnertrack.vppc3y=as_ppc3y[i];
    partnertrack.vppc3z=as_ppc3z[i];
    partnertrack.vthe0=the0;
    partnertrack.vmom=mom;
    partnertrack.vphi0=phi0;
    partnertrack.vdcarm=dcarm;
    partnertrack.vn0=(int)as_n0[i];
    partnertrack.vquality=abs(as_charge[i]);

    //set to good track qualities
    partnertrack.vpc3dphi=0.;
    partnertrack.vpc3dz=0.;
    partnertrack.vemcdphi=0.;
    partnertrack.vemcdz=0.;
    partnertrack.vpc3sdz=0.;
    partnertrack.vpc3sdphi=0.;
    partnertrack.vemcsdz=0.;
    partnertrack.vemcsdphi=0.;
    
    int trackused=LoopThruTrack(trigger, partnertrack,0);

    ngoodtracks+=trackused;
  }
  
  if(total_events%1000==0)cout<<PHWHERE<<": total events = "<<total_events<<" used_events_trig "<<used_events_trig<<endl;
  
  //just a dummy check
  if(fabs(trig_photon_vert-as_zvertex)>5) 
    if(verbosity) 
      cout<<PHWHERE<<": fg vert "<<trig_photon_vert<<" bg vert "<<as_zvertex<<endl;
  
  used_events_trig++;
  used_events_run++;
  
  if(total_events%10000==0) cout<<PHWHERE<<": total  " <<total_events<<endl; 
  if(used_events_run%10000==0) cout<<PHWHERE<<": used  " <<used_events_run<<endl; 
  
  return 1;

}

int combOnetrigbgTrack::End(char *outfile){

  cout<<"combOnetrigbgTrack::End(): total_counter "<<total_events<<endl;
  return 0;
}

void combOnetrigbgTrack::SetTrigProps(float E, float px, float py, float pz, float the, float eta, float phi, float x, float y, float z, float sector, float pc3dr, float emctof,  float vertex, float centrality){
 
  trig_photon_vert=vertex;
  trig_photon_cent=centrality;
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

   return;
}


int combOnetrigbgTrack::IsFinished(){
  return isfinished;

}


void combOnetrigbgTrack::ResetFinished(){
  used_events_trig=0;
  isfinished=0;
  return;
}


void combOnetrigbgTrack::SetMixingChainFile(char *filelistname){
  
  ifstream filelist(filelistname);
  partners=new TChain("singles");
  
  TChain * chpartners = (TChain *) partners;
  
  while(!filelist.eof()){
    char filename[1000];
    filelist>>filename;
    chpartners->Add(filename);
  }
  
  filelist.close();
  
}

void combOnetrigbgTrack::DeleteMixingChainFile(){
  if(partners) delete partners;

}


void combOnetrigbgTrack::Mix(TNtuple * trig, TTree * inpartners)
{
  
  int nentries = trig->GetEntries();
  cout << "combOnetrigbgTrack::Mix(): here we are" << endl;
  float vertex=0;
  float runnumber=0;
  float E=0;
  float px=0;
  float py=0;
  float pz=0;
  float x=0;
  float y=0;
  float z=0;
  float eta=0;
  float theta=0;
  float phi=0;
  //   float id=0;
  float sector=0;
  float pc3dr=0;
  float emctof=0;
  float seqno_fg=0;
  float percent=0;
  
  trig->SetBranchAddress("runno",&runnumber);
  trig->SetBranchAddress("seqno",&seqno_fg);
  trig->SetBranchAddress("E",&E);
  trig->SetBranchAddress("px",&px);
  trig->SetBranchAddress("py",&py);
  trig->SetBranchAddress("pz",&pz);
  trig->SetBranchAddress("x",&x);
  trig->SetBranchAddress("y",&y);
  trig->SetBranchAddress("z",&z);
  trig->SetBranchAddress("eta",&eta);
  trig->SetBranchAddress("the",&theta);

  trig->SetBranchAddress("phi",&phi);
  trig->SetBranchAddress("vertex",&vertex);
  trig->SetBranchAddress("sector",&sector);
  trig->SetBranchAddress("pc3dr",&pc3dr);
  trig->SetBranchAddress("emctof",&emctof);
  if(species>0)trig->SetBranchAddress("percent",&percent);

  int trigs_forthis_rungroup=0;
  
  int nentlo= 0;
  int nenthi= nentries;
  
  if( verbosity ) cout <<PHWHERE<<": nentries "<<nentries<<endl;
  if( verbosity ) cout<<PHWHERE<<":  nentlo "<<nentlo<<" nenthi "<<nenthi<<endl;
  
  int lastrunno=0;
  
  for(int i=nentlo;i<nenthi;i++){
    
    trig->GetEntry(i);
    
    if( verbosity ) cout<<PHWHERE<<": trig # "<<i<<" runnumber "<<runnumber<<" nentlo "<<nentlo<<" nenthi "<<nenthi<<endl;
    
    SetTrigProps(E,px,py,pz,theta,eta,phi,x,y,z,sector,pc3dr,emctof,vertex,percent);
    
    if (inpartners) //then use inpartners
    {
      partners = inpartners;
    }
    else //fill partners based on runnumber and species
    {
      int irunnumber=(int)runnumber;
      
      if(irunnumber!=lastrunno){
        
        if(lastrunno>0){ 
          delete partners;
          DeleteMixingChainFile();
        }
        
        lastrunno=irunnumber;
        
        char filelistname[1000];

        if (_mixFilepath.Contains("XXXX"))
        {
          cout << "FATAL ERROR: Mixing path not set/mixing file not found: must be set in the macro"
          << endl;
          exit(0);
        }

        if (_segFile >= 0)
          sprintf(filelistname, "%s%d_%d.list",  _mixFilepath.Data(),irunnumber,_segFile);
        else
          sprintf(filelistname, "%s%d.list",  _mixFilepath.Data(),irunnumber);
        
        //cout<<" grabbed file for run "<<irunnumber<<endl;
        
        trigs_forthis_rungroup++;
        
        ifstream filelist(filelistname);
        if (filelist.fail()) continue;
        
        vector<string>filenames;
        
        partners= new TChain("singles");
        TChain * chpartners = (TChain *) partners;
        
        while(!filelist.eof()){
          char filename[1000];
          filelist>>filename;
          filenames.push_back(filename);
          chpartners->Add(filename);
        }
        
        filelist.close();
        SetMixingChainFile(filelistname);
      }
      
    }
    
    int istrigfinished=0;
    ResetMultCounter(); 
    ResetFinished();
    
    //find the triggered event/ segment
    int Tnent = partners->GetEntries();
    
    int seqno_bg, evtno_bg;
    int mixby_seqno_evtno=0;
    
    if(mixby_seqno_evtno==0)partners->SetBranchAddress("segno",&seqno_bg);
    if(mixby_seqno_evtno==1)partners->SetBranchAddress("evtno",&evtno_bg);    
    
    if(partners){
      
      int num2skip=0;
      
      partners->GetEntry(1);
      //Use in-segment mixing:
      if(mixby_seqno_evtno==0){
        int trigfound=0;
        while(num2skip<Tnent){
          int nument=++num2skip;
          partners->GetEntry(nument);
          if(seqno_bg>=seqno_fg){ 
            trigfound=1;
            break;
          }
        }
        if(trigfound==0){
          cout<<" segment not found !!! "<<endl;
          num2skip=0;
        }
      }

      int ic=0;
      while(ic<10&&istrigfinished==0){
        
        for(int ip=(int)num2skip;ip<Tnent;ip++){
          if(ip>=Tnent-1){
            ip=0;
            ic++;
          }
          
          partners->GetEntry(ip);
          int returncode=0; 
          returncode = process_event(ip);
          
          // this means we have enough pairs
          if(returncode==-2){
            istrigfinished=1;
            break;
          }
        }
      }
    }
    else{
      cout<<PHWHERE<<": couldn't find mixing file "<<endl;
    }
    
    if(istrigfinished==0) cout<<PHWHERE<<": didn't finish trigger number "<<i<<" from runnumber "<<runnumber<<endl;
  }
  
  if( verbosity ) cout<<PHWHERE<<": ran over "<<trigs_forthis_rungroup<<" trigs "<<endl;
  
}




int combOnetrigbgTrack::process_event(PHCompositeNode * topNode)
{
  //do nothing:  this node should only function in the Init and End
  // when registered with the fun4all mgr
  return 0;
}

int combOnetrigbgTrack::End(PHCompositeNode * topNode)
{
  if (_mixTrig && _mixSing)
    Mix(_mixTrig->get_triggaz(), _mixSing->get_assocNtuple()); 

  return 0;
}

