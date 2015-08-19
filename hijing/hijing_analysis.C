#include <hijing_analysis.h>

#include <THmulf.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <iostream>
#include <Fun4AllHistoManager.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <PHAngle.h>
#include <Fun4AllReturnCodes.h>
#include <PHHepMCGenEvent.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <AParticle.h>
#include <ACluster.h>
#include <ATrack.h>
#include <APiZero.h>
#include <AMixingTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom.h>

using namespace std;
using namespace findNode;

hijing_analysis::hijing_analysis(const char* output, const char* name) 
{
  ThisName  = name;
  verbosity = 0;
  outfile = output;
  nevents = 0;
  Rcut = 0.3;
  Centrality = 0;
  _MaxEta = 1.0;
  _MinAssocPt = 0.5;
  _MinTrigPt = 1.0;
  NMIX = 500;
}

hijing_analysis::~hijing_analysis(){
   delete manager;
}

int hijing_analysis::Init(PHCompositeNode *topNode){

  cout << "Initializing histograms" << endl;
  nevents = 0;
  manager = new Fun4AllHistoManager("hijing_analysis");

  string name;
  ostringstream bin;
  
  bin.str("");
  bin << "_c" << 0;
  name = "hmult" + bin.str();
  Init1DHisto(hmult, name.c_str(),"multiplicity",100,0,100.0);

  name = "h1_mass" + bin.str();
  Init1DHisto(h1_mass, name.c_str(),"mass [GeV/c^2]",500,0,1.0);
  
  name = "h1_trig_pt_all" + bin.str();
  Init1DHisto(h1_trigger_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trig_pt_pi0" + bin.str();
  Init1DHisto(h1_trigger_pi0_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trig_pt_inc" + bin.str();
  Init1DHisto(h1_trigger_dir_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trig_pt_inc_iso" + bin.str();
  Init1DHisto(h1_trigger_iso_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trig_pt_pi0_iso" + bin.str();
  Init1DHisto(h1_trigger_iso_pi0_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trig_pt_iso_dir" + bin.str();
  Init1DHisto(h1_trigger_iso_dir_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trig_pt_dec" + bin.str();
  Init1DHisto(h1_trigger_dec_pt, name.c_str(),"p_{T} bin",5,-0.5,4.5);
  
  name = "h1_trig_pt_dec_iso" + bin.str();
  Init1DHisto(h1_trigger_dec_iso_pt, name.c_str(),"p_{T} bin",5,-0.5,4.5);
  
  name = "h3_dphi_fold" + bin.str();
  Init3DHisto(h3_dphi, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0, PI);
  
  name = "h3_dphi_iso_fold" + bin.str();
  Init3DHisto(h3_dphi_iso, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0, PI);
  
  name = "h3_dphi_dir_fold" + bin.str();
  Init3DHisto(h3_dphi_dir, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0, PI);
  
  name = "h3_dphi_pi0_fold" + bin.str();
  Init3DHisto(h3_dphi_pi0, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0, PI);
  
  name = "h3_dphi_dir_iso_fold" + bin.str();
  Init3DHisto(h3_dphi_dir_iso, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0, PI);
  
  name = "h3_dphi_pi0_iso_fold" + bin.str();
  Init3DHisto(h3_dphi_pi0_iso, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0, PI);
  
  name = "h3_dphi_mix_fold" + bin.str();
  Init3DHisto(h3_dphi_mix, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);

  name = "h3_dphi_mix_iso_fold" + bin.str();
  Init3DHisto(h3_dphi_mix_iso, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);

  name = "h3_dphi_pi0_mix_fold" + bin.str();
  Init3DHisto(h3_dphi_pi0_mix, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);

  name = "h3_dphi_pi0_mix_iso_fold" + bin.str();
  Init3DHisto(h3_dphi_pi0_mix_iso, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);

  TH2F* temp2;
  for(int ipt=0; ipt<5; ipt++){
    bin.str("");
    bin << "_p" << ipt <<"_c"<<0;
    name = "h2_dphi_dec_fold" + bin.str();
    Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
    h2_dphi_dec.push_back(temp2);
    name = "h2_dphi_dec_iso_fold" + bin.str();
    Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
    h2_dphi_dec_iso.push_back(temp2);
    name = "h2_dphi_dec_mix_fold" + bin.str();
    Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
    h2_dphi_dec_mix.push_back(temp2);
    name = "h2_dphi_dec_mix_iso_fold" + bin.str();
    Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
    h2_dphi_dec_mix_iso.push_back(temp2);
  }

  name = "h2_cluster_dphi_dR";
  Init2DHisto(h2_cluster_wdR, name.c_str(), "#Delta R ", 200, 0.0, 2.0, "#Delta #phi [rad]", 60, 0.0, PI);
  
  name = "h2_cluster_pi0_dphi_dR";
  Init2DHisto(h2_cluster_pi0_wdR, name.c_str(), "#Delta R ", 200, 0.0, 2.0, "#Delta #phi [rad]", 60, 0.0, PI);
  
  name = "h2_cluster_dphi_dR";
  Init2DHisto(h2_cluster_dir_wdR, name.c_str(), "#Delta R ", 200, 0.0, 2.0, "#Delta #phi [rad]", 60, 0.0, PI);
  
  name = "h3_cluster_etot_dR";
  Init3DHisto(h3_cluster_etot, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{tot} [GeV]", 40, 0.0, 2.0,"#Delta R",200,0.0,2.0);
  
  name = "h3_cluster_etot_pi0_dR";
  Init3DHisto(h3_cluster_pi0_etot, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{tot} [GeV]", 40, 0.0, 2.0,"#Delta R",200,0.0,2.0);
  
  name = "h3_cluster_etot_dir_dR";
  Init3DHisto(h3_cluster_dir_etot, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{tot} [GeV]", 40, 0.0, 2.0,"#Delta R",200,0.0,2.0);
  
  name = "h3_cluster_dR";
  Init3DHisto(h3_cluster_dR, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{assoc} [GeV]", 40, 0.0, 10.0,"#Delta R",200,0.0,2.0);
  
  name = "h3_cluster_pi0_dR";
  Init3DHisto(h3_cluster_pi0_dR, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{assoc} [GeV]", 40, 0.0, 10.0,"#Delta R",200,0.0,2.0);
  
  name = "h3_cluster_dir_dR";
  Init3DHisto(h3_cluster_dir_dR, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{assoc} [GeV]", 40, 0.0, 10.0,"#Delta R",200,0.0,2.0);
  
  name = "h2_cluster_etot_dR";
  Init2DHisto(h2_cluster_etot, name.c_str(), "E_{tot}/p_{T}", 60, 0.0, 5.0, "p_{T} [GeV/c]", 60, 0.0, 15.0);
  
  name = "h2_cluster_etot_pi0_dR";
  Init2DHisto(h2_cluster_pi0_etot, name.c_str(), "E_{tot}/p_{T}", 60, 0.0, 5.0, "p_{T} [GeV/c]", 60, 0.0, 15.0);

  name = "h2_cluster_etot_dir_dR";
  Init2DHisto(h2_cluster_dir_etot, name.c_str(), "E_{tot}/p_{T}", 60, 0.0, 5.0, "p_{T} [GeV/c]", 60, 0.0, 15.0);
  
  atree = new AMixingTree();
  atree->SetTriggerBranches();
  atree->SetPartnerBranches();

  return 0;
}

int hijing_analysis::process_event(PHCompositeNode* topNode)
{
  if( verbosity )
    cout << "Event: " << nevents << endl;
  if( nevents%500==0 ) cout << "Event: " << nevents << endl;
  HepMC::GenEvent* evt = NULL;
  PHNodeIterator iter(topNode);
  PHHepMCGenEvent *genevent = findNode::getClass<PHHepMCGenEvent>(topNode,"PHHepMCGenEvent");
  evt = genevent->getEvent();
  if( !evt ) {
    cout << "Failed to find HepMC Node!" << endl;
    return -1;
  }
  
  nevents++;
  vector<ACluster*> clusters;
  vector<ATrack*> tracks;
  vector<APiZero*> pizeros;

  int mult = 0;
  for(HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p)
  {
    int id = (*p)->pdg_id();
    if( verbosity > 1 ) cout << "Checking particle with id = " << id << endl;
    double eta = (*p)->momentum().eta();
    if( eta > -4 && eta < -3 ) mult++;

    if( id==111 )
    {
      APiZero piz;
      if( MakePiZero((*p),&piz) ) {
        //dec weighting
        vector<float> mwweight;
        for(int i=0; i<5; i++) mwweight.push_back(0.0);
        EvalDecWeights(&piz,mwweight);
        piz.SetDecayWeights(mwweight);

        pizeros.push_back(piz.clone());
      }
    }
    if( (*p)->status()!=1 )continue;
    if( id==22 )
    {
      ACluster clus;
      if( MakeCluster((*p),&clus) ) {
        /*
        for( HepMC::GenVertex::particles_in_const_iterator ip = (*v)->particles_in_const_begin(); ip != (*v)->particles_in_const_end(); ++ip )
        {
          int pid = (*ip)->pdg_id();
          if( pid==111 || pid==221 || pid==223 ) clus.SetTag(true);
          if( verbosity ) cout << "Photon with E = " << clus.E() << " and parent id = " << pid << endl;
        }
        */
        clusters.push_back(clus.clone());
      }
    }
    if( (fabs(id)==211 || fabs(id)==321 || fabs(id)==2212) )
    {
      ATrack track;
      if( MakeTrack((*p),&track) ) {
        tracks.push_back(track.clone());
        atree->SetPartnerData(track.Pt(),track.Phi(),track.Eta(),track.E(),track.GetPemcx(),track.GetPemcy(),track.GetPemcz(),tracks.size()-1);
      }
    }
  }
  if( verbosity ) cout << "Found " << tracks.size() << " associated tracks" << endl;

  if( clusters.size() == 0 && pizeros.size() == 0 ) {
    cout << "Found no triggers: cleaning up." << endl;
    ClearVector(clusters);
    ClearVector(tracks);
    ClearVector(pizeros);
    return 0; 
  }

  for( unsigned int i = 0; i < clusters.size()-1; i++)
  {
    for( unsigned int j = i+1; j < clusters.size(); j++)
    {
      TLorentzVector pair = TLorentzVector(*clusters[i] + *clusters[j]);
      if( verbosity ) cout << "checking cluster pair with mass = " << pair.M() << endl;
      if( pair.M() > 0.120 && pair.M() < 0.160 )
      {
        cout << "Found pi0-pair with mass = " << pair.M() << endl;
        clusters[i]->SetTag(true);
        clusters[j]->SetTag(true);
      }
    }
  }
  
  if( verbosity ) cout << "Getting event multiplicity" << endl;
  int cent_bin = GetCentrality(mult);
  if( cent_bin != Centrality && Centrality>=0 ) {
    ClearVector(clusters);
    ClearVector(tracks);
    ClearVector(pizeros);
    return 0;
  }

  if( verbosity )
    cout << "Event multiplicity = " << mult << endl;
  hmult->Fill((float)mult);
  if( verbosity ) cout << "Looping over " << pizeros.size() << " pizeros" << endl;
  for( unsigned int i = 0; i < pizeros.size(); i++ )
  {
    if(pizeros[i]->Pt() < 4.0 ) continue;
    h1_trigger_pi0_pt->Fill(pizeros[i]->Pt());
    SetIso(pizeros[i],tracks,clusters,Rcut,h3_cluster_pi0_dR,h3_cluster_pi0_etot,h2_cluster_pi0_wdR,h2_cluster_pi0_etot);
    if( pizeros[i]->IsIso() ) h1_trigger_iso_pi0_pt->Fill(pizeros[i]->Pt());
    //cout << "Set iso<" << Rcut << " cut for piz to " << pizeros[i]->IsIso() << endl;
    //dec trigger counting
    vector<float> mwweight = pizeros[i]->GetDecayWeights();
    for(int ipw=0; ipw<5; ipw++){
      h1_trigger_dec_pt->Fill(ipw,mwweight[ipw]);
      if( pizeros[i]->IsIso() )
        h1_trigger_dec_iso_pt->Fill(ipw,mwweight[ipw]);
    }
    atree->SetTriggerData(pizeros[i]->Pt(),pizeros[i]->Phi(),pizeros[i]->Eta(),pizeros[i]->E(),
      ((ACluster*)pizeros[i]->Daughter1())->GetX(),((ACluster*)pizeros[i]->Daughter1())->GetY(),((ACluster*)pizeros[i]->Daughter1())->GetZ(),
      pizeros[i]->IsIso(),pizeros.size()-1);
    float ph_phi = PHAngle(pizeros[i]->Phi());
    vector<float> weight = pizeros[i]->GetDecayWeights();
    for( unsigned int j = 0; j < tracks.size(); j++ )
    {
      float trk_phi = PHAngle(tracks[j]->Phi());
      float dphifold = CalculateFoldedDphi(trk_phi,ph_phi);

      h3_dphi_pi0->Fill(pizeros[i]->Pt(), tracks[j]->Pt(), dphifold);
      if( pizeros[i]->IsIso() ) h3_dphi_pi0_iso->Fill(pizeros[i]->Pt(), tracks[j]->Pt(), dphifold);
      for(unsigned int ipw=0;ipw<weight.size();ipw++){
        if(weight[ipw]>0) {
          h2_dphi_dec[ipw]->Fill(dphifold,tracks[j]->Pt(),weight[ipw]);
          if( pizeros[i]->IsIso() ) h2_dphi_dec_iso[ipw]->Fill(dphifold,tracks[j]->Pt(),weight[ipw]);
        }
      }
    }
  }
  
  if( verbosity ) cout << "Looping over " << clusters.size() << " clusters" << endl;
  for( unsigned int i = 0; i < clusters.size(); i++ )
  {    
    if( clusters[i]->Pt() < 5.0 ) continue;
    h1_trigger_pt->Fill(clusters[i]->Pt());
    if( clusters[i]->IsTagged() ) continue;

    h1_trigger_dir_pt->Fill(clusters[i]->Pt());
    SetIso(clusters[i],tracks,clusters,Rcut,h3_cluster_dir_dR,h3_cluster_dir_etot,h2_cluster_dir_wdR,h2_cluster_dir_etot);
    if( clusters[i]->IsIso() ) h1_trigger_iso_dir_pt->Fill(clusters[i]->Pt());

    SetIso(clusters[i],tracks,clusters,Rcut,h3_cluster_dR,h3_cluster_etot,h2_cluster_wdR,h2_cluster_etot);
    if( clusters[i]->IsIso() ) h1_trigger_iso_pt->Fill(clusters[i]->Pt());

    atree->SetTriggerData(clusters[i]->Pt(),clusters[i]->Phi(),clusters[i]->Eta(),clusters[i]->E(),
      clusters[i]->GetX(),clusters[i]->GetY(),clusters[i]->GetZ(),clusters[i]->IsIso(),pizeros.size()+clusters.size()-1);
    float ph_phi = PHAngle(clusters[i]->Phi());
    for( unsigned int j = 0; j < tracks.size(); j++ )
    {
      float trk_phi = PHAngle(tracks[j]->Phi());
      float dphifold = CalculateFoldedDphi(trk_phi,ph_phi);

      h3_dphi->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
      if( clusters[i]->IsIso() ) h3_dphi_iso->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
      h3_dphi_dir->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
      if( clusters[i]->IsIso() ) h3_dphi_dir_iso->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
    }
  }

  ClearVector(clusters);
  ClearVector(tracks);
  ClearVector(pizeros);

  atree->SetEventData(nevents,0,mult,(int)clusters.size(),(int)pizeros.size(),(int)tracks.size());
  if(clusters.size() > 0 || pizeros.size() > 0) atree->_ttrig->Fill();
  atree->_tpart->Fill();

  return 0;
}

void hijing_analysis::EvalDecWeights(APiZero* pi0trigger, vector<float>& mwweight)
{
  float pi0trigpt = pi0trigger->Pt();
  
  int trigptbin = hshark[0]->FindBin(pi0trigpt);
  if(trigptbin>400) trigptbin=400;
  
  for(int idecl=0;idecl<5;idecl++){
    float mattshark=hshark[idecl]->GetBinContent(trigptbin);
    if(mattshark>0) mwweight[idecl]=mattshark;
    if( verbosity>1 ) cout << "setting weight = " << mattshark << " for pi0 with pt = " << pi0trigpt << endl;
  }
}

void hijing_analysis::SetSharkFin(const char* filename)
{
  TFile *fshark_exodus=NULL;
  fshark_exodus = new TFile(filename);
  
  cout<<PHWHERE<<" loading sharkfins "<< fshark_exodus->GetName() <<endl;
  for(int idecl=0;idecl<5;idecl++){
    char sharkname[100];
    sprintf(sharkname,"hshark_large_sum_%d",idecl);
    hshark[idecl]=(TH1D*)fshark_exodus->Get(sharkname);
  }  
}

void hijing_analysis::ApplyEnergyResolution(AParticle* mom4, int pbsc_pbgl)
{
  double pt, E, eta, phi;
  double sigma=0;

  E  = mom4->E();
  eta = mom4->Eta();
  phi   = mom4->Phi();

  if(pbsc_pbgl==0)  
    sigma = E*sqrt((0.06)*(0.06)+(0.081/sqrt(E))*(0.081/sqrt(E))); //more realistic
  if(pbsc_pbgl==1)
    sigma = E*sqrt((0.05)*(0.05)+(0.06/sqrt(E))*(0.06/sqrt(E))); 
  
  // cout<<" pbsc_pbgl "<<pbsc_pbgl<<" sigma "<<sigma<<endl;
  E  = gRandom->Gaus(E,sigma);
  pt = E/TMath::CosH(eta);

  //mom4->SetE(E);
  mom4->SetPtEtaPhiE(pt,eta,phi,E);

}

void hijing_analysis::DoMixing(TTree* trig, TTree* assoc, int size)
{
  int pooldepth = 0;
  int nloop = 0;

  int evt_trig;
  float cent_trig;
  int ntrig_photons;
  int ntrig_pi0s;
  int ntrig = 500;
  float trigpt[ntrig];
  float trigphi[ntrig];
  float trigeta[ntrig];
  float trige[ntrig];
  int iso[ntrig];

  int evt_part;
  float cent_part;
  int npart = 500;
  int nphotons;
  float partpt[npart];
  float partphi[npart];
  float parteta[npart];
  float parte[npart];

  trig->SetBranchAddress("evt", &evt_trig);
  trig->SetBranchAddress("cent", &cent_trig);
  trig->SetBranchAddress("nphotons", &ntrig_photons);
  trig->SetBranchAddress("npi0s", &ntrig_pi0s);
  trig->SetBranchAddress("ntrig", &ntrig);
  trig->SetBranchAddress("pt", trigpt);
  trig->SetBranchAddress("phi", trigphi);
  trig->SetBranchAddress("eta", trigeta);
  trig->SetBranchAddress("e",trige);
  trig->SetBranchAddress("iso",iso);

  assoc->SetBranchAddress("evt", &evt_part);
  assoc->SetBranchAddress("cent", &cent_part);
  assoc->SetBranchAddress("ntracks", &npart);
  assoc->SetBranchAddress("nphotons", &nphotons);
  assoc->SetBranchAddress("pt", partpt);
  assoc->SetBranchAddress("phi", partphi);
  assoc->SetBranchAddress("eta", parteta);
  assoc->SetBranchAddress("e", parte);

  int nentrig = trig->GetEntries();
  int nenpart = assoc->GetEntries();
  //cout<<"Mixing: nentrig = " <<nentrig<<"; nenpart = "<<nenpart<<endl;
 
  
  for(int i=0; i<nentrig; i++){
    //cout<<"trigger count = "<<i<<endl;
    pooldepth = 0;
    nloop = 0;
    trig->GetEntry(i);

    vector<ACluster*> clusters;
    vector<APiZero*> pizeros;

    for(int itrig=0; itrig<ntrig_pi0s; itrig++){
      APiZero pi0;
      MakePi0Object(&pi0, trigpt[itrig], trigphi[itrig], trigeta[itrig], trige[itrig]);
      pi0.SetIso(iso[itrig]);
      //dec weighting
      vector<float> mwweight;
      for(int i=0; i<5; i++) mwweight.push_back(0.0);
      EvalDecWeights(&pi0,mwweight);
      pi0.SetDecayWeights(mwweight);

      pizeros.push_back(pi0.clone());
    }

    for(int itrig=ntrig_pi0s; itrig<ntrig; itrig++){
      ACluster pho;
      MakeClusterObject(&pho, trigpt[itrig], trigphi[itrig], trigeta[itrig], trige[itrig]);
      pho.SetIso(iso[itrig]);
      clusters.push_back(pho.clone());
    }

    int ncent_fg = 0.;
    ncent_fg = GetCentrality((int)cent_trig);
    
    for(int j=0; j<nenpart; j++){
      assoc->GetEntry(j);

      //check if trigger and assoc belong to the same event
      if(evt_trig == evt_part){
        j = CheckPool(nenpart,j,pooldepth,size,nloop);
        //cout<<"same event. nloop = "<<nloop<<endl;
        if(nloop > NMIX) break;//to prevent infinite loop in case there is no interesting hadrons available 
        continue;
      }
      
      //check if trigger and partner belongs to the same event class
      int ncent_bg = GetCentrality((int)cent_part);
      if(ncent_fg != ncent_bg) {
        j = CheckPool(nenpart,j,pooldepth,size,nloop);
        //cout<<"not same centrality bin. nloop = "<<nloop<<endl;
        if(nloop > NMIX) break;
        continue;
      }

      vector<ATrack*> tracks;
      for(int ipart=0; ipart<npart; ipart++){
        ATrack trk;
        MakeTrackObject(&trk, partpt[ipart], partphi[ipart], parteta[ipart], parte[ipart]);
        tracks.push_back(trk.clone());
      }

      pooldepth++;
      if( verbosity>3 ) cout<<"pooldepth = "<<pooldepth<<endl;

      for( unsigned int i = 0; i < pizeros.size(); i++ )
      {
        float ph_phi = PHAngle(pizeros[i]->Phi());
        vector<float> weight = pizeros[i]->GetDecayWeights();
        for( unsigned int j = 0; j < tracks.size(); j++ )
        {
          float trk_phi = PHAngle(tracks[j]->Phi());
          float dphifold = CalculateFoldedDphi(trk_phi,ph_phi);

          h3_dphi_pi0_mix->Fill(pizeros[i]->Pt(), tracks[j]->Pt(), dphifold);
          if( pizeros[i]->IsIso() ) h3_dphi_pi0_mix_iso->Fill(pizeros[i]->Pt(), tracks[j]->Pt(), dphifold);
          for(unsigned int ipw=0;ipw<weight.size();ipw++){
            if(weight[ipw]>0) {
              h2_dphi_dec_mix[ipw]->Fill(dphifold,tracks[j]->Pt(),weight[ipw]);
              if( pizeros[i]->IsIso() ) h2_dphi_dec_mix_iso[ipw]->Fill(dphifold,tracks[j]->Pt(),weight[ipw]);
            }
          }
        }
      }
      for( unsigned int i = 0; i < clusters.size(); i++ )
      {
        float ph_phi = PHAngle(clusters[i]->Phi());
        for( unsigned int j = 0; j < tracks.size(); j++ )
        {
          float trk_phi = PHAngle(tracks[j]->Phi());
          float dphifold = CalculateFoldedDphi(trk_phi,ph_phi);

          h3_dphi_mix->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
          if( clusters[i]->IsIso() ) h3_dphi_mix_iso->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
        }
      }

      ClearVector(tracks);
      if(pooldepth == size) {/*cout<<"Mixed enough! pooldepth = "<<pooldepth<<"; nvert_fg = "<<nvert_fg<<"; ncent_fg = "<<ncent_fg<<endl; */break;}

      //make sure making NMIX pairs
      j = CheckPool(nenpart,j,pooldepth,size,nloop);
    }
    ClearVector(clusters);
    ClearVector(pizeros);
  }
}

int hijing_analysis::CheckPool(int nenpart, int j, int pooldepth, int size, int& nloop)
{
  // cout<<"in CheckPool: before nloop = "<<nloop<<endl;
  if(j == nenpart-1){
    if(pooldepth < size){
      //cout<<"pool not deep enough."<<endl;
      j = -1;
      nloop++;
      //cout<<"after nloop = "<<nloop<<endl;
    }
  }
  return j;
}

void hijing_analysis::MakeClusterObject(ACluster* aclus, float pt, float phi, float eta, float e)
{
  aclus->SetPtEtaPhiE(pt, eta, phi, e);
}

void hijing_analysis::MakePi0Object(APiZero* api0, float pt, float phi, float eta, float e)
{
  api0->SetPtEtaPhiE(pt, eta, phi, e);
  api0->Daughter1()->SetPtEtaPhiE(pt, eta, phi, e);
  api0->Daughter2()->SetPxPyPzE(0,0,0,0);
}

void hijing_analysis::MakeTrackObject(ATrack* atrk, float pt, float phi, float eta, float e)
{
  atrk->SetPtEtaPhiE(pt, eta, phi, e);
}

int hijing_analysis::GetCentrality(int mult)
{
  if( mult <= 18 ) return 3; // 60-100%
  if( mult > 18 && mult <= 28 ) return 2; // 40-60%
  if( mult > 28 && mult <=41 ) return 1; // 20-40%
  if( mult > 41 ) return 0; // 0-20%
  return 0;
}

void hijing_analysis::Init1DHisto(TH1F* &h1, string name, string xtitle, int nxbin, double xmin, double xmax)
{
  manager->registerHisto(name.c_str(),h1 = new TH1F(name.c_str(), name.c_str(), nxbin, xmin, xmax));
  h1->Sumw2();
  h1->GetXaxis()->SetTitle(xtitle.c_str());
}

void hijing_analysis::Init2DHisto(TH2F* &h2, string name, string xtitle, int nxbin, double xmin, double xmax, string ytitle, int nybin, double ymin, double ymax)
{
  manager->registerHisto(name.c_str(), h2 = new TH2F(name.c_str(), name.c_str(), nxbin, xmin, xmax, nybin, ymin, ymax));
  h2->Sumw2();
  h2->GetXaxis()->SetTitle(xtitle.c_str());
  h2->GetYaxis()->SetTitle(ytitle.c_str());
}

void hijing_analysis::Init3DHisto(TH3F* &h3, string name, string xtitle, int nxbin, double xmin, double xmax, string ytitle, int nybin, double ymin, double ymax, string ztitle, int nzbin, double zmin, double zmax)
{
  manager->registerHisto(name.c_str(), h3 = new TH3F(name.c_str(), name.c_str(), nxbin, xmin, xmax, nybin, ymin, ymax, nzbin, zmin, zmax));
  h3->Sumw2();
  h3->GetXaxis()->SetTitle(xtitle.c_str());
  h3->GetYaxis()->SetTitle(ytitle.c_str());
  h3->GetZaxis()->SetTitle(ztitle.c_str());
}

bool hijing_analysis::MakeCluster(HepMC::GenParticle* p, ACluster* clus)
{
  const HepMC::FourVector& mom_vector = p->momentum();
  
  clus->SetPxPyPzE(mom_vector.px(),mom_vector.py(),mom_vector.pz(),mom_vector.e());
  
  if( clus->E() < _MinTrigPt/2.0 ) return false;
  if( fabs(clus->Eta()) > _MaxEta ) return false;
  if( OutsideAcceptance(clus->Phi()) ) return false;
  int pbsc_pbgl = 0;
  if( clus->Phi() > 15.0*PI/16.0 ) pbsc_pbgl = 1;
  ApplyEnergyResolution(clus, pbsc_pbgl);
  clus->SetTag(false);
  
  // If cluster comes from hadronic decay set tag to true, assume all other photons are direct
  //std::cout << "Adding cluster with E = " << clus->E() << " with status = " << p->status() << endl;

  return true;

  HepMC::GenVertex* vtx = p->production_vertex();
  for( HepMC::GenVertex::particles_in_const_iterator ip = vtx->particles_in_const_begin(); ip != vtx->particles_in_const_end(); ++ip )
  {
    int pid = (*ip)->pdg_id();
    if( pid==111 || pid==221 || pid==223 ) clus->SetTag(true);
    if( verbosity ) cout << "Photon with E = " << clus->E() << " and parent id = " << pid << endl;
  }
  return true;
}

bool hijing_analysis::MakeTrack(HepMC::GenParticle* p, ATrack* track)
{
  const HepMC::FourVector& mom_vector = p->momentum();
  
  track->SetPxPyPzE(mom_vector.px(),mom_vector.py(),mom_vector.pz(),mom_vector.e());
  if( track->Pt() < _MinAssocPt || track->Pt() > 15.0 ) return false;
  if( fabs(track->Eta()) > _MaxEta ) return false;
  if( OutsideAcceptance(track->Phi()) ) return false;
  
  return true;
}

bool hijing_analysis::MakePiZero(HepMC::GenParticle* p, APiZero* piz)
{
  const HepMC::FourVector& mom_vector = p->momentum();
  
  piz->SetPxPyPzE(mom_vector.px(),mom_vector.py(),mom_vector.pz(),mom_vector.e());
  piz->Daughter1()->SetPxPyPzE(mom_vector.px(),mom_vector.py(),mom_vector.pz(),mom_vector.e());

  if( piz->E() < _MinTrigPt ) return false;
  if( fabs(piz->Eta()) > _MaxEta ) return false;
  if( OutsideAcceptance(piz->Phi()) ) return false;
  int pbsc_pbgl = 0;
  if( piz->Phi() > 15.0*PI/16.0 ) pbsc_pbgl = 1;

  if( verbosity ) std::cout << "Adding piz with E = " << piz->E() << ", M = " << piz->M();
  ApplyEnergyResolution(piz, pbsc_pbgl);
  if( verbosity ) std::cout << "  smeared to E = " << piz->E() << ", M = " << piz->M();
  HepMC::GenVertex* vtx = p->end_vertex();
  if( !vtx ) {if( verbosity ) cout << " - found no decay vertex for this piz!" << endl;}
  else {
    for( HepMC::GenVertex::particles_in_const_iterator ip = vtx->particles_out_const_begin(); ip != vtx->particles_out_const_end(); ++ip )
    {
      int pid = (*ip)->pdg_id();
      if( verbosity ) cout << " - with daughter id = " << pid;
    }
    if( verbosity ) cout << endl;
  }
  
  return true;
}

bool hijing_analysis::OutsideAcceptance(double phi)
{
  // West arm acceptance
  if( -3.0*PI/16.0 < phi && phi < 5.0*PI/16.0 ) return false;
  // East arm acceptance
  if( 11.0*PI/16.0 < phi && phi < 19*PI/16.0 ) return false;
  // otherwise it is outside acceptance
  return true;
  
}

int hijing_analysis::End(PHCompositeNode *topNode){
  //calculate centrality percentiles....
  cout<<"End."<<endl;
  DoMixing(atree->_ttrig, atree->_tpart, NMIX);
  cout<<"finish mixing."<<endl;
  manager->dumpHistos(outfile);
  
  return 0;
}



