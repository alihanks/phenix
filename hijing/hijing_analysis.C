#include <hijing_analysis.h>

#include <THmulf.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <iostream>
#include <Fun4AllHistoManager.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <PHAngle.h>
#include <Fun4AllReturnCodes.h>
#include <PHHepMCGenEvent.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <ACluster.h>
#include <ATrack.h>
#include <APiZero.h>

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
}

hijing_analysis::~hijing_analysis(){
   delete manager;
}
int hijing_analysis::Init(PHCompositeNode *topNode){

  cout << "Initializing histograms" << endl;
  nevents = 0;
  manager = new Fun4AllHistoManager("hijing_analysis");

  string name;
  
  name = "hmult";
  Init1DHisto(hmult, name.c_str(),"multiplicity",100,0,100.0);

  name = "h1_mass";
  Init1DHisto(h1_mass, name.c_str(),"mass [GeV/c^2]",500,0,1.0);
  
  name = "h1_trigger_pt";
  Init1DHisto(h1_trigger_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trigger_pi0_pt";
  Init1DHisto(h1_trigger_pi0_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trigger_dir_pt";
  Init1DHisto(h1_trigger_dir_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trigger_iso_pt";
  Init1DHisto(h1_trigger_iso_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trigger_iso_pi0_pt";
  Init1DHisto(h1_trigger_iso_pi0_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h1_trigger_iso_dir_pt";
  Init1DHisto(h1_trigger_iso_dir_pt, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0);
  
  name = "h3_dphi";
  Init3DHisto(h3_dphi, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
  
  name = "h3_dphi_iso";
  Init3DHisto(h3_dphi_iso, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
  
  name = "h3_dphi_dir";
  Init3DHisto(h3_dphi_dir, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
  
  name = "h3_dphi_pi0";
  Init3DHisto(h3_dphi_pi0, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
  
  name = "h3_dphi_iso_dir";
  Init3DHisto(h3_dphi_iso_dir, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
  
  name = "h3_dphi_iso_pi0";
  Init3DHisto(h3_dphi_iso_pi0, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
  
  name = "h2_cluster_dphi_dR";
  Init2DHisto(h2_cluster_wdR, name.c_str(), "#Delta R ", 200, 0.0, 2.0, "#Delta #phi [rad]", 60, 0.0, PI);
  
  name = "h2_cluster_pi0_dphi_dR";
  Init2DHisto(h2_cluster_pi0_wdR, name.c_str(), "#Delta R ", 200, 0.0, 2.0, "#Delta #phi [rad]", 60, 0.0, PI);
  
  name = "h2_cluster_dir_dphi_dR";
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
  
  return 0;
}

int hijing_analysis::process_event(PHCompositeNode* topNode)
{
  if( verbosity )
    cout << "Event: " << nevents << endl;
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

  int iv = 0;
  int mult = 0;
  for(HepMC::GenEvent::vertex_const_iterator v = evt->vertices_begin(); v != evt->vertices_end(); ++v)
  {
    if( verbosity ) cout << "Checking vertex " << iv << endl;
    for( HepMC::GenVertex::particles_out_const_iterator p = (*v)->particles_out_const_begin(); p != (*v)->particles_out_const_end(); ++p )
    {
      int id = (*p)->pdg_id();
      double eta = (*p)->momentum().eta();
      if( eta > -4 && eta < -3 ) mult++;

      if( id==111 )
      {
        APiZero piz;
        if( MakePiZero((*p),&piz) ) pizeros.push_back(piz.clone());
      }
      if( (*p)->status()!=1 )continue;
      if( id==22 )
      {
        ACluster clus;
        if( MakeCluster((*p),&clus) ) {
          for( HepMC::GenVertex::particles_in_const_iterator ip = (*v)->particles_in_const_begin(); ip != (*v)->particles_in_const_end(); ++ip )
          {
            int pid = (*ip)->pdg_id();
            if( pid==111 || pid==221 || pid==223 ) clus.SetTag(true);
          if( verbosity ) cout << "Photon with E = " << clus.E() << " and parent id = " << pid << endl;
          }
          clusters.push_back(clus.clone());
        }
      }
      if( (fabs(id)==211 || fabs(id)==321 || fabs(id)==2212) )
      {
        ATrack track;
        if( MakeTrack((*p),&track) ) tracks.push_back(track.clone());
      }
    }
    iv++;
  }

  //cout << "Checking particle of status " << (*p)->status() << endl;
  
//    int ntags = 0;
//    if( clusters.size() > 1 ) {
//      for( unsigned int i = 0; i < clusters.size()-1; i++ )
//      {
//        for( unsigned int j = i+1; j < clusters.size(); j++ )
//        {
//          TLorentzVector sum = *clusters[i] + *clusters[j];
//          h1_mass->Fill(sum.M());
//          cout << "Photons tagged as " << clusters[i]->IsTagged() << ", " << clusters[j]->IsTagged() << ": E = " << clusters[i]->E() << " + " << clusters[j]->E() << " = " << sum.E() << ", pair mass = " << sum.M() << endl;
//          if( (sum.M() > 0.133 && sum.M() < 0.137) ) {
//            clusters[i]->SetTag(true);
//            clusters[j]->SetTag(true);
//            ntags++;
//          }
//        }
//      }
//    }
//   cout << "Tagged photons from " << ntags << " pizero candidates!" << endl;
  
  int cent_bin = GetCentrality(mult);
  if( cent_bin != Centrality && Centrality>=0 ) return 0;
  if( verbosity )
    cout << "Event multiplicity = " << mult << endl;
  hmult->Fill((float)mult);
  if( verbosity ) cout << "Looping over " << pizeros.size() << " pizeros" << endl;
  for( unsigned int i = 0; i < pizeros.size(); i++ )
  {
    h1_trigger_pi0_pt->Fill(pizeros[i]->Pt());
    SetIso(pizeros[i],tracks,clusters,Rcut,h3_cluster_pi0_dR,h3_cluster_pi0_etot,h2_cluster_pi0_wdR,h2_cluster_pi0_etot);
    if( pizeros[i]->IsIso() ) h1_trigger_iso_pi0_pt->Fill(pizeros[i]->Pt());
    //cout << "Set iso<" << Rcut << " cut for piz to " << pizeros[i]->IsIso() << endl;
  }
  
  if( verbosity ) cout << "Looping over " << clusters.size() << " clusters" << endl;
  for( unsigned int i = 0; i < clusters.size(); i++ )
  {
    if( !clusters[i]->IsTagged() ) {
      h1_trigger_dir_pt->Fill(clusters[i]->Pt());
      SetIso(clusters[i],tracks,clusters,Rcut,h3_cluster_dir_dR,h3_cluster_dir_etot,h2_cluster_dir_wdR,h2_cluster_dir_etot);
      if( clusters[i]->IsIso() ) h1_trigger_iso_dir_pt->Fill(clusters[i]->Pt());
    }
    h1_trigger_pt->Fill(clusters[i]->Pt());
    SetIso(clusters[i],tracks,clusters,Rcut,h3_cluster_dR,h3_cluster_etot,h2_cluster_wdR,h2_cluster_etot);
    if( clusters[i]->IsIso() ) h1_trigger_iso_pt->Fill(clusters[i]->Pt());
  }

  if( verbosity ) cout << "Looping over " << tracks.size() << " tracks" << endl;
  for( unsigned int i = 0; i < clusters.size(); i++ )
  {
    float ph_phi = PHAngle(clusters[i]->Phi());
    for( unsigned int j = 0; j < tracks.size(); j++ )
    {
      float trk_phi = PHAngle(tracks[j]->Phi());
      float dphifold = CalculateFoldedDphi(trk_phi,ph_phi);

      h3_dphi->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
      if( clusters[i]->IsIso() ) h3_dphi_iso->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
      if( !clusters[i]->IsTagged() ) h3_dphi_dir->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
      if( !clusters[i]->IsTagged() && clusters[i]->IsIso() ) h3_dphi_iso_dir->Fill(clusters[i]->Pt(), tracks[j]->Pt(), dphifold);
    }
  }
  for( unsigned int i = 0; i < pizeros.size(); i++ )
  {
    float ph_phi = PHAngle(pizeros[i]->Phi());
    for( unsigned int j = 0; j < tracks.size(); j++ )
    {
      float trk_phi = PHAngle(tracks[j]->Phi());
      float dphifold = CalculateFoldedDphi(trk_phi,ph_phi);

      h3_dphi_pi0->Fill(pizeros[i]->Pt(), tracks[j]->Pt(), dphifold);
      if( pizeros[i]->IsIso() ) h3_dphi_iso_pi0->Fill(pizeros[i]->Pt(), tracks[j]->Pt(), dphifold);
    }
  }
  
  ClearVector(clusters);
  ClearVector(tracks);
  ClearVector(pizeros);
  return 0;
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
  
  clus->SetPx(mom_vector.px());
  clus->SetPy(mom_vector.py());
  clus->SetPz(mom_vector.pz());
  clus->SetE(mom_vector.e());
  
  if( clus->E() < _MinTrigPt ) return false;
  if( fabs(clus->Eta()) > _MaxEta ) return false;
  if( OutsideAcceptance(clus->Phi()) ) return false;
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
  
  track->SetPx(mom_vector.px());
  track->SetPy(mom_vector.py());
  track->SetPz(mom_vector.pz());
  track->SetE(mom_vector.e());
  if( track->Pt() < _MinAssocPt && track->Pt() < 99 ) return false;
  if( fabs(track->Eta()) > _MaxEta ) return false;
  if( OutsideAcceptance(track->Phi()) ) return false;
  
  return true;
}

bool hijing_analysis::MakePiZero(HepMC::GenParticle* p, APiZero* piz)
{
  const HepMC::FourVector& mom_vector = p->momentum();
  
  piz->SetPx(mom_vector.px());
  piz->SetPy(mom_vector.py());
  piz->SetPz(mom_vector.pz());
  piz->SetE(mom_vector.e());
  piz->Daughter1()->SetPx(mom_vector.px());
  piz->Daughter1()->SetPy(mom_vector.py());
  piz->Daughter1()->SetPz(mom_vector.pz());
  piz->Daughter1()->SetE(mom_vector.e());
  if( piz->E() < _MinTrigPt ) return false;
  if( fabs(piz->Eta()) > _MaxEta ) return false;
  if( OutsideAcceptance(piz->Phi()) ) return false;
  if( verbosity ) std::cout << "Adding piz with E = " << piz->E() << ", M = " << piz->M();
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
  manager->dumpHistos(outfile);
  
  return 0;
}



