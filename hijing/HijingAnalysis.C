#include <HijingAnalysis.h>

#include <THmulf.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <iostream>
#include <math.h>
#include <Fun4AllHistoManager.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <Fun4AllReturnCodes.h>


using namespace std;
using namespace findNode;

const double pi = acos(-1.0);

HijingAnalysis::HijingAnalysis(const char* output, const char* name) 
{
  ThisName  = name;
  verbosity = 0;
  outfile = output;
  embedding = false;
  d_pythia = NULL;
}
 

HijingAnalysis::~HijingAnalysis(){
   delete manager;
}
int HijingAnalysis::Init(PHCompositeNode *topNode){

  nevents = 0;
  manager = new Fun4AllHistoManager("HijingAnalysis");

  return 0;
}

int HijingAnalysis::process_event(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHDataNode<HepMC::GenEvent> *HepMCNode = dynamic_cast<PHDataNode<HepMC::GenEvent> *>(iter.findFirst("PHDataNode","HEPMC"));
  HepMC::GenEvent* evt = NULL;
  if(HepMCNode) evt = HepMCNode->getData();

  for(HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p)
  {
    if( (*p)->end_vertex() || (*p)->status()!=1 ) continue;
    int id = (*p)->pdg_id();
    const HepMC::FourVector& mom_vector = (*p)->momentum();

    if( id == 22 ) {
      ACluster clus;
      MakeParticle(mom_vector,&clus);
      
      all_clusters.push_back(clus.clone());
      
      if( clus.Pt() < 5.0 ) continue;
      bool is_pi0_daughter = false;
      bool is_daughter = false;
      HepMC::GenVertex* vertex = (*p)->production_vertex();
      for(HepMC::GenVertex::articles_in_const_iterator v = vertex->particles_in_const_begin(); v != vertex->particles_in_const_end(); ++v)
      {
        int pid = (*v)->pdg_id();
        if( pid == 111 ) is_pi0_daughter = true;
        if( pid > 100 ) is_daughter = true;
      }
      if( is_pi0_daughter ) pi0_clusters.push_back(clus.clone());
      if( !is_daughter ) direct_clusters.push_back(clus.clone());
    }
    if( (fabs(id)==211 || fabs(id)==321 || fabs(id)==2212) ) {
      ATrack track;
      MakeParticle(mom_vector,&track);
      
      all_tracks.push_back(track.clone());
      
      if( track.Pt() < 0.5 ) continue;
      tracks.push_back(track.clone());
    }
  }
}

void HijingAnalysis::MakeParticle(HepMC::FourVector mom_vector, AParticle* part)
{
  part->SetPx(mom_vector.px());
  part->SetPy(mom_vector.py());
  part->SetPz(mom_vector.pz());
  part->SetE(mom_vector.e());
}

int HijingAnalysis::process_event(PHCompositeNode *topNode){
   d_part = getClass<PHPythiaContainerV2>(topNode,"PHHijing");
   d_global = getClass<PHHijingHeaderV2>(topNode,"PHHijingHeader");
   if(embedding){
      d_pythia = getClass<PHPythiaContainerV2>(topNode,"PHPythia");
      if(!d_pythia){
	 cout << PHWHERE << "pythia node not found!" << endl;
	 return ABORTEVENT;
      }
   }

   if(!d_part)cout << " no pythia particles! " << endl;
   if(!d_global)cout << "no global information" << endl;
   int mult = 0;
   int mult_Au = 0;
   hnpart->Fill(d_global->GetNp()+d_global->GetNt());
   unsigned int size = d_part->size();
   double et_sum = 0.0;
   double sin_sum[N_ORDERS][N_RXPDET], cos_sum[N_ORDERS][N_RXPDET];
   for(int i=0; i<N_ORDERS; i++){
      for(int j =0; j<N_RXPDET; j++){
	 sin_sum[i][j] = 0.0;
	 cos_sum[i][j] = 0.0;
      }
   }

   for(int i=0; i<N_PT; i++){
      pion_phis[i].clear();
   }

   for(unsigned int i=0; i<size; i++){
      TMCParticle *part = d_part->getParticle(i);
      int pid = part->GetKF();
      double pt = part->GetPx()*part->GetPx() + part->GetPy()*part->GetPy();
      pt = sqrt(pt);
      double theta = atan2(pt,part->GetPz());
      double eta = -log(tan(theta/2.0));


      if((fabs(pid) == 211 || pid == 111) && fabs(eta)<0.35 && pt > pt_ranges[0][0]){
	 double phi = atan2(part->GetPy(),part->GetPx());
	 for(int ii=0; ii<N_PT; ii++){
	    if(pt > pt_ranges[0][ii] && pt < pt_ranges[1][ii])
	       pion_phis[ii].push_back(phi);
	 }
      }

      int status = part->GetKS();
      if((fabs(pid)==211 || fabs(pid)==321 || fabs(pid)==2212) && status ==1){
	 mult++;
	 hdNdeta->Fill(eta);
	 double phi = atan2(part->GetPy(),part->GetPx());
	 if(eta > -4 && eta < -3){
	    for(int ii=0; ii<size; ii++){
	       if(ii==i)continue; //not the same particle
	       TMCParticle *part2 = d_part->getParticle(ii);
	       int status2 = part2->GetKS();
	       int pid2 = fabs(part2->GetKF());
	       if(status2 != 1) continue;
	       if(pid2 != 211 && pid2 != 321 && pid2 != 2212)continue;
	       double pt2 = part2->GetPx()*part2->GetPx() + part2->GetPy()*part2->GetPy();
	       pt2 = sqrt(pt2);
	       if(pt2<0.1)continue;
	       double theta2 = atan2(pt2,part2->GetPz());
	       double eta2 = -log(tan(theta2/2.0));
	       if(eta2 > -4 && eta2 < -3){
		  double phi2 = atan2(part2->GetPy(),part2->GetPx());
		  double dphi = phi - phi2;
		  while(dphi > 1.5*pi)dphi -= 2.0*pi;
		  while(dphi < -0.5*pi)dphi += 2.0*pi;
		  hdphi_BBC->Fill(dphi);
	       }
	    }


	    mult_Au++;
   	    for(int ii=0; ii<N_ORDERS; ii++){
   	       sin_sum[ii][0] += sin((double)(ii+1)*phi);
   	       cos_sum[ii][0] += cos((double)(ii+1)*phi);
   	    }
	    for(int ii=0; ii<size; ii++){
	       TMCParticle *part2 = d_part->getParticle(ii);
	       int pid2 = part2->GetKF();
      
	       double pt2 = part2->GetPx()*part2->GetPx() + part2->GetPy()*part2->GetPy();
	       pt2 = sqrt(pt2);
	       if(pt2<0.3)continue;
	       double theta2 = atan2(pt2,part2->GetPz());
	       double eta2 = -log(tan(theta2/2.0));
	       if(fabs(eta2)>0.35)continue;
	       double phi2 = atan2(part2->GetPy(),part2->GetPx());
	       double dphi = phi - phi2;
	       if(dphi > 1.5*pi) dphi -= 2.0*pi;
	       if(dphi < -0.5*pi) dphi += 2.0*pi;
	       h2pc->Fill(dphi);
	    }

	 }
	 if(eta > -2.8 && eta < -1.0){
   	    for(int ii=0; ii<N_ORDERS; ii++){
   	       sin_sum[ii][2] += sin((double)(ii+1)*phi);
   	       cos_sum[ii][2] += cos((double)(ii+1)*phi);
   	    }
	 }
	 if(fabs(eta) < 0.35){
   	    for(int ii=0; ii<N_ORDERS; ii++){
   	       sin_sum[ii][3] += sin((double)(ii+1)*phi);
   	       cos_sum[ii][3] += cos((double)(ii+1)*phi);
   	    }
	 }

      }
      if(pid == 22 && eta > -3.9 && eta < -3.1){
	 double phi = atan2(part->GetPy(),part->GetPx());
	 double et = part->GetEnergy()*sin(theta);
	 for(int ii=0; ii<N_ORDERS; ii++){
   	    sin_sum[ii][1] += sin((double)(ii+1)*phi*et);
   	    cos_sum[ii][1] += cos((double)(ii+1)*phi*et);
	 }
      }
      double et = sqrt(pt*pt + part->GetMass()*part->GetMass());
      if(fabs(eta)<1.0){
	 et_sum += et;
      }
   }

   et_sum /= 63.0;
   et_sum /= 20.0;
   het_npart->Fill(d_global->GetNp() + d_global->GetNt(),et_sum);
   hmult->Fill(mult);
   hmult_Au->Fill(mult_Au);
   nevents++;
   bool isCentral = false;
   bool isPeripheral = false;
   if(mult_Au >= 17)isCentral = true; //corresponds to the top 20%
   if(mult_Au <= 5) isPeripheral = true;
   
   double psi[N_ORDERS][N_RXPDET];
   if(isCentral){
      for(int ii=0; ii<N_ORDERS; ii++){
      	 for(int jj=0; jj<N_RXPDET; jj++){
   	    double order = (double)(ii + 1);
   	    psi[ii][jj] = atan2(sin_sum[ii][jj],cos_sum[ii][jj]) / order;
   	    hpsi[ii][jj]->Fill(psi[ii][jj]);
	    for(int kk=0; kk<N_PT; kk++){
	    
	       vector<double>::iterator iv;
	       for(iv=pion_phis[kk].begin(); iv!=pion_phis[kk].end(); ++iv){
		  double dphi = psi[ii][jj] - *iv;
		  /*
		  if(dphi > pi) dphi -= pi;
		  if(dphi < -1.0*pi) dphi += pi;*/
		  if(fabs(dphi) > pi/2.0){
		     dphi = pi - fabs(dphi);
		  }
		  hdphi[ii][jj][kk]->Fill(fabs(dphi));
		  if(fabs(dphi) > pi/2.0 && ii==1) cout << psi[ii][jj] << " " << *iv << " " << dphi << endl;
	       }
	    }
      	 }
	 for(int jj=0; jj<N_RXPDET; jj++){
	    for(int kk=0; kk<N_RXPDET; kk++){
	       double dpsi = psi[ii][jj] - psi[ii][kk];
	       double n = (double)(ii+1);
	       while(dpsi > pi) dpsi -= 2.0*pi;
	       while(dpsi < -pi) dpsi += 2.0*pi;
	       if(fabs(dpsi) > pi/n) dpsi = pi - fabs(dpsi);
	       if(fabs(dpsi) > pi/n)
		  cout << dpsi << " " << n << " " << psi[ii][jj] << " " << psi[ii][kk] << endl;
	       hdpsi[ii][jj][kk]->Fill(fabs(dpsi));
	       hcosndphi[ii][jj][kk]->Fill(cos(n*dpsi));
	    }
	 }
      }
      hpsi_corr->Fill(fabs(psi[0][0]),fabs(psi[1][0]));
      hpsi_corr_2_rap->Fill(fabs(psi[1][0]),fabs(psi[1][3]));
   }

   return 0;
}

int HijingAnalysis::End(PHCompositeNode *topNode){
   //calculate centrality percentiles....
   double total = hmult_Au->Integral();
   double tmp = 0;
   //0-5%
   for(int i = hmult_Au->GetNbinsX(); i > 0; i--){
      tmp += hmult_Au->GetBinContent(i);
      if( tmp > 0.05*total){
	 cout << "top 5%: " << hmult_Au->GetBinCenter(i) << endl;
	 break;
      }
   }

   /*
   tmp = 0;
   //50-100%
   for(int i=1; i<= hmult_Au->GetNbinsX(); i++){
      tmp += hmult_Au->GetBinContent(i);
      if(tmp > 0.5*total){
	 cout << "bottom 50%: " << hmult_Au->GetBinCenter(i) << endl;
	 break;
      }
   }
   for(int i=0; i<4; i++){
      if(ntrig_cent)hdphi_cent[i]->Scale(1.0/ntrig_cent);
      if(ntrig_per)hdphi_per[i]->Scale(1.0/ntrig_per);
      hdphi_diff[i]->Add(hdphi_cent[i],hdphi_per[i],1.0,-1.0);
   }*/
   hdNdeta->Scale(1.0/nevents);
   manager->dumpHistos(outfile);
      
   return 0;
}



