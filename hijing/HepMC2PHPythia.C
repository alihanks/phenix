#include <HepMC2PHPythia.h>
#include <cstring>
#include <PHString.h>
#include <PHNode.h>
#include <PHPythiaContainer.h>
#include <PHHijingHeaderV4.h>
#include <getClass.h>
#include <PHCompositeNode.h>
#include <TMCParticle.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
#include <getClass.h>

using namespace std;
using namespace findNode;

HepMC2PHPythia::HepMC2PHPythia(const char* _pythianame){
   pythianame = _pythianame;
   embed = false;
}

int HepMC2PHPythia::Init(PHCompositeNode *topNode){
   phhijing = new PHPythiaContainerV2();
   phhijing->Reset();
   phhijingheader = new PHHijingHeaderV4();
   phhijingheader->Reset();
   /*
   PHCompositeNode *dstNode;
   PHNodeIterator iter(topNode);
   dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
   if(!dstNode){
      cout << PHWHERE << "no DST node" << endl;
      return -1;
   }*/
   topNode->print();
   PHIODataNode<PHPythiaContainerV2> *pythianode = new PHIODataNode<PHPythiaContainerV2>(phhijing,"PHHijing","PHPythia");
   topNode->addNode(pythianode);
   PHIODataNode<PHHijingHeaderV4> *headernode = new PHIODataNode<PHHijingHeaderV4>(phhijingheader,"PHHijingHeader","PHHijingHeader");
   topNode->addNode(headernode);
   topNode->print();

   return 0;
}

int HepMC2PHPythia::process_event(PHCompositeNode *topNode){
   if(verbosity) cout << PHWHERE << "new event" << endl;
   phhijing->Reset();
   phhijingheader->Reset();
   PHNodeIterator iter(topNode);
   PHDataNode<HepMC::GenEvent> *HepMCNode = dynamic_cast<PHDataNode<HepMC::GenEvent> *>(iter.findFirst("PHDataNode","HEPMC"));
   HepMC::GenEvent *evt = NULL;
   if(HepMCNode) evt = HepMCNode->getData();
   if(!evt)cout << "no HepMC Event found!" << endl;

   HepMC::HeavyIon *gl = evt->heavy_ion();

   double psi_2 = gl->event_plane_angle();
   psi_2 = atan2(sin(2*psi_2),cos(2*psi_2))/2.0;

   phhijingheader->SetRP(psi_2);

   int npart_proj = gl->Npart_proj();
   int npart_tar = gl->Npart_targ();

   phhijingheader->SetNt(npart_tar);
   phhijingheader->SetNp(npart_proj);

   int ncoll = gl->Ncoll();
   phhijingheader->SetNbinary(ncoll);

   float b = gl->impact_parameter();
   phhijingheader->SetBimpact(b);




   for(HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p){
      if( (*p)->end_vertex() || (*p)->status()!=1 )continue;
      int id = (*p)->pdg_id();
      const HepMC::FourVector& mom_vector = (*p)->momentum();
      TMCParticle *part = new TMCParticle();
      part->SetPx(mom_vector.px());
      part->SetPy(mom_vector.py());
      part->SetPz(mom_vector.pz());
      part->SetEnergy(mom_vector.e());
      part->SetKF(id);
      part->SetKS(1);
      phhijing->addParticle(*part);
      delete part;
   }

   if(!embed)   return 0;
   if(embed)phpythia = getClass<PHPythiaContainer>(topNode,"PHPythia");
   if(embed && !phpythia){
      cout << PHWHERE << "no PHPythia node" << endl;
      return -1;
   }

   int size = phpythia->size();
   for(int i=0; i<size; i++){
      TMCParticle *part = (TMCParticle *)phpythia->getParticle(i);
      phhijing->addParticle(*part);
   }

   if(verbosity>10)cout << PHWHERE <<  "number of particles: " << phhijing->size() << endl;
   return 0;

}
