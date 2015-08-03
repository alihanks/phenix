#ifndef __HepMC2PHPythia_H__
#define __HepMC2PHPythia_H__

#include <PHPythiaContainerV2.h>
#include <SubsysReco.h>

class PHCompositeNode;
class PHPythiaContainerV2;
class PHHijingHeaderV4;
class PHPythiaContainer;

class HepMC2PHPythia : public SubsysReco
{
   public:
      HepMC2PHPythia(const char* _pythianame = "");
      ~HepMC2PHPythia() {}
      int Init(PHCompositeNode *topNode);
      int process_event(PHCompositeNode *topNode);
      void IsEmbedding(bool val) { embed = val; }

   private:
      const char* pythianame;
      PHPythiaContainerV2 *phhijing;
      PHPythiaContainer *phpythia;
      PHHijingHeaderV4 *phhijingheader;
      int events;
      bool embed;
 
};

#endif
