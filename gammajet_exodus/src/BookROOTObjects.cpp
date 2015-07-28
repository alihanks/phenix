//-----------------------------------------------------------------------------
//
//  Book ROOT objects declared in DeclareROOTObjects
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TNtuple.h>

#define  INCLUDEFLAG extern
#include "DeclareROOTObjects.h"

void BookROOTObjects()
{
  char dummy[128];
  double ptmax=30.0;
  int binnum=300;

  //  primaries = new TNtuple("primaries","primary-particle ntuple",
  //		   "mass:weight:p:th:ph:pid");
  //  pairs = new TNtuple("pairs","pi0 photon-pair ntuple",
  //	       "pte:accept1:weight1:ptp:accept2:weight2");
  //  singles = new TNtuple("singles","single-electron ntuple",
  //	       "weight:charge:sece:pe:the:phe:pide:p:th:ph:pid");

  //  hpairinvm = new TH1D("hpairinvm","Mass",150,0.0,1.5);

  ptpivsptgam = new TH2F("ptpivsptgam","ptpivsptgam",400,0,20,400,0,20);

  eta_trueptvsrecopt = new TH2F("eta_trueptvsrecopt","eta_trueptvsrecopt",200,0,20,400,0,20);
  pi0_trueptvsrecopt = new TH2F("pi0_trueptvsrecopt","pi0_trueptvsrecopt",200,0,20,400,0,20);
  realpi0_trueptvsrecopt = new TH2F("realpi0_trueptvsrecopt","realpi0_trueptvsrecopt",200,0,20,400,0,20);
  direct_trueptvsrecopt = new TH2F("direct_trueptvsrecopt","direct_trueptvsrecopt",200,0,20,400,0,20);

  eta_leadptvsrecopt = new TH2F("eta_leadptvsrecopt","eta_leadptvsrecopt",200,0,20,400,0,20);
  pi0_leadptvsrecopt = new TH2F("pi0_leadptvsrecopt","pi0_leadptvsrecopt",200,0,20,400,0,20);
  realpi0_leadptvsrecopt = new TH2F("realpi0_leadptvsrecopt","realpi0_leadptvsrecopt",200,0,20,400,0,20);
  direct_leadptvsrecopt = new TH2F("direct_leadptvsrecopt","direct_leadptvsrecopt",200,0,20,400,0,20);


  ZVERTEX=new TH1F("ZVERTEX","ZVERTEX",100,-50,50);  
  //  PTpi_ZEMCpi_PTgam = new TH3F("PTpi_ZEMCpi_PTgam","PTpi_ZEMCpi_PTgam",400,0,20,360,-180,180,400,0,20);
  //RECONPI_ZEMC_PT = new TH2F("RECONPI_ZEMC_PT","RECONPI_ZEMC_PT",400,-200,200, 20,0,20);

  //reduced #of bins in 3d histos
  
  PTpi_ZEMCpi_PTgam = new TH3F("PTpi_ZEMCpi_PTgam","PTpi_ZEMCpi_PTgam",400,0,20,35,-175,175,400,0,20);
  PTpi_ZEMCpi_PTgam_MISS = new TH3F("PTpi_ZEMCpi_PTgam_MISS","PTpi_ZEMCpi_PTgam_MISS",400,0,20,35,-175,175,400,0,20);

  //Made 2D
  //PTpi_ZEMCpi_PTgam_MISS_PI0 = new TH3F("PTpi_ZEMCpi_PTgam_MISS_PI0","PTpi_ZEMCpi_PTgam_MISS_PI0",400,0,20,35,-175,175,400,0,20);
  //PTpi_ZEMCpi_PTgam_MISS_ETA = new TH3F("PTpi_ZEMCpi_PTgam_MISS_ETA","PTpi_ZEMCpi_PTgam_MISS_ETA",400,0,20,35,-175,175,400,0,20);

  PTpi_PTgam_MISS_PI0 = new TH2F("PTpi_PTgam_MISS_PI0","PTpi_PTgam_MISS_PI0",400,0,20,400,0,20);
  PTpi_PTgam_MISS_ETA = new TH2F("PTpi_PTgam_MISS_ETA","PTpi_PTgam_MISS_ETA",400,0,20,400,0,20);

  RECONPI_ZEMC_PT = new TH2F("RECONPI_ZEMC_PT","RECONPI_ZEMC_PT",35,-175,175, 400,0,20);
  
  //RECONPI_ZEMC_RECONPT_REALPT= new TH3F("ZEMCpi_PTreco_PTreal","ZEMCpi_PTreco_PTreal",35,-175,175,400,0,20,400,0,20);

  //RECONPI_ZEMC_PTTRUE_PTRECO_false= new TH2F("RECONPI_ZEMC_PTTRUE_PTRECO_false","RECONPI_ZEMC_PTTRUE_PTRECO_false",35,-175,175,400,0,20);

  // PTpi_ZEMCpi_PTgam_false= new TH3F("PTpi_ZEMCpi_PTgam_false","PTpi_ZEMCpi_PTgam_false",35,-175,175,400,0,20,400,0,20);

  PTpi_PTgam_false= new TH2F("PTpi_PTgam_false","PTpi_PTgam_false",400,0,20,400,0,20);

  PTpi_PTgam_falsetag= new TH2F("PTpi_PTgam_falsetag","PTpi_PTgam_falsetag",400,0,20,400,0,20);
  PTpi_PTgam_falsetag_eta= new TH2F("PTpi_PTgam_falsetag_eta","PTpi_PTgam_falsetag_eta",400,0,20,400,0,20);
  PTpi_PTgam_falsetag_dir= new TH2F("PTpi_PTgam_falsetag_dir","PTpi_PTgam_falsetag_dir",400,0,20,400,0,20);
  PTpi_PTgam_truetag= new TH2F("PTpi_PTgam_truetag","PTpi_PTgam_truetag",400,0,20,400,0,20);

 RECONPI_ZEMC_PTTRUE_PTRECO_real= new TH3F("RECONPI_ZEMC_PTTRUE_PTRECO_real","RECONPI_ZEMC_PTTRUE_PTRECO_real",35,-175,175,400,0,20,400,0,20);

 /*WRONG
 PTpi_ZEMCpi_PTgam_TRUE= new TH3F("PTpi_ZEMCpi_PTgam_TRUE","PTpi_ZEMCpi_PTgam_TRUE",35,-175,175,400,0,20,400,0,20);
 PTpi_ZEMCpi_PTgam_PI0= new TH3F("PTpi_ZEMCpi_PTgam_PI0","PTpi_ZEMCpi_PTgam_PI0",35,-175,175,400,0,20,400,0,20);
 PTpi_ZEMCpi_PTgam_ETA= new TH3F("PTpi_ZEMCpi_PTgam_ETA","PTpi_ZEMCpi_PTgam_ETA",35,-175,175,400,0,20,400,0,20);
 PTpi_ZEMCpi_PTgam_DIR= new TH3F("PTpi_ZEMCpi_PTgam_DIR","PTpi_ZEMCpi_PTgam_DIR",35,-175,175,400,0,20,400,0,20);
 */

 //made 3D 2D
 /*
 PTpi_ZEMCpi_PTgam_TRUE= new TH3F("PTpi_ZEMCpi_PTgam_TRUE","PTpi_ZEMCpi_PTgam_TRUE",400,0,20,35,-175,175,400,0,20);
 PTpi_ZEMCpi_PTgam_PI0= new TH3F("PTpi_ZEMCpi_PTgam_PI0","PTpi_ZEMCpi_PTgam_PI0",400,0,20,35,-175,175,400,0,20);
 PTpi_ZEMCpi_PTgam_ETA= new TH3F("PTpi_ZEMCpi_PTgam_ETA","PTpi_ZEMCpi_PTgam_ETA",400,0,20,35,-175,175,400,0,20);
 PTpi_ZEMCpi_PTgam_DIR= new TH3F("PTpi_ZEMCpi_PTgam_DIR","PTpi_ZEMCpi_PTgam_DIR",400,0,20,35,-175,175,400,0,20);
 */

 PTpi_PTgam_TRUE= new TH2F("PTpi_PTgam_TRUE","PTpi_PTgam_TRUE",400,0,20,400,0,20);
 PTpi_PTgam_False_PI0pair= new TH2F("PTpi_PTgam_False_PI0pair","PTpi_PTgam_False_PI0pair",400,0,20,400,0,20);
 PTpi_PTgam_PI0= new TH2F("PTpi_PTgam_PI0","PTpi_PTgam_PI0",400,0,20,400,0,20);
 PTpi_PTgam_ETA= new TH2F("PTpi_PTgam_ETA","PTpi_PTgam_ETA",400,0,20,400,0,20);
 PTpi_PTgam_DIR= new TH2F("PTpi_PTgam_DIR","PTpi_PTgam_DIR",400,0,20,400,0,20);


 // RECONPI_ZEMC_PTTRUE_PTRECO_allfalse= new TH3F("RECONPI_ZEMC_PTTRUE_PTRECO_allfalse","RECONPI_ZEMC_PTTRUE_PTRECO_allfalse",35,-175,175,400,0,20,400,0,20);



  ptpivsptgam_tag = new TH2F("ptpivsptgam_tag","ptpivsptgam_tag",400,0,20,400,0,20);
  ppivspgam = new TH2F("ppivspgam","ppivspgam",400,0,20,400,0,20);
  ppivspgam_tag = new TH2F("ppivspgam_tag","ppivspgam_tag",400,0,20,400,0,20);

  realpi0pt = new TH1D("realpi0pt","real pi0 pt",binnum,0,ptmax);
  realetapt = new TH1D("realetapt","real eta pt",binnum,0,ptmax);
  realetaprimept = new TH1D("realetaprimept","real eta prime pt",binnum,0,ptmax);
  realomegapt = new TH1D("realomegapt","real omega prime pt",binnum,0,ptmax);

  gammapt   = new TH1D("gammapt","all gamma in full space",binnum,0,ptmax);
 
  //5/3/2006 for double check
  gamma_acc = new TH1D("gamma_acc","gamma in PbSc acceptance",binnum,0,ptmax);
  gamma_acc_nosm = new TH1D("gamma_acc_nosm","gamma in PbSc acceptance (nosm) ",binnum,0,ptmax);
  gamma_acc_phi = new TH1D("gamma_acc_phi","gamma in PbSc acceptance",200,-2,5);
  gamma_acc_phi_live = new TH1D("gamma_acc_phi_live","gamma in PbSc acceptance",200,-2,5);
  gamma_acc_ert_phi = new TH1D("gamma_acc_ert_phi","gamma in PbSc acc_erteptance",200,-2,5);
  //  gamma_acc_phi->Sumw2();
  gamma_ecore = new TH1D("gamma_ecore","gamma pass ecore cut",binnum,0,ptmax);

  //gamma_pair = new TH1D("gamma_pair","gammas pass various cuts",9,-0.5,8.5);

  gamma_pair[0] = new TH1D("gamma_pair_0","gamma_pair_0",binnum,0,ptmax);
  gamma_pair[1] = new TH1D("gamma_pair_1","gamma_pair_1",binnum,0,ptmax);
  gamma_pair[2] = new TH1D("gamma_pair_2","gamma_pair_2",binnum,0,ptmax);
  gamma_pair[3] = new TH1D("gamma_pair_3","gamma_pair_3",binnum,0,ptmax);
  gamma_pair[4] = new TH1D("gamma_pair_4","gamma_pair_4",binnum,0,ptmax);
  gamma_pair[5] = new TH1D("gamma_pair_5","gamma_pair_5",binnum,0,ptmax);
  gamma_pair[6] = new TH1D("gamma_pair_6","gamma_pair_6",binnum,0,ptmax);
  gamma_pair[7] = new TH1D("gamma_pair_7","gamma_pair_7",binnum,0,ptmax);

  pi_asym[0] = new TH1D("pi_asym_0","pi_asym_0",600,-30,30);
  pi_asym[1] = new TH1D("pi_asym_1","pi_asym_1",600,-30,30);

  pi0gammapt[0] = new TH1D("pi0gammapt_0","pi0gammapt_0",binnum,0,ptmax);
  pi0gammapt[1] = new TH1D("pi0gammapt_1","pi0gammapt_1",binnum,0,ptmax);
  pi0gammapt[2] = new TH1D("pi0gammapt_2","pi0gammapt_2",binnum,0,ptmax);


  invmassvspt=new TH2F("invmassvspt","inmvassvspt",200,0,20,200,0,1);
  invmassvspt->Sumw2();

    INVMASS[0] = new TH2F("INVMASS_0", "INVMASS_0",200,0,20,200,0,1);
    INVMASS[1] = new TH2F("INVMASS_1","INVMASS_1",200,0,20,200,0,1);
    INVMASS[2] = new TH2F("INVMASS_2","INVMASS_2",200,0,20,200,0,1);
    INVMASS[3] = new TH2F("INVMASS_3","INVMASS_3",200,0,20,200,0,1);
    INVMASS[4] = new TH2F("INVMASS_4","INVMASS_4",200,0,20,200,0,1);
    INVMASS[5] = new TH2F("INVMASS_5","INVMASS_5",200,0,20,200,0,1);
    INVMASS[6] = new TH2F("INVMASS_6","INVMASS_6",200,0,20,200,0,1);
    INVMASS[7] = new TH2F("INVMASS_7","INVMASS_7",200,0,20,200,0,1);

    pi_INVMASS[0] = new TH2F("pi_INVMASS_0", "pi_INVMASS_0",200,0,20,200,0,1);
    pi_INVMASS[1] = new TH2F("pi_INVMASS_1","pi_INVMASS_1",200,0,20,200,0,1);
    pi_INVMASS[2] = new TH2F("pi_INVMASS_2","pi_INVMASS_2",200,0,20,200,0,1);


    INVMASS_PBGL[0] = new TH2F("INVMASS_PBGL_0", "INVMASS_PBGL_0",200,0,20,200,0,1);
    INVMASS_PBGL[1] = new TH2F("INVMASS_PBGL_1","INVMASS_PBGL_1",200,0,20,200,0,1);
    INVMASS_PBGL[2] = new TH2F("INVMASS_PBGL_2","INVMASS_PBGL_2",200,0,20,200,0,1);
    INVMASS_PBGL[3] = new TH2F("INVMASS_PBGL_3","INVMASS_PBGL_3",200,0,20,200,0,1);
    INVMASS_PBGL[4] = new TH2F("INVMASS_PBGL_4","INVMASS_PBGL_4",200,0,20,200,0,1);
    INVMASS_PBGL[5] = new TH2F("INVMASS_PBGL_5","INVMASS_PBGL_5",200,0,20,200,0,1);
    INVMASS_PBGL[6] = new TH2F("INVMASS_PBGL_6","INVMASS_PBGL_6",200,0,20,200,0,1);
    INVMASS_PBGL[7] = new TH2F("INVMASS_PBGL_7","INVMASS_PBGL_7",200,0,20,200,0,1);


  invmass=new TH1D("invmass","inmvass",100,0,1);
  pi0gamma_acc = new TH1D("pi0gamma_acc","pi0 gamma in PbSc acceptance both",binnum,0,ptmax);
  pi0gamma_acc_cut = new TH1D("pi0gamma_acc_cut","pi0 gamma in PbSc after sepa cut",binnum,0,ptmax);

  openangle = new TH1D("openangle","open angle between two pi0 photons",1000,0,3.14);
     
  pi0gamma = new TH1D("pi0gamma","gamma from pi0",binnum,0,ptmax);
  etagamma = new TH1D("etagamma","gamma from eta",binnum,0,ptmax);
  etaprimegamma = new TH1D("etaprimegamma","gamma from eta prime",binnum,0,ptmax);
  omegagamma = new TH1D("omegagamma","gamma from omega",binnum,0,ptmax);

  beforecut = new TH1D("beforecut","pi0 open angle",binnum,0,ptmax);
  aftercut = new TH1D("aftercut","pi0 open angle",binnum,0,ptmax);

  accbeforecut = new TH1D("accbeforecut","pi0 open angle",binnum,0,ptmax);
  accaftercut = new TH1D("accaftercut","pi0 open angle",binnum,0,ptmax);

  accvar = new TH1D("accvar","accept",201,-100.5,100.5);

  cocktail = new TH1D("cocktail","cocktail",binnum,0,ptmax);
  
  tightptvspt = new TH2F("tightptvspt","tightptvspt",binnum,0,ptmax,binnum,0,ptmax);

  tightptvspt_wo = new TH2F("tightptvspt_wo","tightptvspt_wo",binnum,0,ptmax,binnum,0,ptmax);
  
  massvspt = new TH2F("massvspt","minv vs pt",100,0,0.5,binnum,0,ptmax);

  loosein = new TH1F("loosein","loosein",binnum,0,ptmax);
  looseall = new TH1F("looseall","looseall",binnum,0,ptmax);

  beforecut_photon = new TH1D("beforecut_photon","photon before open angle",binnum,0,ptmax);
  aftercut_photon = new TH1D("aftercut_photon","photon after open angle",binnum,0,ptmax);

  TRIGPT= new TH1D("TRIGPT","TRIGPT",binnum,0,ptmax);

  pi_TRIGPT= new TH1D("pi_TRIGPT","pi_TRIGPT",binnum,0,ptmax);


  TRIGPT_0= new TH1D("TRIGPT_0","TRIGPT_0",binnum,0,ptmax);
  TRIGPT_1= new TH1D("TRIGPT_1","TRIGPT_1",binnum,0,ptmax);
  TRIGPT_2= new TH1D("TRIGPT_2","TRIGPT_2",binnum,0,ptmax);
  TRIGPT_3= new TH1D("TRIGPT_3","TRIGPT_3",binnum,0,ptmax);
  TRIGPT_4= new TH1D("TRIGPT_4","TRIGPT_4",binnum,0,ptmax);


  DIRPT= new TH1D("DIRPT","DIRPT",binnum,0,ptmax);

  DIRPT_0= new TH1D("DIRPT_0","DIRPT_0",binnum,0,ptmax);
  DIRPT_1= new TH1D("DIRPT_1","DIRPT_1",binnum,0,ptmax);
  DIRPT_2= new TH1D("DIRPT_2","DIRPT_2",binnum,0,ptmax);
  DIRPT_3= new TH1D("DIRPT_3","DIRPT_3",binnum,0,ptmax);
  DIRPT_4= new TH1D("DIRPT_4","DIRPT_4",binnum,0,ptmax);


  DECPT= new TH1D("DECPT","DECPT",binnum,0,ptmax);

  DECPT_0= new TH1D("DECPT_0","DECPT_0",binnum,0,ptmax);
  DECPT_1= new TH1D("DECPT_1","DECPT_1",binnum,0,ptmax);
  DECPT_2= new TH1D("DECPT_2","DECPT_2",binnum,0,ptmax);
  DECPT_3= new TH1D("DECPT_3","DECPT_3",binnum,0,ptmax);
  DECPT_4= new TH1D("DECPT_4","DECPT_4",binnum,0,ptmax);


  TotalTagged_0= new TH1D("TotalTagged_0","TotalTagged_0",binnum,0,ptmax);
  FalseTag_0= new TH1D("FalseTagged_0","FalseTagged_0",binnum,0,ptmax);
  TrueTag_0 = new TH1D("TrueTagged_0","TrueTagged_0",binnum,0,ptmax);

  TotalTagged_1= new TH1D("TotalTagged_1","TotalTagged_1",binnum,0,ptmax);
  FalseTag_1= new TH1D("FalseTagged_1","FalseTagged_1",binnum,0,ptmax);
  TrueTag_1 = new TH1D("TrueTagged_1","TrueTagged_1",binnum,0,ptmax);

  TotalTagged_2= new TH1D("TotalTagged_2","TotalTagged_2",binnum,0,ptmax);
  FalseTag_2= new TH1D("FalseTagged_2","FalseTagged_2",binnum,0,ptmax);
  TrueTag_2 = new TH1D("TrueTagged_2","TrueTagged_2",binnum,0,ptmax);

  TotalTagged_3= new TH1D("TotalTagged_3","TotalTagged_3",binnum,0,ptmax);
  FalseTag_3= new TH1D("FalseTagged_3","FalseTagged_3",binnum,0,ptmax);
  TrueTag_3 = new TH1D("TrueTagged_3","TrueTagged_3",binnum,0,ptmax);

  TotalTagged_4= new TH1D("TotalTagged_4","TotalTagged_4",binnum,0,ptmax);
  FalseTag_4= new TH1D("FalseTagged_4","FalseTagged_4",binnum,0,ptmax);
  TrueTag_4 = new TH1D("TrueTagged_4","TrueTagged_4",binnum,0,ptmax);

  pi0Tagged= new TH1D("True_pi0_Tagged","True_pi0_Tagged",binnum,0,ptmax);
  Ntag = new TH1D("NTag","NTag",binnum,0,ptmax);
  
  DELTAx = new TH1F("DELTAx","DELTAx",200,-10,10);
  DELTAy = new TH1F("DELTAy","DELTAy",200,-10,10);
  DELTAz = new TH1F("DELTAz","DELTAz",200,-10,10);

  for(int i=0;i<6;i++){
    sprintf(dummy,"hitdiag_sec%d",i);
    hitdiag[i] = new TH2F(dummy,dummy,72,-0.5,71.5,36,-0.5,35.5);
  }  
  for(int i=6;i<8;i++){
    sprintf(dummy,"hitdiag_sec%d",i);
    hitdiag[i] = new TH2F(dummy,dummy,94,-0.5,93.5,46,-0.5,45.5);
  }  


  for(int i=0;i<10;i++)
    {
      sprintf(dummy,"gammaprob%d",i);
      gammaprob[i]=new TH1D(dummy,dummy,binnum,0,ptmax);
      sprintf(dummy,"prob%d",i);
      prob[i]=new TH1D(dummy,dummy,binnum,0,ptmax);
      gammaprob[i]->Sumw2();
      prob[i]->Sumw2();
    }

  realpi0pt->Sumw2();
  realetapt->Sumw2();
  realetaprimept->Sumw2();
  realomegapt->Sumw2();
  gammapt->Sumw2();

  pi0gamma->Sumw2();
  etagamma->Sumw2();
  etaprimegamma->Sumw2();
  omegagamma->Sumw2();

  gamma_acc->Sumw2();
  gamma_acc_nosm->Sumw2();
  pi0gamma_acc->Sumw2();
  pi0gamma_acc_cut->Sumw2();
  openangle->Sumw2();

  beforecut->Sumw2();
  aftercut->Sumw2();

  accbeforecut->Sumw2();
  accaftercut->Sumw2();

  cocktail->Sumw2();

  tightptvspt->Sumw2();
  tightptvspt_wo->Sumw2();
  massvspt->Sumw2();

  loosein->Sumw2();
  looseall->Sumw2();

  beforecut_photon->Sumw2();
  aftercut_photon->Sumw2();

  TotalTagged_0->Sumw2();
  FalseTag_0->Sumw2();
  TrueTag_0->Sumw2();

  TotalTagged_1->Sumw2();
  FalseTag_1->Sumw2();
  TrueTag_1->Sumw2();
  TotalTagged_2->Sumw2();
  FalseTag_2->Sumw2();
  TrueTag_2->Sumw2();
  TotalTagged_3->Sumw2();
  FalseTag_3->Sumw2();
  TrueTag_3->Sumw2();
  TotalTagged_4->Sumw2();
  FalseTag_4->Sumw2();
  TrueTag_4->Sumw2();


  pi0Tagged->Sumw2();
  Ntag->Sumw2();


  return;
}


