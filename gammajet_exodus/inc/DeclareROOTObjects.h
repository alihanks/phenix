class TH2F;
class TH3F;
class TH1D;
class TNtuple;
class TF1;

//INCLUDEFLAG TNtuple *primaries; 
//INCLUDEFLAG TNtuple *pairs;
//INCLUDEFLAG TNtuple *singles;

INCLUDEFLAG TH3F    *PTpi_ZEMCpi_PTgam;
INCLUDEFLAG TH3F    *PTpi_ZEMCpi_PTgam_MISS;
INCLUDEFLAG TH2F    *RECONPI_ZEMC_PT;
INCLUDEFLAG TH3F    *RECONPI_ZEMC_RECONPT_REALPT;
INCLUDEFLAG TH3F    *PTpi_ZEMCpi_PTgam_MISS_ETA;
INCLUDEFLAG TH3F    *PTpi_ZEMCpi_PTgam_MISS_PI0;


INCLUDEFLAG TH2F *RECONPI_ZEMC_PTTRUE_PTRECO_false;
INCLUDEFLAG TH3F *PTpi_ZEMCpi_PTgam_false;
INCLUDEFLAG TH3F *RECONPI_ZEMC_PTTRUE_PTRECO_real;
INCLUDEFLAG TH3F *PTpi_ZEMCpi_PTgam_TRUE;

INCLUDEFLAG TH3F *PTpi_ZEMCpi_PTgam_PI0;
INCLUDEFLAG TH3F *PTpi_ZEMCpi_PTgam_ETA;
INCLUDEFLAG TH3F *PTpi_ZEMCpi_PTgam_DIR;

INCLUDEFLAG TH3F *RECONPI_ZEMC_PTTRUE_PTRECO_allfalse;

INCLUDEFLAG TH2F *PTpi_PTgam_truetag;
INCLUDEFLAG TH2F *PTpi_PTgam_falsetag;
INCLUDEFLAG TH2F *PTpi_PTgam_falsetag_dir;
INCLUDEFLAG TH2F *PTpi_PTgam_falsetag_eta;
INCLUDEFLAG TH2F *PTpi_PTgam_MISS_PI0;
INCLUDEFLAG TH2F *PTpi_PTgam_MISS_ETA;

INCLUDEFLAG TH2F *PTpi_PTgam_PI0;
INCLUDEFLAG TH2F *PTpi_PTgam_ETA;
INCLUDEFLAG TH2F *PTpi_PTgam_DIR;
INCLUDEFLAG TH2F *PTpi_PTgam_TRUE;
INCLUDEFLAG TH2F *PTpi_PTgam_false;
INCLUDEFLAG TH2F *PTpi_PTgam_False_PI0pair;

INCLUDEFLAG TH2F    *ptpivsptgam;
INCLUDEFLAG TH2F    *ptpivsptgam_tag;
INCLUDEFLAG TH2F    *ppivspgam;
INCLUDEFLAG TH2F    *ppivspgam_tag;

INCLUDEFLAG TH2F    *eta_trueptvsrecopt;
INCLUDEFLAG TH2F    *direct_trueptvsrecopt;
INCLUDEFLAG TH2F    *pi0_trueptvsrecopt;
INCLUDEFLAG TH2F    *realpi0_trueptvsrecopt;


INCLUDEFLAG TH2F    *eta_leadptvsrecopt;
INCLUDEFLAG TH2F    *direct_leadptvsrecopt;
INCLUDEFLAG TH2F    *pi0_leadptvsrecopt;
INCLUDEFLAG TH2F    *realpi0_leadptvsrecopt;

INCLUDEFLAG TH1D    *pi_asym[2];
INCLUDEFLAG TH1D    *hpairinvm;
INCLUDEFLAG TH1D    *realpi0pt;
INCLUDEFLAG TH1D    *realetapt;
INCLUDEFLAG TH1D    *realetaprimept;
INCLUDEFLAG TH1D    *realomegapt;
INCLUDEFLAG TH1D    *gammapt;
INCLUDEFLAG TH1D    *EVENTS;
INCLUDEFLAG TH1F    *ZVERTEX;
INCLUDEFLAG TH1F    *DELTAx;
INCLUDEFLAG TH1F    *DELTAy;
INCLUDEFLAG TH1F    *DELTAz;

INCLUDEFLAG TH1D    *gammaprob[10];
INCLUDEFLAG TH1D    *prob[10];

INCLUDEFLAG TH2F    *INVMASS[8];
INCLUDEFLAG TH2F    *INVMASS_PBGL[8];
INCLUDEFLAG TH2F    *pi_INVMASS[3];

INCLUDEFLAG TH1D    *invmass;
INCLUDEFLAG TH2F    *invmassvspt;
INCLUDEFLAG TH1D    *gamma_acc;
INCLUDEFLAG TH1D    *gamma_acc_nosm;
INCLUDEFLAG TH1D    *gamma_acc_phi;
INCLUDEFLAG TH1D    *gamma_acc_ert_phi;
INCLUDEFLAG TH1D    *gamma_acc_phi_live;
INCLUDEFLAG TH1D    *gamma_acc_phi_0;
INCLUDEFLAG TH1D    *gamma_acc_phi_1;
INCLUDEFLAG TH1D    *gamma_acc_phi_m1;
INCLUDEFLAG TH1D    *gamma_acc_phi_2;
INCLUDEFLAG TH1D    *pi0gamma_acc;
INCLUDEFLAG TH1D    *pi0gamma_acc_cut;
INCLUDEFLAG TH1D    *pi0gamma;
INCLUDEFLAG TH1D    *etagamma;
INCLUDEFLAG TH1D    *etaprimegamma;
INCLUDEFLAG TH1D    *omegagamma;
INCLUDEFLAG TH1D    *gamma_ecore;
INCLUDEFLAG TH1D    *gamma_pair[8];
INCLUDEFLAG TH1D    *pi0gammapt[3];

INCLUDEFLAG TH1D    *beforecut;
INCLUDEFLAG TH1D    *aftercut;

INCLUDEFLAG TH1D    *accbeforecut;
INCLUDEFLAG TH1D    *accaftercut;
INCLUDEFLAG TH1D    *cocktail;
INCLUDEFLAG TH1D    *accvar;

INCLUDEFLAG TH2F    *hitdiag[8];
INCLUDEFLAG TH2F    *tightptvspt;
INCLUDEFLAG TH2F    *tightptvspt_wo;
INCLUDEFLAG TH2F    *massvspt;

INCLUDEFLAG TH1D    *openangle;

INCLUDEFLAG TH1F    *loosein;
INCLUDEFLAG TH1F    *looseall;

INCLUDEFLAG TH1D    *beforecut_photon;
INCLUDEFLAG TH1D    *aftercut_photon;

INCLUDEFLAG TH1D    *TotalTagged_0;
INCLUDEFLAG TH1D    *FalseTag_0;
INCLUDEFLAG TH1D    *TrueTag_0;

INCLUDEFLAG TH1D    *TotalTagged_1;
INCLUDEFLAG TH1D    *FalseTag_1;
INCLUDEFLAG TH1D    *TrueTag_1;

INCLUDEFLAG TH1D    *TotalTagged_2;
INCLUDEFLAG TH1D    *FalseTag_2;
INCLUDEFLAG TH1D    *TrueTag_2;

INCLUDEFLAG TH1D    *TotalTagged_3;
INCLUDEFLAG TH1D    *FalseTag_3;
INCLUDEFLAG TH1D    *TrueTag_3;

INCLUDEFLAG TH1D    *TotalTagged_4;
INCLUDEFLAG TH1D    *FalseTag_4;
INCLUDEFLAG TH1D    *TrueTag_4;

INCLUDEFLAG TH1D    *TRIGPT;
INCLUDEFLAG TH1D    *pi_TRIGPT;

INCLUDEFLAG TH1D    *TRIGPT_0;
INCLUDEFLAG TH1D    *TRIGPT_1;
INCLUDEFLAG TH1D    *TRIGPT_2;
INCLUDEFLAG TH1D    *TRIGPT_3;
INCLUDEFLAG TH1D    *TRIGPT_4;

INCLUDEFLAG TH1D    *DIRPT;
INCLUDEFLAG TH1D    *DIRPT_0;
INCLUDEFLAG TH1D    *DIRPT_1;
INCLUDEFLAG TH1D    *DIRPT_2;
INCLUDEFLAG TH1D    *DIRPT_3;
INCLUDEFLAG TH1D    *DIRPT_4;

INCLUDEFLAG TH1D    *DECPT;
INCLUDEFLAG TH1D    *DECPT_0;
INCLUDEFLAG TH1D    *DECPT_1;
INCLUDEFLAG TH1D    *DECPT_2;
INCLUDEFLAG TH1D    *DECPT_3;
INCLUDEFLAG TH1D    *DECPT_4;


INCLUDEFLAG TH1D    *pi0Tagged;
INCLUDEFLAG TH1D    *Ntag;
		      
