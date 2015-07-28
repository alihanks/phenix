//-----------------------------------------------------------------------------
//
//  Generate ROOT output file
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>

#include <TRandom3.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TNtuple.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticlePropertyList.h"
#include "THmulf.h"
#include "TLorentzVector.h"

#define PI  3.141592653589793238

#define  INCLUDEFLAG extern
#include "DeclareROOTObjects.h"

bool   ParticleIsLepton(Particle*);
double thetaof(Mom3);
double phiof(Mom3);
Mom4   ApplyResolution(int, Mom4);
Mom4   ApplyEnergyResolution(Mom4, int);
Mom4   ApplyPositionResolution(int, Mom4, float, char*);
int    PHENIXFilter(Particle*, ParticlePropertyList*);
bool   PHENIXPhotonFilter(Particle* kp, int fid = 0);
bool   PairInAcceptance(int, Mom4, Mom4);
float phenix_geom(int,float,float,float,float,float,char*, float & smeared_x, float & smeared_y, float & smeared_z);

void FillROOTObjects(int setup, double dNdy_pi0, double N_coll,
		     char *output_file, 
		     ParticleList *PList, 
		     ParticlePropertyList *PPList)
{

  //  gRandom = new TRandom3();
  //  gRandom->SetSeed(0);
  
  //run10
  //char deadfile[]="/phenix/hhj/jfrantz/gammajet_exodus/livetowers_r10_v4.txt";

  //run 8 dau
  char deadfile[]="/phenix/hhj/ahanks/gammajet_exodus/ttemp_livetowers_r8_test.dat";

  double weight, fillweight, fillweight1, fillweight2;
  double n_pi0_generated;
  double sigma_inelastic, abs_norm;

  double ptcut = 0.2;
  double type,type1,type2,per=16;
  double sepacut=12.0/500.0;
  double opangle;

  // storage for pi0 pair
  Mom4 pairmom[2];
  double pairweight[2];
  int pairaccept[2];
  int pi0pair=-1;
  double pairtype[2];
  double minv;

  //for decay photons
  Mom4   mom1, mom2;
  int accept, accept1, accept2;
  double pxe, pye, pze, pxp, pyp, pzp;
  double pte, ptp, pe, pp;
  double pxe_nosm, pye_nosm, pte_nosm;
  double phie, phip;
  double vtxz=0;
  double theta;
  double thetae, thetap, etae, etap;
  double ye,yp;
  double px1, py1, pz1, px2, py2, pz2;
  double pt1, pt2, p1, p2;
  double phi1, phi2;
  double theta1, theta2, eta1, eta2;
  double ecore1, ecore2;
  double parent1, parent2;
  double pipt1, pipt2;
  int iparent1, iparent2;

  //for prime particle
  int 	 pid = 0;
  Mom4 primmom1;
  Mom4 testmom1;
  double primweight, primfillweight, primpt, primp;

  // z-position in the emcal plane
  double primzemc;

  n_pi0_generated = 10.0e+6;
  sigma_inelastic = 42.2;
  abs_norm        = dNdy_pi0; // for AuAu: multiplicity
  // abs_norm        = dNdy_pi0 * sigma_inelastic; // for pp: cross section

  //cout << "multiplicity: " << abs_norm <<endl;

  PLNode   * CurrentNode     = PList->GetHeadNode();
  PLNode   * NextNode        = 0;
  Particle * CurrentParticle = 0;
  Particle * NextParticle    = 0;
  Particle * CurrentPrimary  = 0;
  
  // here we start with the diphoton analysis
  static int  neventK = 0;


  //variables I added for tagging method
  int nparts=10000;
  // static int nparts=PList->GetLength();

  float gamenergy[nparts];
  float gammom[nparts];
  float gamaccept[nparts];
  float gamphi[nparts]; 
  float gamparent[nparts]; 
  int   gamindex[nparts]; 
  float gampi0pt[nparts]; 
  float gamtheta[nparts];
  float gampt[nparts];
  float gampx[nparts];
  float gampy[nparts];
  float gampz[nparts];
  float gameta[nparts];
  float gamzemc[nparts];

  //int eventcounter=0;

  for(int itest=1; itest<nparts; itest++)
    {
      gamenergy[itest]=0.0;
      gamzemc[itest]=0.0;
      gammom[itest]=0.0;
      gamaccept[itest]=0.0;
      gamphi[itest]=0.0; 
      gamparent[itest]=0.0; 
      gamindex[itest]=0; 
      gampi0pt[itest]=0.0;  
      gamtheta[itest]=0.0;
      gampt[itest]=0.0;
      gampx[itest]=0.0;
      gampy[itest]=0.0;
      gampz[itest]=0.0;
      gameta[itest]=0.0;
    }
  int parentid=0;

  //cout <<"HERE we are in fill.cpp" <<endl;
  //cout << "nparts: " << PList->GetLength() <<endl;
  for ( int ipart=1; ipart<=PList->GetLength(); ipart++) 
    {
      //cout<<" ipart "<<ipart<< " of "<< PList->GetLength() << endl;
    CurrentNode = CurrentNode->GetNextNode();
    CurrentParticle = CurrentNode->Get(0);

    //cout << CurrentParticle->GetID() << ":" << CurrentParticle->GetGeneration() << endl;
    
    //Particle of generation=1 is always immediately
    //followed by it's decay particles so this format works    
    //cout << "for ipart " << ipart << " parentid=" << parentid <<endl;
    if ( CurrentParticle->GetGeneration()==1 )
      {
	CurrentPrimary = CurrentParticle;
	pid            = CurrentParticle->GetID();
	parentid = ipart;
	//cout <<"particle number:" << ipart << " out of " << nparts <<endl;

	primweight     = CurrentParticle->GetWeight();
	primfillweight = abs_norm*primweight/(0.5*n_pi0_generated);
	primfillweight = primfillweight/0.025;

	// pi0 4-vector
	primmom1 = CurrentParticle->Get4mom();
	primpt = sqrt(primmom1.Getp().Getpx() *
		      primmom1.Getp().Getpx() +
		      primmom1.Getp().Getpy() *
		      primmom1.Getp().Getpy() );

	primp  = sqrt(primmom1.Getp().Getpx() *
		      primmom1.Getp().Getpx() +
		      primmom1.Getp().Getpy() *
		      primmom1.Getp().Getpy() +
		      primmom1.Getp().Getpz() *
		      primmom1.Getp().Getpz() );

	vtxz = CurrentParticle->GetzVertex();
	ZVERTEX->Fill(vtxz);
	//double primtheta = acos(primmom1.Getp().Getpz()/primp);
	//cout<<" vtxz pi0 "<<vtxz<<endl;
	//use pbsc as ref radius
	primzemc  = 510.0*primmom1.Getp().Getpz()/primpt+vtxz;
	//primzemc = 510.0*cos(primtheta)+vtxz;


	//Do not weight the input particles
	if (pid == 111)
	  realpi0pt->Fill(primpt);
      
	if (pid == 221)
	  realetapt->Fill(primpt);
	
	if (pid == 331)
	  realetaprimept->Fill(primpt);
	
	//if (pid == 223 ) //Now filling omega histo with photon/direct photons
	if(pid==22 || pid==-111)
	  realomegapt->Fill(primpt);
	
      }//close if generation==1
    
    
    if ( CurrentParticle->GetID()==22 	 || CurrentParticle->GetID()== -111  ){ // -111 = direct photon

      //if ( CurrentParticle->GetID()== 111  ) 
      
	
      
      // correct weighting for pairs (if decay chains are full implemented!) 
      weight   = CurrentParticle->GetWeight();
      fillweight = abs_norm*weight/(0.5*n_pi0_generated);
      fillweight = fillweight/0.025;
      // photon 4-vector
      mom1 = CurrentParticle->Get4mom();
      //vtxz = CurrentParticle->GetzVertex();
      //cout<<" vtxz gamma "<<vtxz<<endl;

      // resolution
      
      pxe_nosm = mom1.Getp().Getpx();
      pye_nosm = mom1.Getp().Getpy();
      pte_nosm = sqrt(pxe_nosm*pxe_nosm+pye_nosm*pye_nosm);
      
            
      float phitemp = atan2(mom1.Getp().Getpy(),mom1.Getp().Getpx());
      //fixed bug below that caused phitemp to always be pi/4 or -3pi/4
      //float phitemp = atan2(mom1.Getp().Getpy(),mom1.Getp().Getpy());
      //phitemp is only used for the energy resolution
      //cout<< "phitemp" << phitemp <<endl;
      if(phitemp<-PI/2.) phitemp+=2*PI;
      if(phitemp>3*PI/2.) phitemp-=2*PI;
      
      if(phitemp<2.95)mom1 = ApplyEnergyResolution(mom1,0);
      else mom1 = ApplyEnergyResolution(mom1,1);
      
      
      pxe = mom1.Getp().Getpx();
      pye = mom1.Getp().Getpy();
      pze = mom1.Getp().Getpz();
      
      float smeared_x=pxe;
      float smeared_y=pye;
      float smeared_z=pze;
      
      //new position resolution smearing in phenix_geom
      
      int dopossmear=1;
      
      accept = -9;
      for(int jj=0;jj<8;jj++){
	//testmom1=ApplyPositionResolution(jj,mom1,vtxz,deadfile);
	
	int smearedjj=jj;
	if(dopossmear==1) smearedjj=jj+10;
	
	int isaccept = phenix_geom(smearedjj,(float)pxe,(float)pye,(float)pze,(float)vtxz,(float)mom1.GetE(),deadfile,smeared_x,smeared_y,smeared_z);
	if(isaccept>0){
	  accept=isaccept;
	  //cout<< "before_smearing " << pxe << " " << pye << " " << pze <<endl;
	  if(dopossmear==1)
	    {
	      pxe=smeared_x;
	      pye=smeared_y;
	      pze=smeared_z;

	      
	      mom1.Setp(pxe,pye,pze);
	      mom1.SetE(sqrt(pxe*pxe+pye*pye+pze*pze));
	      
	    }
	  	  
	}
      }

      pe=mom1.GetE();
      
      pxe = mom1.Getp().Getpx();
      pye = mom1.Getp().Getpy();
      pze = mom1.Getp().Getpz();
      pte = sqrt(pxe*pxe+pye*pye);
      thetae = tan(pte/pze);
      etae = -log(tan(thetae/2.));
      cocktail->Fill(pte,1.0/pte);
 
	
      //complicated.  -9 means did not hit emcal.  all other negative means it hit emcal but near at a bad ert tile (p+p only) -man
	
      
      float eff=0.95;
      //if(pte<20.0) eff=0.0075*pte+0.8; //for centrality 0-10
      //Maybe 80% is too low?
      
      if(gRandom->Rndm()>eff) accept=-9;
      
      //if(accept>0) cout << "eff_passed" <<endl;
      
      accvar->Fill(accept);
       
     if(pte>ptcut){
	//gammapt->Fill(pte, fillweight);
	
	//gammapt->Fill(pte); //Not using fill weight
	
	if(accept!=-9){
	  
	  //Not using fill weights!
	  //gamma_acc->Fill(pte,fillweight);    DAT[i]->Project3D("yx")->Draw("zcol");
	  
	  //gamma_acc_nosm->Fill(pte_nosm,fillweight);
	  
	  //gamma_acc->Fill(pte);
	  gamma_acc_nosm->Fill(pte_nosm);
	  
	  float phi=atan2(pye,pxe);
	  if(phi<-PI/2.)phi+=2*PI;
	  if(phi>3*PI/2.)phi-=2*PI;
	  //if(accept>0)gamma_acc_phi->Fill(phi);
	  //if(accept==0)gamma_acc_phi->Fill(phi);
	  
	  
	}
      }
      
      //cout <<"HERE after filling gamma_acc histos" <<endl;
      
      if(pid==111||pid==221)
	{
	  ////cout<<" pid "<<pid<<endl;
	  if(pte>ptcut){
	    //Not using fill weights
	    //pi0gamma->Fill(pte,fillweight);
	    pi0gamma->Fill(pte);
	    
	  }
	  //store pi0 pair information for pair analysis
	  pi0pair++;
	  pairmom[pi0pair]=mom1;
	  pairaccept[pi0pair]=accept;
	  pairweight[pi0pair]=fillweight;
	  pairtype[pi0pair]=type;
	}
      
      
      
      float phi=atan2(pye,pxe);
      gamenergy[ipart]=pe;
      //gammom[ipart]=mom1;
      gamaccept[ipart]=accept;
      gamphi[ipart]=phi; 
      gamindex[ipart]=parentid;
      gamzemc[ipart]=primzemc;  
      gamparent[ipart]=pid;
      gampi0pt[ipart]=primpt; 
      gamtheta[ipart]=thetae;
      gampt[ipart]=pte;
      gameta[ipart]=etae;
      
      gampx[ipart]=pxe;
      gampy[ipart]=pye;
      gampz[ipart]=pze;
      
      gammapt->Fill(pte);
      
      if(pid==331)
	etaprimegamma->Fill(pte);
      if(pid==221)
	etagamma->Fill(pte);
      if(pid==223)
	omegagamma->Fill(pte);
      
    }//if particle is a gamma
    
    if(pi0pair==1)
      {
	int intagwindow=0;
	int inpiwindow=0;
	//cout<<" pi0pair!!!!! "<<endl;
	pxe = pairmom[0].Getp().Getpx();
	pye = pairmom[0].Getp().Getpy();
	pze = pairmom[0].Getp().Getpz();
	pte = sqrt(pxe*pxe+pye*pye);

	pe = sqrt(pxe*pxe+pye*pye+pze*pze);

	phie = atan2(pye,pxe);
	if(phie<PI/2.) phie+=2*PI;
	if(phie>3*PI/2.) phie-=2*PI;
       
	pxp = pairmom[1].Getp().Getpx();
	pyp = pairmom[1].Getp().Getpy();
	pzp = pairmom[1].Getp().Getpz();
	ptp = sqrt(pxp*pxp+pyp*pyp);

	phip =  atan2(pyp,pxp);
	if(phip<PI/2.) phip+=2*PI;
	if(phip>3*PI/2.) phip-=2*PI;

	thetap = tan(ptp/pzp);
	etap = -log(tan(thetap/2.));
	pp = sqrt(pxp*pxp+pyp*pyp+pzp*pzp);
	
	fillweight1=pairweight[0];
	fillweight2=pairweight[1];
	  
	accept1=pairaccept[0];
	accept2=pairaccept[1];

	type1=pairtype[0];
	type2=pairtype[1];

	float ecoree=sqrt(pxe*pxe+pye*pye+pze*pze);
	float ecorep=sqrt(pxp*pxp+pyp*pyp+pzp*pzp);

	float mass= sqrt(2.0*ecoree*ecorep - 2.0*(pxe*pxp+pye*pyp+pze*pzp)); 


	if(mass>0.12 && mass<0.16) { inpiwindow=1;}
	if(mass>0.118 && mass<0.162) { intagwindow=1;}

	TLorentzVector pho1, pho2, pi0;
	pho1.SetPxPyPzE(pxe,pye,pze,ecoree);
	pho2.SetPxPyPzE(pxp,pyp,pzp,ecorep);
	pi0=pho1+pho2;

	float pi0pt=pi0.Pt();		 
      
	if(accept1>0&&accept2>0) pi0gammapt[2]->Fill(pi0pt);
	else if(accept1>0 || accept2>0) pi0gammapt[1]->Fill(pi0pt);
	else pi0gammapt[0]->Fill(pi0pt);

	if(accept1>-9&&pte>5.0)gamma_acc_phi->Fill(phie);
	if(accept2>-9&&ptp>5.0)gamma_acc_phi->Fill(phip);

	if(accept1>0&&pte>5.0)gamma_acc_phi_live->Fill(phie);
	if(accept2>0&&ptp>5.0)gamma_acc_phi_live->Fill(phip);


	int firedert=0;

	//make fiducial cuts  -- 155 cm for higher pt gamma
	double zemcp= 510.0*pzp/ptp + vtxz;
	double zemce= 510.0*pze/pte + vtxz;
      
      
	//ignore ert for now:
	firedert=1;
	//Is photon 1 included in my inclusive sample?
	if(firedert ==1 && accept1> 0 && fabs(zemce) < 155 && ecoree>1.0){
	  
	  //Not using ert weight from accept value
	  ptpivsptgam->Fill(pte,primpt);
	  PTpi_ZEMCpi_PTgam->Fill(primpt,primzemc,pte);
	  ppivspgam->Fill(pe,primp);

	  if(accept1>-9&&pte>5.0)gamma_acc_ert_phi->Fill(phie);

	  //Do I not tag photon 1 because I miss it's partner?
	  if(accept2==-9||ecorep<1.0) PTpi_ZEMCpi_PTgam_MISS->Fill(primpt,primzemc,pte);
	}
	//Repeat for photon 2
	if(firedert ==1 && accept2> 0 && fabs(zemcp) < 155 && ecorep>1.0){
	  
	  ptpivsptgam->Fill(ptp,primpt);
	  PTpi_ZEMCpi_PTgam->Fill(primpt,primzemc,ptp);	  
	  ppivspgam->Fill(pp,primp);
	  
	  if(accept2>-9&&ptp>5.0)gamma_acc_ert_phi->Fill(phip);

	  //don't use ert weighting and changed pt cut
	  //Do I not tag photon 2 because I miss it's partner?
	  if(accept1==-9||ecoree<1.0) PTpi_ZEMCpi_PTgam_MISS->Fill(primpt,primzemc,ptp);
	}
	//For miss sharkfin do I measure both photons but reconstruct the invmass outside our window?
	if(accept1>0&&accept2>0&&ecoree>1.0&&ecorep>1.0 && !intagwindow && (fabs(zemcp) < 155 || fabs(zemce) < 155)) PTpi_ZEMCpi_PTgam_MISS->Fill(primpt,primzemc,pte);	

	
	if(accept1>-9&&(accept2<-8||ecorep<1.0)) {
	  ppivspgam_tag->Fill(pe,primp);
	  ppivspgam_tag->Fill(pe,primp);
	}
	if(accept2>-9&&(accept1<-8||ecoree<1.0)) {
	  ppivspgam_tag->Fill(pp,primp);
	  ppivspgam_tag->Fill(pp,primp);
	}
	
	//invmass->Fill(mass);
	
	
	//Apply cuts that match data
	float sume12 = ecoree+ecorep;
	//double zemc1= 510.0*pz1/pt1 + vtxz;
	//double zemc2= 510.0*pz2/pt2 + vtxz;
	
	float asym12 =  (ecoree - ecorep)/ sume12; 
	
	const float mina = 0.15;  //this is the smallest asym value used for the cut
	
	int passasym=1;
	if (sume12 < 5.25 && fabs(asym12) > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25))) { passasym=0; }
	
	
	//Apply pi0 cuts to fill the RECON & invmas histos
	
	if(accept1>0&&accept2>0&&ecoree>1.0&&ecorep>1.0&&sume12>=4.0&& passasym==1){ 
	  invmassvspt->Fill(pi0pt,mass);
	  if(inpiwindow){
	    RECONPI_ZEMC_PT->Fill(primzemc,primpt);
	  }
	  
	}
      
	//openning angle of two pion-photons
	opangle = (pairmom[0].Getp()*pairmom[1].Getp()) /
	  (sqrt(pairmom[0].Getp()*pairmom[0].Getp())*
	   sqrt(pairmom[1].Getp()*pairmom[1].Getp()));
	opangle=acos(opangle);
	
	//Not using fill weights
	//openangle->Fill(opangle,primfillweight);
	openangle->Fill(opangle);
	
	//pi0pair process
	pi0pair=-1;
	
      }
    //cout <<"HERE out of the pi0pair loop" <<endl;
    
    }

  //cout <<"HERE out of the particle loop" <<endl;
  
  
  //off for single particle mode
  int tagging=1;

  if(tagging==1)
    {
      //cout <<"HERE in tagging loop with " << nparts <<endl;
      for(int gamtrig=1; gamtrig<nparts; gamtrig++)
	{
	  int taggedtrig=0;
	  int falsepair=1;
	  
	  double zemc1=0.0;
	  double zemc2=0.0;

	  if(gamaccept[gamtrig]>0)
	    {
	      //zemc1=gamzemc[gamtrig];
	      accept1=gamaccept[gamtrig];
	      ecore1=gamenergy[gamtrig];
	      pt1=gampt[gamtrig];
	      px1=gampx[gamtrig];
	      py1=gampy[gamtrig];
	      pz1=gampz[gamtrig];
	      phi1=gamphi[gamtrig];
	      eta1=gameta[gamtrig];
	      theta1=gamtheta[gamtrig];
	      parent1=gamparent[gamtrig];
	      iparent1=gamindex[gamtrig];
	      pipt1=gampi0pt[gamtrig];
	      
	      //if(parent1==0 && ecore1<1) continue;
	      if(ecore1<1) continue;
	      if(accept1<= 0) continue; //shouldn't do anything

	      zemc1= 510.0*pz1/pt1 + vtxz;
	      if(fabs(zemc1) > 155) continue;
	      
	      int ispi0=0;
	      
	      if(pt1>5 && pt1<20){
		TRIGPT->Fill(pt1);
		if(parent1==22 || parent1==-111) DIRPT->Fill(pt1);
		else DECPT->Fill(pt1);
		
		for(int gampart=1; gampart<nparts; gampart++)
		  {
		    if(gampart==gamtrig) continue;		      
		    accept2=gamaccept[gampart];
		    ecore2=gamenergy[gampart];
		    pt2=gampt[gampart];
		    px2=gampx[gampart];
		    py2=gampy[gampart];
		    pz2=gampz[gampart];
		    phi2=gamphi[gampart];
		    eta2=gameta[gampart];
		    theta2=gamtheta[gampart];
		    parent2=gamparent[gampart];
		    iparent2=gamindex[gampart];
		    pipt2=gampi0pt[gampart];
		    

		    //This is only done for pi0 loop
		    //if(ecore2<1.0||ecore1<ecore2|| ecore1 == ecore2) continue;
		    if(ecore2<1.0) continue;
		    
		    float mass= sqrt(2.0*ecore1*ecore2 - 2.0*(px1*px2+py1*py2+pz1*pz2));
		    
		    //invmass->Fill(mass);
		    //cout<< "at mass cut with "<< mass <<endl;
		    
		    float dphi = fabs(phi1-phi2);
		    
		    //if(dphi>PI/2) continue;
		    if(accept2<=0)continue;
		    
		  
		    //For pi0 loop not tagging in data
		    //float sume12 = ecore1+ecore2;
		    //if (sume12 < 4.0) continue;
		    
		    //float asym12 =  (ecore1 - ecore2)/ sume12; // fabs not neccesary because upper > check
		    //const float mina = 0.15;  //this is the smallest asym value used for the cut		      
		    //if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25)))continue;		 
		  
		      
		    invmass->Fill(mass);
		    
		    TLorentzVector pho1, pho2, pi0;
		    pho1.SetPxPyPzE(px1,py1,pz1,ecore1);
		    pho2.SetPxPyPzE(px2,py2,pz2,ecore2);
		    pi0=pho1+pho2;
		    
		    float pi0pt=pi0.Pt();		 
		    
		    
		    //fill for all looped pairs don't stop after tagging
		    if(ecore2>1.0) INVMASS[0]->Fill(pt1,mass);
		    if(ecore2>1.5) INVMASS[1]->Fill(pt1,mass);
		    if(ecore2>2.0) INVMASS[2]->Fill(pt1,mass);
		    if(ecore2>2.5) INVMASS[3]->Fill(pt1,mass);
		    if(ecore2>3.0) INVMASS[4]->Fill(pt1,mass);
		    
		    /*
		    if(ispi0<1 && ecore2>1.0) INVMASS[0]->Fill(pt1,mass);
		    if(ispi0<2 && ecore2>1.5) INVMASS[1]->Fill(pt1,mass);
		    if(ispi0<3 && ecore2>2.0) INVMASS[2]->Fill(pt1,mass);
		    if(ispi0<4 && ecore2>2.5) INVMASS[3]->Fill(pt1,mass);
		    if(ispi0<5 && ecore2>3.0) INVMASS[4]->Fill(pt1,mass);
		    */		    
		  
		    /*
		    //switched pt1 to pi0pt		  
		    if(phi1<2.95){		  
		    if(ecore2>.5) INVMASS[0]->Fill(pi0pt,mass);
		    if(ecore2>.75) INVMASS[1]->Fill(pi0pt,mass);
		    if(ecore2>1) INVMASS[2]->Fill(pi0pt,mass);
		    if(ecore2>1.25) INVMASS[3]->Fill(pi0pt,mass);
		    if(ecore2>1.5) INVMASS[4]->Fill(pi0pt,mass);
		    if(ecore2>1.75) INVMASS[5]->Fill(pi0pt,mass);
		    if(ecore2>2) INVMASS[6]->Fill(pi0pt,mass);
		    }else{
		    if(ecore2>.5) INVMASS_PBGL[0]->Fill(pi0pt,mass);
		    if(ecore2>.75) INVMASS_PBGL[1]->Fill(pi0pt,mass);
		    if(ecore2>1) INVMASS_PBGL[2]->Fill(pi0pt,mass);
		    if(ecore2>1.25) INVMASS_PBGL[3]->Fill(pi0pt,mass);
		    if(ecore2>1.5) INVMASS_PBGL[4]->Fill(pi0pt,mass);
		    if(ecore2>1.75) INVMASS_PBGL[5]->Fill(pi0pt,mass);
		    if(ecore2>2) INVMASS_PBGL[6]->Fill(pi0pt,mass);
		    INVMASS_PBGL[7]->Fill(pi0pt,mass);
		    }
		    
		    INVMASS[7]->Fill(pi0pt,mass);
		    */
		    
		    if(mass>0.118 && mass<0.162){ 
		      //TAGGED!
		      //cout << "TAGGED!!!" <<endl;			
		      if(iparent1==iparent2){
			falsepair=0;
			PTpi_PTgam_truetag->Fill(pipt1,pt1);
		      }

		      if(ispi0<1 && ecore2>1.0) ispi0=1; 
		      if(ispi0<2 && ecore2>1.5) ispi0=2;
		      if(ispi0<3 && ecore2>2.0) ispi0=3;
		      if(ispi0<4 && ecore2>2.5) ispi0=4;
		      if(ispi0<5 && ecore2>3.0) ispi0=5;
		      
		    }
		  }

		if(ispi0<5) TRIGPT_4->Fill(pt1);
		if(ispi0<4) TRIGPT_3->Fill(pt1);
		if(ispi0<3) TRIGPT_2->Fill(pt1);
		if(ispi0<2) TRIGPT_1->Fill(pt1);
		if(ispi0<1) TRIGPT_0->Fill(pt1);
		
		if(parent1==22 || parent1==-111){
		  if(ispi0<5) DIRPT_4->Fill(pt1);
		  if(ispi0<4) DIRPT_3->Fill(pt1);
		  if(ispi0<3) DIRPT_2->Fill(pt1);
		  if(ispi0<2) DIRPT_1->Fill(pt1);
		  if(ispi0<1) DIRPT_0->Fill(pt1);
		}else{
		  if(ispi0<5) DECPT_4->Fill(pt1);
		  if(ispi0<4) DECPT_3->Fill(pt1);
		  if(ispi0<3) DECPT_2->Fill(pt1);
		  if(ispi0<2) DECPT_1->Fill(pt1);
		  if(ispi0<1) DECPT_0->Fill(pt1);
		}
		
		if(ispi0>0)
		  {
		    //replaced pi0pt with pt1

		    TotalTagged_0->Fill(pt1);
		    if(parent1==111 || parent1==221) TrueTag_0->Fill(pt1);
		    if(parent1==22 || parent1==-111) FalseTag_0->Fill(pt1);
		    
		    
		    if(ispi0>1){
		      TotalTagged_1->Fill(pt1);
		      if(parent1==111|| parent1==221) TrueTag_1->Fill(pt1);
		      if(parent1==22 || parent1==-111) FalseTag_1->Fill(pt1);
		    }
		    if(ispi0>2){
		      TotalTagged_2->Fill(pt1);
		      if(parent1==111|| parent1==221) TrueTag_2->Fill(pt1);
		      if(parent1==22 || parent1==-111) FalseTag_2->Fill(pt1);
		      }
		    if(ispi0>3){
		      TotalTagged_3->Fill(pt1);
		      if(parent1==111|| parent1==221) TrueTag_3->Fill(pt1);
		      if(parent1==22 || parent1==-111) FalseTag_3->Fill(pt1);
		    }
		    if(ispi0>4){
		      TotalTagged_4->Fill(pt1);
		      if(parent1==111|| parent1==221) TrueTag_4->Fill(pt1);
		      if(parent1==22 || parent1==-111) FalseTag_4->Fill(pt1);
		    }
		    
		    //adding this 2/25
		    if(falsepair && parent1==111) PTpi_PTgam_falsetag->Fill(pipt1,pt1);
		    if(falsepair && parent1==221) PTpi_PTgam_falsetag_eta->Fill(pipt1,pt1);
		    if(falsepair && (parent1==22 || parent1==-111)) PTpi_PTgam_falsetag_dir->Fill(pipt1,pt1);
		    /*
		    if(!falsepair && parent1==111) PTpi_PTgam_truetag->Fill(pipt1,pt1);
		    if(!falsepair && parent1==221) PTpi_PTgam_truetag_eta->Fill(pipt1,pt1);
		    */
		    
		  }else{
		    //Not tagged  
		    //removed ZEMC dependence 2/25
		    if(parent1==111) PTpi_PTgam_MISS_PI0->Fill(pipt1,pt1);
		    if(parent1==221) PTpi_PTgam_MISS_ETA->Fill(pipt1,pt1);

		    
		  }
		
		
	      }//pt bwtn 5 & 20
	    }//if accept1>0
	  //cout << "end of gamtrig loop" <<endl;
	}//gamtrig loop
      //cout <<"end tagging loop" <<endl;
    }//tagging


  //Now do pi0 loop!
  
  if(tagging==1)//pi0 loop
    {
      //cout <<"HERE in pi0 tagging loop with " << nparts <<endl;
      for(int gamtrig=1; gamtrig<nparts; gamtrig++)
	{
	  int taggedtrig=0;
	  int falsepair=0;
	  
	  //First check if photon1 is in the acceptance
	  if(gamaccept[gamtrig]>0)
	    {
	      accept1=gamaccept[gamtrig];
	      ecore1=gamenergy[gamtrig];
	      pt1=gampt[gamtrig];
	      px1=gampx[gamtrig];
	      py1=gampy[gamtrig];
	      pz1=gampz[gamtrig];
	      phi1=gamphi[gamtrig];
	      eta1=gameta[gamtrig];
	      theta1=gamtheta[gamtrig];
	      parent1=gamparent[gamtrig];
	      iparent1=gamindex[gamtrig];
	      pipt1=gampi0pt[gamtrig];
	      	      
	      gamma_acc->Fill(pt1);
	      
	      //if(parent1==0 && ecore1<1) continue;
	      if(ecore1<1) continue;
	      //cout << gamtrig << "Trig is greater than 1 GeV" <<endl;
	      
	      gamma_ecore->Fill(pt1);
	      
	      int ispi0=0;
	      double zemc1= 510.0*pz1/pt1 + vtxz;
	      
	      for(int gampart=1; gampart<nparts; gampart++)
		{
		  if(gampart==gamtrig) continue;
		    
		  //if(gamparent[gampart]<1) continue;
		  //cout << gampart << " not self nor gamparent < 1" <<endl;
		  
		  accept2=gamaccept[gampart];
		  ecore2=gamenergy[gampart];
		  pt2=gampt[gampart];
		  px2=gampx[gampart];
		  py2=gampy[gampart];
		  pz2=gampz[gampart];
		  phi2=gamphi[gampart];
		  eta2=gameta[gampart];
		  theta2=gamtheta[gampart];
		  parent2=gamparent[gampart];
		  iparent2=gamindex[gampart];
		  pipt2=gampi0pt[gampart];
		  
		  
		  TLorentzVector pho1, pho2, pi0;
		  pho1.SetPxPyPzE(px1,py1,pz1,ecore1);
		  pho2.SetPxPyPzE(px2,py2,pz2,ecore2);
		  pi0=pho1+pho2;
		  
		  float pi0pt=pi0.Pt();		 
		  
		  if(accept2!=0) gamma_pair[0]->Fill(pi0pt);
		  
		  if(accept2<=0) gamma_pair[1]->Fill(pi0pt);
		  if(accept2<=0)continue;
		  
		  if(ecore2<1.0) gamma_pair[2]->Fill(pi0pt);
		  if(ecore1<ecore2|| ecore1 == ecore2) gamma_pair[3]->Fill(pi0pt);
		  if(ecore2<1.0||ecore1<ecore2|| ecore1 == ecore2) continue;   
		  
	  
		  float mass= sqrt(2.0*ecore1*ecore2 - 2.0*(px1*px2+py1*py2+pz1*pz2)); 
		
		  
		  //make fiducial cuts  -- 155 cm for higher pt gamma
		  double zemc2= 510.0*pz2/pt2 + vtxz;
		  
		  float dphi = fabs(phi1-phi2);
		  
		  float sume12 = ecore1+ecore2;
		  if (sume12 < 4.0) gamma_pair[4]->Fill(pi0pt);;
		  
		  if (sume12 < 4.0) continue;
		
		  
		  
		  if(fabs(zemc1) > 155){
		    gamma_pair[5]->Fill(pi0pt); 
		    //continue;
		  }

		  
		  float asym12 =  (ecore1 - ecore2)/ sume12; // fabs not neccesary because upper > check
		  const float mina = 0.15;  //this is the smallest asym value used for the cut
		  
		  pi_asym[0]->Fill(asym12);
		  
		  if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25))) gamma_pair[6]->Fill(pi0pt);
		  if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25)))continue;		 
		
		  
		  gamma_pair[7]->Fill(pi0pt);
		  
		  invmass->Fill(mass);
	     		  
		  //switched pt1 to pi0pt		  
		  if(phi1<2.95){		  
		    INVMASS_PBGL[0]->Fill(pi0pt,mass);
		  }else{
		    INVMASS_PBGL[1]->Fill(pi0pt,mass);
		  }
		  if (iparent1!=iparent2) pi_INVMASS[0]->Fill(pi0pt,mass);
		  if (iparent1==iparent2) pi_INVMASS[1]->Fill(pi0pt,mass);
		  pi_INVMASS[2]->Fill(pi0pt,mass);
		  
		  //if(pi0pt>4&&pi0pt<17){		  
		  if(mass>0.12 && mass<0.165){ 
		    pi_TRIGPT->Fill(pi0pt);
		    if(pi0pt>4) pi_asym[1]->Fill(asym12);

		    //First fill histos for true pi0 pairs
		    if(iparent2==iparent1 && parent1==111){
		      falsepair=1;


		      //Changed back to reco pt on 2/26 at 4:13
		      RECONPI_ZEMC_PTTRUE_PTRECO_real->Fill(zemc1,pipt1,pi0pt);
		      
		      PTpi_PTgam_TRUE->Fill(pipt1,pt1); //changed pipt1 to pi0pt 2/20/09
		      PTpi_PTgam_TRUE->Fill(pipt2,pt2); //changed pipt1 to pi0pt 2/20/09
		      
		    }else{
		      //For all not true
		      PTpi_PTgam_false->Fill(pi0pt,pt1);
		      PTpi_PTgam_false->Fill(pi0pt,pt2);

		      //If both are pi0
		      if(parent1==111 && parent2==111){
			PTpi_PTgam_False_PI0pair->Fill(pi0pt,pt1); 
			PTpi_PTgam_False_PI0pair->Fill(pi0pt,pt2); 
		      }else{
			//do it individually
			if(parent1==111){ 
			  PTpi_PTgam_PI0->Fill(pi0pt,pt1);
			}
			if(parent1==221) PTpi_PTgam_ETA->Fill(pi0pt,pt1);
			if(parent1==-111 || parent1==22) PTpi_PTgam_DIR->Fill(pi0pt,pt1);

			if(parent2==111){ 
			  PTpi_PTgam_PI0->Fill(pi0pt,pt2);
			}
			if(parent2==221) PTpi_PTgam_ETA->Fill(pi0pt,pt2);
			if(parent2==-111 || parent2==22) PTpi_PTgam_DIR->Fill(pi0pt,pt2);
			

			//RECONPI_ZEMC_PTTRUE_PTRECO_allfalse->Fill(zemc1,pipt1,pi0pt);
		      }

		    }
		 






		    if(parent1==221)eta_trueptvsrecopt->Fill(pi0pt,pipt1);

		    if(parent1==111){
		      
		      //if(iparent1==iparent2) RECONPI_ZEMC_RECONPT_REALPT->Fill(zemc1,pi0pt,pipt1);
		      
		      if(parent1==111 && iparent1==iparent2){
			realpi0_trueptvsrecopt->Fill(pi0pt,pipt1);
		      }else{
			pi0_trueptvsrecopt->Fill(pi0pt,pipt1);
		      }
		    }
		    if(parent1==-111 || parent1==22)direct_trueptvsrecopt->Fill(pi0pt,pipt1);
		  
		    //Now fill with second photon's pi0 info
		    if(parent2==221)eta_trueptvsrecopt->Fill(pi0pt,pipt2);
		    if(parent2==111){
		      if(parent1==111 && iparent1==iparent2){
			realpi0_trueptvsrecopt->Fill(pi0pt,pipt2);
		      }else{
			pi0_trueptvsrecopt->Fill(pi0pt,pipt2);
		      }
		    }
		    if(parent2==-111 || parent2==22)direct_trueptvsrecopt->Fill(pi0pt,pipt2);
		  
		    //And for the leading pt
		    double lead_pipt=pipt1;
		    double lead_parent=parent1;
		    if(pipt2>pipt1){
		      lead_pipt=pipt2;
		      lead_parent=parent2;
		    }
		    if(lead_parent==221)eta_leadptvsrecopt->Fill(pi0pt,lead_pipt);
		   
		    if(lead_parent==111 && iparent1==iparent2){
		      realpi0_leadptvsrecopt->Fill(pi0pt,lead_pipt);
		    }else{
		      pi0_leadptvsrecopt->Fill(pi0pt,lead_pipt);
		    }
		  
		    if(lead_parent==-111 || lead_parent==22) direct_leadptvsrecopt->Fill(pi0pt,lead_pipt);

		  
		  }
		  
		  //not in invmass window
		  
		}//gampart

	      //I don't think this made sense...
	      /*  
	      if(!falsepair && parent1==111){
		RECONPI_ZEMC_PTTRUE_PTRECO_false->Fill(zemc1,pipt1);
		PTpi_ZEMCpi_PTgam_false->Fill(pipt1,zemc1,pt1);
	      }
	      */

	    }//if accept1>0
	  //cout << "end of gamtrig loop" <<endl;
	}//gamtrig loop
      //cout <<"end tagging loop" <<endl;
    }//tagging
  
  //cout << "Filled Root Objects!! HERE" <<endl;  
  return;
}

