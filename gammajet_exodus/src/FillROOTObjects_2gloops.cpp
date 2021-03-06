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

#define PI  3.141592653589793238

#define  INCLUDEFLAG extern
#include "DeclareROOTObjects.h"

bool   ParticleIsLepton(Particle*);
double thetaof(Mom3);
double phiof(Mom3);
Mom4   ApplyResolution(int, Mom4);
Mom4   ApplyEnergyResolution(Mom4, int);
int    PHENIXFilter(Particle*, ParticlePropertyList*);
bool   PHENIXPhotonFilter(Particle* kp, int fid = 0);
bool   PairInAcceptance(int, Mom4, Mom4);
int phenix_geom(int,float,float,float,float,char*);

void FillROOTObjects(int setup, double dNdy_pi0, double N_coll,
		     char *output_file, 
		     ParticleList *PList, 
		     ParticlePropertyList *PPList)
{

  //  gRandom = new TRandom3();
  //  gRandom->SetSeed(0);
    
  //char deadfile[]="/direct/phenix+workarea/manguyen/gammajet.exodus/deadmap_test.lst";
  //char deadfile[]="/phenix/workarea/manguyen/offline/analysis/run5pp_photon/offline/macros/jobs5/hotdeadmap.txt"; 
  //  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/hotdeadmap2.list";  
  //  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/hotdeadmap_run7_taxi76.list";  
  //  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/hotdeadmap_run4_hpdsts.list";  
  //char emcalfile[]="/direct/phenix+workarea/manguyen/gammajet.exodus/deadmap.lst.isobe";

  //using live files instead of deadfiles now
  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run5_taxi50.list";
  //      char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run6_taxi71.list";
  //char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run4_hpdsts.list";
  //  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run7_stripe.list";

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
  double vtxz;
  double theta;
  double thetae, thetap, etae, etap;
  double ye,yp;

  //for prime particle
  int 	 pid = 0;
  Mom4 primmom1;
  double primweight, primfillweight, primpt, primp;

  // z-position in the emcal plane
  double primzemc;

  n_pi0_generated = 10.0e+6;
  sigma_inelastic = 42.2;
  abs_norm        = dNdy_pi0; // for AuAu: multiplicity
  // abs_norm        = dNdy_pi0 * sigma_inelastic; // for pp: cross section

  PLNode   * CurrentNode     = PList->GetHeadNode();
  PLNode   * NextNode        = 0;
  Particle * CurrentParticle = 0;
  Particle * NextParticle    = 0;
  Particle * CurrentPrimary  = 0;
  
  // here we start with the diphoton analysis
  static int  neventK = 0;


  for ( int ipart=1; ipart<=PList->GetLength(); ipart++) 
    {
      //cout<<" ipart "<<ipart<<endl;
    CurrentNode = CurrentNode->GetNextNode();
    CurrentParticle = CurrentNode->Get(0);

    //cout << CurrentParticle->GetID() << ":" << CurrentParticle->GetGeneration() << endl;

    //Particle of generation=1 is always immediately
    //followed by it's decay particles so this format works    
    if ( CurrentParticle->GetGeneration()==1 )
      {
	CurrentPrimary = CurrentParticle;
	pid            = CurrentParticle->GetID();

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

	//double primtheta = acos(primmom1.Getp().Getpz()/primp);
	//cout<<" vtxz pi0 "<<vtxz<<endl;
	//use pbsc as ref radius
	primzemc  = 510.0*primmom1.Getp().Getpz()/primpt+vtxz;
	//primzemc = 510.0*cos(primtheta)+vtxz;

	/*
	if (pid == 111)
	  realpi0pt->Fill(primpt, primfillweight);
      
	if (pid == 221)
	  realetapt->Fill(primpt, primfillweight);
	
	if (pid == 331)
	  realetaprimept->Fill(primpt, primfillweight);
	
	if (pid == 223 )
	  realomegapt->Fill(primpt, primfillweight);
	*/


	//Do not weight the input particles
	if (pid == 111)
	  realpi0pt->Fill(primpt);
      
	if (pid == 221)
	  realetapt->Fill(primpt);
	
	if (pid == 331)
	  realetaprimept->Fill(primpt);
	
	if (pid == 223 )
	  realomegapt->Fill(primpt);
	
	//      cout << "found a primary: " << pid << " wt " << endl;
      }//close if generation==1
    
    if ( CurrentParticle->GetID()==22 	 || CurrentParticle->GetID()== -111  ){ // -111 = direct photon

    //if ( CurrentParticle->GetID()== 111  ) 
      
	
	
      // cout << "found a photon" << endl;
	// correct weighting for pairs (if decay chains are full implemented!) 
	weight   = CurrentParticle->GetWeight();
	fillweight = abs_norm*weight/(0.5*n_pi0_generated);
	fillweight = fillweight/0.025;
	// photon 4-vector
	mom1 = CurrentParticle->Get4mom();
	//vtxz = CurrentParticle->GetzVertex();
	//cout<<" vtxz gamma "<<vtxz<<endl;
	type = gRandom->Rndm();
	// resolution

	float phitemp = atan2(mom1.Getp().Getpy(),mom1.Getp().Getpy());
	if(phitemp<-PI/2.) phitemp+=2*PI;
	if(phitemp>3*PI/2.) phitemp-=2*PI;

	pxe_nosm = mom1.Getp().Getpx();
	pye_nosm = mom1.Getp().Getpy();
	pte_nosm = sqrt(pxe_nosm*pxe_nosm+pye_nosm*pye_nosm);

	
	if(phitemp<2.95)mom1 = ApplyEnergyResolution(mom1,0);
	else mom1 = ApplyEnergyResolution(mom1,1);


	
	// acceptance filter
	//accept1 = true;
	//	accept1 = PHENIXPhotonFilter(CurrentParticle);
	//	bool accept1_fid = PHENIXPhotonFilter(CurrentParticle, 1);
      
	pe  = sqrt(mom1.GetE()*mom1.GetE()-mom1*mom1);
	pxe = mom1.Getp().Getpx();
	pye = mom1.Getp().Getpy();
	pze = mom1.Getp().Getpz();
	pte = sqrt(pxe*pxe+pye*pye);
	thetae = tan(pte/pze);
	etae = -log(tan(thetae/2.));
	cocktail->Fill(pte,1.0/pte);

	/*
	accept = -1;
	for(int j=0;j<8;j++){
	  if(phenix_geom(j,(float)pxe,(float)pye,(float)pze,(float)vtxz,deadfile)) accept=j;
	  
	}
	*/
	
	//complicated.  -9 means did not hit emcal.  all other negative means it hit emcal but near at a bad ert tile (p+p only) -man
	

	accept = -9;
	for(int j=0;j<8;j++){
	  int isaccept = phenix_geom(j,(float)pxe,(float)pye,(float)pze,(float)vtxz,deadfile);
	  //if(isaccept==1) accept=j+1;
	  //if(isaccept==-1) accept=-j-1;
	  if(isaccept>0) accept=isaccept;
	  //cout<<" isaccept "<<isaccept<<" accept "<<accept<<endl;
	  
	}
	

	//simulate single particle efficiency
	if(gRandom->Rndm()>0.97) accept=-9;

	if(type<per&&pte>ptcut){
	  //gammapt->Fill(pte, fillweight);
	  gammapt->Fill(pte); //Not using fill weight

	  if(accept!=-9){

	    //Not using fill weights!
	    //gamma_acc->Fill(pte,fillweight);
	    //gamma_acc_nosm->Fill(pte_nosm,fillweight);

	    gamma_acc->Fill(pte);
	    gamma_acc_nosm->Fill(pte_nosm);

	    float phi=atan2(pye,pxe);
	    if(phi<-PI/2.)phi+=2*PI;
	    if(phi>3*PI/2.)phi-=2*PI;
	    //if(accept>0)gamma_acc_phi->Fill(phi);
	    //if(accept==0)gamma_acc_phi->Fill(phi);

	    
	  }
	}
    
	if(pid==111||pid==221)
	  {
	    //cout<<" pid "<<pid<<endl;
	    if(type<per&&pte>ptcut){
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
	    //cout<<" pi0pair "<<pi0pair<<endl;
	  }


	//Not using fill weights
	/*     	
	if(pid==331)
	  etaprimegamma->Fill(pte,fillweight);
	//if(pid==221)
	//etagamma->Fill(pte,fillweight);
	if(pid==223)
	  omegagamma->Fill(pte,fillweight);
	*/

	if(pid==331)
	  etaprimegamma->Fill(pte);
	//if(pid==221)
	//etagamma->Fill(pte,fillweight);
	if(pid==223)
	  omegagamma->Fill(pte);



	///////do the tagging loop here////////////////////
	//match will every other gamma in event
	//count how many times each gamma matches
	//count how many times each gamma makes a false match
	//both as a function of pT
	//get Ntag also in pT bins

	PLNode   * PartnerNode     = PList->GetHeadNode();
	//PLNode   * NextNode        = 0;
	Particle * PartnerParticle = 0;
	//Particle * NextParticle    = 0;
	//Particle * CurrentPrimary  = 0;
	int partpid = 0;

	for ( int ippart=1; ippart<=PList->GetLength(); ippart++) 
	  {
	    PartnerNode = PartnerNode->GetNextNode();
	    PartnerParticle = PartnerNode->Get(0);
	    
	    //cout << CurrentParticle->GetID() << ":" << CurrentParticle->GetGeneration() << endl;
	    
	    partpid = PartnerParticle->GetID();
	    if ( CurrentParticle->GetID()==22 	 || CurrentParticle->GetID()== -111  ){ // -111 = direct photon

	      //if ( CurrentParticle->GetID()== 111  ) 
	      
	      // photon 4-vector
	      mom1 = CurrentParticle->Get4mom();
	      //vtxz = CurrentParticle->GetzVertex();
	      //cout<<" vtxz gamma "<<vtxz<<endl;
	      type = gRandom->Rndm();
	      // resolution
	      
	      float phitemp = atan2(mom1.Getp().Getpy(),mom1.Getp().Getpy());
	      if(phitemp<-PI/2.) phitemp+=2*PI;
	      if(phitemp>3*PI/2.) phitemp-=2*PI;
	      
	      pxe_nosm = mom1.Getp().Getpx();
	      pye_nosm = mom1.Getp().Getpy();
	      pte_nosm = sqrt(pxe_nosm*pxe_nosm+pye_nosm*pye_nosm);
	      
	      
	      if(phitemp<2.95)mom1 = ApplyEnergyResolution(mom1,0);
	      else mom1 = ApplyEnergyResolution(mom1,1);
	      
	      
	      
	      // acceptance filter
	      //accept1 = true;
	      //	accept1 = PHENIXPhotonFilter(CurrentParticle);
	      //	bool accept1_fid = PHENIXPhotonFilter(CurrentParticle, 1);
      
	      pe  = sqrt(mom1.GetE()*mom1.GetE()-mom1*mom1);
	      pxe = mom1.Getp().Getpx();
	      pye = mom1.Getp().Getpy();
	      pze = mom1.Getp().Getpz();
	      pte = sqrt(pxe*pxe+pye*pye);
	      thetae = tan(pte/pze);
	      etae = -log(tan(thetae/2.));
	      cocktail->Fill(pte,1.0/pte);

	      /*
		accept = -1;
		for(int j=0;j<8;j++){
		if(phenix_geom(j,(float)pxe,(float)pye,(float)pze,(float)vtxz,deadfile)) accept=j;
		
		}
	      */
	      
	      //complicated.  -9 means did not hit emcal.  all other negative means it hit emcal but near at a bad ert tile (p+p only) -man
	      
	      
	      accept = -9;
	      for(int j=0;j<8;j++){
		int isaccept = phenix_geom(j,(float)pxe,(float)pye,(float)pze,(float)vtxz,deadfile);
		//if(isaccept==1) accept=j+1;
		//if(isaccept==-1) accept=-j-1;
		if(isaccept>0) accept=isaccept;
		//cout<<" isaccept "<<isaccept<<" accept "<<accept<<endl;
		
	      }
	

	//simulate single particle efficiency
	if(gRandom->Rndm()>0.97) accept=-9;

	if(type<per&&pte>ptcut){
	  //gammapt->Fill(pte, fillweight);
	  gammapt->Fill(pte); //Not using fill weight

	  if(accept!=-9){

	    //Not using fill weights!
	    //gamma_acc->Fill(pte,fillweight);
	    //gamma_acc_nosm->Fill(pte_nosm,fillweight);

	    gamma_acc->Fill(pte);
	    gamma_acc_nosm->Fill(pte_nosm);

	    float phi=atan2(pye,pxe);
	    if(phi<-PI/2.)phi+=2*PI;
	    if(phi>3*PI/2.)phi-=2*PI;
	    //if(accept>0)gamma_acc_phi->Fill(phi);
	    //if(accept==0)gamma_acc_phi->Fill(phi);

	    
	  }
	}
    
	if(pid==111||pid==221)
	  {
	    //cout<<" pid "<<pid<<endl;
	    if(type<per&&pte>ptcut){
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
	    //cout<<" pi0pair "<<pi0pair<<endl;
	  }


	//Not using fill weights
	/*     	
	if(pid==331)
	  etaprimegamma->Fill(pte,fillweight);
	//if(pid==221)
	//etagamma->Fill(pte,fillweight);
	if(pid==223)
	  omegagamma->Fill(pte,fillweight);
	*/

	if(pid==331)
	  etaprimegamma->Fill(pte);
	//if(pid==221)
	//etagamma->Fill(pte,fillweight);
	if(pid==223)
	  omegagamma->Fill(pte);

    }







	  }



    }//if particle is a gamma
    
    
    if(pi0pair==1)
      {
	//cout<<" hello "<<endl;
	pxe = pairmom[0].Getp().Getpx();
	pye = pairmom[0].Getp().Getpy();
	pze = pairmom[0].Getp().Getpz();
	pte = sqrt(pxe*pxe+pye*pye);
	pe = sqrt(pxe*pxe+pye*pye+pze*pze);

	phie =  atan2(pye,pxe);
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
	
	//cout<<" accept1 "<<accept1<<" accept2 "<<accept2<<endl;

	type1=pairtype[0];
	type2=pairtype[1];
	float ecoree=sqrt(pxe*pxe+pye*pye+pze*pze);
	float ecorep=sqrt(pxp*pxp+pyp*pyp+pzp*pzp);

	float mass= sqrt(2.0*ecoree*ecorep - 2.0*(pxe*pxp+pye*pyp+pze*pzp)); 

	//filling these without weights--> flat distribution

	//	ptpivsptgam->Fill(pte,sqrt((pxe+pxp)*(pxe+pxp)+(pye+pyp)*(pye+pyp)));
	//	ptpivsptgam->Fill(ptp,sqrt((pxe+pxp)*(pxe+pxp)+(pye+pyp)*(pye+pyp)));

	//make sure one of the cluster fired the trigger.  requires that one cluster with ecore > 3 hit a live area  ert threshold ~ 3 GeV

	//Not Using fill weights
	/*
	if(accept1>-9&&pte>5.0)gamma_acc_phi->Fill(phie,fillweight1);
	if(accept2>-9&&ptp>5.0)gamma_acc_phi->Fill(phip,fillweight2);

	if(accept1>0&&pte>5.0)gamma_acc_phi_live->Fill(phie,fillweight1);
	if(accept2>0&&ptp>5.0)gamma_acc_phi_live->Fill(phip,fillweight2);
	*/

	if(accept1>-9&&pte>5.0)gamma_acc_phi->Fill(phie);
	if(accept2>-9&&ptp>5.0)gamma_acc_phi->Fill(phip);

	if(accept1>0&&pte>5.0)gamma_acc_phi_live->Fill(phie);
	if(accept2>0&&ptp>5.0)gamma_acc_phi_live->Fill(phip);


	int firedert=0;
	//if(accept1>0&&accept2>0) firedert=1;

	//screw the ert for now
	//if(accept1>0&&ecoree>3)firedert=1;
	//if(accept2>0&&ecorep>3)firedert=1;

	//else if(accept1>-9&&accept1<0&&accept2>0&&ecorep>3) firedert=1;
	//else if(accept2>-9&&accept2<0&&accept1>0&&ecoree>3) firedert=1;
	

	//if(firedert==0) cout<<" happens ? "<<endl;


	//make fiducial cuts  -- 155 cm for higher pt gamma
	double zemcp= 510.0*pzp/ptp + vtxz;
	double zemce= 510.0*pze/pte + vtxz;

		
	//ignore ert for now:
	firedert=1;
	//if(firedert ==1 && accept> 0 && leadzemc < 155&& subzemc< 165){
	if(firedert ==1 && accept1> 0 && fabs(zemce) < 155){
	  
	  //if(primpt> 5.0 && primzemc > 165){
	  //cout<<" PTpi "<<primpt<<" pte "<<pte<<" ptp "<<ptp<<endl;
	  //}
	  
	  ptpivsptgam->Fill(pte,primpt,(float)accept1);
	  PTpi_ZEMCpi_PTgam->Fill(primpt,primzemc,pte,(float)accept1);
	  ppivspgam->Fill(pe,primp,(float)accept1);

	  //Not using fill weight
	  //if(accept1>-9&&pte>5.0)gamma_acc_ert_phi->Fill(phie,fillweight1);
	  if(accept1>-9&&pte>5.0)gamma_acc_ert_phi->Fill(phie);

	  if(accept2==-9||ptp<0.5) PTpi_ZEMCpi_PTgam_MISS->Fill(primpt,primzemc,pte,(float)accept1);
	}
	
	if(firedert ==1 && accept2> 0 && fabs(zemcp) < 155){
	  
	  //if(primpt> 5.0 && primzemc > 165){
	  //cout<<" PTpi "<<primpt<<" pte "<<pte<<" ptp "<<ptp<<endl;
	  //}
	  
	  ptpivsptgam->Fill(ptp,primpt,(float)accept2);
	  PTpi_ZEMCpi_PTgam->Fill(primpt,primzemc,ptp,(float)accept2);	  
	  ppivspgam->Fill(pp,primp,(float)accept2);

	  //Not Using Fill Weight
	  //if(accept2>-9&&ptp>5.0)gamma_acc_ert_phi->Fill(phip,fillweight2);
	  if(accept2>-9&&ptp>5.0)gamma_acc_ert_phi->Fill(phip);

	  if(accept1==-9||pte<0.5) PTpi_ZEMCpi_PTgam_MISS->Fill(primpt,primzemc,pte,(float)accept2);
	}
	
	if(firedert==1&&accept1>0&&accept2>0&&ecoree>1&&ecorep>1){
	  //here use the greater ert efficiency
	  if(accept1>accept2)RECONPI_ZEMC_PT->Fill(primzemc,primpt,(float)accept1);
	  else RECONPI_ZEMC_PT->Fill(primzemc,primpt,(float)accept2);
	}

	

	//if(accept1>-1&&(accept2<0||ecorep<0.5)) ptpivsptgam_tag->Fill(pte,sqrt((pxe+pxp)*(pxe+pxp)+(pye+pyp)*(pye+pyp)));
	//	if(accept2>-1&&(accept1<0||ecoree<0.5)) ptpivsptgam_tag->Fill(ptp,sqrt((pxe+pxp)*(pxe+pxp)+(pye+pyp)*(pye+pyp)));
	
	
	if(accept1>-9&&(accept2<-8||ecorep<1.0)) {
	  ppivspgam_tag->Fill(pe,primp);
	  ppivspgam_tag->Fill(pe,primp);
	}
	if(accept2>-9&&(accept1<-8||ecoree<1.0)) {
	  ppivspgam_tag->Fill(pp,primp);
	  ppivspgam_tag->Fill(pp,primp);
	}
	invmass->Fill(mass);
      
    //cout << pxe << pye << pze << endl;
    //	  cout << NextParticle->GetID() << "Accept1:" << accept1 << "Vtxz:" << vtxz << endl;
    
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
    }

  return;

}

