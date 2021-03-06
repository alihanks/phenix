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

  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run7_stripe.list";

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
  double px1, py1, pz1, px2, py2, pz2;
  double pt1, pt2, p1, p2;
  double phi1, phi2;
  double theta1, theta2, eta1, eta2;
  double ecore1, ecore2;
  double parent1, parent2;
  double pipt1, pipt2;

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
  float gampi0pt[nparts]; 
  float gamtheta[nparts];
  float gampt[nparts];
  float gampx[nparts];
  float gampy[nparts];
  float gampz[nparts];
  float gameta[nparts];

  //int eventcounter=0;

  for(int itest=1; itest<nparts; itest++)
    {
      gamenergy[itest]=0.0;
      gammom[nparts]=0.0;
      gamaccept[nparts]=0.0;
      gamphi[nparts]=0.0; 
      gamparent[nparts]=0.0; 
      gampi0pt[nparts]=0.0;  
      gamtheta[nparts]=0.0;
      gampt[nparts]=0.0;
      gampx[nparts]=0.0;
      gampy[nparts]=0.0;
      gampz[nparts]=0.0;
      gameta[nparts]=0.0;
    }
  int parentid=0;

  //cout <<"HERE we are in fill.cpp" <<endl;
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
	parentid = ipart;

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
	
	if (pid == 223 )
	  realomegapt->Fill(primpt);
	
	//      cout << "found a primary: " << pid << " wt " << endl;
      }//close if generation==1
    

    //cout <<"HERE before gamma loop" <<endl;

    if ( CurrentParticle->GetID()==22 	 || CurrentParticle->GetID()== -111  ){ // -111 = direct photon
      //cout << "HERE we found a photon" << endl;
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

	float phitemp = atan2(mom1.Getp().Getpy(),mom1.Getp().Getpx());
	if(phitemp<-PI/2.) phitemp+=2*PI;
	if(phitemp>3*PI/2.) phitemp-=2*PI;

	pxe_nosm = mom1.Getp().Getpx();
	pye_nosm = mom1.Getp().Getpy();
	pte_nosm = sqrt(pxe_nosm*pxe_nosm+pye_nosm*pye_nosm);
	
	if(phitemp<2.95)mom1 = ApplyEnergyResolution(mom1,0);
	else mom1 = ApplyEnergyResolution(mom1,1);

	float smeared_x=pxe;
	float smeared_y=pye;
	float smeared_z=pze;

	/*
	for(int jj=0;jj<8;jj++){
	  //testmom1=ApplyPositionResolution(jj,mom1,vtxz,deadfile);
	  int isaccept = phenix_geom(jj+10,(float)pxe,(float)pye,(float)pze,(float)vtxz,(float)pe,deadfile,smeared_x,smeared_y,smeared_z);
	  if(isaccept>0){
	    pxe=smeared_x;
	    pye=smeared_y;
	    pze=smeared_z;
	    mom1.Setp(pxe,pye,pze);
	  }
	}
	*/

	pe  = sqrt(mom1.GetE()*mom1.GetE()-mom1*mom1);
	pxe = mom1.Getp().Getpx();
	pye = mom1.Getp().Getpy();
	pze = mom1.Getp().Getpz();
	pte = sqrt(pxe*pxe+pye*pye);
	thetae = tan(pte/pze);
	etae = -log(tan(thetae/2.));
	cocktail->Fill(pte,1.0/pte);

	float junkx=0.0;
	float junky=0.0;
	float junkz=0.0;

	accept = -9;
	for(int j=0;j<8;j++){
	  int isaccept = phenix_geom(j,(float)pxe,(float)pye,(float)pze,(float)vtxz,(float)pe,deadfile,junkx,junky,junkz);
	  
	  if(isaccept>0) accept=isaccept;	  
	}
    
	accvar->Fill(accept);
 
	//simulate single particle efficiency
	if(gRandom->Rndm()>0.97) accept=-9;

	if(type<per&&pte>ptcut){
	  //gammapt->Fill(pte, fillweight);
	  gammapt->Fill(pte); //Not using fill weight

	  if(accept!=-9){
	    //Not using fill weights!
	    gamma_acc->Fill(pte);
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

	float phi=atan2(pye,pxe);
	gamenergy[ipart]=pe;
	//gammom[ipart]=mom1;
	gamaccept[ipart]=accept;
	gamphi[ipart]=phi; 
	gamparent[ipart]=parentid;
	gampi0pt[ipart]=primpt; 
	gamtheta[ipart]=thetae;
	gampt[ipart]=pte;
	gameta[ipart]=etae;

	gampx[ipart]=pxe;
  	gampy[ipart]=pye;
	gampz[ipart]=pze;
  
	if(pid==331)
	  etaprimegamma->Fill(pte);
	//if(pid==221)
	//etagamma->Fill(pte,fillweight);
	if(pid==223)
	  omegagamma->Fill(pte);

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
	float sume12 = ecore1+ecore2;
	//double zemc1= 510.0*pz1/pt1 + vtxz;
	//double zemc2= 510.0*pz2/pt2 + vtxz;

	float asym12 =  (ecore1 - ecore2)/ sume12; // fabs not neccesary because upper > check
	const float mina = 0.15;  //this is the smallest asym value used for the cut

	if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25))) { continue; }

	//cout << "******failed assym cut!!!" <<endl; continue;}	

	//changed accept1>-9&&accept2>-9 to >0
	//if(accept1==0||accept2==0) {cout << "aaaaaaaaaaaaa--accept=0"<<endl;}

	//if(fabs(zemcp) > 155) cout << "zzzzzzzzzzz cut on the zemcp" <<endl; 

      	if(accept1>0&&accept2>0&&ecoree>1.0&&ecorep>1.0&&sume12>=4.0&&fabs(zemcp) < 155) invmassvspt->Fill(primpt,mass);
      
     
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
    }  //cout <<"HERE out of the particle loop" <<endl;


    int tagging=1;

    if(tagging==1)
      {
	//cout <<"HERE in tagging loop"<<endl;
	for(int gamtrig=1; gamtrig<nparts; gamtrig++)
	  {
	    int taggedtrig=0;
	    int falsepair=0;

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
 		pipt1=gampi0pt[gamtrig];
 

		if(parent1==0 && ecore1<1) continue;

	    for(int gampart=1; gampart<nparts; gampart++)
	      {
		if(gampart==gamtrig) continue;
		if(gamparent[gampart]<1) continue;
			  
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
 		pipt2=gampi0pt[gampart];

		if(ecore2<1.0||ecore1<ecore2|| ecore1 == ecore2) continue;   
		  
		float mass= sqrt(2.0*ecore1*ecore2 - 2.0*(px1*px2+py1*py2+pz1*pz2)); 

		//make fiducial cuts  -- 155 cm for higher pt gamma
		double zemc1= 510.0*pz1/pt1 + vtxz;
		double zemc2= 510.0*pz2/pt2 + vtxz;

		float dphi = fabs(phi1-phi2);

		//if(dphi>PI/2) continue;
		if(ecore2<0.5) continue;
		if(accept2<=0)continue;
		
		if(accept1<= 0) continue;
		if(fabs(zemc1) > 155) continue;

		  float sume12 = ecore1+ecore2;
		  if (sume12 < 4.0) continue;
		  
		  float asym12 =  (ecore1 - ecore2)/ sume12; // fabs not neccesary because upper > check
		  const float mina = 0.15;  //this is the smallest asym value used for the cut
		  
		  if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25)))continue;		 


		  //cout << "HERE we're gonna fill INVMASS histos" <<endl;
		  invmass->Fill(mass);

		  TLorentzVector pho1, pho2, pi0;
		  pho1.SetPxPyPzE(px1,py1,pz1,ecore1);
		  pho2.SetPxPyPzE(px2,py2,pz2,ecore2);
		  pi0=pho1+pho2;

		  float pi0pt=pi0.Pt();		 

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
		  

		if(mass>0.108 && mass<0.165){ 
		  //TAGGED!
		  //cout << "TAGGED!!!" <<endl;			

		  TotalTagged->Fill(pi0pt);

		  if(pi0pt>5 && pi0pt<7) 
		    {
		      //Only fill with the higher pT pi0
		      if(pipt1>pipt2)
			{
			  pi0Tagged->Fill(pipt1);
			}else{
			  pi0Tagged->Fill(pipt2);
			}
		  if(parent2==parent1)
		    { 
		      if(pipt1!=pipt2) cout << "parent pt mismatch!!!!" <<endl;
		      //TrueTag->Fill(pt1);
		      TrueTag->Fill(pipt1);
		    }
			}
		  if(parent2!=parent1)
		    { 
		      if(pipt1>pipt2)
			{
			  FalseTag->Fill(pipt1);
			}else{
			  FalseTag->Fill(pipt2);
			}

		  /*
		  if(falsepair<1 && parent2!=parent1)
		    {
		      FlaseTag->Fill(pipt1);
		      falsepair=1;
		    }
		  */
		    }
		  if(taggedtrig<1)
		    {
		      Ntag->Fill(pt1);
		      taggedtrig=1;
		    }
		}
	      }//gampart loop
	      }//if accept1>0
	    //cout << "end of gamtrig loop" <<endl;
	  }//gamtrig loop
	//cout <<"end tagging loop" <<endl;
      }//tagging
    //cout << "Filled Root Objects!! HERE" <<endl;
    return;
}

