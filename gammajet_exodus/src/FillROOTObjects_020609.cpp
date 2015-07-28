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
    
  //char deadfile[]="/direct/phenix+workarea/manguyen/gammajet.exodus/deadmap_test.lst";
  //char deadfile[]="/phenix/workarea/manguyen/offline/analysis/run5pp_photon/offline/macros/jobs5/hotdeadmap.txt"; 
  //  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/hotdeadmap2.list;  
  //  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/hotdeadmap_run7_taxi76.list";  
  //  char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/hotdeadmap_run4_hpdsts.list";  
  //char emcalfile[]="/direct/phenix+workarea/manguyen/gammajet.exodus/deadmap.lst.isobe";

  //using live files instead of deadfiles now
  //char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run5_taxi50.list";
  //      char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run6_taxi71.list";
  //char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run4_hpdsts.list";
  char deadfile[]="/direct/phenix+u/workarea/mjuszkie/gammajet_exodus/livetowers_r4all.txt";
  //char deadfile[]="/phenix/u/workarea/manguyen/gammajet.exodus/livetowers_run7_stripe.list";

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
	
	//if (pid == 223 ) //Now filling omega histo with photon/direct photons
	if(pid==22 || pid==-111)
	  realomegapt->Fill(primpt);
	
	//      cout << "found a primary: " << pid << " wt " << endl;
      }//close if generation==1
    
    
    //cout <<"HERE before gamma loop" <<endl;

    if ( CurrentParticle->GetID()==22 	 || CurrentParticle->GetID()== -111  ){ // -111 = direct photon

    //if ( CurrentParticle->GetID()== 111  ) 
      
	
	
      //cout << "HERE we found a photon" << endl;
	// correct weighting for pairs (if decay chains are full implemented!) 
	weight   = CurrentParticle->GetWeight();
	fillweight = abs_norm*weight/(0.5*n_pi0_generated);
	fillweight = fillweight/0.025;
	// photon 4-vector
	mom1 = CurrentParticle->Get4mom();
	//vtxz = CurrentParticle->GetzVertex();
	//cout<<" vtxz gamma "<<vtxz<<endl;
	//type = gRandom->Rndm();
	type=1;
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
		/*
		  if(fabs(mom1.GetE()-sqrt(pxe*pxe+pye*pye+pze*pze))>.01) 
		  cout << "E_before_smearing " << mom1.GetE() << " after " << sqrt(pxe*pxe+pye*pye+pze*pze) <<endl;
		*/
		
		mom1.Setp(pxe,pye,pze);
		mom1.SetE(sqrt(pxe*pxe+pye*pye+pze*pze));
		
		//These differ as one would guess
		//cout << "E_after_smearing" << mom1.GetE() << " or " << sqrt(pxe*pxe+pye*pye+pze*pze) <<endl;
	      }
	    
	    //cout<< "after_smearing " << pxe << " " << pye << " " << pze <<endl;
	  }
	}

	
	/*
	float phitemp = atan2(mom1.Getp().Getpy(),mom1.Getp().Getpx());
	//fixed bug below that caused phitemp to always be pi/4 or -3pi/4
	//float phitemp = atan2(mom1.Getp().Getpy(),mom1.Getp().Getpy());
	//phitemp is only used for the energy resolution
	//cout<< "phitemp" << phitemp <<endl;
	if(phitemp<-PI/2.) phitemp+=2*PI;
	if(phitemp>3*PI/2.) phitemp-=2*PI;

	if(phitemp<2.95)mom1 = ApplyEnergyResolution(mom1,0);
	else mom1 = ApplyEnergyResolution(mom1,1);
	*/


	/*
	int nochange=0;
	for(int jj=0;jj<8;jj++){
	  testmom1=ApplyPositionResolution(jj,mom1,vtxz,deadfile);
	    double pxe_test = testmom1.Getp().Getpx();
	    double pye_test = testmom1.Getp().Getpy();
	    
	  if(pxe_nosm==pxe_test && pye_nosm==pye_test){
	    nochange++;
	    //if(nochange==8) cout<<"no change to mom"<<endl;
	  }else{
	    //cout << "mom was "<< pxe_nosm << " " << pye_nosm <<endl;
	    //cout << "mom is  "<< pxe_test << " " << pye_test <<endl;
	    mom1=testmom1;
	    break;
	  }
	    
	}	
	*/  

	
	// acceptance filter
	//accept1 = true;
	//	accept1 = PHENIXPhotonFilter(CurrentParticle);
	//	bool accept1_fid = PHENIXPhotonFilter(CurrentParticle, 1);
	//Ee=mom1.GetE();

	//if(mom1.GetE()*mom1.GetE()<mom1*mom1) cout << "nan!!!!" <<endl;
	//pe  = sqrt(mom1.GetE()*mom1.GetE()-mom1*mom1);
	pe=mom1.GetE();

	//These are all equivalent
	//cout << "E= " << mom1.GetE() << "  pe_sub= " << pe << " pe_add= "<< sqrt(pxe*pxe+pye*pye+pze*pze) <<endl;

	//if(mom1.GetE()!=pe || pe!=sqrt(pxe*pxe+pye*pye+pze*pze)) cout << "E= " << mom1.GetE() << "  pe_sub= " << pe << " pe_add= "<< sqrt(pxe*pxe+pye*pye+pze*pze) <<endl;

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
	

	//cout << "HERE we are doing the acceptance & phenix geo" <<endl;

	//Already filled accept for acceptance above with pos res
	/*
	float junkx=0.0;
	float junky=0.0;
	float junkz=0.0;

	accept = -9;
	for(int j=0;j<8;j++){
	  int isaccept = phenix_geom(j,(float)pxe,(float)pye,(float)pze,(float)vtxz,(float)pe,deadfile,junkx,junky,junkz);
	  //if(isaccept==1) accept=j+1;
	  //if(isaccept==-1) accept=-j-1;
	  
	  if(isaccept>0) accept=isaccept;
	  //cout<<" isaccept "<<isaccept<<" accept "<<accept<<endl;
	  //if(isaccept>0) accvar->Fill(isaccept);
	  
	}
	*/
    
	//cout << "HERE we've gotten the acceptance & phenix geo" <<endl;

	//simulate single particle efficiency

	//if(accept>0) cout << "geom_passed" << pid <<endl;
	//else cout << "geom_failed_parent" << pid <<endl;
    
	float eff=0.95;
	//if(pte<20.0) eff=0.0075*pte+0.8; //for centrality 0-10
	//Maybe 80% is too low?

	if(gRandom->Rndm()>eff) accept=-9;

	//if(accept>0) cout << "eff_passed" <<endl;

	accvar->Fill(accept);
	//if(pte<ptcut) cout << "ptcut actually does something" <<endl;
	//if(type>per) cout << "typecut actually does something" <<endl;

	if(type<per&&pte>ptcut){
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

	//cout <<"HERE after pi0 check" <<endl;
	

	float phi=atan2(pye,pxe);
	gamenergy[ipart]=pe;
	//gammom[ipart]=mom1;
	gamaccept[ipart]=accept;
	gamphi[ipart]=phi; 
	gamindex[ipart]=parentid;
	gamzemc[ipart]=primzemc;
	//if(pid==-111) cout << "ipart " << ipart << " from parentid: " << pid <<endl;
	//cout << "ipart " << ipart << " from parentid: " << gamindex[ipart] <<endl;


	if(pe>1.0 && primpt<0.65) cout << "gampt " <<pe << " from pi0pt:" << primpt <<endl;

	gamparent[ipart]=pid;
	gampi0pt[ipart]=primpt; 
	gamtheta[ipart]=thetae;
	gampt[ipart]=pte;
	gameta[ipart]=etae;

	gampx[ipart]=pxe;
  	gampy[ipart]=pye;
	gampz[ipart]=pze;
  
	gammapt->Fill(pte);


	//cout <<"HERE after filling gamma arrays" <<endl;


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


	//cout <<"HERE at the end of the gamma loop" <<endl;

	///////do the tagging loop here////////////////////
	//match will every other gamma in event
	//count how many times each gamma matches
	//count how many times each gamma makes a false match
	//both as a function of pT
	//get Ntag also in pT bins


    }//if particle is a gamma
    
    //cout <<"HERE we are out of the gamma loop pi0pair=" << pi0pair <<endl;

    if(pi0pair==1)
      {
	//cout<<" pi0pair!!!!! "<<endl;
	pxe = pairmom[0].Getp().Getpx();
	pye = pairmom[0].Getp().Getpy();
	pze = pairmom[0].Getp().Getpz();
	pte = sqrt(pxe*pxe+pye*pye);

	pe = sqrt(pxe*pxe+pye*pye+pze*pze);
	//These are equal if we set the E in mom1 after smearing
	//if((pe-pairmom[0].GetE())!=0) cout << "E= " << pairmom[0].GetE() << "  pe= " << pe << endl;

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

	//cout << "filled mass, now filling TLorentzVectors " <<endl;

	TLorentzVector pho1, pho2, pi0;
	pho1.SetPxPyPzE(pxe,pye,pze,ecoree);
	pho2.SetPxPyPzE(pxp,pyp,pzp,ecorep);
	pi0=pho1+pho2;

	float pi0pt=pi0.Pt();		 

      
	if(accept1>0&&accept2>0) pi0gammapt[2]->Fill(pi0pt);
	else if(accept1>0 || accept2>0) pi0gammapt[1]->Fill(pi0pt);
	else pi0gammapt[0]->Fill(pi0pt);



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

	  //Remove 3D histos to save space
	  
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

	  //if(accept1==-9||pte<0.5) PTpi_ZEMCpi_PTgam_MISS->Fill(primpt,primzemc,pte,(float)accept2);
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

	//invmass->Fill(mass);
      

	//Apply cuts that match data
	float sume12 = ecoree+ecorep;
	//double zemc1= 510.0*pz1/pt1 + vtxz;
	//double zemc2= 510.0*pz2/pt2 + vtxz;

	float asym12 =  (ecoree - ecorep)/ sume12; 

	const float mina = 0.15;  //this is the smallest asym value used for the cut

	int passasym=1;
	if (sume12 < 5.25 && fabs(asym12) > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25))) { passasym=0; }
      
	//cout << "******failed assym cut!!!" <<endl; continue;}	

	//changed accept1>-9&&accept2>-9 to >0
	//if(accept1==0||accept2==0) {cout << "aaaaaaaaaaaaa--accept=0"<<endl;}

	//if(fabs(zemcp) > 155) cout << "zzzzzzzzzzz cut on the zemcp" <<endl; 

      	//if(accept1>0&&accept2>0&&ecoree>1.0&&ecorep>1.0&&sume12>=4.0&&fabs(zemcp) < 155) invmassvspt->Fill(primpt,mass);

	/*
	  cout<< "pi0pair_b4_cuts" <<endl;
	  if(accept1>0&&accept2>0) cout << "pi0pair_accept_ok" <<endl;      
	  if(accept1>0&&accept2>0&&ecoree>1.0&&ecorep>1.0) cout << "pi0pair_ecore_ok" <<endl; 
	  if(accept1>0&&accept2>0&&ecoree>1.0&&ecorep>1.0&&sume12>=4.0) cout << "pi0pair_sume12_ok" <<endl; 
	  if(accept1>0&&accept2>0&&ecoree>1.0&&ecorep>1.0&&sume12>=4.0&&fabs(zemcp) < 155) cout << "pi0pair_zemc_ok" <<endl; 
	*/

	if(accept1>0&&accept2>0&&ecoree>1.0&&ecorep>1.0&&sume12>=4.0&&fabs(zemcp) < 155 && passasym==1){ 
	  invmassvspt->Fill(pi0pt,mass);
	  //cout << "mass: " << mass << " pi0pt: " << pi0pt << " from " << ecoree << " and " << ecorep <<endl;
	  //cout<< "pi0pair_ok" <<endl;
	}
      
      
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
    //cout <<"HERE out of the pi0pair loop" <<endl;

    }

  //cout <<"HERE out of the particle loop" <<endl;


    int tagging=1;

    if(tagging==1)
      {
	//cout <<"HERE in tagging loop with " << nparts <<endl;
	for(int gamtrig=1; gamtrig<nparts; gamtrig++)
	  {
	    int taggedtrig=0;
	    int falsepair=0;

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
		      //if(gamparent[gampart]<1) continue;
		      /*
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
		      */
		      
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


		      if(ecore2<1.0||ecore1<ecore2|| ecore1 == ecore2) continue;   
		      
		     
		      //cout<<" accept1 "<<accept1<<" accept2 "<<accept2<<endl;
		      
		      //float ecoree=sqrt(pxe*pxe+pye*pye+pze*pze);
		      //float ecorep=sqrt(pxp*pxp+pyp*pyp+pzp*pzp);
		      
		      //float mass= sqrt(2.0*ecoree*ecorep - 2.0*(pxe*pxp+pye*pyp+pze*pzp)); 
		      float mass= sqrt(2.0*ecore1*ecore2 - 2.0*(px1*px2+py1*py2+pz1*pz2));

		     		      //cout<<"mass: " << mass << " E1: " << ecore1 << " E2: " << ecore2 << " px1: " << px1 << " px2: " << px2 << " second term:" << 2.0*(px1*px2+py1*py2+pz1*pz2) <<endl;
		      
		      //if(accept1>-9&&pte>5.0)gamma_acc_phi->Fill(phie);
		      //if(accept2>-9&&ptp>5.0)gamma_acc_phi->Fill(phip);
		      
		      //if(accept1>0&&pte>5.0)gamma_acc_phi_live->Fill(phie);
		      //if(accept2>0&&ptp>5.0)gamma_acc_phi_live->Fill(phip);
		      
		      
		      
		      //make fiducial cuts  -- 155 cm for higher pt gamma
		      zemc1= 510.0*pz1/pt1 + vtxz;
		      zemc2= 510.0*pz2/pt2 + vtxz;
		    
		      
		      //ignore ert for now:
		      
		      //if(accept1> 0 && fabs(zemc1) < 155){
		      
		      //if(primpt> 5.0 && primzemc > 165){
		      //cout<<" PTpi "<<primpt<<" pte "<<pte<<" ptp "<<ptp<<endl;
		      //}
		      
		      //not fill histos here
		      //ptpivsptgam->Fill(pte,primpt,(float)accept1);
		      //PTpi_ZEMCpi_PTgam->Fill(primpt,primzemc,pte,(float)accept1);
		      //ppivspgam->Fill(pe,primp,(float)accept1);
		      
		      
		      //Not using fill weight
		      //if(accept1>-9&&pte>5.0)gamma_acc_ert_phi->Fill(phie,fillweight1);
		      //if(accept1>-9&&pte>5.0)gamma_acc_ert_phi->Fill(phie);
		      
		      //if(accept2==-9||ptp<0.5) PTpi_ZEMCpi_PTgam_MISS->Fill(primpt,primzemc,pte,(float)accept1);
		      //}
		      /*
			if(accept2> 0 && fabs(zemcp) < 155){
			
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
			
			
			if(accept1>-9&&(accept2<-8||ecore2<1.0)) {
			ppivspgam_tag->Fill(pe,primp);
			ppivspgam_tag->Fill(pe,primp);
			}
			if(accept2>-9&&(accept1<-8||ecore1<1.0)) {
			ppivspgam_tag->Fill(pp,primp);
			ppivspgam_tag->Fill(pp,primp);
			}
		      */
		      
		      //invmass->Fill(mass);
		      //cout<< "at mass cut with "<< mass <<endl;
		      
		      float dphi = fabs(phi1-phi2);
		      
		      //if(dphi>PI/2) continue;
		      if(ecore2<0.5) continue;
		      if(accept2<=0)continue;
		      
		      if(accept1<= 0) continue;
		      if(fabs(zemc1) > 155) continue;
		    
		      float sume12 = ecore1+ecore2;
		      if (sume12 < 4.0) continue;
		      
		      //float asym12 =  (ecore1 - ecore2)/ sume12; // fabs not neccesary because upper > check
		      //const float mina = 0.15;  //this is the smallest asym value used for the cut


		      //pi_asym->Fill(asym12);
		      
		      //if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25)))continue;		 
		      
		      
		      //cout << "HERE we're gonna fill INVMASS histos" <<endl;
		      invmass->Fill(mass);
		      
		      TLorentzVector pho1, pho2, pi0;
		      pho1.SetPxPyPzE(px1,py1,pz1,ecore1);
		      pho2.SetPxPyPzE(px2,py2,pz2,ecore2);
		      pi0=pho1+pho2;
		      
		      float pi0pt=pi0.Pt();		 
		      
		      //cout << "mass: " << mass << " pi0pt: " << pi0pt << " from " << ecore1 << " and " << ecore2 <<endl;
		      
		      
		      if(ispi0<1 && ecore2>1.0) INVMASS[0]->Fill(pt1,mass);
		      if(ispi0<2 && ecore2>1.5) INVMASS[1]->Fill(pt1,mass);
		      if(ispi0<3 && ecore2>2.0) INVMASS[2]->Fill(pt1,mass);
		      if(ispi0<4 && ecore2>2.5) INVMASS[3]->Fill(pt1,mass);
		      if(ispi0<5 && ecore2>3.0) INVMASS[4]->Fill(pt1,mass);
		      
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
		      
		      if(mass>0.108 && mass<0.165){ 
			//TAGGED!
			//cout << "TAGGED!!!" <<endl;			
			
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



		      //cout << "parent = " << parent1 << endl;
		      /*
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
		      */
		      /*
			if(falsepair<1 && parent2!=parent1)
			{
			FlaseTag->Fill(pipt1);
			falsepair=1;
			}
		      */
		    }else{
		      //Not tagged
		      
		      if(parent1==111) PTpi_ZEMCpi_PTgam_MISS->Fill(pipt1,zemc1,pt1);
		    }


		  /*
		    if(taggedtrig<1)
		    {
		    Ntag->Fill(pt1);
		    taggedtrig=1;
		    }
		  */
		  //}//gampart loop
		}//pt bwtn 5 & 20
	      }//if accept1>0
	    //cout << "end of gamtrig loop" <<endl;
	  }//gamtrig loop
	//cout <<"end tagging loop" <<endl;
      }//tagging


    //Now do pi0 loop!

    if(tagging==1)
      {
	//cout <<"HERE in pi0 tagging loop with " << nparts <<endl;
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
 		iparent1=gamindex[gamtrig];
 		pipt1=gampi0pt[gamtrig];
 

		gamma_acc->Fill(pt1);

		//if(parent1==0 && ecore1<1) continue;
		if(ecore1<1) continue;
		//cout << gamtrig << "Trig is greater than 1 GeV" <<endl;

		gamma_ecore->Fill(pt1);

		int ispi0=0;
		
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

		      

		      //cout << gamtrig << " passed ecore cuts with " << gampart <<endl;


		      float mass= sqrt(2.0*ecore1*ecore2 - 2.0*(px1*px2+py1*py2+pz1*pz2)); 
		      
		      //if(iparent1==iparent2) cout << "match " << iparent1  << "with mass " <<mass << " and ecore2 " << ecore2 << " ecore1 " << ecore1 << " accept1 "<< accept1 << " accept2 " << accept2 <<endl; 


		      //make fiducial cuts  -- 155 cm for higher pt gamma
		      double zemc1= 510.0*pz1/pt1 + vtxz;
		      double zemc2= 510.0*pz2/pt2 + vtxz;
		      
		      float dphi = fabs(phi1-phi2);
		      
		      //if(dphi>PI/2) continue;
		      //if(ecore2<0.5) continue;
		      
		      //if(accept1<= 0) continue;
		  
		      //cout << gamtrig << " passed acceptance with " << gampart <<endl;

		      float sume12 = ecore1+ecore2;
		      if (sume12 < 4.0) gamma_pair[4]->Fill(pi0pt);;

		      if (sume12 < 4.0) continue;

		      

		      if(fabs(zemc1) > 155){
			gamma_pair[5]->Fill(pi0pt); 
			//continue;
		      }
		      //cout << gamtrig << " passed zemc with " << gampart <<endl;
		      
		  


		      //cout << gamtrig << " passed sume12 cuts with " << gampart <<endl;
   

		      float asym12 =  (ecore1 - ecore2)/ sume12; // fabs not neccesary because upper > check
		      const float mina = 0.15;  //this is the smallest asym value used for the cut
		      
		      //if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25))) cout << "asymfailed for ecores " <<  asym12 << " > " << (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25)) <<endl;

		      //if (pipt1==pipt2 && sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25))) cout << "real_asymfail" <<  pipt1 <<endl;

		      pi_asym->Fill(asym12);

		      if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25))) gamma_pair[6]->Fill(pi0pt);
		      if (sume12 < 5.25 && asym12 > (mina + (1.0-mina)*((sume12 - 4.0)*(sume12 - 4.0)/1.25/1.25)))continue;		 
		      
		      //cout << gamtrig << " passed assym cuts with " << gampart <<endl;

		      
		      gamma_pair[7]->Fill(pi0pt);


		      //cout << "HERE we're gonna fill INVMASS histos" <<endl;
		      invmass->Fill(mass);
		      //cout << "mass: " << mass << " pi0pt: " << pi0pt << " from " << ecore1 << " and " << ecore2 <<endl;
		      
		      //if(pipt1==pipt2) 
		      //cout << " iparent1 " << iparent1 << " for ipart " << gamtrig <<" iparent2 " << iparent2 << "for ipart " << gampart <<endl;



		      

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

			if(parent1==221)eta_trueptvsrecopt->Fill(pi0pt,pipt1);
			if(parent1==111){
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
			if(lead_parent==111){
			  if(lead_parent==111 && iparent1==iparent2){
			    realpi0_leadptvsrecopt->Fill(pi0pt,lead_pipt);
			  }else{
			    pi0_leadptvsrecopt->Fill(pi0pt,lead_pipt);
			  }
			}
			if(lead_parent==-111 || lead_parent==22) direct_leadptvsrecopt->Fill(pi0pt,lead_pipt);
		      
		      }
		      
		  }//gampart

	      }//if accept1>0
	    //cout << "end of gamtrig loop" <<endl;
	  }//gamtrig loop
	//cout <<"end tagging loop" <<endl;
      }//tagging



    //cout << "Filled Root Objects!! HERE" <<endl;
    
    return;
}

