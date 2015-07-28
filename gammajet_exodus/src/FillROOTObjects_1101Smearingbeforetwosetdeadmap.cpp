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
Mom4   ApplyEnergyResolution(Mom4);
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
    
  char deadfile[]="/direct/phenix+workarea/hgong/gammajet.exodus/deadmap.lst";
  char emcalfile[]="/direct/phenix+workarea/hgong/gammajet.exodus/deadmap.lst.isobe";

  double weight, fillweight, fillweight1, fillweight2;
  double n_pi0_generated;
  double sigma_inelastic, abs_norm;

  double ptcut = 0.2;
  double type,type1,type2,per=0.16;
  double sepacut=12.0/500.0;
  double opangle;

  // storage for pi0 pair
  Mom4 pairmom[2];
  double pairweight[2];
  int pairaccept[2];
  int pi0pair=-1;
  double pairtype[2];

  //for decay photons
  Mom4   mom1, mom2;
  int accept, accept1, accept2;
  double pxe, pye, pze, pxp, pyp, pzp;
  double pte, ptp, pe, pp;
  double vtxz;
  double theta;
  double ye,yp;

  //for prime particle
  int 	 pid = 0;
  Mom4 primmom1;
  double primweight, primfillweight, primpt;

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
    CurrentNode = CurrentNode->GetNextNode();
    CurrentParticle = CurrentNode->Get(0);
    //   cout << CurrentParticle->GetID() << ":" << CurrentParticle->GetGeneration() << endl;

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

	if (pid == 111)
	  realpi0pt->Fill(primpt, primfillweight);
      
	if (pid == 221)
	  realetapt->Fill(primpt, primfillweight);
	
	if (pid == 331)
	  realetaprimept->Fill(primpt, primfillweight);
	
	if (pid == 223 )
	  realomegapt->Fill(primpt, primfillweight);
	
	//      cout << "found a primary: " << pid << " wt " << endl;
      }
    
    if ( CurrentParticle->GetID()==22 )
      //	 || CurrentParticle->GetID()== -111  ) direct photon
      {
	// cout << "found a photon" << endl;
	// correct weighting for pairs (if decay chains are full implemented!) 
	weight   = CurrentParticle->GetWeight();
	fillweight = abs_norm*weight/(0.5*n_pi0_generated);
	fillweight = fillweight/0.025;
	// photon 4-vector
	mom1 = CurrentParticle->Get4mom();
	vtxz = CurrentParticle->GetzVertex();
	type = gRandom->Rndm();
	// resolution
	mom1 = ApplyEnergyResolution(mom1);
	// acceptance filter
	//accept1 = true;
	//	accept1 = PHENIXPhotonFilter(CurrentParticle);
	//	bool accept1_fid = PHENIXPhotonFilter(CurrentParticle, 1);
      
	pe  = sqrt(mom1.GetE()*mom1.GetE()-mom1*mom1);
	pxe = mom1.Getp().Getpx();
	pye = mom1.Getp().Getpy();
	pze = mom1.Getp().Getpz();
	pte = sqrt(pxe*pxe+pye*pye);

	cocktail->Fill(pte,1.0/pte);

	accept = -1;
	for(int j=0;j<6;j++)
	    if(phenix_geom(j,(float)pxe,(float)pye,(float)pze,(float)vtxz,deadfile)) accept=j;

	if(type<per&&pte>ptcut){
	  gammapt->Fill(pte, fillweight);
	  if(accept!=-1)
	    gamma_acc->Fill(pte,fillweight);
	}

	if(pid==111)
	  {
	    if(type<per&&pte>ptcut)
	      pi0gamma->Fill(pte,fillweight);

	    //store pi0 pair information for pair analysis
	    pi0pair++;
	    pairmom[pi0pair]=mom1;
	    pairaccept[pi0pair]=accept;
	    pairweight[pi0pair]=fillweight;
	    pairtype[pi0pair]=type;
	  }
	
	if(pid==331)
	  etaprimegamma->Fill(pte,fillweight);
	if(pid==221)
	  etagamma->Fill(pte,fillweight);
	if(pid==223)
	  omegagamma->Fill(pte,fillweight);

      }

    if(pi0pair==1)
      {
	pxe = pairmom[0].Getp().Getpx();
	pye = pairmom[0].Getp().Getpy();
	pze = pairmom[0].Getp().Getpz();
	pte = sqrt(pxe*pxe+pye*pye);
	
	pxp = pairmom[1].Getp().Getpx();
	pyp = pairmom[1].Getp().Getpy();
	pzp = pairmom[1].Getp().Getpz();
	ptp = sqrt(pxp*pxp+pyp*pyp);
	
	fillweight1=pairweight[0];
	fillweight2=pairweight[1];
	  
	accept1=pairaccept[0];
	accept2=pairaccept[1];
	
	type1=pairtype[0];
	type2=pairtype[1];
	
	//cout << pxe << pye << pze << endl;
	//	  cout << NextParticle->GetID() << "Accept1:" << accept1 << "Vtxz:" << vtxz << endl;
	
	//openning angle of two pion-photons
	opangle = (pairmom[0].Getp()*pairmom[1].Getp()) /
	  (sqrt(pairmom[0].Getp()*pairmom[0].Getp())*
	   sqrt(pairmom[1].Getp()*pairmom[1].Getp()));
	opangle=acos(opangle);

	openangle->Fill(opangle,primfillweight);
	
	//	cout << pxp << pyp << pzp << endl;
	//	  cout << NextParticle->GetID() << "Accept2:" << accept2 << "Vtxz:" << vtxz << endl;
	
	
	if(fillweight1==0 || fillweight2==0)
	  cout << "Something is wrong!" << endl;
	
	//	  pairs->Fill(pte,accept1,fillweight1,ptp,accept2,fillweight2);
	
	if(opangle>sepacut)
	  {

	    //after 10/19/2006 discussion with T.T., f factor without ptcut

	    if(accept1!=-1&&type1<=per&&pte>ptcut)
	      {
		gammaprob[9]->Fill(pte,fillweight1);
		
		if(accept2!=-1)
		//change on 10/25 to make it consistent with Torsten
		//if(accept2!=-1&&ptp>0.2)
		  prob[9]->Fill(pte,fillweight1);
	      }
	    
	    if(accept2!=-1&&type2<=per&&ptp>ptcut)
	      {
		gammaprob[9]->Fill(ptp,fillweight2);
		
		if(accept1!=-1)
		//change on 10/25 make it consistent with Torsten
		//if(accept1!=-1&&pte>ptcut)
		  prob[9]->Fill(ptp,fillweight2);
	      }

	    if(accept1!=-1&&type1<=per&&pte>ptcut)
	      {
		gammaprob[(int)accept1]->Fill(pte,fillweight1);
		if(((int)accept2==(int)accept1)&&ptp>ptcut)
		  prob[(int)accept1]->Fill(pte,fillweight1);
		
		gammaprob[6]->Fill(pte,fillweight1);
		if(accept2!=-1&&ptp>ptcut)
		  prob[6]->Fill(pte,fillweight1);
		
		if(accept1<=3)
		  {
		    gammaprob[7]->Fill(pte,fillweight1);
		    if(accept2<=3&&ptp>ptcut)
		      prob[7]->Fill(pte,fillweight1);
		    }
		else
		  {
		    gammaprob[8]->Fill(pte,fillweight1);
		    if(accept2>3&&ptp>ptcut)
		      prob[8]->Fill(pte,fillweight1);
		  }
	      }
	      
	    if(accept2!=-1&&type2<=per&&ptp>ptcut)
	      {
		gammaprob[6]->Fill(ptp,fillweight2);
		if(accept1!=-1&&pte>ptcut)
		  prob[6]->Fill(ptp,fillweight2);
		
		gammaprob[(int)accept2]->Fill(ptp,fillweight2);
		if(((int)accept1==(int)accept2)&&pte>ptcut)
		  prob[(int)accept2]->Fill(ptp,fillweight2);
		
		if(accept2<=3)
		  {
		    gammaprob[7]->Fill(ptp,fillweight2);
		    if(accept1<=3&&pte>ptcut)
		      prob[7]->Fill(ptp,fillweight2);
		  }
		else
		  {
		    gammaprob[8]->Fill(ptp,fillweight2);
		    if(accept1>3&&pte>ptcut)
		      prob[8]->Fill(ptp,fillweight2);
		  }
	      }
	  }
	
	beforecut->Fill(primpt,primfillweight);
	if(accept1!=-1&&accept2!=-1)
	  accbeforecut->Fill(primpt,primfillweight);

	
	  //for weighting study
	if(pte>ptcut&&ptp>ptcut)
	  {
	    //  if(accept1!=-1&&accept2!=-1)
		{
		  if(type1<=per)
		    {
		      tightptvspt->Fill(pte,ptp,fillweight1);
		      loosein->Fill(ptp,fillweight1);
		    }
		  if(type2<=per)
		    {
		      tightptvspt->Fill(ptp,pte,fillweight1);
		      loosein->Fill(pte,fillweight1);
		    }
		}
	      
	      if(accept1!=-1)
		{
		  if(type1<=per)
		    looseall->Fill(ptp,fillweight1);
		  if(type2<=per)
		    looseall->Fill(pte,fillweight1);
		}
	  }
	
	//  if(pte>ptcut&&ptp>ptcut)
	//if(type1<=per)
	//change on 10/27/2006
	if(accept1!=-1&&accept2!=-1)
	  {
	    if(type1<=per&&pte>ptcut)
	      tightptvspt_wo->Fill(pte,ptp,fillweight1);
	    if(type2<=per&&ptp>ptcut)
	      tightptvspt_wo->Fill(ptp,pte,fillweight1);
	  }
	
	
	//add on 5/3/2006 for double check of formula
	
	if(accept1!=-1&&accept2!=-1&&pte>ptcut&&ptp>ptcut)
	  {
	    if(type1<per)
	      pi0gamma_acc->Fill(pte, fillweight1);
	    if(type2<per)
	      pi0gamma_acc->Fill(ptp, fillweight2);
	    //add on 5/18/2006 for separation cut study
	    if(opangle>sepacut)
	      {
		if(type1<per)
		  pi0gamma_acc_cut->Fill(pte, fillweight1);
		if(type2<per)
		  pi0gamma_acc_cut->Fill(ptp, fillweight2);
	      }
	  }
	
	if(opangle>sepacut)
	  {
	    aftercut->Fill(primpt,primfillweight);
	    if(accept1!=-1&&accept2!=-1)
	      accaftercut->Fill(primpt,primfillweight);
	  }
	
	//pi0pair process
	pi0pair=-1;
      }
  }

  return;

}

