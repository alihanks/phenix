//-----------------------------------------------------------------------------
//
//  Generate ROOT output file
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>

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
bool   PHENIXPhotonFilter(Particle*);
bool   PairInAcceptance(int, Mom4, Mom4);

void FillROOTObjects(int setup, double dNdy_pi0, double N_coll,
		     char *output_file, 
		     ParticleList *PList, 
		     ParticlePropertyList *PPList)
{
  double weight, fillweight;
  double n_pi0_generated;
  double sigma_inelastic, abs_norm;
  Mom4   pold, mom, mom1, mom2;
  Mom4   p_mom1, p_mom2;
  bool   accept1, accept2;
  double E, px, py, pz, pt, y;
  double minv, pxe, pye, pze, pxp, pyp, pzp;
  double pe, pp;
  int 	 pid = 0;

  Mom4 primmom1;
  double primweight, primfillweight, primpt;

  n_pi0_generated = 1.0e+5;
  sigma_inelastic = 42.2;
  abs_norm        = dNdy_pi0; // for AuAu: multiplicity
  // abs_norm        = dNdy_pi0 * sigma_inelastic; // for pp: cross section

  PLNode   * CurrentNode     = PList->GetHeadNode();
  PLNode   * NextNode        = 0;
  Particle * CurrentParticle = 0;
  Particle * NextParticle    = 0;
  Particle * CurrentPrimary  = 0;
  
  // here we start with the diphoton analysis
  static THmulf * gghasym = 0;
  if (gghasym  == 0)
    {
      mgg->GetDirectory()->cd();
      gghasym = new THmulf("gghasym","Pair asymmetry");
      gghasym->AddAxis("mass","Invariant Mass",100,0.,1.0);
      gghasym->AddAxis("pt","Pair Pt",20,0.,4.0);
      //gghasym->AddAxis("mix","Real or mixed",2,-0.5,1.5); 
      gghasym->AddAxis("part","0, pi, 1 eta, 2 etaPrime, 3 rho, 4 omega, 5 phi, ",10,-0.5,9.5);
      //  gghasym->AddAxis("det","detector",1,-0.5,1.5);
      gghasym->AddAxis("asym","Energy Asymmetry Bin",10,0.0,1.0);
      //      se->registerHisto(gghasym->GetName(), gghasym);
      gghasym->Sumw2();
     
    }
  static int  neventK = 0;
  for ( int ipart=1; ipart<=PList->GetLength(); ipart++) 
  {
    CurrentNode = CurrentNode->GetNextNode();
    CurrentParticle = CurrentNode->Get(0);
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
	gghasym->Fill(primfillweight,0.98, primpt, 8. ,0.98);
      if (pid == 221)
	gghasym->Fill(primfillweight,0.98, primpt, 9. ,0.98);      

      //      cout << "found a primary: " << pid << " wt " << endl;
    }
    if ( CurrentParticle->GetID()==22 )
    {
      // cout << "found a photon" << endl;
      // correct weighting for pairs (if decay chains are full implemented!)
      weight     = CurrentParticle->GetWeight();
      fillweight = abs_norm*weight/(0.5*n_pi0_generated);
      fillweight = fillweight/0.025;
      // photon 4-vector
      mom1 = CurrentParticle->Get4mom();
      // resolution
      mom1 = ApplyEnergyResolution(mom1);
      // acceptance filter
      //accept1 = true;
      accept1 = PHENIXPhotonFilter(CurrentParticle);
      NextNode = CurrentNode;
      bool do_pairs = true;
      do
      {
	NextNode     = NextNode->GetNextNode();
	NextParticle = NextNode->Get(0);
	if ( NextParticle->GetID()==22 )
	{
	  // cout << "pair analysis" << endl;
	  // photon 4-vector
	  mom2 = NextParticle->Get4mom();
	  // resolution
	  mom2 = ApplyEnergyResolution(mom2);
	  // acceptance filter
	  // accept2 = true;
	  accept2 = PHENIXPhotonFilter(NextParticle);
	  // single kinematics
	  pe  = sqrt(mom1.GetE()*mom1.GetE()-mom1*mom1);
	  pxe = mom1.Getp().Getpx();
	  pye = mom1.Getp().Getpy();
	  pze = mom1.Getp().Getpz();
	  pp  = sqrt(mom2.GetE()*mom2.GetE()-mom2*mom2);
	  pxp = mom2.Getp().Getpx();
	  pyp = mom2.Getp().Getpy();
	  pzp = mom2.Getp().Getpz();
	  // pair kinematics
	  E    = mom1.GetE() + mom2.GetE();
	  px   = pxe + pxp;
	  py   = pye + pyp;
	  pz   = pze + pzp;
	  pt   = sqrt(px*px+py*py);
	  y    = log((E+pz)/(E-pz))/2.0;	  
	  minv = sqrt((mom1+mom2)*(mom1+mom2));
	  
	  float Easym = ( fabs( mom1.GetE() - mom2.GetE() )/( mom1.GetE() + mom2.GetE() ) );
	  
	  if ( accept1 && accept2 && fabs(y)<=0.5 && pt > 0.4 && pt < 1.0 && Easym > 0.8)
	    //if (true)
	  {
	    // fill histograms
	    if ( pid==111 )
	    {
	      mgg->Fill(minv,fillweight);
	      mggPion->Fill(minv,fillweight);
	      gghasym->Fill(fillweight,minv, pt, 0, Easym);
	    }
	    if ( pid==221 )
	    {
	      mgg->Fill(minv,fillweight);
	      mggEta->Fill(minv,fillweight);
	      gghasym->Fill(fillweight,minv, pt, 1, Easym);
	    }
	    if ( pid==331 )
	    {
	      mgg->Fill(minv,fillweight);
	      mggEtaprime->Fill(minv,fillweight);
	      gghasym->Fill(fillweight,minv, pt, 2, Easym);
	    }
	    if ( pid==113 )
	    {
	      mgg->Fill(minv,fillweight);
	      mggRho->Fill(minv,fillweight);
	      gghasym->Fill(fillweight,minv, pt, 3, Easym);
	    }
	    if ( pid==223 )
	    {
	      mgg->Fill(minv,fillweight);
	      mggOmega->Fill(minv,fillweight);
	      gghasym->Fill(fillweight,minv, pt, 4, Easym);
	    }
	    if ( pid==333 )
	    {
	      mgg->Fill(minv,fillweight);
	      mggPhi->Fill(minv,fillweight);
	      gghasym->Fill(fillweight,minv, pt, 5, Easym);
	    }
	  }
	}
	if ( NextParticle->GetGeneration()==1 ) do_pairs = false;
      } while ( do_pairs );
    }
  }

  return;

}

