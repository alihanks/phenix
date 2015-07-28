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

#define PI  3.141592653589793238

#define  INCLUDEFLAG extern
#include "DeclareROOTObjects.h"

bool   ParticleIsLepton(Particle*);
double thetaof(Mom3);
double phiof(Mom3);
Mom4   ApplyResolution(int, Mom4);
int    PHENIXFilter(Particle*, ParticlePropertyList*);
bool   PairInAcceptance(int, Mom4, Mom4);

void FillROOTObjects(int setup, double dNdy_pi0, double N_coll,
		     char *output_file, 
		     ParticleList *PList, 
		     ParticlePropertyList *PPList)
{
  double weight, fillweight, fillweight_inv;
  double binwidthcorr, binwidthcorrT, convprob;
  double n_pi0_generated;
  double sigma_inelastic, abs_norm;
  Mom4   pold, mom, mom1, mom2;
  Mom4   p_mom1, p_mom2;
  double sece, secp;
  double mass, E, px, py, pz, pt, theta, phi, y;
  double minv, pxe, pye, pze, pxp, pyp, pzp;
  double pe, pp, thp, php, p;
  double th = 0.0;
  double ph = 0.0;
  double pt_parent = 0.0;
  int    charge, pide;
  int 	 pid = 0;

  n_pi0_generated = 0.1e+06;
  sigma_inelastic = 42.2;
  abs_norm        = dNdy_pi0; // for AuAu: multiplicity
  // abs_norm        = dNdy_pi0 * sigma_inelastic; // for pp: cross section

  PLNode   * CurrentNode     = PList->GetHeadNode();
  PLNode   * NextNode        = 0;
  Particle * CurrentParticle = 0;
  Particle * NextParticle    = 0;
  Particle * CurrentPrimary  = 0;

  CurrentNode     = CurrentNode->GetNextNode();
  CurrentParticle = CurrentNode->Get(0);

  if ( PList->GetLength()>1 )
  {
    NextNode        = CurrentNode->GetNextNode();
    NextParticle    = NextNode->Get(0);
  }
  CurrentPrimary  = CurrentParticle;

  for ( int i=1; i<=PList->GetLength(); i++) 
  {
    mom1   = CurrentParticle->Get4mom();
    weight = CurrentParticle->GetWeight();
    mass   = sqrt(mom1*mom1);
    E      = mom1.GetE();
    p      = sqrt(E*E-mass*mass);
    px     = mom1.Getp().Getpx();
    py     = mom1.Getp().Getpy();
    pz     = mom1.Getp().Getpz();
    pt     = sqrt(px*px+py*py);
    y      = log((E+pz)/(E-pz))/2.0;
    theta  = thetaof(mom1.Getp());
    phi    = phiof(mom1.Getp());

    if ( CurrentParticle->GetGeneration()==1 )
    {
      CurrentPrimary = CurrentParticle;
      pid            = CurrentParticle->GetID();
      // primaries->Fill(mass,weight,p,theta,phi,pid);
    }

    bool final_state = false;
    if ( i != PList->GetLength() ) {
      if ( (CurrentParticle->GetGeneration()>=NextParticle->GetGeneration()) ) 
	final_state = true;
    }
    else {
      final_state = true;
    }

    if ( final_state ) {
      // parent:
      mom    = CurrentPrimary->Get4mom();
      pid    = CurrentPrimary->GetID();
      p      = sqrt(mom.GetE()*mom.GetE()-mom*mom);
      th     = thetaof(mom.Getp());
      ph     = phiof(mom.Getp());
      pt_parent = p*sin(thetaof(mom.Getp()));
      // particle:
      weight = CurrentParticle->GetWeight();
      mom1   = CurrentParticle->Get4mom();
      //      mom1   = ApplyResolution(setup,mom1);
      pide   = CurrentParticle->GetID();
      pp     = sqrt(mom1.GetE()*mom1.GetE()-mom1*mom1);
      thp    = thetaof(mom1.Getp());
      php    = phiof(mom1.Getp());
      mass   = PPList->GetByID(CurrentParticle->GetID())->GetMass();
      charge = PPList->GetByID(CurrentParticle->GetID())->GetCharge();
      secp   = PHENIXFilter(CurrentParticle,PPList);
      CurrentParticle->SetAccept((int)secp);
      // if ( secp>=-10. ) 
      //   singles->Fill(weight,charge,secp,pp,thp,php,pide,p,th,ph,pid);

      if ( abs(pide)==11 || pide==21 ) {
      // if ( pide==22 ) {
        fillweight = abs_norm*weight/(0.5*n_pi0_generated);
        fillweight = fillweight/2.0;
        fillweight = fillweight/0.1;
        fillweight = fillweight/(2.0*PI);
	binwidthcorr  = 0.1/pteR->GetBinWidth(pteR->FindBin(pt));
	binwidthcorrT = 0.1/pteRT->GetBinWidth(pteRT->FindBin(pt));

        fillweight_inv = fillweight/pt;

	// AuAu Run-2
	// convprob = 1.573; // first guess from photonic spectra
	// convprob = 1.17; // update from pizero simulation
	convprob = 1.25; // update from pizero simulation, second try

	// AuAu Run-4
	// convprob = 0.65;

	// pp
	// convprob = 0.730;

	if ( convprob<0.0 ) convprob = 0.0; 

        if ( fabs(y) <= 0.5 )
        {
          if (pid==111) 
	  {
	    fillweight     = fillweight     * (1.+0.636/exp(7.74*pt_parent));
	    fillweight_inv = fillweight_inv * (1.+0.636/exp(7.74*pt_parent));
	    ptePion->Fill(pt,fillweight_inv);
	    pteCPion->Fill(pt,convprob*fillweight_inv);
            pte->Fill(pt,fillweight_inv);
	    pte->Fill(pt,convprob*fillweight_inv);
	    pteC->Fill(pt,convprob*fillweight_inv);
	    pte2Pion->Fill(pt,fillweight);
	    pte2CPion->Fill(pt,convprob*fillweight);
            pte2->Fill(pt,fillweight);
	    pte2->Fill(pt,convprob*fillweight);
	    pte2C->Fill(pt,convprob*fillweight);
	    pteRPion->Fill(pt,binwidthcorr*fillweight_inv);
	    pteRCPion->Fill(pt,binwidthcorr*convprob*fillweight_inv);
            pteR->Fill(pt,binwidthcorr*fillweight_inv);
	    pteR->Fill(pt,binwidthcorr*convprob*fillweight_inv);
	    pteRC->Fill(pt,binwidthcorr*convprob*fillweight_inv);
	    pteR2Pion->Fill(pt,binwidthcorr*fillweight);
	    pteR2CPion->Fill(pt,binwidthcorr*convprob*fillweight);
            pteR2->Fill(pt,binwidthcorr*fillweight);
	    pteR2->Fill(pt,binwidthcorr*convprob*fillweight);
	    pteR2C->Fill(pt,binwidthcorr*convprob*fillweight);
	    pteRTPion->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRTCPion->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
            pteRT->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRT->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
	    pteRTC->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
	    pteRT2Pion->Fill(pt,binwidthcorrT*fillweight);
	    pteRT2CPion->Fill(pt,binwidthcorrT*convprob*fillweight);
            pteRT2->Fill(pt,binwidthcorrT*fillweight);
	    pteRT2->Fill(pt,binwidthcorrT*convprob*fillweight);
	    pteRT2C->Fill(pt,binwidthcorrT*convprob*fillweight);
	  }
          if (pid==221)
	  {
	    convprob = convprob/1.255;
	    pteEta->Fill(pt,fillweight_inv);
	    pteCEta->Fill(pt,convprob*fillweight_inv);
            pte->Fill(pt,fillweight_inv);
	    pte->Fill(pt,convprob*fillweight_inv);
	    pteC->Fill(pt,convprob*fillweight_inv);
	    pte2Eta->Fill(pt,fillweight);
	    pte2CEta->Fill(pt,convprob*fillweight);
            pte2->Fill(pt,fillweight);
	    pte2->Fill(pt,convprob*fillweight);
	    pte2C->Fill(pt,convprob*fillweight);
	    pteREta->Fill(pt,binwidthcorr*fillweight_inv);
	    pteRCEta->Fill(pt,binwidthcorr*convprob*fillweight_inv);
            pteR->Fill(pt,binwidthcorr*fillweight_inv);
	    pteR->Fill(pt,binwidthcorr*convprob*fillweight_inv);
	    pteRC->Fill(pt,binwidthcorr*convprob*fillweight_inv);
	    pteR2Eta->Fill(pt,binwidthcorr*fillweight);
	    pteR2CEta->Fill(pt,binwidthcorr*convprob*fillweight);
            pteR2->Fill(pt,binwidthcorr*fillweight);
	    pteR2->Fill(pt,binwidthcorr*convprob*fillweight);
	    pteR2C->Fill(pt,binwidthcorr*convprob*fillweight);
	    pteRTEta->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRTCEta->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
            pteRT->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRT->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
	    pteRTC->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
	    pteRT2Eta->Fill(pt,binwidthcorrT*fillweight);
	    pteRT2CEta->Fill(pt,binwidthcorrT*convprob*fillweight);
            pteRT2->Fill(pt,binwidthcorrT*fillweight);
	    pteRT2->Fill(pt,binwidthcorrT*convprob*fillweight);
	    pteRT2C->Fill(pt,binwidthcorrT*convprob*fillweight);
	  }
          if (pid==331)
	  {
	    convprob = convprob/3.50;
	    pteEtaprime->Fill(pt,fillweight_inv);
	    pteCEtaprime->Fill(pt,convprob*fillweight_inv);
            pte->Fill(pt,fillweight_inv);
	    pte->Fill(pt,convprob*fillweight_inv);
	    pteC->Fill(pt,convprob*fillweight_inv);
	    pte2Etaprime->Fill(pt,fillweight);
	    pte2CEtaprime->Fill(pt,convprob*fillweight);
            pte2->Fill(pt,fillweight);
	    pte2->Fill(pt,convprob*fillweight);
	    pte2C->Fill(pt,convprob*fillweight);
	    pteREtaprime->Fill(pt,binwidthcorr*fillweight_inv);
	    pteRCEtaprime->Fill(pt,binwidthcorr*convprob*fillweight_inv);
            pteR->Fill(pt,binwidthcorr*fillweight_inv);
	    pteR->Fill(pt,binwidthcorr*convprob*fillweight_inv);
	    pteRC->Fill(pt,binwidthcorr*convprob*fillweight_inv);
	    pteR2Etaprime->Fill(pt,binwidthcorr*fillweight);
	    pteR2CEtaprime->Fill(pt,binwidthcorr*convprob*fillweight);
            pteR2->Fill(pt,binwidthcorr*fillweight);
	    pteR2->Fill(pt,binwidthcorr*convprob*fillweight);
	    pteR2C->Fill(pt,binwidthcorr*convprob*fillweight);
	    pteRTEtaprime->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRTCEtaprime->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
            pteRT->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRT->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
	    pteRTC->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
	    pteRT2Etaprime->Fill(pt,binwidthcorrT*fillweight);
	    pteRT2CEtaprime->Fill(pt,binwidthcorrT*convprob*fillweight);
            pteRT2->Fill(pt,binwidthcorrT*fillweight);
	    pteRT2->Fill(pt,binwidthcorrT*convprob*fillweight);
	    pteRT2C->Fill(pt,binwidthcorrT*convprob*fillweight);
	  }
          if (pid==113)
	  { 
	    pteRho->Fill(pt,fillweight_inv);
            pte->Fill(pt,fillweight_inv);
	    pte2Rho->Fill(pt,fillweight);
            pte2->Fill(pt,fillweight);
	    pteRRho->Fill(pt,binwidthcorr*fillweight_inv);
            pteR->Fill(pt,binwidthcorr*fillweight_inv);
	    pteR2Rho->Fill(pt,binwidthcorr*fillweight);
            pteR2->Fill(pt,binwidthcorr*fillweight);
	    pteRTRho->Fill(pt,binwidthcorrT*fillweight_inv);
            pteRT->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRT2Rho->Fill(pt,binwidthcorrT*fillweight);
            pteRT2->Fill(pt,binwidthcorrT*fillweight);
	  }
          if (pid==223)
	  {
	    pteOmega->Fill(pt,fillweight_inv);
            pte->Fill(pt,fillweight_inv);
	    pte2Omega->Fill(pt,fillweight);
            pte2->Fill(pt,fillweight);
	    pteROmega->Fill(pt,binwidthcorr*fillweight_inv);
            pteR->Fill(pt,binwidthcorr*fillweight_inv);
	    pteR2Omega->Fill(pt,binwidthcorr*fillweight);
            pteR2->Fill(pt,binwidthcorr*fillweight);
	    pteRTOmega->Fill(pt,binwidthcorrT*fillweight_inv);
            pteRT->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRT2Omega->Fill(pt,binwidthcorrT*fillweight);
            pteRT2->Fill(pt,binwidthcorrT*fillweight);
	  }
          if (pid==333)
	  {
            ptePhi->Fill(pt,fillweight_inv);
            pte->Fill(pt,fillweight_inv);
            pte2Phi->Fill(pt,fillweight);
            pte2->Fill(pt,fillweight);
            pteRPhi->Fill(pt,binwidthcorr*fillweight_inv);
            pteR->Fill(pt,binwidthcorr*fillweight_inv);
            pteR2Phi->Fill(pt,binwidthcorr*fillweight);
            pteR2->Fill(pt,binwidthcorr*fillweight);
            pteRTPhi->Fill(pt,binwidthcorrT*fillweight_inv);
            pteRT->Fill(pt,binwidthcorrT*fillweight_inv);
            pteRT2Phi->Fill(pt,binwidthcorrT*fillweight);
            pteRT2->Fill(pt,binwidthcorrT*fillweight);
	  }
          if (pid==-111) 
	  {
	    fillweight     = fillweight*N_coll/dNdy_pi0;
	    fillweight_inv = fillweight_inv*N_coll/dNdy_pi0;
	    pteGamma->Fill(pt,1.255*fillweight_inv);
	    pteCGamma->Fill(pt,convprob*fillweight_inv);
            pte->Fill(pt,1.255*fillweight_inv);
	    pte->Fill(pt,convprob*fillweight_inv);
	    pteC->Fill(pt,convprob*fillweight_inv);
	    pte2Gamma->Fill(pt,1.255*fillweight);
	    pte2CGamma->Fill(pt,convprob*fillweight);
            pte2->Fill(pt,1.255*fillweight);
	    pte2->Fill(pt,convprob*fillweight);
	    pte2C->Fill(pt,convprob*fillweight);
	    pteRGamma->Fill(pt,binwidthcorr*1.255*fillweight_inv);
	    pteRCGamma->Fill(pt,binwidthcorr*convprob*fillweight_inv);
            pteR->Fill(pt,binwidthcorr*1.255*fillweight_inv);
	    pteR->Fill(pt,binwidthcorr*convprob*fillweight_inv);
	    pteRC->Fill(pt,binwidthcorr*convprob*fillweight_inv);
	    pteR2Gamma->Fill(pt,binwidthcorr*1.255*fillweight);
	    pteR2CGamma->Fill(pt,binwidthcorr*convprob*fillweight);
            pteR2->Fill(pt,binwidthcorr*1.255*fillweight);
	    pteR2->Fill(pt,binwidthcorr*convprob*fillweight);
	    pteR2C->Fill(pt,binwidthcorr*convprob*fillweight);
	    pteRTGamma->Fill(pt,binwidthcorrT*1.255*fillweight_inv);
	    pteRTCGamma->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
            pteRT->Fill(pt,binwidthcorrT*1.255*fillweight_inv);
	    pteRT->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
	    pteRTC->Fill(pt,binwidthcorrT*convprob*fillweight_inv);
	    pteRT2Gamma->Fill(pt,binwidthcorrT*1.255*fillweight);
	    pteRT2CGamma->Fill(pt,binwidthcorrT*convprob*fillweight);
            pteRT2->Fill(pt,binwidthcorrT*1.255*fillweight);
	    pteRT2->Fill(pt,binwidthcorrT*convprob*fillweight);
	    pteRT2C->Fill(pt,binwidthcorrT*convprob*fillweight);
	  }
          if (pid==21) 
	  {
	    pteKe3->Fill(pt,fillweight_inv);
            pte->Fill(pt,fillweight_inv);
	    pte2Ke3->Fill(pt,fillweight);
            pte2->Fill(pt,fillweight);
	    pteRKe3->Fill(pt,binwidthcorr*fillweight_inv);
            pteR->Fill(pt,binwidthcorr*fillweight_inv);
	    pteR2Ke3->Fill(pt,binwidthcorr*fillweight);
            pteR2->Fill(pt,binwidthcorr*fillweight);
	    pteRTKe3->Fill(pt,binwidthcorrT*fillweight_inv);
            pteRT->Fill(pt,binwidthcorrT*fillweight_inv);
	    pteRT2Ke3->Fill(pt,binwidthcorrT*fillweight);
            pteRT2->Fill(pt,binwidthcorrT*fillweight);
	  }
        }
      }

      if ( pide==22 ) {
        fillweight = abs_norm*weight/(0.5*n_pi0_generated);
        fillweight = fillweight/0.1;
        fillweight = fillweight/(2.0*PI);
        fillweight_inv = fillweight/pt;

        if ( fabs(y) <= 0.5 )
        {
          if (pid==111) 
	  {
	    ptg->Fill(pt,fillweight_inv);
	    ptgPion->Fill(pt,fillweight_inv);
	  }
          if (pid==221) 
	  {
	    ptg->Fill(pt,fillweight_inv);
	    ptgEta->Fill(pt,fillweight_inv);
	  }
          if (pid==331) 
	  {
	    ptg->Fill(pt,fillweight_inv);
	    ptgEtaprime->Fill(pt,fillweight_inv);
	  }
          if (pid==223) 
	  {
	    ptg->Fill(pt,fillweight_inv);
	    ptgOmega->Fill(pt,fillweight_inv);
	  }
          if (pid==333) 
	  {
	    ptg->Fill(pt,fillweight_inv);
	    ptgPhi->Fill(pt,fillweight_inv);
	  }
          if (pid==-111) 
	  {
	    fillweight     = fillweight*N_coll/dNdy_pi0;
	    fillweight_inv = fillweight_inv*N_coll/dNdy_pi0;
	    ptg->Fill(pt,fillweight_inv);
	    ptgGamma->Fill(pt,fillweight_inv);
	  }
	}
      }
    }

    if ( i==PList->GetLength() ) break;

    if ( ParticleIsLepton(CurrentParticle) &&
         ParticleIsLepton(NextParticle) &&
         CurrentParticle->GetGeneration()==NextParticle->GetGeneration() &&
      	 PPList->GetByID(CurrentParticle->GetID())->GetCharge() !=
      	 PPList->GetByID(NextParticle->GetID())->GetCharge() ) 
    {
      p_mom1 = mom1;
      p_mom2 = NextParticle->Get4mom();
      p_mom1 = ApplyResolution(setup,p_mom1);
      p_mom2 = ApplyResolution(setup,p_mom2);

      if ( PPList->GetByID(CurrentParticle->GetID())->GetCharge()==1 )
      {
	pe  = sqrt(p_mom2.GetE()*p_mom2.GetE()-p_mom2*p_mom2);
        pxe = p_mom2.Getp().Getpx();
        pye = p_mom2.Getp().Getpy();
        pze = p_mom2.Getp().Getpz();
	pp  = sqrt(p_mom1.GetE()*p_mom1.GetE()-p_mom1*p_mom1);
        pxp = p_mom1.Getp().Getpx();
        pyp = p_mom1.Getp().Getpy();
        pzp = p_mom1.Getp().Getpz();
      }
      else
      {
	pe  = sqrt(p_mom1.GetE()*p_mom1.GetE()-p_mom1*p_mom1);
        pxe = p_mom1.Getp().Getpx();
        pye = p_mom1.Getp().Getpy();
        pze = p_mom1.Getp().Getpz();
 	pp  = sqrt(p_mom2.GetE()*p_mom2.GetE()-p_mom2*p_mom2);
        pxp = p_mom2.Getp().Getpx();
        pyp = p_mom2.Getp().Getpy();
        pzp = p_mom2.Getp().Getpz();
      }
      E  = p_mom1.GetE() + p_mom2.GetE();
      px = pxe + pxp;
      py = pye + pyp;
      pz = pze + pzp;
      pt = sqrt(px*px+py*py);
      y  = log((E+pz)/(E-pz))/2.0;
      
      // sece = (double)ElectronInPHENIXAcceptance(pxe,pye,pze,0.0,-1.,
      //                                          0.15,35.);
      // secp = (double)ElectronInPHENIXAcceptance(pxp,pyp,pzp,0.0,1.,
      //                                          0.15,35.);
      sece = 0.;
      secp = 0.;

      if ( secp>=-10. && sece>=-10. && fabs(y)<=0.5 ) 
      {
        minv = sqrt((p_mom1+p_mom2)*(p_mom1+p_mom2));

        weight     = CurrentParticle->GetWeight();
        fillweight = abs_norm*weight/(0.5*n_pi0_generated);
        fillweight = fillweight/0.025;

	if ( pid==111 )
	{
	  mee->Fill(minv,fillweight);
	  meePion->Fill(minv,fillweight);
	}
	if ( pid==221 )
	{
	  mee->Fill(minv,fillweight);
	  meeEta->Fill(minv,fillweight);
	}
	if ( pid==331 )
	{
	  mee->Fill(minv,fillweight);
	  meeEtaprime->Fill(minv,fillweight);
	}
	if ( pid==113 )
	{
	  mee->Fill(minv,fillweight);
	  meeRho->Fill(minv,fillweight);
	}
	if ( pid==223 )
	{
	  mee->Fill(minv,fillweight);
	  meeOmega->Fill(minv,fillweight);
	}
	if ( pid==333 )
	{
	  mee->Fill(minv,fillweight);
	  meePhi->Fill(minv,fillweight);
	}
	if ( pid==443 )
	{
	  mee->Fill(minv,fillweight);
	  meeJPsi->Fill(minv,fillweight);
	}
      }
    }

    CurrentNode     = NextNode;
    CurrentParticle = NextParticle;
    if ( i<PList->GetLength()-1 )
    {
      NextNode = CurrentNode->GetNextNode();
      NextParticle    = NextNode->Get(0);
    }  

  }

  return;

}

