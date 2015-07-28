//-----------------------------------------------------------------------------
//
// Generate lepton-pair mass distributions for Dalitz decays according
// to the Kroll-Wada parametrization 
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>

#include <TH1.h>
#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticlePropertyList.h"
#include "DecayList.h"

TH1F * InitializeRandomHist(TH1F *);
double FormFactor(double, ParticleProperty*);

TH1F * CreateDalitzLeptonPairMass(Decay *DDalitz, ParticlePropertyList *PPList)
{
  int ibody, ID;
  int ibin, nbins=1000;
  double min, max, binwidth;
  double pmass, lmass, omass;
  double epsilon, delta, m_ll, q, kw_help, kroll_wada, form_factor, weight;

  ParticleProperty * PParent=0;
  ParticleProperty * PLepton=0;
  ParticleProperty * POther=0;

  PParent = PPList->GetByID(DDalitz->GetParentID());
  for ( ibody=1; ibody<=3; ibody++ )
  {
    ID = DDalitz->GetChildID(ibody);
    if ( abs(ID)==11 || abs(ID)==13 )
      PLepton = PPList->GetByID(ID);
    else
      POther = PPList->GetByID(ID);
  }  
  pmass = PParent->GetMass();
  lmass = PLepton->GetMass();
  omass = POther->GetMass();

  min = 2.0*lmass;
  max = pmass-omass;
  binwidth = (max-min)/(double)nbins;
  TH1F * hdal = new TH1F("hdal","Dalitz",nbins,min,max);

  epsilon = (lmass/pmass)*(lmass/pmass);
  delta   = (omass/pmass)*(omass/pmass);

  for ( ibin=1; ibin<=nbins; ibin++ )
  {
    m_ll = min + (double)(ibin-1)*binwidth+binwidth/2.0;
    q    = (m_ll/pmass)*(m_ll/pmass);
    if ( q<=4.0*epsilon )
    {
      cout << "Error in calculating Dalitz mass histogram binning! q" << endl;
      hdal = 0;
      return hdal;
    }	
    kw_help = (1.0+q/(1.0-delta))*(1.0+q/(1.0-delta))
            - 4.0*q/((1.0-delta)*(1.0-delta));    
    if ( kw_help<=0.0 )
    {
      cout << "Error in calculating Dalitz mass histogram binning! kw" << endl;
      hdal = 0;
      return hdal;
    }	
    kroll_wada = (2.0/m_ll) * exp(1.5*log(kw_help))
                            * sqrt(1.0-4.0*epsilon/q)   
                            * (1.0+2.0*epsilon/q);
    form_factor = FormFactor(m_ll*m_ll, PParent);
    weight = kroll_wada * form_factor;    
    hdal->AddBinContent(ibin,weight);
  }

  hdal = InitializeRandomHist(hdal);

  return hdal;

}





