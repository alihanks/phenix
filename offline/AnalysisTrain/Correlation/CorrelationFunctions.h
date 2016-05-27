// CorrelationFunction.h includes utility functions for Correlation code package
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH3.h>
#include <TH2.h>
#include <PHAngle.h>

#include <APiZero.h>
#include <ACluster.h>
#include <ATrack.h>

const float PI = acos(-1.0);
const float WEST_LOW_EDGE = -3.0*PI/16.0;
const float WEST_HIGH_EDGE = 5.0*PI/16.0;
const float EAST_LOW_EDGE = 11.0*PI/16.0;
const float EAST_HIGH_EDGE = 19*PI/16.0;

inline float CalculateFoldedDphi(float phi1, float phi2)//0-pi
{
  float dphi = PHAngle(phi1-phi2);
  // if ( dphi < 0 )    dphi+=2*PI;
  // if ( dphi > 2*PI ) dphi-=2*PI;
  // if ( dphi > PI )   dphi=2*PI-dphi;

  if(dphi < -PI)             dphi += 2*PI;
  if(dphi < 0 && dphi >= -PI) dphi = fabs(dphi);
  if(dphi > PI)              dphi = 2*PI - dphi;
  return dphi;
}

inline float CalculateDphi(float phi1, float phi2)//-0.5pi - 1.5pi
{
  
  float dphi = PHAngle(phi1-phi2);
  
  if(dphi < -0.5*PI) {/*std::cout <<"in CalculateDphi: phi1 = " << phi1 <<"; phi2 = " << phi2 << std::endl; std::cout <<"dphi = " << dphi <<std::endl; std::cout << "dphi<-0.5PI! +=2PI!" << std::endl; */dphi += 2*PI;}
  if(dphi > 1.5*PI ) {/*std::cout <<"in CalculateDphi: phi1 = " << phi1 <<"; phi2 = " << phi2 << std::endl; std::cout <<"dphi = " << dphi <<std::endl; std::cout << "dphi>1.5PI! -=2PI!" << std::endl;*/ dphi -= 2*PI;}
  if(fabs(dphi+0.5*PI)<0.00001) {/*std::cout << "out of bound: dphi = " << dphi << std::endl; */dphi += 0.00001;}
  if(fabs(dphi-1.5*PI)<0.00001) {/*std::cout << "out of bound: dphi = " << dphi << std::endl; */dphi -= 0.00001;}
  return dphi;
}

template<class T> void SetIso(T* trigger,
                              std::vector<ATrack*> lessqualtrk_vec,
                              std::vector<ACluster*> all_clus_vec,
                              float Rcut = 0.3,
                              float Emin = 0,
                              TH3F* h3_dR = NULL,
                              TH3F* h3_etot = NULL,
                              TH2F* h2_dR = NULL,
                              TH2F* h2_iso = NULL,
                              TH3F* h3_iso_acc = NULL
                              )
{
  float pt = trigger->E();
  float phi = trigger->Phi();
  float eta = trigger->Eta();

  float etot = 0;
  for ( unsigned int iclus = 0; iclus < all_clus_vec.size(); iclus++ )
  {
    float pt1 = all_clus_vec[iclus]->E();
    if ( pt1 < 0.2 ) continue;
    float dR = all_clus_vec[iclus]->DeltaR(*(TLorentzVector*)trigger);
    float dphi = CalculateFoldedDphi(all_clus_vec[iclus]->Phi(), phi);
    // if ( dR < 0.0001 ) continue;

    if ( dR < Rcut ) etot += pt1;

    if ( h3_dR ) h3_dR->Fill(pt, pt1, dR);
    if ( h3_etot ) h3_etot->Fill(pt, pt1/pt, dR);
    if ( h2_dR ) h2_dR->Fill(dR, dphi);
  }

  for ( unsigned int itrk = 0; itrk < lessqualtrk_vec.size(); itrk++ ) {
    float pt1 = lessqualtrk_vec[itrk]->Pt();
    if ( pt1 < 0.2 ) continue;
    if ( pt1 > 7.0 ) continue;
    float dR = lessqualtrk_vec[itrk]->DeltaR(*(TLorentzVector*)trigger);
    float aeta = lessqualtrk_vec[itrk]->Eta();
    float deta = eta - aeta;
    float dphi = CalculateFoldedDphi(lessqualtrk_vec[itrk]->Phi(), phi);
    dR = sqrt(dphi*dphi + deta*deta);

    if ( dR < Rcut ) etot += pt1;

    if ( h3_dR ) h3_dR->Fill(pt, pt1, dR);
    if ( h3_etot ) h3_etot->Fill(pt, pt1/pt, dR);
    if ( h2_dR ) h2_dR->Fill(dR, dphi);
  }
  if ( h2_iso ) h2_iso->Fill(etot, pt);
  float emin_frac = Emin/pt;
  // Trigger energy IS included in total: check against 110% + background energy in % of trigger pt
  if ( etot < 1*(.1+emin_frac)*pt ) {
    trigger->SetIso(true);
    if ( h3_iso_acc ) h3_iso_acc->Fill(phi, eta, pt);
  }
  else trigger->SetIso(false);
}

template<class T> void ClearVector(std::vector<T*> vec)
{
  for ( unsigned int i=0; i < vec.size(); i++ ) delete vec[i];
  vec.clear();
}


