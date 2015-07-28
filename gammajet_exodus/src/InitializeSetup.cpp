//-----------------------------------------------------------------------------
//
//  Initialize the ParticleGeneratorList according to the chosen setup
//
//-----------------------------------------------------------------------------

#include "Momentum.h"
#include "Particle.h"
#include "ParticlePropertyList.h"
#include "ParticleGeneratorList.h"

TF2  * CreateTestFunction(int, ParticlePropertyList*);

TH1F * InitializeWhitePt(double, double);
TH1F * InitializeWhiteY(double, double);
TH1F * InitializeExpPt(double, double, double, double);
TH1F * InitializePtCERES(int, ParticlePropertyList*);
TH1F * InitializeYCERES(int, ParticlePropertyList*);
TH1F * InitializePtISR(int, ParticlePropertyList*);
TH1F * InitializeYISR(int, ParticlePropertyList*);
TH1F * InitializePtPHENIX(int, int, double, double, double, double, double,
			  ParticlePropertyList*);
TH1F * InitializeYPHENIX(int, ParticlePropertyList*);
TH1F * InitializeM(int, float, float, ParticlePropertyList*);

ParticleGeneratorList * InitializeSetup(int setup, 
					ParticlePropertyList* PPList, float dNdypi)
{
  ParticleGeneratorList * PGList = new ParticleGeneratorList;
  ParticleGenerator * Pion       = new ParticleGenerator; 
  Pion->SetID(111);
  ParticleGenerator * directPion = new ParticleGenerator; 
  //directPion->SetID(-111);
  directPion->SetID(22);
  ParticleGenerator * Ke3        = new ParticleGenerator; 
  Ke3->SetID(21);
  ParticleGenerator * Eta        = new ParticleGenerator; 
  Eta->SetID(221);
  ParticleGenerator * Etaprime   = new ParticleGenerator; 
  Etaprime->SetID(331);
  ParticleGenerator * Rho        = new ParticleGenerator; 
  Rho->SetID(113);
  ParticleGenerator * Omega      = new ParticleGenerator; 
  Omega->SetID(223);
  ParticleGenerator * Phi        = new ParticleGenerator; 
  Phi->SetID(333);
  ParticleGenerator * Piplus     = new ParticleGenerator; 
  Piplus->SetID(211);
  ParticleGenerator * Piminus    = new ParticleGenerator; 
  Piminus->SetID(-211);
  ParticleGenerator * Kplus      = new ParticleGenerator; 
  Kplus->SetID(321);
  ParticleGenerator * Kminus     = new ParticleGenerator; 
  Kminus->SetID(-321);
  ParticleGenerator * Proton     = new ParticleGenerator; 
  Proton->SetID(2212);
  ParticleGenerator * Antiproton = new ParticleGenerator; 
  Antiproton->SetID(-2212);
  ParticleGenerator * JPsi       = new ParticleGenerator; 
  JPsi->SetID(443);
  ParticleGenerator * EtaC       = new ParticleGenerator; 
  EtaC->SetID(433);
  ParticleGenerator * Upsilon    = new ParticleGenerator; 
  Upsilon->SetID(553);


  TH1F * Histogram = 0;
  TF2  * Function  = 0;

  int    pt_setup  = 1;
  double weight = 1.0;
  double mass, inv_slope;
  double f_c  = 1.;
  double f_p0 = 1.;
  double f_a = 1.;
  double f_b = 1.;
  double f_n = 1.;

  switch(setup)
  {
    case 1:  cout << "Initializing CERES" << endl; 

             Pion->SetWeight(1.0);
             Eta->SetWeight(0.053);
             Etaprime->SetWeight(0.009);
             Rho->SetWeight(0.065);
             Omega->SetWeight(0.065);
             Phi->SetWeight(0.0033);

             Function = CreateTestFunction(Pion->GetID(),PPList);
	     Function = 0;
             Pion->SetPtYFunction(Function);
             Eta->SetPtYFunction(Function);
             Etaprime->SetPtYFunction(Function);
             Rho->SetPtYFunction(Function);
             Omega->SetPtYFunction(Function);
             Phi->SetPtYFunction(Function);

             Histogram = InitializePtCERES(Pion->GetID(),PPList);
             Pion->SetPtHistogram(Histogram);
             Histogram = InitializePtCERES(Eta->GetID(),PPList);
             Eta->SetPtHistogram(Histogram);
             Histogram = InitializePtCERES(Etaprime->GetID(),PPList);
             Etaprime->SetPtHistogram(Histogram);
             Histogram = InitializePtCERES(Rho->GetID(),PPList);
             Rho->SetPtHistogram(Histogram);
             Histogram = InitializePtCERES(Omega->GetID(),PPList);
             Omega->SetPtHistogram(Histogram);
             Histogram = InitializePtCERES(Phi->GetID(),PPList);
             Phi->SetPtHistogram(Histogram);

             Histogram = InitializeYCERES(Pion->GetID(),PPList);
             Pion->SetYHistogram(Histogram);
             Histogram = InitializeYCERES(Eta->GetID(),PPList);
             Eta->SetYHistogram(Histogram);
             Histogram = InitializeYCERES(Etaprime->GetID(),PPList);
             Etaprime->SetYHistogram(Histogram);
             Histogram = InitializeYCERES(Rho->GetID(),PPList);
             Rho->SetYHistogram(Histogram);
             Histogram = InitializeYCERES(Omega->GetID(),PPList);
             Omega->SetYHistogram(Histogram);
             Histogram = InitializeYCERES(Phi->GetID(),PPList);
             Phi->SetYHistogram(Histogram);

             Histogram = 0;
             Pion->SetMHistogram(Histogram);
             Eta->SetMHistogram(Histogram);
             Etaprime->SetMHistogram(Histogram);
	     Histogram = InitializeM(Rho->GetID(),0.,0.,PPList);
             Rho->SetMHistogram(Histogram);
	     Histogram = InitializeM(Omega->GetID(),0.,0.,PPList);
             Omega->SetMHistogram(Histogram);
	     Histogram = InitializeM(Phi->GetID(),0.,0.,PPList);
             Phi->SetMHistogram(Histogram);

	     PGList->Insert(Pion);
	     PGList->Insert(Eta);
	     PGList->Insert(Etaprime);
	     PGList->Insert(Rho);
             PGList->Insert(Omega);
	     PGList->Insert(Phi);

             break;

    case 2:  cout << "Initializing ISR" << endl;

             Pion->SetWeight(1.0);
             Eta->SetWeight(0.053);
             Etaprime->SetWeight(0.009);
             Rho->SetWeight(0.065);
             Omega->SetWeight(0.065);
             Phi->SetWeight(0.0033);

             Histogram = InitializePtISR(Pion->GetID(),PPList);
             Pion->SetPtHistogram(Histogram);
             Histogram = InitializePtISR(Eta->GetID(),PPList);
             Eta->SetPtHistogram(Histogram);
             Histogram = InitializePtISR(Etaprime->GetID(),PPList);
             Etaprime->SetPtHistogram(Histogram);
             Histogram = InitializePtISR(Rho->GetID(),PPList);
             Rho->SetPtHistogram(Histogram);
             Histogram = InitializePtISR(Omega->GetID(),PPList);
             Omega->SetPtHistogram(Histogram);
             Histogram = InitializePtISR(Phi->GetID(),PPList);
             Phi->SetPtHistogram(Histogram);

             Histogram = InitializeYISR(Pion->GetID(),PPList);
             Pion->SetYHistogram(Histogram);
             Histogram = InitializeYISR(Eta->GetID(),PPList);
             Eta->SetYHistogram(Histogram);
             Histogram = InitializeYISR(Etaprime->GetID(),PPList);
             Etaprime->SetYHistogram(Histogram);
             Histogram = InitializeYISR(Rho->GetID(),PPList);
             Rho->SetYHistogram(Histogram);
             Histogram = InitializeYISR(Omega->GetID(),PPList);
             Omega->SetYHistogram(Histogram);
             Histogram = InitializeYISR(Phi->GetID(),PPList);
             Phi->SetYHistogram(Histogram);

             Histogram = 0;
             Pion->SetMHistogram(Histogram);
             Eta->SetMHistogram(Histogram);
             Etaprime->SetMHistogram(Histogram);
	     Histogram = InitializeM(Rho->GetID(),0.,0.,PPList);
             Rho->SetMHistogram(Histogram);
	     Histogram = InitializeM(Omega->GetID(),0.,0.,PPList);
             Omega->SetMHistogram(Histogram);
	     Histogram = InitializeM(Phi->GetID(),0.,0.,PPList);
             Phi->SetMHistogram(Histogram);

	     PGList->Insert(Pion);
	     PGList->Insert(Eta);
	     PGList->Insert(Etaprime);
	     PGList->Insert(Rho);
             PGList->Insert(Omega);
	     PGList->Insert(Phi);

             break;

    case 3:  cout << "Initializing PHENIX" << endl << endl;

             do
             {
	       cout << "Choose the shape of the pt distribution:" << endl;
	       cout << "----------------------------------------" << endl;
	       cout << endl;
	       cout << "1) power law" << endl;
	       cout << "2) exponential (with flow)" << endl;
	       cout << "3) hydrodynamical parameterization" << endl;
	       cout << "4) Run-2 AuAu" << endl;
	       cout << "5) Run-2 pp" << endl;
	       cout << "6) straight power law" << endl;
	       cout << endl;
	       cout << "Your choice (1-6): "; 
	       cin  >> pt_setup; 
	       cout << endl;
	     } while ( pt_setup<1 || pt_setup>6 );


// 	     cin  >> weight;
// 	     cout << endl;
	     Pion->SetWeight(1);
 	     cout << "Pion weight: =1"<<endl;
// 	     cout << "Eta weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
	     // Eta->SetWeight(1);
	     //cout << "Eta weight: 1"<<endl;;
// 	     cout << "Etaprime weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
// 	     Etaprime->SetWeight(weight);
// 	     cout << "Rho weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
// 	     Rho->SetWeight(weight);
// 	     cout << "Omega weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
// 	     Omega->SetWeight(weight);
// 	     cout << "Phi weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
// 	     Phi->SetWeight(weight);
// 	     cout << "J/Psi weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
// 	     JPsi->SetWeight(weight);
// 	     cout << "EtaC weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
// 	     EtaC->SetWeight(weight);
// 	     cout << "Direct Photon weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
 	     //directPion->SetWeight(weight);
	     // 	     directPion->SetWeight(1);
// 	     cout << "Ke3 weight: ";
// 	     cin  >> weight;
// 	     cout << endl;
// 	     Ke3->SetWeight(weight);


	     cout << "f_c: ";
	     cin  >> f_c;
	     cout << endl;
	     cout << "f_p0: ";
	     cin  >> f_p0;
	     cout << endl;
	     cout << "f_a: ";
	     cin  >> f_a;
	     cout << endl;
	     cout << "f_b: ";
	     cin  >> f_b;
	     cout << endl;
	     cout << "f_n: ";
	     cin  >> f_n;
	     cout << endl;
	     
             Histogram = InitializePtPHENIX(pt_setup,Pion->GetID(),
	     f_c, f_p0, f_a, f_b, f_n,
	     PPList);
             Pion->SetPtHistogram(Histogram);
	     /*
             Histogram = InitializePtPHENIX(pt_setup,Eta->GetID(),
	     f_c, f_p0, f_a, f_b, f_n,
	     PPList);
             Eta->SetPtHistogram(Histogram);
	     */
	     /*
             Histogram = InitializePtPHENIX(pt_setup,directPion->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             directPion->SetPtHistogram(Histogram);
	     */     
	     /*
             Histogram = InitializePtPHENIX(pt_setup,Ke3->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Ke3->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Etaprime->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Etaprime->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Rho->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Rho->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Omega->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Omega->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Phi->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Phi->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,JPsi->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             JPsi->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,EtaC->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             EtaC->SetPtHistogram(Histogram);
	     */
	     
             Histogram = InitializeYPHENIX(Pion->GetID(),PPList);
	     Pion->SetYHistogram(Histogram);
	     
	     /*
	     Histogram = InitializeYPHENIX(Eta->GetID(),PPList);
	     Eta->SetYHistogram(Histogram);
	     */
	     /*
             Histogram = InitializeYPHENIX(directPion->GetID(),PPList);
             directPion->SetYHistogram(Histogram);
	     */
	     /*
             Histogram = InitializeYPHENIX(Ke3->GetID(),PPList);
             Ke3->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Etaprime->GetID(),PPList);
             Etaprime->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Rho->GetID(),PPList);
             Rho->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Omega->GetID(),PPList);
             Omega->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Phi->GetID(),PPList);
             Phi->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(JPsi->GetID(),PPList);
             JPsi->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(EtaC->GetID(),PPList);
             EtaC->SetYHistogram(Histogram);
	     */
	     
             Histogram = 0;

	     Pion->SetMHistogram(Histogram);
	     //Eta->SetMHistogram(Histogram);
	     //directPion->SetMHistogram(Histogram);
	     /*
             Ke3->SetMHistogram(Histogram);
             Etaprime->SetMHistogram(Histogram);
	     Histogram = InitializeM(Rho->GetID(),0.,0.,PPList);
             Rho->SetMHistogram(Histogram);
	     Histogram = InitializeM(Omega->GetID(),0.,0.,PPList);
             Omega->SetMHistogram(Histogram);
	     Histogram = InitializeM(Phi->GetID(),0.,0.,PPList);
             Phi->SetMHistogram(Histogram);
	     Histogram = InitializeM(JPsi->GetID(),2.5,3.6,PPList);
             JPsi->SetMHistogram(Histogram);
	     Histogram = InitializeM(EtaC->GetID(),2.4,3.5,PPList);
             EtaC->SetMHistogram(Histogram);
	     */

	     PGList->Insert(Pion);
	     //PGList->Insert(Eta);
	     //PGList->Insert(Etaprime);
	     //	     PGList->Insert(Rho);
	     //PGList->Insert(Omega);
	     //	     PGList->Insert(Phi);
	     //	     PGList->Insert(JPsi);
	     //	     PGList->Insert(EtaC);
	     //	     PGList->Insert(Ke3);
	     //PGList->Insert(directPion);

	     //cout<<" mass of direct photons "<<PPList->GetByID(directPion->GetID())->GetMass()<<endl;
             break;

    case 4:  cout << "Initializing Phi->KK" << endl; 

             Phi->SetWeight(1.);

	     mass = PPList->GetByID(Phi->GetID())->GetMass();
	     inv_slope = 0.240;
             Histogram = InitializeExpPt(0.,4.,mass,inv_slope);
             Phi->SetPtHistogram(Histogram);
             Histogram = InitializeWhiteY(-0.5,0.5);
             Phi->SetYHistogram(Histogram);
	     Histogram = InitializeM(Phi->GetID(),0.,0.,PPList);
             Phi->SetMHistogram(Histogram);

	     PGList->Insert(Phi);

             break;

    case 5:  cout << "Initializing single-particle generator" << endl; 

             break;

    case 6:  cout << "Initializing PHENIX: complete events" << endl << endl;

             do
             {
	       cout << "Choose the shape of the pt distribution:" << endl;
	       cout << "----------------------------------------" << endl;
	       cout << endl;
	       cout << "1) power law" << endl;
	       cout << "2) exponential (with flow)" << endl;
	       cout << "3) hydrodynamical parameterization" << endl;
	       cout << "6) Takao's parameterization" << endl;
	       cout << endl;
	       cout << "Your choice (1-3): "; 
	       cin  >> pt_setup; 
	       cout << endl;
             } while ( pt_setup<1 || pt_setup>6 );
	     /*
             Piplus->SetWeight(0.401);
             Piminus->SetWeight(0.401);
             Kplus->SetWeight(0.062);
             Kminus->SetWeight(0.062);
             Proton->SetWeight(0.039);
             Antiproton->SetWeight(0.036);
	     */
	     //cout << "f_n was " << f_n << " now filled with " << dNdypi <<endl;
	     //f_n=dNdypi;

             Pion->SetWeight(1.0); //changed 1 back to .445 on 1/16/09
	     //weight was originally set to 0.445 then I thought 1
	     //might make sense but for Takao dN/dy is really 2*dN/dy
	     //so I'll leave the weight at this and continue to use
	     //2*dN/dy as my dN_ch/dy input? 
	     //change back to 1 after talking to takao

	     directPion->SetWeight(0.0138);
	     //Found that by multiplying Ndir by 1.2 I get the right Rgamma
	     //Therefore I mulitiplied my 0.0115 estimate from the spectra
	     //by 1.2 to get above 0.0138 --MEC 02/02/09
	     //changed .008 to .018 (=.008/.445) when .445 went to 1.0 
	     //0.018 is too high so I redid the estimate from the direct photon spectrum --> integral*2pi--> 1 direct and 88 pi0s--> 1.15437478150859554e-02
	     

	     
             //Eta->SetWeight(0.043);
	     //UPDATED 1/29/09MEC
             Eta->SetWeight(0.19);//should be .48 but only doing eta->gg 
	     
	     /*
             Etaprime->SetWeight(0.0038);
             Rho->SetWeight(1.0);
             Omega->SetWeight(1.0); //0.9 from hep-ex/0510017
             Phi->SetWeight(1.0);
	     JPsi->SetWeight(1.0);
	     Upsilon->SetWeight(1.0);
	     */

	     /* 
             Histogram = InitializePtPHENIX(pt_setup,Piplus->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Piplus->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Piminus->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Piminus->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Kplus->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Kplus->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Kminus->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Kminus->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Proton->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Proton->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Antiproton->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Antiproton->SetPtHistogram(Histogram);
	     */
             Histogram = InitializePtPHENIX(pt_setup,Pion->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Pion->SetPtHistogram(Histogram);

	     
             Histogram = InitializePtPHENIX(pt_setup,Eta->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
	     
             Eta->SetPtHistogram(Histogram);


	     Histogram = InitializePtPHENIX(pt_setup,directPion->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
	     directPion->SetPtHistogram(Histogram);
	     cout << "direct ID " << directPion->GetID() <<endl;


	     /*
             Histogram = InitializePtPHENIX(pt_setup,Etaprime->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
	     
             Etaprime->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Rho->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Rho->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Omega->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Omega->SetPtHistogram(Histogram);
             Histogram = InitializePtPHENIX(pt_setup,Phi->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Phi->SetPtHistogram(Histogram);
	     Histogram = InitializePtPHENIX(pt_setup,JPsi->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             JPsi->SetPtHistogram(Histogram);
	     Histogram = InitializePtPHENIX(pt_setup,Upsilon->GetID(),
					    f_c, f_p0, f_a, f_b, f_n,
					    PPList);
             Upsilon->SetPtHistogram(Histogram);
	     */
	     /*
             Histogram = InitializeYPHENIX(Piplus->GetID(),PPList);
             Piplus->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Piminus->GetID(),PPList);
             Piminus->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Kplus->GetID(),PPList);
             Kplus->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Kminus->GetID(),PPList);
             Kminus->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Proton->GetID(),PPList);
             Proton->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Antiproton->GetID(),PPList);
             Antiproton->SetYHistogram(Histogram);
	     */
             Histogram = InitializeYPHENIX(Pion->GetID(),PPList);
             Pion->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(directPion->GetID(),PPList);
             directPion->SetYHistogram(Histogram);
	     
             Histogram = InitializeYPHENIX(Eta->GetID(),PPList);
             Eta->SetYHistogram(Histogram);
	     /*
             Histogram = InitializeYPHENIX(Etaprime->GetID(),PPList);
             Etaprime->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Rho->GetID(),PPList);
             Rho->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Omega->GetID(),PPList);
             Omega->SetYHistogram(Histogram);
             Histogram = InitializeYPHENIX(Phi->GetID(),PPList);
             Phi->SetYHistogram(Histogram);
	     Histogram = InitializeYPHENIX(JPsi->GetID(),PPList);
             JPsi->SetYHistogram(Histogram);
	     Histogram = InitializeYPHENIX(Upsilon->GetID(),PPList);
             Upsilon->SetYHistogram(Histogram);
	     */

             Histogram = 0;
	     /*
             Piplus->SetMHistogram(Histogram);
             Piminus->SetMHistogram(Histogram);
             Kplus->SetMHistogram(Histogram);
             Kminus->SetMHistogram(Histogram);
             Proton->SetMHistogram(Histogram);
             Antiproton->SetMHistogram(Histogram);
	     */
             Pion->SetMHistogram(Histogram);
             directPion->SetMHistogram(Histogram);
	     
             Eta->SetMHistogram(Histogram);
	     /*
             Etaprime->SetMHistogram(Histogram);
	     Histogram = InitializeM(Rho->GetID(),0.,0.,PPList);
             Rho->SetMHistogram(Histogram);
	     Histogram = InitializeM(Omega->GetID(),0.,0.,PPList);
             Omega->SetMHistogram(Histogram);
	     Histogram = InitializeM(Phi->GetID(),0.,0.,PPList);
             Phi->SetMHistogram(Histogram);
	     Histogram = InitializeM(JPsi->GetID(),3.0927,3.1014,PPList);
             JPsi->SetMHistogram(Histogram);
	     Histogram = InitializeM(Upsilon->GetID(),9.4574,9.4626,PPList);
             Upsilon->SetMHistogram(Histogram);
	     
	     PGList->Insert(Piplus);
	     PGList->Insert(Piminus);
	     PGList->Insert(Kplus);
	     PGList->Insert(Kminus);
	     PGList->Insert(Proton);
	     PGList->Insert(Antiproton);
	     */
	     PGList->Insert(Pion);
	     PGList->Insert(directPion);
	     
	     PGList->Insert(Eta);
	     /*
	     PGList->Insert(Etaprime);
	     PGList->Insert(Rho);
             PGList->Insert(Omega);
	     PGList->Insert(Phi);
	     PGList->Insert(JPsi);
	     PGList->Insert(Upsilon);
	     */

             break;

    default: cout << "Error: Setup " << setup << " not predefined" << endl;

             break;
  }

  cout << endl;

  return PGList;
}
