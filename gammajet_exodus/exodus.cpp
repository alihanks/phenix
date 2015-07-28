//-----------------------------------------------------------------------------
//
//  Decay machine of the EXODUS package
//
//  Main program
//
//-----------------------------------------------------------------------------

#include <stdlib.h>

#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TString.h>
#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticlePropertyList.h"
#include "ParticleGeneratorList.h"
#include "DecayList.h"

TROOT exodus("exodus","Initialize ROOT for exodus");
#define  INCLUDEFLAG
#include "DeclareROOTObjects.h"

TFile        * OpenROOTFile(char *);
TFile        * CloseROOTFile(TFile *);
ofstream     * OpenFullEventFile(char *);
void           CloseFullEventFile(ofstream *);
void           BookROOTObjects();
ParticleGeneratorList * InitializeSetup(int, ParticlePropertyList *,float);
void                    GenerateParticles(ParticleList *, int, 
					  ParticleGeneratorList *,
					  ParticlePropertyList *);
void                    GenerateFullEvent(ParticleList *, int, 
					  ParticleGeneratorList *,
					  ParticlePropertyList *,
					  int pflag = 0);
ParticleList * DoAllDecays(ParticleList *, ParticlePropertyList *,
			   DecayList *);
ParticleList * AdjustDecaySum(ParticleList *, DecayList *);
ParticlePropertyList * DefineParticleProperties(void);
DecayList            * DefineDecayProperties(ParticlePropertyList *);
void FillROOTObjects(int, double, double, char *, 
		     ParticleList *,ParticlePropertyList *);

void WriteFullEventFile(ofstream *, int, ParticleList *, 
			ParticlePropertyList *);

void InitializeRandom();
void CleanupRandom();

int main()
{
  int            setup, events;
  ParticleList * OneEvent = 0; 
  char           output_file[100];
  char           event_file[100];
  int            ievent;
  int            dnch_dy = 100;
  double         dNdy_pi0 = 1.;
  double         N_coll = 1.;
  ofstream     * output_stream = 0;
  TFile        * root_file = 0;

  cout << endl << endl;
  cout << "**********************************" << endl;
  cout << "*                                *" << endl;
  cout << "*  W E L C O M E to E X O D U S  *" << endl;
  cout << "*                                *" << endl;
  cout << "*        GENESIS version         *" << endl;
  cout << "*                                *" << endl;
  cout << "*                                *" << endl;
  cout << "*    Wir kriegen alles klein!    *" << endl;
  cout << "*                                *" << endl;
  cout << "**********************************" << endl;
  cout << endl << endl;

  do
  {
    cout << "Choose one of the predefined setups:" << endl;
    cout << "------------------------------------" << endl;
    cout << endl;
    cout << "1) CERES" << endl;
    cout << "2) ISR" << endl;
    cout << "3) PHENIX" << endl;
    cout << "4) Phi->KK" << endl;
    cout << "5) single-particle generator" << endl;
    cout << "6) PHENIX: complete events" << endl;
    cout << endl;
    cout << "Your choice (1-5): "; 
    cin  >> setup; 
    cout << endl;
  } while ( setup<1 || setup>6 );

  cout << "How many events? ";
  cin  >> events;
  cout << endl;

  cout << "Output file containing ROOT objects: ";
  cin >> output_file;
  cout << endl;

  TString myout(output_file);
  int ptflag = 0;
  if (myout.Contains("flat"))
    ptflag = 1;
  
  cout << "setting up ptflag" << endl;

  if ( setup==3 )
  {
    cout << "dN/dy of pizero: ";
    cin  >> dNdy_pi0;
    cout << endl;
    cout << "N_coll: ";
    cin  >> N_coll; 
    cout << endl;
  }

  if ( setup==6 )
  {
    cout << "Output file containing full events: ";
    cin  >> event_file;
    cout << endl;
    //output_stream = OpenFullEventFile(event_file);
    cout << "Select dN_ch/dy at y=0 (integer value): ";
    cin  >> dnch_dy; 
    cout << endl;
  }

  root_file = OpenROOTFile(output_file);
  BookROOTObjects();
  InitializeRandom();

  if(setup==6) dNdy_pi0=dnch_dy;

  cout<<" define decay "<<endl;
  ParticlePropertyList  * Species       = DefineParticleProperties();
  cout<<" init setup "<<endl;
  ParticleGeneratorList * GeneratorList = InitializeSetup(setup,Species,dNdy_pi0);
  cout<<" decay list "<<endl;
  DecayList             * Decays        = DefineDecayProperties(Species);

  if ( Decays==0 ) return 0;
  
  for ( ievent=1; ievent<=events; ievent++ )
  {
    OneEvent = new ParticleList;
    //if ( fmod((double)ievent,100000.0)==0.0 )
    if ( ievent%100000==0 )
      cout << ievent << " events done" << endl;

    if ( setup!=6 )
    {
      GenerateParticles(OneEvent,1,GeneratorList,Species);
    }
    else
    {
      //cout << "Generating EVENT " << ievent <<endl; 
      GenerateFullEvent(OneEvent,dnch_dy,GeneratorList,Species,ptflag);

      //cout << "Generated full EVENT " << ievent <<endl; 

    }
    //cout << "going to adjust decay sum for EVENT" <<endl;
    OneEvent = AdjustDecaySum(OneEvent,Decays);
    //cout << "now do all decays for EVENT" <<endl;
    OneEvent = DoAllDecays(OneEvent,Species,Decays);
    //cout << "Filling root objects for EVENT" <<endl;
    FillROOTObjects(setup,dNdy_pi0,N_coll,output_file,OneEvent,Species);
    //cout << "exodus: filled EVENT " << ievent <<endl;
    if ( setup==6 )  
      //WriteFullEventFile(output_stream,ievent,OneEvent,Species);
      delete OneEvent;
      OneEvent = 0;
  }
  //cout << "exodus: close file" <<endl;
  root_file = CloseROOTFile(root_file);
  if ( setup==6 )
  {
    //CloseFullEventFile(output_stream);
  }
  CleanupRandom();

  return 0;
}

  







