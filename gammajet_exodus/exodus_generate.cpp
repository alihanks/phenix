//-----------------------------------------------------------------------------
//
//  First-generation generator of the EXODUS package
//
//  Main program
//
//-----------------------------------------------------------------------------

#include <stdlib.h>

#include <TROOT.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticleList.h"
#include "ParticlePropertyList.h"
#include "ParticleGeneratorList.h"

ParticlePropertyList  * DefineParticleProperties(void);
ParticleGeneratorList * InitializeSetup(int, ParticlePropertyList *,float);
void                  * GenerateParticles(ParticleList *, int, 
					  ParticleGeneratorList *,
					  ParticlePropertyList *);
void                  * GenerateSingleParticles(ParticleList *, int, 
						ParticlePropertyList *);
void WriteParticleListToFile(char *, ParticleList *);
void WriteParticleListToOscarFile(char *, ParticleList *, ParticlePropertyList *);
void InitializeRandom();
void CleanupRandom();

int main()
{
  ParticleList * FirstGeneration = new ParticleList; 
  int            setup, events;
  char           output_file[100];
  char           oscar_file[100];

  cout << endl << endl;
  cout << "**********************************" << endl;
  cout << "*                                *" << endl;
  cout << "*  W E L C O M E to E X O D U S  *" << endl;
  cout << "*                                *" << endl;
  cout << "*        FIRST GENERATION        *" << endl;
  cout << "*                                *" << endl;
  cout << "*           GENERATOR            *" << endl;
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
    cout << "5) Single-particle generator" << endl;
    cout << endl;
    cout << "Your choice (1-5): "; 
    cin  >> setup; 
    cout << endl;
  } while ( setup<1 || setup>5 );

  cout << "How many events? ";
  cin  >> events;
  cout << endl;

  if ( setup!=5 )
    {
      cout << "Output file containing first-generation particles: ";
      cin  >> output_file;
      cout << endl;
    }
  else
    {
      cout << "Output file containing first-generation particles (OSCAR compliant): ";
      cin  >> oscar_file;
      cout << endl;
    }

  TROOT exodus("exodus","Initialize ROOT for exodus");
  InitializeRandom();
  cout<<" defining particle properties "<<endl;
  ParticlePropertyList * Species = DefineParticleProperties();
   cout<<" defined particle properties "<<endl;
  ParticleGeneratorList * GeneratorList = InitializeSetup(setup,Species,1.0);

  if ( setup!=5 )
    {
      GenerateParticles(FirstGeneration,events,GeneratorList,Species);
    }
  else
    {
      GenerateSingleParticles(FirstGeneration,events,Species);
    }

  cout << "Lenght of first-generation particle list: "
       << FirstGeneration->GetLength() << " particles" << endl << endl;
  
  if ( setup!=5 )  
    WriteParticleListToFile(output_file,FirstGeneration);

  if ( setup==5 )  
    WriteParticleListToOscarFile(oscar_file,FirstGeneration,Species);

  CleanupRandom();

  return 0;
 
}


