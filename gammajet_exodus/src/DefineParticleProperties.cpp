//-----------------------------------------------------------------------------
//
//  Read particle properties from file 
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <stdlib.h>
#include <string>

#include "Momentum.h"
#include "Particle.h"
#include "ParticlePropertyList.h"

ParticlePropertyList * DefineParticleProperties(void)
{
  double mass,width;
  int    id;
  int    charge, spin;
  char   ch;
  string label, name;

  ParticlePropertyList * DefinedParticles = new ParticlePropertyList;

  cout << "Reading particle definitions from file: defined_particles.txt"
       << endl;

  ifstream fin("defined_particles.txt");
  while ( fin.get(ch) ) 
  {
    fin >> label;
    if ( label=="@{" )
    {
      fin >> name;
      for(;;)
      {
        fin >> label;
	if ( label=="ID=" )     fin >> id;
	if ( label=="Mass=" )   fin >> mass;
	if ( label=="Width=" )  fin >> width;
	if ( label=="Charge=" ) fin >> charge;
	if ( label=="Spin=" )   fin >> spin;
	if ( label=="}" )       break;
      }
    }
    else continue;
    ParticleProperty * pPart = new ParticleProperty;
    pPart->SetName(name);
    pPart->SetID(id);
    pPart->SetMass(mass);
    pPart->SetWidth(width);
    pPart->SetCharge(charge);
    pPart->SetSpin(spin);
    DefinedParticles->Insert(pPart);
    cout << "Defined: " << name << endl;
  }
  fin.close();
  
  cout << endl;
  
  return DefinedParticles;
}
























