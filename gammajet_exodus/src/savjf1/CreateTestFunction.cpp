//-----------------------------------------------------------------------------
//
// Create 2D function
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>

#include <TFormula.h>
#include <TF2.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticlePropertyList.h"

TF2 * CreateTestFunction(int ParticleID, ParticlePropertyList *PPList)
{

  TFormula * formula = new TFormula("formula","y*y");
  formula->Update();
  //  TF2 * pty = new TF2("pty","x+formula",0.,4.,-2.,2.);
  TF2 * pty = new TF2("pty","x+formula",0.2,20.,-2.,2.);
  pty->SetNpx(10000);
  pty->SetNpy(10000);

  return pty;

}





