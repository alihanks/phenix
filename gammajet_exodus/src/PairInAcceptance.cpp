//-----------------------------------------------------------------------------
//
//  Call pair acceptance routines for different setups
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "Momentum.h"

using namespace std;

bool InCERESAcceptance(Mom4, Mom4);

bool PairInAcceptance(int setup, Mom4 mom1, Mom4 mom2)
{
  bool accept = 0;

  switch(setup)
  {
    case 1:  accept = InCERESAcceptance(mom1,mom2);

             break;

    case 2:  accept = InCERESAcceptance(mom1,mom2);

             break;

    case 3:  accept = InCERESAcceptance(mom1,mom2);

             break;

    default: cout << "Error: Setup " << setup << " not predefined" << endl;

             break;
  }

  return (accept);

}

