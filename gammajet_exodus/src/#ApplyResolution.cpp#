//-----------------------------------------------------------------------------
//
//  Call resolution routines for different setups
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <iostream>
#include <TH1.h>
#include "Momentum.h"

using namespace std;

Mom4 ApplyCERESResolution(Mom4);
Mom4 ApplyPHENIXResolution(Mom4);

Mom4 ApplyResolution(int setup, Mom4 mom4)
{
  Mom4   result;

  switch(setup)
  {
    case 1:  result = ApplyCERESResolution(mom4);

             break;

    case 2:  result = ApplyCERESResolution(mom4);

             break;

    case 3:  result = ApplyPHENIXResolution(mom4);

             break;

    case 6:  result = mom4;

             break;

    default: cout << "Error: Setup " << setup << " not predefined" << endl;

             break;
  }

  return (result);

}















