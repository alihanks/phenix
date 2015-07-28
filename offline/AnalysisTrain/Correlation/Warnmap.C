#include "Warnmap.h"
	
	
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

Warnmap::Warnmap(const unsigned int bins)
{
  ResetMap(bins);
}

void Warnmap::ResetMap(const unsigned int bins)
{
  for( int i = 0; i < N_ARMSECT; i++ ) {
    for( int j = 0; j < N_YPOS_PBGL; j++ ) {
      for( int k = 0; k < N_ZPOS_PBGL; k++ ) {
        for( unsigned int z = 0; z < bins; z++ ) {
          m_warnmap[i][j][k].push_back(0);
        }
      }
    }
  }
  //memset(m_warnmap, vector<int>(bins), sizeof(m_warnmap));
}

void Warnmap::SetPtRange(vector<double> maxpT)
{
  for( unsigned int i = 0; i < maxpT.size(); i++ )
  {
    binCut.push_back(maxpT[i]);
  }
}

// Warnmap class keeps track of if/how many pt ranges are looked at separately
void Warnmap::ReadMap(const char* fname, int bin)
{
  std::ifstream if_status(fname);
  if (!(if_status.is_open())) {
    cout << "Warnmap::ReadMap():" << endl
    	 << "  File '" << fname << "' does not exist.  Abort." << endl;
    exit(0);

  }
  
  int armsect, ypos, zpos, status;
  while (if_status >> armsect >> ypos >> zpos >> status) {
    if (IsValidYZ(armsect, ypos, zpos)) {
      if( CheckBin(bin) ) {
        m_warnmap[armsect][ypos][zpos][bin] = status;
      }
    }
  }
  if_status.close();
  return;
}

bool Warnmap::CheckBin(int bin)
{
  if( bin >= 0 && bin < (int)binCut.size() ) return true;
  return false;
}

int Warnmap::IsValidYZ(int as, int y, int z)
{
  int ret = 0;
  if (as == 4 || as == 5) {
    if (y >= 0 && y < N_YPOS_PBGL && z >= 0 && z < N_ZPOS_PBGL) ret = 1;
  } else {
    if (y >=0 && y < N_YPOS_PBSC && z >=0 && z < N_ZPOS_PBSC)   ret = 1;
  }
  return ret;
}
