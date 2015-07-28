#ifndef __WARNMAP_H__
#define __WARNMAP_H__

#include <vector>

class Warnmap {

  static const int N_ARM       = 2;
  static const int N_SECTOR    = 4;
  static const int N_ARMSECT   = 8;
  static const int N_YPOS_PBGL = 48;
  static const int N_YPOS_PBSC = 36;
  static const int N_ZPOS_PBGL = 96;
  static const int N_ZPOS_PBSC = 72;
  static const int N_TOWER     = 24768;
  static const int N_SUPERMOD  = 32;

  static const int MASK_HOT         = 0x10;
  static const int MASK_DEAD        = 0x08;
  static const int MASK_AROUND_HOT  = 0x04;
  static const int MASK_AROUND_DEAD = 0x02;
  static const int MASK_EDGE        = 0x01;

  std::vector<int> m_warnmap[N_ARMSECT][N_YPOS_PBGL][N_ZPOS_PBGL];
  std::vector<double> binCut;
  
public:
  Warnmap(const unsigned int bins = 3);
  virtual ~Warnmap(){}
  
  void ResetMap(const unsigned int bins);
  void ReadMap(const char* fname, int bin);
  void SetPtRange(std::vector<double> maxPt);
	
  int IsBad(int as, int y, int z, double pt);
  int IsHot(int as, int ys, int z, int pt);
  int IsDead(int as, int y, int z, int pt);
  int IsAroundHot(int as, int y, int z, int pt);
  int IsAroundDead(int as, int y, int z, int pt);
  int IsEdge(int as, int y, int z, int pt);

  int IsValidYZ(int as, int y, int z);
  
private:
  bool CheckBin(int ptbin);
};

inline int Warnmap::IsBad(int as, int y, int z, double pt)
{
  bool isbad = false;
  for( unsigned int ipt = 0; ipt < binCut.size(); ipt++ )
  {
    if( pt < binCut[ipt] )
      isbad = isbad || (IsHot(as, y, z, ipt)||IsDead(as, y, z, ipt)||IsAroundHot(as, y, z, ipt)||IsAroundDead(as, y, z, ipt)||IsEdge(as, y, z, ipt));
  }
  return isbad;
};
	
inline int Warnmap::IsHot(int as, int y, int z, int pt) 
{
  return (m_warnmap[as][y][z][pt] & MASK_HOT);
};
        
inline int Warnmap::IsDead(int as, int y, int z, int pt)
{
  return (m_warnmap[as][y][z][pt] & MASK_DEAD);
};
	
inline int Warnmap::IsAroundHot(int as, int y, int z, int pt)
{
  return (m_warnmap[as][y][z][pt] & MASK_AROUND_HOT);
};
	
inline int Warnmap::IsAroundDead(int as, int y, int z, int pt)
{
  return (m_warnmap[as][y][z][pt] & MASK_AROUND_DEAD);
};

inline int Warnmap::IsEdge(int as, int y, int z, int pt)
{
  return (m_warnmap[as][y][z][pt] & MASK_EDGE);
};
	
#endif // __WARNMAP_H__
