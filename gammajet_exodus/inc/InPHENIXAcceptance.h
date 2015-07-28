#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#define MOM_CUT     .2
#define JMASS	3.07
#define EMASS   .000511
#define PHIMASS 1.02
#define KMASS   .49365 
#define OMASS   .782
#define PI     acos(-1.)
#define MODE   0 
#if MODE

#define PAI     3.1415926
#define DPAI    6.2831853
#define THETA_TOP 1.9198621
#define THETA_BOT 1.2217305
#define PHI_TOP 2.6179939
#define PHI_BOT .87266463

#else
#define PAI  180
#define DPAI 360
#define THETA_TOP 110
#define THETA_BOT 70
#define PHI_TOP 150
#define PHI_BOT 50

#endif
#define NUMBER     1000000

using namespace std;

/******************** in file accept_particle.cc ***********************/
//this class only exsist in this routine.
class CURPAR
{
public: 
  float p0;
  float p1;
  float p2;
};
typedef struct{
  short arm;
  short half_arm;
  short sector;
} ACCEPT;
class Particletype
{
public:
  float m;
  float p;
  float theta;
  float phi;
  int   charge;
};

ACCEPT  accept_particle(const Particletype&);
ACCEPT  accept_particle2(const Particletype&);

int InPHENIXAcceptance(double p,double theta, double phi,int charge, double m);
/*************************** end here ****************************/

/******************** in file decay.cc  *********************************/


void decay(const Particletype &, Particletype &,Particletype &);
/******************** end *********************************************/







