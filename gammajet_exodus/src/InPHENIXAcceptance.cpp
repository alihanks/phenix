/****************************************************************************
 ****************************************************************************

 Accept_particle
 ---------------

 DESCRIPTION: For given three momentum px,py,pz and charge of a particle
 determine whether it will be accepted or not, if accepted, px,py,pz will 
 be smeared to the values that are actually measured when it passes phenix 
 detector. Return value is 0,1,2,3,4 corresponding to not accepted, accepted by
 one of the four half arm.

 CAUTION: some particle may go pass two half arms so are discarded by this filter,
 I will rewrite this filter so that the return value will tell the sector number
 the half arm number, the arm number by which the particle is accepted.

 ALGORITHM: For px,py,pz, first we change them into sphere varables 
 p,theta,phi. Then we see if p theta phi pass our filter. Filter for p and 
 theta are direct. Filter for phi is eight  simple function line1 - line8.
 they behave as exponential function at low momentum, and p goes to infinity 
 when phi approching some value.
   
 AUTHOR/CONTACT: J. Jiangyong,Stephe Johnson, StonyBrook
 jjia@skipper.physics.sunysb.edu
 snoopy@skipper.physics.sunysb.edu

 REVISIONS:
       Date            Author          Description

       7/17/98         J. Jiangyong     Original
		
 INPUT VARIABLES: px = x component of the momentum at the vertex
                  py = y component of the momentum at the vertex
		  pz = z component of the momentum at the vertex
		  charge = the number of charge of the particle, 1 for positron
		  and -1 for electron.
		  
 OUTPUT VARIABLES: px,py,pz
 we use reference for px,py,pz, they will store the particle momentum components leaving the detector.
 ***************************************************************************
 ***************************************************************************/

 /***************************************************************************
 ***************************************************************************
 Smear
 -----
 
 DESCRIPTION: Smear input value by 1% according to Gaussian distribution.
              rewrite from some code of J. T. Mitchel.

 AUTHOR/CONTACT: J. Jiangyong,Stephe Johnson, StonyBrook

 REVISIONS:
       Date            Author          Description

       7/10/98         J. Jiangyong     Originally type in

 INPUT VARIABLES: p = input value .
                  return value is smeared one.
		  
 ***************************************************************************
 ***************************************************************************/

#include "InPHENIXAcceptance.h"
#include "Momentum.h"

/***********************************************************************/
/******************* description of acceptance filter for phi************/
/************************************************************************ 
   We use eight curve to define our filter, they will take curve parameters
   from Struct array curpar.
   for angles we use: MODE=1  radian; MODE=0   degree
   ***********************************************************************/
#if MODE  
static CURPAR curpar[20]={
{  2.864789,       -.58904862,        .14},\
{  2.4064227,      -.19634954,        .15},\
{  2.4064227,       .19634954,        .16},\
{  2.4064227,       .58904862,        .17},\
{  2.6356059,       .9817477,         .18},\
{  2.864789,       2.1598449,         .14},\
{  2.4064227,      2.552544,          .15},\
{  2.4064227,      2.9452431,         .16},\
{  2.4064227,      3.3379422,         .17},\
{  2.6356059,      3.7306413,         .18},\
{ -2.6356059,      -.58904862,        .14},\
{ -2.4064227,      -.19634954,        .15},\
{ -2.4064227,       .19634954,        .16},\
{ -2.4064227,       .58904862,        .17},\
{ -2.864789,        .9817477,         .18},\
{ -2.6356059,      2.1598449,         .14},\
{ -2.4064227,      2.552544,          .15},\
{ -2.4064227,      2.9452431,         .16},\
{ -2.4064227,      3.3379422,         .17},\
{ -2.864789,       3.7306413,         .18},\
};
#else
static CURPAR curpar[20]={
  {  .050,    -33.75,      .14},\
  {  .042,    -11.25,      .15},\
  {  .042,     11.25,      .16},\
  {  .042,     33.75,      .17},\
  {  .046,     56.25,      .18},\
  {  .050,     123.75,     .14},\
  {  .042,     146.25,     .15},\
  {  .042,     168.75,     .16},\
  {  .042,     191.25,     .17},\
  {  .046,     213.75,     .18},\
  { -.046,    -33.75,      .14},\
  { -.042,    -11.25,      .15},\
  { -.042,     11.25,      .16},\
  { -.042,     33.75,      .17},\
  { -.050,     56.25,      .18},\
  { -.046,     123.75,     .14},\
  { -.042,     146.25,     .15},\
  { -.042,     168.75,     .16},\
  { -.042,     191.25,     .17},\
  { -.050,     213.75,     .18},\
};
#endif

static float line1(const float &phi)
{
  if(phi>=curpar[0].p1) return NUMBER; 
  return .5*(1/(1-exp((curpar[0].p0)*(phi-curpar[0].p1)))-1)+curpar[0].p2;
}
static float line2(const float &phi)
{
  if(phi>=curpar[1].p1) return NUMBER;
  return .5*(1/(1-exp((curpar[1].p0)*(phi-curpar[1].p1)))-1)+curpar[1].p2;
}
static float line3(const float &phi)
{
  if(phi>=curpar[2].p1) return NUMBER;
  return .5*(1/(1-exp((curpar[2].p0)*(phi-curpar[2].p1)))-1)+curpar[2].p2;
}
static float line4(const float &phi)
{
  if(phi>=curpar[3].p1) return NUMBER;
  return .5*(1/(1-exp((curpar[3].p0)*(phi-curpar[3].p1)))-1)+curpar[3].p2;
}
static float line5(const float &phi)
{
  if(phi>=curpar[4].p1) return NUMBER; 
  return .5*(1/(1-exp((curpar[4].p0)*(phi-curpar[4].p1)))-1)+curpar[4].p2;
}
static float line6(const float &phi)
{
  if(phi>=curpar[5].p1) return NUMBER; 
  return .5*(1/(1-exp((curpar[5].p0)*(phi-curpar[5].p1)))-1)+curpar[5].p2;
}
static float line7(const float &phi)
{
  if(phi>=curpar[6].p1) return NUMBER; 
  return .5*(1/(1-exp((curpar[6].p0)*(phi-curpar[6].p1)))-1)+curpar[6].p2;
}
static float line8(const float &phi)
{
  if(phi>=curpar[7].p1) return NUMBER; 
  return .5*(1/(1-exp((curpar[7].p0)*(phi-curpar[7].p1)))-1)+curpar[7].p2;
}
static float line9(const float &phi)
{
  if(phi>=curpar[8].p1) return NUMBER; 
  return .5*(1/(1-exp((curpar[8].p0)*(phi-curpar[8].p1)))-1)+curpar[8].p2;
}
static float line10(const float &phi)
{
  if(phi>=curpar[9].p1) return NUMBER; 
  return .5*(1/(1-exp((curpar[9].p0)*(phi-curpar[9].p1)))-1)+curpar[9].p2;
}

static float line11(const float &phi)
{
  if(phi<=curpar[10].p1) return -NUMBER;
  return  -.5*(1/(1-exp((curpar[10].p0)*(phi-curpar[10].p1)))-1)-curpar[10].p2;
}
static float line12(const float &phi)
{
  if(phi<=curpar[11].p1) return -NUMBER;
  return -.5*(1/(1-exp((curpar[11].p0)*(phi-curpar[11].p1)))-1)-curpar[11].p2;
}
static float line13(const float &phi)
{
  if(phi<=curpar[12].p1) return -NUMBER;
  return -.5*(1/(1-exp((curpar[12].p0)*(phi-curpar[12].p1)))-1)-curpar[12].p2;
}
static float line14(const float &phi)
{
  if(phi<=curpar[13].p1) return -NUMBER;
  return -.5*(1/(1-exp((curpar[13].p0)*(phi-curpar[13].p1)))-1)-curpar[13].p2;
}
static float line15(const float &phi)
{
  if(phi<=curpar[14].p1) return -NUMBER;
  return -.5*(1/(1-exp((curpar[14].p0)*(phi-curpar[14].p1)))-1)-curpar[14].p2;
}
static float line16(const float &phi)
{
  if(phi<=curpar[15].p1) return -NUMBER;
  return -.5*(1/(1-exp((curpar[15].p0)*(phi-curpar[15].p1)))-1)-curpar[15].p2;
}
static float line17(const float &phi)
{
  if(phi<=curpar[16].p1) return -NUMBER;
  return -.5*(1/(1-exp((curpar[16].p0)*(phi-curpar[16].p1)))-1)-curpar[16].p2;
}
static float line18(const float &phi)
{
  if(phi<=curpar[17].p1) return -NUMBER;
  return  -.5*(1/(1-exp((curpar[17].p0)*(phi-curpar[17].p1)))-1)-curpar[17].p2;
}
static float line19(const float &phi)
{
  if(phi<=curpar[18].p1) return -NUMBER;
  return -.5*(1/(1-exp((curpar[18].p0)*(phi-curpar[18].p1)))-1)-curpar[18].p2;
}
static float line20(const float &phi)
{
  if(phi<=curpar[19].p1) return -NUMBER;
  return  -.5*(1/(1-exp((curpar[19].p0)*(phi-curpar[19].p1)))-1)-curpar[19].p2;
}
/*************** defination of phi & momentum filter end here***************/


/**********************************************************************
********* function smear will randomize the initial momentum according**
********* to momentum    resolution   *******************************/ 
float smear(float p)
{ 
  float v1, v2,f, rsq,sigma;
  sigma=.01*p;

  /* Generate Gaussian random value.  
   algorithm by J.T.MITCHEL*/
  do {
    v1 = 2.0*drand48()-1.0;
    v2 = 2.0*drand48()-1.0;
    rsq = v1*v1 + v2*v2;
  } while (rsq >=1.0 || rsq == 0);
  f = sqrt(-2.0*log(rsq)/rsq);

  p += v2*f*sigma;
  return p;
}

ACCEPT accept_particle2(const Particletype& x)
{
  ACCEPT Revalue={-1,-1,-1};
  float mom=x.p,theta=x.theta,phi=x.phi,charge=x.charge;

  // define function array: left and right,  
  static float (*cur[20])(const float &)={line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12,line13,line14,line15,line16,line17,line18,line19,line20};
  // static float (*right[8])(const float &)={line2,line4,line6,line8,line10,line12,line14,line16};
  
  //deal with neutral particle.
  if(charge==0)
    { 
      if(mom<0)mom=-mom;
      if(mom<MOM_CUT) return Revalue;//momentum cut
      if(theta>THETA_TOP||theta<THETA_BOT) return Revalue;//theta cut
      if(phi<-PHI_BOT)phi+=DPAI;
      if(phi>curpar[0].p1&&phi<curpar[4].p1)
	{
	  Revalue.arm=1;
	  if(phi<curpar[1].p1)Revalue.sector=4;
	  else if(phi<curpar[2].p1)Revalue.sector=5;
	  else if(phi<curpar[3].p1)Revalue.sector=6;
	  else Revalue.sector=7;
	}
      else if(phi>curpar[5].p1&&phi<curpar[9].p1)
	{
	  Revalue.arm=0;
	  if(phi<curpar[6].p1) Revalue.sector=0;
	  else if(phi<curpar[7].p1)Revalue.sector=1;
	  else if(phi<curpar[8].p1)Revalue.sector=2;
	  else Revalue.sector=3;
	}
      return Revalue;
    }
  /*now when |charge|>1 we need to adjust mom, so we can still use the same filter**/
  mom/=charge;
  if(mom<MOM_CUT&&mom>-MOM_CUT)return Revalue;//momentum cut
  if(theta>THETA_TOP||theta<THETA_BOT)return Revalue;//theta cut
  
  //phi and momentum filter
  if(charge>0)
    {
      if(mom<0)mom=-mom;
      if(phi<-PHI_TOP)phi+=DPAI;
      if(mom>cur[4](phi)&&mom<cur[0](phi))
      {
	Revalue.arm=1;
	if(mom>cur[1](phi)) Revalue.sector=4;
	else if(mom>cur[2](phi)) Revalue.sector=5;
	else if(mom>cur[3](phi)) Revalue.sector=6;
	else Revalue.sector=7;
      }
      else if(mom>cur[9](phi)&&mom<cur[5](phi))
	{
	  Revalue.arm=0;
	  if(mom>cur[6](phi)) Revalue.sector=0;
	  else if(mom>cur[7](phi)) Revalue.sector=1;
	  else if(mom>cur[8](phi)) Revalue.sector=2;
	  else Revalue.sector=3;
      }
      return Revalue;
    }
  else if(charge<0)
    {
      if(mom>0)mom=-mom;
      if(phi<-PHI_BOT)phi+=DPAI;
      if(mom>cur[14](phi)&&mom<cur[10](phi))
	{
	  Revalue.arm=1;
	  if(mom>cur[11](phi)) Revalue.sector=4;
	  else if(mom>cur[12](phi)) Revalue.sector=5;
	  else if(mom>cur[13](phi)) Revalue.sector=6;
	  else Revalue.sector=7;
	}
      if(mom>cur[19](phi)&&mom<cur[15](phi))
	{
	  Revalue.arm=0;
	  if(mom>cur[16](phi)) Revalue.sector=0;
	  else if(mom>cur[17](phi)) Revalue.sector=1;
	  else if(mom>cur[18](phi)) Revalue.sector=2;
	  else Revalue.sector=3;
	}
      return Revalue;
    }

  cout << "this line should never be reached!" << endl;
  return Revalue;

}


double arg(double x, double y)	{
double u=180./3.1415926535897932384626;
	if ((x==0)&&(y==0)) return 0;
	if ((x==0)&&(y>0)) return 90;
	if ((x==0)&&(y<0)) return 270;
	if (x>0) {
		if (y>=0) return u*atan(y/x);
		return 360-u*atan(-y/x);
	}
	if (x<0) {
		if (y>=0) return 180-u*atan(y/-x);
		return 180+u*atan(y/x);
	}
	cout<<"ERROR IN arg()";
	return -1;
}	


int InPHENIXAcceptance(double p, double theta, double phi, 
		       int charge, double m=0)
{
  Particletype teil;
  teil.m=m;
  teil.p=p;
  teil.theta=theta;
  teil.phi=phi-90.;
  //  teil.phi=phi;
  teil.charge=charge;
  return accept_particle2(teil).sector;
}

bool InPHENIXAcceptance(Mom4 mom, int charge)
{

Particletype teil;
Mom3 p=mom.Getp();

teil.p=sqrt(p.Getpx()*p.Getpx()+p.Getpy()*p.Getpy()+p.Getpz()*p.Getpz());
teil.m=sqrt(mom.GetE()*mom.GetE()-teil.p*teil.p);
teil.theta=arg(p.Getpz(),sqrt(p.Getpx()*p.Getpx()+p.Getpy()*p.Getpy()));
teil.phi=arg(p.Getpx(),p.Getpy());
teil.charge=charge;

 return (accept_particle2(teil).arm==1);
}
