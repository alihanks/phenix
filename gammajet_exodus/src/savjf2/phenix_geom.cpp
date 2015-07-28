/* Standard includes */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <istream>
#include <fstream>
#include <TRandom.h>
#include <iostream>

#include <TH2.h>

#define  INCLUDEFLAG extern
#include "DeclareROOTObjects.h"

using namespace std;


int output_num[16];
int evaluate_max(int,float,float,float);
int Calculate_cross_point(float,float,float,int,float*,float*,float*);
int convert_to_local_and_check(float,float,float,int, float & smeared_x, float & smeared_y, float & smeared_z, float, int);
void Construct_EMCal(char);
int check_bad_tile(int,int,int);

static int checkhit_Sc[6][72][36];
static int checkhit_Gl[2][96][48];
//new for live tower weighting
static int checkhit_Sc_weight[6][72][36];
static int checkhit_Gl_weight[2][96][48];
float deg = M_PI/180.;
/* Sector dimensions */
int NxSc = 72;
int NySc = 36;
int NxGl = 96;
int NyGl = 48;
/* Zero tower position in sector frame defined by geodesists (cm) */
float ZeroTowerSc[3]={2.859, 2.859, -48.};
float ZeroTowerGl[3]={0., 0.,   0.};
/* Tower Size (cm) */
float TowerSizeSc = 5.562;
float TowerSizeGl = 4.068;
/* Gaps between sectors (cm) */
float GapXSc = 0.256;
float GapYSc = 0.156;
float GapXGl = 0.150;
float GapYGl = 0.163;
static int first_flag=1;

/* Rotation angles (deg) */
float emcRA[8] =
{
     22.47,                   /*  W0 */
     -2.3e-3,                 /*  W1 */
     -22.486,                 /*  W2 */
     -45.010,                 /*  W3 */

     -45.016,                 /*  E3 */
     -22.502,                 /*  E2 */
     0.189,                   /*  E1 */
     22.5                     /*  E0 */
};
/* Translation vector */
float emcTR1[8][3] =
{
    {4.745e+2, -3.053e+2, 2.012e+2},         /*  W0 */
    {5.552e+2, -1.005e+2, 2.010e+2},         /*  W1 */
    {5.513e+2,  1.196e+2, 2.011e+2},         /*  W2 */
    {4.636e+2,  3.215e+2, 2.011e+2},         /*  W3 */

    {-4.640e+2,  3.215e+2, -2.009e+2},       /*  E3 */
    {-5.516e+2,  1.198e+2, -2.012e+2},       /*  E2 */
    {-5.398e+2, -0.978e+2, -1.938e+2},       /*  E1 */
    {-4.606e+2, -2.965e+2, -1.934e+2}        /*  E0 */
};
float emcTR2[8][3];
float emcTR3[8][3];
float emcTR4[8][3];
float emcTR5[8][3];

float sxSc,sySc,sxGl,syGl;
float posx,posy,posz;
static float eta_max[3]={-5.0,-5.0,-5.0};
static float eta_min[3]={5.0,5.0,5.0};
static float phi_max[3]={-100.0,-100.0,-100.0};
static float phi_min[3]={450.0,450.0,450.0};
static float x_max[3]={-10000.0,-10000.0,-10000.0};
static float y_max[3]={-10000.0,-10000.0,-10000.0};
static float z_max[3]={-10000.0,-10000.0,-10000.0};
static float x_min[3]={10000.0,10000.0,10000.0};
static float y_min[3]={10000.0,10000.0,10000.0};
static float z_min[3]={10000.0,10000.0,10000.0};
float theta,eta,phi;
float A[8],B[8],C[8];

int i,flag;

void Construct_EMCal(char *deadfile)
{
  int iS,j,k;
  float Edge_end_x,Edge_end_y,Edge_end_z;
  float Total_Sc_z,Total_Sc_y,Total_Gl_z,Total_Gl_y;
  float Offset_Sc_y,Offset_Sc_z,Offset_Gl_y,Offset_Gl_z;
  float X1,X2,Y1,Y2,Z1,Z2;
  int   Xpos,Ypos,sector,arm,keycent,iz,iy;

  // Calculate Towersize with implmentation of small gaps between modules
  sxSc = ( (NxSc-1)*TowerSizeSc + 5*GapXSc ) / (NxSc-1);
  sySc = ( (NySc-1)*TowerSizeSc + 2*GapYSc ) / (NySc-1);
  sxGl = ( (NxGl-1)*TowerSizeGl + 5*GapXGl ) / (NxGl-1);
  syGl = ( (NyGl-1)*TowerSizeGl + 2*GapYGl ) / (NyGl-1);

  cout << "sxSc= " << sxSc << endl;
  cout << "sySc= " << sySc << endl;

  // They are the total size of Sector
  Total_Sc_z=sxSc*72.0;
  Total_Sc_y=sySc*36.0;
  Total_Gl_z=sxGl*96.0;
  Total_Gl_y=syGl*48.0;

  Offset_Sc_y= Total_Sc_y;
  Offset_Sc_z= Total_Sc_z;

  Offset_Gl_y= Total_Gl_y;
  Offset_Gl_z= Total_Gl_z;

  ifstream fin(deadfile);
  if(!fin){cout << "No deadmap file of" << deadfile << endl; exit(100);}

  // At first there should not be any dead channels in it
  //j: iz in emcClusterContainer, k: iy in emcClusterContainer
  //  for(iS=0;iS<6;iS++)for(j=0;j<72;j++)for(k=0;k<36;k++)checkhit_Sc[iS][j][k]=1;
  //  for(iS=0;iS<2;iS++)for(j=0;j<96;j++)for(k=0;k<48;k++)checkhit_Gl[iS][j][k]=1;

  //switched to live towers -- default is 0 instead
    for(iS=0;iS<6;iS++)for(j=0;j<72;j++)for(k=0;k<36;k++)checkhit_Sc[iS][j][k]=0;
    for(iS=0;iS<2;iS++)for(j=0;j<96;j++)for(k=0;k<48;k++)checkhit_Gl[iS][j][k]=0;
    for(iS=0;iS<6;iS++)for(j=0;j<72;j++)for(k=0;k<36;k++)checkhit_Sc_weight[iS][j][k]=0;
    for(iS=0;iS<2;iS++)for(j=0;j<96;j++)for(k=0;k<48;k++)checkhit_Gl_weight[iS][j][k]=0;



  /* for Takao's Run2 definition
  while(fin >> keycent)
    {
      arm = keycent/100000;
      sector = (keycent%100000)/10000;
      if (arm==1) sector = 7 - sector;
      Ypos = ((keycent % 100000) % 10000) / 100;
      Xpos = ((keycent % 100000) % 10000) % 100;
      checkhit_Sc[sector][Xpos][Ypos]=0;
    }
   */

  /* for Run4 Fast Track
  while (fin >> sector >> iz >> iy)
    {
      if(sector<6){
	checkhit_Sc[sector][iz][iy] = 0;
      }
      else{
	checkhit_Gl[sector-6][iz][iy] = 0;
      }
    }
  */

  /* for Run4   */
  /*
  int ishotdead, twrid;
    cout<<" loading dead map "<<endl;
  while(fin>>twrid>>ishotdead){
    // cout<<" twrid "<<twrid<<" ishotdead "<<ishotdead<<endl;
    if(ishotdead!=1){
      if (twrid<15552){
	//  First 15552 towers are PbSc in section 0,1,2,3,6,7
	sector=(int)twrid/2592;  // get sector
	if (sector > 3) {sector = sector+2;}  // correct for sector shift 
	iy=(int)(twrid%2592)/72;  //get iy
	iz=(int)(twrid%2592)%72;  //get iz
	checkhit_Sc[sector][iz][iy] = 0;
      }else{
	//  Rest are from PbGl sectors 4,5
	sector=(int)((twrid-15552)/4608) + 4;  // get sector (+4 to get PbGl)
	iy=(int)((twrid-15552)%4608)/96;  // get iy
	iz=(int)((twrid-15552)%4608)%96;  // get iz
	checkhit_Gl[sector-6][iz][iy] = 0;
      }
    }
  }

  cout<<" loaded dead map "<<endl;
  */

  //beware, for matt's tower number convention

  int ishotdead, twrid, effweight;
    cout<<" loading dead map "<<endl;
    while(fin>>twrid>>effweight){
      

      cout << twrid << " " << effweight << endl;
    if (twrid<15552){
      //  First 15552 towers are PbSc in section 0,1,2,3,6,7
      sector=(int)twrid/2592;  // get sector

      //if (sector > 3) {sector = sector+2;}  // correct for sector shift 

      //need sector swap here
      if(sector==4) sector=-1;
      if(sector==5) sector=4;
      if(sector==-1) sector=5;

      iy=(int)(twrid%2592)/72;  //get iy
      iz=(int)(twrid%2592)%72;  //get iz
      //checkhit_Sc[sector][iz][iy] = 0;
      //switch to live towers
      checkhit_Sc[sector][iz][iy] = 1;
      checkhit_Sc_weight[sector][iz][iy] = effweight;
    }else{
      //  Rest are from PbGl sectors 4,5
      sector=(int)((twrid-15552)/4608);  // get sector (+4 to get PbGl)

      //need sector swap here
      if(sector==1) sector=-1;
      if(sector==0) sector=1;
      if(sector==-1) sector=0;
      
      iy=(int)((twrid-15552)%4608)/96;  // get iy
      iz=(int)((twrid-15552)%4608)%96;  // get iz
      //checkhit_Gl[sector][iz][iy] = 0;
      //switch to live towers
      checkhit_Gl[sector][iz][iy] = 1;
      checkhit_Gl_weight[sector][iz][iy] = effweight;
    }

    }


  cout<<" loaded dead map "<<endl;

  /*  
  while (fin >> arm >> sector >> iz >> iy)
    {
      sector = arm == 0 ? sector : 7-sector;
      if(sector<6){
	checkhit_Sc[sector][iz][iy] = 0;
      }
      else{
	checkhit_Gl[sector-6][iz][iy] = 0;
      }
    }
  */    
  // For Lead Scintillator (West) Geometry
  for(iS=0;iS<4;iS++){
        /* This is the very point of the sector edge after rotation */
        emcTR5[iS][0]=emcTR1[iS][0]+cos(emcRA[iS]*deg)*ZeroTowerSc[2]-sin(emcRA[iS]*deg)*ZeroTowerSc[1];
        emcTR5[iS][1]=emcTR1[iS][1]-sin(emcRA[iS]*deg)*ZeroTowerSc[2]-cos(emcRA[iS]*deg)*ZeroTowerSc[1];
        emcTR5[iS][2]=emcTR1[iS][2]+ZeroTowerSc[0];

        Edge_end_x=emcTR5[iS][0]+sin(emcRA[iS]*deg)*Offset_Sc_y;
        Edge_end_y=emcTR5[iS][1]+cos(emcRA[iS]*deg)*Offset_Sc_y;
        Edge_end_z=emcTR5[iS][2]-Offset_Sc_z;

        emcTR2[iS][0]=emcTR5[iS][0];
        emcTR2[iS][1]=emcTR5[iS][1];
        emcTR2[iS][2]=Edge_end_z;

        emcTR3[iS][0]=Edge_end_x;
        emcTR3[iS][1]=Edge_end_y;
        emcTR3[iS][2]=emcTR5[iS][2];

        emcTR4[iS][0]=Edge_end_x;
        emcTR4[iS][1]=Edge_end_y;
        emcTR4[iS][2]=Edge_end_z;

        // Calculate the equation of plane in 3D space 
	X1=emcTR2[iS][0]-emcTR5[iS][0];
	X2=emcTR3[iS][0]-emcTR5[iS][0];

	Y1=emcTR2[iS][1]-emcTR5[iS][1];
	Y2=emcTR3[iS][1]-emcTR5[iS][1];

	Z1=emcTR2[iS][2]-emcTR5[iS][2];
	Z2=emcTR3[iS][2]-emcTR5[iS][2];

	B[iS]=1; C[iS]=0.0;
	A[iS]=-1*(B[iS]*Y2+C[iS]*Z2)/X2;

/*	printf("X1=%f,X2=%f,Y1=%f,Y2=%f,Z1=%f,Z2=%f\n",X1,X2,Y1,Y2,Z1,Z2); */

/*        printf("A=%f,B=%f,C=%f\n",A,B,C);
	printf("TR1=%f,TR2=%f,TR3=%f\n",emcTR1[iS][0],emcTR1[iS][1],emcTR1[iS][2]); */

//        for(j=0;j<1;j++) for(k=0;k<36;k++) checkhit_Sc[iS][j][k]=1;
//        for(j=71;j<72;j++) for(k=0;k<36;k++) checkhit_Sc[iS][j][k]=1;
//        for(j=35;j<40;j++) for(k=20;k<30;k++) checkhit_Sc[iS][j][k]=0;
  }

  // For Lead Scintillator (East) Geometry
  for(iS=4;iS<6;iS++){
        /* This is the very point of the sector edge after rotation */
        emcTR5[iS][0]=emcTR1[iS][0]-cos(emcRA[iS]*deg)*ZeroTowerSc[2]-sin(emcRA[iS]*deg)*ZeroTowerSc[1];
        emcTR5[iS][1]=emcTR1[iS][1]-sin(emcRA[iS]*deg)*ZeroTowerSc[2]+cos(emcRA[iS]*deg)*ZeroTowerSc[1];
        emcTR5[iS][2]=emcTR1[iS][2]-ZeroTowerSc[0];

        Edge_end_x=emcTR5[iS][0]-sin(emcRA[iS]*deg)*Offset_Sc_y;
        Edge_end_y=emcTR5[iS][1]+cos(emcRA[iS]*deg)*Offset_Sc_y;
        Edge_end_z=emcTR5[iS][2]+Offset_Sc_z;

        emcTR2[iS][0]=emcTR5[iS][0];
        emcTR2[iS][1]=emcTR5[iS][1];
        emcTR2[iS][2]=Edge_end_z;

        emcTR3[iS][0]=Edge_end_x;
        emcTR3[iS][1]=Edge_end_y;
        emcTR3[iS][2]=emcTR5[iS][2];

        emcTR4[iS][0]=Edge_end_x;
        emcTR4[iS][1]=Edge_end_y;
        emcTR4[iS][2]=Edge_end_z;

        // Calculate the equation of plane in 3D space 
	X1=emcTR2[iS][0]-emcTR5[iS][0];
	X2=emcTR3[iS][0]-emcTR5[iS][0];

	Y1=emcTR2[iS][1]-emcTR5[iS][1];
	Y2=emcTR3[iS][1]-emcTR5[iS][1];

	Z1=emcTR2[iS][2]-emcTR5[iS][2];
	Z2=emcTR3[iS][2]-emcTR5[iS][2];

	B[iS]=1; C[iS]=0.0;
	A[iS]=-1*(B[iS]*Y2+C[iS]*Z2)/X2;

/*	printf("X1=%f,X2=%f,Y1=%f,Y2=%f,Z1=%f,Z2=%f\n",X1,X2,Y1,Y2,Z1,Z2); */

/*        printf("A=%f,B=%f,C=%f\n",A,B,C);
	printf("TR1=%f,TR2=%f,TR3=%f\n",emcTR1[iS][0],emcTR1[iS][1],emcTR1[iS][2]); */

//        for(j=0;j<1;j++) for(k=0;k<36;k++) checkhit_Sc[iS][j][k]=1;
//        for(j=71;j<72;j++) for(k=0;k<36;k++) checkhit_Sc[iS][j][k]=1;
//        for(j=35;j<40;j++) for(k=20;k<30;k++) checkhit_Sc[iS][j][k]=0;
  }

  // For Lead Glass (East) Geometry
  for(iS=6;iS<8;iS++){
        /* This is the very point of the sector edge after rotation */
        emcTR5[iS][0]=emcTR1[iS][0]-cos(emcRA[iS]*deg)*ZeroTowerGl[2]-sin(emcRA[iS]*deg)*ZeroTowerGl[1];
        emcTR5[iS][1]=emcTR1[iS][1]-sin(emcRA[iS]*deg)*ZeroTowerGl[2]+cos(emcRA[iS]*deg)*ZeroTowerGl[1];
        emcTR5[iS][2]=emcTR1[iS][2]-ZeroTowerGl[0];

        Edge_end_x=emcTR5[iS][0]-sin(emcRA[iS]*deg)*Offset_Gl_y;
        Edge_end_y=emcTR5[iS][1]+cos(emcRA[iS]*deg)*Offset_Gl_y;
        Edge_end_z=emcTR5[iS][2]+Offset_Gl_z;

        emcTR2[iS][0]=emcTR5[iS][0];
        emcTR2[iS][1]=emcTR5[iS][1];
        emcTR2[iS][2]=Edge_end_z;

        emcTR3[iS][0]=Edge_end_x;
        emcTR3[iS][1]=Edge_end_y;
        emcTR3[iS][2]=emcTR5[iS][2];

        emcTR4[iS][0]=Edge_end_x;
        emcTR4[iS][1]=Edge_end_y;
        emcTR4[iS][2]=Edge_end_z;

        // Calculate the equation of plane in 3D space 
	X1=emcTR2[iS][0]-emcTR5[iS][0];
	X2=emcTR3[iS][0]-emcTR5[iS][0];

	Y1=emcTR2[iS][1]-emcTR5[iS][1];
	Y2=emcTR3[iS][1]-emcTR5[iS][1];

	Z1=emcTR2[iS][2]-emcTR5[iS][2];
	Z2=emcTR3[iS][2]-emcTR5[iS][2];

	B[iS]=1; C[iS]=0.0;
	A[iS]=-1*(B[iS]*Y2+C[iS]*Z2)/X2;
  }
 
}


int convert_to_local_and_checkhit(float cross_x,float cross_y,float cross_z,int iS, float & smeared_x, float & smeared_y, float & smeared_z, float E, int gn)
{

  //cout << "gn: " << gn <<endl;

   float local_x,local_y,local_z;
   int num_ptr_x,num_ptr_y;
   if(iS<4){
       if(cross_x<=0) return(0);
       //cout <<"crossx: " << cross_x <<endl;
       local_x= cross_x-emcTR1[iS][0]-cos(emcRA[iS]*deg)*ZeroTowerSc[2]+sin(emcRA[iS]*deg)*ZeroTowerSc[1];
       local_y= cross_y-emcTR1[iS][1]+sin(emcRA[iS]*deg)*ZeroTowerSc[2]+cos(emcRA[iS]*deg)*ZeroTowerSc[1];
       local_z= cross_z-emcTR1[iS][2]-ZeroTowerSc[0];

       //cout <<"crossx_mid: " << local_x <<endl;

       local_x/= sin(emcRA[iS]*deg);
       local_y/= cos(emcRA[iS]*deg);
       local_z/= -1;

       //cout <<"localx: " << local_x <<endl;


       //
       // Here comes Fiducial Area Cuts(old version)
       //
       //if(local_z<8 || local_z>389 || local_y<2 || local_y>191) return(0);

       //
       // Get Pointer to the position of calorimeter. Also check if the
       //  pointer points within fiducial area
       //

       if(local_z<0) return(0);
       if(local_y<0) return(0);

       num_ptr_x=(int)(local_z/sxSc); num_ptr_y=(int)(local_y/sySc);
       //if(num_ptr_x>=72||num_ptr_y>=36||num_ptr_x<0||num_ptr_y<0) return(0);

       //use edge cut:
       if(num_ptr_x>=70||num_ptr_y>=34||num_ptr_x<2||num_ptr_y<2) return(0);//2004.10.01

       //cout <<local_z << endl;

       //
       // Check the dead channel map.
       //   "checkhit_Sc[][][]==1" means that it is live channel.
       //

	 /*
       if(checkhit_Sc[iS][num_ptr_x+0][num_ptr_y+0]==1 &&
	  checkhit_Sc[iS][num_ptr_x+0][num_ptr_y+1]==1 &&
	  checkhit_Sc[iS][num_ptr_x+0][num_ptr_y-1]==1 &&
	  checkhit_Sc[iS][num_ptr_x+1][num_ptr_y+0]==1 &&
	  checkhit_Sc[iS][num_ptr_x+1][num_ptr_y+1]==1 &&
	  checkhit_Sc[iS][num_ptr_x+1][num_ptr_y-1]==1 &&
	  checkhit_Sc[iS][num_ptr_x-1][num_ptr_y+0]==1 &&
	  checkhit_Sc[iS][num_ptr_x-1][num_ptr_y+1]==1 &&
	  checkhit_Sc[iS][num_ptr_x-1][num_ptr_y-1]==1 ) 
	 */
       if(checkhit_Sc[iS][num_ptr_x][num_ptr_y]==1) 
	 {
	   //hitdiag[iS]->Fill(local_z/sxSc,local_y/sySc);
	   int effweight = checkhit_Sc_weight[iS][num_ptr_x][num_ptr_y];

	   //MAKE HIT MAP WITHOUT ERT WEIGHT
	   //hitdiag[iS]->Fill(num_ptr_x,num_ptr_y,(float)effweight);
	   if(E>1.0) hitdiag[iS]->Fill(num_ptr_x,num_ptr_y);

	   if(gn>9){
	     //float sigmax=sqrt(0.155*0.155+(0.57/sqrt(E))*(0.57/sqrt(E)));
	     //float sigmax=sqrt(0.155*0.155+(0.67/sqrt(E))*(0.67/sqrt(E)));
	     float sigmax=sqrt(0.16*0.16+(0.67/sqrt(E))*(0.67/sqrt(E)));

	     //cout << "sigmax " << sigmax << "  E " << E <<endl;
	     //cout << "nonsmeared_localy" << local_y <<endl;
	     //cout << "nonsmeared_localx" << local_x <<endl;

	     local_y  = gRandom->Gaus(local_y,sigmax);
	     local_z  = gRandom->Gaus(local_z,sigmax);
	     local_x  = gRandom->Gaus(local_x,sigmax);

	     //cout << "smeared_localy" << local_y <<endl;
	     //cout << "smeared_localx" << local_x <<endl;

	     smeared_x = local_x*sin(emcRA[iS]*deg);
	     smeared_y = local_y*cos(emcRA[iS]*deg);
	     smeared_z = local_z*-1;

	     //cout << "smeared_midx" << smeared_x <<endl;
	     
	     smeared_x= smeared_x+emcTR1[iS][0]+cos(emcRA[iS]*deg)*ZeroTowerSc[2]-sin(emcRA[iS]*deg)*ZeroTowerSc[1];
	     smeared_y= smeared_y+emcTR1[iS][1]-sin(emcRA[iS]*deg)*ZeroTowerSc[2]-cos(emcRA[iS]*deg)*ZeroTowerSc[1];
	     smeared_z= smeared_z+emcTR1[iS][2]+ZeroTowerSc[0];

	     //cout << "smeared_cross_x " << smeared_x <<endl;
	   }


	   /*  -- p+p stuff
	   int ertbad = check_bad_tile(iS,num_ptr_y,num_ptr_x);
	   if(ertbad==1)hitdiag[iS]->Fill(num_ptr_x,num_ptr_y);
	   return(ertbad);
	   */
	   //right return code?
	   //return(1);
	   return(effweight);
	 }
       else return(0);
   }
   else if(iS<6){
       if(cross_x>0) return(0);
       local_x= cross_x-emcTR1[iS][0]+cos(emcRA[iS]*deg)*ZeroTowerSc[2]+sin(emcRA[iS]*deg)*ZeroTowerSc[1];
       local_y= cross_y-emcTR1[iS][1]+sin(emcRA[iS]*deg)*ZeroTowerSc[2]-cos(emcRA[iS]*deg)*ZeroTowerSc[1];
       local_z= cross_z-emcTR1[iS][2]+ZeroTowerSc[0];
       local_x/= -sin(emcRA[iS]*deg);
       local_y/= cos(emcRA[iS]*deg);
       local_z/= 1;

       //
       // Here comes Fiducial Area Cuts(old version)
       //
       //if(local_z<8 || local_z>389 || local_y<2 || local_y>191) return(0);

       //
       // Get Pointer to the position of calorimeter. Also check if the
       //  pointer points within fiducial area
       //
       if(local_z<0) return(0);
       if(local_y<0) return(0);
       num_ptr_x=(int)(local_z/sxSc); num_ptr_y=(int)(local_y/sySc);
       //if(num_ptr_x>=72||num_ptr_y>=36||num_ptr_x<0||num_ptr_y<0) return(0);
       //use edge cut
       if(num_ptr_x>=70||num_ptr_y>=34||num_ptr_x<2||num_ptr_y<2) return(0);//2004.10.01
       //
       // Check the dead channel map.
       //   "checkhit_Sc[][][]==1" means that it is live channel.
       //

	 /*
       if(checkhit_Sc[iS][num_ptr_x+0][num_ptr_y+0]==1 &&
	  checkhit_Sc[iS][num_ptr_x+0][num_ptr_y+1]==1 &&
	  checkhit_Sc[iS][num_ptr_x+0][num_ptr_y-1]==1 &&
	  checkhit_Sc[iS][num_ptr_x+1][num_ptr_y+0]==1 &&
	  checkhit_Sc[iS][num_ptr_x+1][num_ptr_y+1]==1 &&
	  checkhit_Sc[iS][num_ptr_x+1][num_ptr_y-1]==1 &&
	  checkhit_Sc[iS][num_ptr_x-1][num_ptr_y+0]==1 &&
	  checkhit_Sc[iS][num_ptr_x-1][num_ptr_y+1]==1 &&
	  checkhit_Sc[iS][num_ptr_x-1][num_ptr_y-1]==1 ) 
	 */
       if(checkhit_Sc[iS][num_ptr_x][num_ptr_y]==1) 
	 {
	   //hitdiag[iS]->Fill(local_z/sxSc,local_y/sySc);
	   int effweight=checkhit_Sc_weight[iS][num_ptr_x][num_ptr_y];
	   //Fill without weight
	   //hitdiag[iS]->Fill(num_ptr_x,num_ptr_y,(float) effweight);
	   if(E>1.0) hitdiag[iS]->Fill(num_ptr_x,num_ptr_y);
	   /*
	   int ertbad = check_bad_tile(iS,num_ptr_y,num_ptr_x);
	   if(ertbad==1)hitdiag[iS]->Fill(num_ptr_x,num_ptr_y);
	   return(ertbad);
	   */

	   if(gn>9){
	     //float sigmax=sqrt(0.155*0.155+(0.57/sqrt(E))*(0.57/sqrt(E)));
	     //float sigmax=sqrt(0.155*0.155+(0.67/sqrt(E))*(0.67/sqrt(E)));
	     float sigmax=sqrt(0.16*0.16+(0.67/sqrt(E))*(0.67/sqrt(E)));

	     local_y  = gRandom->Gaus(local_y,sigmax);
	     local_z  = gRandom->Gaus(local_z,sigmax);
	     local_x  = gRandom->Gaus(local_x,sigmax);

	     //cout << "smeared_localx" << local_x <<endl;

	     smeared_z=local_z*1;
	     smeared_y=local_y*cos(emcRA[iS]*deg);
	     smeared_x=local_x*-sin(emcRA[iS]*deg);

	     //cout << "smeared_midx" << smeared_x <<endl;

	     smeared_x= smeared_x+emcTR1[iS][0]-cos(emcRA[iS]*deg)*ZeroTowerSc[2]-sin(emcRA[iS]*deg)*ZeroTowerSc[1];
	     smeared_y= smeared_y+emcTR1[iS][1]-sin(emcRA[iS]*deg)*ZeroTowerSc[2]+cos(emcRA[iS]*deg)*ZeroTowerSc[1];
	     smeared_z= smeared_z+emcTR1[iS][2]-ZeroTowerSc[0];

	     //cout << "smeared_cross_x " << smeared_x <<endl;

	   }

	   return(effweight);
	   //return(1);
	 }
       else return(0);
   }
   else{
     if(cross_x>0) return(0);
       local_x= cross_x-emcTR1[iS][0]+cos(emcRA[iS]*deg)*ZeroTowerGl[2]+sin(emcRA[iS]*deg)*ZeroTowerGl[1];
       local_y= cross_y-emcTR1[iS][1]+sin(emcRA[iS]*deg)*ZeroTowerGl[2]-cos(emcRA[iS]*deg)*ZeroTowerGl[1];
       local_z= cross_z-emcTR1[iS][2]+ZeroTowerGl[0];
       local_x/= -sin(emcRA[iS]*deg);
       local_y/= cos(emcRA[iS]*deg);
       local_z/= 1;
       


       //       cout << iS << ", " << local_y << " " << local_x << endl; 
       //
       // Here comes Fiducial Area Cuts
       //  Temporary not applied to PbGl.
       //       if(local_z<20 || local_z>380 || local_y<10 || local_y>190) return(0);

       //
       // Get Pointer to the position of calorimeter. Also check if the
       //  pointer points within fiducial area
       //
       if(local_z<0) return(0);
       if(local_y<0) return(0);
       num_ptr_x=(int)(local_z/sxGl); num_ptr_y=(int)(local_y/syGl);
       //if(num_ptr_x>=96||num_ptr_y>=48||num_ptr_x<0||num_ptr_y<0) return(0);
       //use edge cut:
       if(num_ptr_x>=94||num_ptr_y>=46||num_ptr_x<2||num_ptr_y<2) return(0);

       //
       // Check the dead channel map.
       //   "checkhit_Gl[][][]==1" means that it is live channel.
       //

	 /*
       if(checkhit_Gl[iS-6][num_ptr_x+0][num_ptr_y+0]==1 &&
	  checkhit_Gl[iS-6][num_ptr_x+0][num_ptr_y+1]==1 &&
	  checkhit_Gl[iS-6][num_ptr_x+0][num_ptr_y-1]==1 &&
	  checkhit_Gl[iS-6][num_ptr_x+1][num_ptr_y+0]==1 &&
	  checkhit_Gl[iS-6][num_ptr_x+1][num_ptr_y+1]==1 &&
	  checkhit_Gl[iS-6][num_ptr_x+1][num_ptr_y-1]==1 &&
	  checkhit_Gl[iS-6][num_ptr_x-1][num_ptr_y+0]==1 &&
	  checkhit_Gl[iS-6][num_ptr_x-1][num_ptr_y+1]==1 &&
	  checkhit_Gl[iS-6][num_ptr_x-1][num_ptr_y-1]==1 ) 
	 */
       //if(checkhit_Gl[iS-6][num_ptr_x][num_ptr_y]==0)cout<<" pbgl hit ? "<<checkhit_Gl[iS-6][num_ptr_x][num_ptr_y]<<endl;
       if(checkhit_Gl[iS-6][num_ptr_x][num_ptr_y]==1) 
	 {	 
	   
	   int effweight  = checkhit_Gl_weight[iS-6][num_ptr_x][num_ptr_y];
	   //hitdiag[iS]->Fill(num_ptr_x,num_ptr_y,(float) effweight);

	   //Fill without ERT weight
	   if(E>1.0) hitdiag[iS]->Fill(num_ptr_x,num_ptr_y);
	   /*
	   int ertbad = check_bad_tile(iS,num_ptr_y,num_ptr_x);
	   if(ertbad==1)hitdiag[iS]->Fill(num_ptr_x,num_ptr_y);
	   return(ertbad);
	   */

	   if(gn>9){
	     //float sigmax=sqrt(0.03*0.03+(0.87/sqrt(E))*(0.87/sqrt(E)));
	     //float sigmax=sqrt(0.1*0.1+(0.87/sqrt(E))*(0.87/sqrt(E)));
	     float sigmax=sqrt(0.16*0.16+(0.67/sqrt(E))*(0.67/sqrt(E)));

	     local_y  = gRandom->Gaus(local_y,sigmax);
	     local_z  = gRandom->Gaus(local_z,sigmax);
	     local_x  = gRandom->Gaus(local_x,sigmax);

	     //cout << "smeared_localx" << local_x <<endl;

	     smeared_z=local_z*1;
	     smeared_y=local_y*cos(emcRA[iS]*deg);
	     smeared_x=local_x*-sin(emcRA[iS]*deg);

	     //cout << "smeared_midx" << smeared_x <<endl;

	     smeared_x= smeared_x+emcTR1[iS][0]-cos(emcRA[iS]*deg)*ZeroTowerGl[2]-sin(emcRA[iS]*deg)*ZeroTowerGl[1];
	     smeared_y= smeared_y+emcTR1[iS][1]-sin(emcRA[iS]*deg)*ZeroTowerGl[2]+cos(emcRA[iS]*deg)*ZeroTowerGl[1];
	     smeared_z= smeared_z-emcTR1[iS][2]-ZeroTowerGl[0];

	     //cout << "smeared_cross_x " << smeared_x <<endl;

	   }


	   return(effweight);
	   //	   return(1);
	 }
       else return(0);
       //else return(0);
       
   }

//   cout << iS << ", " << local_y/local_x << endl; 
//  cout << "Till here4 OK! iS=" << iS << "x= " << num_ptr_x << "y= " << num_ptr_y << endl;
}


float phenix_geom(int gnm,float Px, float Py, float Pz, float vtxz, float E, char *deadfile, float & smeared_x, float & smeared_y, float & smeared_z)
{

  float cross_x,cross_y,cross_z,K;

  if(first_flag==1){Construct_EMCal(deadfile); first_flag=0;}

  posx=Px; posy=Py; posz=Pz;

  if((gnm>=0&&gnm<8)||(gnm>9&&gnm<18)){
    int gm,result;
    gm=gnm;

    if(gnm>9) gm=gnm-10;
    
    // Geometry number 0 is entire W0 and W1

    //
    // Get the crossing point of detector surface plane
    // and particle derection vector 
    //
    K= (A[gm]*emcTR5[gm][0]+B[gm]*emcTR5[gm][1])/(A[gm]*posx+B[gm]*posy);
    //    if(K<0) cout << "output";
    if(K<0) return(0);
    cross_x= K*posx; cross_y= K*posy; cross_z= K*posz+vtxz;
    //cout << "emcalposition: " << cross_x << " " << cross_y << " "<< cross_z << " " << gm << " " << K <<endl;
    //cout << "emcalposition: " << posx << " " << posy << " "<< posz << " " << gm << " " << K <<endl;

    //cout << "emcalposition: " << cross_x << " " << cross_y << " "<< cross_z <<endl;
//    if (gm < 3 || gm==4 || gm==5)
/*
    float smeared_x;
    float smeared_y;
    float smeared_z;
*/

    if (gm < 6)
      { 
	//cout << "calling convert_to_local_and_checkhit " << cross_x << "," << cross_y << "," << cross_z << "," << gm << "," << smeared_x << "," << smeared_y << "," << smeared_z << "," <<  E << "," <<  gnm <<endl;
	result=convert_to_local_and_checkhit(cross_x,cross_y,cross_z,gm,smeared_x,smeared_y,smeared_z, E, gnm);
	//cout << "result " << result <<endl;
	//if(result>0) cout << "cross_x " << cross_x << " smeared_cross_x " << smeared_x << endl;

	float deltax=smeared_x-cross_x;
	float deltay=smeared_y-cross_y;
	float deltaz=smeared_z-cross_z;

	float r=sqrt(cross_x*cross_x+cross_y*cross_y+cross_z*cross_z);

	float deltaposx=deltax*E/r;
	float deltaposy=deltay*E/r;
	float deltaposz=deltaz*E/r;

	float MattE = sqrt((posx+deltaposx)*(posx+deltaposx)+(posy+deltaposy)*(posy+deltaposy)+(posz+deltaposz)*(posz+deltaposz));

	smeared_x=posx+deltaposx;
	smeared_y=posy+deltaposy;
	smeared_z=posz+deltaposz;


	//if(result>0) cout << "pos_x " << posx << " Matt_smeared_pos_x " << posx+deltaposx <<endl;
	//if(result>0) cout << "pos_E " << E << " Matt_smeared_pos_E " << posx+deltaposx <<endl;

	//smeared_x= smeared_x/K; smeared_y= smeared_y/K; smeared_z=(smeared_z-vtxz)/K;
	//if(result>0) cout << "pos_x " << posx << " smeared_pos_x " << smeared_x << endl;


	if(result>0){
	  DELTAx->Fill(deltax);
	  DELTAy->Fill(deltay);
	  DELTAz->Fill(deltaz);
	}
	//if(deltax==0){cout << "XXXXX" << gnm << " deltay=" << deltay << " deltaz=" << deltaz <<endl;}

	//if(result>0) cout << "smeared_x " << smeared_x << endl;

      }
    else 
      {
	//result = 0;
	//cout << "calling convert_to_local_and_checkhit " << cross_x << "," << cross_y << "," << cross_z << "," << gm << "," << smeared_x << "," << smeared_y << "," << smeared_z << "," <<  E << "," <<  gnm <<endl;

	result=convert_to_local_and_checkhit(cross_x,cross_y,cross_z,gm,smeared_x,smeared_y,smeared_z, E, gnm);

	//cout << "result " << result <<endl;

	float deltax=smeared_x-cross_x;
	float deltay=smeared_y-cross_y;
	float deltaz=smeared_z-cross_z;

	float r=sqrt(cross_x*cross_x+cross_y*cross_y+cross_z*cross_z);

	float deltaposx=deltax*E/r;
	float deltaposy=deltay*E/r;
	float deltaposz=deltaz*E/r;

	if(result>0){
	  DELTAx->Fill(deltax);
	  DELTAy->Fill(deltay);
	  DELTAz->Fill(deltaz);
	}

	//if(deltax==0){cout << "XXXXX" << gnm << " deltay=" << deltay << " deltaz=" << deltaz <<endl;}

	//smeared_x= smeared_x/K; smeared_y= smeared_y/K; smeared_z=(smeared_z-vtxz)/K;

	smeared_x=posx+deltaposx;
	smeared_y=posy+deltaposy;
	smeared_z=posz+deltaposz;


      }
    //    if(gnm>=7&&gnm<=8&&result)evaluate_max(gnm-7,cross_x,cross_y,cross_z);
    //    cout << "results_" << gm << "= " << result << endl;
    if(gnm>9&&result>0) return(K);
    else return(result);
  }
  else if(gnm==8){
    //should never get here, what is it for? -man

    int result1,result2;

    //
    // Do the geometry check both for W0 and W1
    //
    K= (A[0]*emcTR5[0][0]+B[0]*emcTR5[0][1])/(A[0]*posx+B[0]*posy);
    if(K>0){
      cross_x= K*posx; cross_y= K*posy; cross_z= K*posz+vtxz;

      //cout << "calling convert_to_local_and_checkhit " << cross_x << "," << cross_y << "," << cross_z << "," << 0 << "," << smeared_x << "," << smeared_y << "," << smeared_z << "," <<  E << "," <<  gnm <<endl;

      result1=convert_to_local_and_checkhit(cross_x,cross_y,cross_z,0,smeared_x,smeared_y,smeared_z,E,gnm);
      //cout << "result " << result1 <<endl;

    }
    else{result1=0;}

//    if(result1)evaluate_max(0,cross_x,cross_y,cross_z);

    K= (A[1]*emcTR5[1][0]+B[1]*emcTR5[1][1])/(A[1]*posx+B[1]*posy);
    if(K>0){
      cross_x= K*posx; cross_y= K*posy; cross_z= K*posz+vtxz;

      //cout << "calling convert_to_local_and_checkhit " << cross_x << "," << cross_y << "," << cross_z << "," << 1 << "," << smeared_x << "," << smeared_y << "," << smeared_z << "," <<  E << "," <<  gnm <<endl;

      result2=convert_to_local_and_checkhit(cross_x,cross_y,cross_z,1,smeared_x,smeared_y,smeared_z,E,gnm);
      //cout << "result " << result2 <<endl;

    }
    else{result2=0;}

//    if(result2)evaluate_max(1,cross_x,cross_y,cross_z);

//    cout << "result1= " << result1 << ", result2= " << result2 << ", OR= " << (result1|result2) << endl;
    return(result1|result2);
  }
  else return(1);

//printf("K=%f,cross_x=%f,cross_y=%f,cross_z=%f\n",K,cross_x,cross_y,cross_z);

}


int evaluate_max(int iS,float cross_x,float cross_y,float cross_z)
{
	phi=180.0/M_PI*atan2(posy,posx);
	theta=atan2(sqrt(posx*posx+posy*posy),posz);

	eta=-log(tan(theta/2));

	flag=0;
	if(phi>phi_max[iS+1]){phi_max[iS+1]=phi; flag=1;}
	if(phi<phi_min[iS+1]){phi_min[iS+1]=phi; flag=1;}
	if(eta>eta_max[iS+1]){eta_max[iS+1]=eta; flag=1;}
	if(eta<eta_min[iS+1]){eta_min[iS+1]=eta; flag=1;}
	if(cross_x<x_min[iS+1]){x_min[iS+1]=cross_x; flag=1;}
	if(cross_x>x_max[iS+1]){x_max[iS+1]=cross_x; flag=1;}
	if(cross_y<y_min[iS+1]){y_min[iS+1]=cross_y; flag=1;}
	if(cross_y>y_max[iS+1]){y_max[iS+1]=cross_y; flag=1;}
	if(cross_z<z_min[iS+1]){z_min[iS+1]=cross_z; flag=1;}
	if(cross_z>z_max[iS+1]){z_max[iS+1]=cross_z; flag=1;}

	if(flag==1){
		printf("phi_max1= %3.4f, phi_min1= %3.4f, eta_max1= %3.4f, eta_min1= %3.4f\n",phi_max[1],phi_min[1],eta_max[1],eta_min[1]);
		printf("phi_max2= %3.4f, phi_min2= %3.4f, eta_max2= %3.4f, eta_min2= %3.4f\n",phi_max[2],phi_min[2],eta_max[2],eta_min[2]);
		printf("x_max1= %3.4f, x_min1= %3.4f, ",x_max[1],x_min[1]);
		printf("x_max2= %3.4f, x_min2= %3.4f\n",x_max[2],x_min[2]);
		printf("y_max1= %3.4f, y_min1= %3.4f, ",y_max[1],y_min[1]);
		printf("y_max2= %3.4f, y_min2= %3.4f\n",y_max[2],y_min[2]);
		printf("z_max1= %3.4f, z_min1= %3.4f, ",z_max[1],z_min[1]);
		printf("z_max2= %3.4f, z_min2= %3.4f\n",z_max[2],z_min[2]);
		printf("posx= %3.4f, posy= %3.4f, posz= %3.4f\n",posx,posy,posz);
		printf("theta= %3.4f, phi= %3.4f\n\n",180.0/M_PI*theta,phi);
	}
        return flag;
}


int check_bad_tile(int sector, int iy, int iz){
  //need to exclude photons which landed on a bad tile
  //some are only ~ 1/2 bad

  if(sector==0){
    if(iy>=0&&iy<=11&&iz>=60&&iz<=72) return -1;
    if(iy>=24&&iy<=35&&iz>=48&&iz<=49&&iy%2==0) return -1;
  }
  if(sector==1){
    if(iy>=24&&iy<=35&&iz>=48&&iz<=59) return -1;
  }
  if(sector==2){
    if(iy>=12&&iy<=23&&iz>=12&&iz<=23) return -1;
    if(iy>=12&&iy<=23&&iz>=36&&iz<=47&&iy%2==0) return -1;
  }
  if(sector==3){
    if(iy>=24&&iy<=35&&iz>=12&&iz<=23) return -1;
    if(iy>=12&&iy<=23&&iz>=48&&iz<=59) return -1;
    if(iy>=0&&iy<=11&&iz>=48&&iz<=59&&iy%2==0) return -1;
  }
  if(sector==4){
    if(iy>=12&&iy<=23&&iz>=12&&iz<=23) return -1;
    if(iy>=0&&iy<=11&&iz>=24&&iz<=35) return -1;
    if(iy>=24&&iy<=35&&iz>=24&&iz<=35) return -1;
    if(iy>=12&&iy<=23&&iz>=36&&iz<=47) return -1;
  }
  if(sector==5){
   if(iy>=0&&iy<=11&&iz>=12&&iz<=23) return -1;
   if(iy>=12&&iy<=23&&iz>=60&&iz<=72) return -1;
  }
  if(sector==7){
   if(iy>=36&&iy<=45&&iz>=24&&iz<=35) return -1;
  }


  return 1;

}
