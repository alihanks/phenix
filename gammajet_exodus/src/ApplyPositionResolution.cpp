/* Standard includes */
#include <stdio.h>
#include <math.h>
#include <istream>
#include <fstream>

#include <TH2.h>

#define  INCLUDEFLAG extern
#include "DeclareROOTObjects.h"
#include <TH1.h>
#include <TRandom.h>
#include "Momentum.h"
#define PI  3.141592653589793238


using namespace std;

double thetaof(Mom3);
double phiof(Mom3);
float phenix_geom(int,float,float,float,float,char*);

Mom4 ApplyPositionResolution(int gnm, Mom4 mom4, float vtxz, char *deadfile)
{
  double p, px, py, pz, E, theta, phi;
  double sigmax=0.0;
  Mom4   result;
  double P, Px, Py, Pz, E1, theta1, Phi, mass;

  Px = mom4.Getp().Getpx();
  Py = mom4.Getp().Getpy();
  Pz = mom4.Getp().Getpz();
  P  = sqrt(Px*Px+Py*Py+Pz*Pz);
  E  = mom4.GetE();

  theta = thetaof(mom4.Getp());
  Phi   = phiof(mom4.Getp());

  //cout << "vector values " << Px << " " << Py << " " << Pz << " " << P << " " << E << " " << theta << " " << Phi << " " << atan2(Py,Px) <<endl;

  mass = mom4*mom4;
  if ( mass<0.0001 )
  {
    mass = 0.51099906e-3;
  }
  else
  {
    mass = sqrt(mass);
  }

  int gm;
  float accept;
  float cross_x,cross_y,cross_z,K;
  cross_x=0;
  cross_y=0;
  cross_z=0;



  accept=phenix_geom(gnm,Px,Py,Pz,vtxz,deadfile);
  //if(first_flag==1){Construct_EMCal(deadfile); first_flag=0;}

  
  double posx=Px; 
  double posy=Py; 
  double posz=Pz;
  

  if(accept>0){
    K=accept;
    cross_x= K*posx; cross_y= K*posy; cross_z= K*posz+vtxz;

  /*
  if(gnm>=0&&gnm<8){

    // Geometry number 0 is entire W0 and W1

    gm=gnm;
    //
    // Get the crossing point of detector surface plane
    // and particle derection vector 
    //
    K= (A[gm]*emcTR5[gm][0]+B[gm]*emcTR5[gm][1])/(A[gm]*posx+B[gm]*posy);
    //    if(K<0) cout << "output";
    if(K<0) return(mom4);
    //cross_x= K*posx; cross_y= K*posy; cross_z= K*posz+vtxz;

    
    sigmax=1.4+5.9/sqrt(E);
 
    cross_x= K*posx; cross_y= K*posy; cross_z= K*posz+vtxz;


    if (gm < 6)
      { 
	accept=convert_to_local_and_checkhit(cross_x,cross_y,cross_z,gm);
      }
    else 
      {
	//accept = 0;
	accept=convert_to_local_and_checkhit(cross_x,cross_y,cross_z,gm);
      }
    //    if(gnm>=7&&gnm<=8&&accept)evaluate_max(gnm-7,cross_x,cross_y,cross_z);
    //    cout << "accepts_" << gm << "= " << accept << endl;
    //return(accept);
  }
  */  
    float phitemp=Phi; //changed phitemp=phi to =Phi on 12/08/10
    if(phitemp<-PI/2.) phitemp+=2*PI;
    if(phitemp>3*PI/2.) phitemp-=2*PI;

    if(phitemp<2.95) sigmax=sqrt(0.155*0.155+(0.57/sqrt(E))*(0.57/sqrt(E)));
    else sigmax=sqrt(0.03*0.03+(0.87/sqrt(E))*(0.87/sqrt(E)));

    //if(phitemp<2.95) sigmax=0.155+0.57/sqrt(E);
    //else sigmax=0.03+0.87/sqrt(E);

  //if(accept>0){
    //sigmax=1.4+5.9/sqrt(E);//I tried this first in mm
    //maybe I should use cm
    //sigmax=.14+.59/sqrt(E);

    //sigmax=0.155+0.57/sqrt(E); //NIM says this for pbsc
    //sigmax=0.03+0.87/sqrt(E); //NIM says this for pbgl

    //cout << "crossing: " << cross_x << " " << cross_y <<endl;
 
   //cross_x  = gRandom->Gaus(cross_x,sigmax);
    cross_y  = gRandom->Gaus(cross_y,sigmax);
    cross_z  = gRandom->Gaus(cross_z,sigmax);
    //cout << "smeared_crossing: " << cross_x << " " << cross_y <<endl;

    //double PI=3.141592653589793238;
    phi=atan2(cross_y,cross_x);
    //cout << "New_Old_phis: " << phi << " " << Phi <<endl;
    if(phi<0) phi+=2*PI;
    //if(phi>3*PI/2.) phi-=2*PI;
    //cout << "New_Old_phis: " << phi << " " << Phi <<endl;

    //theta=tan(sqrt()/crossz)


    px    = E*sin(theta)*cos(phi);
    py    = E*sin(theta)*sin(phi);
    pz    = E*cos(theta);

    //E     = sqrt(p*p+mass*mass);
    
    //  px    = px*E/p;;
    //  py    = py*E/p;
    //  pz    = pz*E/p;
    
    result.SetE(E);
    result.Setp(px,py,pz);

    return (result);
  }else{
    return(mom4);
  }

}



