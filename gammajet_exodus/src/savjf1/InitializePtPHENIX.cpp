//-----------------------------------------------------------------------------
//
// Generate transverse-momentum distributions for the PHENIX generator
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <math.h>

#include <TH1.h>
#include <TRandom.h>

#include "Momentum.h"
#include "Particle.h"
#include "ParticlePropertyList.h"
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"


TH1F * InitializeRandomHist(TH1F*);
float  InitializePtHydro(double, double, double, float);

TH1F * InitializePtPHENIX(int setup, int ParticleID,
			  double f_c, double f_p0, double f_a, 
			  double f_b, double f_n,			 
			  ParticlePropertyList *PPList)
{
  int    nbins, ibin;
  float  binwidth;
  float  pt, ptmin, ptmax, mt, weight;
  double mass, pimass;
  double t_fo = 0.1;
  double beta = 0.0;

      char title[100];
       ifstream fin("FitResult.txt", ifstream::in );
       ifstream fdirin("ppg042dirfit.dat", ifstream::in );

       Float_t AVal,p0Val,nVal,chi2Val,NDFVal;
       Float_t A1Val,p1Val;
       TF1 *fitfuncPN[10];
       TF1 *dirfit; 
       Float_t dNdyPN[10];
       int centralitybin=0;
       TF1 *f3;
 

  ParticleProperty * PParticle = PPList->GetByID(111);
  pimass = PParticle->GetMass();
  PParticle = PPList->GetByID(ParticleID);
  mass = PParticle->GetMass();

  if ( ParticleID == -111 ) {
    nbins = 30000;
    //ptmin = 0.0;
    ptmin = 0.5;
    ptmax = 30.0;
  }
  else {
    //    nbins = 15000;
    nbins = 30000;
    ptmin = 0.5;
    //ptmin = 3.0;
    //    ptmax = 15.0;
    ptmax = 30.0;
  }

  //  ptmin = 4.5;
  //ptmax = 10.0;
  binwidth = (ptmax-ptmin)/(double)nbins;
  TH1F * pthistogram = new TH1F("pt","pt",nbins,ptmin,ptmax);

  weight = 1.0;
  for ( ibin=1; ibin<=nbins; ibin++ )
  {

    pt = ptmin+(double)(ibin-1)*binwidth+binwidth/2.0;
    mt = sqrt(pt*pt+mass*mass);

    
    //cout<<" pt "<<pt<<" mass "<<mass<<" setup "<<setup<<endl;
    switch(setup)
    {
      case 1:  if ( ibin==1 )
	         cout << "Initializing power law for particle " 
		      << ParticleID << endl; 
      
               weight = pt*exp(9.99*log(1.219/
			(1.219+sqrt(pt*pt+mass*mass-pimass*pimass))));

               break;

      case 2:  if ( ibin==1 )
	         cout << "Initializing exponential with flow for particle " 
		      << ParticleID << endl; 

               if ( ParticleID==111 )
               {
                 t_fo = 0.157;
	         beta = 0.40;
	       }
	       else
	       {
                 t_fo = 0.179;
	         beta = 0.0;
               }
	       t_fo = 0.157;
	       beta = 0.40;
	       weight = pt*exp(-1.0*mt/(t_fo+beta*beta*mass));

               break;

      case 3:  if ( ibin==1 )



	         cout << "Initializing hydrodynamical parameterization " 
		      << "for particle " << ParticleID << endl; 
      
               t_fo = 0.100;
               beta = 0.73;
               weight = InitializePtHydro(mass,t_fo,beta,pt);

               break;

      case 4:  if ( ibin==1 )
	         cout << "AuAu: Initializing most realistic parameterization " 
		      << "for particle " << ParticleID << endl; 

              weight = pt*exp(9.99*log(1.219/
			(1.219+sqrt(pt*pt+mass*mass-pimass*pimass))));


	       weight = pt * f_c / 
		 pow(exp(-f_a*sqrt(mt*mt-pimass*pimass)
			 -f_b*sqrt(mt*mt-pimass*pimass)*
			 sqrt(mt*mt-pimass*pimass)) + 
		     sqrt(mt*mt-pimass*pimass)/f_p0,f_n);

	       // direct photons
	       if ( ParticleID == -111 || 22) {
		 weight = pt/pow(exp(0.1*pt)+pt/0.34,5.89);
	       }

	       // Ke3
	       if ( ParticleID == 21 ) {
		 weight = pt/pow(1.8+pt,17.0);
	       }

	       if(weight==0)
		 cout << "Something is Wrong!" << endl;
	       break;

    case 5:  if ( ibin==1 )
      cout << "pp: Initializing most realistic parameterization " 
	   << "for particle " << ParticleID << endl; 
      
      //from run 3 parameterization 
      //weight = pt * 14.43/pow(pt,8.1028);
      //weight = pt/pow(1.8+pt,17.0);
      weight = pt * 14.43/pow(1.0*pt,8.1);

      //cout<<" pt "<<pt<<" weight "<<weight<<endl;

      /*      
	     weight = pt * f_c / 
	       pow(exp(-f_a*sqrt(mt*mt-pimass*pimass)
		       -f_b*sqrt(mt*mt-pimass*pimass)*
		       sqrt(mt*mt-pimass*pimass)) + 
		   sqrt(mt*mt-pimass*pimass)/f_p0,f_n);

	     // direct photons
	     if ( ParticleID == -111 ) {
	       weight = pt/pow(exp(0.1*pt)+pt/0.34,5.89);
	     }

	     // Ke3
	     if ( ParticleID == 21 ) {
	       weight = pt*exp(-1.0*pt/0.453)/pow(0.408+pt,7.249);
	     }
      */
      break;
      /*
      case 6:  if ( ibin==1 )
	//straight power law -- man

	       weight = f_a / pow(f_b+pt,f_n);

               break;
      */
     case 6:  if ( ibin==1 )
       {
	 cout << "Initializing Takao's parameterization for particle " 
	      << ParticleID << endl;
 
	 //from Takao's MC -- meg

	 for(int i=0;i<9;i++){
	   fin>>AVal>>p0Val>>A1Val>>p1Val>>nVal>>chi2Val>>NDFVal;
	   sprintf(title,"fitfuncPN%d",i);
	   fitfuncPN[i] = new TF1(title,"([0]*x^(-1*[1])*(1-1/(1+exp((x-3.75)/0.1)))+[2]/((1+x/[3])^[4])*(1/(1+exp((x-3.75)/0.1))))*x",ptmin,ptmax);
	   
	   //fitfuncPN[i]->SetLineWidth(0.8);
	   fitfuncPN[i]->SetLineStyle(2);
	   fitfuncPN[i]->SetParameter(0,AVal);
	   fitfuncPN[i]->SetParameter(1,p0Val);
	   fitfuncPN[i]->SetParameter(2,A1Val);
	   fitfuncPN[i]->SetParameter(3,p1Val);
	   fitfuncPN[i]->SetParameter(4,nVal);
	   fitfuncPN[i]->SetLineColor(1);
	   dNdyPN[i]=fitfuncPN[i]->Integral(ptmin,ptmax)*2.0*3.141592653589793238; //change 6 to 2
	   if(ParticleID==111) cout << "Use " << dNdyPN[i] << " for centralitybin " << i << endl;
	   if(ParticleID==111) cout << "Use " << fitfuncPN[i]->Integral(0.5,30.0)*2.0*3.141592653589793238 << " for 0.5-30 " << i << endl;

	   //if(dNdyPN[i]<f_n+1.0 && dNdyPN[i]>f_n-1.0 ) centralitybin=i;

	 }
	 fin.close();              
	 //bins: 0-10%:10-20%:20-30%:30-40%:40-50%:50-60%:70-80%:80-92%:00-92%

	 //}//if ibin==1


       if(ParticleID==-111 || ParticleID==22) {
	 
	 ifstream fin;
	 char line[500];
	 fin.open("TadaAki_R.dat");
	 char c[20];
	 
	 double A[9][30][2]={{{0.0}}}; 
	 double B[9][30]={{0.0}}, Be[9][30]={{0.0}};
	 double x[9][30]={{0.0}};
	 float pt,value, err1, err2,staterr1,staterr2, syserr1,syserr2;
	 float staterr,syserr,pt_uncorr_err,pt_corr_err;
	 double comb_value, comb_err;
	 int cent_low=0, cent_high=0;
	 int cent_index = -1;
	 int pt_index = 0;
	 int point[9] = {0};
	 TGraphErrors *gr[4];
	 
	 
	 
	 while(fin.getline(line,500))
	   {
	     if( !(line && strlen(line)) ) continue;
	     if( strncmp( line,"//",2) == 0 ) continue;
	     if( strncmp( line, "Centrality:", 11) == 0 ) {
	       cent_index++;
	       
	       sscanf( line, "%s%i%i", &c, &cent_low, &cent_high);
	       pt_index=0;
	       //cout<<"centrality bin "<<cent_index<< "  within  "<<cent_low<<" and "<<cent_high<<endl;
	       continue;
	     }
	     if(sscanf( line, "%f%f%f%f%f%f%f", &pt, &value, &staterr,&syserr, &err1,&pt_uncorr_err, &pt_corr_err) !=7) continue;
	     x[cent_index][pt_index] = pt;
	     B[cent_index][pt_index] = value;
	     //A[cent_index][pt_index][1] = err1; //used total error
	     Be[cent_index][pt_index] = staterr; //now using stat error only
	     pt_index++; 
	   }
	  fin.close();
       
	 for(int i = 0; i < 1; i++) {
	   for(int k = 0; k < 30; k++){
	     if(B[i][k]==0) {  point[i] = k; break; }
	   }
	   //cout<<"for centrality "<<i<<"  we have "<<point[i]<<" points."<<endl;
	 }
	 for (int i=0; i<1; i++){
	   //cout << "centrality bin: " << i <<endl;
	   for (int j=0; j<point[i]; j++){
	     cout << j << " " << x[i][j] << " " << B[i][j] << " " << Be[i][j] <<endl;
	   }
	 }
	 
	 char name[512];
	 double xe[5][30] = {{0.0}};
	 
	 gr[0] = new TGraphErrors(point[0],x[0],B[0],xe[0],Be[0]);
	 gr[0]->SetName("gr0"); gr[0]->SetTitle("gr0");
	 
	 //TF1 *f3 = new TF1("fp3","pol3",0.2,14);  
	 f3= new TF1("f3","1.0+[0]*x^[1]",0.2,16);
	 
	 gr[0]->Fit(f3,"I","",1,14);

	 // cout << "param0: " << f3->GetParameter(0) << "param0: " << f3->GetParameter(1) << endl;

	 float totalgammas=0.0;
	 for(int nents=0; nents<100; nents++){
	   float sizebin=(ptmax-1.5)/100;
	   float npt=sizebin*nents+1.5;
	   //float Rg=f3->Integral(npt,npt+sizebin);
	   float Rg=f3->Eval(npt+sizebin/2);
	   float pidN=fitfuncPN[centralitybin]->Integral(npt,npt+sizebin);
	   totalgammas=totalgammas+pidN*(1.0-1.0/Rg);
	   //cout << npt << " to " << npt+sizebin << " rgamma= " << Rg << " Npi0 " << pidN << " gammas= " << totalgammas <<endl;
	 }
	 totalgammas=totalgammas*2.0*3.14; //switch 6 to 2
	 cout << "Use " <<  totalgammas << " for centralitybin 0 " << endl;
	 

	 //float temp=fitfuncPN[0]->Integral(0.2,30)*6.0*3.14;
	 //cout << "Use " << temp*(1.0-1.0/rgamma) << " for centralitybin 0 " << endl;
	 // cout << "Use " <<  f3->Integral(0.2,30)*6.0*3.14 << " for centralitybin 0 " << endl;
	 

	 /*
	 dirfit= new TF1("dirfit","[0]*x^(-1*[1])",ptmin,ptmax);
	 for(int i=0;i<1;i++){
	   fdirin>>AVal>>p0Val>>A1Val>>p1Val>>nVal;
	   //fdirin>>AVal>>p0Val;
	   if(i==centralitybin){
	     dirfit->SetParameter(0,AVal);
	     dirfit->SetParameter(1,p0Val);
	     dirfit->SetParameter(2,A1Val);
	     dirfit->SetParameter(3,p1Val);
	     dirfit->SetParameter(4,nVal);

	     cout << "params used: " << AVal << "  " << p0Val << endl;

	   }
	 }
	 fdirin.close();
	 */

	 dirfit = new TF1("dirfit","([0]*x^(-1*[1])*(1-1/(1+exp((x-3.75)/0.1)))+[2]/((1+x/[3])^[4])*(1/(1+exp((x-3.75)/0.1))))*x",ptmin,ptmax);
	 
	   dirfit->SetParameter(0,1.22);
	   dirfit->SetParameter(1,6.52);
	   dirfit->SetParameter(2,1.88);
	   dirfit->SetParameter(3,3.81e+04);
	   dirfit->SetParameter(4,1.02e+05);
	 
	   float dirdN=dirfit->Integral(ptmin,ptmax)*2.0*3.14;
	   cout << "Use " << dirdN << " from AN353 " << endl;
	   //cout << "Use " << dirdN << " from ppg042 " << endl;
       }
       
       }//if ibin==1

       //cout << "using centrality: " << centralitybin <<endl;

       if(ParticleID==-111 || ParticleID==22) {
	 
	 float rgamma=1.0;
	 rgamma=f3->Eval(pt);
	 
	 //weight=(fitfuncPN[centralitybin]->Eval(pt))*(1.0-1.0/rgamma);
	 //if(weight<0) weight=0;
	 //if(weight>5.0) weight=0;
	 weight=dirfit->Eval(pt);
	 
       } else{ 
	 weight=fitfuncPN[centralitybin]->Eval(pt);
       }


       
       //weight = f_a / pow(f_b+pt,f_n);

       break;

      
      default: if ( ibin==1 ) 
	       {
		 cout << "Error: pt setup " << setup 
		      << " not predefined" << endl;
		 cout << "using flat pt distribution" << endl;
	       }

               weight = 1.0;

               break;
    }

    pthistogram->AddBinContent(ibin,weight);
  }

  pthistogram->Scale(nbins/pthistogram->Integral());
  
  return pthistogram;

}



