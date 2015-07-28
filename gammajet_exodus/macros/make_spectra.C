void make_spectra(){
  ifstream fin("FitResult.txt", ifstream::in );
  cout << "Initializing Takao's parameterization for particle " <<endl;

      char title[100];

       Float_t AVal,p0Val,nVal,chi2Val,NDFVal;
       Float_t A1Val,p1Val;
       TF1 *fitfuncPN[10];
       TF1 *dirfit; 
       Float_t dNdyPN[10];
       int centralitybin=0;
       //TF1 *f3;
       TH1F *dirgammas= new TH1F("dirgammas","dirgammas",300,0,30);
  int    nbins, ibin;
  float  binwidth;
  float  pt, ptmin, ptmax, mt, weight;
  double mass, pimass;

   nbins = 30000;
    ptmin = 0.5;
    //ptmin = 3.0;
    //    ptmax = 15.0;
    ptmax = 30.0;


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
	   fitfuncPN[i]->SetLineColor(2);
	   dNdyPN[i]=fitfuncPN[i]->Integral(0.2,30)*6.0*3.14;
	   cout << "Use " << dNdyPN[i] << " for centralitybin " << i << endl;
	   
	 }
	 fin.close();              
	 //bins: 0-10%:10-20%:20-30%:30-40%:40-50%:50-60%:70-80%:80-92%:00-92%

	 //}//if ibin==1

	 TCanvas * c1 = new TCanvas("c1","c1",1);
	 fitfuncPN[0]->Draw();

	   //if(ParticleID==-111 || ParticleID==22) {

	 ifstream fin;
	 char line[500];
	 fin.open("TadaAki_R.dat");
	 char c[20];
	 
	 double A[9][30][2]={0.0}; 
	 double B[9][30]={0.0}, Be[9][30]={0.0};
	 double x[9][30]={0.0};
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
	 double xe[5][30] = {0.0};
	 
	 gr[0] = new TGraphErrors(point[0],x[0],B[0],xe[0],Be[0]);
	 gr[0]->SetName("gr0"); gr[0]->SetTitle("gr0");
	 
	 //TF1 *f3 = new TF1("fp3","pol3",0.2,14);  
	 TF1 *f3= new TF1("f3","1.0+[0]*x^[1]",0.2,16);
	 
	 gr[0]->Fit(f3,"I","",1,14);

	 cout << "param0: " << f3->GetParameter(0) << "param0: " << f3->GetParameter(1) << endl;

	 float totalgammas=0.0;
	 for(int nents=0; nents<100; nents++){
	   float sizebin=(ptmax)/100;
	   float npt=sizebin*nents;  //+1.5;
	   //float Rg=f3->Integral(npt,npt+sizebin);
	   float Rg=f3->Eval(npt+sizebin/2);
	   if(Rg<1) Rg=1;
	   float pidN=fitfuncPN[centralitybin]->Integral(npt,npt+sizebin);
	   totalgammas=totalgammas+pidN*(1.0-1.0/Rg);
	   cout << npt << " to " << npt+sizebin << " rgamma= " << Rg << " Npi0 " << pidN << " gammas= " << totalgammas <<endl;

	   //loop over all higher pTs that may contribute to the decay
	   float Ndec=0.0;
	   for(int j=nents+1; j<101; j++){
	     float jpt=sizebin*j;
	     Ndec=Ndec+sizebin/(sizebin*j)*(fitfuncPN[centralitybin]->Integral(jpt,jpt+sizebin));

	   }

	  
	   dirgammas->SetBinContent(nents+1,Ndec*(Rg-1));

	   //dirgammas->SetBinContent(nents+1,pidN*(1.0-1.0/Rg));
	   //dirgammas->Fill(pidN*(1.0-1.0/Rg));
	 }
	 TCanvas * c2 = new TCanvas("c2","c2",1);
	 dirgammas->Draw();

	 cout << "Use " << dirgammas->Integral()  << " from histo " << endl;

	 totalgammas=totalgammas*6.0*3.14;
	 cout << "Use " <<  totalgammas << " for centralitybin 0 " << endl;


	 //float temp=fitfuncPN[0]->Integral(0.2,30)*6.0*3.14;
	 //cout << "Use " << temp*(1.0-1.0/rgamma) << " for centralitybin 0 " << endl;
	 // cout << "Use " <<  f3->Integral(0.2,30)*6.0*3.14 << " for centralitybin 0 " << endl;

       
       /* 
       }//if ibin==1

       if(ParticleID==-111 || ParticleID==22) {
	 
	 float rgamma=1.0;
	 rgamma=f3->Eval(pt);
	 
	 weight=(fitfuncPN[centralitybin]->Eval(pt))*(1.0-1.0/rgamma);
	 if(weight<0) weight=0;
	 if(weight>5.0) weight=0;
	 
       } else{ 
	 weight=fitfuncPN[centralitybin]->Eval(pt);
       }

       */



	   dirfit = new TF1("dirfit","([0]*x^(-1*[1])*(1-1/(1+exp((x-3.75)/0.1)))+[2]/((1+x/[3])^[4])*(1/(1+exp((x-3.75)/0.1))))*x",ptmin,ptmax);

	   dirfit->SetParameter(0,1.22);
	   dirfit->SetParameter(1,6.52);
	   dirfit->SetParameter(2,1.88);
	   dirfit->SetParameter(3,3.81e+04);
	   dirfit->SetParameter(4,1.02e+05);

	   float dirdN=dirfit->Integral(0.2,30)*6.0*3.14;

	   cout << "Use " << dirdN << " from AN353 " << endl;




}
