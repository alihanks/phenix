void fit_spectra(){

  int ispi0=0;

  ifstream fin;
  char line[500];
  //if(ispi0==0)fin.open("AN353dir.dat");
  if(ispi0==0)fin.open("ppg042.dat");

  if(ispi0==1) fin.open("pi0spectra.dat");
  char c[20];

  double A[9][30][2]={0.0}; 
  double B[4][30]={0.0}, Be[4][30]={0.0};
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
      if( strncmp( line, "Cent:", 5) == 0 ) {
	cent_index++;

	sscanf( line, "%s%i%i", &c, &cent_low, &cent_high);
	pt_index=0;
	cout<<"centrality bin "<<cent_index<< "  within  "<<cent_low<<" and "<<cent_high<<endl;
	continue;
      }
      //cent_index=0;

      if( ispi0==1 && (sscanf( line, "%f%f%f%f%f%f%f%f%f", &pt, &value, &staterr,&staterr1,&pt_uncorr_err, &err1,&syserr,&pt_corr_err,&err2) !=9)) continue;
      if( ispi0==0 && (sscanf( line, "%f%f%f%f%f%f%f%f", &pt, &value,  &err1, &err2, &staterr,&staterr1,&syserr,&pt_uncorr_err) !=8)) continue;
      //if(ispi0==0 && sscanf( line, "%f%f%f", &pt, &value, &staterr) !=3) continue;

	//cout <<"filling arrays"<< pt << value <<endl;
	x[cent_index][pt_index] = pt;
	A[cent_index][pt_index][0] = value; 
	B[cent_index][pt_index]=value*pt*2*3.14;
	//B[cent_index][pt_index]=value;
	//A[cent_index][pt_index][1] = err1; //used total error
	A[cent_index][pt_index][1] = staterr; //now using stat error only
	Be[cent_index][pt_index]=staterr;

	//if( ispi0==0) Be[cent_index][pt_index]=err1; //total error
	//else Be[cent_index][pt_index]=staterr;
	pt_index++; 
	cout << "B= " << B[cent_index][pt_index] <<endl;

    }

  int cent=0;

  point[cent]=23;
  TGraphErrors *gr[4];
  double xe[5][30] = {0.0};
  gr[cent] = new TGraphErrors(point[cent],x[cent],B[cent],xe[cent],Be[cent]);
  gr[cent]->SetName("gr0"); gr[cent]->SetTitle("gr0");
  gr[cent]->SetMarkerColor(2);
  gr[cent]->SetMarkerStyle(8);

  /*
  point[1]=23;

  gr[1] = new TGraphErrors(point[1],x[1],B[1],xe[1],Be[1]);
  gr[1]->SetName("gr1"); gr[1]->SetTitle("gr1");
  gr[1]->SetMarkerColor(4);
  gr[1]->SetMarkerStyle(8);
  */

  TCanvas *C1=new TCanvas("c1","c1",1);
  //C1->SetLogy();


  gr[cent]->Draw("AP");
  //gr[1]->Draw("P");

  //TF1 *fa = new TF1("fa","[0]/pow(x+[2],[1])",2,15); 
  TF1 *fa = new TF1("fa","([0]*x^(-1*[1])*(1-1/(1+exp((x-3.75)/0.1)))+[2]/((1+x/[3])^[4])*(1/(1+exp((x-3.75)/0.1))))*x",1,15);
  //TF1 *fa = new TF1("fa","[0]*x*exp(9.99*log(1.219/(1.219+x)))",2,15); 
  //fa->SetParameters(1,8,1);
  //fa->SetParameters(76,8,1187,2.74,15.67);
  //75.9661 8.09399 1187.04 2.74135 15.6747 5.21176 29
  //68.1934 8.15446 867.284 2.42962 14.2742 3.8447 28
  fa->SetParameters(68.1934,8.15446,867.284,2.42962,14.2742);
  //fa->Multiply(2.0*3.14);

  //gr[0]->Fit("fa","R");
  
  //TF1 *fb = new TF1("fb","[0]/pow(x+[2],[1])",8,20); 
  //fb->SetParameters(1,6,1);
  //gr[0]->Fit("fb","R");
  TF1 *fb = new TF1("fb","([0]*x^(-1*[1])*(1-1/(1+exp((x-3.75)/0.1)))+[2]/((1+x/[3])^[4])*(1/(1+exp((x-3.75)/0.1))))*x",1,15);
  //TF1 *fb = new TF1("fb","[0]*x^(-1*[1])",1,15);

  //fb->SetParameters(75.9661, 8.09399, 1187.04, 2.74135, 15.6747);
  //fb->SetParameters(76,8,1187,2.74,15.67);
  //fb->SetParameters(1,7,1,34273,97032);
  //params below are what were used in the exodus code for direct photons
  fb->SetParameter(0,1.22);
  fb->SetParameter(1,6.52);
  
  fb->SetParameter(2,1.88);
  fb->SetParameter(3,3.81e+04);
  fb->SetParameter(4,1.02e+05);
  

  fb->SetLineColor(2);
  
  cout << "integral of fb" << fb->Integral(4,15) <<endl;

  if(ispi0==1) fa->Draw("same");

  if(ispi0==0) gr[cent]->Fit("fb","R");
 
  if(ispi0==0) fb->Draw("same");
  
  cout << fb->GetParameter(0) << " " <<
    fb->GetParameter(1) << " " <<
    fb->GetParameter(2) << " " <<
    fb->GetParameter(3) << " " <<
    fb->GetParameter(4) << endl;
 
 /*
 cout << fb->GetParameter(0) << " " <<
    fb->GetParameter(1) << " " << endl;
  */

  //TFile *fsim=new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/takao_all.root");
  //TFile *fsim=new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_cent0_etas.root");
  //TFile *fsim=new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_hieff.root");
  TFile *fsim=new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/meandmatt1.root");


  TH1D *simpi0pt;
  if(ispi0==0)simpi0pt=(TH1D*)fsim->Get("realomegapt"); //direct photons
  else simpi0pt=(TH1D*)fsim->Get("realpi0pt");
  //float simcounts= simpi0pt->GetBinContent(13);
  //float simcounts=simpi0pt->GetEntries()/88.0;
  float simcounts=2.0525e+08;
  cout << "nevents in sim: " << simcounts <<endl;
  //simpi0pt->Scale(1.0/498500000.0);
  //simpi0pt->Scale(1.0/83500000.0);
  simpi0pt->Scale(1.0/simpi0pt->GetBinWidth(1));
  simpi0pt->Scale(1.0/simcounts);
  cout << "integral of sim" << simpi0pt->Integral(41,150) <<endl;
  simpi0pt->SetMarkerColor(4);
  simpi0pt->SetMarkerStyle(4);
  
  simpi0pt->Draw("same");

  //gr[cent]->Draw("P,same");
  //TF1 *fa = new TF1("fa","[0]*x*exp(9.99*log(1.219/(1.219+x)))",1,15); 
  //simpi0pt->Fit("fa","R");


}
