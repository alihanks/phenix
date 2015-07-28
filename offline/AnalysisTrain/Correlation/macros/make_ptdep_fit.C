/*
  Matt's code to generate xi corrections for gamma-jet event mixing.
  Adopted from Justin, Mike M., Sarah.

no pt binning

  8-23-05 - man

  Adapted again 10-27-09 - meg
Removed C0 prefix 8-19-10
*/
using namespace std;

void make_ptdep_fit(int type,int pttbin=0, int ptabin=0)
{
  //type == 0: inclusive photon
  //type == 1: pi0
  //type == 2: dec
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize(0.06,"x");
  gStyle->SetLabelSize(0.06,"y");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"y");
  
  int pttlo, ptthi, ptalo, ptahi;
  //5-7
  if(pttbin==0){
    pttlo=51;
    ptthi=70;
  }
  //7-9
  if(pttbin==1){
    pttlo=71;
    ptthi=90;
  }
  //9-12
  if(pttbin==2){
    pttlo=91;
    ptthi=120;
  }
  //12-15
  if(pttbin==3){
    pttlo=121;
    ptthi=150;
  }
  
  //0.5-1
  if(ptabin==0){
    ptalo=6;
    ptahi=10;
  }
  //1-2
  if(ptabin==1){
    ptalo=11;
    ptahi=20;
  }
  //2-3
  if(ptabin==2){
    ptalo=21;
    ptahi=30;
  }
  //3-5
  if(ptabin==3){
    ptalo=31;
    ptahi=50;
  }
  //5-7
  if(ptabin==4){
    ptalo=51;
    ptahi=70;
  }
  
    /*
  //0.4-1
  if(ptabin==4){
    ptalo=5;
    ptahi=10;
  }
  //0.4-7
  if(ptabin==5){
    ptalo=5;
    ptahi=70;
  }
  
  //0.5-1
  if(ptabin==4){
    ptalo=6;
    ptahi=10;
  }
  if(ptabin==5){
    ptalo=6;
    ptahi=70;
  }
  */
  /*
  //0.5-1
  if(ptabin==0){
    ptalo=11;
    ptahi=20;
  }
  //1-2
  if(ptabin==1){
    ptalo=21;
    ptahi=40;
  }
  //2-3
  if(ptabin==2){
    ptalo=41;
    ptahi=60;
  }
  //3-5
  if(ptabin==3){
    ptalo=61;
    ptahi=100;
  }
  //5-7
  if(ptabin==4){
    ptalo=101;
    ptahi=140;
  }
  */
  TCanvas *c1=new TCanvas("c1","Ntrigg/Nassoc vs Ncoll/Npart");
  c1->Divide(2,2);
 
  ofstream trigfileout;
  ofstream partfileout;

  //trigger - satexp, arctan, satexp, arctan

  TGraphErrors * gTrigsatexpnpart;
  TGraphErrors * gTrigarctannpart;
  TGraphErrors * gTrigsatexpncoll;
  TGraphErrors * gTrigarctanncoll;
  TGraphErrors * gPartsatexpnpart;
  TGraphErrors * gPartarctannpart;
  TGraphErrors * gPartsatexpncoll;
  TGraphErrors * gPartarctanncoll;


  gTrigsatexpnpart=new TGraphErrors(17);
  gTrigarctannpart=new TGraphErrors(17);    
  gTrigsatexpncoll=new TGraphErrors(17);
  gTrigarctanncoll=new TGraphErrors(17);    
  gPartsatexpnpart=new TGraphErrors(17);
  gPartarctannpart=new TGraphErrors(17);    
  gPartsatexpncoll=new TGraphErrors(17);
  gPartarctanncoll=new TGraphErrors(17);    
    
  gTrigsatexpnpart->SetTitle("Sat Exp, Npart");
  gTrigarctannpart->SetTitle("Arc Tan, Npart");
  gTrigsatexpncoll->SetTitle("Sat Exp, Ncoll");
  gTrigarctanncoll->SetTitle("Arc Tan, Ncoll");
  gPartsatexpnpart->SetTitle("Sat Exp, Npart");
  gPartarctannpart->SetTitle("Arc Tan, Npart");
  gPartsatexpncoll->SetTitle("Sat Exp, Ncoll");
  gPartarctanncoll->SetTitle("Arc Tan, Ncoll");


  //https://www.phenix.bnl.gov/phenix/WWW/p/draft/reygers/glauber/tables_auau_200gev.html

  float ncoll[17] = {1065.4, 845.4, 672.4, 532.7, 421.8, 325.6, 251.0, 188.6/*178.6*/, 139.4, 101.3, 72.1, 49.9, 34.4, 22.6, 14.8, 9.9, 4.9};
  float ncolerr[17] = {105.3, 82.1, 66.8, 52.1, 46.8, 32.4, 25.9, 20.6, 15.4, 12.1, 10.5, 9.6, 8.7, 6.6, 5.1, 3.3, 1.2};
  float npart[17] = {351.4, 299.0, 253.9, 215.3, 181.6/*171.6*/, 151.5, 125.7, 102.7, 82.9, 65.9, 51.6, 39.4, 29.8, 21.5, 15.5, 11.3, 6.3};
  float nparterr[17] = {2.9, 3.8, 4.3, 5.3, 5.6, 4.9, 4.9, 4.3, 4.3, 3.4, 3.2, 3.5, 4.1, 3.8, 3.4, 2.6, 1.2};

  // const int src_bin[6] = {0,20,40,60,80,92};
  // const int mult_suffix[5] = {10,20,50,160,440};

  //Open source files
  //TFile *filein = new TFile ("/direct/phenix+hhj/hge/taxi/Run11AuAu200MinBias/4509/data/correlation_taxi4509.root");
  //TFile *filein = new TFile ("/direct/phenix+hhj/hge/taxi/Run11AuAu200MinBias/4702/data/correlation_taxi4702.root");
  //TFile *filein = new TFile ("/direct/phenix+hhj/hge/taxi/Run11AuAu200MinBias/4734/data/correlation_taxi4734.root");
  TFile *filein = new TFile ("/direct/phenix+hhj/hge/taxi/Run11AuAu200MinBias/4755/data/correlation_taxi4755.root");

  char txtfilename[500];
  if(type == 1) sprintf(txtfilename,"trigfits_pi0_%d.txt",pttbin);
  else sprintf(txtfilename,"trigfits_inc_%d.txt",pttbin);
  trigfileout.open(txtfilename);
  sprintf(txtfilename,"partfits_%d.txt",ptabin);
  partfileout.open(txtfilename);

  TH2D *ZVTXCENTTR;
  if(type == 1) ZVTXCENTTR =  new TH2D(*(TH2D*)filein->Get("h2_ptvscent_trig_pi0"));
  else ZVTXCENTTR =  new TH2D(*(TH2D*)filein->Get("h2_ptvscent_trig_inc"));
  TH1D *TRIGGERCENT = new TH1D(*(TH1D*)ZVTXCENTTR->ProjectionY("TRIGGERCENT",pttlo,ptthi));
  
  TH2D *ZVTXCENTPA = new TH2D(*(TH2D*)filein->Get("h2_ptvscent_part"));
  TH1D *PARTNERCENT = new TH1D(*(TH1D*) ZVTXCENTPA->ProjectionY("PARTNERCENT",ptalo,ptahi));
  
  cout<<" trig int "<<TRIGGERCENT->Integral()<<endl;
  cout<<" part int "<<PARTNERCENT->Integral()<<endl;
  
  TRIGGERCENT->Sumw2();
  PARTNERCENT->Sumw2();
  
 //loop over centrality bins
 
  for (int icent = 0; icent < 17; icent++)
    {
      float ntriggers =0;
      float nassocs =0;
      float ntrigerr =0;
      float nassocerr = 0;
      
      int icshi = (icent+1)*5+1;
      if(icent==16) icshi=100;
      
      for(int ics = icent*5+1; ics<icshi; ics++){
	
	ntriggers += TRIGGERCENT->GetBinContent(ics+1); 
	nassocs += PARTNERCENT->GetBinContent(ics+1);
	
	ntrigerr += TRIGGERCENT->GetBinError(ics+1)*TRIGGERCENT->GetBinError(ics+1); 
	nassocerr += PARTNERCENT->GetBinError(ics+1)*PARTNERCENT->GetBinError(ics+1);
      }

     ntrigerr=sqrt(ntriggers);
     nassocerr=sqrt(nassocs);

     gTrigsatexpnpart->SetPoint(icent, npart[icent], ntriggers);
     gTrigarctannpart->SetPoint(icent, npart[icent], ntriggers);
     gTrigsatexpncoll->SetPoint(icent, ncoll[icent], ntriggers);
     gTrigarctanncoll->SetPoint(icent, ncoll[icent], ntriggers);    

     gTrigsatexpnpart->SetPointError(icent, 0., ntrigerr);
     gTrigarctannpart->SetPointError(icent, 0., ntrigerr);
     gTrigsatexpncoll->SetPointError(icent, 0., ntrigerr);       
     gTrigarctanncoll->SetPointError(icent, 0., ntrigerr);       

     cout<<" ntriggers "<<ntriggers<<" associated per trigger "<<nassocs<<endl;
     cout<<" data: "<<npart[icent]<<" "<<nassocs<<" "<<ntriggers<<endl;
     
     gPartsatexpnpart->SetPoint(icent, npart[icent], nassocs);    
     gPartarctannpart->SetPoint(icent, npart[icent], nassocs);    
     gPartsatexpncoll->SetPoint(icent, ncoll[icent], nassocs);	  
     gPartarctanncoll->SetPoint(icent, ncoll[icent], nassocs);	     

     gPartsatexpnpart->SetPointError(icent, 0., nassocerr);
     gPartarctannpart->SetPointError(icent, 0., nassocerr);
     gPartsatexpncoll->SetPointError(icent, 0., nassocerr);
     gPartarctanncoll->SetPointError(icent, 0., nassocerr);
    }
  
  TF1 * f1 = new TF1("fTrigatannc","[0]*atan(-abs([1])*x^[2])", 0, 1200);
  //f1->SetParameter(0,-1e7);
  f1->SetParameter(0,200);
  f1->SetParameter(1,0.02);
  //f1->SetParameter(1,0.01);//pi0 bin 1
  f1->SetParameter(2,0.088);
  //f1->SetParameter(2,1.5);//pi0 bin 1
  cout << "=======> n1 arctan nc vals==========" << endl;
  //if(type == 1)gTrigarctanncoll->Fit("fTrigatannc","","",170,1200);
  if(type == 1)gTrigarctanncoll->Fit("fTrigatannc","","",10,1200);
  else gTrigarctanncoll->Fit("fTrigatannc","","",1,1200);
  float chi2 = f1->GetChisquare();
  float ndf = f1->GetNDF();
  cout<<" chi2 "<<chi2<<" ndf "<<ndf<<" chi2/ndf "<<chi2/ndf<<endl;
  cout << "end----------n1 arctan nc vals" << endl;
 
  c1->cd(1);
  gTrigarctanncoll->GetXaxis()->SetTitle("N_{coll}");
  gTrigarctanncoll->GetYaxis()->SetTitle("N_{trig}");
  //gTrigarctanncoll->GetYaxis()->SetTitleOffset(1.15);
  gTrigarctanncoll->GetYaxis()->SetTitleSize(0.07);
  gTrigarctanncoll->Draw("AP");
  
  trigfileout<<f1->GetParameter(0)<<" "<<f1->GetParameter(1)<<" "<<f1->GetParameter(2)<<"\n";
  cout<<f1->GetParameter(0)<<" "<<f1->GetParameter(1)<<" "<<f1->GetParameter(2)<<endl;;
 
  TF1 * f2 = new TF1("fTrigsexpnc","[0]*(1 - exp(-1.0*abs([1])*x^[2]))", 0, 1200);
  f2->SetLineColor(2);
 
  if(pttbin<2)f2->SetParameter(0,200);
  else f2->SetParameter(0,1e5);
  f2->SetParameter(1,0.02);
  //f2->SetParameter(1,0.05);//pi0 bin 0
  f2->SetParameter(2,0.88);

  cout<<f2->GetParameter(0)<<" "<<f2->GetParameter(1)<<" "<<f2->GetParameter(2)<<endl;
  
  cout << "=======> n1 satexp nc vals==========" << endl;
  //if(type == 1) gTrigsatexpncoll->Fit("fTrigsexpnc","","",170,1200);
  if(type == 1) gTrigsatexpncoll->Fit("fTrigsexpnc","","",1,1200);
  else gTrigsatexpncoll->Fit("fTrigsexpnc","","",1,1200);
  float chi2 = f2->GetChisquare();
  float ndf = f2->GetNDF();
  cout<<" chi2 "<<chi2<<" ndf "<<ndf<<" chi2/ndf "<<chi2/ndf<<endl;
  cout << "end----------n1 satexp nc vals" << endl;
  
  gTrigsatexpncoll->Draw("P");
  
  trigfileout<<f2->GetParameter(0)<<" "<<f2->GetParameter(1)<<" "<<f2->GetParameter(2)<<"\n";
  cout<<f2->GetParameter(0)<<" "<<f2->GetParameter(1)<<" "<<f2->GetParameter(2)<<endl;;
  
  TF1 * f3 = new TF1("fTrigatannp","[0]*atan(-abs([1])*x^[2])", 0, 400);
  if(pttbin<2)f3->SetParameter(0,300);
  //else f3->SetParameter(0,-1e5);
  else f3->SetParameter(0,300);
  
  f3->SetParameter(1,0.0005);
  //f3->SetParameter(1,0.05);//pi0 bin 1
  f3->SetParameter(2,0.7);
  //f3->SetParameter(2,0.02);//pi0 bin 1
  cout << "=======> n1 arctan vals==========" << endl;
  //if(type == 1) gTrigarctannpart->Fit("fTrigatannp","","",100,400);
  if(type == 1) gTrigarctannpart->Fit("fTrigatannp","","",10,400);
  else gTrigarctannpart->Fit("fTrigatannp","","",1,400);
  float chi2 = f3->GetChisquare();
  float ndf = f3->GetNDF();
  cout<<" chi2 "<<chi2<<" ndf "<<ndf<<" chi2/ndf "<<chi2/ndf<<endl;
  cout << "end----------n1 arctan np vals" << endl;
  
  c1->cd(2);
  gTrigarctannpart->GetXaxis()->SetTitle("N_{part}");
  gTrigarctannpart->GetYaxis()->SetTitle("N_{trig}");
  gTrigarctannpart->GetYaxis()->SetTitleSize(0.07);
  //gTrigarctannpart->GetYaxis()->SetTitleOffset(1.15);
  gTrigarctannpart->Draw("AP");
  
  trigfileout<<f3->GetParameter(0)<<" "<<f3->GetParameter(1)<<" "<<f3->GetParameter(2)<<"\n";
  cout<<f3->GetParameter(0)<<" "<<f3->GetParameter(1)<<" "<<f3->GetParameter(2)<<endl;
  
  TF1 * f4 = new TF1("fTrigsexpnp","[0]*(1 - exp(-1.0*abs([1])*x^[2]))", 0, 400);
  f4->SetLineColor(2);
  
  if(pttbin<2)  f4->SetParameter(0,500);
  // else f4->SetParameter(0,1e5);
   else f4->SetParameter(0,500);
  //f4->SetParameter(1,0.00008);
  f4->SetParameter(1,0.00005);
  f4->SetParameter(2,1.77);
  
  cout<<f4->GetParameter(0)<<" "<<f4->GetParameter(1)<<" "<<f4->GetParameter(2)<<endl;;
  
 
  cout << "=======> n1 satexp vals==========" << endl;
  //  if(type == 1) gTrigsatexpnpart->Fit("fTrigsexpnp","","",10,400);
  if(type == 1) gTrigsatexpnpart->Fit("fTrigsexpnp","","",10,400);//pi0 bin 3
  else gTrigsatexpnpart->Fit("fTrigsexpnp","","",1,400);
  float chi2 = f4->GetChisquare();
  float ndf = f4->GetNDF();
  cout<<" chi2 "<<chi2<<" ndf "<<ndf<<" chi2/ndf "<<chi2/ndf<<endl;
  cout << "end----------n1 satexp np vals" << endl;
  
  gTrigsatexpnpart->Draw("P");
  
  trigfileout<<f4->GetParameter(0)<<" "<<f4->GetParameter(1)<<" "<<f4->GetParameter(2)<<"\n";
  cout<<f4->GetParameter(0)<<" "<<f4->GetParameter(1)<<" "<<f4->GetParameter(2)<<endl;;
  
  TF1 * f5 = new TF1("fPartatannc","[0]*atan(-abs([1])*x^[2])", 0, 1200);
  if(ptabin<1){
    f5->SetParameter(0,-1e9);
    f5->SetParameter(1,0.003);
    f5->SetParameter(2,0.9);
  }
  else{
    f5->SetParameter(0,-1e7);
    f5->SetParameter(1,0.151);
    f5->SetParameter(2,0.9);
  }
  
  cout<<f5->GetParameter(0)<<" "<<f5->GetParameter(1)<<" "<<f5->GetParameter(2)<<endl;;
  
  
  cout << "=======> n2 arctan nc vals==========" << endl;
  gPartarctanncoll->Fit("fPartatannc","","",1,1200);
  float chi2 = f5->GetChisquare();
  float ndf = f5->GetNDF();
  cout<<" chi2 "<<chi2<<" ndf "<<ndf<<" chi2/ndf "<<chi2/ndf<<endl;
  cout << "end----------n2 arctan nc vals" << endl;
  
  c1->cd(3);
  gPartarctanncoll->GetXaxis()->SetTitle("N_{coll}");
  gPartarctanncoll->GetYaxis()->SetTitle("N_{asso}");
  gPartarctanncoll->GetYaxis()->SetTitleOffset(1.15);
  gPartarctanncoll->GetYaxis()->SetTitleSize(0.07);
  gPartarctanncoll->Draw("AP");
  
  partfileout<<f5->GetParameter(0)<<" "<<f5->GetParameter(1)<<" "<<f5->GetParameter(2)<<"\n";
  cout<<f5->GetParameter(0)<<" "<<f5->GetParameter(1)<<" "<<f5->GetParameter(2)<<endl;;
  
  TF1 * f6 = new TF1("fPartsexpnc","[0]*(1 - exp(-1.0*abs([1])*x^[2]))", 0, 1200);
  f6->SetLineColor(2);
  
  f6->SetParameter(0,100000.0);//bin 0
  //f6->SetParameter(0,1500.0);
  f6->SetParameter(1,0.000005); //bin 0
  //f6->SetParameter(1,0.05);
  //f6->SetParameter(2,1.1);
  f6->SetParameter(2,0.005);
  
  
  cout << "=======> n2 satexp nc vals==========" << endl;
  gPartsatexpncoll->Fit("fPartsexpnc","","",10,1200);
  float chi2 = f6->GetChisquare();
  float ndf = f6->GetNDF();
  cout<<" chi2 "<<chi2<<" ndf "<<ndf<<" chi2/ndf "<<chi2/ndf<<endl;
  cout << "end----------n2 satexp nc vals" << endl;
  
  gPartsatexpncoll->Draw("P");
  
  partfileout<<f6->GetParameter(0)<<" "<<f6->GetParameter(1)<<" "<<f6->GetParameter(2)<<"\n";
  cout<<f6->GetParameter(0)<<" "<<f6->GetParameter(1)<<" "<<f6->GetParameter(2)<<endl;;
  
  TF1 * f7 = new TF1("fPartatannp","[0]*atan(-abs([1])*x^[2])", 0, 400);
  if(ptabin<1)  f7->SetParameter(0,-1.5e9);
  //else  f7->SetParameter(0,-1e5);
  else  f7->SetParameter(0,-1.0e10);
  
  //if(ptabin<1)  f7->SetParameter(1,3.5e-3);
  //else   f7->SetParameter(1,0.02);
  f7->SetParameter(1,3.5e-3);
  
  if(ptabin<1)  f7->SetParameter(2,0.9);
  else  f7->SetParameter(2,1.39);
  
  cout<<f7->GetParameter(0)<<" "<<f7->GetParameter(1)<<" "<<f7->GetParameter(2)<<endl;  
  
  
  cout << "=======> n2 arctan vals==========" << endl;
  gPartarctannpart->Fit("fPartatannp","","",1,400);
  float chi2 = f7->GetChisquare();
  float ndf = f7->GetNDF();
  cout<<" chi2 "<<chi2<<" ndf "<<ndf<<" chi2/ndf "<<chi2/ndf<<endl;
  cout << "end----------n2 arctan np vals" << endl;
  
  
  c1->cd(4);
  gPartarctannpart->GetXaxis()->SetTitle("N_{part}");
  gPartarctannpart->GetYaxis()->SetTitle("N_{asso}");
  gPartarctannpart->GetYaxis()->SetTitleOffset(1.15);
  gPartarctannpart->GetYaxis()->SetTitleSize(0.07);
  
  gPartarctannpart->Draw("AP");
  
  partfileout<<f7->GetParameter(0)<<" "<<f7->GetParameter(1)<<" "<<f7->GetParameter(2)<<"\n";
  cout<<f7->GetParameter(0)<<" "<<f7->GetParameter(1)<<" "<<f7->GetParameter(2)<<endl;  
 
  TF1 * f8 = new TF1("fPartsexpnp","[0]*(1 - exp(-1.0*abs([1])*x^[2]))", 0, 400);
  f8->SetLineColor(2);

  // if(ptabin<1)f8->SetParameter(0,-2e9.0);
  // else f8->SetParameter(0,1e3.0);
  if(ptabin == 3) f8->SetParameter(0,1e3.0);
  else f8->SetParameter(0,-2e9.0);

  if(ptabin<1)f8->SetParameter(1,0.007);
  else  f8->SetParameter(1,0.3);
  
  //f8->SetParameter(2,1.05);
  f8->SetParameter(2,0.5);
  
  cout<<f8->GetParameter(0)<<" "<<f8->GetParameter(1)<<" "<<f8->GetParameter(2)<<endl;;
  
  
  cout << "=======> n2 satexp vals==========" << endl;
  gPartsatexpnpart->Fit("fPartsexpnp","","",1,400);
  float chi2 = f8->GetChisquare();
  float ndf = f8->GetNDF();
  cout<<" chi2 "<<chi2<<" ndf "<<ndf<<" chi2/ndf "<<chi2/ndf<<endl;
  cout << "end----------n2 satexp np vals" << endl;
  
  gPartsatexpnpart->Draw("P");
  
  partfileout<<f8->GetParameter(0)<<" "<<f8->GetParameter(1)<<" "<<f8->GetParameter(2)<<"\n";
  cout<<f8->GetParameter(0)<<" "<<f8->GetParameter(1)<<" "<<f8->GetParameter(2)<<endl;;
  
  gTrigsatexpnpart->SetMarkerStyle(8);
  gTrigarctannpart->SetMarkerStyle(8);
  gTrigsatexpncoll->SetMarkerStyle(8);
  gTrigarctanncoll->SetMarkerStyle(8);
  gPartsatexpnpart->SetMarkerStyle(8);
  gPartarctannpart->SetMarkerStyle(8);
  gPartsatexpncoll->SetMarkerStyle(8);
  gPartarctanncoll->SetMarkerStyle(8);
  
  TLegend *t=new TLegend(0.6,0.6,0.9,0.9);
  
  t->AddEntry(f1,"Arctangent","l");
  t->AddEntry(f2,"Saturated Exponential","l");
  t->SetTextSize(0.05);
  t->Draw();
  /*
  c2->cd(1);
  gTrigsatexpnpart->Draw("AP");
  gTrigarctannpart->Draw("P");
  c2->cd(2);
  gTrigsatexpncoll->Draw("AP");
  gTrigarctanncoll->Draw("P");
  
  TLegend *t=new TLegend(0.6,0.6,0.9,0.9);
  
  t->AddEntry(f1,"Arctangent","l");
  t->AddEntry(f2,"Saturated Exponential","l");
  
  t->SetBorderSize(0);
  t->SetFillStyle(0);
  t->Draw();
  
  c2->cd(3);
  gPartsatexpnpart->Draw("AP");
  gPartarctannpart->Draw("P");
  c2->cd(4);
  gPartsatexpncoll->Draw("AP");
  gPartarctanncoll->Draw("P");
  */
  
  trigfileout.close();
  partfileout.close();
  /*
  char plotfilename[500];
  if(type == 1) sprintf(plotfilename,"plots/xi_fits_pi0_%d_%d.png",pttbin,ptabin);
  else sprintf(plotfilename,"plots/xi_fits_inc_%d_%d.png",pttbin,ptabin);
  
  c1->SaveAs(plotfilename);
  */
}



 

       
 

					


