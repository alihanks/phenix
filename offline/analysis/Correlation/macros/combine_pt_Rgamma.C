#include <fstream>
#include <iostream>
#include <sstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <cmath>

using namespace std;

double StatError(double n1, double sn1, double x1, double sx1, double n2, double sn2, double x2, double sx2, double n3=0, double sn3=0, double x3=0, double sx3=0, double n4=0, double sn4=0, double x4=0, double sx4=0)
{
  //double err = pow(n1/(n1+n2)*sx1,2)+pow(n2/(n1+n2)*sx2,2)+pow((x2-x1)*n2/(n1+n2)/(n1+n2)*sn1,2)+pow((x1-x2)*n1/(n1+n2)/(n1+n2)*sn2,2);
  double err = pow(n1/(n1+n2+n3+n4)*sx1,2) + pow(n2/(n1+n2+n3+n4)*sx2,2) + pow(n3/(n1+n2+n3+n4)*sx3,2) + pow(n4/(n1+n2+n3+n4)*sx4,2)
  + pow(((x2-x1)*n2+(x3-x1)*n3+(x4-x1)*n4)/(n1+n2+n3+n4)/(n1+n2+n3+n4)*sn1,2)
  + pow(((x1-x2)*n1+(x3-x2)*n3+(x4-x2)*n4)/(n1+n2+n3+n4)/(n1+n2+n3+n4)*sn2,2)
  + pow(((x1-x3)*n1+(x2-x3)*n2+(x4-x3)*n4)/(n1+n2+n3+n4)/(n1+n2+n3+n4)*sn3,2)
  + pow(((x1-x4)*n1+(x2-x4)*n2+(x3-x4)*n3)/(n1+n2+n3+n4)/(n1+n2+n3+n4)*sn4,2);
  
  return sqrt(err);
}

double CorrError(double n1, double x1, double sx1, double n2, double x2, double sx2, double n3=0, double x3=0, double sx3=0, double n4=0, double x4=0, double sx4=0)
{
  double Rup = (n1*(x1+sx1)+n2*(x2+sx2)+n3*(x3+sx3)+n4*(x4+sx4))/(n1+n2+n3+n4);
  double Rdown = (n1*(x1-sx1)+n2*(x2-sx2)+n3*(x3-sx3)+n4*(x4-sx4))/(n1+n2+n3+n4);
  double sys = fabs(Rup-Rdown)/2.0;
  //double sys = (n1*sx1 + n2*sx2 + n3*sx3 + n4*sx4)/(n1+n2+n3+n4);
  return sys;
}

void combine_pt_Rgamma(const char* output = "Rgamma_final_test.root", const char* histo_input = "all_histos.root", const char* rgamma_input = "Rgamma_cent_test.root")
{
  TFile* fhistos = new TFile(histo_input,"read");

  const int NCENT = 4;
  const int NPT = 13;
  const int NPT_TRIG = 4;
  float pt_range[NPT+1] = {5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,12.0,14.0,16.0};

  TGraphErrors *grnew[NCENT];
  TGraphErrors *gstatnew[NCENT];
  TGraphErrors *gsysnew[NCENT];
  float meanincpt[NPT];
  float ntrig[NPT];
  
  TCanvas* cfit = new TCanvas("cfit","cfit");
  cfit->Divide(2,2,0.001,0.001);
  TH1F *h = new TH1F("h","h",100,0.0,20.0);
  h->SetMaximum(6.0);
  h->SetMinimum(0.0);
  
  for(int icent = 0; icent<NCENT; icent++){
    std::ostringstream name;
    name << "C" << icent << "_TRIGPT";
    TH1F* TRIGPT_INC= (TH1F*)fhistos->Get(name.str().c_str());
    
    for(int rbin=0; rbin<NPT; rbin++)
    {
      meanincpt[rbin]=0;
      int lowbin = TRIGPT_INC->FindBin(pt_range[rbin]);
      int highbin = TRIGPT_INC->FindBin(pt_range[rbin+1]);

      ntrig[rbin] = TRIGPT_INC->Integral(lowbin,highbin);
      for( int it = lowbin; it <= highbin; it++)
	meanincpt[rbin] += TRIGPT_INC->GetBinContent(it) * TRIGPT_INC->GetBinCenter(it);
      meanincpt[rbin] /= ntrig[rbin];
    }
    
    TFile* myRfile = new TFile(rgamma_input,"read");
    char grname[100];
    sprintf(grname,"Rgamma_stat_%d",icent);
    TGraphErrors* gRgamma = (TGraphErrors*)myRfile->Get(grname);
    sprintf(grname,"Rgamma_sys_%d",icent);
    TGraphErrors* gRsys = (TGraphErrors*)myRfile->Get(grname);

    cfit->cd(icent+1);
    h->Draw();
    TF1 *f3 = new TF1("fp3","pol3",5,16);
    gRgamma->Fit(f3,"I","",5,16);
    gRgamma->Draw("sameP");
    gRsys->Draw("sameP");

    double bincenter, fitratio;
    double newrgamma[NPT];
    double oldrgamma, dum;
    for( int rbin = 0; rbin < NPT; rbin++ )
    {
      gRgamma->GetPoint(rbin+2,dum,oldrgamma);
      cout << rbin+2 << " pt: " << dum << " Rgamma: " << oldrgamma <<endl;
      bincenter = (pt_range[rbin+1]+pt_range[rbin])/2.0;
      cout << "max: " << pt_range[rbin+1] << " min: " << pt_range[rbin] << " bincenter " << meanincpt[rbin] <<endl;
      if ((f3->Eval(bincenter))==0){
	cout << "error with bincenter" <<endl;
	fitratio=0;
      }else if((f3->Eval(meanincpt[rbin]))==0){
	cout << "error with meanpt" <<endl;
	fitratio=0;
      }else{
	fitratio = (f3->Eval(meanincpt[rbin]))/(f3->Eval(bincenter));
      }
      newrgamma[rbin]=fitratio*oldrgamma;
      cout << "new Rgamma: " << newrgamma[rbin] << " fitratio: " << fitratio <<endl;
    }
    
    cout << "filling the Rerr" << endl;
    double Rerr[NPT];
    double Rsys[NPT];
    for (int i=0; i<NPT; i++)
    {
      Rerr[i] = gRgamma->GetErrorY(i+2);
      Rsys[i] = gRsys->GetErrorY(i+2);
      if( Rerr[i] < 0 ) Rerr[i] = Rerr[i-1];
      if( Rsys[i] < 0 ) Rsys[i] = Rsys[i-1];
      cout << i << " " << Rerr[i] <<endl;
    }
    
    cout << "calculating new rgammas" <<endl;
    double Rgamma[NPT_TRIG];
    double Rerrcomb[NPT_TRIG];
    double Rstatcomb[NPT_TRIG];
    double Rsyscomb[NPT_TRIG];
    double pt[NPT_TRIG];
    double pterr[NPT_TRIG];
    //float pt_range_trig[NPT_TRIG+1] = {5.0,7.0,9.0,12.0,15.0};

    for( int ipt = 0; ipt < NPT_TRIG; ipt++ )
    {
      pterr[ipt] = 0;
      if( ipt < 2 ){
	int bin1 = 4*ipt+0; int bin2 = 4*ipt+1; int bin3 = 4*ipt+2; int bin4 = 4*ipt+3;
	Rgamma[ipt] = (ntrig[bin1]*newrgamma[bin1] + ntrig[bin2]*newrgamma[bin2] + ntrig[bin3]*newrgamma[bin3] + ntrig[bin4]*newrgamma[bin4])/(ntrig[bin1]+ntrig[bin2]+ntrig[bin3]+ntrig[bin4]);

	pt[ipt] = (ntrig[bin1]*meanincpt[bin1] + ntrig[bin2]*meanincpt[bin2] + ntrig[bin3]*meanincpt[bin3] + ntrig[bin4]*meanincpt[bin4])/(ntrig[bin1]+ntrig[bin2]+ntrig[bin3]+ntrig[bin4]);

	Rstatcomb[ipt]= StatError(ntrig[bin1],0,newrgamma[bin1],Rerr[bin1],ntrig[bin2],0,newrgamma[bin2],Rerr[bin2],ntrig[bin3],0,newrgamma[bin3],Rerr[bin3],ntrig[bin4],0,newrgamma[bin4],Rerr[bin4]);

	Rsyscomb[ipt] = CorrError(ntrig[bin1],newrgamma[bin1],Rsys[bin1],ntrig[bin2],newrgamma[bin2],Rsys[bin2],ntrig[bin3],newrgamma[bin3],Rsys[bin3],ntrig[bin4],newrgamma[bin4],Rsys[bin4]);
      }
      else if( ipt == 2 ) {
	int bin1 = 4*ipt+0; int bin2 = 4*ipt+1; int bin3 = 4*ipt+2;
	Rgamma[ipt]= (ntrig[bin1]*newrgamma[bin1] + ntrig[bin2]*newrgamma[bin2] + ntrig[bin3]*newrgamma[bin3])/(ntrig[bin1]+ntrig[bin2]+ntrig[bin3]);

	pt[ipt] = (ntrig[bin1]*meanincpt[bin1] + ntrig[bin2]*meanincpt[bin2] + ntrig[bin3]*meanincpt[bin3])/(ntrig[bin1]+ntrig[bin2]+ntrig[bin3]);

	Rstatcomb[ipt]= StatError(ntrig[bin1],0,newrgamma[bin1],Rerr[bin1],ntrig[bin2],0,newrgamma[bin2],Rerr[bin2],ntrig[bin3],0,newrgamma[bin3],Rerr[bin3]);

	Rsyscomb[ipt] = CorrError(ntrig[bin1],newrgamma[bin1],Rsys[bin1],ntrig[bin2],newrgamma[bin2],Rsys[bin2],ntrig[bin3],newrgamma[bin3],Rsys[bin3]);
      }
      else {
	int bin = 4*ipt+0-1;
	Rgamma[ipt] = newrgamma[bin];
	pt[ipt] = meanincpt[bin];
	Rstatcomb[ipt] = Rerr[bin];
	Rsyscomb[ipt] = Rsys[bin];
      }
      Rerrcomb[ipt] = sqrt(pow(Rstatcomb[ipt],2)+pow(Rsyscomb[ipt],2));
    }

    cout << "filling graph" <<endl;
    sprintf(grname,"gr%d",icent);
    grnew[icent] = new TGraphErrors(NPT_TRIG,pt,Rgamma,pterr,Rerrcomb);
    grnew[icent]->SetName(grname); grnew[icent]->SetTitle(grname);
    sprintf(grname,"stat%d",icent);
    gstatnew[icent] = new TGraphErrors(NPT_TRIG,pt,Rgamma,pterr,Rstatcomb);
    gstatnew[icent]->SetName(grname); gstatnew[icent]->SetTitle(grname);
    sprintf(grname,"sys%d",icent);
    gsysnew[icent] = new TGraphErrors(NPT_TRIG,pt,Rgamma,pterr,Rsyscomb);
    gsysnew[icent]->SetName(grname); gsysnew[icent]->SetTitle(grname);
  }

  TFile* fout = new TFile(output,"RECREATE");
  cout << "outfile is opened" <<endl;
  TCanvas* can = new TCanvas("can","can");
  can->Divide(2,2,0.001,0.001);
  for(int i = 0; i < NPT_TRIG; i++) {
    cfit->cd(i+1);
    gstatnew[i]->Draw("sameP");
    gsysnew[i]->SetLineColor(kRed-9);
    gsysnew[i]->SetLineWidth(3);
    gsysnew[i]->Draw("sameP");

    can->cd(i+1);
    h->Draw();
    grnew[i]->Draw("sameP");
    gstatnew[i]->Write();
    gsysnew[i]->Write();
    grnew[i]->Write();
  }
  cfit->Write();
  can->Write();
  h->Write();
  cout << "wrote graph to file" <<endl;
}
