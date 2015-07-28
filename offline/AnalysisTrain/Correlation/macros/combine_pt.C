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

double GetNTriggers(TH1D* trigpt, double trigptmin, double trigptmax)
{
  double ntrig = 0;
  for(int i = 0; i < trigpt->GetNbinsX(); i++)
  {
    double bincenter = trigpt->GetBinCenter(i+1);
    if(bincenter >=trigptmin && bincenter <= trigptmax)
    {
      double bincontent = trigpt->GetBinContent(i+1);
      ntrig += bincontent;
    }
  }
  return ntrig;
}

void CombinePtHistos(double itrig, double ntrig1, double ntrig2, TH1D* h1, TH1D* h2, TH1D* comb)
{
  double y[4][5][30];//[NPT_TRIG][NPT_PART][NDPHI]
  double yerr[4][5][30];
  double y_comb[2][5][30];
  double yerr_comb[2][5][30];
  for(int ipart=0; ipart<5; ipart++){
    for(int idphi=0; idphi<30; idphi++){
      y[itrig][ipart][idphi] = h1->GetBinContent(idphi+1);
      y[itrig+1][ipart][idphi] = h2->GetBinContent(idphi+1);
      cout<<"Y["<<itrig<<"]["<<ipart<<"]["<<idphi<<"] = "<<y[itrig][ipart][idphi]<<endl;
      yerr[itrig][ipart][idphi] = h1->GetBinError(idphi+1);
      yerr[itrig+1][ipart][idphi] = h2->GetBinError(idphi+1);
      cout<<"YERR["<<itrig<<"]["<<ipart<<"]["<<idphi<<"] = "<<yerr[itrig][ipart][idphi]<<endl;

      int it;
      if(itrig<2) it = 0;
      else it = 1;
      y_comb[it][ipart][idphi] = (ntrig1*y[itrig][ipart][idphi]+ntrig2*y[itrig+1][ipart][idphi])/(ntrig1+ntrig2);
      yerr_comb[it][ipart][idphi] = StatError(ntrig1,0,y[itrig][ipart][idphi],yerr[itrig][ipart][idphi],ntrig2,0,y[itrig+1][ipart][idphi],yerr[itrig+1][ipart][idphi]);
      comb->SetBinContent(idphi+1,y_comb[it][ipart][idphi]);
      comb->SetBinError(idphi+1,yerr_comb[it][ipart][idphi]);
    }
  }
  cout<<"leaving CombinePtHistos function."<<endl;
}

void GetTrigPartBin(string file_string, int& itrig, int& ipart)
{
  int p1 = file_string.find_last_of("/");
  cout<<"p1 = "<<p1<<endl;
  string  photon_bin_string = file_string.substr(p1+6,1);
  int photon_bin = atoi(photon_bin_string.c_str());
  int p2 = file_string.find_last_of(".");
  cout<<"p2 = "<<p2<<endl;
  string hadron_bin_string = file_string.substr(p2-1,1);
  int hadron_bin = atoi(hadron_bin_string.c_str());
  cout<<"photon_bin: "<<photon_bin<<" hadron_bin: "<<hadron_bin<<endl;
  switch(photon_bin){
  case 5:
    itrig = 0;
    break;
  case 7:
    itrig = 1;
    break;
  case 9:
    itrig = 2;
    break;
  case 1:
    itrig = 3;
    break;
  default:
    cout << "no such photon bin!" << endl;
  }

  switch(hadron_bin){
  case 1:
    ipart = 0;
    break;
  case 2:
    ipart = 1;
    break;
  case 3:
    ipart = 2;
    break;
  case 5:
    ipart = 3;
    break;
  case 0:
    ipart = 4;
    break;
  default:
    cout << "no such hadron bin!" << endl;
  }
  cout<<"leaving GetTrigPartBin function..."<<endl;
}

double StatError(double n1, double sn1, double x1, double sx1, double n2, double sn2, double x2, double sx2)
{
  double sum = n1+n2;
  double err = pow(n1/sum*sx1,2) + pow(n2/sum*sx2,2) + pow(((x2-x1)*n2)/sum/sum*sn1,2) + pow(((x1-x2)*n1)/sum/sum*sn2,2);
  
  return sqrt(err);
}

void Combine1DHistos(TH1D* h1, TH1D* h2, TH1D* comb)
{
  for( int ix = 1; ix <= comb->GetNbinsX(); ix++ ){
    double n1 = h1->GetBinContent(ix);
    double n2 = h2->GetBinContent(ix);
    cout<<"n1 = "<<n1<<"; n2 = "<<n2<<endl;
    if( n1!=0 && n2!=0 ){
      cout<<"n1&&n2 case"<<endl;
      double w1 = 1/pow(h1->GetBinError(ix),2);
      double w2 = 1/pow(h2->GetBinError(ix),2);
      double sumw = w1+w2;
      cout<<"w1 = "<<w1<<"; w2 = "<<w2<<"; sumw = "<< sumw<<endl;
      double ntot = w1*n1/sumw + w2*n2/sumw;
      double errtot = 1/sqrt(sumw);
      cout<<"ntot = "<<ntot<<"; errtot = "<<errtot<<endl;

      comb->SetBinContent(ix,ntot);
      comb->SetBinError(ix,errtot);
    }
    else if(n1){
      cout<<"n1 case"<<endl;
      comb->SetBinContent(ix,n1);
      comb->SetBinError(ix,h1->GetBinError(ix));
    }
    else if(n2){
      cout<<"n2 case"<<endl;
      comb->SetBinContent(ix,n2);
      comb->SetBinError(ix,h2->GetBinError(ix));
    }
    else {
      comb->SetBinContent(ix,0);
      comb->SetBinError(ix,0);
    }
  }
}

void combine_pt(int type, const string filein, const string filetaxi, const string fileout)
{
  const int NCENT = 4;
  const int NPT_TRIG = 4;
  const int NPT_PART = 5;
  const int NDPHI = 30;
  const double PI = TMath::Pi();

  ostringstream bin;
  string name;
  TH1D* h1_trig_pt[NCENT];
  //TH1D* h1_part_pt[NCENT];
  TH1D* h1_JF[NCENT][NPT_TRIG][NPT_PART];
  TH1D* h1_JF_COMB[NCENT][2][NPT_PART];
  TH1D* h1_JF_COMBTOT[2][2][NPT_PART];
  int Ntrigs[NCENT][NPT_TRIG];
  double trig_pt_range[NPT_TRIG+1] = {5.0,7.0,9.0,12.0,15.0};
  double part_pt_range[NPT_PART+1] = {0.5,1.0,2.0,3.0,5.0,7.0};
  double Y[NCENT][NPT_TRIG][NPT_PART][NDPHI];
  double YERR[NCENT][NPT_TRIG][NPT_PART][NDPHI];
  double Y_COMB[NCENT][2][NPT_PART][NDPHI];//[5-9,9-12]
  double YERR_COMB[NCENT][2][NPT_PART][NDPHI];

  TFile* fin = new TFile(filein.c_str());
  TFile* ftaxi = new TFile(filetaxi.c_str());
  TFile* fout = new TFile(fileout.c_str(),"recreate");
  fout->cd();

  for(int ic=0; ic<NCENT; ic++){
    bin.str("");
    bin << ic;
    if(type == 0) name = "h1_trig_pt_inc_c" + bin.str();
    if(type == 1) name = "h1_trig_pt_pi0_c" + bin.str();
    if(type == 2) name = "h1_trig_pt_dec_c" + bin.str();
    h1_trig_pt[ic] = new TH1D(*(TH1D*)ftaxi->Get(name.c_str()));   
    // name = "h1_part_pt_c" + bin.str();
    // h1_trig_pt[ic] = new TH1D(*(TH1D*)filein->Get(name.c_str()));

    for(int itrig=0; itrig<NPT_TRIG; itrig++){
      int binmin = h1_trig_pt[ic]->FindBin(trig_pt_range[itrig]);
      int binmax = h1_trig_pt[ic]->FindBin(trig_pt_range[itrig+1]);
      if(type == 0 || type == 1) {
	Ntrigs[ic][itrig] = h1_trig_pt[ic]->Integral(binmin, binmax-1);
	//cout<<"Ntrigs = "<<Ntrigs[ic][itrig]<<" = "<<GetNTriggers(h1_trig_pt[ic], trig_pt_range[itrig], trig_pt_range[itrig+1])<<endl;
      }
      if(type == 2) {
	if(itrig<3){
	  Ntrigs[ic][itrig] = h1_trig_pt[ic]->GetBinContent(itrig+1);
	}
	else{
	  Ntrigs[ic][itrig] = h1_trig_pt[ic]->GetBinContent(itrig+2);
	}
      }
      for(int ipart=0; ipart<NPT_PART; ipart++){
	bin.str("");	
	bin <<ic<<"_p"<<itrig<<"_h"<<ipart;
	name = "JF_c" + bin.str();
	cout<<"name = "<<name.c_str()<<endl;
	h1_JF[ic][itrig][ipart] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
	for(int idphi=0; idphi<NDPHI; idphi++){
	  Y[ic][itrig][ipart][idphi] = h1_JF[ic][itrig][ipart]->GetBinContent(idphi+1);
	  cout<<"Y["<<ic<<"]["<<itrig<<"]["<<ipart<<"]["<<idphi<<"] = "<<Y[ic][itrig][ipart][idphi]<<endl;
	  YERR[ic][itrig][ipart][idphi] = h1_JF[ic][itrig][ipart]->GetBinError(idphi+1);
	  cout<<"YERR["<<ic<<"]["<<itrig<<"]["<<ipart<<"]["<<idphi<<"] = "<<YERR[ic][itrig][ipart][idphi]<<endl;
	}
      }
    }
  }
  cout<<"part 1 done"<<endl;
  for(int ic=0; ic<NCENT; ic++){
    for(int itrig=0; itrig<NPT_TRIG; itrig+=2){
      for(int ipart=0; ipart<NPT_PART; ipart++){
	for(int idphi=0; idphi<NDPHI; idphi++){
	  if(itrig<2) {
	    cout<<"itrig = "<<itrig<<endl;
	    // cout<<"Y["<<ic<<"]["<<itrig<<"]["<<ipart<<"]["<<idphi<<"] = "<<Y[ic][itrig][ipart][idphi]<<endl;
	    // cout<<"Y["<<ic<<"]["<<itrig+1<<"]["<<ipart<<"]["<<idphi<<"] = "<<Y[ic][itrig+1][ipart][idphi]<<endl;
	    // cout<<"YERR["<<ic<<"]["<<itrig<<"]["<<ipart<<"]["<<idphi<<"] = "<<YERR[ic][itrig][ipart][idphi]<<endl;
	    // cout<<"YERR["<<ic<<"]["<<itrig+1<<"]["<<ipart<<"]["<<idphi<<"] = "<<YERR[ic][itrig+1][ipart][idphi]<<endl;
	  Y_COMB[ic][0][ipart][idphi] = (Ntrigs[ic][itrig]*Y[ic][itrig][ipart][idphi]+Ntrigs[ic][itrig+1]*Y[ic][itrig+1][ipart][idphi])/(Ntrigs[ic][itrig]+Ntrigs[ic][itrig+1]);
	  YERR_COMB[ic][0][ipart][idphi] = StatError(Ntrigs[ic][itrig],0,Y[ic][itrig][ipart][idphi],YERR[ic][itrig][ipart][idphi],Ntrigs[ic][itrig+1],0,Y[ic][itrig+1][ipart][idphi],YERR[ic][itrig+1][ipart][idphi]);
	  }
	  else {
	    cout<<"itrig = "<<itrig<<endl;
	    //cout<<"Y["<<ic<<"]["<<itrig<<"]["<<ipart<<"]["<<idphi<<"] = "<<Y[ic][itrig][ipart][idphi]<<endl;
	    //cout<<"YERR["<<ic<<"]["<<itrig<<"]["<<ipart<<"]["<<idphi<<"] = "<<YERR[ic][itrig][ipart][idphi]<<endl;
	    Y_COMB[ic][1][ipart][idphi] = (Ntrigs[ic][itrig]*Y[ic][itrig][ipart][idphi]+Ntrigs[ic][itrig+1]*Y[ic][itrig+1][ipart][idphi])/(Ntrigs[ic][itrig]+Ntrigs[ic][itrig+1]);
	    YERR_COMB[ic][1][ipart][idphi] = StatError(Ntrigs[ic][itrig],0,Y[ic][itrig][ipart][idphi],YERR[ic][itrig][ipart][idphi],Ntrigs[ic][itrig+1],0,Y[ic][itrig+1][ipart][idphi],YERR[ic][itrig+1][ipart][idphi]);
	  }
	  // cout<<"Y_COMB["<<ic<<"][0]["<<ipart<<"]["<<idphi<<"] = "<<Y_COMB[ic][0][ipart][idphi]<<endl;
	  // cout<<"YERR_COMB["<<ic<<"][0]["<<ipart<<"]["<<idphi<<"] = "<<YERR_COMB[ic][0][ipart][idphi]<<endl;
	  // cout<<"Y_COMB["<<ic<<"][1]["<<ipart<<"]["<<idphi<<"] = "<<Y_COMB[ic][1][ipart][idphi]<<endl;
	  // cout<<"YERR_COMB["<<ic<<"][1]["<<ipart<<"]["<<idphi<<"] = "<<YERR_COMB[ic][1][ipart][idphi]<<endl;

	}
      }
    }
  }

  for(int ic=0; ic<NCENT; ic++){
    for(int itrig=0; itrig<2; itrig++){
      for(int ipart=0; ipart<NPT_PART; ipart++){
	bin.str("");
	bin << ic<<"_p"<<itrig<<"_h"<<ipart;
	name = "JF_COMB_c" + bin.str();
	cout<<"name = "<<name.c_str()<<endl;
	h1_JF_COMB[ic][itrig][ipart] = new TH1D(name.c_str(),name.c_str(),30,0.0,PI);
	for(int idphi=0; idphi<NDPHI; idphi++){
	  h1_JF_COMB[ic][itrig][ipart]->SetBinContent(idphi+1,Y_COMB[ic][itrig][ipart][idphi]);
	  h1_JF_COMB[ic][itrig][ipart]->SetBinError(idphi+1,YERR_COMB[ic][itrig][ipart][idphi]);
	}
	
	h1_JF_COMB[ic][itrig][ipart]->Write();
      }
    }
  }

  for(int ic=0; ic<NCENT; ic+=2){
    for(int itrig=0; itrig<2; itrig++){
      for(int ipart=0; ipart<NPT_PART; ipart++){
	bin.str("");
	if(ic<2){
	  bin << ic << "_p"<<itrig<<"_h"<<ipart;
	  name = "JF_COMBTOT_c" + bin.str();
	  cout<<"name = "<<name.c_str()<<endl;
	  h1_JF_COMBTOT[0][itrig][ipart] = new TH1D(name.c_str(),name.c_str(),30,0.0,PI);
	  //cout<<"t1: "<<h1_JF_COMB[ic][itrig][ipart]->GetName()<<"; t2: "<<h1_JF_COMB[ic+1][itrig][ipart]->GetName()<<endl;
	  Combine1DHistos(h1_JF_COMB[ic][itrig][ipart],h1_JF_COMB[ic+1][itrig][ipart],h1_JF_COMBTOT[0][itrig][ipart]);
	  h1_JF_COMBTOT[0][itrig][ipart]->Write();
	}
	else {
	  bin << ic-1 << "_p"<<itrig<<"_h"<<ipart;
	  name = "JF_COMBTOT_c" + bin.str();
	  cout<<"name = "<<name.c_str()<<endl;
	  h1_JF_COMBTOT[1][itrig][ipart] = new TH1D(name.c_str(),name.c_str(),30,0.0,PI);
	  Combine1DHistos(h1_JF_COMB[ic][itrig][ipart],h1_JF_COMB[ic+1][itrig][ipart],h1_JF_COMBTOT[1][itrig][ipart]);
	  h1_JF_COMBTOT[1][itrig][ipart]->Write();
	}
      }
    }
  }
  
  cout<<"finish writing out."<<endl;

}

//(0,"run10_output_020.lst","run10_output_2040.lst","/phenix/hhj/ahanks/run10AuAu/combinedsimple_run10AuAu_020_merged_taxi313_total.root","/phenix/hhj/ahanks/run10AuAu/combinedsimple_run10AuAu_2040_merged_taxi313_total.root","fileout.root")
void combine_pt_run10(int type, const string flist020, const string flist2040, const string ftaxi020, const string ftaxi2040, const string fileout)
{
  const int NCENT = 2;//only using 020 and 2040 files
  const int NPT_TRIG = 4;
  const int NPT_PART = 5;
  const int NDPHI = 30;
  const double PI = TMath::Pi();

  ostringstream bin;
  string name;
  string file_string;
  double Ntrigs[NCENT][NPT_TRIG];
  TH1D* h1_trig_pt[NCENT];
  TH1D* h1_JF[NCENT][NPT_TRIG][NPT_PART];
  TH1D* h1_JF_COMB[NCENT][2][NPT_PART];
  TH1D* h1_JF_COMBTOT[2][NPT_PART];//no centrality bin, 2 trig pt bins

  double trig_pt_range[NPT_TRIG+1] = {5.0,7.0,9.0,12.0,15.0};

  TFile* fin020 = new TFile(ftaxi020.c_str());
  TFile* fin2040 = new TFile(ftaxi2040.c_str());

  TFile* fout = new TFile(fileout.c_str(),"recreate");

  vector<string> Run10FileList020;
  char buffer[256];
  cout<<"opening run10 020 input file list: " << flist020.c_str() << endl;
  ifstream flst020(flist020.c_str());
  if( !flst020.is_open()) return -1;

  while (flst020.getline(buffer, 256, '\n')){
    cout<<"in the while loop"<<endl;
    char this_file[512];
    strcpy(this_file, buffer);
    cout<<"this_file: "<<this_file<<endl;
    string s (this_file);
    cout<<"s: "<<s.c_str()<<endl;
    Run10FileList020.push_back(s);
  }
  const int Nfiles_020 = Run10FileList020.size();
  cout<<"Nfiles_020 = "<<Nfiles_020<<endl;

  cout<<"start looping run10 020 files."<<endl;
  for(int ifile = 0; ifile < Nfiles_020; ifile++){
    cout<<"Run10FileList020 ["<<ifile<<"] = "<<Run10FileList020[ifile]<<endl;
    file_string = Run10FileList020[ifile];

    TFile* file = new TFile(file_string.c_str());
    if (! file->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }

    int cbin = 0;
    int trigbin, partbin;    
    GetTrigPartBin(file_string,trigbin,partbin);
    if(type == 0) h1_JF[cbin][trigbin][partbin] = new TH1D(*(TH1D*)file->Get("JF_INC"));
    if(type == 1) h1_JF[cbin][trigbin][partbin] = new TH1D(*(TH1D*)file->Get("JF_PI0"));
    if(type == 2) h1_JF[cbin][trigbin][partbin] = new TH1D(*(TH1D*)file->Get("JF_DEC"));
    if(type == 3) h1_JF[cbin][trigbin][partbin] = new TH1D(*(TH1D*)file->Get("JF_DIR"));
    bin.str("");
    bin << cbin << "_p"<<trigbin<<"_h"<<partbin;
    name = "JFrun7_c"+bin.str();
    h1_JF[cbin][trigbin][partbin]->SetName(name.c_str());
  }
  cout<<"part 1 done."<<endl;

  if(type == 0) name = "C0_TRIGPT";
  if(type == 1) name = "piC0_TRIGPT";
  if(type == 2) name = "piC0_MWTRIGPT";
  h1_trig_pt[0] = new TH1D(*(TH1D*)fin020->Get(name.c_str()));   
  
  for(int itrig=0; itrig<NPT_TRIG; itrig++){
    int binmin = h1_trig_pt[0]->FindBin(trig_pt_range[itrig]);
    int binmax = h1_trig_pt[0]->FindBin(trig_pt_range[itrig+1]);
    if(type == 0 || type == 1) {
      Ntrigs[0][itrig] = h1_trig_pt[0]->Integral(binmin, binmax-1);
      //cout<<"Ntrigs = "<<Ntrigs[ic][itrig]<<" = "<<GetNTriggers(h1_trig_pt[ic], trig_pt_range[itrig], trig_pt_range[itrig+1])<<endl;
    }
    if(type == 2) {
      if(itrig<3){
	Ntrigs[0][itrig] = h1_trig_pt[0]->GetBinContent(itrig+1);
      }
      else{
	Ntrigs[0][itrig] = h1_trig_pt[0]->GetBinContent(itrig+2);
      }
    }
  }
  cout<<"Ntrigs calculated."<<endl;
  for(int ic=0; ic<NCENT; ic++){
    for(int itrig=0; itrig<2; itrig++){
      for(int ipart=0; ipart<NPT_PART; ipart++){
	bin.str("");
	bin << ic <<"_p"<<itrig<<"_h"<<ipart;
	name = "JFrun7_COMB_c" + bin.str();
	h1_JF_COMB[ic][itrig][ipart] = new TH1D(name.c_str(),name.c_str(),30,0.0,PI);
      }
    }
  }
  cout<<"h1_JF_COMB defined."<<endl;
  for(int ipart=0; ipart<NPT_PART; ipart++){
    CombinePtHistos(0,Ntrigs[0][0],Ntrigs[0][1],h1_JF[0][0][ipart],h1_JF[0][1][ipart],h1_JF_COMB[0][0][ipart]);
    CombinePtHistos(2,Ntrigs[0][2],Ntrigs[0][3],h1_JF[0][2][ipart],h1_JF[0][3][ipart],h1_JF_COMB[0][1][ipart]);
    fout->cd();
    h1_JF_COMB[0][0][ipart]->Write();
    h1_JF_COMB[0][1][ipart]->Write();
  }


  //dealing with run10 2040 files
  vector<string> Run10FileList2040;
  char buffer[256];
  cout<<"opening run10 2040 input file list: " << flist2040.c_str() << endl;
  ifstream flst2040(flist2040.c_str());
  if( !flst2040.is_open()) return -1;

  while (flst2040.getline(buffer, 256, '\n')){
    cout<<"in the while loop"<<endl;
    char this_file[512];
    strcpy(this_file, buffer);
    cout<<"this_file: "<<this_file<<endl;
    string s (this_file);
    cout<<"s: "<<s.c_str()<<endl;
    Run10FileList2040.push_back(s);
  }
  const int Nfiles_2040 = Run10FileList2040.size();
  cout<<"Nfiles_2040 = "<<Nfiles_2040<<endl;

  cout<<"start looping run10 2040 files."<<endl;
  for(int ifile = 0; ifile < Nfiles_2040; ifile++){
    cout<<"Run10FileList2040["<<ifile<<"] = "<<Run10FileList2040[ifile]<<endl;
    file_string = Run10FileList2040[ifile];

    TFile* file = new TFile(file_string.c_str());
    if (! file->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }

    cbin = 1;
    int trigbin, partbin;    
    GetTrigPartBin(file_string,trigbin,partbin);
    if(type == 0) h1_JF[cbin][trigbin][partbin] = new TH1D(*(TH1D*)file->Get("JF_INC"));
    if(type == 1) h1_JF[cbin][trigbin][partbin] = new TH1D(*(TH1D*)file->Get("JF_PI0"));
    if(type == 2) h1_JF[cbin][trigbin][partbin] = new TH1D(*(TH1D*)file->Get("JF_DEC"));
    if(type == 3) h1_JF[cbin][trigbin][partbin] = new TH1D(*(TH1D*)file->Get("JF_DIR"));
    bin.str("");
    bin << cbin << "_p"<<trigbin<<"_h"<<partbin;
    name = "JFrun7_c"+bin.str();
    h1_JF[cbin][trigbin][partbin]->SetName(name.c_str());
  }

  if(type == 0) name = "C1_TRIGPT";
  if(type == 1) name = "piC1_TRIGPT";
  if(type == 2) name = "piC1_MWTRIGPT";
  h1_trig_pt[1] = new TH1D(*(TH1D*)fin2040->Get(name.c_str()));   
  
  for(int itrig=0; itrig<NPT_TRIG; itrig++){
    int binmin = h1_trig_pt[1]->FindBin(trig_pt_range[itrig]);
    int binmax = h1_trig_pt[1]->FindBin(trig_pt_range[itrig+1]);
    if(type == 0 || type == 1) {
      Ntrigs[1][itrig] = h1_trig_pt[1]->Integral(binmin, binmax-1);
      //cout<<"Ntrigs = "<<Ntrigs[ic][itrig]<<" = "<<GetNTriggers(h1_trig_pt[ic], trig_pt_range[itrig], trig_pt_range[itrig+1])<<endl;
    }
    if(type == 2) {
      if(itrig<3){
	Ntrigs[1][itrig] = h1_trig_pt[1]->GetBinContent(itrig+1);
      }
      else{
	Ntrigs[1][itrig] = h1_trig_pt[1]->GetBinContent(itrig+2);
      }
    }
  }
  for(int ipart=0; ipart<NPT_PART; ipart++){
    CombinePtHistos(0, Ntrigs[1][0],Ntrigs[1][1],h1_JF[1][0][ipart], h1_JF[1][1][ipart], h1_JF_COMB[1][0][ipart]);
    CombinePtHistos(2, Ntrigs[1][2],Ntrigs[1][3],h1_JF[1][2][ipart], h1_JF[1][3][ipart], h1_JF_COMB[1][1][ipart]);
    fout->cd();
    h1_JF_COMB[1][0][ipart]->Write();
    h1_JF_COMB[1][1][ipart]->Write();
  }
  cout<<"finish writing out pt combined jfs."<<endl;

  //Centrality combining.
  for(int itrig=0; itrig<2; itrig++){
    for(int ipart=0; ipart<NPT_PART; ipart++){
      bin.str("");
      bin << "p"<<itrig<<"_h"<<ipart;
      name = "JFrun7_COMBTOT_c0_" + bin.str();
      cout<<"name = "<<name.c_str()<<endl;
      h1_JF_COMBTOT[itrig][ipart] = new TH1D(name.c_str(),name.c_str(),30,0.0,PI);
      //cout<<"t1: "<<h1_JF_COMB[ic][itrig][ipart]->GetName()<<"; t2: "<<h1_JF_COMB[ic+1][itrig][ipart]->GetName()<<endl;
      Combine1DHistos(h1_JF_COMB[0][itrig][ipart],h1_JF_COMB[1][itrig][ipart],h1_JF_COMBTOT[itrig][ipart]);
      fout->cd();
      h1_JF_COMBTOT[itrig][ipart]->Write();
    }
  }
  
  cout<<"finish writing out."<<endl;
}
