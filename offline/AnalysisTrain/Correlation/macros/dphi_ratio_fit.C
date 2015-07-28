#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <vector>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>

using namespace std;
const int N_RUNS = 538;

int dphi_ratio_fit(const char* flist, const string fin, const string fout, char* chi2_off_run)
{
  vector<string> FileList;
  char buffer[256];
  cout << "opening input file: " << flist << endl;
  ifstream run_file_list(flist);
  if( !run_file_list.is_open())
    return 1;

  while( run_file_list.getline(buffer, 256, '\n') ) {
    char this_file[512];
    strcpy( this_file, buffer );
    string s(this_file);
    FileList.push_back(s);
  }
  cout << "Done making FileList vector." <<endl;
  const int N_RUNS = FileList.size();
  cout<<"# of run files: "<<N_RUNS<<endl;

  string file_string;
  string run_string[N_RUNS];
  double run_nbr[N_RUNS];

  for(int ifile = 0; ifile < N_RUNS; ifile++){
    file_string = FileList[ifile];
    cout<<"file_string: "<<file_string<<endl;
    int p1 = file_string.find_last_of(".");
    run_string[ifile] = file_string.substr(p1-6,6);
    cout<<"run_string: "<<run_string[ifile].c_str()<<endl;
   
    run_nbr[ifile] = atof(run_string[ifile].c_str());
    cout<<"run_nbr[ifile] = "<<run_nbr[ifile]<<endl;
  }

  ostringstream bin;
  string name;

  TH1D* dphi_ratio[N_RUNS];
  TCanvas* can_dphi_ratio[27];
  TF1* ratio_fit[N_RUNS];
  //double fit_val[N_RUNS];
  double chi2_per_ndf[N_RUNS];
  TGraph* g_chi2;
  int nBad = 0;

  for(int i=0; i<27; i++){
    bin.str("");
    bin<<i;
    name = "can_dphi_ratio_" + bin.str();
    can_dphi_ratio[i] = new TCanvas(name.c_str(),name.c_str());
    can_dphi_ratio[i]->Divide(5,4);
  }


  FILE* bad_chi2_file;
  if ((bad_chi2_file = fopen(chi2_off_run, "w")) == 0){
    fprintf(stderr, "File open error: %s\n", chi2_off_run);
    return 2;
  }

  TFile* fdphi = new TFile(fin.c_str());
  TFile* foutput = new TFile(fout.c_str(), "recreate");
  foutput->cd();
  for(int ifile=0; ifile<N_RUNS; ifile++){
    bin.str("");
    bin << run_string[ifile];
    name = "dphi_ratio_" + bin.str();
    //cout<<"name = "<<name.c_str()<<endl;
    dphi_ratio[ifile] = (TH1D*)fdphi->Get(name.c_str());
    // dphi_ratio[ifile]->SetMarkerSize(0.6);
    // dphi_ratio[ifile]->SetMarkerStyle(20);
   
   
    int can_nbr = (int)(ifile/20);
    //    cout<<"can_nbr = "<<can_nbr<<endl;
    int pad_nbr = (int)(ifile%20);
    //cout<<"pad_nbr = "<<pad_nbr<<endl;

    can_dphi_ratio[can_nbr]->cd(pad_nbr+1);
    dphi_ratio[ifile]->Draw();
    

    name = "fit_dphi_ratio_" + bin.str();
    ratio_fit[ifile] = new TF1(name.c_str(),"pol0(0)");
    dphi_ratio[ifile]->Fit(name.c_str(),"BQ");


    double fit_chi2 = ratio_fit[ifile]->GetChisquare();
    int fit_ndf = ratio_fit[ifile]->GetNDF();
    chi2_per_ndf[ifile] = fit_chi2/fit_ndf;
    //cout<<"chi2/ndf = "<<chi2_per_ndf[ifile]<<endl;
    if(chi2_per_ndf[ifile] > 1.7){
      nBad++;
      fprintf(bad_chi2_file, "%f %f\n", run_nbr[ifile], chi2_per_ndf[ifile]);
      dphi_ratio[ifile]->Write();
    }
  }
  cout<<"# of bad runs: "<<nBad<<endl;

  g_chi2 = new TGraph(N_RUNS,run_nbr,chi2_per_ndf);
  g_chi2->GetXaxis()->SetTitle("run #");
  g_chi2->GetYaxis()->SetTitle("#chi^{2}/NDF");
  g_chi2->SetName("chi2_vs_run");

  foutput->cd();
  for(int i=0; i<27; i++){
    can_dphi_ratio[i]->Write();
  }
  g_chi2->Write();

  return 0;
}
