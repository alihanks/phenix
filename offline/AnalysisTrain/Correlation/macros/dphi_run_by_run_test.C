#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <vector>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>

using namespace std;

const double PI = TMath::Pi();

double Round(double number)
{
  number*=100;
  number+=0.5;
  number=(int)number;
  number/=100;
  return number;
}

// int dphi_run_by_run(const char* flist, const string ftaxi, const string fout)
// {
//   vector<string> FileList;
//   char buffer[256];
//   cout << "opening input file: " << flist << endl;
//   ifstream run_file_list(flist);
//   if( !run_file_list.is_open())
//     return 1;

//   while( run_file_list.getline(buffer, 256, '\n') ) {
//     // Push all paths that exist onto list
//     char this_file[512];
//     strcpy( this_file, buffer );
//     string s(this_file);
//     FileList.push_back(s);
//   }
//   cout << "Done making FileList vector." <<endl;
//   const int N_RUNS = FileList.size();
//   cout<<"# of run files: "<<N_RUNS<<endl;

//   string file_string;
//   string run_string[N_RUNS];
//   int run_nbr[N_RUNS];

//   string name;
//   string name_mix;
//   ostringstream bin;

//   TH1D* dphi_1d[N_RUNS];
//   TH1D* dphi_1d_merged;

//   TH1D* dphi_ratio[N_RUNS];
  
//   TFile* foutput = new TFile(fout.c_str(), "recreate");

//   //statistics merged
//   TFile* fTaxi = new TFile (ftaxi.c_str());
//   cout<<"dealing with Merged taxi input."<<endl;
//   for(int ic = 0; ic < 10; ic++){
//     cout<<"merged: ic = "<<ic<<endl;
//     for(int iv = 0; iv < 10; iv++){
//       cout<<"merged: iv = "<<iv<<endl;
//       bin.str("");
//       bin << iv <<"_" << ic;
//       name = "h3_deltaphi_" + bin.str();
//       TH3D* temp3D = (TH3D*)fTaxi->Get(name.c_str());
//       temp3D->Sumw2();
      
//       name = "h1_dphi_" + bin.str();

//       TH1D* temp = (TH1D*)temp3D->Project3D("z");
//       if( iv==0 && ic == 0 ){
// 	dphi_1d_merged = new TH1D (*temp);
// 	dphi_1d_merged->SetNameTitle(name.c_str(),"");
//       }
//       else {
// 	dphi_1d_merged->Add(temp);
//       }
//       cout<<"merged: finish combining zvtx bins"<<endl;
//     }
//   }
//   cout<<"merged: finish combining centrality bins."<<endl;
//   foutput->cd();
//   dphi_1d_merged->Write();

//   //individual run

//   for(int ifile = 0; ifile < N_RUNS; ifile++){
//     file_string = FileList[ifile];
//     cout<<"file_string: "<<file_string<<endl;
//     int p1 = file_string.find_last_of(".");
//     run_string[ifile] = file_string.substr(p1-6,6);
//     cout<<"run_string: "<<run_string[ifile].c_str()<<endl;
   
//     run_nbr[ifile] = atoi(run_string[ifile].c_str());
//     cout<<"Run_nbr "<<ifile<<" is "<<run_nbr[ifile]<<endl;

//     TFile* fin = new TFile(file_string.c_str());
//     if (! fin->IsOpen() ) {
//       cout << "cannot open file '" << file_string << "'\n";
//       continue;
//     }
//     for(int ic = 0; ic < 10; ic++){
//       cout<<"ic = "<<ic<<endl;
//       for (int iv = 0; iv < 10; iv++){
// 	bin.str("");
// 	bin << iv <<"_" << ic;
// 	name = "h3_deltaphi_" + bin.str();
// 	TH3D* temp3D = (TH3D*)fin->Get(name.c_str());
// 	temp3D->Sumw2();

// 	name = "h1_delphi_" + run_string[ifile];
// 	TH1D* temp = (TH1D*)temp3D->Project3D("z");
// 	if( iv==0&&ic==0 ) {
// 	  dphi_1d[ifile] = new TH1D(*temp);
// 	  dphi_1d[ifile]->SetNameTitle(name.c_str(),"");
// 	}
// 	else{
// 	  dphi_1d[ifile]->Add(temp);
// 	}
//       }
//       cout<<"finish combining zvtx bins."<<endl;
//     }
//     cout<<"Done combining centrality bins!"<<endl;
//     foutput->cd();
//     dphi_1d[ifile]->Write();

//     bin.str("");
//     bin<<ifile;
//     string ratio_run = "dphi_ratio_" + run_string[ifile];
//     dphi_ratio[ifile] = (TH1D*)dphi_1d[ifile]->Clone();
//     dphi_ratio[ifile]->SetNameTitle(ratio_run.c_str(),ratio_run.c_str());
//     dphi_ratio[ifile]->Divide(dphi_1d_merged);
//     dphi_ratio[ifile]->SetMarkerSize(0.6);
//     dphi_ratio[ifile]->SetMarkerStyle(20);
//     dphi_ratio[ifile]->Write();

//     fin->Close();
//     delete fin;
//   }

//   cout<<"finish writing out dphi for each run."<<endl;
  
//   return 0;
// }

//10/21/14
//dphi_run_by_run_v2("runfile.lst","/direct/phenix+hhj/hge/taxi/Run11AuAu200MinBias/4759/data/correlation_taxi4759.root","dphi_run_by_run_v2.root")
int dphi_run_by_run_v2(const char* flist, const string ftaxi, const string fout)
{
  vector<string> FileList;
  char buffer[256];
  cout << "opening input file: " << flist << endl;
  ifstream run_file_list(flist);
  if( !run_file_list.is_open())
    return 1;

  while( run_file_list.getline(buffer, 256, '\n') ) {
    // Push all paths that exist onto list
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
  int run_nbr[N_RUNS];

  string name;
  string name_mix;
  ostringstream bin;

  TH1D* dphi_1d[N_RUNS];
  TH1D* dphi_1d_mix[N_RUNS];
  TH1D* dphi_1d_merged;
  TH1D* dphi_1d_mix_merged;

  TH1D* dphi_ratio[N_RUNS];
  TH1D* dphi_mix_ratio[N_RUNS];
  TF1* ratio_fit[N_RUNS];
  TF1* ratio_fit_mix[N_RUNS];
  double chi2[N_RUNS];
  double chi2_mix[N_RUNS];

  TCanvas* can[27];//538 runs, 20 pads per canvas.
  TCanvas* can_mix[27];
  for(int i=0; i<27; i++){
    bin.str("");
    bin<<i;
    name = "can_dphi_ratio_" + bin.str();
    can[i] = new TCanvas(name.c_str(),name.c_str());
    can[i]->Divide(5,4);
    name_mix = "can_dphi_mix_ratio_" + bin.str();
    can_mix[i] = new TCanvas(name_mix.c_str(),name_mix.c_str());
    can_mix[i]->Divide(5,4);
  }

  TFile* foutput = new TFile(fout.c_str(), "recreate");

  //statistics merged
  TFile* fTaxi = new TFile (ftaxi.c_str());
  cout<<"dealing with Merged taxi input."<<endl;
  for(int ic = 0; ic < 4; ic++){
    cout<<"merged: ic = "<<ic<<endl;
    bin.str("");
    bin << ic;
    name = "h3_dphi_fold_c" + bin.str();
    TH3D* temp3D = (TH3D*)fTaxi->Get(name.c_str());
    name_mix = "h3_dphi_mix_fold_c" + bin.str();
    TH3D* temp3D_mix = (TH3D*)fTaxi->Get(name_mix.c_str()); 
      
    name = "h1_dphi_c" + bin.str(); 
    name_mix = "h1_dphi_mix_c" + bin.str();   
    TH1D* temp = (TH1D*)temp3D->Project3D("z");
    TH1D* temp_mix = (TH1D*)temp3D_mix->Project3D("z");
    
    if( ic == 0 ){
      dphi_1d_merged = new TH1D (*temp);
      dphi_1d_merged->SetName(name.c_str());
      dphi_1d_mix_merged = new TH1D (*temp_mix);
      dphi_1d_mix_merged->SetName(name_mix.c_str());
    }
    else {
      dphi_1d_merged->Add(temp);
      dphi_1d_mix_merged->Add(temp_mix);
    }
  }
  cout<<"merged: finish combining centrality bins."<<endl;
  foutput->cd();
  dphi_1d_merged->Write();
  dphi_1d_mix_merged->Write();

  //individual run

  for(int ifile = 0; ifile < N_RUNS; ifile++){
    file_string = FileList[ifile];
    cout<<"file_string: "<<file_string<<endl;
    int p1 = file_string.find_last_of(".");
    run_string[ifile] = file_string.substr(p1-6,6);
    cout<<"run_string: "<<run_string[ifile].c_str()<<endl;
   
    run_nbr[ifile] = atoi(run_string[ifile].c_str());
    cout<<"Run_nbr "<<ifile<<" is "<<run_nbr[ifile]<<endl;

    TFile* fin = new TFile(file_string.c_str());
    if (! fin->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }
    for(int ic = 0; ic < 4; ic++){
      bin.str("");
      bin << ic;
      name = "h3_dphi_fold_c" + bin.str();
      TH3D* temp3D = (TH3D*)fin->Get(name.c_str());
      name_mix = "h3_dphi_mix_fold_c" + bin.str();
      TH3D* temp3D_mix = (TH3D*)fin->Get(name_mix.c_str());
                 
      name = "h1_delphi_" + run_string[ifile];
      TH1D* temp = (TH1D*)temp3D->Project3D("z");
      name_mix = "h1_delphi_mix_" + run_string[ifile];
      TH1D* temp_mix = (TH1D*)temp3D_mix->Project3D("z");
      
      if( ic==0 ) {
	dphi_1d[ifile] = new TH1D(*temp);
	dphi_1d[ifile]->SetName(name.c_str());
	dphi_1d_mix[ifile] = new TH1D(*temp_mix);
	dphi_1d_mix[ifile]->SetName(name_mix.c_str());
      }
      else{
	dphi_1d[ifile]->Add(temp);
	dphi_1d_mix[ifile]->Add(temp_mix);
      }
    }
    cout<<"Done combining centrality bins!"<<endl;
    foutput->cd();
    dphi_1d[ifile]->Write();
    dphi_1d_mix[ifile]->Write();
    
    string ratio_run = "dphi_ratio_" + run_string[ifile];
    cout<<"ratio_run string: "<<ratio_run.c_str()<<endl;
    dphi_ratio[ifile] = (TH1D*)dphi_1d[ifile]->Clone();
    dphi_ratio[ifile]->SetNameTitle(ratio_run.c_str(),ratio_run.c_str());
    dphi_ratio[ifile]->Divide(dphi_1d_merged);
    dphi_ratio[ifile]->SetMarkerSize(0.6);
    dphi_ratio[ifile]->SetMarkerStyle(20);

    name = "fit_dphi_ratio_" + run_string[ifile];
    ratio_fit[ifile] = new TF1(name.c_str(),"pol0(0)");
    dphi_ratio[ifile]->Fit(name.c_str(),"BQ");
    chi2[ifile] = Round(ratio_fit[ifile]->GetChisquare()/ratio_fit[ifile]->GetNDF());
    cout<<"chi2["<<ifile<<"] = "<<chi2[ifile]<<endl;
    dphi_ratio[ifile]->Write();

    ratio_run = "dphi_mix_ratio_" + run_string[ifile];
    cout<<"ratio_run(mix) string: "<<ratio_run.c_str()<<endl;
    dphi_mix_ratio[ifile] = (TH1D*)dphi_1d_mix[ifile]->Clone();
    dphi_mix_ratio[ifile]->SetNameTitle(ratio_run.c_str(),ratio_run.c_str());
    dphi_mix_ratio[ifile]->Divide(dphi_1d_mix_merged);
    dphi_mix_ratio[ifile]->SetMarkerSize(0.6);
    dphi_mix_ratio[ifile]->SetMarkerStyle(20);

    name = "fit_dphi_mix_ratio_" + run_string[ifile];
    ratio_fit_mix[ifile] = new TF1(name.c_str(),"pol0(0)");
    dphi_mix_ratio[ifile]->Fit(name.c_str(),"BQ");
    chi2_mix[ifile] = Round(ratio_fit_mix[ifile]->GetChisquare()/ratio_fit_mix[ifile]->GetNDF());
    cout<<"chi2_mix["<<ifile<<"] = "<<chi2_mix[ifile]<<endl;
    dphi_mix_ratio[ifile]->Write();

    int can_nbr = (int)(ifile/20);
    int pad_nbr = (int)(ifile%20);
    cout<<"can_nbr: "<<can_nbr<<"; pad_nbr: "<<pad_nbr<<endl;

    can[can_nbr]->cd(pad_nbr+1);
    cout<<"can_nbr: "<<can_nbr<<"; pad_nbr: "<<pad_nbr<<endl;
    cout<<"dphi_ratio["<<ifile<<"] name: "<<dphi_ratio[ifile]->GetName()<<endl;
    dphi_ratio[ifile]->Draw();
    bin.str("");
    bin << "#chi^{2} = "<<chi2[ifile];
    TLatex *la_chi2 = new TLatex(0.7, 0.2, bin.str().c_str());
    la_chi2->SetNDC();
    la_chi2->Draw("same");

    can_mix[can_nbr]->cd(pad_nbr+1);
    cout<<"can_nbr: "<<can_nbr<<"; pad_nbr: "<<pad_nbr<<endl;
    cout<<"dphi_mix_ratio["<<ifile<<"] name: "<<dphi_mix_ratio[ifile]->GetName()<<endl;
    dphi_mix_ratio[ifile]->Draw();
    bin.str("");
    bin << "#chi^{2} = "<<chi2_mix[ifile];
    cout<<"chi2_mix: "<<bin.str().c_str()<<endl;
    TLatex *la_chi2_mix = new TLatex(0.7, 0.2, bin.str().c_str());
    la_chi2_mix->SetNDC();
    la_chi2_mix->Draw("same");

    cout<<"finish drawing"<<endl;
    fin->Close();
    delete fin;
  }

  foutput->cd();
  for(int i=0; i<27; i++){
    can[i]->Write();
    can_mix[i]->Write();
  }

  cout<<"finish writing out dphi for each run."<<endl;
  
  return 0;
}

int dphi_run_by_run_single(const char* flist, const string ftaxi, const string fout, const string name = "h3_dphi_fold_c0" )
{
  vector<string> FileList;
  char buffer[256];
  cout << "opening input file: " << flist << endl;
  ifstream run_file_list(flist);
  if( !run_file_list.is_open())
    return 1;

  while( run_file_list.getline(buffer, 256, '\n') ) {
    // Push all paths that exist onto list
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
  int run_nbr[N_RUNS];

  ostringstream bin;

  TH1D* dphi_1d[N_RUNS];
  TH1D* dphi_1d_merged;

  TH1D* dphi_ratio[N_RUNS];
  TF1* ratio_fit[N_RUNS];
  double chi2[N_RUNS];

  string cname;
  TCanvas* can[27];//538 runs, 20 pads per canvas.
  for(int i=0; i<27; i++){
    bin.str("");
    bin<<i;
    cname = "can_dphi_ratio_" + bin.str();
    can[i] = new TCanvas(cname.c_str(),cname.c_str());
    can[i]->Divide(5,4);
  }

  TFile* foutput = new TFile(fout.c_str(), "recreate");

  //statistics merged
  TFile* fTaxi = new TFile (ftaxi.c_str());
  cout<<"dealing with Merged taxi input."<<endl;
  
  TH3D* temp3D = (TH3D*)fTaxi->Get(name.c_str());

  string pname;
  pname = "h1_dphi_c0";
  TH1D* temp = (TH1D*)temp3D->Project3D("z");
  
  dphi_1d_merged = new TH1D (*temp);
  dphi_1d_merged->SetName(pname.c_str());

  foutput->cd();
  dphi_1d_merged->Write();

  //individual run

  for(int ifile = 0; ifile < N_RUNS; ifile++){
    file_string = FileList[ifile];
    cout<<"file_string: "<<file_string<<endl;
    int p1 = file_string.find_last_of(".");
    run_string[ifile] = file_string.substr(p1-6,6);
    //run_string[ifile] = file_string.substr(p1-4,4);
    cout<<"run_string: "<<run_string[ifile].c_str()<<endl;
   
    run_nbr[ifile] = atoi(run_string[ifile].c_str());
    cout<<"Run_nbr "<<ifile<<" is "<<run_nbr[ifile]<<endl;

    TFile* fin = new TFile(file_string.c_str());
    if (! fin->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }

    TH3D* temp3D = (TH3D*)fin->Get(name.c_str());
    
    pname = "h1_delphi_" + run_string[ifile];
    TH1D* temp = (TH1D*)temp3D->Project3D("z");
    
    dphi_1d[ifile] = new TH1D(*temp);
    dphi_1d[ifile]->SetName(pname.c_str());
   
    cout<<"Done combining centrality bins!"<<endl;
    foutput->cd();
    dphi_1d[ifile]->Write();
    
    string ratio_run = "dphi_ratio_" + run_string[ifile];
    cout<<"ratio_run string: "<<ratio_run.c_str()<<endl;
    dphi_ratio[ifile] = (TH1D*)dphi_1d[ifile]->Clone();
    dphi_ratio[ifile]->SetNameTitle(ratio_run.c_str(),ratio_run.c_str());
    dphi_ratio[ifile]->Divide(dphi_1d_merged);
    dphi_ratio[ifile]->GetYaxis()->SetLabelSize(0.06);
    dphi_ratio[ifile]->GetXaxis()->SetLabelSize(0.06);
    dphi_ratio[ifile]->GetXaxis()->SetTitleSize(0.05);
    dphi_ratio[ifile]->SetMarkerSize(0.6);
    dphi_ratio[ifile]->SetMarkerStyle(20);
    double R = dphi_ratio[ifile]->Integral("width")/PI;
    dphi_ratio[ifile]->Scale(1/R);

    int can_nbr = (int)(ifile/20);
    int pad_nbr = (int)(ifile%20);
    cout<<"can_nbr: "<<can_nbr<<"; pad_nbr: "<<pad_nbr<<endl;

    can[can_nbr]->cd(pad_nbr+1);

    pname = "fit_dphi_ratio_" + run_string[ifile];
    ratio_fit[ifile] = new TF1(pname.c_str(),"pol0(0)");
    dphi_ratio[ifile]->Fit(pname.c_str(),"BQ");
    chi2[ifile] = Round(ratio_fit[ifile]->GetChisquare()/ratio_fit[ifile]->GetNDF());
    cout<<"chi2["<<ifile<<"] = "<<chi2[ifile]<<endl;
    dphi_ratio[ifile]->Write();

    cout<<"can_nbr: "<<can_nbr<<"; pad_nbr: "<<pad_nbr<<endl;
    cout<<"dphi_ratio["<<ifile<<"] name: "<<dphi_ratio[ifile]->GetName()<<endl;
    dphi_ratio[ifile]->Draw();
    bin.str("");
    bin << "#chi^{2} = "<<chi2[ifile];
    TLatex *la_chi2 = new TLatex(0.7, 0.2, bin.str().c_str());
    la_chi2->SetTextSize(0.06);
    la_chi2->SetNDC();
    la_chi2->Draw("same");

    cout<<"finish drawing"<<endl;
    fin->Close();
    delete fin;
  }

  foutput->cd();
  for(int i=0; i<27; i++){
    can[i]->Write();
  }

  cout<<"finish writing out dphi for each run."<<endl;
  
  return 0;
}

//check the effect of increasing pool depth etc by comparing two taxi output
int dphi_double_ratio(const char* flist, const string file1, const string file2, const string fout)
{
  vector<string> FileList;
  char buffer[256];
  cout << "opening input file: " << flist << endl;
  ifstream run_file_list(flist);
  if( !run_file_list.is_open())
    return 1;

  while( run_file_list.getline(buffer, 256, '\n') ) {
    // Push all paths that exist onto list
    char this_file[512];
    strcpy( this_file, buffer );
    string s(this_file);
    FileList.push_back(s);
  }
  cout << "Done making FileList vector." <<endl;
  const int N_RUNS = FileList.size();
  cout<<"# of run files: "<<N_RUNS<<endl;

  ostringstream bin;
  string name;
  string file_string;
  string run_string[N_RUNS];
  int run_nbr[N_RUNS];

  TH1D* h1_fin1[N_RUNS];
  TH1D* h1_fin2[N_RUNS];
  TH1D* h1_double_ratio[N_RUNS];

  TFile* fin1 = new TFile(file1.c_str());
  TFile* fin2 = new TFile(file2.c_str()); 
  TFile* foutput = new TFile(fout.c_str(), "recreate");

  TCanvas* can[27];//538 runs, 20 pads per canvas.
  TCanvas* can_both[27];
  for(int i=0; i<27; i++){
    bin.str("");
    bin<<i;
    name = "can_dphi_double_ratio_" + bin.str();
    can[i] = new TCanvas(name.c_str(),name.c_str());
    can[i]->Divide(5,4);
    name = "can_dphi_ratio_" + bin.str();
    can_both[i] = new TCanvas(name.c_str(),name.c_str());
    can_both[i]->Divide(5,4);
  }

  for(int ifile = 0; ifile < N_RUNS; ifile++){
  //for(int ifile = 0; ifile < 20; ifile++){
    file_string = FileList[ifile];
    cout<<"file_string: "<<file_string<<endl;
    int p1 = file_string.find_last_of(".");
    run_string[ifile] = file_string.substr(p1-6,6);
    cout<<"run_string: "<<run_string[ifile].c_str()<<endl;
   
    run_nbr[ifile] = atoi(run_string[ifile].c_str());
    
    name = "dphi_ratio_" + run_string[ifile];
    cout<<"name: "<<name.c_str()<<endl;
    TH1D* temp1D = (TH1D*)fin1->Get(name.c_str());
    h1_fin1[ifile] = new TH1D(*(TH1D*)temp1D);
    //name += "_taxi5286";//|zvtx|<20cm
    //name += "_taxi4755";//pd=100
    name += "_taxi5132";//pd=500
    cout<<"name: "<<name.c_str()<<endl;
    h1_fin1[ifile]->SetName(name.c_str());

    name = "dphi_ratio_" + run_string[ifile];
    cout<<"name: "<<name.c_str()<<endl;
    TH1D* temp1D_2 = (TH1D*)fin2->Get(name.c_str());
    h1_fin2[ifile] = new TH1D(*(TH1D*)temp1D_2);
    //name += "_taxi5132";//pd=500
    name += "_taxi5152";//pd=1000, 9 runs missing though
    //name += "_taxi5190";//pd=25
    //name += "_taxi5288";//|zvtx|<30cm
    // name += "_taxi5292";//|zvtx|<20cm, more restricted zed
    cout<<"name: "<<name.c_str()<<endl;
    h1_fin2[ifile]->SetName(name.c_str());

    h1_double_ratio[ifile] = new TH1D(*(TH1D*)h1_fin1[ifile]);
    h1_double_ratio[ifile]->Divide(h1_fin2[ifile]);
    
    int can_nbr = (int)(ifile/20);
    int pad_nbr = (int)(ifile%20); 

    can_both[can_nbr]->cd(pad_nbr+1);
    h1_fin1[ifile]->SetLineColor(4);
    h1_fin1[ifile]->SetMarkerColor(4);
    h1_fin1[ifile]->Draw();
    h1_fin2[ifile]->SetLineColor(2);
    h1_fin2[ifile]->SetMarkerColor(2);
    h1_fin2[ifile]->Draw("same");
    can[can_nbr]->cd(pad_nbr+1);
    h1_double_ratio[ifile]->Draw();
  }
  foutput->cd();
  for(int i=0; i<27; i++){
    can_both[i]->Write();
    can[i]->Write();
  }
}

//look at run by run phi distributions, zed distributions, pool vertex counter for 0-10% events
int phi_run_by_run(const char* flist, const string fout)//for some reason, canvases are empty
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
  int run_nbr[N_RUNS];

  string name;
  ostringstream bin;

  TH1D* phi[N_RUNS];
  TH1D* phi_mix[N_RUNS];
  TH2D* zed_vs_phid[N_RUNS];
  TH2D* zed_vs_phid_mix[N_RUNS];

  TH2D* pcount[N_RUNS];
  TH2D* pcount_real[N_RUNS];

  TH2D* trig_pt_zvtx[N_RUNS];
  TH2D* part_pt_zvtx[N_RUNS];

  // TCanvas* can_phi[27];
  // TCanvas* can_pool_counter[27];
  // for(int i=0; i<27; i++){
  //    bin.str("");
  //   bin<<i;
  //   name = "can_phi_" + bin.str();
  //   can_phi[i] = new TCanvas(name.c_str(),name.c_str());
  //   can_phi[i]->Divide(5,4);
  //   name = "can_pool_counter_" + bin.str();
  //   can_pool_counter[i] = new TCanvas(name.c_str(),name.c_str());
  //   can_pool_counter[i]->Divide(5,4);
  // }

  TFile* foutput = new TFile(fout.c_str(), "recreate");

  for(int ifile=0; ifile<N_RUNS; ifile++){
    file_string = FileList[ifile];
    cout<<"file_string: "<<file_string<<endl;
    int p1 = file_string.find_last_of(".");
    run_string[ifile] = file_string.substr(p1-6,6);
    cout<<"run_string: "<<run_string[ifile].c_str()<<endl;
   
    run_nbr[ifile] = atoi(run_string[ifile].c_str());
    cout<<"Run_nbr "<<ifile<<" is "<<run_nbr[ifile]<<endl;

    TFile* fin = new TFile(file_string.c_str());
    if (! fin->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }

    TH1D* temp = (TH1D*)fin->Get("h1_phi");
    phi[ifile] = new TH1D(*(TH1D*)temp);
    name = "h1_phi_" + run_string[ifile];
    cout<<"name: "<<name.c_str()<<endl;
    phi[ifile]->SetNameTitle(name.c_str(),name.c_str());
    phi[ifile]->SetLineColor(4);
    phi[ifile]->SetMarkerColor(4);
    phi[ifile]->SetMarkerStyle(20);
    phi[ifile]->SetMarkerSize(0.7);

    TH1D* temp_mix = (TH1D*)fin->Get("h1_phi_mix");
    phi_mix[ifile] = new TH1D(*(TH1D*)temp_mix);
    name = "h1_phi_mix_" + run_string[ifile];
    cout<<"name_mix: "<<name.c_str()<<endl;
    phi_mix[ifile]->SetNameTitle(name.c_str(),name.c_str());
    phi_mix[ifile]->Scale(0.01);
    phi_mix[ifile]->SetLineColor(2);
    phi_mix[ifile]->SetMarkerColor(2);
    phi_mix[ifile]->SetMarkerStyle(20);
    phi_mix[ifile]->SetMarkerSize(0.7);
    
    zed_vs_phid[ifile] = new TH2D(*(TH2D*)fin->Get("h2_phi_zed_aft"));
    name = "h2_phi_zed_" + run_string[ifile];
    zed_vs_phid[ifile]->SetNameTitle(name.c_str(),name.c_str());

    zed_vs_phid_mix[ifile] = new TH2D(*(TH2D*)fin->Get("h2_phi_zed_mix"));
    name = "h2_phi_zed_mix_" + run_string[ifile];
    zed_vs_phid_mix[ifile]->SetNameTitle(name.c_str(),name.c_str());

    pcount[ifile] = new TH2D(*(TH2D*)fin->Get("h2_pool_counter"));
    name = "h2_pool_counter_" + run_string[ifile];
    pcount[ifile]->SetNameTitle(name.c_str(),name.c_str());

    pcount_real[ifile] = new TH2D(*(TH2D*)fin->Get("h2_event_counter"));
    name = "h2_event_counter_" + run_string[ifile];
    pcount_real[ifile]->SetNameTitle(name.c_str(),name.c_str());

    trig_pt_zvtx[ifile] = new TH2D(*(TH2D*)fin->Get("h2_trig_pt_zvtx_inc"));
    name = "h2_trig_pt_zvtx_" + run_string[ifile];
    trig_pt_zvtx[ifile]->SetNameTitle(name.c_str(),name.c_str());

    part_pt_zvtx[ifile] = new TH2D(*(TH2D*)fin->Get("h2_part_pt_zvtx"));
    name = "h2_part_pt_zvtx_" + run_string[ifile];
    part_pt_zvtx[ifile]->SetNameTitle(name.c_str(),name.c_str());

    foutput->cd();
    phi[ifile]->Write();
    phi_mix[ifile]->Write();
    zed_vs_phid[ifile]->Write();
    zed_vs_phid_mix[ifile]->Write();
    pcount[ifile]->Write();
    pcount_real[ifile]->Write();
    trig_pt_zvtx[ifile]->Write();
    part_pt_zvtx[ifile]->Write();

    // int can_nbr = (int)(ifile/20);
    // int pad_nbr = (int)(ifile%20);
    // cout<<"can_nbr = "<<can_nbr<<"; pad_nbr = "<<pad_nbr<<endl;
    // cout<<"name for phi["<<ifile<<"] is: "<<phi[ifile]->GetName()<<endl;
    // cout<<"name for phi_mix["<<ifile<<"] is: "<<phi_mix[ifile]->GetName()<<endl;

    // can_phi[can_nbr]->cd(pad_nbr+1);
    // phi[ifile]->Draw();    
    // phi_mix[ifile]->Draw("same");

    // can_pool_counter[can_nbr]->cd(pad_nbr+1);
    // pcount[ifile]->Draw("colz");

    fin->Close();
    delete fin;
  }

  // foutput->cd();
  // for(int i=0; i<27; i++){
  //   can_phi[i]->Write();
  //   can_pool_counter[i]->Write();
  // }
  cout<<"finish writing out phi and pcount for each run."<<endl;
  return 0;
}
//plot # of triggers, pairs as a function of run #
/*
int triggers_run_by_run(const char* flist, const string fout)
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

  double trigpt[5] = {5.0,7.0,9.0,12.0,15.0};
  string file_string;
  string run_string[N_RUNS];
  double run_nbr[N_RUNS];
  ostringstream bin;
  string name;

  TH1D* h1_zvertex[N_RUNS];
  TH1D* h1_num_trig_inc[4][N_RUNS];//[cent]
  TH1D* h1_num_trig_pi0[4][N_RUNS];
  TH1D* h1_num_trig_dec[4][N_RUNS];
  TH3D* h3_num_pair_inc[4][N_RUNS];
  TH3D* h3_num_pair_pi0[4][N_RUNS];
  TH2D* h2_num_pair_dec[4][4][N_RUNS];//[cent][itrig]

  TGraph* gNtrigInc[4][4];//[cent][itrig]
  TGraph* gNtrigPi0[4][4];
  TGraph* gNtrigDec[4][4];
  TGraph* gNpairInc[4][4];//[cent][itrig]
  TGraph* gNpairPi0[4][4];
  TGraph* gNpairDec[4][4];

  double ntrig_inc[4][4][N_RUNS];//[cent][itrig]
  double ntrig_pi0[4][4][N_RUNS];
  double ntrig_dec[4][4][N_RUNS];
  double npair_inc[4][4][N_RUNS];
  double npair_pi0[4][4][N_RUNS];
  double npair_dec[4][4][N_RUNS];

  TFile* foutput = new TFile (fout.c_str(),"recreate");
  for(int ifile = 0; ifile < N_RUNS; ifile++){
    file_string = FileList[ifile];
    //cout<<"file_string: "<<file_string<<endl;
    int p1 = file_string.find_last_of(".");
    run_string[ifile] = file_string.substr(p1-6,6);
   
    run_nbr[ifile] = atoi(run_string[ifile].c_str());
    //cout<<"Run_nbr "<<ifile<<" is "<<run_nbr[ifile]<<endl;

    TFile* fin = new TFile(file_string.c_str());
    if (! fin->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }

    h1_zvertex[ifile] = new TH1D(*(TH1D*)fin->Get("h1_zvertex"));
    double nevt = h1_zvertex[ifile]->Integral();
    //cout<<"h1_zvertex["<<ifile<<"] = "<<nevt<<endl;

    for(int ic=0; ic<4; ic++){
      bin.str("");
      bin << ic;
      name = "h1_trig_pt_inc_c" + bin.str();
      h1_num_trig_inc[ic][ifile] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
      name = "h1_trig_pt_pi0_c" + bin.str();
      h1_num_trig_pi0[ic][ifile] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
      name = "h1_trig_pt_dec_c" + bin.str();
      h1_num_trig_dec[ic][ifile] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
      name = "h3_dphi_fold_c" + bin.str();
      h3_num_pair_inc[ic][ifile] = new TH3D(*(TH3D*)fin->Get(name.c_str()));
      name = "h3_dphi_pi0_fold_c" + bin.str();
      h3_num_pair_pi0[ic][ifile] = new TH3D(*(TH3D*)fin->Get(name.c_str()));

      for(int ip=0; ip<4; ip++){
      int x1 = h1_num_trig_inc[ic][ifile]->GetXaxis()->FindBin(trigpt[ip]);
      int x2 = h1_num_trig_inc[ic][ifile]->GetXaxis()->FindBin(trigpt[ip+1]);
      ntrig_inc[ic][ip][ifile] = h1_num_trig_inc[ic][ifile]->Integral(x1,x2-1);
      ntrig_inc[ic][ip][ifile]/=nevt;
      //cout<<"ntrig_inc["<<ic<<"]["<<ip<<"]["<<ifile<<"] = "<<ntrig_inc[ic][ip][ifile]<<endl;

      x1 = h1_num_trig_pi0[ic][ifile]->GetXaxis()->FindBin(trigpt[ip]);
      x2 = h1_num_trig_pi0[ic][ifile]->GetXaxis()->FindBin(trigpt[ip+1]);
      ntrig_pi0[ic][ip][ifile] = h1_num_trig_pi0[ic][ifile]->Integral(x1,x2-1);
      ntrig_pi0[ic][ip][ifile]/=nevt;
      //cout<<"ntrig_pi0["<<ic<<"]["<<ip<<"]["<<ifile<<"] = "<<ntrig_pi0[ic][ip][ifile]<<endl;

      if(ip<3) ntrig_dec[ic][ip][ifile] = h1_num_trig_dec[ic][ifile]->Integral(ip+1,ip+1);
      else ntrig_dec[ic][ip][ifile] = h1_num_trig_dec[ic][ifile]->Integral(ip+2,ip+2);
      ntrig_dec[ic][ip][ifile]/=nevt;
      //cout<<"ntrig_dec["<<ic<<"]["<<ip<<"]["<<ifile<<"] = "<<ntrig_dec[ic][ip][ifile]<<endl;
      
      x1 = h3_num_pair_inc[ic][ifile]->GetXaxis()->FindBin(trigpt[ip]);
      x2 = h3_num_pair_inc[ic][ifile]->GetXaxis()->FindBin(trigpt[ip+1]);
      npair_inc[ic][ip][ifile] = h3_num_pair_inc[ic][ifile]->Integral(x1,x2-1);
      npair_inc[ic][ip][ifile]/=nevt;
      //cout<<"npair_inc["<<ic<<"]["<<ip<<"]["<<ifile<<"] = "<<npair_inc[ic][ip][ifile]<<endl;

      x1 = h3_num_pair_pi0[ic][ifile]->GetXaxis()->FindBin(trigpt[ip]);
      x2 = h3_num_pair_pi0[ic][ifile]->GetXaxis()->FindBin(trigpt[ip+1]);
      npair_pi0[ic][ip][ifile] = h3_num_pair_pi0[ic][ifile]->Integral(x1,x2-1);
      npair_pi0[ic][ip][ifile]/=nevt;
      //cout<<"npair_pi0["<<ic<<"]["<<ip<<"]["<<ifile<<"] = "<<npair_pi0[ic][ip][ifile]<<endl;

      bin.str("");
      if(ip<3) bin << ip <<"_c"<<ic;
      else bin << ip+1 <<"_c"<<ic;
      name = "h2_dphi_dec_fold_p" + bin.str();
      //cout<<"name: "<< name.c_str()<<endl;
      h2_num_pair_dec[ic][ip][ifile] = new TH2D(*(TH2D*)fin->Get(name.c_str()));
      npair_dec[ic][ip][ifile] = h2_num_pair_dec[ic][ip][ifile]->Integral();
      npair_dec[ic][ip][ifile]/=nevt;
      //cout<<"npair_dec["<<ic<<"]["<<ip<<"]["<<ifile<<"] = "<<npair_dec[ic][ip][ifile]<<endl;
      }
    }
    fin->Close();
    delete fin;
  }
  for(int ic=0; ic<4; ic++){
    for(int ip=0; ip<4; ip++){
      bin.str("");
      bin << ic << "_p"<<ip;
      gNtrigInc[ic][ip] = new TGraph(N_RUNS,run_nbr,ntrig_inc[ic][ip]);
      name = "gNtrigInc_c" + bin.str();
      gNtrigInc[ic][ip]->SetName(name.c_str());
      gNtrigPi0[ic][ip] = new TGraph(N_RUNS,run_nbr,ntrig_pi0[ic][ip]);
      name = "gNtrigPi0_c" + bin.str();
      gNtrigPi0[ic][ip]->SetName(name.c_str());
      gNtrigDec[ic][ip] = new TGraph(N_RUNS,run_nbr,ntrig_dec[ic][ip]);
      name = "gNtrigDec_c" + bin.str();
      gNtrigDec[ic][ip]->SetName(name.c_str());

      gNpairInc[ic][ip] = new TGraph(N_RUNS,run_nbr,npair_inc[ic][ip]);
      name = "gNpairInc_c" + bin.str();
      gNpairInc[ic][ip]->SetName(name.c_str());
      gNpairPi0[ic][ip] = new TGraph(N_RUNS,run_nbr,npair_pi0[ic][ip]);
      name = "gNpairPi0_c" + bin.str();
      gNpairPi0[ic][ip]->SetName(name.c_str());
      gNpairDec[ic][ip] = new TGraph(N_RUNS,run_nbr,npair_dec[ic][ip]);
      name = "gNtrigDec_c" + bin.str();
      gNpairDec[ic][ip]->SetName(name.c_str());

      foutput->cd();
      gNtrigInc[ic][ip]->Write();
      gNtrigPi0[ic][ip]->Write();
      gNtrigDec[ic][ip]->Write();
      gNpairInc[ic][ip]->Write();
      gNpairPi0[ic][ip]->Write();
      gNpairDec[ic][ip]->Write();

    }
  }
  return 0;
}
*/
//plot # of triggers, pairs as a function of run #
int triggers_run_by_run_v2(const char* flist, const string fout)
{
  const int NCENT = 2;
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
  ostringstream bin;
  string name;

  TH1D* h1_zvertex[N_RUNS];
  TH1D* h1_num_trig_inc[2][N_RUNS];//[cent]
  TH1D* h1_num_trig_pi0[2][N_RUNS];
  TH1D* h1_num_trig_dec[2][N_RUNS];
  TH3D* h3_num_pair_inc[2][N_RUNS];
  TH3D* h3_num_pair_pi0[2][N_RUNS];
  TH2D* h2_num_pair_dec[2][4][N_RUNS];//[cent][itrig]

  TH1D* hNtrigInc[2];
  TH1D* hNtrigPi0[2];
  TH1D* hNtrigDec[2];
  TH1D* hNpairInc[2];
  TH1D* hNpairPi0[2];
  TH1D* hNpairDec[2];

  TGraph* gNtrigInc[2];//[cent]
  TGraph* gNtrigPi0[2];
  TGraph* gNtrigDec[2];
  TGraph* gNpairInc[2];//[cent]
  TGraph* gNpairPi0[2];
  TGraph* gNpairDec[2];
  
  TF1* fNtrigIncMean[2];
  TF1* fNtrigIncRMS1[2];
  TF1* fNtrigIncRMS2[2];
  TF1* fNtrigPi0Mean[2];
  TF1* fNtrigPi0RMS1[2];
  TF1* fNtrigPi0RMS2[2];
  TF1* fNtrigDecMean[2];
  TF1* fNtrigDecRMS1[2];
  TF1* fNtrigDecRMS2[2];
  TF1* fNpairIncMean[2];
  TF1* fNpairIncRMS1[2];
  TF1* fNpairIncRMS2[2];
  TF1* fNpairPi0Mean[2];
  TF1* fNpairPi0RMS1[2];
  TF1* fNpairPi0RMS2[2];
  TF1* fNpairDecMean[2];
  TF1* fNpairDecRMS1[2];
  TF1* fNpairDecRMS2[2];
  
  double ntrig_inc[2][N_RUNS];//[cent]
  double ntrig_pi0[2][N_RUNS];
  double ntrig_dec[2][N_RUNS];
  double npair_inc[2][N_RUNS];
  double npair_pi0[2][N_RUNS];
  double npair_dec[2][N_RUNS];

  double mean_ntrig_inc[2];
  double rms_ntrig_inc[2];
  double mean_ntrig_pi0[2];
  double rms_ntrig_pi0[2];
  double mean_ntrig_dec[2];
  double rms_ntrig_dec[2];
  double mean_npair_inc[2];
  double rms_npair_inc[2];
  double mean_npair_pi0[2];
  double rms_npair_pi0[2];
  double mean_npair_dec[2];
  double rms_npair_dec[2];

  double avg_ntrig_inc[2];//[cent]
  double avg_ntrig_pi0[2];
  double avg_ntrig_dec[2];
  double avg_npair_inc[2];
  double avg_npair_pi0[2];
  double avg_npair_dec[2];

  for(int i=0; i<2; i++){
    avg_ntrig_inc[i] = 0.0;
    avg_ntrig_pi0[i] = 0.0;
    avg_ntrig_dec[i] = 0.0;
    avg_npair_inc[i] = 0.0;
    avg_npair_pi0[i] = 0.0;
    avg_npair_dec[i] = 0.0;
  }

  TFile* foutput = new TFile (fout.c_str(),"recreate");
  for(int ifile = 0; ifile < N_RUNS; ifile++){
    file_string = FileList[ifile];
    cout<<"file_string: "<<file_string<<endl;
    int p1 = file_string.find_last_of(".");
    run_string[ifile] = file_string.substr(p1-6,6);
   
    run_nbr[ifile] = atoi(run_string[ifile].c_str());
    cout<<"Run_nbr "<<ifile<<" is "<<run_nbr[ifile]<<endl;

    TFile* fin = new TFile(file_string.c_str());
    if (! fin->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }

    h1_zvertex[ifile] = new TH1D(*(TH1D*)fin->Get("h1_zvertex"));
    double nevt = h1_zvertex[ifile]->Integral();
    //cout<<"h1_zvertex["<<ifile<<"] = "<<nevt<<endl;

    for(int ic=0; ic<2; ic++){
      bin.str("");
      bin << ic;
      name = "h1_trig_pt_inc_c" + bin.str();
      h1_num_trig_inc[ic][ifile] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
      ntrig_inc[ic][ifile] = h1_num_trig_inc[ic][ifile]->Integral();
      ntrig_inc[ic][ifile]/=nevt;
      //cout<<"ntrig_inc["<<ic<<"]["<<ifile<<"] = "<<ntrig_inc[ic][ifile]<<endl;

      name = "h1_trig_pt_pi0_c" + bin.str();
      h1_num_trig_pi0[ic][ifile] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
      ntrig_pi0[ic][ifile] = h1_num_trig_pi0[ic][ifile]->Integral(51,150);
      ntrig_pi0[ic][ifile]/=nevt;
      //cout<<"ntrig_pi0["<<ic<<"]["<<ifile<<"] = "<<ntrig_pi0[ic][ifile]<<endl;

      name = "h1_trig_pt_dec_c" + bin.str();
      h1_num_trig_dec[ic][ifile] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
      ntrig_dec[ic][ifile] = h1_num_trig_dec[ic][ifile]->Integral(1,3)+h1_num_trig_dec[ic][ifile]->Integral(4,4);
      ntrig_dec[ic][ifile]/=nevt;
      //cout<<"ntrig_dec["<<ic<<"]["<<ifile<<"] = "<<ntrig_dec[ic][ifile]<<endl;

      name = "h3_dphi_fold_c" + bin.str();
      h3_num_pair_inc[ic][ifile] = new TH3D(*(TH3D*)fin->Get(name.c_str()));
      npair_inc[ic][ifile] = h3_num_pair_inc[ic][ifile]->Integral();
      npair_inc[ic][ifile]/=nevt;
      //cout<<"npair_inc["<<ic<<"]["<<ifile<<"] = "<<npair_inc[ic][ifile]<<endl;

      name = "h3_dphi_pi0_fold_c" + bin.str();
      h3_num_pair_pi0[ic][ifile] = new TH3D(*(TH3D*)fin->Get(name.c_str()));
      npair_pi0[ic][ifile] = h3_num_pair_pi0[ic][ifile]->Integral(51,150);
      npair_pi0[ic][ifile]/=nevt;
      //cout<<"npair_pi0["<<ic<<"]["<<ifile<<"] = "<<npair_pi0[ic][ifile]<<endl;

      for(int ip=0; ip<4; ip++){
	bin.str("");
	if(ip<3) bin << ip <<"_c"<<ic;
	else bin << ip+1 <<"_c"<<ic;
	name = "h2_dphi_dec_fold_p" + bin.str();
	h2_num_pair_dec[ic][ip][ifile] = new TH2D(*(TH2D*)fin->Get(name.c_str()));
	npair_dec[ic][ifile] += h2_num_pair_dec[ic][ip][ifile]->Integral();
	//cout<<"npair_dec["<<ic<<"]["<<ip<<"]["<<ifile<<"] = "<<npair_dec[ic][ip][ifile]<<endl;
      }
      npair_dec[ic][ifile]/=nevt;
      //cout<<"npair_dec["<<ic<<"]["<<ifile<<"] = "<<npair_dec[ic][ifile]<<endl;
      avg_ntrig_inc[ic]+=ntrig_inc[ic][ifile];
      avg_ntrig_pi0[ic]+=ntrig_pi0[ic][ifile];
      avg_ntrig_dec[ic]+=ntrig_dec[ic][ifile];
      avg_npair_inc[ic]+=npair_inc[ic][ifile];
      avg_npair_pi0[ic]+=npair_pi0[ic][ifile];
      avg_npair_dec[ic]+=npair_dec[ic][ifile];
    }
    
    fin->Close();
    delete fin;
  }

  //filling histograms to fit and find means.
  for(int ic=0; ic<2; ic++){
    avg_ntrig_inc[ic]/=N_RUNS;
    avg_ntrig_pi0[ic]/=N_RUNS;
    avg_ntrig_dec[ic]/=N_RUNS;
    avg_npair_inc[ic]/=N_RUNS;
    avg_npair_pi0[ic]/=N_RUNS;
    avg_npair_dec[ic]/=N_RUNS;
    // cout<<"avg_ntrig_inc["<<ic<<"] = "<<avg_ntrig_inc[ic]<<endl;
    // cout<<"avg_ntrig_pi0["<<ic<<"] = "<<avg_ntrig_pi0[ic]<<endl;
    // cout<<"avg_ntrig_dec["<<ic<<"] = "<<avg_ntrig_dec[ic]<<endl;
    // cout<<"avg_npair_inc["<<ic<<"] = "<<avg_npair_inc[ic]<<endl;
    // cout<<"avg_npair_pi0["<<ic<<"] = "<<avg_npair_pi0[ic]<<endl;
    // cout<<"avg_npair_dec["<<ic<<"] = "<<avg_npair_dec[ic]<<endl;

    bin.str("");
    bin << ic;
    name = "hNtrigInc_c" + bin.str();
    hNtrigInc[ic] = new TH1D(name.c_str(),name.c_str(),50,avg_ntrig_inc[ic]*0.8,avg_ntrig_inc[ic]*1.2);
    name = "hNtrigPi0_c" + bin.str();
    hNtrigPi0[ic] = new TH1D(name.c_str(),name.c_str(),50,avg_ntrig_pi0[ic]*0.8,avg_ntrig_pi0[ic]*1.2);
    name = "hNtrigDec_c" + bin.str();
    hNtrigDec[ic] = new TH1D(name.c_str(),name.c_str(),50,avg_ntrig_dec[ic]*0.8,avg_ntrig_dec[ic]*1.2);
    name = "hNpairInc_c" + bin.str();
    hNpairInc[ic] = new TH1D(name.c_str(),name.c_str(),50,avg_npair_inc[ic]*0.8,avg_npair_inc[ic]*1.2);
    name = "hNpairPi0_c" + bin.str();
    hNpairPi0[ic] = new TH1D(name.c_str(),name.c_str(),50,avg_npair_pi0[ic]*0.8,avg_npair_pi0[ic]*1.2);
    name = "hNpairDec_c" + bin.str();
    hNpairDec[ic] = new TH1D(name.c_str(),name.c_str(),50,avg_npair_dec[ic]*0.8,avg_npair_dec[ic]*1.2);

    for(int ifile=0; ifile<N_RUNS; ifile++){
      hNtrigInc[ic]->Fill(ntrig_inc[ic][ifile]);
      hNtrigPi0[ic]->Fill(ntrig_pi0[ic][ifile]);
      hNtrigDec[ic]->Fill(ntrig_dec[ic][ifile]);
      hNpairInc[ic]->Fill(npair_inc[ic][ifile]);
      hNpairPi0[ic]->Fill(npair_pi0[ic][ifile]);
      hNpairDec[ic]->Fill(npair_dec[ic][ifile]);
    }
    mean_ntrig_inc[ic] = hNtrigInc[ic]->GetMean();
    rms_ntrig_inc[ic] = hNtrigInc[ic]->GetRMS();
    mean_ntrig_pi0[ic] = hNtrigPi0[ic]->GetMean();
    rms_ntrig_pi0[ic] = hNtrigPi0[ic]->GetRMS();
    mean_ntrig_dec[ic] = hNtrigDec[ic]->GetMean();
    rms_ntrig_dec[ic] = hNtrigDec[ic]->GetRMS();
    mean_npair_inc[ic] = hNpairInc[ic]->GetMean();
    rms_npair_inc[ic] = hNpairInc[ic]->GetRMS();
    mean_npair_pi0[ic] = hNpairPi0[ic]->GetMean();
    rms_npair_pi0[ic] = hNpairPi0[ic]->GetRMS();
    mean_npair_dec[ic] = hNpairDec[ic]->GetMean();
    rms_npair_dec[ic] = hNpairDec[ic]->GetRMS();

  }

  cout<<"making graphs"<<endl;
  for(int ic=0; ic<2; ic++){
    bin.str("");
    bin << ic;
    gNtrigInc[ic] = new TGraph(N_RUNS,run_nbr,ntrig_inc[ic]);
    name = "gNtrigInc_c" + bin.str();
    gNtrigInc[ic]->SetName(name.c_str());
    
    name = "fNtrigIncMean_c" + bin.str();
    fNtrigIncMean[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    fNtrigIncMean[ic]->SetParameter(0,mean_ntrig_inc[ic]);
    name = "fNtrigIncRMS1_c" + bin.str();
    fNtrigIncRMS1[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNtrigIncRMS1[ic]->SetParameter(0,mean_ntrig_inc[ic]-3*rms_ntrig_inc[ic]);
    fNtrigIncRMS1[ic]->SetParameter(0,mean_ntrig_inc[ic]-1.5*rms_ntrig_inc[ic]);
    name = "fNtrigIncRMS2_c" + bin.str();
    fNtrigIncRMS2[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNtrigIncRMS2[ic]->SetParameter(0,mean_ntrig_inc[ic]+3*rms_ntrig_inc[ic]);
    fNtrigIncRMS2[ic]->SetParameter(0,mean_ntrig_inc[ic]+1.5*rms_ntrig_inc[ic]);

    gNtrigPi0[ic] = new TGraph(N_RUNS,run_nbr,ntrig_pi0[ic]);
    name = "gNtrigPi0_c" + bin.str();
    gNtrigPi0[ic]->SetName(name.c_str());
    
    name = "fNtrigPi0Mean_c" + bin.str();
    fNtrigPi0Mean[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    fNtrigPi0Mean[ic]->SetParameter(0,mean_ntrig_pi0[ic]);
    name = "fNtrigPi0RMS1_c" + bin.str();
    fNtrigPi0RMS1[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNtrigPi0RMS1[ic]->SetParameter(0,mean_ntrig_pi0[ic]-3*rms_ntrig_pi0[ic]);
    fNtrigPi0RMS1[ic]->SetParameter(0,mean_ntrig_pi0[ic]-1.5*rms_ntrig_pi0[ic]);
    name = "fNtrigPi0RMS2_c" + bin.str();
    fNtrigPi0RMS2[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNtrigPi0RMS2[ic]->SetParameter(0,mean_ntrig_pi0[ic]+3*rms_ntrig_pi0[ic]);
    fNtrigPi0RMS2[ic]->SetParameter(0,mean_ntrig_pi0[ic]+1.5*rms_ntrig_pi0[ic]);
    
    gNtrigDec[ic] = new TGraph(N_RUNS,run_nbr,ntrig_dec[ic]);
    name = "gNtrigDec_c" + bin.str();
    gNtrigDec[ic]->SetName(name.c_str());
    
    name = "fNtrigDecMean_c" + bin.str();
    fNtrigDecMean[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    fNtrigDecMean[ic]->SetParameter(0,mean_ntrig_dec[ic]);
    name = "fNtrigDecRMS1_c" + bin.str();
    fNtrigDecRMS1[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNtrigDecRMS1[ic]->SetParameter(0,mean_ntrig_dec[ic]-3*rms_ntrig_dec[ic]);
    fNtrigDecRMS1[ic]->SetParameter(0,mean_ntrig_dec[ic]-1.5*rms_ntrig_dec[ic]);
    name = "fNtrigDecRMS2_c" + bin.str();
    fNtrigDecRMS2[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNtrigDecRMS2[ic]->SetParameter(0,mean_ntrig_dec[ic]+3*rms_ntrig_dec[ic]);
    fNtrigDecRMS2[ic]->SetParameter(0,mean_ntrig_dec[ic]+1.5*rms_ntrig_dec[ic]);
    
    // lNtrigPi0Mean[ic] = new TLine(343786,mean_ntrig_dec[ic],349678,mean_ntrig_dec[ic]);
    // lNtrigPi0RMS1[ic] = new TLine(343786,mean_ntrig_dec[ic]+3*rms_ntrig_dec[ic],349678,mean_ntrig_dec[ic]+3*rms_ntrig_dec[ic]);
    // lNtrigPi0RMS2[ic] = new TLine(343786,mean_ntrig_dec[ic]-3*rms_ntrig_dec[ic],349678,mean_ntrig_dec[ic]-3*rms_ntrig_dec[ic]);

    gNpairInc[ic] = new TGraph(N_RUNS,run_nbr,npair_inc[ic]);
    name = "gNpairInc_c" + bin.str();
    gNpairInc[ic]->SetName(name.c_str());
    
    name = "fNpairIncMean_c" + bin.str();
    fNpairIncMean[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    fNpairIncMean[ic]->SetParameter(0,mean_npair_inc[ic]);
    name = "fNpairIncRMS1_c" + bin.str();
    fNpairIncRMS1[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNpairIncRMS1[ic]->SetParameter(0,mean_npair_inc[ic]-3*rms_npair_inc[ic]);
    fNpairIncRMS1[ic]->SetParameter(0,mean_npair_inc[ic]-1.5*rms_npair_inc[ic]);
    name = "fNpairIncRMS2_c" + bin.str();
    fNpairIncRMS2[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNpairIncRMS2[ic]->SetParameter(0,mean_npair_inc[ic]+3*rms_npair_inc[ic]);
    fNpairIncRMS2[ic]->SetParameter(0,mean_npair_inc[ic]+1.5*rms_npair_inc[ic]);
   
    gNpairPi0[ic] = new TGraph(N_RUNS,run_nbr,npair_pi0[ic]);
    name = "gNpairPi0_c" + bin.str();
    gNpairPi0[ic]->SetName(name.c_str());
    
    name = "fNpairPi0Mean_c" + bin.str();
    fNpairPi0Mean[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    fNpairPi0Mean[ic]->SetParameter(0,mean_npair_pi0[ic]);
    name = "fNpairPi0RMS1_c" + bin.str();
    fNpairPi0RMS1[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNpairPi0RMS1[ic]->SetParameter(0,mean_npair_pi0[ic]-3*rms_npair_pi0[ic]);
    fNpairPi0RMS1[ic]->SetParameter(0,mean_npair_pi0[ic]-1.5*rms_npair_pi0[ic]);
    name = "fNpairPi0RMS2_c" + bin.str();
    fNpairPi0RMS2[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNpairPi0RMS2[ic]->SetParameter(0,mean_npair_pi0[ic]+3*rms_npair_pi0[ic]);
    fNpairPi0RMS2[ic]->SetParameter(0,mean_npair_pi0[ic]+1.5*rms_npair_pi0[ic]);

    gNpairDec[ic] = new TGraph(N_RUNS,run_nbr,npair_dec[ic]);
    name = "gNpairDec_c" + bin.str();
    gNpairDec[ic]->SetName(name.c_str());
    
    name = "fNpairDecMean_c" + bin.str();
    fNpairDecMean[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    fNpairDecMean[ic]->SetParameter(0,mean_npair_dec[ic]);
    name = "fNpairDecRMS1_c" + bin.str();
    fNpairDecRMS1[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNpairDecRMS1[ic]->SetParameter(0,mean_npair_dec[ic]-3*rms_npair_dec[ic]);
    fNpairDecRMS1[ic]->SetParameter(0,mean_npair_dec[ic]-1.5*rms_npair_dec[ic]);
    name = "fNpairDecRMS2_c" + bin.str();
    fNpairDecRMS2[ic] = new TF1(name.c_str(),"[0]",343786.,349678.);
    //fNpairDecRMS2[ic]->SetParameter(0,mean_npair_dec[ic]+3*rms_npair_dec[ic]);
    fNpairDecRMS2[ic]->SetParameter(0,mean_npair_dec[ic]+1.5*rms_npair_dec[ic]);

    foutput->cd();
    gNtrigInc[ic]->Write();
    gNtrigPi0[ic]->Write();
    gNtrigDec[ic]->Write();
    gNpairInc[ic]->Write();
    gNpairPi0[ic]->Write();
    gNpairDec[ic]->Write();
    hNtrigInc[ic]->Write();
    hNtrigPi0[ic]->Write();
    hNtrigDec[ic]->Write();
    hNpairInc[ic]->Write();
    hNpairPi0[ic]->Write();
    hNpairDec[ic]->Write();
    
    fNtrigIncMean[ic]->Write();
    fNtrigIncRMS1[ic]->Write();
    fNtrigIncRMS2[ic]->Write();
    fNtrigPi0Mean[ic]->Write();
    fNtrigPi0RMS1[ic]->Write();
    fNtrigPi0RMS2[ic]->Write();
    fNtrigDecMean[ic]->Write();
    fNtrigDecRMS1[ic]->Write();
    fNtrigDecRMS2[ic]->Write();
    fNpairIncMean[ic]->Write();
    fNpairIncRMS1[ic]->Write();
    fNpairIncRMS2[ic]->Write();
    fNpairPi0Mean[ic]->Write();
    fNpairPi0RMS1[ic]->Write();
    fNpairPi0RMS2[ic]->Write();
    fNpairDecMean[ic]->Write();
    fNpairDecRMS1[ic]->Write();
    fNpairDecRMS2[ic]->Write();
    
  }
  return 0;
}

//plot # of triggers and partners as a function of run #, centrality integrated
int triggers_run_by_run_v3(const char* flist, const string fout)
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
  ostringstream bin;
  string name;

  TH1D* h1_zvertex[N_RUNS];
  TH1D* h1_num_trig[N_RUNS];
  TH1D* h1_num_part[N_RUNS];

  TH1D* hNtrig;
  TH1D* hNpart;
  TH1D* hNpart_rg1;
  TH1D* hNpart_rg2;
  TH1D* hNpart_rg3;
  //TH1D* hNpart_rg4;

  TGraph* gNtrig;
  TGraph* gNpart; 

  TF1* fNtrigMean;
  TF1* fNtrigRMS1;
  TF1* fNtrigRMS2;

  TF1* fNpartMean;
  TF1* fNpartRMS1;
  TF1* fNpartRMS2;
  TF1* fNpartMeanRG1_1;
  TF1* fNpartRMS1RG1_1;
  TF1* fNpartRMS2RG1_1;
  TF1* fNpartMeanRG1_2;
  TF1* fNpartRMS1RG1_2;
  TF1* fNpartRMS2RG1_2;
  TF1* fNpartMeanRG2_1;
  TF1* fNpartRMS1RG2_1;
  TF1* fNpartRMS2RG2_1;
  TF1* fNpartMeanRG2_2;
  TF1* fNpartRMS1RG2_2;
  TF1* fNpartRMS2RG2_2;
  TF1* fNpartMeanRG3;
  TF1* fNpartRMS1RG3;
  TF1* fNpartRMS2RG3;
  
  double ntrig[N_RUNS];
  double npart[N_RUNS];

  double mean_ntrig;
  double rms_ntrig;
  double mean_npart;
  double rms_npart;
  double mean_npart_rg1;
  double rms_npart_rg1;
  double mean_npart_rg2;
  double rms_npart_rg2;
  double mean_npart_rg3;
  double rms_npart_rg3;

  double avg_ntrig = 0.0;
  double avg_npart = 0.0;

  TFile* foutput = new TFile (fout.c_str(),"recreate");
  for(int ifile = 0; ifile < N_RUNS; ifile++){
    file_string = FileList[ifile];
    cout<<"file_string: "<<file_string<<endl;
    int p1 = file_string.find_last_of(".");
    run_string[ifile] = file_string.substr(p1-6,6);
   
    run_nbr[ifile] = atoi(run_string[ifile].c_str());
    cout<<"Run_nbr "<<ifile<<" is "<<run_nbr[ifile]<<endl;

    TFile* fin = new TFile(file_string.c_str());
    if (! fin->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }

    h1_zvertex[ifile] = new TH1D(*(TH1D*)fin->Get("h1_zvertex"));
    double nevt = h1_zvertex[ifile]->Integral();
    cout<<"h1_zvertex["<<ifile<<"] = "<<nevt<<endl;

    name = "h1_trig_pt_inc_tot";
    h1_num_trig[ifile] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
    ntrig[ifile] = h1_num_trig[ifile]->Integral();
    ntrig[ifile]/=nevt;
    cout<<"ntrig["<<ifile<<"] = "<<ntrig[ifile]<<endl;

    name = "h1_part_pt_tot";
    h1_num_part[ifile] = new TH1D(*(TH1D*)fin->Get(name.c_str()));
    npart[ifile] = h1_num_part[ifile]->Integral();
    npart[ifile]/=nevt;
    cout<<"npart["<<ifile<<"] = "<<npart[ifile]<<endl;
    
    avg_ntrig+=ntrig[ifile];
    avg_npart+=npart[ifile];
  }
  fin->Close();
  delete fin;

  //filling Npart histograms according tonrun groups
  // hNpart_rg1 = new TH1D("hNpart_rg1","hNpart_rg1",50,8.5,10.5);//run group 1: up to 345875, 152 runs
  // hNpart_rg2 = new TH1D("hNpart_rg2","hNpart_rg2",50,7.0,10.5);//run group 2: 346009-348875, 320 runs
  // hNpart_rg3 = new TH1D("hNpart_rg3","hNpart_rg3",50,8.0,10.5);//run group 3: 348984 and above 66 runs

  //up to 343966(19 runs), 344126-344709(55 runs), 344894-346670(165 runs), 346973-348875(233 runs), 348984 and above(66 runs)
  hNpart_rg1 = new TH1D("hNpart_rg1","hNpart_rg1",50,7.0,10.5);
  hNpart_rg2 = new TH1D("hNpart_rg2","hNpart_rg2",50,7.0,10.5);
  hNpart_rg3 = new TH1D("hNpart_rg3","hNpart_rg3",50,7.0,10.5);

  for(int ifile=0; ifile<19; ifile++) hNpart_rg1->Fill(npart[ifile]);
  for(int ifile=19; ifile<74; ifile++) hNpart_rg2->Fill(npart[ifile]);
  for(int ifile=74; ifile<239; ifile++) hNpart_rg1->Fill(npart[ifile]);
  for(int ifile=239; ifile<472; ifile++) hNpart_rg3->Fill(npart[ifile]);
  for(int ifile=472; ifile<538; ifile++) hNpart_rg2->Fill(npart[ifile]);
  mean_npart_rg1 = hNpart_rg1->GetMean();
  rms_npart_rg1 = hNpart_rg1->GetRMS();
  mean_npart_rg2 = hNpart_rg2->GetMean();
  rms_npart_rg2 = hNpart_rg2->GetRMS();
  mean_npart_rg3 = hNpart_rg3->GetMean();
  rms_npart_rg3 = hNpart_rg3->GetRMS();

  fNpartMeanRG1_1 = new TF1("fNpartMeanRG1_1","[0]",343786.,343966.);
  fNpartMeanRG1_1->SetParameter(0,mean_npart_rg1);
  fNpartRMS1RG1_1 = new TF1("fNpartRMS1RG1_1","[0]",343786.,343966.);
  fNpartRMS1RG1_1->SetParameter(0,mean_npart_rg1-1.5*rms_npart_rg1);
  fNpartRMS2RG1_1 = new TF1("fNpartRMS2RG1_1","[0]",343786.,343966.);
  fNpartRMS2RG1_1->SetParameter(0,mean_npart_rg1+1.5*rms_npart_rg1);

  fNpartMeanRG1_2 = new TF1("fNpartMeanRG1_2","[0]",344894.,346670.);
  fNpartMeanRG1_2->SetParameter(0,mean_npart_rg1);
  fNpartRMS1RG1_2 = new TF1("fNpartRMS1RG1_2","[0]",344894.,346670.);
  fNpartRMS1RG1_2->SetParameter(0,mean_npart_rg1-1.5*rms_npart_rg1);
  fNpartRMS2RG1_2 = new TF1("fNpartRMS2RG1_2","[0]",344894.,346670.);
  fNpartRMS2RG1_2->SetParameter(0,mean_npart_rg1+1.5*rms_npart_rg1);

  fNpartMeanRG2_1 = new TF1("fNpartMeanRG2_1","[0]",344126.,344709.);
  fNpartMeanRG2_1->SetParameter(0,mean_npart_rg2);
  fNpartRMS1RG2_1 = new TF1("fNpartRMS1RG2_1","[0]",344126.,344709.);
  fNpartRMS1RG2_1->SetParameter(0,mean_npart_rg2-1.5*rms_npart_rg2);
  fNpartRMS2RG2_1 = new TF1("fNpartRMS2RG2_1","[0]",344126.,344709.);
  fNpartRMS2RG2_1->SetParameter(0,mean_npart_rg2+1.5*rms_npart_rg2);

  fNpartMeanRG2_2 = new TF1("fNpartMeanRG2_2","[0]",348984.,349678.);
  fNpartMeanRG2_2->SetParameter(0,mean_npart_rg2);
  fNpartRMS1RG2_2 = new TF1("fNpartRMS1RG2_2","[0]",348984.,349678.);
  fNpartRMS1RG2_2->SetParameter(0,mean_npart_rg2-1.5*rms_npart_rg2);
  fNpartRMS2RG2_2 = new TF1("fNpartRMS2RG2_2","[0]",348984.,349678.);
  fNpartRMS2RG2_2->SetParameter(0,mean_npart_rg2+1.5*rms_npart_rg2);

  fNpartMeanRG3 = new TF1("fNpartMeanRG3","[0]",346973.,348875.);
  fNpartMeanRG3->SetParameter(0,mean_npart_rg3);
  fNpartRMS1RG3 = new TF1("fNpartRMS1RG3","[0]",346973.,348875.);
  fNpartRMS1RG3->SetParameter(0,mean_npart_rg3-1.5*rms_npart_rg3);
  fNpartRMS2RG3 = new TF1("fNpartRMS2RG3","[0]",346973.,348875.);
  fNpartRMS2RG3->SetParameter(0,mean_npart_rg3+1.5*rms_npart_rg3);

  //filling histograms to fit and find means.
  avg_ntrig/=N_RUNS;
  avg_npart/=N_RUNS;
  name = "hNtrig";
  hNtrig = new TH1D(name.c_str(),name.c_str(),50,avg_ntrig*0.8,avg_ntrig*1.2);
  name = "hNpart";
  hNpart = new TH1D(name.c_str(),name.c_str(),50,avg_npart*0.8,avg_npart*1.2);

  for(int ifile=0; ifile<N_RUNS; ifile++){
    hNtrig->Fill(ntrig[ifile]);
    hNpart->Fill(npart[ifile]);
  }
  mean_ntrig = hNtrig->GetMean();
  rms_ntrig = hNtrig->GetRMS();
  mean_npart = hNpart->GetMean();
  rms_npart = hNpart->GetRMS();

  cout<<"making graphs"<<endl;
  gNtrig = new TGraph(N_RUNS,run_nbr,ntrig);
  name = "gNtrig";
  gNtrig->SetName(name.c_str());
  gNpart = new TGraph(N_RUNS,run_nbr,npart);
  name = "gNpart";
  gNpart->SetName(name.c_str());
    
  name = "fNtrigMean";
  fNtrigMean = new TF1(name.c_str(),"[0]",343786.,349678.);
  fNtrigMean->SetParameter(0,mean_ntrig);
  name = "fNtrigRMS1";
  fNtrigRMS1 = new TF1(name.c_str(),"[0]",343786.,349678.);
  fNtrigRMS1->SetParameter(0,mean_ntrig-1.5*rms_ntrig);
  name = "fNtrigRMS2";
  fNtrigRMS2 = new TF1(name.c_str(),"[0]",343786.,349678.);
  fNtrigRMS2->SetParameter(0,mean_ntrig+1.5*rms_ntrig);

  name = "fNpartMean";
  fNpartMean = new TF1(name.c_str(),"[0]",343786.,349678.);
  fNpartMean->SetParameter(0,mean_npart);
  name = "fNpartRMS1";
  fNpartRMS1 = new TF1(name.c_str(),"[0]",343786.,349678.);
  fNpartRMS1->SetParameter(0,mean_npart-1.5*rms_npart);
  name = "fNpartRMS2";
  fNpartRMS2 = new TF1(name.c_str(),"[0]",343786.,349678.);
  fNpartRMS2->SetParameter(0,mean_npart+1.5*rms_npart);

  foutput->cd();
  gNtrig->Write();
  hNtrig->Write();
  fNtrigMean->Write();
  fNtrigRMS1->Write();
  fNtrigRMS2->Write();

  fNpartMeanRG1_1->Write();
  fNpartRMS1RG1_1->Write();
  fNpartRMS2RG1_1->Write();
  fNpartMeanRG1_2->Write();
  fNpartRMS1RG1_2->Write();
  fNpartRMS2RG1_2->Write();
  fNpartMeanRG2_1->Write();
  fNpartRMS1RG2_1->Write();
  fNpartRMS2RG2_1->Write();
  fNpartMeanRG2_2->Write();
  fNpartRMS1RG2_2->Write();
  fNpartRMS2RG2_2->Write();
  fNpartMeanRG3->Write();
  fNpartRMS1RG3->Write();
  fNpartRMS2RG3->Write();

  gNpart->Write();
  hNpart->Write();
  hNpart_rg1->Write();
  hNpart_rg2->Write();
  hNpart_rg3->Write();
  //hNpart_rg4->Write();
  fNpartMean->Write();
  fNpartRMS1->Write();
  fNpartRMS2->Write();
  
  return 0;
}
