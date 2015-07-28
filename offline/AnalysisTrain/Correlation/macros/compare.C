#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>

using namespace std;

const double PI = TMath::Pi();

void DoTtest(double n1, double n2, double mean1, double mean2, double var1, double var2, double &t, double& prob)
{
  double sig1 = var1*var1/n1;
  double sig2 = var2*var2/n2;
  t = (mean1-mean2)/sqrt(sig1 + sig2);
  std::cout << "t-stat = " << t << std::endl;

  TF1* gaus = new TF1("gaus","TMath::Gaus(x,[0],[1],1)",-20,20);
  gaus->SetParameter(0,0);
  gaus->SetParameter(1,1);
  double integral = gaus->Integral(-1*fabs(t),fabs(t));
  prob = 1-integral;
  //prob = prob/(1-1/pow(df,2));
  std::cout << "prob = " << prob << std::endl;
}

double GetChiSquare(double mean1, double mean2, double var1, double var2)
{
  double chi2 = pow((mean1 - mean2),2)/(var1*var1 + var2*var2);
  return chi2;
}

double GetChiSquareProb(double ndf, double chi2)
{
  TF1* func = new TF1("func","1/(pow(2,[0]/2.0)*TMath::Gamma([0]/2.0))*pow(x,[0]/2.0-1.0)*TMath::Exp(-x/2.0)",0,1000);//pdf of chi-squared distribution
  func->SetParameter(0,ndf);
  double integral = func->Integral(0,chi2);
  double prob = 1-integral;
  std::cout << "chi2/df = " << chi2/ndf << " with p-value = " << prob << std::endl;
  return prob;
}

double Round(double number)
{
  number*=100;
  number+=0.5;
  number=(int)number;
  number/=100;
  return number;
}

int compare(int centbin, int type, int isjf, const string frun10, const string fout)//frun10:"run10_output_old.lst"
{
  string frun11; 

  if(isjf){
    if(type == 0) frun11 = "/direct/phenix+u/hge/offline/AnalysisTrain/Correlation/macros/makejfs_inc_fold_taxi5286.root";
    if(type == 1) frun11 = "/direct/phenix+u/hge/offline/AnalysisTrain/Correlation/macros/makejfs_pi0_fold_taxi5286.root";
    if(type == 2) frun11 = "/direct/phenix+u/hge/offline/AnalysisTrain/Correlation/macros/makejfs_dec_fold_taxi5286.root";
    if(type == 3) frun11 = "/direct/phenix+u/hge/offline/AnalysisTrain/Correlation/macros/makedirs_taxi5286.root";
  }
  else{
    if(type == 0) frun11 = "/direct/phenix+u/hge/offline/AnalysisTrain/Correlation/macros/makecfs_inc_fold_taxi5286.root";
    if(type == 1) frun11 = "/direct/phenix+u/hge/offline/AnalysisTrain/Correlation/macros/makecfs_pi0_fold_taxi5286.root";
    if(type == 2) frun11 = "/direct/phenix+u/hge/offline/AnalysisTrain/Correlation/macros/makecfs_dec_fold_taxi5286.root";
  }

  TFile* f11 = new TFile(frun11.c_str());
  
  vector<string> Run10FileList;
  char buffer[256];
  cout<<"opening run10 input file list: " << frun10.c_str() << endl;
  ifstream flist(frun10.c_str());
  if( !flist.is_open()) return -1;

  while (flist.getline(buffer, 256, '\n')){
    cout<<"in the while loop"<<endl;
    char this_file[512];
    strcpy(this_file, buffer);
    cout<<"this_file: "<<this_file<<endl;
    string s (this_file);
    cout<<"s: "<<s.c_str()<<endl;
    Run10FileList.push_back(s);
  }
  cout << "Done making Run10FileList vector." << endl;
  const int Nfiles = Run10FileList.size();
  cout<<"Nfiles = "<<Nfiles<<endl;

  TFile* foutput = new TFile(fout.c_str(), "recreate");
  foutput->cd();

  TH1D *my_corr;
  TH1D *corr;
  TH1D *ratio;
  string file_string;
  string photon_bin_string;
  string hadron_bin_string;
  int photon_bin;
  int hadron_bin;
  ostringstream bin;
  ostringstream cbin;
  ostringstream cfname;
  ostringstream jfname;
  ostringstream dirname;
  ostringstream can_name;
  ostringstream latex_name;
  string ratio_name;

  cout<<"define canvases."<<endl;
  TCanvas *can[4];
  TCanvas *can_ratio[4];
  for(int i=0; i<4; i++){
    can_name.str("");
    can_name << "can_p"<<i;
    can[i] = new TCanvas(can_name.str().c_str(),can_name.str().c_str());
    can[i]->Divide(3,2);
    can_name.str("");
    can_name << "can_ratio_p"<<i;
    can_ratio[i] = new TCanvas(can_name.str().c_str(),can_name.str().c_str());
    can_ratio[i]->Divide(3,2);
  }

  TH1D* ttest_CF = new TH1D("ttest_CF","ttest_CF",50,-5,5);
  ttest_CF->GetXaxis()->SetTitle("t");
  TH1D* tprob_CF = new TH1D("tprob_CF","tprob_CF",50,0,1);
  tprob_CF->GetXaxis()->SetTitle("prob");
  TH1D* chitest_CF = new TH1D("chitest_CF","chitest_CF",50,0,5);
  chitest_CF->GetXaxis()->SetTitle("#chi^{2}/NDF");
  TH1D* chiprob_CF = new TH1D("chiprob_CF","chiprob_CF",50,0,1);
  chiprob_CF->GetXaxis()->SetTitle("prob");

  cout<<"start looping run10 files."<<endl;
  for(int ifile = 0; ifile < Nfiles; ifile++){
    cout<<"Run10FileList ["<<ifile<<"] = "<<Run10FileList[ifile]<<endl;
    file_string = Run10FileList[ifile];
    int p1 = file_string.find_last_of("/");
    cout<<"p1 = "<<p1<<endl;
    photon_bin_string = file_string.substr(p1+6,1);
    photon_bin = atoi(photon_bin_string.c_str());
    int p2 = file_string.find_last_of(".");
    cout<<"p2 = "<<p2<<endl;
    hadron_bin_string = file_string.substr(p2-1,1);
    hadron_bin = atoi(hadron_bin_string.c_str());

    cout<<"photon_bin: "<<photon_bin<<" hadron_bin: "<<hadron_bin<<endl;
    cbin.str("");
    cbin << centbin;
    cfname.str("");
    jfname.str("");
    dirname.str("");
    can_name.str("");
    switch(photon_bin){
    case 5:
      cout<<"photon case 5"<<endl;
      latex_name.str("");
      latex_name << "5-7 #times ";
      switch(hadron_bin){
      case 1:
	latex_name << "0.5-1 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p0_h0";
	jfname << "JF_c"<<cbin.str()<<"_p0_h0";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p0_h0";
	can[0]->cd(1);
	break;
      case 2:
	latex_name << "1-2 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p0_h1";
	jfname << "JF_c"<<cbin.str()<<"_p0_h1";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p0_h1";
	can[0]->cd(2);
	break;
      case 3:
	latex_name << "2-3 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p0_h2";
	jfname << "JF_c"<<cbin.str()<<"_p0_h2";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p0_h2";
	can[0]->cd(3);
	break;
      case 5:
	latex_name << "3-5 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p0_h3";
	jfname << "JF_c"<<cbin.str()<<"_p0_h3";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p0_h3";
	can[0]->cd(4);
	break;
	
      case 0:
	latex_name << "5-7 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p0_h4";
	jfname << "JF_c"<<cbin.str()<<"_p0_h4";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p0_h4";
	can[0]->cd(5);
	break;
	
      }
      break;

    case 7:
      cout<<"photon case 7"<<endl;
      latex_name.str("");
      latex_name << "7-9 #times ";
      switch(hadron_bin){
      case 1:
	latex_name << "0.5-1 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p1_h0";
	jfname << "JF_c"<<cbin.str()<<"_p1_h0";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p1_h0";
	can[1]->cd(1);
	break;
      case 2:
	latex_name << "1-2 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p1_h1";
	jfname << "JF_c"<<cbin.str()<<"_p1_h1";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p1_h1";
	can[1]->cd(2);
	break;
      case 3:
	latex_name << "2-3 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p1_h2";
	jfname << "JF_c"<<cbin.str()<<"_p1_h2";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p1_h2";
	can[1]->cd(3);
	break;
      case 5:
	latex_name << "3-5 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p1_h3";
	jfname << "JF_c"<<cbin.str()<<"_p1_h3";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p1_h3";
	can[1]->cd(4);
	break;
	
      case 0:
	latex_name << "5-7 GeV/c";
        cfname << "CF_c"<<cbin.str()<<"_p1_h4";
	jfname << "JF_c"<<cbin.str()<<"_p1_h4";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p1_h4";
	can[1]->cd(5);
	break;
	
      }
      break;

    case 9:
      cout<<"photon case 9"<<endl;
      latex_name.str("");
      latex_name << "9-12 #times ";
      switch(hadron_bin){
      case 1:
	latex_name << "0.5-1 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p2_h0";
	jfname << "JF_c"<<cbin.str()<<"_p2_h0";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p2_h0";
	can[2]->cd(1);
	break;
      case 2:
	latex_name << "1-2 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p2_h1";
	jfname << "JF_c"<<cbin.str()<<"_p2_h1";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p2_h1";
	can[2]->cd(2);
	break;
      case 3:
	latex_name << "2-3 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p2_h2";
	jfname << "JF_c"<<cbin.str()<<"_p2_h2";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p2_h2";
	can[2]->cd(3);
	break;
      case 5:
	latex_name << "3-5 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p2_h3";
	jfname << "JF_c"<<cbin.str()<<"_p2_h3";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p2_h3";
	can[2]->cd(4);
	break;
	
      case 0:
	latex_name << "5-7 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p2_h4";
	jfname << "JF_c"<<cbin.str()<<"_p2_h4";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p2_h4";
	can[2]->cd(5);
	break;
	
      }
      break;

    case 1:
      cout<<"photon case 12"<<endl;
      latex_name.str("");
      latex_name << "12-15 #times ";
      switch(hadron_bin){
      case 1:
	latex_name << "0.5-1 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p3_h0";
	jfname << "JF_c"<<cbin.str()<<"_p3_h0";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p3_h0";
	can[3]->cd(1);
	break;
      case 2:
	latex_name << "1-2 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p3_h1";
	jfname << "JF_c"<<cbin.str()<<"_p3_h1";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p3_h1";
	can[3]->cd(2);
	break;
      case 3:
	latex_name << "2-3 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p3_h2";
	jfname << "JF_c"<<cbin.str()<<"_p3_h2";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p3_h2";
	can[3]->cd(3);
	break;
      case 5:
	latex_name << "3-5 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p3_h3";
	jfname << "JF_c"<<cbin.str()<<"_p3_h3";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p3_h3";
	can[3]->cd(4);
	break;
	
      case 0:
	latex_name << "5-7 GeV/c";
	cfname << "CF_c"<<cbin.str()<<"_p3_h4";
	jfname << "JF_c"<<cbin.str()<<"_p3_h4";
	dirname << "DIR_JF_c"<<cbin.str()<<"_p3_h4";
        can[3]->cd(5);
	break;
	
      }
      break;
    default:
      cout << "no such photon bin!" << endl;
    }

    TFile* file = new TFile(file_string.c_str());
    if (! file->IsOpen() ) {
      cout << "cannot open file '" << file_string << "'\n";
      continue;
    }
    cout<<"f11 CF: "<<cfname.str().c_str()<<endl;
    cout<<"f11 JF: "<<jfname.str().c_str()<<endl;
    cout<<"f11 DIR: "<<dirname.str().c_str()<<endl;
    cout<<"f10: "<<file_string.c_str()<<endl;

    if(isjf){
      if(type==3) my_corr = new TH1D (*(TH1D*)f11->Get(dirname.str().c_str()));
      else my_corr = new TH1D (*(TH1D*)f11->Get(jfname.str().c_str()));
      if(type == 0) corr = new TH1D (*(TH1D*)file->Get("JF_INC"));
      if(type == 1) corr = new TH1D (*(TH1D*)file->Get("JF_PI0"));
      if(type == 2) corr = new TH1D (*(TH1D*)file->Get("JF_DEC"));
      if(type == 3) corr = new TH1D (*(TH1D*)file->Get("JF_DIR"));
    }
    else{
       my_corr = new TH1D (*(TH1D*)f11->Get(cfname.str().c_str()));
      if(type == 0) corr = new TH1D (*(TH1D*)file->Get("hDPHI_INC"));
      if(type == 1) corr = new TH1D (*(TH1D*)file->Get("hDPHI_PI0"));
      if(type == 2) corr = new TH1D (*(TH1D*)file->Get("hDPHI_DEC"));
    }
    
    int nReb;
    if(photon_bin == 5) nReb = 3;
    else if(photon_bin == 7) nReb = 5;
    else nReb = 6;
    corr->Rebin(nReb);
    my_corr->Rebin(nReb);
    
    // if(!isjf){
      double my_R = my_corr->Integral("width")/PI;
      cout<<"my_R = "<<my_R<<endl;
      my_corr->Scale(1./my_R);
      double R = corr->Integral("width")/PI;
      cout<<"R = "<<R<<endl;
      //corr->Sumw2();
      cout<<"my_R/R = "<<my_R/R<<endl;
      corr->Scale(1./R);
      //   }
    
    corr->SetMarkerStyle(20);
    corr->SetMarkerSize(0.8);
    my_corr->SetMarkerColor(1);
    my_corr->SetLineColor(1);
    if(isjf){
      double ymin = corr->GetMinimum();
      double ymax = corr->GetMaximum();
      if(my_corr->GetMaximum() > ymax) ymax = my_corr->GetMaximum();
      if(my_corr->GetMinimum() < ymin) ymin = my_corr->GetMinimum();
      corr->SetAxisRange(ymin-0.01, 1.2*ymax,"Y");
    }
    corr->GetYaxis()->SetTitleOffset(1.9);
    corr->Draw();
    my_corr->Draw("same");

    TLegend *l1 = new TLegend(0.7,0.7,0.9,0.9,"","brNDC");
    l1->AddEntry(my_corr,"Run11","lpf");
    l1->AddEntry(corr,"Run10","lpf");
    if(hadron_bin==0) l1->Draw("same");

    TLatex *la = new TLatex(0.3, 0.75, latex_name.str().c_str());
    la->SetNDC();
    la->Draw("same"); 
      
    string pbin_string = cfname.str().substr(7,1);
    string hbin_string = cfname.str().substr(10,1);
    if(pbin_string.compare("0") == 0){
      if(hbin_string.compare("0") == 0) can_ratio[0]->cd(1);
      if(hbin_string.compare("1") == 0) can_ratio[0]->cd(2);
      if(hbin_string.compare("2") == 0) can_ratio[0]->cd(3);
      if(hbin_string.compare("3") == 0) can_ratio[0]->cd(4);
      if(hbin_string.compare("4") == 0) can_ratio[0]->cd(5);
    }
    if(pbin_string.compare("1") == 0){
      if(hbin_string.compare("0") == 0) can_ratio[1]->cd(1);
      if(hbin_string.compare("1") == 0) can_ratio[1]->cd(2);
      if(hbin_string.compare("2") == 0) can_ratio[1]->cd(3);
      if(hbin_string.compare("3") == 0) can_ratio[1]->cd(4);
      if(hbin_string.compare("4") == 0) can_ratio[1]->cd(5);
    }
    if(pbin_string.compare("2") == 0){
      if(hbin_string.compare("0") == 0) can_ratio[2]->cd(1);
      if(hbin_string.compare("1") == 0) can_ratio[2]->cd(2);
      if(hbin_string.compare("2") == 0) can_ratio[2]->cd(3);
      if(hbin_string.compare("3") == 0) can_ratio[2]->cd(4);
      if(hbin_string.compare("4") == 0) can_ratio[2]->cd(5);
    }
    if(pbin_string.compare("3") == 0){
      if(hbin_string.compare("0") == 0) can_ratio[3]->cd(1);
      if(hbin_string.compare("1") == 0) can_ratio[3]->cd(2);
      if(hbin_string.compare("2") == 0) can_ratio[3]->cd(3);
      if(hbin_string.compare("3") == 0) can_ratio[3]->cd(4);
      if(hbin_string.compare("4") == 0) can_ratio[3]->cd(5);
    }
    if(isjf) ratio_name = "ratio_" + jfname.str();
    else ratio_name = "ratio_" + cfname.str();
    ratio = new TH1D (*(TH1D*)my_corr);
    ratio->SetName(ratio_name.c_str());
    cout<<"before dividing"<<endl;
    ratio->Divide(corr);
    if(isjf) ratio->SetAxisRange(-3.0, 5.0, "Y");
    ratio->Draw();
    la->Draw("same"); 

    if(isjf){ 
      ratio_name = "ratio_fit_" + jfname.str();
      /*  
      ratio_fit = new TF1(ratio_name.c_str(),"pol0(0)");
      ratio->Fit(ratio_name.c_str(),"BQ");
      double fitconst = Round(ratio_fit->GetParameter(0));
      bin.str("");
      bin << "ratio #approx "<<fitconst;
      TLatex *la_const = new TLatex(0.6, 0.25, bin.str().c_str());
      la_const->SetNDC();
      la_const->Draw("same");
      double fit_chi2 = ratio_fit->GetChisquare();
      int fit_ndf = ratio_fit->GetNDF();
      double chi2 = Round(fit_chi2/fit_ndf);      
      bin.str("");
      bin << "#chi^{2} = "<<chi2;
      TLatex *la_chi2 = new TLatex(0.6, 0.2, bin.str().c_str());
      la_chi2->SetNDC();
      la_chi2->Draw("same");
      */
    }
       
    //        if(!isjf){         
      //Do a chisquare test for CF comparison.
    double prob = 0; double t = 0; double chi2 = 0;
    double nbins = (double)corr->GetNbinsX();
    double ncounts1 = my_corr->GetEntries();
    double ncounts2 = corr->GetEntries();
    double var1 = 0; double var2 = 0;
    double mean1 = 0; double mean2 = 0;
    double sumw7 = 0;  double sumw10 = 0;
    
    my_corr->IntegralAndError(1,nbins,sumw7);
    my_corr->IntegralAndError(1,nbins,sumw10);
    for( int ib = 1; ib <= nbins; ib++ )
      {
	double err1 = my_corr->GetBinError(ib);
	mean1 += my_corr->GetBinContent(ib);
	double err2 = corr->GetBinError(ib);
	mean2 += corr->GetBinContent(ib);
	chi2 += GetChiSquare(my_corr->GetBinContent(ib),corr->GetBinContent(ib),err1,err2);
	var1 += err1*err1;
	var2 += err2*err2;
      }
    std::cout << std::endl;
    std::cout << "t-test for " << ncounts1 << "," << ncounts2 << " counts with means of " << mean1/nbins << "," << mean2/nbins << " and uncertainties of " <<  sqrt(var1) << "," << sqrt(var2) << std::endl;
    //var1 += pow(JF_INC_run7->GetMeanError(),2);
    //var2 += pow(JF_INC_run10->GetMeanError(),2);

    DoTtest(nbins,nbins,mean1/nbins,mean2/nbins,sqrt(var1),sqrt(var2),t,prob);
    ttest_CF->Fill(t);
    tprob_CF->Fill(prob);

    double cprob = GetChiSquareProb(nbins,chi2);
    double reducedchi2 = chi2/nbins;
    cout<<"reducedchi2 = "<<reducedchi2<<endl;

    reducedchi2 = Round(reducedchi2);
    cout<<reducedchi2<<endl;
    cprob = Round(cprob);

    latex_name.str("");
    latex_name << "#chi^{2} = "<<reducedchi2;
    TLatex *la_chi2 = new TLatex(0.6, 0.2, latex_name.str().c_str());
    la_chi2->SetNDC();
    la_chi2->Draw("same");
    latex_name.str("");
    latex_name << "p = "<<cprob;
    TLatex *la_cprob = new TLatex(0.6, 0.15, latex_name.str().c_str());
    la_cprob->SetNDC();
    la_cprob->Draw("same");
    
    chitest_CF->Fill(chi2/nbins);
    chiprob_CF->Fill(cprob);
    //        }
    
  }
  
  
 
  foutput->cd();
  for(int i=0; i<4; i++){
    can[i]->Write();
    can_ratio[i]->Write();
  }
    
  if(!isjf){
    ttest_CF->Write();
    tprob_CF->Write();
    chitest_CF->Write();
    chiprob_CF->Write();
  }
  
  
  cout<<"finish writing out."<<endl;
  return 0;
}

//for 0-40% combined stats comparison
void compare_combined(string f11, string f10, string fout)
{
  ostringstream bin;
  string name;
  TH1D* h10;
  TH1D* h11;
  TFile* fin10 = new TFile(f10.c_str());
  TFile* fin11 = new TFile(f11.c_str());
  TFile* foutput = new TFile(fout.c_str(),"recreate");
  foutput->cd();

  TCanvas* c1[2];//for 5-9,9-15 trig pt bins.
  for(int itrig=0; itrig<2; itrig++){
    bin.str("");
    bin << itrig;
    name = "can_p" + bin.str();
    cout<<"can name = "<<name.c_str()<<endl;
    c1[itrig] = new TCanvas(name.c_str(),name.c_str());
    c1[itrig]->Divide(3,2);

    for(int ipart=0; ipart<5; ipart++){
      bin.str("");
      bin << itrig <<"_h"<<ipart;
      name = "JFrun7_COMBTOT_c0_p"+bin.str();
      h10 = new TH1D(*(TH1D*)fin10->Get(name.c_str()));
      name = "JF_COMBTOT_c0_p"+bin.str();
      h11 = new TH1D(*(TH1D*)fin11->Get(name.c_str()));
      c1[itrig]->cd(ipart+1);
      h10->SetMarkerColor(2);
      h10->SetLineColor(2);
      h10->Rebin(3);
      double R10 = h10->Integral("width")/PI;
      //h10->Scale(1/R10);
      h10->Draw();
      h11->Rebin(3);
      double R11 = h11->Integral("width")/PI;
      h11->Scale(1/0.8);
      h11->Draw("same");
    }
    c1[itrig]->Write();
  }
}

 

