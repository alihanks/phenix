#include "TH1.h"
#include "TH3F.h"
#include "TFile.h"


#include <fstream> 
#include <iostream> 

void fast_merge(char *input_file_list = "inc_bg_229.txt", char*outFile ="test.root", int tagflag = 0)
//void fast_merge(char *input_file_list = "temp_list.txt", char*outFile ="test.root")
{

  const int numhist = 4;

  string histoname[numhist] = {
    "C0_ZVERTEX",
    "C0_MODMAP0"
    "C0_TRIGPT",

    "bgC0_TRIGPT",
    "Cent_0__MWTRIGPT",
    "Cent_0_MWPHIPT0",     

    "Cent_0_MWPHIPT1",
    "Cent_0_MWPHIPT2",
    "Cent_0_MWPHIPT3",

    "Cent_0_MWPHIFOLDPT0",     
    "Cent_0_MWPHIFOLDPT1",
    "Cent_0_MWPHIFOLDPT2",

    "Cent_0_MWPHIFOLDPT3",
    " ",
    " " };
  
    
  string histoname1 = "Cent_0__PT1PT2DPHI";

  
  TH3F *hPT1PT2DPHI_BG_INC = new TH3F(); 
  TH1 * hclone[numhist] = {0};
  
  
  char pwgfile[1000];
  ifstream runlist1(input_file_list);
  int counter = 0;

  TH3F *temp_hPT1PT2DPHI_BG_INC;
  TH1 * temp_h;
  TFile *inFile;
  TFile *inFile0;

  int maxrun = numhist;
  if (tagflag == 0) maxrun = 3;

  while (runlist1.getline(pwgfile,500))
    {
      //cout << pwgfile << " " << counter << endl;
      if (counter % 25 == 0) cout << pwgfile << " " << counter << endl;
      if(counter==0) inFile0 = new TFile(pwgfile);
      else inFile = new TFile(pwgfile);
      
      for (int i=0; i < maxrun; i++)
	{
	  //TH3F *temp_hPT1PT2DPHI_BG_INC = (TH3F*) inFile->Get(histoname1.c_str())->Clone();
	  if(counter==0) temp_h = (TH1*) inFile0->Get(histoname[i].c_str());
	  else temp_h = (TH1*) inFile->Get(histoname[i].c_str());

	  if (!temp_h) continue;
 
	  if (counter == 0) // fg file
	    {
	      cout << "merging " << histoname[i] << endl;
	      hclone[i] = (TH1*) temp_h->Clone();

	    }
	  else
	    {
	      temp_h->SetName("temp");
	      //      cout << "adding here" << endl;
	      hclone[i]->Add(temp_h);
	      //cout << "here" << endl;
	      //	      delete temp_h;
	      //delete inFile;
	      //temp_h->Delete();
	      //inFile->Delete();
	      delete temp_h;
	    }
	}

      delete inFile;

      counter++;
      if(counter %100 ==0) cout<<counter<<" of files have been merged."<<endl;
    }
  
  cout<<"total "<<counter<<" files have been merged"<<endl;
  
  TFile *hfile = TFile::Open(outFile, "RECREATE", "fast_merge");
  
  //  hPT1PT2DPHI_BG_INC->Write();
  
  for (int i = 0; i < numhist; i++)
    {
      if (hclone[i]) hclone[i]->Write();
    }
  cout << "done writing " << endl;

  hfile->Close();
  
  cout<< "the output file is "<<outFile<<endl;
}
