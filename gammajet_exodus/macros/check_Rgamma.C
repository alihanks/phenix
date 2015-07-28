void check_Rgamma(){

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>

  TH1F *TRIGPT = new TH1F("TRIGPT","TRIGPT",200,0,20);



  char infilelist[1000];
  //sprintf(infilelist,"/phenix/hp/data10/manguyen/run6/compressed_files_stripe/x%d",ifile);

  sprintf(infilelist,"/phenix/hp/data10/mjuszkie/exodus/piaa");

  char goodfname[500];
  ifstream goodf(infilelist);
  int nfiles=0;
  while (goodf.getline(goodfname, 500)){


    TFile *fin=new TFile(goodfname);
 

  //TFile *oscarin = new TFile(oscarfilename.c_str());
  //TFile *oscarin = new TFile("/phenix/scratch/mjuszkie/run4/fb_hpdst/ntuples/ntuple_0000122223_0001.root.0_20.trigpi0.root");
  TTree *ch1 = (TTree*)fin->Get("triggaz");
  int n1ent = ch1->GetEntries();

  float runno;   
  float seqno;   
  float evtno;   
  float E;	 
  float Esub;    
  float eta;	 
  float the;  	 
  float phi; 	 
  float vertex;  
  float ngoodpi0s;
  float px;	 
  float py;   	 
  float pz;  	 
  float x;  	 
  float y;   	 
  float z;   	 
  float sector;  
  float pc3dr;	 
  float emctof;	 
  float percent; 
  float thetaDEG;
  float towerid1;
  float towerid2;

 
  ch1->SetBranchAddress("runno",		&runno);   
  ch1->SetBranchAddress("seqno",		&seqno);   
  ch1->SetBranchAddress("evtno",		&evtno);   
  ch1->SetBranchAddress("E",		&E);	 
  ch1->SetBranchAddress("Esub",		&Esub);    
  ch1->SetBranchAddress("eta",		&eta);	 
  ch1->SetBranchAddress("the",		&the);  	 
  ch1->SetBranchAddress("phi",		&phi); 	 
  ch1->SetBranchAddress("vertex",	&vertex);  
  ch1->SetBranchAddress("ngoodpi0s",	&ngoodpi0s);
  ch1->SetBranchAddress("px",		&px);	 
  ch1->SetBranchAddress("py",		&py);   	 
  ch1->SetBranchAddress("pz",		&pz);  	 
  ch1->SetBranchAddress("x",		&x);  	 
  ch1->SetBranchAddress("y",		&y);   	 
  ch1->SetBranchAddress("z",		&z);   	 
  ch1->SetBranchAddress("sector",		&sector);  
  ch1->SetBranchAddress("pc3dr",		&pc3dr);	 
  ch1->SetBranchAddress("emctof",		&emctof);	 
  ch1->SetBranchAddress("percent",	&percent); 
  ch1->SetBranchAddress("thetaDEG",	&thetaDEG);
  ch1->SetBranchAddress("towerid1",	&towerid1);
  ch1->SetBranchAddress("towerid2",	&towerid2);
 


  for(int i=0; i<n1ent; i++)
    {
      //if(as_eta[i]<0.35 && eta[i]>-0.35)
      //{
   //float pt=sqrt(as_px[i]*as_px[i]+as_py[i]*as_py[i]);
   //if(i%100==0) cout << "pt= " << as_pt << " py= " << as_py[i] << " px= " << as_px[i] <<endl;

   ch1->GetEntry(i);
   if(percent<20 && percent>0)
     {               
       //float as_pt=sqrt(px*px+py*py);

       TRIGPT->Fill(E-Esub);

       /*
       CF->Fill(as_pt);
       float wt = dNdpt->Eval(as_pt,0,0);
       cout <<"as_pt" << as_pt << " wt " << wt<<endl;
       wTRIGPT->Fill(as_pt, wt);
       //wTRIGPT->Fill(wt);
       CFw->Fill(as_pt, wt);
       */

     }
   //}
    }
  }
  TRIGPT->Draw();
  TFile *fout=new TFile("piaa.root","RECREATE");
  TRIGPT->Write();
}
