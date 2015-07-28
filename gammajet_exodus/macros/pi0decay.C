void pi0decay(char* inputfile)
{

  //Reset ROOT and connect tree file
  gROOT->Reset();
  gStyle->SetPalette(1);
  TFile *f = new TFile(inputfile,"READ");
  
  TH1D *gamma[9];
  TH1D *prob[9];
  TH2D *hitdiag[6];

  char dummy[128];

  for(int i=0;i<9;i++){
    sprintf(dummy,"gammaprob%d",i);
    gamma[i]=(TH1D*)(gDirectory->Get(dummy))->Clone();

    sprintf(dummy,"prob%d",i);
    prob[i]=(TH1D*)(gDirectory->Get(dummy))->Clone();

    prob[i]->SetMarkerStyle(20+i);
    prob[i]->SetMarkerSize(1);
    prob[i]->SetMarkerColor(kBlack+i);
    prob[i]->Rebin(10);
    gamma[i]->Rebin(10);
  }  

  leg = new TLegend(0.7,0.2,0.9,0.6);
  
  for(int i=0;i<6;i++)
    {
      sprintf(dummy,"hitdiag_sec%d",i);
      hitdiag[i]=(TH2D*)gDirectory->Get(dummy);
    }
  
  for(i=0; i<6; i++)
    {
      sprintf(dummy,"sec%d",i);
      leg->AddEntry(prob[i],dummy,"p");
    }
  
  sprintf(dummy,"PbSc");
  leg->AddEntry(prob[6],dummy,"p");
  //  sprintf(dummy,"West");
  //  leg->AddEntry(prob[7],dummy,"p");
  //  sprintf(dummy,"East");
  //  leg->AddEntry(prob[8],dummy,"p");
 
  TCanvas *c1 = new TCanvas("c1","factor study");
   
  for(i=0;i<6;i++)
    {
      prob[i]->Divide(gamma[i]);
      if(i==0)
	prob[i]->Draw("p");
      else
	prob[i]->Draw("psame");
    }
  
  leg->Draw();
  
  gStyle->SetOptStat(0);
  TCanvas *c2 = new TCanvas("c2","hitdiag");
  c2->Divide(2,3);
  for(i=0;i<6;i++)
    {
      c2->cd(i+1);
      hitdiag[i]->Draw();
    }

  c1->cd();  
}
