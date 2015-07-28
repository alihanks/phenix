void make_weight_funcs(){

  TFile *fin=new TFile("100m_pi0_smeared_seff_tag_momw.root");

  // TH2F *pivsgampi=(TH2F*)fin->Get("ptpivsptgam_tag");
 TH2F *pivsgampi=(TH2F*)fin->Get("ptpivsptgam");


  TH1D *wfpi57=(TH1D*)pivsgampi->ProjectionY("wfpi57",101,140);
  TH1D *wfpi79=(TH1D*)pivsgampi->ProjectionY("wfpi79",141,180);
  TH1D *wfpi912=(TH1D*)pivsgampi->ProjectionY("wfpi912",181,240);
  TH1D *wfpi510=(TH1D*)pivsgampi->ProjectionY("wfpi510",101,200);

  TH1D *wfpi56=(TH1D*)pivsgampi->ProjectionY("wfpi56",101,120);
  TH1D *wfpi67=(TH1D*)pivsgampi->ProjectionY("wfpi67",121,140);
  TH1D *wfpi78=(TH1D*)pivsgampi->ProjectionY("wfpi78",141,160);
  TH1D *wfpi89=(TH1D*)pivsgampi->ProjectionY("wfpi89",161,180);
  TH1D *wfpi910=(TH1D*)pivsgampi->ProjectionY("wfpi910",181,200);
  TH1D *wfpi1012=(TH1D*)pivsgampi->ProjectionY("wfpi1012",201,240);
  TH1D *wfpi1220=(TH1D*)pivsgampi->ProjectionY("wfpi1220",241,400);

  TCanvas *c1=new TCanvas("c1","c1",1);
  c1->Divide(2,2);
  c1->cd(1);
  wfpi57->Draw();
  c1->cd(2);
  wfpi79->Draw();
  c1->cd(3);
  wfpi912->Draw();
  c1->cd(4);
  wfpi510->Draw();




  TCanvas *cg=new TCanvas("cg","cg",1);
  int write=1;
  
  if(write){
    TGraph *wf57=new TGraph();
    TGraph *wf79=new TGraph();
    TGraph *wf912=new TGraph();
    TGraph *wf510=new TGraph();


    TGraph *wf56=new TGraph();
    TGraph *wf67=new TGraph();
    TGraph *wf78=new TGraph();
    TGraph *wf89=new TGraph();
    TGraph *wf910=new TGraph();
    TGraph *wf1012=new TGraph();
    TGraph *wf1220=new TGraph();

    
    wf57->SetName("wf57");
    wf79->SetName("wf79");
    wf912->SetName("wf912");
    wf510->SetName("wf510");

    wf56->SetName("wf56");
    wf67->SetName("wf67");
    wf78->SetName("wf78");
    wf89->SetName("wf89");
    wf910->SetName("wf910");
    wf1012->SetName("wf1012");
    wf1220->SetName("wf1220");


   
  for(int i=0;i<wfpi57->GetNbinsX();i++){
    float pt=wfpi57->GetBinCenter(i+1);  

    float bc57=wfpi57->GetBinContent(i+1);  
    float bc79=wfpi79->GetBinContent(i+1);  
    float bc912=wfpi912->GetBinContent(i+1);  
    float bc510=wfpi510->GetBinContent(i+1);  

    float bc56=wfpi56->GetBinContent(i+1);  
    float bc67=wfpi67->GetBinContent(i+1);  
    float bc78=wfpi78->GetBinContent(i+1);  
    float bc89=wfpi89->GetBinContent(i+1);  
    float bc910=wfpi910->GetBinContent(i+1);  
    float bc1012=wfpi1012->GetBinContent(i+1);  
    float bc1220=wfpi1220->GetBinContent(i+1);  

    wf57->SetPoint(i,pt,bc57);
    wf79->SetPoint(i,pt,bc79);
    wf912->SetPoint(i,pt,bc912);
    wf510->SetPoint(i,pt,bc510);

    wf56->SetPoint(i,pt,bc56);
    wf67->SetPoint(i,pt,bc67);
    wf78->SetPoint(i,pt,bc78);
    wf89->SetPoint(i,pt,bc89);
    wf910->SetPoint(i,pt,bc910);
    wf1012->SetPoint(i,pt,bc1012);
    wf1220->SetPoint(i,pt,bc1220);


   }

  cg->Divide(2,2);
  cg->cd(1);
  wf57->Draw("AP");
   cg->cd(2);
  wf79->Draw("AP");
  cg->cd(3);
  wf912->Draw("AP");
  cg->cd(4);
  wf510->Draw("AP");
  
  TFile *fout=new TFile("smearedwfuncsnorm3.root","RECREATE");
  wf57->Write();
  wf79->Write();
  wf912->Write();
  wf510->Write();

  wf56->Write();
  wf67->Write();
  wf78->Write();
  wf89->Write();
  wf910->Write();
  wf1012->Write();
  wf1220->Write();

  
  fout->Close();
  
}
}
