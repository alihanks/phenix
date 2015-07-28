void conv_corrs(){

  TFile *fin=new TFile("shark57.root");
  TH1D *py=(TH1D *)fin->Get("py");

  TCanvas *c1=new TCanvas("c1","c1",1);
  c1->Divide(1,3);
  
  c1->cd(1);
  py->SetTitle("5-7 shark fin");
  py->Draw();
  TF1 *fpower=new TF1("fpower","10000./pow(x,7.1)",3,20);

  c1->cd(2);
  fpower->Draw();

  TH1D *conv=py->Clone();
  conv->Reset();
  conv->SetName("conv");

  for(int i=0;i<conv->GetNbinsX();i++){

    float xval  = conv->GetBinCenter(i+1);
    
    float bc = py->GetBinContent(i+1);

    float product = fpower->Eval(xval)*bc;

    conv->SetBinContent(i+1,product);




  }
  c1->cd(3);
  conv->SetTitle("shark fin * pi0 spectrum");
  conv->Draw();


}
