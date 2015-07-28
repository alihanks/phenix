void make_v2_graphs(const char* output = "v2_inputs.root")
{
  TFile* fout = new TFile(output,"recreate");
  fout->cd();
  TGraphErrors* trig_inc_v2[4];  TGraphErrors* trig_inc_sys[4];
  TGraphErrors* trig_dec_v2[4];  TGraphErrors* trig_dec_sys[4];
  TGraphErrors* trig_pi0_v2[4];  TGraphErrors* trig_pi0_sys[4];
  TGraphErrors* hassoc_v2[4];  TGraphErrors* hassoc_sys[4];
  
  //[centrality][pt]
  double inc_v2[4][4]; double inc_v2_err[4][4]; double inc_v2_sys[4][4];
  double dec_v2[4][4]; double dec_v2_err[4][4]; double dec_v2_sys[4][4];
  double pi0_v2[4][4]; double pi0_v2_err[4][4]; double pi0_v2_sys[4][4];
  double hadron_v2[4][5]; double hadron_v2_err[4][5]; double hadron_v2_sys[4][5];

  inc_v2[0][0] = 0.0423226; inc_v2_err[0][0] = 0.00175307; inc_v2_sys[0][0] = 0.0076029;
  inc_v2[1][0] = 0.0852710; inc_v2_err[1][0] = 0.00207489; inc_v2_sys[1][0] = 0.0179611;
  inc_v2[2][0] = 0.1379880; inc_v2_err[2][0] = 0.00432101; inc_v2_sys[2][0] = 0.0183026;
  inc_v2[3][0] = 0.4418120; inc_v2_err[3][0] = 0.01932520; inc_v2_sys[3][0] = 0.0197302;
  pi0_v2[0][0] = 0.0835297; pi0_v2_err[0][0] = 0.00287215; pi0_v2_sys[0][0] = 0.0118235;
  pi0_v2[1][0] = 0.1335370; pi0_v2_err[1][0] = 0.00330527; pi0_v2_sys[1][0] = 0.0229443;
  pi0_v2[2][0] = 0.1744460; pi0_v2_err[2][0] = 0.00671865; pi0_v2_sys[2][0] = 0.0244753;
  pi0_v2[3][0] = 0.4214010; pi0_v2_err[3][0] = 0.02954310; pi0_v2_sys[3][0] = 0.0485074;
  dec_v2[0][0] = 0.0802729; dec_v2_err[0][0] = 0.00259655; dec_v2_sys[0][0] = 0.0106477;
  dec_v2[1][0] = 0.1270690; dec_v2_err[1][0] = 0.00292113; dec_v2_sys[1][0] = 0.0236515;
  dec_v2[2][0] = 0.1708000; dec_v2_err[2][0] = 0.00588472; dec_v2_sys[2][0] = 0.0221000;
  dec_v2[3][0] = 0.4716660; dec_v2_err[3][0] = 0.02824160; dec_v2_sys[3][0] = 0.0515539;

  inc_v2[0][1] = 0.0461296; inc_v2_err[0][1] = 0.00454517; inc_v2_sys[0][1] = 0.00432369;
  inc_v2[1][1] = 0.0738878; inc_v2_err[1][1] = 0.00558803; inc_v2_sys[1][1] = 0.01509860;
  inc_v2[2][1] = 0.1432810; inc_v2_err[2][1] = 0.01209420; inc_v2_sys[2][1] = 0.01267580;
  inc_v2[3][1] = 0.5100750; inc_v2_err[3][1] = 0.05529200; inc_v2_sys[3][1] = 0.04872020;
  pi0_v2[0][1] = 0.0755782; pi0_v2_err[0][1] = 0.00531345; pi0_v2_sys[0][1] = 0.00960588;
  pi0_v2[1][1] = 0.1162020; pi0_v2_err[1][1] = 0.00587208; pi0_v2_sys[1][1] = 0.02417620;
  pi0_v2[2][1] = 0.1548150; pi0_v2_err[2][1] = 0.01164770; pi0_v2_sys[2][1] = 0.01545110;
  pi0_v2[3][1] = 0.5469130; pi0_v2_err[3][1] = 0.05027790; pi0_v2_sys[3][1] = 0.03882330;
  dec_v2[0][1] = 0.0769908; dec_v2_err[0][1] = 0.00440387; dec_v2_sys[0][1] = 0.01112830;
  dec_v2[1][1] = 0.1142120; dec_v2_err[1][1] = 0.00505774; dec_v2_sys[1][1] = 0.02610160;
  dec_v2[2][1] = 0.1723320; dec_v2_err[2][1] = 0.01048860; dec_v2_sys[2][1] = 0.01933210;
  dec_v2[3][1] = 0.5762520; dec_v2_err[3][1] = 0.05424370; dec_v2_sys[3][1] = 0.04244810;

  inc_v2[0][2] = 0.0398415; inc_v2_err[0][2] = 0.00836215; inc_v2_sys[0][2] = 0.00344175;
  inc_v2[1][2] = 0.0520371; inc_v2_err[1][2] = 0.01062420; inc_v2_sys[1][2] = 0.02247170;
  inc_v2[2][2] = 0.1245420; inc_v2_err[2][2] = 0.02340030; inc_v2_sys[2][2] = 0.01237760;
  inc_v2[3][2] = 0.4952110; inc_v2_err[3][2] = 0.11218300; inc_v2_sys[3][2] = 0.15661100;
  pi0_v2[0][2] = 0.0839048; pi0_v2_err[0][2] = 0.01013930; pi0_v2_sys[0][2] = 0.01805810;
  pi0_v2[1][2] = 0.1030550; pi0_v2_err[1][2] = 0.01122400; pi0_v2_sys[1][2] = 0.03859360;
  pi0_v2[2][2] = 0.0801250; pi0_v2_err[2][2] = 0.05710570; pi0_v2_sys[2][2] = 0.02456170;
  pi0_v2[3][2] = 0.5638220; pi0_v2_err[3][2] = 0.10233000; pi0_v2_sys[3][2] = 0.02597060;
  //testing old v2 because 9-12 seems high
  //v2[2][2] = 0.242351; err_v2_pi0[2][2] = 0.0226531; sys_err_v2_pi0[2][2] = 0.0129369;
  dec_v2[0][2] = 0.0831338; dec_v2_err[0][2] = 0.00951379; dec_v2_sys[0][2] = 0.02080890;
  dec_v2[1][2] = 0.1005180; dec_v2_err[1][2] = 0.01043080; dec_v2_sys[1][2] = 0.03353270;
  dec_v2[2][2] = 0.2019460; dec_v2_err[2][2] = 0.02221120; dec_v2_sys[2][2] = 0.01217420;
  dec_v2[3][2] = 0.5698930; dec_v2_err[3][2] = 0.11125100; dec_v2_sys[3][2] = 0.03953840;

  inc_v2[0][3] = 0.0169588; inc_v2_err[0][3] = 0.0189458; inc_v2_sys[0][3] = 0.00465717;
  inc_v2[1][3] = 0.0644017; inc_v2_err[1][3] = 0.0252076; inc_v2_sys[1][3] = 0.02599590;
  inc_v2[2][3] = 0.0260897; inc_v2_err[2][3] = 0.0570124; inc_v2_sys[2][3] = 0.05901420;
  inc_v2[3][3] = -0.299766; inc_v2_err[3][3] = 0.2958920; inc_v2_sys[3][3] = 0.79791400;
  pi0_v2[0][3] = 0.1055550; pi0_v2_err[0][3] = 0.0258798; pi0_v2_sys[0][3] = 0.04751440;
  pi0_v2[1][3] = 0.0861201; pi0_v2_err[1][3] = 0.0290550; pi0_v2_sys[1][3] = 0.02462110;
  pi0_v2[2][3] = 0.1169110; pi0_v2_err[2][3] = 0.0579355; pi0_v2_sys[2][3] = 0.06628670;
  pi0_v2[3][3] = 0.4231040; pi0_v2_err[3][3] = 0.3135600; pi0_v2_sys[3][3] = 0.25851400;
  dec_v2[0][3] = 0.1026670; dec_v2_err[0][3] = 0.0238731; dec_v2_sys[0][3] = 0.03616490;
  dec_v2[1][3] = 0.0893677; dec_v2_err[1][3] = 0.0271222; dec_v2_sys[1][3] = 0.03356940;
  dec_v2[2][3] = 0.1307690; dec_v2_err[2][3] = 0.0609939; dec_v2_sys[2][3] = 0.07152490;
  dec_v2[3][3] = 0.3373860; dec_v2_err[3][3] = 0.2751690; dec_v2_sys[3][3] = 0.37330200;

  hadron_v2[0][0]=0.0440219; hadron_v2_err[0][0]=0.000145508; hadron_v2_sys[0][0]=0.00137609;
  hadron_v2[1][0]=0.0826719; hadron_v2_err[1][0]=0.000171826; hadron_v2_sys[1][0]=0.00296920;
  hadron_v2[2][0]=0.0980676; hadron_v2_err[2][0]=0.000408317; hadron_v2_sys[2][0]=0.00221491;
  hadron_v2[3][0]=0.0868451; hadron_v2_err[3][0]=0.001215900; hadron_v2_sys[3][0]=0.00435710;

  hadron_v2[0][1]=0.0745276; hadron_v2_err[0][1]=7.62922e-05; hadron_v2_sys[0][1]=0.00202200;
  hadron_v2[1][1]=0.1381520; hadron_v2_err[1][1]=8.20737e-05; hadron_v2_sys[1][1]=0.00419830;
  hadron_v2[2][1]=0.1646640; hadron_v2_err[2][1]=0.000219383; hadron_v2_sys[2][1]=0.00260620;
  hadron_v2[3][1]=0.1693740; hadron_v2_err[3][1]=0.001336240; hadron_v2_sys[3][1]=0.00693613;

  hadron_v2[0][2]=0.1095070; hadron_v2_err[0][2]=0.000216781; hadron_v2_sys[0][2]=0.00307190;
  hadron_v2[1][2]=0.1933760; hadron_v2_err[1][2]=0.000249992; hadron_v2_sys[1][2]=0.00579030;
  hadron_v2[2][2]=0.2177190; hadron_v2_err[2][2]=0.000623148; hadron_v2_sys[2][2]=0.00315090;
  hadron_v2[3][2]=0.2150790; hadron_v2_err[3][2]=0.003462490; hadron_v2_sys[3][2]=0.00682347;

  hadron_v2[0][3]=0.1103720; hadron_v2_err[0][3]=0.000591522; hadron_v2_sys[0][3]=0.00322900;
  hadron_v2[1][3]=0.1907190; hadron_v2_err[1][3]=0.000680041; hadron_v2_sys[1][3]=0.00616680;
  hadron_v2[2][3]=0.2077700; hadron_v2_err[2][3]=0.001614320; hadron_v2_sys[2][3]=0.00375637;
  hadron_v2[3][3]=0.2052480; hadron_v2_err[3][3]=0.008518280; hadron_v2_sys[3][3]=0.01208620;

  hadron_v2[0][4]=0.0827500; hadron_v2_err[0][4]=0.00844; hadron_v2_sys[0][4]=0.00500;
  hadron_v2[1][4]=0.1297000; hadron_v2_err[1][4]=0.01350; hadron_v2_sys[1][4]=0.00830;
  hadron_v2[2][4]=0.1718000; hadron_v2_err[2][4]=0.01718; hadron_v2_sys[2][4]=0.00858;
  hadron_v2[3][4]=0.2052480; hadron_v2_err[3][4]=8.52e-4; hadron_v2_sys[3][4]=0.05*hadron_v2[3][4];

  double pt_trig[4] = {5.5,7.5,9.5,12.5};
  double pt_trig_err[4] = {0,0,0,0};
  double pt_assoc[5] = {0.6,1.3,2.3,3.5,5.5};
  double pt_assoc_err[5] = {0,0,0,0,0};
  std::ostringstream name;
  for( int i = 0; i < 4; i++ )
  {
    name.str("");
    name << "gamma_inc_v2_" << i;
    trig_inc_v2[i] = new TGraphErrors(4,pt_trig,inc_v2[i],pt_trig_err,inc_v2_err[i]);
    trig_inc_v2[i]->SetName(name.str().c_str());
    trig_inc_v2[i]->SetTitle(name.str().c_str());
    trig_inc_v2[i]->Write();

    name.str("");
    name << "gamma_inc_v2sys_" << i;
    trig_inc_sys[i] = new TGraphErrors(4,pt_trig,inc_v2[i],pt_trig_err,inc_v2_sys[i]);
    trig_inc_sys[i]->SetName(name.str().c_str());
    trig_inc_sys[i]->SetTitle(name.str().c_str());
    trig_inc_sys[i]->Write();

    name.str("");
    name << "gamma_dec_v2_" << i;
    trig_dec_v2[i] = new TGraphErrors(4,pt_trig,dec_v2[i],pt_trig_err,dec_v2_err[i]);
    trig_dec_v2[i]->SetName(name.str().c_str());
    trig_dec_v2[i]->SetTitle(name.str().c_str());
    trig_dec_v2[i]->Write();
 
    name.str("");
    name << "gamma_dec_v2sys_" << i;
    trig_dec_sys[i] = new TGraphErrors(4,pt_trig,dec_v2[i],pt_trig_err,dec_v2_sys[i]);
    trig_dec_sys[i]->SetName(name.str().c_str());
    trig_dec_sys[i]->SetTitle(name.str().c_str());
    trig_dec_sys[i]->Write();
 
    name.str("");
    name << "pi0_v2_" << i;
    trig_pi0_v2[i] = new TGraphErrors(4,pt_trig,pi0_v2[i],pt_trig_err,pi0_v2_err[i]);
    trig_pi0_v2[i]->SetName(name.str().c_str());
    trig_pi0_v2[i]->SetTitle(name.str().c_str());
    trig_pi0_v2[i]->Write();

    name.str("");
    name << "pi0_v2sys_" << i;
    trig_pi0_sys[i] = new TGraphErrors(4,pt_trig,pi0_v2[i],pt_trig_err,pi0_v2_sys[i]);
    trig_pi0_sys[i]->SetName(name.str().c_str());
    trig_pi0_sys[i]->SetTitle(name.str().c_str());
    trig_pi0_sys[i]->Write();

    name.str("");
    name << "hadron_v2_" << i;
    hassoc_v2[i] = new TGraphErrors(5,pt_assoc,hadron_v2[i],pt_assoc_err,hadron_v2_err[i]);
    hassoc_v2[i]->SetName(name.str().c_str());
    hassoc_v2[i]->SetTitle(name.str().c_str());
    hassoc_v2[i]->Write();

    name.str("");
    name << "hadron_v2sys_" << i;
    hassoc_sys[i] = new TGraphErrors(5,pt_assoc,hadron_v2[i],pt_assoc_err,hadron_v2_sys[i]);
    hassoc_sys[i]->SetName(name.str().c_str());
    hassoc_sys[i]->SetTitle(name.str().c_str());
    hassoc_sys[i]->Write();

  }

  int cent[5] = {0,20,40,60,90};
  int color[4] = {kBlack,kRed,kBlue,kViolet-7};
  int sys_color[4] = {kGray,kRed-9,kBlue-9,kViolet-9};
  TCanvas* can = new TCanvas("can","can");
  can->Divide(2,2,0.001,0.001);
  TH1D* thisto = new TH1D("thisto",";p^{#gamma}_{T}   ;v_{2}   ",100,0.0,15.0);
  thisto->SetAxisRange(0.0,0.75,"Y");
  thisto->SetAxisRange(4.5,14.0,"X");
  TH1D* ahisto = new TH1D("ahisto",";p^{h}_{T}    ;v_{2}   ",100,0.0,15.0);
  ahisto->SetAxisRange(0.0,0.75,"Y");
  ahisto->SetAxisRange(0.0,6.0,"X");
  for( int i = 0; i < 4; i++ )
  {
    TVirtualPad* pad = can->cd(i+1);
    pad->SetRightMargin(0.01);
    pad->SetTopMargin(0.01);
    if( i < 3 ) thisto->Draw();
    else ahisto->Draw();
  }
  TLegend* leg = new TLegend(0.5,0.5,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  for( int ic = 0; ic < 4; ic++ )
  {
    can->cd(1);
    trig_inc_sys[ic]->SetMarkerSize(0);
    trig_inc_sys[ic]->SetLineWidth(10);
    trig_inc_sys[ic]->SetLineColor(sys_color[ic]);
    trig_inc_sys[ic]->Draw("E1,Psame");
    trig_inc_v2[ic]->SetMarkerColor(color[ic]);
    trig_inc_v2[ic]->SetLineColor(color[ic]);
    trig_inc_v2[ic]->Draw("Psame");
    // trig_inc_v2[ic]->Write();
    // trig_inc_sys[ic]->Write();
    can->cd(2);
    trig_dec_sys[ic]->SetMarkerSize(0);
    trig_dec_sys[ic]->SetLineWidth(10);
    trig_dec_sys[ic]->SetLineColor(sys_color[ic]);
    trig_dec_sys[ic]->Draw("E1,Psame");
    trig_dec_v2[ic]->SetMarkerColor(color[ic]);
    trig_dec_v2[ic]->SetLineColor(color[ic]);
    trig_dec_v2[ic]->Draw("Psame");
    can->cd(3);
    trig_pi0_sys[ic]->SetMarkerSize(0);
    trig_pi0_sys[ic]->SetLineWidth(10);
    trig_pi0_sys[ic]->SetLineColor(sys_color[ic]);
    trig_pi0_sys[ic]->Draw("E1,Psame");
    trig_pi0_v2[ic]->SetMarkerColor(color[ic]);
    trig_pi0_v2[ic]->SetLineColor(color[ic]);
    trig_pi0_v2[ic]->Draw("Psame");
    can->cd(4);
    hassoc_sys[ic]->SetMarkerSize(0);
    hassoc_sys[ic]->SetLineWidth(10);
    hassoc_sys[ic]->SetLineColor(sys_color[ic]);
    hassoc_sys[ic]->Draw("E1,Psame");
    hassoc_v2[ic]->SetMarkerColor(color[ic]);
    hassoc_v2[ic]->SetLineColor(color[ic]);
    hassoc_v2[ic]->Draw("Psame");
    name.str("");
    name << cent[ic] << " - " << cent[ic+1] << "%";
    leg->AddEntry(trig_inc_v2[ic],name.str().c_str(),"P");
  }
  leg->Draw();
  can->Write();
}
