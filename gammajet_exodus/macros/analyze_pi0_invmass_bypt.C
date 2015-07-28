void analyze_pi0_invmass_bypt(){
  gROOT->SetStyle("Plain");

  //TFile *fin=new TFile("/phenix/hp/data06/jfrantz/merge_qm08/bf/fullsall1.root");
  //TFile *fin=new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/lowpt_all.root");
  //TFile *fin=new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent1/pi.root");

  //TFile *fin=new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/new_calibs/some.root");

  //TFile *fin=new TFile("all973.root");
  //TFile *fin=new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent0/cent0_posres_all.root");

  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4_b2w.root");

  //TFile *fin =  new TFile("/direct/phenix+u/workarea/manguyen/AllCode/awayside/wrk/megan_MC/cent0_newrange_merge.root");

  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_cent0_etas.root");
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_nopos.root");
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/all_posres.root");
  //TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent0/direct/cent0_newdir_a2p.root");
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/new/all_posresmymap.root");
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/new/set1/all_set1.root");
  //data
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4all.root");
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/nov172008/mergedhistos/r4alla.root");

  //TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/fe_sharks/good_ak.root");
  //TFile *fin =  new TFile("/phenix/scratch/mjuszkie/gammajet.exodus/fe_sharks/cent0_bigIII.root");

  //TFile *fin =  new TFile("output_merged/fe/new/all_fenew.root");
  // TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/exodus/cent0/meandmatt1.root");
  //TFile *fin =  new TFile("/phenix/u/workarea/mjuszkie/run7_QM09/taxi138/q2/first_rg6.root");
  //TFile *fin =  new TFile("/phenix/u/workarea/mjuszkie/run7_QM09/taxi136/rg6.root");

  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/QM09/taxi141/020/all_rg6.root");

  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/QM09/taxi141/2040/pwg_rg6.root");
  //  TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/QM09/taxi141/4060/pwg_rg6.root");
  //TFile *fin =  new TFile("/phenix/hp/data10/mjuszkie/run7/postQM09/taxi174/merged/oct292009/C2_rg6.root");


  //  TFile *fin =  new TFile("run10_ppg113/cs.root");
  TFile *fin =  new TFile("test/mergelarge1.root");
 
  


  TH2D *INVMASS[8]; TH1D *INVMASS_pro[4][8];
  
  TCanvas *c[8];

  char name[100];
  

  TGraph *peak=new TGraph(8);
  TGraph *width=new TGraph(8);


  for(int icent=0;icent<8;icent++){
    sprintf(name,"c%d",icent);
    c[icent]=new TCanvas(name,name,1);
    c[icent]->Divide(2,2);
  }
  

    

 

      //for(int isect=0;isect<8;isect++){
      //sprintf(name,"Cent_%d__INVMASS_%dorg",isect,icent);
      //sprintf(name,"piC%d_INVMASS_%d",isect,icent);
      sprintf(name,"INVMASS_%d",icent);
      //sprintf(name,"itagC%d_INVMASS_%d",isect,icent);
	//sprintf(name,"INVMASS_PBGL_%d",icent);
      //if(icent<3)sprintf(name,"pi_INVMASS_%d",icent);
      //else sprintf(name,"INVMASS_%d",icent-3);
    for(int icent=0;icent<8;icent++){


      int isect=icent;

	cout<<name<<endl;



	//INVMASS[icent]=(TH2D*)fin->Get(name);


	if(icent<7) INVMASS[icent]=(TH2D*)fin->Get(name);
	cout << "got it" <<endl;
	
	if(icent==0)INVMASS[7]=(TH2D*)fin->Get(name);
	else {
	  TH2F *htemp=(TH2F*)fin->Get(name);
	  INVMASS[7]->Add(htemp);
	}
	
	//}
      
      float par0=0;
      float par1=0;
      float par2=0;
      float par3=0;
      float par4=0;
    
      for(int ipt=0;ipt<4;ipt++){
	cout<<" ipt "<<ipt<<endl;
	c[icent]->cd(ipt+1);      
	char name[100];
	sprintf(name,"invmass_%d_%d",ipt,icent);
	
	if(ipt==0)INVMASS_pro[ipt][icent]=(TH1D*)INVMASS[icent]->ProjectionY(name,51,70);
	if(ipt==1)INVMASS_pro[ipt][icent]=(TH1D*)INVMASS[icent]->ProjectionY(name,71,90);
	if(ipt==2)INVMASS_pro[ipt][icent]=(TH1D*)INVMASS[icent]->ProjectionY(name,91,120);
	if(ipt==3)INVMASS_pro[ipt][icent]=(TH1D*)INVMASS[icent]->ProjectionY(name,121,150);

	INVMASS_pro[ipt][icent]->Sumw2();

	if(icent==7){
	  for(icomb=0; icomb<7; icomb++){
	    INVMASS_pro[ipt][icent]->Add(INVMASS_pro[ipt][icomb]);
	  }
	}

      INVMASS_pro[ipt][icent]->GetXaxis()->SetRangeUser(0.,0.3);
      //if(icent==3){
      
      INVMASS_pro[ipt][icent]->Fit("gaus","","",0.09,0.21);
      
      TF1 *fgaustemp=INVMASS_pro[ipt][icent]->GetFunction("gaus");

      
      //float max =      INVMASS_pro[ipt][icent]->GetMaximum();

	TF1 *gauspol = NULL;
	if(isect==0&&ipt==4){
	  gauspol=new TF1("gauspol","gaus(0)+pol0(3)",0.07,0.21);
	  gauspol->SetParameters(100.,0.14,0.13,1.0);
	}
	else{
	  gauspol=new TF1("gauspol","gaus(0)+pol1(3)",0.07,0.21);
	  //gauspol->SetParameters(100.,0.14,0.13,1.0,0.0);
	  gauspol->SetParameters(fgaustemp->GetParameter(0),fgaustemp->GetParameter(1),fgaustemp->GetParameter(2),1.0,0.0);
	}
	
	par0=gauspol->GetParameter(0);
	par1=gauspol->GetParameter(1);
	par2=gauspol->GetParameter(2);
	par3=gauspol->GetParameter(3);
	par4=gauspol->GetParameter(4);

	       

	INVMASS_pro[ipt][icent]->Fit("gauspol","R");


	if(isect==3){
	  cout<<" peak position "<<gauspol->GetParameter(2)<<" width "<<gauspol->GetParameter(3)<<endl;
	  peak[isect]->SetPoint(ipt,gauspol->GetParameter(1));
	  width[isect]->SetPoint(ipt,gauspol->GetParameter(2));

	}

	TF1 *gauspart = new TF1("gauspart","gaus",0.09,0.22);
	TF1 *pol1part = NULL;
	if(isect==0&&ipt==4){
	  pol1part= new TF1("pol1part","pol0",0.07,0.22);
	}
	else{
	  pol1part= new TF1("pol1part","pol1",0.07,0.22);
	}
	gauspart->SetParameters(gauspol->GetParameter(0),gauspol->GetParameter(1),gauspol->GetParameter(2));

	if(isect==0&&ipt==4){
	pol1part->SetParameter(0,gauspol->GetParameter(3));

	}
	else{
	pol1part->SetParameters(gauspol->GetParameter(3),gauspol->GetParameter(4));
	}


	
	if(ipt==0)sprintf(name,"5-7 GeV",ipt);
	if(ipt==1)sprintf(name,"7-9 GeV",ipt);
	if(ipt==2)sprintf(name,"9-12 GeV",ipt);
	if(ipt==3)sprintf(name,"12-15 GeV",ipt);
	INVMASS_pro[ipt][icent]->SetTitle(name);

	
	//INVMASS_pro[ipt][icent]->SetLabelSize(0.08,"X");
	//INVMASS_pro[ipt][icent]->SetLabelOffset(0.0,"X");
	//INVMASS_pro[ipt][icent]->SetLabelSize(0.08,"Y");
	//INVMASS_pro[ipt][icent]->SetTitleSize(0.08,"X");
	//INVMASS_pro[ipt][icent]->SetTitleSize(0.08,"Y");
	INVMASS_pro[ipt][icent]->SetXTitle("m_{#gamma#gamma} [GeV]");
	INVMASS_pro[ipt][icent]->SetXTitle("m_{#gamma#gamma} [GeV]");
	INVMASS_pro[ipt][icent]->Draw();

	cout <<"drawing " << icent <<ipt <<endl;

	pol1part->SetLineColor(2);
	pol1part->Draw("same");
	
	cout << "sec " << icent << "ipt " << ipt <<" peak position "
	     <<gauspol->GetParameter(1)<<" width "<<gauspol->GetParameter(2)<<endl;
	
	float window1=gauspol->GetParameter(1)-2*(gauspol->GetParameter(2));
	float window2=gauspol->GetParameter(1)+2*(gauspol->GetParameter(2));

	cout << "peak_window " << window1 << " to " << window2 <<endl;

	window1=0.12;
	window2=0.16;

	float s2b=0.0;


	
	if(pol1part->Integral(window1,window2)>0) s2b=gauspart->Integral(window1,window2)/pol1part->Integral(window1,window2);


	//float s2b=gauspart->Integral(0.12,0.16)/pol1part->Integral(0.12,0.16);

	TLatex t;
	t.SetNDC();
	char s2bname[100];
	sprintf(s2bname,"S/B = %2.1f",s2b);
	t.SetTextSize(0.12);
	t.DrawLatex(0.5,0.7,s2bname);

      }
    
    }

    /*
  TCanvas *c2=new TCanvas("c2","c2",1);
  c2->Divide(1,2);
  c2->cd(1);
  peak->SetMarkerStyle(8);
  peak->SetTitle("Peak Position");
  peak->GetXaxis()->SetTitle("sector");
  peak->GetYaxis()->SetTitle("Peak pos [cm]");
  peak->Draw("AP");
  c2->cd(2);
  width->SetTitle("Width");
  width->GetXaxis()->SetTitle("sector");
  width->GetYaxis()->SetTitle("width [GeV]");
  width->SetMarkerStyle(8);
  width->Draw("AP");

    */

}
