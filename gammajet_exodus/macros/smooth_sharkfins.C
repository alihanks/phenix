void smooth_sharkfins(char * reg, char * flat, char * out)
{
  TH1D *hshark_large[5][33];
  TH1D *hshark_small[7][33];
  TH1D *hshark_alt[7][33];

  TFile *fin=new TFile(reg);
  TFile *fflat=new TFile(flat);

  TFile *fout=new TFile(out,"UPDATE");
  
    for(int izemc=0;izemc<33;izemc++){ //33
    cout<<" izemc "<<izemc<<endl;
      for(int idecl=0;idecl<5;idecl++){
      char sharkname[100]; 
      sprintf(sharkname,"hshark_large_%d_%d",idecl,izemc);
      TH1D * hcop = fout->Get(sharkname);
      TH1D * horig = fin->Get(sharkname);
      //      TH1D * hnflat = fflat->Get(sharkname);
      
      horig->Rebin(5);
      horig->Scale(1.0/5.0);
      
      horig->Smooth(5);

      //      float nflt = hnflat->Integral(160,200);
      // float flt = hcop->Integral(160,200);
      
      //      hnflat->Scale(flt/nflt);


      for (int i =0; i<400; i++)
	{
	  if ( (idecl == 0 && i > 100) ||
	       (idecl == 1 && i > 130 ) ||
	       (idecl == 2 && i > 170 ) ||
	       (idecl == 4 && i > 230 ) )
	    {
	      //	      hcop->SetBinContent(i+1,horig->GetBinContent(i+1));	       	      
	      hcop->SetBinContent(i+1,horig->Interpolate(hcop->GetBinCenter(i+1)));
	    }
	  
	}
      
      }
    }
    
    fout->Write();
    
}
