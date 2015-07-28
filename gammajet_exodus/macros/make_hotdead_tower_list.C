
void make_hotdead_tower_list(const char* input = "correlations_dAu_total.root", const char* output = "ttemp_livetowers_r8_test.dat", const char* map_name = "nhit_mywarn_as")
{

  TFile *fin=new TFile(input);

  TH3F *MODMAP[8];
  TH2D *MODMAP_pxy[8];

  for(int i=0;i<8;i++){
    
    char name[100];
    sprintf(name,"%s%d",map_name,i);
    cout << "getting map: " << name << endl;
    MODMAP[i]=(TH3F*)fin->Get(name);
    MODMAP[i]->SetAxisRange(5.0,7.0,"X");
    MODMAP_pxy[i]=(TH2D*)MODMAP[i]->Project3D("yz");

    //**********I don't think writing directly to ff1 works 
    //I just so like root -b -q 'make_hotdead_tower_list.C' > taxi230.log &
    FILE *ff1;
 
    ff1 = fopen(output,"w");

    for(int j=0;j<MODMAP_pxy[i]->GetNbinsX();j++){
      for(int k=0;k<MODMAP_pxy[i]->GetNbinsY();k++){
	

	if(i<6){

          if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)> 5 && MODMAP_pxy[i]->GetBinContent(j+1,k+1)< 1e80){ 
            //cout<<i*2592+j+72*k<<" "<<MODMAP_pxy[i]->GetBinContent(j+1,k+1)<<endl;
            fprintf(ff1,"%d %d \n",i*2592+j+72*k,MODMAP_pxy[i]->GetBinContent(j+1,k+1));
          }
	}
	else{
	  if(MODMAP_pxy[i]->GetBinContent(j+1,k+1)>5  && MODMAP_pxy[i]->GetBinContent(j+1,k+1)< 1e80){
	     //cout<<15552+(i-6)*4608+j+96*k<<" "<<MODMAP_pxy[i]->GetBinContent(j+1,k+1)<<endl;
	     fprintf(ff1,"%d %d \n",15552+(i-6)*4608+j+96*k,MODMAP_pxy[i]->GetBinContent(j+1,k+1));
          }

	}
      }
    }

  }

  fclose(ff1); 

}
