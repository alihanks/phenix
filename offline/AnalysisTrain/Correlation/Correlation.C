#include <Correlation.h>
#include <Warnmap.h>
#include <AParticle.h>
#include <ACluster.h>
#include <APiZero.h>
#include <ATrack.h>
#include <AEvent.h>
#include <AMixingPool.h>
#include <AMixingTree.h>

#include <TFile.h>
#include <THmulf.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h>
#include <TTree.h>

#include <Fun4AllServer.h>
#include <Fun4AllReturnCodes.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <PHAngle.h>
#include <PHCentralTrack.h>
#include <PHGlobal.h>
#include <emcClusterContainer.h>
#include <emcClusterContent.h>
// #include <SvxCentralTrackList.h>
// #include <SvxCentralTrack.h>
// #include <SvxClusterList.h>
// #include <SvxCluster.h>
// #include <ConversionVeto.h>
#include <RunHeader.h>
#include <EventHeader.h>
#include <ErtOut.h>

#include <TOAD.h>

#include <string>
#include <sstream>
#include <iostream>

using namespace std;

//some constants for the pbgl
static const float vnx[2][4] = { { 0.92619,     0.999986,     0.921968,     0.702254 }, {-0.92627,     -0.999996,    -0.921187,    -0.702855 } };
static const float vny[2][4] = { {-0.377058,    0.00528539,   0.387266,     0.711927 }, {-0.37686,      0.00268898,   0.389118,     0.711332 } };
static const float vnz[2][4] = { { 4.61055e-05, 7.46704e-05, -0.000252431, -0.000332109}, {-0.000400881, -0.000738336, -0.000926837, -0.000577316} };

Correlation::Correlation(const char* outfile)
{
  NZVTX = 12;
  NCENT = 20;
  NCBINS = 1;
  NMIX = 100;
  Rcut = 0.3;
  data_set = INVALID;
  zVertexCut = 20.0;
  fieldPolarity = 1.0;
  //nsvxpart = 0;
  output = outfile;
  evt = 0;
  event = 0;
  event_z = -9999.;
  cbin = -1;
  photon_pt_min = 5.0; photon_pt_max = 15.0;
  hadron_pt_min = 0.5; hadron_pt_max = 7.0;
  econe_min[0] = 0; econe_min[1] = 0; econe_min[2] = 0; econe_min[3] = 0;
  photon_ecore_min = 1.0;
  pi0_pt_min = 4.0; pi0_pt_max = 17.0;
  //useVtx = 0;
  DiagFlag = 0;
  RecalFlag = 1;
  PC3_NSIGMA = 2.0;
  EMC_NSIGMA = 2.0;
  vetoPtCut = 1.0;
  minAsym = 0.15;
  dofilltime = 1;
  //  fhadeff = NULL;
  fexemb = NULL;
  fhadroneff = NULL;
  //v2file = NULL;
  fpi0eff_0 = NULL;
  fpi0eff_1 = NULL;
  fpi0eff_2 = NULL;
  fpi0eff_3 = NULL;
  grpi0eff_0 = NULL;
  grpi0eff_1 = NULL;
  grpi0eff_2 = NULL;
  grpi0eff_3 = NULL;
  manager = NULL;
  global = NULL;
  emcclustercontainer = NULL;
  particle = NULL;
  // svxcntlist = NULL;
  // svxcluslist = NULL;
  // svxcnttrk = NULL;
  toad_loader = NULL;
  warnmap = NULL;
  atree = NULL;
  
  // for(int i=0; i<4; i++){
  //   gr_inc_v2[i] = NULL;
  //   gr_dec_v2[i] = NULL;
  //   gr_pi0_v2[i] = NULL;
  //   gr_had_v2[i] = NULL;
  //   gr_inc_v2sys[i] = NULL;
  //   gr_dec_v2sys[i] = NULL;
  //   gr_pi0_v2sys[i] = NULL;
  //   gr_had_v2sys[i] = NULL;
  // }
  InitHistos();
}

Correlation::~Correlation()
{
  if(warnmap) delete warnmap;
  
  if(fhadroneff) delete fhadroneff;

  if(grpi0eff_0) delete grpi0eff_0;
  if( fpi0eff_0 ){
    fpi0eff_0->Close();
    delete fpi0eff_0;
  }
  if( grpi0eff_1 ) delete grpi0eff_1;
  if( fpi0eff_1 ){
    fpi0eff_1->Close();
    delete fpi0eff_1;
  }
  if( grpi0eff_2 ) delete grpi0eff_2;
  if( fpi0eff_2 ){
    fpi0eff_2->Close();
    delete fpi0eff_2;
  }
  if( grpi0eff_3 ) delete grpi0eff_3;
  if( fpi0eff_3 ){
    fpi0eff_3->Close();
    delete fpi0eff_3;
  }
  delete toad_loader;
  delete manager;
  if(atree) delete atree;

  for(int ic=0; ic<4; ic++){
    for(int it=0; it<4; it++){
      for(int ip=0; ip<5; ip++){
        delete IncAcc[ic][it][ip];
        delete Pi0Acc[ic][it][ip];
        delete DecAcc[ic][it][ip];
        delete IncAccIso[ic][it][ip];
        delete Pi0AccIso[ic][it][ip];
        delete DecAccIso[ic][it][ip];
      }
      for(int ix=0; ix<6; ix++){
        delete IncAccXi[ic][it][ix];
        delete Pi0AccXi[ic][it][ix];
        delete DecAccXi[ic][it][ix];
        delete IncAccXiIso[ic][it][ix];
        delete Pi0AccXiIso[ic][it][ix];
        delete DecAccXiIso[ic][it][ix];
      }
    }
  }
}

void Correlation::SetWarnmaps(vector<string> filename, vector<double> ptcut)
{
  if( filename.size() != ptcut.size() )
  {
    cout << "ERROR: Correlation::SetWarnmaps() - must have a pt cutoff for each input file!" << endl;
    return;
  }
  for( unsigned int i = 0; i < filename.size(); i++ )
  {
    warnmap_filenames.push_back(filename[i]);
    warnmap_cuts.push_back(ptcut[i]);
  }
}

int Correlation::Init(PHCompositeNode* topNode)
{
  manager = new Fun4AllHistoManager("MyHistos");
  toad_loader = new TOAD("Correlation");
  warnmap = new Warnmap(warnmap_filenames.size());
  warnmap->SetPtRange(warnmap_cuts);
  
  for( unsigned int i = 0; i < warnmap_filenames.size(); i++ )
  {
    cout << "finding WARNMAP input " << i << ": " << warnmap_filenames[i] << "..." << endl;
    string file_location = toad_loader->location(warnmap_filenames[i].c_str());
    warnmap->ReadMap(file_location.c_str(),i);
    cout << "WARNMAP input " << i << ": " << file_location << " loaded" << endl;
  }
  
  string fhadroneff_location = toad_loader->location(_hadeffFilename.c_str());
  SetHadronEfficiency(fhadroneff_location.c_str());
  
  string feff_0_location = toad_loader->location(_pi0effFilename_0.c_str());
  string feff_1_location = toad_loader->location(_pi0effFilename_1.c_str());
  string feff_2_location = toad_loader->location(_pi0effFilename_2.c_str());
  string feff_3_location = toad_loader->location(_pi0effFilename_3.c_str());
  SetTriggerEfficiency(feff_0_location.c_str(), feff_1_location.c_str(), feff_2_location.c_str(), feff_3_location.c_str());
  cout << "pi0 trigger efficiency loaded" << endl;

  string sharkfin_file = toad_loader->location(_sharkfinname.c_str());
  SetSharkFin(sharkfin_file.c_str());
  cout << "sharkfin input loaded" << sharkfin_file<<endl;

  string flow_filename = toad_loader->location(_flowfilename.c_str());
  cout<<"flow_filename: "<<flow_filename.c_str()<<endl;
  SetV2(flow_filename.c_str());
  
  ostringstream bin;
  string name, name_mix;
  
  TH1F* temp;
  TH2F* temp2;
  TH3F* temp3;


  for(int ic = 0; ic < NCBINS; ic++){
    bin.str("");
    bin << ic;

    name = "h1_centrality_c" + bin.str();
    Init1DHisto(temp, name.c_str(), "centrality", 2, 0, 2);
    h1_centrality.push_back(temp);

    name = "h3_dphi_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_dphi.push_back(temp3);
    
    name = "h3_dphi_iso_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_dphi_iso.push_back(temp3);
    
    name = "h3_dphi_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_dphi_mix.push_back(temp3);
    
    name = "h3_dphi_iso_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_dphi_iso_mix.push_back(temp3);
    
    name = "h3_dphi_pi0_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_dphi_pi0.push_back(temp3);
    
    name = "h3_dphi_pi0_iso_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_dphi_pi0_iso.push_back(temp3);
    
    name = "h3_dphi_pi0_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_dphi_pi0_mix.push_back(temp3);

    name = "h3_dphi_pi0_iso_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_dphi_pi0_iso_mix.push_back(temp3);

    name = "h3_ptxidphi_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptxidphi.push_back(temp3);

    name = "h3_ptxidphi_iso_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptxidphi_iso.push_back(temp3);

    name = "h3_ptztdphi_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "zt", 40, 0.0, 2.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptztdphi.push_back(temp3);

    name = "h3_ptxidphi_pi0_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptxidphi_pi0.push_back(temp3);

    name = "h3_ptxidphi_pi0_iso_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptxidphi_pi0_iso.push_back(temp3);

    name = "h3_ptztdphi_pi0_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "zt", 40, 0.0, 2.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptztdphi_pi0.push_back(temp3);

    name = "h3_ptxidphi_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptxidphi_mix.push_back(temp3);

    name = "h3_ptxidphi_iso_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptxidphi_iso_mix.push_back(temp3);

    name = "h3_ptztdphi_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "zt", 40, 0.0, 2.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptztdphi_mix.push_back(temp3);

    name = "h3_ptxidphi_pi0_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptxidphi_pi0_mix.push_back(temp3);

    name = "h3_ptxidphi_pi0_iso_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptxidphi_pi0_iso_mix.push_back(temp3);

    name = "h3_ptztdphi_pi0_mix_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "zt", 40, 0.0, 2.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
    h3_ptztdphi_pi0_mix.push_back(temp3);
    
    name = "h3_dphi_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_dphi_fold.push_back(temp3);

    name = "h3_dphi_iso_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_dphi_iso_fold.push_back(temp3);

    name = "h3_ptxidphi_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptxidphi_fold.push_back(temp3);

    name = "h3_ptxidphi_iso_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptxidphi_iso_fold.push_back(temp3);

    name = "h3_ptztdphi_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "zt", 40, 0.0, 2.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptztdphi_fold.push_back(temp3);
    
    name = "h3_dphi_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_dphi_mix_fold.push_back(temp3);

    name = "h3_dphi_iso_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_dphi_iso_mix_fold.push_back(temp3);

    name = "h3_ptxidphi_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptxidphi_mix_fold.push_back(temp3);

    name = "h3_ptxidphi_iso_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptxidphi_iso_mix_fold.push_back(temp3);

    name = "h3_ptztdphi_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "zt", 40, 0.0, 2.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptztdphi_mix_fold.push_back(temp3);

    name = "h3_dphi_pi0_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_dphi_pi0_fold.push_back(temp3);

    name = "h3_dphi_pi0_iso_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_dphi_pi0_iso_fold.push_back(temp3);

    name = "h3_ptxidphi_pi0_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptxidphi_pi0_fold.push_back(temp3);

    name = "h3_ptxidphi_pi0_iso_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptxidphi_pi0_iso_fold.push_back(temp3);

    name = "h3_ptztdphi_pi0_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "zt", 40, 0.0, 2.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptztdphi_pi0_fold.push_back(temp3);
    
    name = "h3_dphi_pi0_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_dphi_pi0_mix_fold.push_back(temp3);

    name = "h3_dphi_pi0_iso_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_dphi_pi0_iso_mix_fold.push_back(temp3);

    name = "h3_ptxidphi_pi0_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptxidphi_pi0_mix_fold.push_back(temp3);

    name = "h3_ptxidphi_pi0_iso_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "xi", 25, -1.0, 4.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptxidphi_pi0_iso_mix_fold.push_back(temp3);

    name = "h3_ptztdphi_pi0_mix_fold_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 20, 0.0, 20.0, "zt", 40, 0.0, 2.0, "#Delta#phi [rad]", 30, 0., PI);
    h3_ptztdphi_pi0_mix_fold.push_back(temp3);
    
    name = "h1_trig_pt_all_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_all.push_back(temp);
    
    name = "h1_trig_pt_all_iso_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_all_iso.push_back(temp);
    
    name = "h1_trig_pt_inc_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_inc.push_back(temp);
    
    name = "h1_trig_pt_inc_iso_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_inc_iso.push_back(temp);
    
    name = "h1_trig_pt_pi0_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_pi0.push_back(temp);
    
    name = "h1_trig_pt_pi0_iso_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_pi0_iso.push_back(temp);
    
    name = "h1_trig_pt_dec_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} bin",5,-0.5,4.5);
    h1_trig_pt_dec.push_back(temp);
    
    name = "h1_trig_pt_dec_iso_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} bin",5,-0.5,4.5);
    h1_trig_pt_dec_iso.push_back(temp);
    
    name = "h1_trig_pt_dec_niso_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} bin",5,-0.5,4.5);
    h1_trig_pt_dec_niso.push_back(temp);
    
    name = "h1_trig_pt_inc_mix_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_inc_mix.push_back(temp);
    
    name = "h1_trig_pt_inc_iso_mix_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_inc_iso_mix.push_back(temp);
    
    name = "h1_trig_pt_pi0_mix_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_pi0_mix.push_back(temp);
    
    name = "h1_trig_pt_pi0_iso_mix_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",200,0.,20.);
    h1_trig_pt_pi0_iso_mix.push_back(temp);
    
    name = "h1_trig_pt_dec_mix_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} bin",5,-0.5,4.5);
    h1_trig_pt_dec_mix.push_back(temp);
    
    name = "h1_trig_pt_dec_iso_mix_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} bin",5,-0.5,4.5);
    h1_trig_pt_dec_iso_mix.push_back(temp);
    
    name = "h1_part_pt_c" + bin.str();
    Init1DHisto(temp, name.c_str(),"p_{T} [GeV/c]",100,0.,10.);
    h1_part_pt.push_back(temp);
    
    name = "h2_pi0mass_c" + bin.str();
    Init2DHisto(temp2, name.c_str(), "p_{T} [GeV/c]", 100, 0.0, 20.0, "m [GeV/c^{2}]", 200, 0.0, 0.6);
    h2_pi0mass.push_back(temp2);
    
    name = "h2_cluster_dphi_dR_c" + bin.str();
    Init2DHisto(temp2, name.c_str(), "#Delta R [GeV]", 200, 0.0, 2.0, "#Delta #phi", 60, 0.0, PI);
    h2_cluster_wdR.push_back(temp2);
    
    name = "h2_cluster_pi0_dphi_dR_c" + bin.str();
    Init2DHisto(temp2, name.c_str(), "#Delta R [GeV]", 200, 0.0, 2.0, "#Delta #phi", 60, 0.0, PI);
    h2_cluster_pi0_wdR.push_back(temp2);
    
    name = "h2_cluster_mix_dphi_dR_c" + bin.str();
    Init2DHisto(temp2, name.c_str(), "#Delta R [GeV]", 200, 0.0, 2.0, "#Delta #phi", 60, 0.0, PI);
    h2_cluster_mix_wdR.push_back(temp2);
    
    name = "h3_cluster_etot_dR_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{tot} [GeV]", 40, 0.0, 2.0,"#Delta R",200,0.0,2.0);
    h3_cluster_etot.push_back(temp3);
    
    name = "h3_cluster_etot_pi0_dR_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{tot} [GeV]", 40, 0.0, 2.0,"#Delta R",200,0.0,2.0);
    h3_cluster_pi0_etot.push_back(temp3);
    
    name = "h3_cluster_etot_mix_dR_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{tot} [GeV]", 40, 0.0, 2.0,"#Delta R",200,0.0,2.0);
    h3_cluster_mix_etot.push_back(temp3);
    
    name = "h3_cluster_dR_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{assoc} [GeV]", 40, 0.0, 10.0,"#Delta R",200,0.0,2.0);
    h3_cluster_dR.push_back(temp3);
    
    name = "h3_cluster_pi0_dR_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{assoc} [GeV]", 40, 0.0, 10.0,"#Delta R",200,0.0,2.0);
    h3_cluster_pi0_dR.push_back(temp3);
    
    name = "h3_cluster_mix_dR_c" + bin.str();
    Init3DHisto(temp3, name.c_str(), "p_{T} [GeV/c]", 60, 0.0, 15.0, "E_{assoc} [GeV]", 40, 0.0, 10.0,"#Delta R",200,0.0,2.0);
    h3_cluster_mix_dR.push_back(temp3);
    
    name = "h2_cluster_etot_dR_c" + bin.str();
    Init2DHisto(temp2, name.c_str(), "E_{tot}/p_{T}", 60, 0.0, 5.0, "p_{T} [GeV/c]", 60, 0.0, 15.0);
    h2_cluster_etot.push_back(temp2);
    
    name = "h2_cluster_etot_pi0_dR_c" + bin.str();
    Init2DHisto(temp2, name.c_str(), "E_{tot}/p_{T}", 60, 0.0, 5.0, "p_{T} [GeV/c]", 60, 0.0, 15.0);
    h2_cluster_pi0_etot.push_back(temp2);
    
    name = "h2_cluster_etot_mix_dR_c" + bin.str();
    Init2DHisto(temp2, name.c_str(), "E_{tot}/p_{T}", 60, 0.0, 5.0, "p_{T} [GeV/c]", 60, 0.0, 15.0);
    h2_cluster_mix_etot.push_back(temp2);
    
    if(DiagFlag){
      name = "h3_dphi_1sig_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
      h3_dphi_1sig.push_back(temp3);
      
      name = "h3_dphi_2sig_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
      h3_dphi_2sig.push_back(temp3);
      
      name = "h3_dphi_fold_1sig_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
      h3_dphi_fold_1sig.push_back(temp3);
      
      name = "h3_dphi_fold_2sig_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "p_{T, #gamma} [GeV/c]", 20, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 30, 0., PI);
      h3_dphi_fold_2sig.push_back(temp3);
      
      name = "h3_iso_acc_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "#phi [rad]", 60, -1*PI, PI, "#eta", 80, -0.4, 0.4,"p_{T} [GeV/c]",60,0.0,15.0);
      h3_iso_acc.push_back(temp3);
      
      name = "h3_iso_pi0_acc_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "#phi [rad]", 60, -1*PI, PI, "#eta", 80, -0.4, 0.4,"p_{T} [GeV/c]",60,0.0,15.0);
      h3_iso_pi0_acc.push_back(temp3);
      
      name = "h3_iso_mix_acc_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "#phi [rad]", 60, -1*PI, PI, "#eta", 80, -0.4, 0.4,"p_{T} [GeV/c]",60,0.0,15.0);
      h3_iso_mix_acc.push_back(temp3);
      
      name = "h3_dphi_pi0_bg_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 200, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
      h3_dphi_pi0_bg.push_back(temp3);
      name = "h3_dphi_pi0_bg_mix_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "p_{T, #pi^{0}} [GeV/c]", 200, 0.0, 20.0, "p_{T, h} [GeV/c]", 100, 0.0, 10.0, "#Delta#phi [rad]", 60, -PI/2, 3*PI/2);
      h3_dphi_pi0_bg_mix.push_back(temp3);
      name = "h3_mintrackdist_fg_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "track_dist[cm]", 100,0.,100., "trigger pt [GeV/c]", 20, 0.,20, "veto track p_{T} [GeV/c]", 40,0.,20.);
      h3_mintrackdist_fg.push_back(temp3);
      name = "h3_mintrackdist_bg_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "track_dist[cm]", 100,0.,100., "trigger pt [GeV/c]", 20, 0.,20, "veto track p_{T} [GeV/c]", 40,0.,20.);
      h3_mintrackdist_bg.push_back(temp3);
      name = "h3_EoverP_c" + bin.str();
      Init3DHisto(temp3, name.c_str(), "p_{T} [GeV/c]", 100, 0., 10., "E/p", 50, 0., 1., "#Delta#phi [rad]", 30, 0.0, PI);
      h3_EoverP.push_back(temp3);
    }
  }
  
  name = "h2_pi0mass_PbGl";
  Init2DHisto(h2_pi0mass_PbGl, name.c_str(), "p_{T} [GeV/c]", 100, 0.0, 20.0, "m [GeV/c^{2}]", 200, 0.0, 0.6);
  name = "h2_pi0mass_PbSc";
  Init2DHisto(h2_pi0mass_PbSc, name.c_str(), "p_{T} [GeV/c]", 100, 0.0, 20.0, "m [GeV/c^{2}]", 200, 0.0, 0.6);

  double pos_array[100];
  for (int i = 0; i < 100; i++) {pos_array[i] = i;}

  for(int ias=0; ias<N_ARMSECT; ias++){
    bin.str("");
    bin << ias;
    name = "h2_pi0mass_as" + bin.str();
    Init2DHisto(h2_pi0mass_as[ias], name.c_str(), "p_{T} [GeV/c]", 100, 0.0, 20.0, "m [GeV/c^{2}]", 200, 0.0, 0.6);
    if(DiagFlag){
      name = "h2_pi0_bg_as" + bin.str();
      Init2DHisto(h2_pi0_bg_as[ias], name.c_str(), "p_{T} [GeV/c]", 100, 0.0, 20.0, "m [GeV/c^{2}]", 200, 0.0, 0.6);
    }

    name = "nhit_as" + bin.str();
    int ny_max, nz_max;
    if (IsPbGl(ias)) { ny_max = N_YPOS_PBGL;  nz_max = N_ZPOS_PBGL; }
    else             { ny_max = N_YPOS_PBSC;  nz_max = N_ZPOS_PBSC; }

    int pt_max=25;

    Init3DHisto(h3_nhit[ias], name.c_str(), "ecore [GeV]", pt_max, pos_array, "ypos", ny_max, pos_array, "zpos", nz_max, pos_array);
    name = "nhit_mywarn_as" + bin.str();
    Init3DHisto(h3_nhit_mywarn[ias], name.c_str(), "ecore [GeV]", pt_max, pos_array, "ypos", ny_max, pos_array, "zpos", nz_max, pos_array);

    name = "h1_chi2_as" + bin.str();
    Init1DHisto(h1_chi2[ias], name.c_str(), "#chi^{2}", 20, 0., 10.);
    name = "h3_disp_as" + bin.str();
    Init3DHisto(h3_disp[ias], name.c_str(), "dispMax", 20, 0., 2., "dispCut", 20, 0., 2., "pass", 2,0,2);
    name = "h1_trigger_as" + bin.str();
    Init1DHisto(h1_trigger[ias], name.c_str(), "ERT Supermodule", 50, 0., 50.);
  }

  if(DiagFlag){
    name = "h2_pi0_bg";
    Init2DHisto(h2_pi0_bg, name.c_str(), "p_{T} [GeV/c]", 100, 0.0, 20.0, "m [GeV/c^{2}]", 200, 0.0, 0.6);
    name = "h2_pi0_bg_PbGl";
    Init2DHisto(h2_pi0_bg_PbGl, name.c_str(), "p_{T} [GeV/c]", 100, 0.0, 20.0, "m [GeV/c^{2}]", 200, 0.0, 0.6);
    name = "h2_pi0_bg_PbSc";
    Init2DHisto(h2_pi0_bg_PbSc, name.c_str(), "p_{T} [GeV/c]", 100, 0.0, 20.0, "m [GeV/c^{2}]", 200, 0.0, 0.6);

    Init2DHisto(h2_emctrkdphi_bf, "h2_emctrkdphi_bf", "p_{T} [GeV/c]", 30, 0.0, 15.0, "emctrkdphi [rad]", 60, -0.2, 0.2);
    Init2DHisto(h2_emctrkdz_bf, "h2_emctrkdz_bf", "p_{T} [GeV/c]", 30, 0.0, 15.0, "emctrkdz [cm]", 70, -70.0, 70.0);
    Init2DHisto(h2_emcpc3dphi_bf, "h2_emcpc3dphi_bf", "p_{T} [GeV/c]", 30, 0.0, 15.0, "emcpc3dphi [rad]", 60, -0.2, 0.2);
    Init2DHisto(h2_emcpc3dz_bf, "h2_emcpc3dz_bf", "p_{T} [GeV/c]", 30, 0.0, 15.0, "emcpc3dz [cm]", 70, -70.0, 70.0);
    
    Init2DHisto(h2_emctrkdphi_aft, "h2_emctrkdphi_aft", "p_{T} [GeV/c]", 30, 0.0, 15.0, "emctrkdphi [rad]", 60, -0.2, 0.2);
    Init2DHisto(h2_emctrkdz_aft, "h2_emctrkdz_aft", "p_{T} [GeV/c]", 30, 0.0, 15.0, "emctrkdz [cm]", 70, -70.0, 70.0);
    Init2DHisto(h2_emcpc3dphi_aft, "h2_emcpc3dphi_aft", "p_{T} [GeV/c]", 30, 0.0, 15.0, "emcpc3dphi [rad]", 60, -0.2, 0.2);
    Init2DHisto(h2_emcpc3dz_aft, "h2_emcpc3dz_aft", "p_{T} [GeV/c]", 30, 0.0, 15.0, "emcpc3dz [cm]", 70, -70.0, 70.0);
    
    //track information-before cut
    Init2DHisto(h2_pc3sdphi_bf, "h2_pc3sdphi_bf", "p_{T} [GeV/c]", 35, 0.0, 7.0, "pc3sdphi [rad]", 40, -10.0, 10.0);
    Init2DHisto(h2_pc3sdz_bf, "h2_pc3sdz_bf", "p_{T} [GeV/c]", 35, 0.0, 7.0, "pc3sdz [cm]", 40, -10.0, 10.0);
    //-after cut
    Init2DHisto(h2_pc3sdphi_aft, "h2_pc3sdphi_aft", "p_{T} [GeV/c]", 35, 0.0, 7.0, "pc3sdphi [rad]", 40, -10.0, 10.0);
    Init2DHisto(h2_pc3sdz_aft, "h2_pc3sdz_aft", "p_{T} [GeV/c]", 35, 0.0, 7.0, "pc3sdz [cm]", 40, -10.0, 10.0);
    
    //DC acceptance
    Init2DHisto(h2_phi_zed_bf, "h2_phi_zed_bf", "zed [cm]", 160, -80, 80,"#phi [rad]", 100, -1.0, 4.0);
    Init2DHisto(h2_phi_zed_aft, "h2_phi_zed_aft", "zed [cm]", 160, -80, 80,"#phi [rad]", 100, -1.0, 4.0);
    
    //PC3 acceptance
    Init2DHisto(h2_ppc3_west_bf, "h2_ppc3_west_bf", "ppc3z [cm]", 400, -200, 200, "ppc3y [cm]", 750, -300, 450);
    Init2DHisto(h2_ppc3_east_bf, "h2_ppc3_east_bf", "ppc3z [cm]", 400, -200, 200, "ppc3y [cm]", 750, -300, 450);
    Init2DHisto(h2_ppc3_west_aft, "h2_ppc3_west_aft", "ppc3z [cm]", 400, -200, 200, "ppc3y [cm]", 750, -300, 450);
    Init2DHisto(h2_ppc3_east_aft, "h2_ppc3_east_aft", "ppc3z [cm]", 400, -200, 200, "ppc3y [cm]", 750, -300, 450);
    //PC1 acceptance
    Init2DHisto(h2_ppc1_west_bf, "h2_ppc1_west_bf", "ppc1z [cm]", 200, -100, 100, "ppc1y [cm]", 380, -150, 230);
    Init2DHisto(h2_ppc1_east_bf, "h2_ppc1_east_bf", "ppc1z [cm]", 200, -100, 100, "ppc1y [cm]", 380, -150, 230);
    Init2DHisto(h2_ppc1_west_aft, "h2_ppc1_west_aft", "ppc1z [cm]", 200, -100, 100, "ppc1y [cm]", 380, -150, 230);
    Init2DHisto(h2_ppc1_east_aft, "h2_ppc1_east_aft", "ppc1z [cm]", 200, -100, 100, "ppc1y [cm]", 380, -150, 230);
    Init2DHisto(h2_EoverPvspt, "h2_EoverPvspt", "p_{T} [GeV/c]", 100, 0., 10., "E/p", 50, 0., 1.);
    
    InitPhotonCutChecker("photon_cut_check", photon_cut_check);
    InitHadronCutChecker("hadron_cut_check", hadron_cut_check);

    Init3DHisto(h3_mintrackdist_fg_allcent, "h3_mintrackdist_fg_allcent", "track_dist[cm]", 100,0.,100., "trigger pt [GeV/c]", 20, 0.,20, "veto track p_{T} [GeV/c]", 40,0.,20.);
    Init3DHisto(h3_mintrackdist_bg_allcent, "h3_mintrackdist_bg_allcent", "track_dist[cm]", 100,0.,100., "trigger pt [GeV/c]", 20, 0.,20, "veto track p_{T} [GeV/c]", 40,0.,20.);
  }

  Init1DHisto(h1_zvertex, "h1_zvertex","z [cm]",60,-30,30);
  Init1DHisto(h1_trig_pt_inc_tot, "h1_trig_pt_inc_tot","p_{T} [GeV/c]",200,0.,20.);
  Init1DHisto(h1_trig_pt_pi0_tot, "h1_trig_pt_pi0_tot","p_{T} [GeV/c]",200,0.,20.);
  Init1DHisto(h1_trig_pt_dec_tot, "h1_trig_pt_dec_tot","p_{T} bin",5,-0.5,4.5);
  Init1DHisto(h1_trig_pt_inc_mix_tot, "h1_trig_pt_inc_mix_tot","p_{T} [GeV/c]",200,0.,20.);
  Init1DHisto(h1_trig_pt_pi0_mix_tot, "h1_trig_pt_pi0_mix_tot","p_{T} [GeV/c]",200,0.,20.);
  Init1DHisto(h1_trig_pt_dec_mix_tot, "h1_trig_pt_dec_mix_tot","p_{T} bin",5,-0.5,4.5);
  Init1DHisto(h1_part_pt_tot, "h1_part_pt_tot","p_{T} [GeV/c]",100,0.,10.);
  Init1DHisto(h1_part_pt_05sig, "h1_part_pt_05sig","p_{T} [GeV/c]",100,0.,10.);
  Init1DHisto(h1_part_pt_1sig, "h1_part_pt_1sig","p_{T} [GeV/c]",100,0.,10.);
  Init1DHisto(h1_part_pt_15sig, "h1_part_pt_15sig","p_{T} [GeV/c]",100,0.,10.);
  Init1DHisto(h1_part_pt_2sig, "h1_part_pt_2sig","p_{T} [GeV/c]",100,0.,10.);
  Init1DHisto(h1_part_pt_3sig, "h1_part_pt_3sig","p_{T} [GeV/c]",100,0.,10.);
  Init2DHisto(h2_ptvscent_trig_inc, "h2_ptvscent_trig_inc","p_{T} [GeV/c]",200,0.,20., "centrality [%]",101,0,100);
  Init2DHisto(h2_ptvscent_part, "h2_ptvscent_part","p_{T} [GeV/c]",100,0.,10., "centrality [%]",101,0,100);
  Init2DHisto(h2_ptvscent_trig_pi0, "h2_ptvscent_trig_pi0","p_{T} [GeV/c]",200,0.,20., "centrality [%]",101,0,100);
  Init2DHisto(h2_bfpaircut_inc,"h2_bfpaircut_inc","p_{T,#gamma} bin",4,-0.5,3.5,"p_{T,h} bin",5,-0.5,4.5);
  Init2DHisto(h2_aftpaircut_inc,"h2_aftpaircut_inc","p_{T,#gamma} bin",4,-0.5,3.5,"p_{T,h} bin",5,-0.5,4.5);
  Init2DHisto(h2_bfpaircut_pi0,"h2_bfpaircut_pi0","p_{T,#gamma} bin",4,-0.5,3.5,"p_{T,h} bin",5,-0.5,4.5);
  Init2DHisto(h2_aftpaircut_pi0,"h2_aftpaircut_pi0","p_{T,#gamma} bin",4,-0.5,3.5,"p_{T,h} bin",5,-0.5,4.5);

  h2_dphi_dec = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphi_dec_fold = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphixi_dec = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphixi_dec_fold = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphizt_dec = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphizt_dec_fold = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphi_dec_iso_fold = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphi_dec_mix = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphi_dec_mix_fold = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphixi_dec_mix = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphixi_dec_mix_fold = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphizt_dec_mix = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphizt_dec_mix_fold = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );
  h2_dphi_dec_mix_iso_fold = vector<vector<TH2F*> > (NCBINS, vector<TH2F*>() );

  for(int ic=0; ic<NCBINS; ic++){
    for(int ipt=0; ipt<5; ipt++){
      bin.str("");
      bin << ipt <<"_c"<<ic;
      name = "h2_dphi_dec_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec[ic].push_back(temp2);

      name = "h2_dphi_dec_iso_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec_iso[ic].push_back(temp2);

      name = "h2_dphi_dec_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec_fold[ic].push_back(temp2);

      name = "h2_dphi_dec_iso_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec_iso_fold[ic].push_back(temp2);

      name = "h2_dphixi_dec_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "xi",25,-1.0,4.0);
      h2_dphixi_dec[ic].push_back(temp2);

      name = "h2_dphixi_dec_iso_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "xi",25,-1.0,4.0);
      h2_dphixi_dec_iso[ic].push_back(temp2);

      name = "h2_dphixi_dec_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "xi",25,-1.0,4.0);
      h2_dphixi_dec_fold[ic].push_back(temp2);

      name = "h2_dphixi_dec_iso_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "xi",25,-1.0,4.0);
      h2_dphixi_dec_iso_fold[ic].push_back(temp2);

      name = "h2_dphizt_dec_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "zt",40,0.0,2.0);
      h2_dphizt_dec[ic].push_back(temp2);

      name = "h2_dphizt_dec_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "zt",40,0.0,2.0);
      h2_dphizt_dec_fold[ic].push_back(temp2);

      name = "h2_dphi_dec_iso_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec_iso_fold[ic].push_back(temp2);

      name = "h2_dphi_dec_mix_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec_mix[ic].push_back(temp2);

      name = "h2_dphi_dec_iso_mix_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec_iso_mix[ic].push_back(temp2);

      name = "h2_dphi_dec_mix_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec_mix_fold[ic].push_back(temp2);

      name = "h2_dphi_dec_iso_mix_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "p_{T, h} [GeV/c]",100,0.0,10.0);
      h2_dphi_dec_iso_mix_fold[ic].push_back(temp2);

      name = "h2_dphixi_dec_mix_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "xi",25,-1.0,4.0);
      h2_dphixi_dec_mix[ic].push_back(temp2);

      name = "h2_dphixi_dec_iso_mix_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "xi",25,-1.0,4.0);
      h2_dphixi_dec_iso_mix[ic].push_back(temp2);

      name = "h2_dphixi_dec_mix_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "xi",25,-1.0,4.0);
      h2_dphixi_dec_mix_fold[ic].push_back(temp2);

      name = "h2_dphixi_dec_iso_mix_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "xi",25,-1.0,4.0);
      h2_dphixi_dec_iso_mix_fold[ic].push_back(temp2);

      name = "h2_dphizt_dec_mix_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 60, -PI/2, 3*PI/2, "zt",40,0.0,2.0);
      h2_dphizt_dec_mix[ic].push_back(temp2);

      name = "h2_dphizt_dec_mix_fold_p" + bin.str();
      Init2DHisto(temp2, name.c_str(), "#Delta#phi [rad]", 30, 0., PI, "zt",40,0.0,2.0);
      h2_dphizt_dec_mix_fold[ic].push_back(temp2);
    }
  }

  Init1DHisto(h1_n0,"h1_n0","n0",11,-1,10);

  //filltime debugging histos
  if( DiagFlag ) {
    cout<<"init new histos~~~~"<<endl;
    //REAL pairs
    name = "h2_dphi";
    Init2DHisto(h2_dphi,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI);
    name = "h2_dphi_accw";
    Init2DHisto(h2_dphi_accw,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI);
    name = "h3_dphi_accw";
    Init3DHisto(h3_dphi_accw,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI, "weight", 80, 0., 80.);

    name = "h2_dphi_xi";
    Init2DHisto(h2_dphi_xi,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI);
    name = "h2_dphi_accw_xi";
    Init2DHisto(h2_dphi_accw_xi,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI);
    name = "h3_dphi_accw_xi";
    Init3DHisto(h3_dphi_accw_xi,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI, "weight", 80, 0., 80. );

    //MIX pairs
    name = "h2_dphi_mix";
    Init2DHisto(h2_dphi_mix,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI);
    name = "h2_dphi_accw_mix";
    Init2DHisto(h2_dphi_accw_mix,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI);
    name = "h3_dphi_accw_mix";
    Init3DHisto(h3_dphi_accw_mix,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI, "weight", 80, 0., 80.);

    name = "h2_dphi_xi_mix";
    Init2DHisto(h2_dphi_xi_mix,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI);
    name = "h2_dphi_accw_xi_mix";
    Init2DHisto(h2_dphi_accw_xi_mix,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI);
    name = "h3_dphi_accw_xi_mix";
    Init3DHisto(h3_dphi_accw_xi_mix,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI, "weight", 80, 0., 80. );

    //REALPI pairs
    name = "h2_dphi_pi0";
    Init2DHisto(h2_dphi_pi0,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI);
    name = "h2_dphi_accw_pi0";
    Init2DHisto(h2_dphi_accw_pi0,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI);
    name = "h3_dphi_accw_pi0";
    Init3DHisto(h3_dphi_accw_pi0,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI, "weight", 80, 0., 80.);

    name = "h2_dphi_xi_pi0";
    Init2DHisto(h2_dphi_xi_pi0,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI);
    name = "h2_dphi_accw_xi_pi0";
    Init2DHisto(h2_dphi_accw_xi_pi0,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI);
    name = "h3_dphi_accw_xi_pi0";
    Init3DHisto(h3_dphi_accw_xi_pi0,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI, "weight", 80, 0., 80. );

    //MIXPI pairs
    name = "h2_dphi_pi0_mix";
    Init2DHisto(h2_dphi_pi0_mix,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI);
    name = "h2_dphi_accw_pi0_mix";
    Init2DHisto(h2_dphi_accw_pi0_mix,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI);
    name = "h3_dphi_accw_pi0_mix";
    Init3DHisto(h3_dphi_accw_pi0_mix,name.c_str(),"p_{T}^{h}", 100, 0., 10, "#Delta#phi",30, 0., PI, "weight", 80, 0., 80.);

    name = "h2_dphi_xi_pi0_mix";
    Init2DHisto(h2_dphi_xi_pi0_mix,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI);
    name = "h2_dphi_accw_xi_pi0_mix";
    Init2DHisto(h2_dphi_accw_xi_pi0_mix,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI);
    name = "h3_dphi_accw_xi_pi0_mix";
    Init3DHisto(h3_dphi_accw_xi_pi0_mix,name.c_str(),"#xi", 25, -1.0, 4.0, "#Delta#phi",30, 0., PI, "weight", 80, 0., 80. );

    name = "h2_partpt_xi";
    Init2DHisto(h2_partpt_xi, name.c_str(),"p_{T}^{h}", 100, 0., 10, "#xi", 25, -1.0, 4.0);
    name = "h2_partpt_xi_mix";
    Init2DHisto(h2_partpt_xi_mix, name.c_str(),"p_{T}^{h}", 100, 0., 10, "#xi", 25, -1.0, 4.0);
    name = "h2_partpt_xi_pi0";
    Init2DHisto(h2_partpt_xi_pi0, name.c_str(),"p_{T}^{h}", 100, 0., 10, "#xi", 25, -1.0, 4.0);
    name = "h2_partpt_xi_pi0_mix";
    Init2DHisto(h2_partpt_xi_pi0_mix, name.c_str(),"p_{T}^{h}", 100, 0., 10, "#xi", 25, -1.0, 4.0);
  }

  name = "h3_pt_phi_eta_clus";
  Init3DHisto(h3_pt_phi_eta_clus, name.c_str(), "p_{T}^{#gamma} [GeV/c]", 200, 0., 20., "#phi [rad]", 100, -1.0, 4.0, "#eta [rad]", 80, -0.4, 0.4);
  name = "h3_pt_phi_eta_trk";
  Init3DHisto(h3_pt_phi_eta_trk, name.c_str(), "p_{T}^{h} [GeV/c]",100, 0., 10., "#phi [rad]", 100, -1.0, 4.0, "#eta [rad]", 80, -0.4, 0.4);

  string acc_filename = toad_loader->location(_accfilename.c_str());
  cout<<"get acc_filename: "<<acc_filename.c_str()<<endl;
  GetAcceptanceWeights(acc_filename.c_str());

  atree = new AMixingTree();
  atree->SetTriggerBranches();
  atree->SetPartnerBranches();

  return 0;
}

int Correlation::ResetEvent(PHCompositeNode* topNode)
{
  emcclustercontainer = 0;
  particle = 0;
  return 0;
}

int Correlation::process_event(PHCompositeNode* topNode)
{

  if( evt%5000==0 )
    cout << "Correlation: at event " << evt << endl;
  evt++;
  
  global = findNode::getClass<PHGlobal>(topNode,"PHGlobal");
  if(!global)
  {
      cout << PHWHERE << "findNode::GetClass: PHGlobal not in Node Tree" << endl;
      //return ABORTEVENT;
  }
  
  emcclustercontainer = findNode::getClass<emcClusterContainer>(topNode, "emcClusterContainer");
  if(!emcclustercontainer)
    {
      cout << PHWHERE << "emcclustercontainter not in Node Tree" << endl;
      //return EVENT_OK;
      return -2;
    }

  particle = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if(!particle) {
    cout << PHWHERE << "PHCentralTrack not in Node Tree" << endl;
    //return EVENT_OK;
    return -2;
  }

  RunHeader* runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if(!runheader){
    cout << PHWHERE << " RunHeader not found"  << endl;
    //return EVENT_OK;
    return 0;
  }

  int RunNumber = runheader->get_RunNumber();
  data_set = GetDataSet(RunNumber);
  if( data_set == INVALID ) {
    cout << PHWHERE << "Warning: run# outside valid range! Add this Run to DataSet options" << endl;
    return EVENT_OK;
  }

  EventHeader *evtsync = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if(!evtsync){
    cout << PHWHERE << "CombinedSimple:: EventHeader not in Node Tree" << endl;
    return 0;
  }
  
  event_z = global->getBbcZVertex();
  event_c = global->getCentrality();
  
  if (verbosity > 2) cout << "BBC vertex: "<< event_z <<" Centrality: "<< event_c << endl;

  if (event_z < -1*zVertexCut || event_z > zVertexCut) {
    if(verbosity > 0) cout<<"Outside vertex range: z = "<<event_z<<endl;
    //return EVENT_OK;
    return 0;
  }

  if(runheader->get_currentCentral()>0) fieldPolarity = 1.0;
  else fieldPolarity = -1.0;

  // For looking at triggered data
  ErtOut* ertout = findNode::getClass<ErtOut>(topNode, "ErtOut");

  AEvent aevent;
  MakeEventObject(global, &aevent);
  int vtxbin = aevent.GetVertexBin();
  int centbin = aevent.GetCentralityBin();
  int centrality = (int)event_c;
  cbin = GetCentBin(centrality);
  if(verbosity > 2) cout<<"vtxbin = "<<vtxbin<<"; centbin = "<<centbin<< "; cbin = "<<cbin<<endl;
  
  int nclus = emcclustercontainer->size();
  int npart = particle->get_npart();
  if (verbosity > 2)   cout<<"At event #" <<evt<<": nclus = "<<nclus<<" npart = "<<npart<<endl;

  h1_zvertex->Fill(event_z);

  if ( (vtxbin < 0) || (cbin < 0) || (cbin > (NCBINS-1))) {
    if(verbosity > 0) cout<<"Invalid ZVertex/centrality bin! vtxbin = "<<vtxbin<<"; cbin = "<<cbin<<"; event_c = "<<event_c<<endl;
    return 0;
  }
  h1_centrality[cbin]->Fill(1);
  
  //*****************************************************************************
  //*  Select track sample (quality > 7) to compute track min dist later        *
  //*****************************************************************************
  for(int ipart = 0; ipart < npart; ipart++){
    if(particle->get_the0(ipart) < -99) continue;
    if(particle->get_quality(ipart) <= 7) continue;
    ATrack atrack;
    MakeTrackObject(particle, ipart, &atrack);
    lessqualtrk_vector.push_back(atrack.clone());
  }
  if(verbosity > 1) cout<<"NEW: lessqualtrk_vector size: "<<lessqualtrk_vector.size()<<endl;

  //*************************************************
  //*  Selecting inclusive photon sample            *
  //*************************************************
  for(int iclus = 0; iclus < nclus; iclus++){
    emcClusterContent* clus = emcclustercontainer->getCluster(iclus);
    ACluster acluster;
    MakeClusterObject(clus, &acluster);
    double ecore = acluster.GetEcore();
    if( ecore < 0.1 ) continue;
    int armsect = acluster.GetArmSect();
    int ypos = acluster.GetIypos();
    int zpos = acluster.GetIzpos();
    h3_nhit[armsect]->Fill(ecore, ypos, zpos);
    if( !IsGoodTower(&acluster) ) continue;
    if( !IsGoodCluster(&acluster,data_set) ) continue;
    // Check trigger bit - for triggered data
    if( data_set == Run8dAu ) {
      int ert_sm = (acluster.GetArm()==1&&acluster.GetSec()<2) ? acluster.GetIypos()/12*8 + acluster.GetIzpos()/12 : acluster.GetIypos()/12*6 + acluster.GetIzpos()/12;
      if( verbosity > 1 ) cout << "Photon in SM " << ert_sm << endl;
      int sm_fired = 0;
      // set trigger bit if fired ERT4x4 a,b, or c
      if( ertout ) 
        if( ertout->get_ERTbit(0, acluster.GetArm(), acluster.GetSec(), ert_sm) ||
            ertout->get_ERTbit(1, acluster.GetArm(), acluster.GetSec(), ert_sm) )
          //  ertout->get_ERTbit(2, acluster.GetArm(), acluster.GetSec(), ert_sm) ) 
          sm_fired = 1;
      //if( (armsect==2 && ert_sm==17) || (armsect==3 && (ert_sm==4||ert_sm==10||ert_sm==15||ert_sm==16)) || (armsect==4 && ert_sm==30) || (armsect==7 && ert_sm==15) ) sm_fired = 0;
      if( sm_fired ) {
        acluster.SetTrigger(true);
        if( verbosity > 3 ) cout << "Photon in sector " << armsect << ", SM " << ert_sm << " fired ERT trigger!" << endl;
      }
      else acluster.SetTrigger(false);
    }
    else acluster.SetTrigger(true); // For other data sets just set trigger status to true by default
    if( ecore<5.0 || acluster.IsTrigger() )
      h3_nhit_mywarn[armsect]->Fill(ecore, ypos, zpos);
    
    //filling a eta phi 2d histogram
    float pt = acluster.Pt();
    float phi = acluster.Phi();
    float theta = acluster.Theta();
    float eta = -log(tan(theta/2));
    h3_pt_phi_eta_clus->Fill(pt,phi,eta);

    all_clus_vector.push_back(acluster.clone());
  }
  if( verbosity > 2 ) cout << "Found " << all_clus_vector.size() << " raw clusters" << endl;

  //******************************************************
  //*   Select decay photon sample and construct pi0s.   *
  //******************************************************
  if( all_clus_vector.size() > 0 ) MakePi0s(all_clus_vector,centrality,event_z,data_set);
  if(verbosity > 1)/*if (pi0_vector.size())*/ cout<<"evt#: "<<evt<<"; NEW: pi0_vector size: "<<pi0_vector.size()<<endl;
  
  for(unsigned int iclus = 0; iclus < all_clus_vector.size(); iclus++){
    double ecore = all_clus_vector[iclus]->GetEcore();
    
    if( ecore < 0.3 ) continue;

    float cluster_pt = all_clus_vector[iclus]->Pt();
    if((cluster_pt <= photon_pt_min) || (cluster_pt >= photon_pt_max)) continue;
    //cout<<"vtxbin = "<<vtxbin<<"; centbin = "<<centbin<< "; cbin = "<<cbin<<endl;

    if( verbosity > 1 ) cout << "Found cluster with ecore = " << ecore << endl;

    if(DiagFlag) FillClusterQAHistos(0,cluster_pt,all_clus_vector[iclus]->GetEmctrkdphi(),all_clus_vector[iclus]->GetEmctrkdz(),all_clus_vector[iclus]->GetEmcpc3dphi(),all_clus_vector[iclus]->GetEmcpc3dz());
    
    if( verbosity > 1 ) cout << "Checking charge veto cut...    " << endl;
    float mindist = 99.;
    float ptofvetotrack = -1.0;
    VetoTracks(all_clus_vector[iclus],lessqualtrk_vector,mindist,ptofvetotrack);
    if(DiagFlag){
      h3_mintrackdist_fg_allcent->Fill(mindist,all_clus_vector[iclus]->Pt(), ptofvetotrack);
      h3_mintrackdist_fg[cbin]->Fill(mindist,all_clus_vector[iclus]->Pt(), ptofvetotrack);
    }
    if(mindist<8.0) continue;
    
    if( verbosity > 1 ) cout << "Cluster passed!" << endl;
    if(ecore < photon_ecore_min) continue;
    if(DiagFlag) FillClusterQAHistos(1,cluster_pt,all_clus_vector[iclus]->GetEmctrkdphi(),all_clus_vector[iclus]->GetEmctrkdz(),all_clus_vector[iclus]->GetEmcpc3dphi(),all_clus_vector[iclus]->GetEmcpc3dz());
 
    int armsect = all_clus_vector[iclus]->GetArmSect();
    emcClusterContent* clus = all_clus_vector[iclus]->GetemcClusterContent(); 

    if(!Chi2Cut(clus,event_z,h1_chi2[armsect],h3_disp[armsect])) { continue; } //chi2 cut on inclusive photons

    if( !all_clus_vector[iclus]->IsTrigger() ) { continue; }

    int ert_sm = (all_clus_vector[iclus]->GetArm()==1&&all_clus_vector[iclus]->GetSec()<2) ? all_clus_vector[iclus]->GetIypos()/12*8 + all_clus_vector[iclus]->GetIzpos()/12 : all_clus_vector[iclus]->GetIypos()/12*6 + all_clus_vector[iclus]->GetIzpos()/12;
    h1_trigger[armsect]->Fill(ert_sm);

    //fiducial cut on inclusive photon
    float x1 = all_clus_vector[iclus]->GetX();
    float y1 = all_clus_vector[iclus]->GetY();
    float z1 = all_clus_vector[iclus]->GetZ();

    float zemcsub = z1*510.0/sqrt(x1*x1+y1*y1);
    if(fabs(zemcsub)>155.0) continue;
    if( all_clus_vector[iclus]->IsTagged() )
      SetIso(all_clus_vector[iclus],lessqualtrk_vector,all_clus_vector,Rcut,econe_min[cbin],h3_cluster_pi0_dR[cbin],h3_cluster_pi0_etot[cbin],h2_cluster_pi0_wdR[cbin],h2_cluster_pi0_etot[cbin],h3_iso_pi0_acc[cbin]);
    SetIso(all_clus_vector[iclus],lessqualtrk_vector,all_clus_vector,Rcut,econe_min[cbin],h3_cluster_dR[cbin],h3_cluster_etot[cbin],h2_cluster_wdR[cbin],h2_cluster_etot[cbin],h3_iso_acc[cbin]);
    if( verbosity > 0 ) cout << "Photon isTagged = " << all_clus_vector[iclus]->IsTagged() << " and IsIso = " << all_clus_vector[iclus]->IsIso() << endl;
    
    h1_trig_pt_all[cbin]->Fill(cluster_pt); // Keep track of all photons before tagging rejection (for dAu or p+p)
    if( all_clus_vector[iclus]->IsIso() ) h1_trig_pt_all_iso[cbin]->Fill(cluster_pt); // Keep track of all photons before tagging rejection (for dAu or p+p)
    // In Run8 reject all tagged photons prior to making pairs...
    if( data_set == Run8dAu && all_clus_vector[iclus]->IsTagged() ) continue;
    
    h1_trig_pt_inc[cbin]->Fill(cluster_pt);
    if( all_clus_vector[iclus]->IsIso() ) h1_trig_pt_inc_iso[cbin]->Fill(cluster_pt);
    h1_trig_pt_inc_tot->Fill(cluster_pt);
    h2_ptvscent_trig_inc->Fill(cluster_pt, centrality);
    //h2_trig_pt_zvtx_inc->Fill(cluster_pt,vtxbin);
    
    clus_vector.push_back(all_clus_vector[iclus]->clone());
    atree->SetTriggerData(cluster_pt,all_clus_vector[iclus]->Phi(),all_clus_vector[iclus]->Eta(),all_clus_vector[iclus]->E(),all_clus_vector[iclus]->GetX(),all_clus_vector[iclus]->GetY(),all_clus_vector[iclus]->GetZ(),all_clus_vector[iclus]->IsIso(),pi0_vector.size()+clus_vector.size()-1);
  }
  if(verbosity > 1) cout<<"NEW: clus_vector size: "<<clus_vector.size()<<endl;
  
  //*****************************************
  //*   Associate Hadrons                   *
  //*****************************************
  int ntrack = 0;
  //if( data_set == Run11AuAu && useVtx ) ntrack = nsvxpart;
  /* else*/ ntrack = npart;

  for(int ipart = 0; ipart < ntrack; ipart++){
    ATrack atrack;
    // if(useVtx){
    //   svxcnttrk = svxcntlist->getCentralTrack(ipart);  
    //   int icnttrk = svxcnttrk->getDchIndex();
    //   MakeTrackObject(particle, icnttrk, &atrack);
    // }
    /*else*/ MakeTrackObject(particle, ipart, &atrack);
    float trk_pt = atrack.Pt();
    //    float charge = atrack.GetCharge();

    // if(useVtx){
    //   ConversionVeto vetofunc;
    //   if(!vetofunc.calculate(fieldPolarity, trk_pt, charge, svxcnttrk, svxcluslist)) continue;
    // }
    
    if (atrack.GetPhi() <= -99) continue;    
    if(trk_pt < hadron_pt_min || trk_pt > hadron_pt_max) continue;
    
    if(DiagFlag){
      FillTrackQAHistos(0,atrack.GetPpc3x(),atrack.GetPpc3z(),atrack.GetPpc3y(),atrack.GetPpc1z(),atrack.GetPpc1y(),atrack.GetZed(),atrack.GetPhiD(),trk_pt,atrack.GetQuality(),atrack.GetN0(),atrack.GetPc3sdphi(),atrack.GetPc3sdz());
    }
    if((atrack.GetQuality()!=31) && (atrack.GetQuality()!=63)) continue;

    int trk_n0 = atrack.GetN0();
    if(trk_n0 < 0) trk_n0 = -1;
    h1_n0->Fill(trk_n0);
    if(trk_pt < 5.0 && trk_n0 >= 0) continue;    
    
    if(DiagFlag){
      FillTrackQAHistos(1,atrack.GetPpc3x(),atrack.GetPpc3z(),atrack.GetPpc3y(),atrack.GetPpc1z(),atrack.GetPpc1y(),atrack.GetZed(),atrack.GetPhiD(),trk_pt,atrack.GetQuality(),atrack.GetN0(),atrack.GetPc3sdphi(),atrack.GetPc3sdz());
    }
    //********make different matching cuts on hadron samples*************
    if(data_set == Run11AuAu){
      if(PassMatchingCut(&atrack,0.5,0.5)) trk_vector_05sig.push_back(atrack.clone());
      if(PassMatchingCut(&atrack,1.0,1.0)) trk_vector_1sig.push_back(atrack.clone());
      if(PassMatchingCut(&atrack,1.5,1.5)) trk_vector_15sig.push_back(atrack.clone());
      if(PassMatchingCut(&atrack,2.0,2.0)) trk_vector_2sig.push_back(atrack.clone());
      if(PassMatchingCut(&atrack,3.0,3.0)) trk_vector_3sig.push_back(atrack.clone());
    }
    else{
      if(PassMatchingCut(&atrack,0.5)) trk_vector_05sig.push_back(atrack.clone());
      if(PassMatchingCut(&atrack,1.0)) trk_vector_1sig.push_back(atrack.clone());
      if(PassMatchingCut(&atrack,1.5)) trk_vector_15sig.push_back(atrack.clone());
      if(PassMatchingCut(&atrack,2.0)) trk_vector_2sig.push_back(atrack.clone());
      if(PassMatchingCut(&atrack,3.0)) trk_vector_3sig.push_back(atrack.clone());
    }
    if(!IsGoodTrack(&atrack,data_set)) continue; 
    if(DiagFlag){
      h2_pc3sdphi_aft->Fill(trk_pt, atrack.GetPc3sdphi());
      h2_pc3sdz_aft->Fill(trk_pt, atrack.GetPc3sdz());
    } 
    h1_part_pt[cbin]->Fill(trk_pt);
    h1_part_pt_tot->Fill(trk_pt); 
    h2_ptvscent_part->Fill(trk_pt, centrality);   
    if(DiagFlag){      
      h2_EoverPvspt->Fill(trk_pt, atrack.GetEcore()/trk_pt);
    }
    float theta = atrack.GetTheta();
    float eta = -log(tan(theta/2));
    //float phi = atrack.GetPhi();
    // if(trk_pt>5.0)
    //   if((phi > 2.15 && phi < 2.4) && (eta > -0.38 && eta < 0.02)) continue;
    h3_pt_phi_eta_trk->Fill(trk_pt,atrack.GetPhi(),eta);

    trk_vector.push_back(atrack.clone());
    if( data_set != Run8dAu )
      atree->SetPartnerData(atrack.Pt(),atrack.Phi(),atrack.Eta(),atrack.E(),atrack.GetPemcx(),atrack.GetPemcy(),atrack.GetPemcz(),trk_vector.size()-1);

  }
  if(verbosity > 1) cout<<"NEW: trk_vector size: "<<trk_vector.size()<<endl;
   
  for(unsigned int itrk=0; itrk<trk_vector_05sig.size(); itrk++) h1_part_pt_05sig->Fill(trk_vector_05sig[itrk]->Pt());
  for(unsigned int itrk=0; itrk<trk_vector_1sig.size(); itrk++) h1_part_pt_1sig->Fill(trk_vector_1sig[itrk]->Pt());
  for(unsigned int itrk=0; itrk<trk_vector_15sig.size(); itrk++) h1_part_pt_15sig->Fill(trk_vector_15sig[itrk]->Pt());
  for(unsigned int itrk=0; itrk<trk_vector_2sig.size(); itrk++) h1_part_pt_2sig->Fill(trk_vector_2sig[itrk]->Pt());
  for(unsigned int itrk=0; itrk<trk_vector_3sig.size(); itrk++) h1_part_pt_3sig->Fill(trk_vector_3sig[itrk]->Pt());

  //***********************************************
  //*   Make inclusive photon-h foreground pairs  *
  //***********************************************

  //inclusive photon-hadron delta phi dists.
  MakePairs(clus_vector,trk_vector,REAL,data_set,h3_dphi[cbin],h3_dphi_fold[cbin],h3_ptxidphi[cbin],h3_ptxidphi_fold[cbin],h3_ptztdphi[cbin],h3_ptztdphi_fold[cbin]);
  if( DiagFlag ) MakePairs(clus_vector,trk_vector,REAL,data_set,0,h3_dphi[cbin],h3_dphi_fold[cbin],h3_ptxidphi[cbin],h3_ptxidphi_fold[cbin],h3_ptztdphi[cbin],h3_ptztdphi_fold[cbin],vector<TH2F*>(),vector<TH2F*>(),vector<TH2F*>(),vector<TH2F*>(),vector<TH2F*>(),vector<TH2F*>(),NULL,NULL,h2_dphi,h2_dphi_xi,h2_dphi_accw,h2_dphi_accw_xi,h3_dphi_accw,h3_dphi_accw_xi,h2_partpt_xi);
  MakePairs(clus_vector,trk_vector,REAL,data_set,1,h3_dphi_iso[cbin],h3_dphi_iso_fold[cbin],h3_ptxidphi_iso[cbin],h3_ptxidphi_iso_fold[cbin],NULL,NULL);
    
  //****************************************
  //*   Make pi0-h foreground pairs        *
  //****************************************  
  MakePairs(pi0_vector,trk_vector,REALPI,data_set,h3_dphi_pi0[cbin],h3_dphi_pi0_fold[cbin],h3_ptxidphi_pi0[cbin],h3_ptxidphi_pi0_fold[cbin],h3_ptztdphi_pi0[cbin],h3_ptztdphi_pi0_fold[cbin],h2_dphi_dec[cbin],h2_dphi_dec_fold[cbin],h2_dphixi_dec[cbin],h2_dphixi_dec_fold[cbin],h2_dphizt_dec[cbin],h2_dphizt_dec_fold[cbin]);
  if( DiagFlag ) MakePairs(pi0_vector,trk_vector,REALPI,data_set,0,h3_dphi_pi0[cbin],h3_dphi_pi0_fold[cbin],h3_ptxidphi_pi0[cbin],h3_ptxidphi_pi0_fold[cbin],h3_ptztdphi_pi0[cbin],h3_ptztdphi_pi0_fold[cbin],h2_dphi_dec[cbin],h2_dphi_dec_fold[cbin],h2_dphixi_dec[cbin],h2_dphixi_dec_fold[cbin],h2_dphizt_dec[cbin],h2_dphizt_dec_fold[cbin],NULL,NULL,h2_dphi_pi0,h2_dphi_xi_pi0,h2_dphi_accw_pi0,h2_dphi_accw_xi_pi0,h3_dphi_accw_pi0,h3_dphi_accw_xi_pi0,h2_partpt_xi_pi0);
  MakePairs(pi0_vector,trk_vector,REALPI,data_set,1,h3_dphi_pi0_iso[cbin],h3_dphi_pi0_iso_fold[cbin],h3_ptxidphi_pi0_iso[cbin],h3_ptxidphi_pi0_iso_fold[cbin],NULL,NULL,h2_dphi_dec_iso[cbin],h2_dphi_dec_iso_fold[cbin],h2_dphixi_dec_iso[cbin],h2_dphixi_dec_iso_fold[cbin],NULL,NULL);

  atree->SetEventData(evt,event_z,event_c,(int)clus_vector.size(),(int)pi0_vector.size(),(int)trk_vector.size());
  if( data_set == Run8dAu ) {
    AddMBEvent(data_set);
  }

  if(clus_vector.size() > 0 || pi0_vector.size() > 0) atree->_ttrig->Fill();
  atree->_tpart->Fill();

  Clear();
  event++;
  return EVENT_OK;
}

void Correlation::AddMBEvent(DataSet data_set)
{
  Fun4AllServer* se = Fun4AllServer::instance();
  PHCompositeNode* pwg_mixNode = se->topNode("MIXPWG");
  PHCompositeNode* cnt_mixNode = se->topNode("MIXCNT");
  float bbcz_mix;
  PHGlobal* mix_global = findNode::getClass<PHGlobal>(pwg_mixNode,"PHGlobal");
  if(!mix_global)
  {
    cout << PHWHERE << "findNode::GetClass: PHGlobal not in Node Tree" << endl;
    return;
  }
  bbcz_mix = mix_global->getBbcZVertex();
  // cout<<"bbczvertex = "<<z<<endl;
  float c = mix_global->getCentrality();
  if (c == -9999) {
    if(verbosity > 1) cout<<"Invalid Mixing Centrality!"<<endl;
    return;
  }
  if (bbcz_mix < -1*zVertexCut || bbcz_mix > zVertexCut) {
    if(verbosity > 1) cout<<"Outside vertex range: z = "<<event_z<<endl;
    return;
  }
  AEvent aevent;
  MakeEventObject(mix_global,&aevent);
  PHCentralTrack* mix_tracks = findNode::getClass<PHCentralTrack>(cnt_mixNode, "PHCentralTrack");
  if( !mix_tracks ) {
    cout << PHWHERE << "PHCentralTrack not in mixing Node Tree!" << endl;
    return;
  }
  emcClusterContainer* mix_clusters = findNode::getClass<emcClusterContainer>(pwg_mixNode, "emcClusterContainer");
  if(!mix_clusters)
  {
    cout << PHWHERE << "emcclustercontainter not in mixing Node Tree!" << endl;
    return;
  }

  int nclust = 0;
  int nclus = mix_clusters->size();
  for(int iclus = 0; iclus < nclus; iclus++){
    emcClusterContent* clus = mix_clusters->getCluster(iclus);
    ACluster acluster;
    MakeClusterObject(clus, &acluster);
    double ecore = acluster.GetEcore();
    if( ecore < 0.1 ) continue;
    if( !IsGoodTower(&acluster) ) continue;
    if( !IsGoodCluster(&acluster,data_set) ) continue;
    
    atree->SetPhotonData(acluster.Pt(),acluster.Phi(),acluster.Eta(),acluster.E(),nclust);
    nclust++;
  }

  int ntrack = 0;
  int npart = mix_tracks->get_npart();
  for(int ipart = 0; ipart < npart; ipart++){
    ATrack atrack;
    MakeTrackObject(mix_tracks, ipart, &atrack);
    double trk_pt = atrack.Pt();
    if(trk_pt < hadron_pt_min || trk_pt > hadron_pt_max) continue;
    if((atrack.GetZed() < -75.0) || (atrack.GetZed() > 75.0)) continue;
    if((atrack.GetQuality()!=31) && (atrack.GetQuality()!=63)) continue;
    if(trk_pt < 5.0 && atrack.GetN0() > 0) continue;
    if(!IsGoodTrack(&atrack,data_set)) continue;
    atree->SetPartnerData(atrack.Pt(),atrack.Phi(),atrack.Eta(),atrack.E(),atrack.GetPemcx(),atrack.GetPemcy(),atrack.GetPemcz(),ntrack);
    ntrack++;
  }
  atree->SetMBEventData(evt,aevent.GetVertex(),aevent.GetCentrality(),nclust,ntrack);
}

void Correlation::MakePi0s(vector<ACluster*> all_clusters, int cent, float zvertex, DataSet dataset)
{
  //  for( unsigned int iclus = 0; iclus < clusters.size()-1; iclus++ ){
  for( unsigned int iclus = 0; iclus < all_clusters.size(); iclus++ ){
    emcClusterContent* photon1 = all_clusters[iclus]->GetemcClusterContent(); 
    float ecore1 = all_clusters[iclus]->GetEcore();

    if( !IsGoodCluster(all_clusters[iclus],dataset) ) continue;
    if(ecore1 < photon_ecore_min) continue;
    if (!Chi2Cut(photon1,zvertex)) continue;

    // For ERT data leading photon must have fired the trigger! (bool set to true by default for MB data)
    if( !all_clusters[iclus]->IsTrigger() ) continue;
    
    float mindist = 99.;
    float ptofvetotrack = -1.0;
    VetoTracks(all_clusters[iclus],lessqualtrk_vector,mindist,ptofvetotrack);
    
    if(mindist<8.0) continue;
    
    //for( unsigned int jclus = iclus+1; jclus < clusters.size(); jclus++ ){
    for( unsigned int jclus = 0; jclus < all_clusters.size(); jclus++ ){
      if(jclus == iclus) continue;     
      float ecore2 = all_clusters[jclus]->GetEcore();
      if(ecore2<1.0 || ecore1<ecore2 || ecore1==ecore2) continue;
      
      float sume12 = ecore1 + ecore2;
      if(sume12 < 4.0) continue;
            
      //use asymmetry cut for low pT pi0s. only in 0-40% centralities
      float asym12 = (ecore1 - ecore2)/sume12;      
      if(cbin<2 && sume12<5.25 && asym12 > (minAsym + (1.0-minAsym)*((sume12-4.0)*(sume12-4.0)/1.25/1.25))) continue;
      
      int as = all_clusters[jclus]->GetArmSect();
      // if( as != all_clusters[iclus]->GetArmSect()) continue;
      int arm = all_clusters[jclus]->GetArm();
      if( arm != all_clusters[iclus]->GetArm()) continue;
      
      if( !IsGoodCluster(all_clusters[jclus],dataset) ) continue;

      APiZero api0(all_clusters[iclus], all_clusters[jclus]);
      api0.SetCent(event_c);
      api0.SetZvtx(event_z);

      
      h2_pi0mass[cbin]->Fill(api0.Pt(), api0.M());
      h2_pi0mass_as[as]->Fill(api0.Pt(), api0.M());
      if (as==4 || as==5) h2_pi0mass_PbGl->Fill(api0.Pt(), api0.M());
      else h2_pi0mass_PbSc->Fill(api0.Pt(), api0.M());
      
      if(api0.M() < 0.12 || api0.M() > 0.16) continue;
      all_clusters[iclus]->SetTag(true);

      if(api0.Pt() < pi0_pt_min || api0.Pt() > pi0_pt_max) continue;

      SetIso(&api0,lessqualtrk_vector,all_clusters,Rcut,econe_min[cbin]);
      if( verbosity > 1 )
        cout << "Event " << evt << " - Found pi0 with isolation = " << api0.IsIso() << endl;
      
      pi0_vector.push_back(api0.clone());

      //cout<<"setting pi0 trigger data: pt = "<<api0.Pt()<<"; phi = "<<api0.Phi()<<"; eta = "<<api0.Eta()<<"; e = "<<api0.E()<<"; x = "<<((ACluster*)api0.Daughter1())->GetX()<<"; y = "<<((ACluster*)api0.Daughter1())->GetY()<<"; z = "<<((ACluster*)api0.Daughter1())->GetZ()<<endl;
      atree->SetTriggerData(api0.Pt(),api0.Phi(),api0.Eta(),api0.E(),((ACluster*)api0.Daughter1())->GetX(),((ACluster*)api0.Daughter1())->GetY(),((ACluster*)api0.Daughter1())->GetZ(),api0.IsIso(),pi0_vector.size()-1);

      if( verbosity > 1 ) cout << "Event " << evt << " - Added pi0 with isolation = " << pi0_vector.at(pi0_vector.size()-1)->IsIso() << endl;
    }
  }
  for( unsigned int i = 0; i < pi0_vector.size(); i++)
  {
    h1_trig_pt_pi0[cbin]->Fill(pi0_vector[i]->Pt());

    if( pi0_vector[i]->IsIso() )
      h1_trig_pt_pi0_iso[cbin]->Fill(pi0_vector[i]->Pt());
    h1_trig_pt_pi0_tot->Fill(pi0_vector[i]->Pt());
    h2_ptvscent_trig_pi0->Fill(pi0_vector[i]->Pt(),cent);

    //dec weighting
    vector<float> mwweight;
    for(int n=0; n<5; n++) mwweight.push_back(0.0);
    EvalDecWeights(pi0_vector[i],event_z,cbin,mwweight);
    pi0_vector[i]->SetDecayWeights(mwweight);

    //dec trigger counting
    for(int ipw=0; ipw<5; ipw++){
      h1_trig_pt_dec[cbin]->Fill(ipw,mwweight[ipw]);
      if( pi0_vector[i]->IsIso() )
        h1_trig_pt_dec_iso[cbin]->Fill(ipw,mwweight[ipw]);
      else h1_trig_pt_dec_niso[cbin]->Fill(ipw,mwweight[ipw]);
      h1_trig_pt_dec_tot->Fill(ipw,mwweight[ipw]);
    }
  }
}

void Correlation::GetAcceptanceWeights(string filename)//input file is a previous taxi output, use h3_dphi_mix_c* to do acc. corr.
{
  cout << "In GetAcceptanceWeights" <<endl;
  float trig_pt_range[NTRIGBINS+1] = {5.0,7.0,9.0,12.0,15.0};
  float part_pt_range[NPARTBINS+1] = {0.5,1.0,2.0,3.0,5.0,7.0};
  float xi_range[NXIBINS+1] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};
  ostringstream bin;
  string name;
  
  TH1::AddDirectory(kFALSE);
  TFile* fin = TFile::Open(filename.c_str(), "READ");
  TH1F* temp;
  float norm = 2*PI;

  cout << "acceptance weight input file is opened." <<endl;
  for( int ic = 0; ic < 4; ic++ ){
    bin.str(""); bin << "h3_dphi_mix_c" << ic;
    TH3F* bginc = (TH3F*)fin->Get(bin.str().c_str());
    bin.str(""); bin << "h3_dphi_iso_mix_c" << ic;
    TH3F* bginc_iso = (TH3F*)fin->Get(bin.str().c_str());
    
    bin.str(""); bin << "h3_dphi_pi0_mix_c" << ic;
    TH3F* bgpi0 = (TH3F*)fin->Get(bin.str().c_str());
    bin.str(""); bin << "h3_dphi_pi0_iso_mix_c" << ic;
    TH3F* bgpi0_iso = (TH3F*)fin->Get(bin.str().c_str());
    
    bin.str(""); bin << "h3_ptxidphi_mix_c" << ic;
    TH3F* bginc_xi = (TH3F*)fin->Get(bin.str().c_str());
    bin.str(""); bin << "h3_ptxidphi_iso_mix_c" << ic;
    TH3F* bginc_xi_iso = (TH3F*)fin->Get(bin.str().c_str());
    
    bin.str(""); bin << "h3_ptxidphi_pi0_mix_c" << ic;
    TH3F* bgpi0_xi = (TH3F*)fin->Get(bin.str().c_str());
    bin.str(""); bin << "h3_ptxidphi_pi0_iso_mix_c" << ic;
    TH3F* bgpi0_xi_iso = (TH3F*)fin->Get(bin.str().c_str());

    for( int it = 0; it < NTRIGBINS; it++ ){
      int tbin = it;
      if( it>=3) tbin = it+1;
      bin.str(""); bin << tbin << "_c" << ic;
      name = "h2_dphi_dec_mix_p" + bin.str();
      TH2F* bgdec = (TH2F*)fin->Get(name.c_str());
      name = "h2_dphi_dec_iso_mix_p" + bin.str();
      TH2F* bgdec_iso = (TH2F*)fin->Get(name.c_str());

      name = "h2_dphixi_dec_mix_p" + bin.str();
      TH2F* bgdec_xi = (TH2F*)fin->Get(name.c_str());
      name = "h2_dphixi_dec_iso_mix_p" + bin.str();
      TH2F* bgdec_xi_iso = (TH2F*)fin->Get(name.c_str());
      
      for( int ip = 0; ip < NPARTBINS; ip++ ){
        bin.str("");
        bin << ic << "_p" << it << "_h" << ip;

        name = "h1_IncAcc_c" + bin.str();
        temp = MakeDphiProjection(bginc,trig_pt_range[it],trig_pt_range[it+1],part_pt_range[ip],part_pt_range[ip+1],norm);        
        IncAcc[ic][it][ip] = new TH1F(*temp);
        IncAcc[ic][it][ip]->SetName(name.c_str());
        name = "h1_IncAccIso_c" + bin.str();
        temp = MakeDphiProjection(bginc_iso,trig_pt_range[it],trig_pt_range[it+1],part_pt_range[ip],part_pt_range[ip+1],norm);        
        IncAccIso[ic][it][ip] = new TH1F(*temp);
        IncAccIso[ic][it][ip]->SetName(name.c_str());

        name = "h1_Pi0Acc_c" + bin.str();
        temp = MakeDphiProjection(bgpi0,trig_pt_range[it],trig_pt_range[it+1],part_pt_range[ip],part_pt_range[ip+1],norm);
        Pi0Acc[ic][it][ip] = new TH1F(*temp);
        Pi0Acc[ic][it][ip]->SetName(name.c_str());
        name = "h1_Pi0AccIso_c" + bin.str();
        temp = MakeDphiProjection(bgpi0_iso,trig_pt_range[it],trig_pt_range[it+1],part_pt_range[ip],part_pt_range[ip+1],norm);
        Pi0AccIso[ic][it][ip] = new TH1F(*temp);
        Pi0AccIso[ic][it][ip]->SetName(name.c_str());
        
        name = "h1_DecAcc_c" + bin.str();
        temp = MakeDphiProjection2D(bgdec,part_pt_range[ip],part_pt_range[ip+1],norm);
        DecAcc[ic][it][ip] = new TH1F(*temp);
        DecAcc[ic][it][ip]->SetName(name.c_str());
        name = "h1_DecAccIso_c" + bin.str();
        temp = MakeDphiProjection2D(bgdec_iso,part_pt_range[ip],part_pt_range[ip+1],norm);
        DecAccIso[ic][it][ip] = new TH1F(*temp);
        DecAccIso[ic][it][ip]->SetName(name.c_str());
      }

      for( int ix = 0; ix < NXIBINS; ix++ ){
        bin.str("");
        bin << ic << "_p" << it << "_x" << ix;
        name = "h1_IncAccXi_c" + bin.str();
        temp = MakeDphiProjection(bginc_xi,trig_pt_range[it],trig_pt_range[it+1],xi_range[ix],xi_range[ix+1],norm);
        IncAccXi[ic][it][ix] = new TH1F(*temp);
        IncAccXi[ic][it][ix]->SetName(name.c_str());
        name = "h1_IncAccXiIso_c" + bin.str();
        TH1F* temp = MakeDphiProjection(bginc_xi_iso,trig_pt_range[it],trig_pt_range[it+1],xi_range[ix],xi_range[ix+1],norm);
        IncAccXiIso[ic][it][ix] = new TH1F(*temp);
        IncAccXiIso[ic][it][ix]->SetName(name.c_str());
        
        name = "h1_Pi0AccXi_c" + bin.str();
        temp = MakeDphiProjection(bgpi0_xi,trig_pt_range[it],trig_pt_range[it+1],xi_range[ix],xi_range[ix+1],norm);        
        Pi0AccXi[ic][it][ix] = new TH1F(*temp);
        Pi0AccXi[ic][it][ix]->SetName(name.c_str());
        name = "h1_Pi0AccXiIso_c" + bin.str();
        temp = MakeDphiProjection(bgpi0_xi_iso,trig_pt_range[it],trig_pt_range[it+1],xi_range[ix],xi_range[ix+1],norm);        
        Pi0AccXiIso[ic][it][ix] = new TH1F(*temp);
        Pi0AccXiIso[ic][it][ix]->SetName(name.c_str());

        name = "h1_DecAccXi_c" + bin.str();
        temp = MakeDphiProjection2D(bgdec_xi,xi_range[ix],xi_range[ix+1],norm);
        DecAccXi[ic][it][ix] = new TH1F(*temp);
        DecAccXi[ic][it][ix]->SetName(name.c_str());
        name = "h1_DecAccXiIso_c" + bin.str();
        temp = MakeDphiProjection2D(bgdec_xi_iso,xi_range[ix],xi_range[ix+1],norm);
        DecAccXiIso[ic][it][ix] = new TH1F(*temp);
        DecAccXiIso[ic][it][ix]->SetName(name.c_str());
      }

      delete bgdec; delete bgdec_iso; delete bgdec_xi; delete bgdec_xi_iso;
    }
    delete bginc; delete bginc_xi; delete bginc_iso; delete bginc_xi_iso;
    delete bgpi0; delete bgpi0_xi; delete bgpi0_iso; delete bgpi0_xi_iso;
  }

  delete temp;
  fin->Close();
  delete fin;
}

double Correlation::GetAcc(TH1F* hist, float dphi)
{
  double acc = 0;
  int phibin = hist->FindBin(dphi);
  if(phibin < 1 || phibin > 60) acc = 0.;
  else acc = hist->GetBinContent(phibin);
  return acc;  
}

double Correlation::GetAcceptance(PairType type, int cbin, int tbin, int pbin, float dphi, int isxi, int isiso)
{
  //cout<<"In GetAcceptance..."<<endl;
  double acc = 1.;
  if(tbin < 0 || pbin < 0) return 0.;
  
  if(type == REAL || type == MIX) {
    if( isxi ) {
      if( isiso )
        acc = GetAcc(IncAccXiIso[cbin][tbin][pbin],dphi);
      else
        acc = GetAcc(IncAccXi[cbin][tbin][pbin],dphi);
    }
    else {
      if( isiso )
        acc = GetAcc(IncAccIso[cbin][tbin][pbin],dphi);
      else
        acc = GetAcc(IncAcc[cbin][tbin][pbin],dphi);
    }
  }

  if(type == REALPI || type == MIXPI) {
    if( isxi ) {
      if( isiso )
        acc = GetAcc(Pi0AccXiIso[cbin][tbin][pbin],dphi);
      else
        acc = GetAcc(Pi0AccXi[cbin][tbin][pbin],dphi);
    }
    else {
      if( isiso )
        acc = GetAcc(Pi0AccIso[cbin][tbin][pbin],dphi);
      else
        acc = GetAcc(Pi0Acc[cbin][tbin][pbin],dphi);
    }
  }

  if(type == DEC || type == MIXDEC) {
    if( isxi ) {
      if( isiso )
        acc = GetAcc(DecAccXiIso[cbin][tbin][pbin],dphi);
      else
        acc = GetAcc(DecAccXi[cbin][tbin][pbin],dphi);
    }
    else {
      if( isiso )
        acc = GetAcc(DecAccIso[cbin][tbin][pbin],dphi);
      else
        acc = GetAcc(DecAcc[cbin][tbin][pbin],dphi);
    }
  }
  
  return acc;
}

float Correlation::GetFlowWeights(PairType type, int tbin, int pbin, float dphifold)
{
  if(tbin<0 || pbin<0) return 0.;
  
  int typebin = -1;
  if( type == MIX ) typebin = 0;
  if( type == MIXPI ) typebin = 1;
  if( type == MIXDEC) typebin = 2;
  if( typebin < 0 ) return 1.0;  // only return flow modulation for mixed pairs
  
  float flowweight = (1+trig_v2[typebin][cbin][tbin]*part_v2[cbin][pbin]*cos(2*dphifold));
  //cout<<"flowweight = "<< flowweight << endl;
  return flowweight;
}

TH1F* Corelation::MakeDphiProjection2D(TH2F* h2, float ymin, float ymax, float norm) 
{
  int min = bgdec->GetYaxis()->FindBin(ymin);
  int max = bgdec->GetYaxis()->FindBin(ymax);
  TH1F* proj_hist = (TH1F*)bgdec->ProjectionX("proj",min,max-1);
  proj_hist->Scale(norm/proj_hist->Integral("width"));

  return proj_hist;
}

TH1F* Correlation::MakeDphiProjection(TH3F* h3, float xmin, float xmax, float ymin, float ymax, float norm)
{
  TH1F* proj_x = (TH1F*)h3->ProjectionX("px");
  TH1F* proj_y = (TH1F*)h3->ProjectionY("py");
  int xbinlo = proj_x->FindBin(xmin);
  int xbinhi = proj_x->FindBin(xmax);
  int ybinlo = proj_y->FindBin(ymin);
  int ybinhi = proj_y->FindBin(ymax);

  TH1F* proj_hist = (TH1F*)(h3->ProjectionZ("proj",xbinlo,xbinhi-1,ybinlo,ybinhi-1));
  proj_hist->Scale(norm/proj_hist->Integral("width"));

  delete proj_x;
  delete proj_y;
  
  return proj_hist;
}

void Correlation::MakeAccHistos(TH1F* h1in, TH1F* h1out)//take the projected 1d histo and feed the content into IncAcc, Pi0Acc
{
  int xbins = h1in->GetNbinsX();
  int xbinsout = h1out->GetNbinsX();
  if(xbins != xbinsout) cout<<"trouble!"<<endl;
  
  for (int ibin=0; ibin<=xbins+1; ibin++){
    double bincont = h1in->GetBinContent(ibin);
    double binerr = h1in->GetBinError(ibin);
    h1out->SetBinContent(ibin,bincont);
    h1out->SetBinError(ibin,binerr);
  }
}

void Correlation::MakeEventObject(PHGlobal* glob, AEvent* aevt)
{
  aevt->SetVertex(NZVTX, 60.0, glob->getBbcZVertex());
  aevt->SetCentrality(NCENT, 100.0, glob->getCentrality());
}

void Correlation::MakeClusterObject(emcClusterContent* clus, ACluster* aclus)
{
  aclus->SetCent(event_c);
  aclus->SetZvtx(event_z);
  aclus->SetemcClusterContent(clus);
  aclus->SetWarnDead(clus->warnmap(),clus->deadmap());
  aclus->SetEcore(clus->ecore());
  aclus->SetTheta(clus->theta());
  aclus->SetPhi(clus->phi());
  aclus->SetProb(clus->prob_photon());
  aclus->SetArm(clus->arm());
  aclus->SetSec(clus->sector());
  aclus->SetIypos(clus->iypos());
  aclus->SetIzpos(clus->izpos());
  aclus->SetEmcX(clus->x());
  aclus->SetEmcY(clus->y());
  aclus->SetEmcZ(clus->z());
  aclus->SetEmctrkdphi(clus->emctrkdphi());
  aclus->SetEmctrkdz(clus->emctrkdz());
  aclus->SetEmcpc3dphi(clus->emcpc3dphi());
  aclus->SetEmcpc3dz(clus->emcpc3dz());
  aclus->SetPtEtaPhiE(clus->ecore()*sin(clus->theta()),
                      -1*log(tan(clus->theta()/2)),
                      clus->phi(),clus->ecore());
  aclus->SetFiducial(CheckPhiFiducial(clus->phi()));
  
  if(verbosity > 3){
    cout<<"clus ecore = "<<aclus->GetEcore()<<endl;
    cout<<"clus theta = "<<aclus->Theta()<<endl;
    cout<<"clus phi = "<<aclus->Phi()<<endl;
    cout<<"clus prob = "<<aclus->GetProb()<<endl;
    cout<<"clus sector = "<<aclus->GetArmSect()<<endl;
    cout<<"clus ypos = "<<aclus->GetIypos()<<endl;
    cout<<"clus zpos = "<<aclus->GetIzpos()<<endl;
    cout<<"clus pt = "<<aclus->Pt()<<endl;
    cout<<"clus emcpc3dphi = "<<aclus->GetEmcpc3dphi()<<endl;
    cout<<"clus emcpc3dz = "<<aclus->GetEmcpc3dz()<<endl;
    cout<<"clus emctrkdphi = "<<aclus->GetEmctrkdphi()<<endl;
    cout<<"clus emctrkdz = "<<aclus->GetEmctrkdz()<<endl;
  }
}

void Correlation::MakeTrackObject(PHCentralTrack* trk, int itrk, ATrack* atrk)
{
  atrk->SetZvtx(event_z);
  if(RecalFlag) {
    atrk->SetPc3Match(trk->get_pc3sdphi(itrk), trk->get_pc3sdz(itrk));
    atrk->SetEmcMatch(trk->get_emcsdphi(itrk), trk->get_emcsdz(itrk));
  }
  else {
    atrk->SetPc3Match(trk->get_pc3dphi(itrk), trk->get_pc3dz(itrk));
    atrk->SetEmcMatch(trk->get_emcdphi(itrk), trk->get_emcdz (itrk));
  }
  atrk->SetPHCentralTrack(trk);
  atrk->SetPxPyPz(trk->get_px(itrk), trk->get_py(itrk), trk->get_pz(itrk));
  atrk->SetPc3dz(trk->get_pc3dz(itrk));
  atrk->SetPc3dphi(trk->get_pc3dphi(itrk));
  atrk->SetPpc3(trk->get_ppc3x(itrk), trk->get_ppc3y(itrk), trk->get_ppc3z(itrk));
  atrk->SetPpc1(trk->get_ppc1x(itrk), trk->get_ppc1y(itrk), trk->get_ppc1z(itrk));
  atrk->SetPemc(trk->get_pemcx(itrk), trk->get_pemcy(itrk), trk->get_pemcz(itrk));
  atrk->SetN0(trk->get_n0(itrk));
  atrk->SetQuality(trk->get_quality(itrk));
  atrk->SetAlpha(trk->get_alpha(itrk));
  atrk->SetCharge(trk->get_charge(itrk));
  atrk->SetTheta(trk->get_the0(itrk));
  atrk->SetPhiD(trk->get_phi(itrk));
  atrk->SetPhi(trk->get_phi0(itrk));
  atrk->SetZed(trk->get_zed(itrk));
  atrk->SetDCArm(trk->get_dcarm(itrk));
  atrk->SetBoard(trk->get_phi(itrk), trk->get_dcarm(itrk));
  atrk->SetPtEtaPhiE(sqrt(trk->get_px(itrk)*trk->get_px(itrk) + trk->get_py(itrk)*trk->get_py(itrk)),-1*log(tan(trk->get_the0(itrk)/2)),trk->get_phi0(itrk),trk->get_ecore(itrk));
  atrk->SetEcore(trk->get_ecore(itrk));
  if(verbosity > 3){
    cout<<"trk ppc3x = "<<atrk->GetPpc3x()<<endl;
    cout<<"trk ppc3y = "<<atrk->GetPpc3y()<<endl;
    cout<<"trk ppc3z = "<<atrk->GetPpc3z()<<endl;
    cout<<"trk pc3sdz = "<<atrk->GetPc3sdz()<<endl;
    cout<<"trk pc3sdphi = "<<atrk->GetPc3sdphi()<<endl;
    cout<<"trk n0 = "<<atrk->GetN0()<<endl;
    cout<<"trk qual = "<<atrk->GetQuality()<<endl;
    cout<<"trk alpha = "<<atrk->GetAlpha()<<endl;
    cout<<"trk phi = "<<atrk->Phi()<<endl;
    cout<<"trk phid = "<<atrk->GetPhiD()<<endl;
    cout<<"trk zed = "<<atrk->GetZed()<<endl;
    cout<<"trk pt = "<<atrk->Pt()<<endl;
  }
}

bool Correlation::CheckPhiFiducial(float phi)
{
  float phphi = PHAngle(phi);
  if((phphi>(WEST_LOW_EDGE+Rcut*0.25)&&phphi<(WEST_HIGH_EDGE-Rcut*0.25))||(phphi>(EAST_LOW_EDGE+Rcut*0.25)&&phphi<(EAST_HIGH_EDGE-Rcut*0.25)))
    return true;
  else
    return false;
}

void Correlation::MakeClusterObject(ACluster* aclus, float pt, float phi, float eta, float e, float x, float y, float z, float zvtx)
{
  aclus->SetPtEtaPhiE(pt, eta, phi, e);
  aclus->SetEcore(e);
  aclus->SetEmcX(x);
  aclus->SetEmcY(y);
  aclus->SetEmcZ(z);
  aclus->SetZvtx(zvtx);
  aclus->SetFiducial(CheckPhiFiducial(aclus->Phi()));
}

void Correlation::MakePi0Object(APiZero* api0, float pt, float phi, float eta, float e, float x, float y, float z, float zvtx)
{
  api0->SetPtEtaPhiE(pt, eta, phi, e);
  api0->Daughter1()->SetPtEtaPhiE(pt, eta, phi, e);
  api0->Daughter2()->SetPxPyPzE(0,0,0,0);
  api0->Daughter1()->SetEmcX(x);
  api0->Daughter1()->SetEmcY(y);
  api0->Daughter1()->SetEmcZ(z);
  api0->Daughter1()->SetZvtx(zvtx);
}

void Correlation::MakeTrackObject(ATrack* atrk, float pt, float phi, float eta, float e, float pemcx, float pemcy, float pemcz, float zvtx)
{
  atrk->SetPtEtaPhiE(pt, eta, phi, e);
  atrk->SetPemc(pemcx, pemcy, pemcz);
  atrk->SetZvtx(zvtx);
}

// function Chi2Cut, CalcPbglDisp, IncidentAngle are directly adapted from CombineSimple code -031015
bool Correlation::Chi2Cut(emcClusterContent* sngl_emc, float zvertex, TH1F* h1_chi2, TH3F* h3_disp)
{
  int arm = sngl_emc->arm();
  int sector = sngl_emc->sector();
  int iypos = sngl_emc->iypos();
  int izpos = sngl_emc->izpos();
  
  int index = arm * 100000 + sector * 10000 + iypos * 100 + izpos;
  
  if (index < 100000 || index >=120000) // pbsc
    {
      if(h1_chi2) h1_chi2->Fill(sngl_emc->chi2());
      if(sngl_emc->chi2()<3.0) return true;
    }      
  else // pbgl...dispersion cut
    {
      
      float yz_cg[2], disp[2], pbgl_disp_z_y[2];
      
      yz_cg[0] = sngl_emc->zcg();
      yz_cg[1] = sngl_emc->ycg();
      disp[0] = sngl_emc->dispz(); 
      disp[1] = sngl_emc->dispy(); 
      
      CalcPbglDisp(yz_cg[0], // /cg_z
                   yz_cg[1], //cg_y,
                   disp[0],
                   disp[1],
                   pbgl_disp_z_y);
      
      float fDispYCor = pbgl_disp_z_y[1];
      float fDispZCor = pbgl_disp_z_y[0];
      
      double fX = sngl_emc->x();
      double fY = sngl_emc->y();
      double fZ = sngl_emc->z();
      double fBbcZEvt = zvertex;
      
      // taken straight from offline/ana/PbGl/gampi0 
      // ---> Flag 3 : dispersion cut
      // angle dependent dispersion cuts
      // the cut values are taken from PbGl and are not tuned for PbSc
      static const double p1 =  0.270;
      static const double p2 = -0.0145;
      static const double p3 =  0.00218;
      
      int arm = sngl_emc->arm();
      
      // Incident Angle in degrees (0 = perpendicular to sector surface)
      double theta = IncidentAngle(arm, (index % 100000 / 10000), fX,fY,fZ,fBbcZEvt);
      
      // angle dependent dispersion cut
      double dispCut = p1 + p2*theta + p3*theta*theta;
      
      // maximal dispersion
      double dispMax = fDispZCor;
      if (dispMax < fDispYCor) dispMax = fDispYCor;
      int pass = 1;
      if (dispMax >= dispCut) {
        pass = 0;
        if(h3_disp) h3_disp->Fill(dispMax,dispCut,pass);
      }
      else {// apply dispersion cut
        if(h3_disp) h3_disp->Fill(dispMax,dispCut,pass);
        return true;
      }
      // cout<< "dispMax "<<dispMax<< " dispCut "<<dispCut<<endl;
    } //if PbGl
     
  return false;  
}

void Correlation::CalcPbglDisp(const float& m1z, const float& m1y,
                                      const float& m2z, const float& m2y,
                                            float* returnVals_z_y
                                      ) 
{

  double zModSize, yModSize;

  // modules sizes extracted from CVS
  // offline/packages/emc/mEmcGeometryModule.C (revision 2.8)
  // -> mEmcGeometryModule::BuildGeometry()

  //   if ((iSecId == idSecE0) || (iSecId == idSecE1)) {
  zModSize = 4.09168;
  yModSize = 4.10615;
  //   }
  //   else {
  //     zModSize = 5.58003;
  //     yModSize = 5.57091;
  //   }

  //=====> dispersion stuff
  // correct simple periodic structure of the dispersion
 float fPosZ  = m1z / zModSize;
 float fPosY  = m1y / yModSize;
 
 float fDispZ = m2z / (zModSize * zModSize);
 float fDispY = m2y / (yModSize * yModSize);
 
  double zPosMod = fPosZ - floor(fPosZ);
  double yPosMod = fPosY - floor(fPosY);
  
  float fDispZCor = fDispZ - (zPosMod - zPosMod*zPosMod);
  float fDispYCor = fDispY - (yPosMod - yPosMod*yPosMod);

  returnVals_z_y[0] = fDispZCor;
  returnVals_z_y[1] = fDispYCor;
}

double Correlation::IncidentAngle(const int arm,
                                  const int sect,
                                  const double& fX, 
                                  const double& fY,
                                  const double& fZ,
                                  const double& fBbcZEvt)
{
  
  double CosPhiNum =   vnx[arm][sect] * fX 
    + vny[arm][sect] * fY 
    + vnz[arm][sect] * (fZ-fBbcZEvt); 
  
  double CosPhiDen = sqrt(fX*fX + fY*fY + (fZ-fBbcZEvt)*(fZ-fBbcZEvt));
  
  double Phi;
  if (CosPhiDen != 0) Phi = 180./3.14159 * acos(CosPhiNum / CosPhiDen);
  else                Phi = -999.;
  
  return Phi;  // in degrees 
}

void Correlation::Init1DHisto(TH1F* &h1, string name, string xtitle, int nxbin, double xmin, double xmax)
{
  manager->registerHisto(name.c_str(),h1 = new TH1F(name.c_str(), name.c_str(), nxbin, xmin, xmax));
  h1->Sumw2();
  h1->GetXaxis()->SetTitle(xtitle.c_str());
}

void Correlation::Init2DHisto(TH2F* &h2, string name, string xtitle, int nxbin, double xmin, double xmax, string ytitle, int nybin, double ymin, double ymax)
{
  manager->registerHisto(name.c_str(), h2 = new TH2F(name.c_str(), name.c_str(), nxbin, xmin, xmax, nybin, ymin, ymax));
  h2->Sumw2();
  h2->GetXaxis()->SetTitle(xtitle.c_str());
  h2->GetYaxis()->SetTitle(ytitle.c_str());
}

void Correlation::Init3DHisto(TH3F* &h3, string name, string xtitle, int nxbin, double xmin, double xmax, string ytitle, int nybin, double ymin, double ymax, string ztitle, int nzbin, double zmin, double zmax)
{
  manager->registerHisto(name.c_str(), h3 = new TH3F(name.c_str(), name.c_str(), nxbin, xmin, xmax, nybin, ymin, ymax, nzbin, zmin, zmax));
  h3->Sumw2();
  h3->GetXaxis()->SetTitle(xtitle.c_str());
  h3->GetYaxis()->SetTitle(ytitle.c_str());
  h3->GetZaxis()->SetTitle(ztitle.c_str());
}

void Correlation::Init3DHisto(TH3D*& h3, string name, string xtitle, int nxbins, double* xbins, string ytitle, int nybins, double* ybins, string ztitle, int nzbins, double* zbins)
{
  manager->registerHisto(name.c_str(), h3 = new TH3D(name.c_str(), name.c_str(), nxbins, xbins, nybins, ybins, nzbins, zbins));
  h3->GetXaxis()->SetTitle(xtitle.c_str());
  h3->GetYaxis()->SetTitle(ytitle.c_str());
  h3->GetZaxis()->SetTitle(ztitle.c_str());
}

void Correlation::InitPhotonCutChecker(string name, THmulf*& histo)
{
  manager->registerHisto(name.c_str(), histo = new THmulf(name.c_str(),name.c_str()));
  histo->AddAxis("pt", "p_{T}", 20, photon_pt_min, photon_pt_max);
  histo->AddAxis("prob", "photon_prob_cut", 2, 0, 2);
  // histo->AddAxis("emctrkdphi", "emctrkdphi", 30, -0.2, 0.2);
  // histo->AddAxis("emctrkdz", "emctrkdz", 35, -70., 70.);
  // histo->AddAxis("emcpc3dphi", "emcpc3dphi", 20, -0.2, 0.2);
  // histo->AddAxis("emcpc3dz", "emcpc3dz", 35, -70., 70.);
  histo->AddAxis("emctrk", "emctrk_cut", 2, 0, 2);
  histo->AddAxis("emcpc3", "emcpc3_cut", 2, 0, 2);
  histo->Sumw2();
}

void Correlation::InitHadronCutChecker(string name, THmulf*& histo)
{
  manager->registerHisto(name.c_str(), histo = new THmulf(name.c_str(),name.c_str()));
  histo->AddAxis("pt", "p_{T}", 10, hadron_pt_min, hadron_pt_max);
  histo->AddAxis("qual", "quality", 2, 0, 2);
  histo->AddAxis("n0", "n0", 2, 0, 2);
  histo->AddAxis("pc3sdphi", "pc3sdphi", 40, -10., 10.);
  histo->AddAxis("pc3sdz", "pc3sdz", 40, -10., 10.);
  histo->Sumw2();
}

int Correlation::VetoTracks(ACluster* aclus, vector<ATrack*> lessqualtrk_vec)
{
  float mindist = 99.0;
  for(unsigned int it=0; it<lessqualtrk_vec.size(); it++){
    float dist = FindTrackDistance(aclus, lessqualtrk_vec[it]);
    if(dist < mindist && lessqualtrk_vec[it]->Pt() > vetoPtCut) mindist = dist;
  }
  if (mindist < 8.0) return 1;
  
  return 0;
}

int Correlation::VetoTracks(ACluster* aclus, vector<ATrack*> lessqualtrk_vec, float& mindist, float& vetopt)
{
  for(unsigned int it=0; it<lessqualtrk_vec.size(); it++){
    float dist = FindTrackDistance(aclus, lessqualtrk_vec[it]);
    if(dist < mindist && lessqualtrk_vec[it]->Pt() > vetoPtCut) {
      mindist = dist;
      vetopt = lessqualtrk_vec[it]->Pt();
    }
  }
  if (mindist < 8.0) return 1;
  
  return 0;
}

void Correlation::FillClusterQAHistos(int isafter, float pt, float emctrkdphi, float emctrkdz, float emcpc3dphi, float emcpc3dz)
{
  if(!isafter){
    h2_emctrkdphi_bf->Fill(pt,emctrkdphi);
    h2_emctrkdz_bf->Fill(pt,emctrkdz);
    h2_emcpc3dphi_bf->Fill(pt,emcpc3dphi);
    h2_emcpc3dz_bf->Fill(pt,emcpc3dz);
  }
  else {
    h2_emctrkdphi_aft->Fill(pt,emctrkdphi);
    h2_emctrkdz_aft->Fill(pt,emctrkdz);
    h2_emcpc3dphi_aft->Fill(pt,emcpc3dphi);
    h2_emcpc3dz_aft->Fill(pt,emcpc3dz);
  }
}

void Correlation::FillTrackQAHistos(int isafter, float ppc3x, float ppc3z, float ppc3y, float ppc1z, float ppc1y, float zed, float phid, float pt, float quality, int n0, float pc3sdphi, float pc3sdz)
{
  if(!isafter){
    if(ppc3x > 0) {
      h2_ppc3_west_bf->Fill(ppc3z, ppc3y);
      h2_ppc1_west_bf->Fill(ppc1z, ppc1y);
    }
    if(ppc3x < 0) {
      h2_ppc3_east_bf->Fill(ppc3z, ppc3y);
      h2_ppc1_east_bf->Fill(ppc1z, ppc1y);
    }
    h2_phi_zed_bf->Fill(zed, phid);
    hadron_cut_check->Fill(1, pt, (quality==63||quality==31), n0<=0, pc3sdphi, pc3sdz);
  }
  else{
    if(ppc3x > 0) {
      h2_ppc3_west_aft->Fill(ppc3z, ppc3y);
      h2_ppc1_west_aft->Fill(ppc1z, ppc1y);
    }
    if(ppc3x < 0) {
      h2_ppc3_east_aft->Fill(ppc3z, ppc3y);
      h2_ppc1_east_aft->Fill(ppc1z, ppc1y);
    }
    h2_phi_zed_aft->Fill(zed, phid);
    if(pc3sdphi>-9000 && pc3sdz > -9000){
      h2_pc3sdphi_bf->Fill(pt, pc3sdphi);
      h2_pc3sdz_bf->Fill(pt, pc3sdz);
    }
  }
}

bool Correlation::IsGoodTower(ACluster* aclus)
{
  bool good_tower = true;

  double ecore = aclus->GetEcore();
  int armsect = aclus->GetArmSect();
  int ypos = aclus->GetIypos();
  int zpos = aclus->GetIzpos();

  if(warnmap->IsBad(armsect, ypos, zpos, ecore)){
    if(verbosity > 1) cout<<"Cluster falls in a bad tower! Sect: "<< armsect << ", ypos: "<< ypos << ", zpos: "<< zpos << endl;
    good_tower = false;
  }
  return good_tower;
}

bool Correlation::IsGoodCluster(ACluster* aclus, DataSet dataset)
{
  bool good_clus = true;
  
  //if(((aclus->GetWarnStat() & Cut3x3Map)!= 0) || (aclus->GetProb() < 0.02) || (aclus->GetEcore() < 0.04)) good_clus = false;
  if(dataset == Run10AuAu) Cut3x3Map = 0x7de1ce70;
  else Cut3x3Map = 0x1ce70;
  if(((aclus->GetWarnStat() & Cut3x3Map)!= 0) || ((aclus->GetDeadStat() & Cut3x3Map)!= 0) || (aclus->GetEcore() < 0.04)) good_clus = false;
  return good_clus;
}

bool Correlation::IsGoodTrack(ATrack* atrk, DataSet dataset)
{
  bool good_trk = true;  
  //check and cut on pc3 matching
  if(RecalFlag == 0){ 
    if (((atrk->GetPc3sdphi() < -0.014) || (atrk->GetPc3sdphi() > 0.015)) ||
        ((atrk->GetPc3sdz() < -5.0) || (atrk->GetPc3sdz() > 7.0))){
      if(verbosity > 3) cout<<"bad track! pc3dphi = "<<atrk->GetPc3sdphi()<<"; pc3dz = "<<atrk->GetPc3sdz()<<endl;
      good_trk = false;
    }
  }
  if(RecalFlag == 1){  
    if(dataset == Run11AuAu){
      if(!PassMatchingCut(atrk,PC3_NSIGMA,EMC_NSIGMA)) good_trk = false; 
    }
    else {
      if(!PassMatchingCut(atrk,PC3_NSIGMA)) good_trk = false;
      //if (!(sqrt(atrk->GetPc3sdz*atrk->GetPc3sdz+atrk->GetPc3sdphi*atrk->GetPc3sdphi)<2.0) ) continue;
    } 
  }
  if(dataset==Run8dAu){
    if( atrk->Phi() > 0.25 && atrk->Phi() < 0.4 ) good_trk = false;
  }
  return good_trk;
}

bool Correlation::PassMatchingCut(ATrack* atrk, float pc3nsig)
{
  bool pass = true;
  double circular_pc3cut = sqrt(atrk->GetPc3sdphi()*atrk->GetPc3sdphi()+atrk->GetPc3sdz()*atrk->GetPc3sdz());
  
  if((circular_pc3cut >= pc3nsig) || isnan(circular_pc3cut)) pass = false;
  //if(circular_pc3cut >= pc3nsig) pass = false;
  return pass;
}

bool Correlation::PassMatchingCut(ATrack* atrk, float pc3nsig, float emcnsig)
{
  bool pass = true;
  double circular_pc3cut = sqrt(atrk->GetPc3sdphi()*atrk->GetPc3sdphi()+atrk->GetPc3sdz()*atrk->GetPc3sdz());
  double circular_emccut = sqrt(atrk->GetEmcsdphi()*atrk->GetEmcsdphi()+atrk->GetEmcsdz()*atrk->GetEmcsdz());
  
  if(fabs(atrk->GetPc3dz())<9999 && fabs(atrk->GetPc3dphi())<9999){
    if((circular_pc3cut >= pc3nsig) || isnan(circular_pc3cut)) pass = false;
  }
  else{
    if ((circular_emccut >= emcnsig) || isnan(circular_emccut)) pass = false;
  }
  return pass;
}

int Correlation::GetPtBin(float pt, int istrig)
{
  int ipt = -1;
  if(istrig){
    if(pt>=5.0 && pt<7.0) ipt = 0;
    if(pt>=7.0 && pt<9.0) ipt = 1;
    if(pt>=9.0 && pt<12.0) ipt = 2;
    if(pt>=12.0 && pt<15.0) ipt = 3;
  }
  else {
    if(pt>=0.5 && pt<1.0) ipt = 0;
    if(pt>=1.0 && pt<2.0) ipt = 1;
    if(pt>=2.0 && pt<3.0) ipt = 2;
    if(pt>=3.0 && pt<5.0) ipt = 3;
    if(pt>=5.0 && pt<7.0) ipt = 4;
  }
  return ipt;
}

int Correlation::GetXiBin(float xi)
{
  int ixi = -1;
  if( xi < 0.4 ) ixi = 0;
  if( xi >= 0.4 && xi < 0.8 ) ixi = 1;
  if( xi >= 0.8 && xi < 1.2 ) ixi = 2;
  if( xi >= 1.2 && xi < 1.6 ) ixi = 3;
  if( xi >= 1.6 && xi < 2.0 ) ixi = 4;
  if( xi >= 2.0 ) ixi = 5;
  return ixi;
}

int Correlation::GetCentBin(int cent)
{
  int cbin = -1;
  if(cent > 0 && cent <= 20) cbin = 0;
  else if(cent > 20 && cent <= 40) cbin = 1;
  else if(cent > 40 && cent <= 60) cbin = 2;
  else if(cent > 60 ) cbin = 3;
  return cbin;
}

float Correlation::FindTrackDistance(ACluster* clus, ATrack* trk)
{
  float x1 = clus->GetX();
  float y1 = clus->GetY();
  float z1 = clus->GetZ();
  float zvtx1 = clus->GetZvtx();
  TVector3 vCluster(x1, y1, z1-zvtx1);
  
  float x2 = trk->GetPemcx();
  float y2 = trk->GetPemcy();
  float z2 = trk->GetPemcz();
  float zvtx2 = trk->GetZvtx();
  
  double trackdist = 99.;
  if(x2>-9999 && y2>-9999 && z2>-9999){
    TVector3 vTrack(x2, y2, z2-zvtx2);
    TVector3 v3 = vCluster - vTrack;
    trackdist = v3.Mag();
  }
  return trackdist;
}

void Correlation::EvalDecWeights(APiZero* pi0trigger, float zvertex, int cbin, vector<float>& mwweight)
{
  float pi0trigpt = pi0trigger->Pt();
  //  float pi0trigpz = pi0trigger->Pz();
  
  float pi0trigeff = 0.;
  if(cbin == 0) pi0trigeff = grpi0eff_0->Eval(pi0trigpt); 
  if(cbin == 1) pi0trigeff = grpi0eff_1->Eval(pi0trigpt);
  if(cbin == 2) pi0trigeff = grpi0eff_2->Eval(pi0trigpt);
  if(cbin == 3) pi0trigeff = grpi0eff_3->Eval(pi0trigpt);
  
  int trigptbin = hshark_large[0][0]->FindBin(pi0trigpt);
  if(trigptbin>400) trigptbin=400;
  
  int ipi0zemc = GetPi0ZEMCBin(pi0trigger);
  
  for(int idecl=0;idecl<5;idecl++){
    
    float mattshark=hshark_large[idecl][ipi0zemc]->GetBinContent(trigptbin);
    if(mattshark>0) {
      //cout << "mattshark = " << mattshark << "; pi0trigeff = " << pi0trigeff << endl;
      mwweight[idecl]=mattshark*pi0trigeff;
      //cout << "mwweight[" << idecl << "] = " << mwweight[idecl] << endl;
    }
  }
}

float Correlation::GetFilltimeWeight(PairType type, float dphi, float partpt, int pbin, int tbin, int isxi, int isiso)
{
  //cout << "GetFilltimeWeight: " << endl;
  if(tbin<0 || pbin<0) return 0.;

  float filltimeweight = 1.;
  float filltimeflow = 1.;

  float seffcorr = GetHadronEfficiencyCorr(partpt);
  if( verbosity > 1 ) cout << PHWHERE << "seffcorr = " << seffcorr << endl;
  
  float accw = 1.0;
  accw = GetAcceptance(type, cbin, tbin, pbin, dphi, isxi, isiso);  
  if( verbosity > 1 ) cout << PHWHERE << "accw at dphi = " << dphi << " for decay: " << accw << endl;
  if( accw > 0 ) filltimeweight = seffcorr/accw;
  else filltimeweight = 0.;
  if( verbosity > 1 ) cout << PHWHERE << "filltimeweight = " << filltimeweight << endl;

  // GetFlowWeights returns 1.0 if these are real pairs
  // Don't apply flow modulation for non Au+Au runs (like dAu)
  if( isxi ) pbin = GetPtBin(partpt, 0); //for xi case need to calculate partner pT bin for getting v2
  if( data_set != Run8dAu ) filltimeflow = GetFlowWeights(type,tbin,pbin,dphi)*filltimeweight;
  else filltimeflow = filltimeweight;
  if( verbosity > 1 ) cout << PHWHERE << "filltimeweight = " << filltimeflow << endl;

  return filltimeflow;
}

void Correlation::MakeDecays(PairType type, int isiso, float dphi, float dphifold, float partpt, float trigpt, APiZero* pizero, std::vector<float> weight,
                             std::vector<TH2F*> hdphi, std::vector<TH2F*> hdphi_fold, 
                             std::vector<TH2F*> hdphixi, std::vector<TH2F*> hdphixi_fold,
                             std::vector<TH2F*> hdphizt, std::vector<TH2F*> hdphizt_fold)                     
{
  if( verbosity > 1) cout<<"MakeDecays: type = "<<type<<endl;

  //calculate xi weights for filltime 
    
  float zt = partpt/trigpt;
  //  float xi = log(1.0/zt);

  int pbin = GetPtBin(partpt, 0);
  //int xbin = GetXiBin(xi);
  float filltimeflow = 1.0;
  float filltimeflowxi = 1.0;

  for(unsigned int ipw=0;ipw<hdphi.size();ipw++){
    int tbin = ipw;
    if(ipw == 3) continue;
    if(ipw > 3) tbin = ipw - 1;
    float seffcorr = GetHadronEfficiencyCorr(partpt);
    filltimeflow = GetFilltimeWeight(type,dphi,partpt,pbin,tbin,0,isiso);

    if(weight[ipw]>0) {
      //cout << "weight["<< ipw << "] = " << weight[ipw] << endl;
      hdphi[ipw]->Fill(dphi,partpt,weight[ipw]*seffcorr);
      if( hdphi_fold.size()>ipw ){
        hdphi_fold[ipw]->Fill(dphifold,partpt,weight[ipw]*filltimeflow);
      }
      if( hdphizt.size()>ipw ){         
        hdphizt[ipw]->Fill(dphi,zt,weight[ipw]*filltimeflow);           
      }         
      if( hdphizt_fold.size()>ipw ){            
        hdphizt_fold[ipw]->Fill(dphifold,zt,weight[ipw]*filltimeflow);          
      }

      //getting weights for xi-binned histograms
      int nXibins = hdphixi[ipw]->GetNbinsY();
      float sumtrigweight = 0.0;
      int ipi0zemc = GetPi0ZEMCBin(pizero);
      float decxibw = hdphixi[ipw]->GetYaxis()->GetBinWidth(1);
      for(int ixidecbin=0; ixidecbin<nXibins; ixidecbin++){
        float xidec = hdphixi[ipw]->GetYaxis()->GetBinCenter(ixidecbin+1);
        sumtrigweight += GetDecayXiWeights(decxibw,tbin,xidec,ipi0zemc,trigpt,partpt);
      }

      float fineweightave = 0.0;
      for(int ixidecbin=0; ixidecbin<nXibins; ixidecbin++){
        float xidec = hdphixi[ipw]->GetYaxis()->GetBinCenter(ixidecbin+1);
        fineweightave = GetDecayXiWeights(decxibw,tbin,xidec,ipi0zemc,trigpt,partpt);
        if(sumtrigweight>0) fineweightave *= weight[ipw]/sumtrigweight;

        int ixibin = GetXiBin(xidec);
        filltimeflowxi = GetFilltimeWeight(type,dphi,partpt,ixibin,tbin,1,isiso);
        if (hdphixi.size()>ipw) hdphixi[ipw]->Fill(dphi,xidec,fineweightave*seffcorr);
        if (hdphixi_fold.size()>ipw) hdphixi_fold[ipw]->Fill(dphifold,xidec,fineweightave*filltimeflowxi);
      }
    }
  }
}

void Correlation::SetHadronEfficiency(const char* filename)
{
  TFile* fhadeff = new TFile(filename);
  if( verbosity > 1 ) cout<<PHWHERE<<"loading hadron efficiency"<<fhadeff->GetName()<<endl;
  fhadeff->GetObject("feff",fhadroneff);
  //fhadroneff = (TH1D*)fhadeff->Get("heff2");
  fexemb = new TF1("fexemb","[0]+[1]*exp([2]*x)",5.0,10.0);
  fhadeff->Close();
}

float Correlation::GetHadronEfficiencyCorr(float pt)
{
  float richcorr[4]={0.680,0.835,0.925,0.975};
  //From Andrew's thesis for 0-20% without Rich embedding
  //TF1* fexemb = new TF1("fexemb","[0]+[1]*exp([2]*x)",5.0,10.0);
  fexemb->SetParameters(0.761,1.640,-4.734);
  
  float embcorr[4]={(0.779+0.851)/2.,(0.887+.938)/2.,(0.968+0.988)/2.,(.986+.993+.998)/3.};
  if(pt>5.0 && cbin==0) embcorr[0]=fexemb->Eval(pt);
  
  float seffcorr = 1.0;
  if( verbosity > 1 ) cout<<"fhadroneff->Eval(pt) = "<<fhadroneff->Eval(pt)<<endl;
  seffcorr = 2.0/fhadroneff->Eval(pt);
  if( data_set != Run8dAu ) {
    if(pt>5.0) seffcorr = seffcorr/embcorr[cbin];
    else seffcorr = seffcorr/richcorr[cbin];
  }
  return seffcorr;
}

void Correlation::SetTriggerEfficiency(const char* filename_0, const char* filename_1, const char* filename_2, const char* filename_3)
{
  fpi0eff_0=new TFile(filename_0);
  fpi0eff_1=new TFile(filename_1);
  fpi0eff_2=new TFile(filename_2);
  fpi0eff_3=new TFile(filename_3);
  
  char grname[100];
  sprintf(grname,"ratio_graph_AA_0");
  fpi0eff_0->GetObject(grname,grpi0eff_0);
  sprintf(grname,"ratio_graph_AA_1");
  fpi0eff_1->GetObject(grname,grpi0eff_1);
  sprintf(grname,"ratio_graph_AA_2");
  fpi0eff_2->GetObject(grname,grpi0eff_2);
  sprintf(grname,"ratio_graph_AA_3");
  fpi0eff_3->GetObject(grname,grpi0eff_3);
}

void Correlation::SetV2(const char* v2_inputs)
{
  TFile v2file(v2_inputs);
  cout<<PHWHERE<<" loading v2 file: "<< v2file.GetName() <<endl;
  
  TGraphErrors* gr_inc_v2[4];   
  TGraphErrors* gr_dec_v2[4];   
  TGraphErrors* gr_pi0_v2[4];   
  TGraphErrors* gr_had_v2[4];   
  TGraphErrors* gr_inc_v2sys[4];
  TGraphErrors* gr_dec_v2sys[4];
  TGraphErrors* gr_pi0_v2sys[4];
  TGraphErrors* gr_had_v2sys[4];

  for(int i=0; i<4; i++){
    gr_inc_v2[i] = NULL;
    gr_dec_v2[i] = NULL;
    gr_pi0_v2[i] = NULL;
    gr_had_v2[i] = NULL;
    gr_inc_v2sys[i] = NULL;
    gr_dec_v2sys[i] = NULL;
    gr_pi0_v2sys[i] = NULL;
    gr_had_v2sys[i] = NULL;
  }

  ostringstream name;

  for(int i=0; i<NCBINS; i++){//centrality
    name.str("");
    name << "gamma_inc_v2_"<<i;
    gr_inc_v2[i] = (TGraphErrors*)v2file.Get(name.str().c_str());
    
    double *inc = gr_inc_v2[i]->GetY();
    double *inc_err = gr_inc_v2[i]->GetEY();
    for(int j=0; j<4; j++){
      trig_v2[0][i][j] = inc[j];
      trig_v2_err[0][i][j] = inc_err[j];
    }
  
    name.str("");
    name << "gamma_inc_v2sys_"<<i;
    gr_inc_v2sys[i] = (TGraphErrors*)v2file.Get(name.str().c_str());

    double *inc_sys = gr_inc_v2sys[i]->GetEY();
    for(int j=0; j<4; j++) trig_v2_sys[0][i][j] = inc_sys[j];
    
    name.str("");
    name << "gamma_dec_v2_" << i;
    gr_dec_v2[i] = (TGraphErrors*)v2file.Get(name.str().c_str());
    
    double *dec = gr_dec_v2[i]->GetY();
    double *dec_err = gr_dec_v2[i]->GetEY();
    for(int j=0; j<4; j++){
      trig_v2[2][i][j] = dec[j];
      trig_v2_err[2][i][j] = dec_err[j];
    }

    name.str("");
    name << "gamma_dec_v2sys_" << i;
    gr_dec_v2sys[i] = (TGraphErrors*)v2file.Get(name.str().c_str());
    
    double *dec_sys = gr_dec_v2sys[i]->GetEY();
    for(int j=0; j<4; j++) trig_v2_sys[2][i][j] = dec_sys[j];
    
    name.str("");
    name << "pi0_v2_" << i;
    gr_pi0_v2[i] = (TGraphErrors*)v2file.Get(name.str().c_str());
    
    double *pi0 = gr_pi0_v2[i]->GetY();
    double *pi0_err = gr_pi0_v2[i]->GetEY();
    for(int j=0; j<4; j++){
      trig_v2[1][i][j] = pi0[j];
      trig_v2_err[1][i][j] = pi0_err[j];
    }

    name.str("");
    name << "pi0_v2sys_" << i;
    gr_pi0_v2sys[i] = (TGraphErrors*)v2file.Get(name.str().c_str());
   
    double *pi0_sys = gr_pi0_v2sys[i]->GetEY();
    for(int j=0; j<4; j++) trig_v2_sys[1][i][j] = pi0_sys[j];
    
    name.str("");
    name << "hadron_v2_" << i;
    gr_had_v2[i] = (TGraphErrors*)v2file.Get(name.str().c_str());

    double *had = gr_had_v2[i]->GetY();
    double *had_err = gr_had_v2[i]->GetEY();
    for(int j=0; j<5; j++){
      part_v2[i][j] = had[j];
      part_v2_err[i][j] = had_err[j];
    }

    name.str("");
    name << "hadron_v2sys_" << i;
    gr_had_v2sys[i] = (TGraphErrors*)v2file.Get(name.str().c_str());

    double *had_sys = gr_had_v2sys[i]->GetEY();
    for(int j=0; j<5; j++) part_v2_sys[i][j] = had_sys[j];
  }

  for(int i=0; i<4; i++){
    delete gr_inc_v2[i];
    delete gr_dec_v2[i];
    delete gr_pi0_v2[i];
    delete gr_had_v2[i];
    delete gr_inc_v2sys[i];
    delete gr_dec_v2sys[i];
    delete gr_pi0_v2sys[i];
    delete gr_had_v2sys[i];
  }
  v2file.Close(); 
}

void Correlation::SetSharkFin(const char* filename)
{
  TFile *fshark_exodus=NULL;
  fshark_exodus = new TFile(filename);
  
  //now take 10 cm projections out to 165 cm
  //right now use bins out to 170 since we have two cm bins instead of 1.
  cout<<PHWHERE<<" loading sharkfins "<< fshark_exodus->GetName() <<endl;
  
  for(int izemc=0;izemc<33;izemc++){
    char sharkname[100];
    sprintf(sharkname,"PTpi_PTgam_%d",izemc);
    ptpivsptgam[izemc]=(TH2D*)fshark_exodus->Get(sharkname);
    
    for(int idecl=0;idecl<5;idecl++){
      sprintf(sharkname,"hshark_large_%d_%d",idecl,izemc);
      hshark_large[idecl][izemc]=(TH1D*)fshark_exodus->Get(sharkname);
    }
  }
  
  cout << "applying pisacorr " << endl;
  TF1 * pisacorr_large[5];
  pisacorr_large[0]=new TF1("pisacorr_large0","-7.85e-2*x+1.85",4,20);
  pisacorr_large[1]=new TF1("pisacorr_large1","-8.7e-2*x+2.0",4,20);
  pisacorr_large[2]=new TF1("pisacorr_large2","-9.65e-2*x+2.2",4,20);
  pisacorr_large[3]=new TF1("pisacorr_large3","-12.5e-2*x+2.8",4,20);
  pisacorr_large[4]=new TF1("pisacorr_large4","-12.5e-2*x+2.8",4,20);
  
  for (int ipzemc = 0; ipzemc < 33; ipzemc++){
    for (int ipdecs = 0; ipdecs < 5; ipdecs++){
      for (int ipbin = 1; ipbin < hshark_large[ipdecs][ipzemc]->GetNbinsX()+1; ipbin++){
        double mattshark = hshark_large[ipdecs][ipzemc]->GetBinContent(ipbin);
        double pi0pt =  hshark_large[ipdecs][ipzemc]->GetBinCenter(ipbin);
        double pisacorr=pisacorr_large[ipdecs]->Eval(pi0pt);
        if(pisacorr>1.0) pisacorr=1.0;
        double mwweightfine=mattshark*pisacorr;
        //cout<<"ipzemc: "<<ipzemc<<"; ipdecs: "<<ipdecs<<"; mwweightfine: "<<mwweightfine<<endl;
        hshark_large[ipdecs][ipzemc]->SetBinContent(ipbin,mwweightfine);
      }
    }
  }
  for(int ij=0; ij<5; ij++)
    delete pisacorr_large[ij];
}

Correlation::DataSet Correlation::GetDataSet(int RunNumber)
{
  DataSet dataset = INVALID;
  if(RunNumber>=246000 && RunNumber<=254000)
    dataset = Run8dAu;
  if(RunNumber>=300475 && RunNumber<=310454)
    dataset = Run10AuAu;
  if(RunNumber>=343031 && RunNumber<=349680)
    dataset = Run11AuAu;

  return dataset;
}

float Correlation::GetNTriggers(TH1F* trigpt, float trigptmin, float trigptmax)
{
  float ntrig = 0;
  for(int i = 0; i < trigpt->GetNbinsX(); i++)
  {
    float bincenter = trigpt->GetBinCenter(i+1);
    if(bincenter >=trigptmin && bincenter <= trigptmax)
    {
      float bincontent = trigpt->GetBinContent(i+1);
      ntrig += bincontent;
    }
  }
  return ntrig;
}

void Correlation::Clear()
{
  
  for (unsigned int i = 0; i < clus_everything.size(); i++){
    delete clus_everything[i];
  }
  clus_everything.clear();
  //   for (unsigned int i = 0; i < clus_vector_novetotracks.size(); i++){
  //     delete clus_vector_novetotracks[i];
  //   }
  //   clus_vector_novetotracks.clear();
  for (unsigned int i = 0; i < clus_vector.size(); i++){
    delete clus_vector[i];
  }
  clus_vector.clear();
  for(unsigned int i = 0; i < trk_vector.size(); i++){
    delete trk_vector[i];
  }
  trk_vector.clear();
  for(unsigned int i = 0; i < trk_vector_05sig.size(); i++){
    delete trk_vector_05sig[i];
  }
  trk_vector_05sig.clear();
  for(unsigned int i = 0; i < trk_vector_1sig.size(); i++){
    delete trk_vector_1sig[i];
  }
  trk_vector_1sig.clear();
  for(unsigned int i = 0; i < trk_vector_15sig.size(); i++){
    delete trk_vector_15sig[i];
  }
  trk_vector_15sig.clear();
  for(unsigned int i = 0; i < trk_vector_2sig.size(); i++){
    delete trk_vector_2sig[i];
  }
  trk_vector_2sig.clear();
  for(unsigned int i = 0; i < trk_vector_3sig.size(); i++){
    delete trk_vector_3sig[i];
  }
  trk_vector_3sig.clear();
  for(unsigned int i = 0; i < pi0_vector.size(); i++){
    delete pi0_vector[i];
  }
  pi0_vector.clear();
  for(unsigned int i = 0; i < bgpi0_vector.size(); i++){
    delete bgpi0_vector[i];
  }
  bgpi0_vector.clear();
  for(unsigned int i = 0; i < bgclus_vector.size(); i++){
    delete bgclus_vector[i];
  }
  bgclus_vector.clear();
  for(unsigned int i = 0; i < all_clus_vector.size(); i++){
    delete all_clus_vector[i];
  }
  all_clus_vector.clear();
  for(unsigned int i = 0; i < lessqualtrk_vector.size(); i++){
    delete lessqualtrk_vector[i];
  }
  lessqualtrk_vector.clear();
}

int Correlation::End(PHCompositeNode* topNode)
{
  cout<<"End."<<endl;
  DoMixing(atree->_ttrig, atree->_tpart, NMIX);
  cout<<"finish mixing."<<endl;
  manager->dumpHistos(output.c_str());
  cout<<"total # of events: "<<evt<<endl;
  cout<<"total # of events processed: "<<event<<endl;
  return 0;
}

void Correlation::InitHistos()
{
  h2_dphi = NULL;
  h2_dphi_xi = NULL;
  h2_dphi_accw = NULL;
  h2_dphi_accw_xi = NULL;
  h3_dphi_accw = NULL;
  h3_dphi_accw_xi = NULL;
  h2_partpt_xi = NULL;

  h2_dphi_mix = NULL;
  h2_dphi_xi_mix = NULL;
  h2_dphi_accw_mix = NULL;
  h2_dphi_accw_xi_mix = NULL;
  h3_dphi_accw_mix = NULL;
  h3_dphi_accw_xi_mix = NULL;
  h2_partpt_xi_mix = NULL;

  h2_dphi_pi0 = NULL;
  h2_dphi_xi_pi0 = NULL;
  h2_dphi_accw_pi0 = NULL;
  h2_dphi_accw_xi_pi0 = NULL;
  h3_dphi_accw_pi0 = NULL;
  h3_dphi_accw_xi_pi0 = NULL;
  h2_partpt_xi_pi0 = NULL;

  h2_dphi_pi0_mix = NULL;
  h2_dphi_xi_pi0_mix = NULL;
  h2_dphi_accw_pi0_mix = NULL;
  h2_dphi_accw_xi_pi0_mix = NULL;
  h3_dphi_accw_pi0_mix = NULL;
  h3_dphi_accw_xi_pi0_mix = NULL;
  h2_partpt_xi_pi0_mix = NULL;

  for(int ic=0; ic<4; ic++){
    for(int it=0; it<4; it++){
      for(int ip=0; ip<5; ip++){
        IncAcc[ic][it][ip] = NULL;
        Pi0Acc[ic][it][ip] = NULL;
        DecAcc[ic][it][ip] = NULL;
        IncAccIso[ic][it][ip] = NULL;
        Pi0AccIso[ic][it][ip] = NULL;
        DecAccIso[ic][it][ip] = NULL;
      }
      for(int ix=0; ix<6; ix++){
        IncAccXi[ic][it][ix] = NULL;
        Pi0AccXi[ic][it][ix] = NULL;
        DecAccXi[ic][it][ix] = NULL;
        IncAccXiIso[ic][it][ix] = NULL;
        Pi0AccXiIso[ic][it][ix] = NULL;
        DecAccXiIso[ic][it][ix] = NULL;
      }
    }
  }
  h2_pi0mass_PbGl = NULL;
  h2_pi0mass_PbSc = NULL;

  h2_pi0_bg = NULL;
  h2_pi0_bg_PbGl = NULL;
  h2_pi0_bg_PbSc = NULL;

  for(int ias=0; ias<N_ARMSECT; ias++){
    h2_pi0mass_as[ias] = NULL;
    h2_pi0_bg_as[ias] = NULL;
    h3_nhit[ias] = NULL;
    h3_nhit_mywarn[ias] = NULL;
    h1_chi2[ias] = NULL;
    h3_disp[ias] = NULL;
  }

  h2_emctrkdphi_bf = NULL;
  h2_emctrkdz_bf = NULL;
  h2_emcpc3dphi_bf = NULL;
  h2_emcpc3dz_bf = NULL;
  h2_emctrkdphi_aft = NULL;
  h2_emctrkdz_aft = NULL;
  h2_emcpc3dphi_aft = NULL;
  h2_emcpc3dz_aft = NULL;
  h2_pc3sdphi_bf = NULL;
  h2_pc3sdz_bf = NULL;
  h2_pc3sdphi_aft = NULL;
  h2_pc3sdz_aft = NULL;
  h2_phi_zed_bf = NULL;
  h2_phi_zed_aft = NULL;
  h2_ppc3_west_bf = NULL;
  h2_ppc3_east_bf = NULL;
  h2_ppc3_west_aft = NULL;
  h2_ppc3_east_aft = NULL;
  h2_ppc1_west_bf = NULL;
  h2_ppc1_east_bf = NULL;
  h2_ppc1_west_aft = NULL;
  h2_ppc1_east_aft = NULL;
  h2_EoverPvspt = NULL;
  h1_zvertex = NULL;
  // h1_centrality = NULL;
  h1_trig_pt_inc_tot = NULL;
  h1_trig_pt_pi0_tot = NULL;
  h1_trig_pt_dec_tot = NULL;
  h1_trig_pt_inc_mix_tot = NULL;
  h1_trig_pt_pi0_mix_tot = NULL;
  h1_trig_pt_dec_mix_tot = NULL;
  h1_part_pt_tot = NULL;
  h1_part_pt_05sig = NULL;
  h1_part_pt_1sig = NULL;
  h1_part_pt_15sig = NULL;
  h1_part_pt_2sig = NULL;
  h1_part_pt_3sig = NULL;
  h2_ptvscent_trig_inc = NULL;
  h2_ptvscent_part = NULL;
  h2_ptvscent_trig_pi0 = NULL;
  h3_mintrackdist_fg_allcent = NULL;
  h3_mintrackdist_bg_allcent = NULL;
  h2_bfpaircut_inc = NULL;
  h2_aftpaircut_inc = NULL;
  h2_bfpaircut_pi0 = NULL;
  h2_aftpaircut_pi0 = NULL;
  h3_pt_phi_eta_clus = NULL;
  h3_pt_phi_eta_trk = NULL;
}

void Correlation::DoMixing(TTree* trig, TTree* assoc, int size)
{
  int pooldepth = 0;
  int nloop = 0;

  int evt_trig;
  float zvtx_trig;
  float cent_trig;
  int ntrig_photons;
  int ntrig_pi0s;
  int ntrig = DIM;
  float trigpt[ntrig];
  float trigphi[ntrig];
  float trigeta[ntrig];
  float trige[ntrig];
  float x[ntrig];
  float y[ntrig];
  float z[ntrig];
  int iso[ntrig];

  int evt_part;
  float zvtx_part;
  float cent_part;
  int npart = DIM;
  int nphotons;
  float partpt[npart];
  float partphi[npart];
  float parteta[npart];
  float parte[npart];
  float pemcx[npart];
  float pemcy[npart];
  float pemcz[npart];
  float cluspt[npart];
  float clusphi[npart];
  float cluseta[npart];
  float cluse[npart];

  trig->SetBranchAddress("evt", &evt_trig);
  trig->SetBranchAddress("zvtx", &zvtx_trig);
  trig->SetBranchAddress("cent", &cent_trig);
  trig->SetBranchAddress("nphotons", &ntrig_photons);
  trig->SetBranchAddress("npi0s", &ntrig_pi0s);
  trig->SetBranchAddress("ntrig", &ntrig);
  trig->SetBranchAddress("pt", trigpt);
  trig->SetBranchAddress("phi", trigphi);
  trig->SetBranchAddress("eta", trigeta);
  trig->SetBranchAddress("e",trige);
  trig->SetBranchAddress("x",x);
  trig->SetBranchAddress("y",y);
  trig->SetBranchAddress("z",z);
  trig->SetBranchAddress("iso",iso);

  assoc->SetBranchAddress("evt", &evt_part);
  assoc->SetBranchAddress("zvtx", &zvtx_part);
  assoc->SetBranchAddress("cent", &cent_part);
  assoc->SetBranchAddress("ntracks", &npart);
  assoc->SetBranchAddress("nphotons", &nphotons);
  assoc->SetBranchAddress("pt", partpt);
  assoc->SetBranchAddress("phi", partphi);
  assoc->SetBranchAddress("eta", parteta);
  assoc->SetBranchAddress("e", parte);
  assoc->SetBranchAddress("pemcx", pemcx);
  assoc->SetBranchAddress("pemcy", pemcy);
  assoc->SetBranchAddress("pemcz", pemcz);
  assoc->SetBranchAddress("pt_clus", cluspt);
  assoc->SetBranchAddress("phi_clus", clusphi);
  assoc->SetBranchAddress("eta_clus", cluseta);
  assoc->SetBranchAddress("e_clus", cluse);

  int nentrig = trig->GetEntries();
  int nenpart = assoc->GetEntries();
  //cout<<"Mixing: nentrig = " <<nentrig<<"; nenpart = "<<nenpart<<endl;
 
  
  for(int i=0; i<nentrig; i++){
    //cout<<"trigger count = "<<i<<endl;
    pooldepth = 0;
    nloop = 0;
    trig->GetEntry(i);
    //cout <<"evt_trig = "<<evt_trig<<"; zvtx_trig = "<<zvtx_trig<<"; cent_trig = "<<cent_trig<<"; ntrig_photons = "<<ntrig_photons<<"; ntrig_pi0s = "<<ntrig_pi0s<<"; ntrig = "<<ntrig<<endl;
    cbin = GetCentBin((int)cent_trig);
    //cout<<"Mixed event cbin = "<<cbin<<endl;

    vector<ACluster*> photons;
    vector<APiZero*> pi0s; 

    for(int itrig=0; itrig<ntrig_pi0s; itrig++){
      //cout<<"pi0 itrig = "<<itrig<<"; trigpt["<<itrig<<"] = "<<trigpt[itrig]<<"; trigphi["<<itrig<<"] = "<<trigphi[itrig]<<endl;
      APiZero pi0;
      MakePi0Object(&pi0, trigpt[itrig], trigphi[itrig], trigeta[itrig], trige[itrig], x[itrig], y[itrig], z[itrig], zvtx_trig);
      if( verbosity ) { cout<<"Mixing: pt = "<<pi0.Pt()<<"; phi = "<<pi0.Phi()<<"; eta = "<<pi0.Eta()<<"; e = "<<pi0.E()<<"; x = "<<((ACluster*)pi0.Daughter1())->GetX()<<"; y = "<<((ACluster*)pi0.Daughter1())->GetY()<<"; z = "<<((ACluster*)pi0.Daughter1())->GetZ()<<endl;}
      pi0.SetIso(iso[itrig]);

      //dec weighting
      vector<float> mwweight;
      for(int i=0; i<5; i++) mwweight.push_back(0.0);
      EvalDecWeights(&pi0,zvtx_trig,cbin,mwweight);
      pi0.SetDecayWeights(mwweight);

      pi0s.push_back(pi0.clone());
      if( verbosity ) { cout<<"Mixing: after making pi0 vector. pt = "<<pi0.Pt()<<"; phi = "<<pi0.Phi()<<"; eta = "<<pi0.Eta()<<"; e = "<<pi0.E()<<"; x = "<<((ACluster*)pi0.Daughter1())->GetX()<<"; y = "<<((ACluster*)pi0.Daughter1())->GetY()<<"; z = "<<((ACluster*)pi0.Daughter1())->GetZ()<<endl; }
    }
    //cout<<"DoMixing: made pi0s. pi0 vector size: "<<pi0s.size()<<endl; 

    for(int itrig=ntrig_pi0s; itrig<ntrig; itrig++){
      //cout<<"photon itrig = "<<itrig<<"; trigpt["<<itrig<<"] = "<<trigpt[itrig]<<"; trigphi["<<itrig<<"] = "<<trigphi[itrig]<<endl;
      ACluster pho;
      MakeClusterObject(&pho, trigpt[itrig], trigphi[itrig], trigeta[itrig], trige[itrig], x[itrig], y[itrig], z[itrig], zvtx_trig);
      pho.SetIso(iso[itrig]);
      //cout<<"Mixing: pho.Pt() = "<<pho.Pt()<<"; pho.Phi() = "<<pho.Phi()<<endl;
      photons.push_back(pho.clone());
    }
    //cout<<"DoMixing: made photons. photon vector size: "<<photons.size()<<endl;

    float nvert_fg = 0.;
    float ncent_fg = 0.;
    nvert_fg = TMath::Floor(zvtx_trig/5.0);
    ncent_fg = TMath::Floor((cent_trig-1)/5);
    //cout<<"nvert_fg = "<<nvert_fg<<"; ncent_fg = "<<ncent_fg<<endl;
    
    for(int j=0; j<nenpart; j++){
      
      assoc->GetEntry(j);
      //cout <<"evt_part = "<<evt_part<<"; zvtx_part = "<<zvtx_part<<"; cent_part = "<<cent_part<<"; nphotons = "<<nphotons<<"; npart = "<<npart<<endl;

      //cout<<"evt_trig = "<<evt_trig<<"; evt_part = "<<evt_part<<endl;
      //check if trigger and assoc belong to the same event
      if(evt_trig == evt_part){
        j = CheckPool(nenpart,j,pooldepth,size,nloop);
        //cout<<"same event. nloop = "<<nloop<<endl;
        if(nloop > NMIX) break;//to prevent infinite loop in case there is no interesting hadrons available 
        continue;
      }
      
      //check if trigger and partner belongs to the same event class
      float nvert_bg = TMath::Floor(zvtx_part/5.0);
      float ncent_bg = TMath::Floor((cent_part-1)/5);
      //cout<<"nvert_bg = "<<nvert_bg<<"; ncent_bg = "<<ncent_bg<<endl;
      
      if(nvert_fg != nvert_bg) {
        j = CheckPool(nenpart,j,pooldepth,size,nloop);
        //cout<<"not same vertex bin. nloop = "<<nloop<<endl;
        if(nloop > NMIX) break;
        continue;
      }
      
      if(ncent_fg != ncent_bg) {
        j = CheckPool(nenpart,j,pooldepth,size,nloop);
        //cout<<"not same centrality bin. nloop = "<<nloop<<endl;
        if(nloop > NMIX) break;
        continue;
      }
      
      vector<ATrack*> hadrons;
      for(int ipart=0; ipart<npart; ipart++){
        ATrack trk;
        MakeTrackObject(&trk, partpt[ipart], partphi[ipart], parteta[ipart], parte[ipart], pemcx[ipart], pemcy[ipart], pemcz[ipart], zvtx_part);
        hadrons.push_back(trk.clone());
      }
     
      vector<ACluster*> clusters;
      for( int iclust=0; iclust<nphotons; iclust++ )
      {
        ACluster pho;
        MakeClusterObject(&pho, cluspt[iclust], clusphi[iclust], cluseta[iclust], cluse[iclust], 0, 0, 0, zvtx_part);
        clusters.push_back(pho.clone());
      }
      
      for( unsigned int itrig = 0; itrig < photons.size(); itrig++ ) {
        // Want to use fg isolated acceptance but apply additional isolation based on uncorrelated (mixed) particles
        if( photons[itrig]->IsIso() )
          SetIso(photons[itrig],hadrons,clusters,Rcut,econe_min[cbin],h3_cluster_mix_dR[cbin],h3_cluster_mix_etot[cbin],h2_cluster_mix_wdR[cbin],h2_cluster_mix_etot[cbin],h3_iso_mix_acc[cbin]);
        h1_trig_pt_inc_mix[cbin]->Fill(photons[itrig]->Pt());
        if( photons[itrig]->IsIso() )
          h1_trig_pt_inc_iso_mix[cbin]->Fill(photons[itrig]->Pt());
        h1_trig_pt_inc_mix_tot->Fill(photons[itrig]->Pt()); 
      }
      
      for( unsigned int itrig = 0; itrig < pi0s.size(); itrig++ ) {
        // Want to use fg isolated acceptance but apply additional isolation based on uncorrelated (mixed) particles
        if( pi0s[itrig]->IsIso() )
          SetIso(pi0s[itrig],hadrons,clusters,Rcut,econe_min[cbin]);
        h1_trig_pt_pi0_mix[cbin]->Fill(pi0s[itrig]->Pt());
        if( pi0s[itrig]->IsIso() )
          h1_trig_pt_pi0_iso_mix[cbin]->Fill(pi0s[itrig]->Pt());
        h1_trig_pt_pi0_mix_tot->Fill(pi0s[itrig]->Pt());
        
        //counting dec triggers
        for(int ipw=0; ipw<5; ipw++){
          double weight = pi0s[itrig]->GetDecayWeights()[ipw];
          h1_trig_pt_dec_mix[cbin]->Fill(ipw,weight);
          if( pi0s[itrig]->IsIso() ) 
            h1_trig_pt_dec_iso_mix[cbin]->Fill(ipw,weight);
          h1_trig_pt_dec_mix_tot->Fill(ipw,weight);
        }
      }
      
      if( verbosity>3 ) cout <<"DoMixing: made hadrons. hadron vector size: "<<hadrons.size()<<endl;

      pooldepth++;
      if( verbosity>3 ) cout<<"pooldepth = "<<pooldepth<<endl;

      if( DiagFlag ) {
        MakePairs(photons,hadrons,MIX,data_set,0,h3_dphi_mix[cbin],h3_dphi_mix_fold[cbin],h3_ptxidphi_mix[cbin],h3_ptxidphi_mix_fold[cbin],h3_ptztdphi_mix[cbin],h3_ptztdphi_mix_fold[cbin],vector<TH2F*>(),vector<TH2F*>(),vector<TH2F*>(),vector<TH2F*>(),vector<TH2F*>(),vector<TH2F*>(),h2_bfpaircut_inc,h2_aftpaircut_inc,h2_dphi_mix,h2_dphi_xi_mix,h2_dphi_accw_mix,h2_dphi_accw_xi_mix,h3_dphi_accw_mix,h3_dphi_accw_xi_mix,h2_partpt_xi_mix);
        MakePairs(pi0s,hadrons,MIXPI,data_set,0,h3_dphi_pi0_mix[cbin],h3_dphi_pi0_mix_fold[cbin],h3_ptxidphi_pi0_mix[cbin],h3_ptxidphi_pi0_mix_fold[cbin],h3_ptztdphi_pi0_mix[cbin],h3_ptztdphi_pi0_mix_fold[cbin],h2_dphi_dec_mix[cbin],h2_dphi_dec_mix_fold[cbin],h2_dphixi_dec_mix[cbin],h2_dphixi_dec_mix_fold[cbin],h2_dphizt_dec_mix[cbin],h2_dphizt_dec_mix_fold[cbin],h2_bfpaircut_pi0,h2_aftpaircut_pi0,h2_dphi_pi0_mix,h2_dphi_xi_pi0_mix,h2_dphi_accw_pi0_mix,h2_dphi_accw_xi_pi0_mix,h3_dphi_accw_pi0_mix,h3_dphi_accw_xi_pi0_mix,h2_partpt_xi_pi0_mix);
      }
      else {
        MakePairs(photons,hadrons,MIX,data_set,0,h3_dphi_mix[cbin],h3_dphi_mix_fold[cbin],h3_ptxidphi_mix[cbin],h3_ptxidphi_mix_fold[cbin],h3_ptztdphi_mix[cbin],h3_ptztdphi_mix_fold[cbin]);
        MakePairs(pi0s,hadrons,MIXPI,data_set,0,h3_dphi_pi0_mix[cbin],h3_dphi_pi0_mix_fold[cbin],h3_ptxidphi_pi0_mix[cbin],h3_ptxidphi_pi0_mix_fold[cbin],h3_ptztdphi_pi0_mix[cbin],h3_ptztdphi_pi0_mix_fold[cbin],h2_dphi_dec_mix[cbin],h2_dphi_dec_mix_fold[cbin],h2_dphixi_dec_mix[cbin],h2_dphixi_dec_mix_fold[cbin],h2_dphizt_dec_mix[cbin],h2_dphizt_dec_mix_fold[cbin]);
      }
      MakePairs(photons,hadrons,MIX,data_set,1,h3_dphi_iso_mix[cbin],h3_dphi_iso_mix_fold[cbin],h3_ptxidphi_iso_mix[cbin],h3_ptxidphi_iso_mix_fold[cbin],NULL,NULL);     
      MakePairs(pi0s,hadrons,MIXPI,data_set,1,h3_dphi_pi0_iso_mix[cbin],h3_dphi_pi0_iso_mix_fold[cbin],h3_ptxidphi_pi0_iso_mix[cbin],h3_ptxidphi_pi0_iso_mix_fold[cbin],NULL,NULL,h2_dphi_dec_iso_mix[cbin],h2_dphi_dec_iso_mix_fold[cbin],h2_dphixi_dec_iso_mix[cbin],h2_dphixi_dec_iso_mix_fold[cbin],NULL,NULL);

      //for(unsigned int i=0; i<hadrons.size(); i++) delete hadrons[i];
      //hadrons.clear();
      //for(unsigned int i=0; i<clusters.size(); i++) delete hadrons[i];
      //clusters.clear();
      
      ClearVector(hadrons);
      ClearVector(clusters);
      
      if(pooldepth == size) {/*cout<<"Mixed enough! pooldepth = "<<pooldepth<<"; nvert_fg = "<<nvert_fg<<"; ncent_fg = "<<ncent_fg<<endl; */break;}

      //make sure making NMIX pairs
      j = CheckPool(nenpart,j,pooldepth,size,nloop);
    }
    ClearVector(photons);
    ClearVector(pi0s);
    
    //for(unsigned int i=0; i<photons.size(); i++) delete photons[i];
    //photons.clear();
    //for(unsigned int i=0; i<pi0s.size(); i++) delete pi0s[i];
    //pi0s.clear();
  }
}

int Correlation::CheckPool(int nenpart, int j, int pooldepth, int size, int& nloop)
{
  // cout<<"in CheckPool: before nloop = "<<nloop<<endl;
  if(j == nenpart-1){
    if(pooldepth < size){
      //cout<<"pool not deep enough."<<endl;
      j = -1;
      nloop++;
      //cout<<"after nloop = "<<nloop<<endl;
    }
  }
  return j;
}

float Correlation::GetDecayXiWeights(int decxibw, int itdec, float xidec, int ipi0zemc, float trigpt, float partpt)//to get weights in xi-binned histograms for decays
{
  //cout << "In GetDecayXiWeights: " << endl;

  int trigbinbounds[5] = {5, 7, 9, 12, 15};

  float fineweightave=0.;

  float xideclo = xidec - decxibw/2.;
  float xidechi = xidec + decxibw/2.;

  float ptdeclo = fabs(partpt*exp(xideclo));
  float ptdechi = fabs(partpt*exp(xidechi));

  if(ptdechi < trigbinbounds[itdec] || ptdeclo > trigbinbounds[itdec+1]) return 0.0;
  if(ptdeclo < trigbinbounds[itdec]) ptdeclo = trigbinbounds[itdec];
  if(ptdechi > trigbinbounds[itdec+1]) ptdechi = trigbinbounds[itdec+1];
  
  int ipdecbinlo = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdeclo);       
  int ipdecbinhi = ptpivsptgam[ipi0zemc]->GetYaxis()->FindBin(ptdechi);       

  int trigptbin_shark = hshark_large[0][0]->FindBin(trigpt);
  for(int idecbin=ipdecbinlo;idecbin<ipdecbinhi;idecbin++)
    fineweightave+=ptpivsptgam[ipi0zemc]->GetBinContent(trigptbin_shark,idecbin);
  
  //take into account the phase space for decay
  fineweightave*=(ptdechi-ptdeclo);
  return fineweightave;
}

int Correlation::GetPi0ZEMCBin(APiZero* pizero)
{
  float pt = pizero->Pt();
  float pz = pizero->Pz();
  float zvertex = pizero->GetZvtx();

  float pi0zemc = 510.0*pz/pt+zvertex;
  int ipi0zemc = (int)TMath::Floor((pi0zemc+165.0)/10.);

  if(ipi0zemc<0)ipi0zemc=0;
  if(ipi0zemc>32)ipi0zemc=32; 
  
  return ipi0zemc;   
}
