void plot_integrated_raa_2_5()
{
  int   num_pi0 = 8;
  float n_part_pi0[8]  = {325.2, 234.6, 166.6, 114.2, 
			  74.40, 45.50,  25.7,   9.5};
  float en_part_pi0[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
  float raa_pi0[8]     = {0.294671, 0.369870, 0.439591, 0.521614, 
			  0.611563, 0.661897, 0.705451, 0.811798};
  float stat_pi0[8]    = {0.0115225, 0.0137082, 0.0161476, 0.0189241,
			  0.0215431, 0.0244435, 0.0273755, 0.0304532};
  float sys_pi0[8]     = {0.0282595, 0.0348328, 0.0414679, 0.0492902,
			  0.0565685, 0.0610058, 0.0652277, 0.0747985};
  float sys_taa_pi0[8] = {0.0202060, 0.0257749, 0.0355624, 0.0438834,
			  0.0598733, 0.1049800, 0.1867370, 0.2435390};
  float e_raa_pi0_cent[8];
  for (int ibin=0; ibin<num_pi0; ibin++){
    e_raa_pi0_cent[ibin] = sqrt( sys_pi0[ibin]**2 
				 + sys_taa_pi0[ibin]**2);
  }
  float e_raa_pi0_nocent[8];
  for (int ibin=0; ibin<num_pi0; ibin++){
    e_raa_pi0_nocent[ibin] = 0.10 * raa_pi0[ibin];
  }
  
  int   num_e = 4;
  float n_part_e[4]       = {325.2, 234.6, 140.4, 60.0};
  float en_part_e[4]      = {0., 0., 0., 0.};
  float raa_e[4]          = {0.342, 0.610, 0.512, 0.654};
  float stat_e[4]         = {0.046, 0.068, 0.066, 0.127};
  float e_raa_e_cent[4]   = {0.054, 0.093, 0.084, 0.116};
  float e_raa_e_nocent[4] = {0.142, 0.244, 0.215, 0.284};
  /*
  for (int ibin=0; ibin<4; ibin++){
    stat_e[ibin] = sqrt(stat_e[ibin]*stat_e[ibin] + 
			e_raa_e_cent[ibin]*e_raa_e_cent[ibin]);
  }
  */
  TGraphErrors *raa_pizero_stat = new TGraphErrors(num_pi0,n_part_pi0,raa_pi0,en_part_pi0,stat_pi0);
  TGraphErrors *raa_pizero_cent = new TGraphErrors(num_pi0,n_part_pi0,raa_pi0,en_part_pi0,e_raa_pi0_cent);
  TGraphErrors *raa_pizero_nocent = new TGraphErrors(num_pi0,n_part_pi0,raa_pi0,en_part_pi0,e_raa_pi0_nocent);
  raa_pizero_stat->SetMarkerColor(2);
  raa_pizero_stat->SetMarkerStyle(21);
  raa_pizero_cent->SetMarkerColor(2);
  raa_pizero_cent->SetMarkerStyle(21);
  raa_pizero_nocent->SetMarkerColor(2);
  raa_pizero_nocent->SetMarkerStyle(21);

  TGraphErrors *raa_electron_stat = new TGraphErrors(num_e,n_part_e,raa_e,en_part_e,stat_e);
  TGraphErrors *raa_electron_cent = new TGraphErrors(num_e,n_part_e,raa_e,en_part_e,e_raa_e_cent);
  TGraphErrors *raa_electron_nocent = new TGraphErrors(num_e,n_part_e,raa_e,en_part_e,e_raa_e_nocent);
  raa_electron_stat->SetMarkerColor(4);
  raa_electron_stat->SetMarkerStyle(21);
  raa_electron_cent->SetMarkerColor(4);
  raa_electron_cent->SetMarkerStyle(21);
  raa_electron_nocent->SetMarkerColor(4);
  raa_electron_nocent->SetMarkerStyle(21);

  c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  c1->SetFillColor(10);
  //c1->SetGrid();
  c1->SetTicks();
  c1->SetLogy(0);

  TH1F *h1=new TH1F("h1","",350,0.,350);
  h1->SetMaximum(1.5);
  h1->SetMinimum(0.);
  h1->SetYTitle("R_{AA}^{2.5 - 5.0}");
  h1->SetXTitle("N_{part}");
  h1->Draw();
  TF1 *one = new TF1("one","1",0.,400.);
  one->Draw("same");

  draw_box(raa_electron_nocent);
  draw_box(raa_pizero_nocent);

  raa_pizero_stat->Draw("P");
  draw_bracket(raa_pizero_cent);

  raa_electron_stat->Draw("P");
  draw_bracket(raa_electron_cent);

  //  raa_electron_stat->Fit("pol1");

  TLatex *electron = new TLatex(200.,1.3,"(e^{+}+e^{-})/2: R_{AA}^{2.5 - 5.0}");
  electron->SetTextColor(4);
  electron->SetTextSize(0.06);
  electron->Draw();

  TLatex *pizero = new TLatex(200.,1.15,"#pi^{0}: R_{AA}^{2.5 - 5.0}");
  pizero->SetTextColor(2);
  pizero->SetTextSize(0.06);
  pizero->Draw();

  TLatex *bracket = new TLatex(170.,0.90,"bracket: point-by-point sys. error");
  bracket->SetTextColor(1);
  bracket->SetTextSize(0.04);
  bracket->Draw();
  TLatex *box = new TLatex(170.,0.80,"box: centrality indep. sys. error");
  box->SetTextColor(1);
  box->SetTextSize(0.04);
  box->Draw();

  return;
}

void draw_bracket(TGraphErrors *graph)
{
  int n = graph->GetN();
  TPolyLine *upper_bracket;
  TPolyLine *lower_bracket;
  TArrow *arrows;

  for(int i=0;i<n;i++) {
    Double_t dy0=0.02;
    Double_t dx0=5.0;
    Double_t xp[4];
    Double_t yp[4];
    Double_t yerr = graph->GetErrorY(i);
    Double_t xval;
    Double_t yval;

    graph->GetPoint(i,xval,yval);
    xp[0]= xval - dx0;
    yp[0]= (yval+yerr) - dy0;
    xp[1]= xval - dx0;
    yp[1]= (yval+yerr);
    xp[2]= xval + dx0;
    yp[2]= (yval+yerr);
    xp[3]= xval + dx0;
    yp[3]= (yval+yerr) - dy0;

    upper_bracket = new TPolyLine(4,xp,yp);
    upper_bracket->Draw();

    Double_t xp2[4];
    Double_t yp2[4];
    if(yval - yerr > 0) {
      xp2[0]= xval - dx0;
      yp2[0]= (yval-yerr) + dy0;
      xp2[1]= xval - dx0;
      yp2[1]= (yval-yerr);
      xp2[2]= xval + dx0;
      yp2[2]= (yval-yerr);
      xp2[3]= xval + dx0;
      yp2[3]= (yval-yerr) + dy0;

      lower_bracket = new TPolyLine(4,xp2,yp2);
      lower_bracket->Draw();
    } else {
      Double_t ya1=(yval+yerr);
      Double_t ya2= ya1 - 1.0;
      arrows=new TArrow(xval,ya1,xval,ya2);
      arrows->SetArrowSize(0.03);
      arrows->Draw();
    }
  }
}

void draw_box(TGraphErrors *graph)
{
  int n = graph->GetN();
  TBox *errorbox;

  for(int i=0;i<n;i++) {
    Double_t dy0=0.015;
    Double_t dx0=3.0;
    Double_t xp[2];
    Double_t yp[2];
    Double_t yerr = graph->GetErrorY(i);
    Double_t xval;
    Double_t yval;

    graph->GetPoint(i,xval,yval);
    xp[0]= xval - dx0;
    yp[0]= yval + yerr;
    xp[1]= xval + dx0;
    yp[1]= yval - yerr;

    errorbox = new TBox(xp[0],yp[0],xp[1],yp[1]);
    if (graph->GetMarkerColor() == 2) errorbox->SetFillColor(8);
    if (graph->GetMarkerColor() == 4) errorbox->SetFillColor(5);
    errorbox->Draw();
  }
}
