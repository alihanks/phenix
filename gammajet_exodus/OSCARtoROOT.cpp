//-----------------------------------------------------------------------------
//
//  Translate OSCAR compliant ASCII file into a ROOT Ntuple
//
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>

TROOT OSCARtoROOT("OSCARtoROOT","Converter OSCAR to ROOT");

using namespace std;

int main()
{
  char   input_file[100], output_file[100];
  char   character;
  string line;
  bool   file_end = false;
  int    event_id, max_events;
  int    particles_per_event;
  int    zero;
  int    pnum, pid;
  float  px, py, pz, E, mass, xvtx, yvtx, zvtx, opt;
  float  Pi, theta, phi, rap, eta, mt, mom;

  cout << endl << endl;
  cout << "**********************************" << endl;
  cout << "*                                *" << endl;
  cout << "*  W E L C O M E to OSCARtoROOT  *" << endl;
  cout << "*                                *" << endl;
  cout << "*   Translate OSCAR compliant    *" << endl;
  cout << "* ASCII file into a ROOT Ntuple  *" << endl;
  cout << "*                                *" << endl;
  cout << "**********************************" << endl;
  cout << endl << endl;

  cout << "OSCAR compliant input file: ";
  cin >> input_file;
  cout << endl;

  cout << "Output file containing ROOT Ntuple: ";
  cin >> output_file;
  cout << endl;

  cout << "Number of events to convert (0=all): ";
  cin >> max_events;
  cout << endl;

  cout << "Opening ROOT output file: "  << output_file << endl;
  TFile * hfile = new TFile(output_file,"RECREATE","ROOT file");

  TNtuple *particle = new TNtuple("particle","primary particle ntuple",
	  "event:pnum:pid:px:py:pz:E:mass:xvtx:yvtx:zvtx:theta:phi:rap:eta");

  cout << "Opening input file: " << input_file << endl;
  cout << endl;
  ifstream * input_stream = new ifstream;;
  input_stream->open(input_file);

  if ( input_stream->get(character) )
  {
    *input_stream >> line;
    if ( line!="OSC1999A" ) cout << "error in OSCAR file!" << endl;
    for ( int idummy=0; idummy<11; idummy++) *input_stream >> line;
  }

  Pi = acos(-1.);
  event_id = 0;
  do
  {
    *input_stream >> particles_per_event;
    *input_stream >> zero;
    if ( particles_per_event!=0 )
    {
      event_id++;
      for (int iparticle=0; iparticle<particles_per_event; iparticle++)
      {
        *input_stream >> pnum;
        *input_stream >> pid;
        *input_stream >> zero;
        *input_stream >> px;
        *input_stream >> py;
        *input_stream >> pz;
        *input_stream >> E;
        *input_stream >> mass;
        *input_stream >> xvtx;
        *input_stream >> yvtx;
        *input_stream >> zvtx;
        *input_stream >> opt;

	phi = atan2(py,px);
	if ( py<0.0 ) phi=phi+2.0*Pi;
	mom = sqrt(px*px+py*py+pz*pz);
        if ( mom!= 0.0 )
	{
	  theta = acos(pz/mom);
	}
	else theta = 0.0;
	eta = -log(tan(theta/2.0));
	mt  = sqrt(mass*mass+px*px+py*py);
	if ( E!=pz )
	  rap = log((E+pz)/(E-pz))/2.0;
	else
	  rap = 0.0;
	particle->Fill(event_id,pnum,pid,px,py,pz,E,mass,xvtx,yvtx,zvtx,
		       theta,phi,rap,eta);
      }
      *input_stream >> particles_per_event;
      *input_stream >> zero;
      if ( particles_per_event!=0 || zero!=0 )
        cout << "error in OSCAR file!" << endl;
    }
    else file_end=true;
    if ( event_id==max_events && max_events!=0 ) file_end=true;
  }  
  while ( !file_end );

  cout << "Saving ROOT Ntuple and closing file" << endl;
  hfile->Write();
  hfile->Close();
  hfile = 0;

  return 0;
}






