//-----------------------------------------------------------------------------
//
//  Open ROOT output file
//
//-----------------------------------------------------------------------------

#include <iostream>

#include <TFile.h>

using namespace std;

TFile * OpenROOTFile(char *output_file)
{
  cout << "Opening ROOT output file: "  << output_file << endl;

  TFile * hfile = new TFile(output_file,"RECREATE","Demo ROOT file");

  return hfile;
}
