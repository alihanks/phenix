//-----------------------------------------------------------------------------
//
//  Save ROOT objects in ROOT output file and close it
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <TFile.h>

using namespace std;

TFile * CloseROOTFile(TFile * root_file)
{
  cout << "Saving ROOT objects and closing output file" << endl;

  root_file->Write();
  root_file->Close();

  root_file = 0;

  return root_file;
}
