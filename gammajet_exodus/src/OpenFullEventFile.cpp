//-----------------------------------------------------------------------------
//
//  Open an ASCII output file and return stream
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <iostream>

using namespace std;

ofstream * OpenFullEventFile(char * file)
{
  cout << "Opening output file: " << file << endl;
  cout << endl;

  ofstream * output_file = new ofstream;;
  output_file->open(file);

  *output_file << "# OSC1999A" << endl;
  *output_file << "# final_id_p_x" << endl;
  *output_file << "# EXODUS event generator in full event mode" << endl;
  *output_file << "#" << endl;

  output_file->precision(5);

  return (output_file);
}





