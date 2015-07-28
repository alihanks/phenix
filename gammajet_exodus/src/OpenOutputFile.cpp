//-----------------------------------------------------------------------------
//
//  Open an ASCII output file and return stream
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <iostream>

using namespace std;

ofstream * OpenOutputFile(char * file)
{
  cout << "Opening output file: " << file << endl;
  cout << endl;

  ofstream * output_file = new ofstream;;
  output_file->open(file);

  return (output_file);
}





