//-----------------------------------------------------------------------------
//
//  Open an ASCII input file and return stream
//
//-----------------------------------------------------------------------------

#include <fstream>
#include <iostream>

using namespace std;

ifstream * OpenInputFile(char * file)
{
  cout << "Opening input file: " << file << endl;
  cout << endl;

  ifstream * input_file = new ifstream;;
  input_file->open(file);

  return (input_file);
}





