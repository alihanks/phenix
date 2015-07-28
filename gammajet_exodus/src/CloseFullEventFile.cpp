//-----------------------------------------------------------------------------
//
//  Close the output stream
//
//-----------------------------------------------------------------------------

#include <fstream>

using namespace std;

void CloseFullEventFile(ofstream * output_file)
{
  output_file->close();
  return ;
}





