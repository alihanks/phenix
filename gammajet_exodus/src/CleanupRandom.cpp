#include <fstream>
#include <iostream>
#include <TRandom.h>

using namespace std;

void CleanupRandom()
{
  UInt_t seed;

  seed = gRandom->GetSeed();
  cout << endl;
  cout << "Updating random-generator seed file" << endl;
  ofstream fout("seeds.txt");
  fout << endl;
  fout << seed << endl;
  fout.close();
  cout << "seed: " << seed << endl;
  cout << endl;

  return;
}











