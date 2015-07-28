#include <fstream>
#include <iostream>
#include <TRandom3.h>

using namespace std;

void InitializeRandom()
{
  char character;
  UInt_t seed;

  gRandom = new TRandom3;

  seed = 0;
  cout << endl;
  cout << "Reading random-generator seed from file" << endl;
  ifstream fin("seeds.txt");
  while (fin.get(character)) 
  {
    fin >> seed;
  }
  fin.close();
  cout << "seed read: " << seed << endl;
  cout << endl;
  gRandom->SetSeed(0);
  cout << "seed randomly set to: " << gRandom->GetSeed() << endl;

  return;
}











