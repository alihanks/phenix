//-----------------------------------------------------------------------------
//
//  Declaration of the classes Particle, ParticleProperty and
//  ParticleGenerator
//
//-----------------------------------------------------------------------------

#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <TH1.h>
#include <TF2.h>

using namespace std;

class Particle
{
public:
  Particle();
  Particle(int inittype, Mom4 initmom4, double initdecaysum, 
           int initGeneration);
  Particle(int inittype, double initE, 
           double initpx, double initpy, double initpz,
           double initdecaysum, int initGeneration);
  ~Particle();
  int GetID() const {return itsType;}
  Mom4 Get4mom() const {return itsMom4;}
  void SetID(int type){itsType = type;}
  void Set4mom(Mom4 mom4) {itsMom4 = mom4;}
  void Set4mom(double E, double px, double py, double pz) 
    {
      itsMom4.SetE(E);
      itsMom4.Setp(px,py,pz);
    }
  void Set4mom(double E, Mom3 p)
    {
      itsMom4.SetE(E);
      itsMom4.Setp(p);
    }
  void SetVertex(double xVTX, double yVTX, double zVTX)
    {
      itsVertex[0] = xVTX;
      itsVertex[1] = yVTX;
      itsVertex[2] = zVTX;
    }
  double GetxVertex() const {return itsVertex[0];}
  double GetyVertex() const {return itsVertex[1];}
  double GetzVertex() const {return itsVertex[2];}
  void Show(); 
  double GetDecaysum() const {return itsDecaysum;}
  void SetDecaysum(double decaysum){itsDecaysum = decaysum;}
  int GetGeneration() const {return itsGeneration;}
  void SetGeneration(int Generation){itsGeneration = Generation;}
  void AddGeneration(){itsGeneration++;}
  long double GetMass();
  void SetWeight(float weight){itsWeight = weight;}
  float GetWeight(){return itsWeight;}
  void SetAccept(int accept){itsAccept = accept;}
  int GetAccept(){return itsAccept;}
private:
  int    itsType;
  Mom4   itsMom4;
  double itsVertex[3];
  double itsDecaysum;
  int    itsGeneration;
  float  itsWeight;
  int   itsAccept;
};

class ParticleProperty
{
public:
  ParticleProperty(){}
  ~ParticleProperty(){}
  void Set(int ID, double mass, double width, int spin, int charge);
  int GetID() const {return itsID;}
  void SetID(int ID) {itsID = ID;}
  string GetName() const {return itsName;}
  void SetName(string name){itsName = name;}
  double GetMass() const {return itsMass;}
  void SetMass(double mass){itsMass = mass;}
  double GetWidth() const {return itsWidth;}
  void SetWidth(double width){itsWidth = width;}
  int GetSpin() const {return itsSpin;}
  void SetSpin(int spin){itsSpin = spin;}
  int GetCharge() const {return itsCharge;}
  void SetCharge(int charge){itsCharge = charge;}
private:
  double itsMass, itsWidth;
  string itsName;
  int    itsID, itsSpin, itsCharge;
};

class ParticleGenerator
{
public: 
  ParticleGenerator(){}
  ~ParticleGenerator(){}
  int GetID() const {return itsID;}
  void SetID(int ID) {itsID = ID;}
  double GetWeight() const {return itsWeight;}
  void SetWeight(double weight){itsWeight = weight;}
  void SetPtYFunction(TF2 *function){itsPtYFunction=function;}
  TF2* GetPtYFunction(){return itsPtYFunction;}
  void SetPtHistogram(TH1F *histogram){itsPtHistogram=histogram;}
  TH1F* GetPtHistogram(){return itsPtHistogram;}
  void SetYHistogram(TH1F *histogram){itsYHistogram=histogram;}
  TH1F* GetYHistogram(){return itsYHistogram;}
  void SetMHistogram(TH1F *histogram){itsMHistogram=histogram;}
  TH1F* GetMHistogram(){return itsMHistogram;}
private:
  int    itsID;
  double itsWeight;
  TF2*   itsPtYFunction;
  TH1F*  itsPtHistogram;
  TH1F*  itsYHistogram;
  TH1F*  itsMHistogram;
};
   
#endif /* PARTICLE_H */











