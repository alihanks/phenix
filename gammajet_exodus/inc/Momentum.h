//-----------------------------------------------------------------------------
//
//  Declaration of the classes Mom3 and Mom4
//
//-----------------------------------------------------------------------------

#ifndef MOMENTUM_H
#define MOMENTUM_H

class Mom3
{
public:
  Mom3();
  Mom3(double initpx, double initpy, double initpz);
  ~Mom3();
  double Getpx() const {return itsPx;}
  double Getpy() const {return itsPy;}
  double Getpz() const {return itsPz;}
  void Set(double px,double py,double pz);
  void Setpx(double px){itsPx = px;}
  void Setpy(double py){itsPy = py;}
  void Setpz(double pz){itsPz = pz;}
  Mom3 operator+ (const Mom3 &);
  Mom3 operator- (const Mom3 &);
  double operator* (const Mom3 &);
  Mom3 operator* (const double &);
private:
  double itsPx, itsPy, itsPz;
};

class Mom4
{
public:
  Mom4();
  Mom4(double initE, Mom3 initp);
  Mom4(double initE, double initpx,double initpy, double initpz);
  ~Mom4();
  double GetE() const {return itsE;}
  Mom3 Getp() const {return itsP;}
  void SetE(double E){itsE=E;}
  void Setp(Mom3 p){itsP=p;}
  void Setp(double px, double py, double pz)
    {
      itsP.Setpx(px); itsP.Setpy(py); itsP.Setpz(pz);
    } 
  Mom4 operator+ (const Mom4 &);
  double operator* (const Mom4 &);
private:
  double itsE;
  Mom3   itsP;
};

#endif /* MOMENTUM_H */










