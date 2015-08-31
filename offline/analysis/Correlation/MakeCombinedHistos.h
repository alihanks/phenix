#ifndef __MAKECOMBINEDHISTOS_H__
#define __MAKECOMBINEDHISTOS_H__

#include <string>
#include <sstream>
#include <iostream>

class TH1D;
class TFile;


class MakeCombinedHistos
{
public:
  MakeCombinedHistos(const std::string fin, const std::string fout, const std::string trig_name);
  virtual ~MakeCombinedHistos(){};

private:
  double GetNTrig(TH1D* trigpt, double trigptmin, double trigptmax);
};

#endif