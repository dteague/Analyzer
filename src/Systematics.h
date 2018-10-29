#ifndef Systematics_hh
#define Systematics_hh

#include <string>
#include <iostream>
#include <functional>
#include <unordered_map>
#include "Particle.h"
//#include <boost/unordered_map.hpp>
#include "json/include/nlohmann/json.hpp"


// we will put stuff from the main analyser here once we know what


class TRandom3;

class Systematics {

public:
  Systematics();
  Systematics(json const &distats);
  ~Systematics();

  void init();

  void shiftParticle(Particle& jet, TLorentzVector recJet, double const& ratio, double& dPx, double& dPy, int syst);
  TLorentzVector shiftParticle(TLorentzVector recJet, double const& ratio, double& dPx, double& dPy);
  void shiftLepton(Lepton& lepton, TLorentzVector recoLep, TLorentzVector genLep, double& dPx, double& dPy, int syst);
  TLorentzVector shiftLepton(TLorentzVector recoLep, TLorentzVector genLep, double& dPx, double& dPy);
  void loadScaleRes(const json& smear, const json& syst, std::string syst_name);

private:
  double scale;
  double resolution;
};
#endif /*Systematics_hh*/
