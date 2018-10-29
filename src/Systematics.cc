#include "Systematics.h"

using json = nlohmann::json;

Systematics::Systematics(){}

Systematics::Systematics(json const &distats){

}
Systematics::~Systematics(){}

void Systematics::init(){

}


void Systematics::shiftParticle(Particle& jet, TLorentzVector recJet, double const& ratio, double& dPx, double& dPy, int syst){

   //add the shifted part up
   dPx+=recJet.Px()*(ratio-1);
   dPy+=recJet.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recJet*=ratio;
   jet.addP4Syst(recJet, syst);
   return;
}

TLorentzVector Systematics::shiftParticle(TLorentzVector recJet, double const& ratio, double& dPx, double& dPy){

   //add the shifted part up
   dPx+=recJet.Px()*(ratio-1);
   dPy+=recJet.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recJet*=ratio;
   return recJet;
}


void Systematics::shiftLepton(Lepton& lepton, TLorentzVector recoLep, TLorentzVector genLep, double& dPx, double& dPy, int syst){
  if (genLep == TLorentzVector(0,0,0,0)) {
    lepton.addP4Syst(recoLep, syst);
    return;
  }
  double ratio = ((genLep.Pt()*scale) + (recoLep.Pt() - genLep.Pt())*resolution)/recoLep.Pt();
  //cout<<"ratio  "<<ratio<<"  "<<scale<<"  "<<resolution    <<std::endl;
   //add the shifted part up
   dPx+=recoLep.Px()*(ratio-1);
   dPy+=recoLep.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recoLep*=ratio;
   lepton.addP4Syst(recoLep, syst);
   return;
}


TLorentzVector Systematics::shiftLepton(TLorentzVector recoLep, TLorentzVector genLep, double& dPx, double& dPy){
  if (genLep == TLorentzVector(0,0,0,0)) {
    return recoLep;
  }
  double ratio = ((genLep.Pt()*scale) + (recoLep.Pt() - genLep.Pt())*resolution)/recoLep.Pt();
  //cout<<"ratio  "<<ratio<<"  "<<scale<<"  "<<resolution    <<std::endl;
   //add the shifted part up
   dPx+=recoLep.Px()*(ratio-1);
   dPy+=recoLep.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recoLep*=ratio;
   return recoLep;
}





void Systematics::loadScaleRes(const json& smear, const json& syst, std::string syst_name) {
  scale = 1;
  resolution = 1;
  if(!smear["SmearTheParticle"].empty()) {
    scale = smear["PtScaleOffset"];
    resolution = smear["PtResolutionOffset"];
  } 
  if(syst_name.find("_Res_")) {
    double res = syst["res"];
    resolution = syst_name.find("_Up") ? 1 + res : 1 - res;
    scale=1;
  } else if(syst_name.find("_Scale_")) {
    double sc = syst["scale"];
    scale = syst_name.find("_Up") ? 1+sc : 1- sc;
    resolution=1;
  }
}

