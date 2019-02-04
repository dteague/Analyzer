#ifndef LEPTON_PAIR_H
#define LEPTON_PAIR_H

#include "Particle.h"


struct LepPair {
  int charge1, charge2;
  int nl1, nl2;
  TLorentzVector lep1, lep2;
  bool filled = false;

  LepPair();
  LepPair(Lepton& Muon, std::vector<int>& muonN, Lepton& Electron, std::vector<int>& elecN) {
    for(auto i: muonN) {
      fillup(i, Muon);
    }
    for(auto i: elecN) {
      fillup(i, Electron);
    }
    if(lep1.Pt() < lep2.Pt()) {
      std::swap(lep1, lep2);
      std::swap(charge1, charge2);
    }

  }
  void fillup(int num, Lepton& lep) {
    if(!filled) {
      lep1 = lep.p4(num);
      charge1 = lep.charge(num);
      nl1 = num;
      filled = true;
    } else {
      lep2 = lep.p4(num);
      charge2 = lep.charge(num);
      nl2 = num;
    }
  }
};



#endif
