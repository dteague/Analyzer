#ifndef Analyzer_h
#define Analyzer_h

//// BIG_NUM = sqrt(sizeof(int)) so can use diparticle convention of
//// index = BIG_NUM * i1 + i2
//// This insures easy way to extract indices
//// Needs to be changed if go to size_t instead (if want to play safe
#define BIG_NUM 46340


// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>

#include <TDirectory.h>
#include <TEnv.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1.h>

#include "Particle.h"
#include "MET.h"
#include "Histo.h"

#include "json/include/nlohmann/json.hpp"

/////fix
#include "./btagging/BTagCalibrationStandalone.h"

#include "Cut_enum.h"
#include "FillInfo.h"
//#include "CRTest.h"
#include "Systematics.h"
#include "JetScaleResolution.h"
#include "DepGraph.h"
#include "LeptonPair.h"


double normPhi(double phi);
double absnormPhi(double phi);

//#define const
//using namespace std;

static const int nTrigReq = 2;

class Analyzer {
  //  friend class CRTester;
public:
  Analyzer(std::vector<std::string>, std::string, std::string configFolder="PartDet", bool useMeta = false);
  ~Analyzer();
  void add_metadata(std::vector<std::string> infiles);
  void clear_values();
  bool preprocess(int);
  bool fillCuts(bool);
  void printCuts();
  int nentries;
  void fill_histogram();
  void fill_Tree();

  ///// Functions /////

  void fill_Folder(std::string, const int, Histogramer& ihisto, bool issyst);

  void initializePileupInfo(std::string, std::string, std::string, std::string);
  void read_info(std::string);
  void setupGeneral();
  void branchException(std::string);
  void setCutNeeds();

  void smearLepton(Lepton&, CUTS, const json&, const json&, int syst=0);
  void smearJet(Particle&, CUTS, const json&, int syst=0);

  bool JetMatchesLepton(const Lepton&, const TLorentzVector&, double, CUTS);
  TLorentzVector matchLeptonToGen(const TLorentzVector&, const json&, CUTS);
  TLorentzVector matchJetToGen(const TLorentzVector&, const json&, CUTS);

  int matchToGenPdg(const TLorentzVector& lvec, double minDR);


  void getGoodParticles(int);
  void getGoodGen(const json&);
  void getGoodRecoLeptons(const Lepton&, const CUTS, const CUTS, const json&, const int);
  void getGoodRecoJets(CUTS, const json&, const int);
  void getGoodRecoFatJets(CUTS, const json&, const int);
  int getLooseLepton(const Lepton&);

  void getGoodLeptonPair(const json&, const int, const LepPair&);
  bool isZdecay(const TLorentzVector&, const TLorentzVector&);
  
  void TriggerCuts(CUTS);


  double calculateLeptonMetMt(const TLorentzVector&);
  bool isZdecay(const TLorentzVector&, const Lepton&);
  bool isOverlaping(const TLorentzVector&, Lepton&, CUTS, double);
  bool passProng(std::string, int);
  bool isInTheCracks(float);

  void create_fillInfo();

  inline bool passCutRange(double, const json&);

  void metCuts(int syst=0);
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> getArray();
  void setup_Tree();
  std::vector<std::string> bset(const json& j);

  inline int DiNum(int one, int two) { return BIG_NUM * one + two;}
  inline int poneNum(int big) { return big/BIG_NUM;}
  inline int ptwoNum(int big) { return big%BIG_NUM;}
  
  ///// values /////

  TChain* BOOM;
  TFile* infoFile;
  TFile* routfile;
  std::string filespace = "";
  double hPU[200];
  double hPU_up[200];
  double hPU_down[200];
  
  Generated* _Gen;
  Electron* _Electron;
  Muon* _Muon;
  Jet* _Jet;
  FatJet* _FatJet;
  Met* _MET;

  Histogramer histo;
  Histogramer syst_histo;

  std::unordered_map<CUTS, std::vector<int>*, EnumHash>* active_part;
  static const std::unordered_map<std::string, CUTS> cut_num;

  Systematics systematics;
  JetScaleResolution jetScaleRes;
  json genStat;

  json distats;
  std::unordered_map<std::string, FillVals*> fillInfo;
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> goodParts;
  std::vector<std::unordered_map<CUTS, std::vector<int>*, EnumHash>> syst_parts;

  std::vector<Particle*> allParticles;
  std::vector<std::string> syst_names;
  DepGraph neededCuts;

  std::vector<std::string> trigNames;
  std::vector<bool*> trig_decision;
  std::vector<int> cuts_per, cuts_cumul;

  std::unordered_map< std::string,float > zBoostTree;

  int leadIndex, maxCut, crbins=1;
  bool isData, CalculatePUSystematics, doSystematics;

  float nTruePU = 0;
  int bestVertices = 0;
  float gen_weight = 0;

  BTagCalibration btagCalib = BTagCalibration("csvv1", "Pileup/btagging.csv");
  BTagCalibrationReader btagReader = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central");

  double rho =20.;

  const static std::vector<CUTS> jetCuts;
  const static std::vector<CUTS> nonParticleCuts;
  double pu_weight, wgt, backup_wgt;
  std::unordered_map<int, GenFill*> genMaper;

  clock_t start_time;
  std::chrono::time_point<std::chrono::system_clock> start;

  /* double plus = 0; */
  /* double minus =0; */

  
};


#endif
