#include "Analyzer.h"
#include "Compression.h"
#include <regex>


//// Used to convert Enums to integers
#define ival(x) static_cast<int>(x)
//// BIG_NUM = sqrt(sizeof(int)) so can use diparticle convention of
//// index = BIG_NUM * i1 + i2
//// This insures easy way to extract indices
//// Needs to be changed if go to size_t instead (if want to play safe
#define BIG_NUM 46340

///// Macros defined to shorten code.  Made since lines used A LOT and repeative.  May change to inlines
///// if tests show no loss in speed
#define histAddVal2(val1, val2, name) ihisto.addVal(val1, val2, group, max, name, wgt);
#define histAddVal(val, name) ihisto.addVal(val, group, max, name, wgt);
#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

typedef std::vector<int>::iterator vec_iter;



//////////////////////////////////////////////////////////////////
///////////////////CONSTANTS DEFINITONS///////////////////////////
//////////////////////////////////////////////////////////////////

//Filespace that has all of the .in files
const std::string PUSPACE = "Pileup/";


//////////PUBLIC FUNCTIONS////////////////////

const std::vector<CUTS> Analyzer::genCuts = {
  CUTS::eGTau, CUTS::eNuTau, CUTS::eGTop,
  CUTS::eGElec, CUTS::eGMuon, CUTS::eGZ,
  CUTS::eGW, CUTS::eGHiggs, CUTS::eGJet
};

const std::vector<CUTS> Analyzer::jetCuts = {
  CUTS::eRJet1,  CUTS::eRJet2,  CUTS::eRBJet
};

const std::vector<CUTS> Analyzer::nonParticleCuts = {
  CUTS::eRVertex,CUTS::eRTrig1, CUTS::eRTrig2,
};

const std::unordered_map<std::string, CUTS> Analyzer::cut_num = {
  {"NGenTau", CUTS::eGTau},                             {"NGenTop", CUTS::eGTop},
  {"NGenElectron", CUTS::eGElec},                       {"NGenMuon", CUTS::eGMuon},
  {"NGenZ", CUTS::eGZ},                                 {"NGenW", CUTS::eGW},
  {"NGenHiggs", CUTS::eGHiggs},                         {"NGenJet", CUTS::eGJet},
  {"NRecoMuon1", CUTS::eRMuon1},                        {"NRecoMuon2", CUTS::eRMuon2},
  {"NRecoElectron1", CUTS::eRElec1},                    {"NRecoElectron2",CUTS::eRElec2},
  {"NRecoTau1", CUTS::eRTau1},                          {"NRecoTau2", CUTS::eRTau2},
  {"NRecoJet1", CUTS::eRJet1},                          {"NRecoJet2", CUTS::eRJet2},
  {"NRecoBJet", CUTS::eRBJet},
  {"NRecoTriggers1", CUTS::eRTrig1},                    {"NRecoTriggers2", CUTS::eRTrig2},
  {"NDiMuonCombinations", CUTS::eDiMuon},               {"NDiElectronCombinations", CUTS::eDiElec},
  {"NDiTauCombinations", CUTS::eDiTau},                 {"NDiJetCombinations", CUTS::eDiJet},
  {"NMuon1Tau1Combinations", CUTS::eMuon1Tau1},         {"NMuon1Tau2Combinations", CUTS::eMuon1Tau2},
  {"NMuon2Tau1Combinations", CUTS::eMuon2Tau1},         {"NMuon2Tau2Combinations", CUTS::eMuon2Tau2},
  {"NElectron1Tau1Combinations", CUTS::eElec1Tau1},     {"NElectron1Tau2Combinations", CUTS::eElec1Tau2},
  {"NElectron2Tau1Combinations", CUTS::eElec2Tau1},     {"NElectron2Tau2Combinations", CUTS::eElec2Tau2},
  {"NMuon1Electron1Combinations", CUTS::eMuon1Elec1},   {"NMuon1Electron2Combinations", CUTS::eMuon1Elec2},
  {"NMuon2Electron1Combinations", CUTS::eMuon2Elec1},   {"NMuon2Electron2Combinations", CUTS::eMuon2Elec2},
  {"NElectron1Jet1Combinations", CUTS::eElec1Jet1},     {"NElectron1Jet2Combinations", CUTS::eElec1Jet2},
  {"NElectron2Jet1Combinations", CUTS::eElec2Jet1},     {"NElectron2Jet2Combinations", CUTS::eElec2Jet2},
  {"NLeadJetCombinations", CUTS::eSusyCom},             {"METCut", CUTS::eMET},
  {"NRecoWJet", CUTS::eRWjet},                          {"NRecoVertex", CUTS::eRVertex}
};



//////////////////////////////////////////////////////
//////////////////PUBLIC FUNCTIONS////////////////////
//////////////////////////////////////////////////////

///Constructor
Analyzer::Analyzer(std::vector<std::string> infiles, std::string outfile, bool setCR, std::string configFolder) :
  /****************/ goodParts(getArray()), genName_regex(".*([A-Z][^[:space:]]+)") {

  std::cout << "setup start" << std::endl;

  routfile = new TFile(outfile.c_str(), "RECREATE", outfile.c_str(), ROOT::CompressionSettings(ROOT::kLZMA, 9));
  add_metadata(infiles); /// possibly change so routfile seen to change

  BOOM= new TChain("Events"); // TChain might not be good...?
  infoFile=0;

  for( std::string infile: infiles){
    BOOM->AddFile(infile.c_str());
  }


  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "TOTAL EVENTS: " << nentries << std::endl;

  srand(0);  // do i need to reseed this or?
  
  filespace=configFolder + "/";
  std::cout<<"setupGeneral();"<< std::endl;
  setupGeneral();

  std::cout << std::setw(2) << distats << std::endl;

  
  //isData = distats["Run"]["isData"];
  std::cout<<"reader.load(calib, BTagEntry::FLAV_B, comb);"<<std::endl;
  btagReader.load(btagCalib, BTagEntry::FLAV_B, "comb");

  
  std::cout<<"CalculatePUSystematics..."<<std::endl;
  CalculatePUSystematics = distats["Run"]["CalculatePUSystematics"];
  std::cout<<"initializePileupInfo..."<<std::endl;
  initializePileupInfo(distats["Run"]["MCHistos"], distats["Run"]["DataHistos"],distats["Run"]["DataPUHistName"],distats["Run"]["MCPUHistName"]);
  std::cout<<"syst_names.push_back(orig)"<<std::endl;
  syst_names.push_back("orig");
  std::cout<<"unordered_map<CUTS tmp;"<<std::endl;
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> tmp;
  std::cout<<"syst_parts.push_back(tmp)"<<std::endl;
  syst_parts.push_back(tmp);
  std::cout<<"if isData Systematics useSystematics..."<<std::endl;
  if(!isData && distats["Systematics"]["useSystematics"]) {
    for(auto systname : bset(distats["Systematics"])) {
      if( systname == "useSystematics")
        doSystematics= true;
      else {
        syst_names.push_back(systname);
        syst_parts.push_back(getArray());
      }
    }
  }else {
    doSystematics=false;
  }
  std::cout<<"end of that if.. now start with _Electron = new..."<<std::endl;
  _Electron = new Electron(BOOM, filespace + "Electron_info.json", syst_names);
  _Muon     = new Muon(BOOM, filespace + "Muon_info.json", syst_names);
  _Tau      = new Taus(BOOM, filespace + "Tau_info.json", syst_names);
  _Jet      = new Jet(BOOM, filespace + "Jet_info.json", syst_names);
  _FatJet   = new FatJet(BOOM, filespace + "FatJet_info.json", syst_names);
  _MET      = new Met(BOOM, "MET" , syst_names, distats["Run"]["MT2Mass"]);

  std::cout<<"---------------------------------------------------"<<std::endl;

  if(!isData) {
    std::cout<<"This is MC if not, change the flag!"<<std::endl;
    _Gen = new Generated(BOOM, filespace + "Gen_info.json", syst_names);
    allParticles= {_Gen,_Electron,_Muon,_Tau,_Jet,_FatJet};
  } else {
    std::cout<<"This is Data if not, change the flag!"<<std::endl;
    allParticles= {_Electron,_Muon,_Tau,_Jet,_FatJet};
  }
  std::cout << "particleCutMap set here" << std::endl;
  particleCutMap[CUTS::eGElec]=_Electron;
  particleCutMap[CUTS::eGMuon]=_Muon;
  particleCutMap[CUTS::eGTau]=_Tau;

  std::vector<std::string> cr_variables;
  std::cout << "Histogramer" << std::endl;
  histo = Histogramer(1, filespace+"Hist_entries.in", filespace+"Cuts.in", outfile, isData, cr_variables);
  std::cout << "systematics" << std::endl;
  if(doSystematics)
    syst_histo=Histogramer(1, filespace+"Hist_syst_entries.in", filespace+"Cuts.in", outfile, isData, cr_variables,syst_names);
  //  systematics = Systematics(distats);
  jetScaleRes = JetScaleResolution("Pileup/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt", "",  "Pileup/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt", "Pileup/Spring16_25nsV6_MC_SF_AK4PFchs.txt");
  

  std::cout << "set up ttree" << std::endl;
  ///this can be done nicer
  //put the variables that you use here:
  zBoostTree["tau1_pt"] =0;
  zBoostTree["tau1_eta"]=0;
  zBoostTree["tau1_phi"]=0;
  zBoostTree["tau2_pt"] =0;
  zBoostTree["tau2_eta"]=0;
  zBoostTree["tau2_phi"]=0;
  zBoostTree["met"]     =0;
  zBoostTree["mt_tau1"] =0;
  zBoostTree["mt_tau2"] =0;
  zBoostTree["mt2"]     =0;
  zBoostTree["cosDphi1"]=0;
  zBoostTree["cosDphi2"]=0;
  zBoostTree["jet1_pt"] =0;
  zBoostTree["jet1_eta"]=0;
  zBoostTree["jet1_phi"]=0;
  zBoostTree["jet2_pt"] =0;
  zBoostTree["jet2_eta"]=0;
  zBoostTree["jet2_phi"]=0;
  zBoostTree["jet_mass"]=0;


  histo.createTree(&zBoostTree,"TauTauTree");

  cuts_per.resize(histo.get_cuts()->size());
  cuts_cumul.resize(histo.get_cuts()->size());

  std::cout << "  create_fillInfo();" << std::endl;
  create_fillInfo();
  // for(auto maper: distats["Control_Region"].dmap) {

  //   setupCR(maper.first, maper.second);
  // }
  // check if we need to make gen level cuts to cross clean the samples:
  std::cout <<   "for(auto iselect : gen_selection)" << std::endl;
  for(auto iselect : gen_selection){
    if(iselect.second){
      std::cout<<"Waning: The selection "<< iselect.first<< " is active!"<<std::endl;
    }
  }

  std::cout << "  setCutNeeds();" << std::endl;
  setCutNeeds();

  std::cout << "setup complete" << std::endl << std::endl;
  start = std::chrono::system_clock::now();
}

/*
   Add info from nano file that isn't all of the data
   this is info that might be useful for plotting or 
   tracking info. 
 */
void Analyzer::add_metadata(std::vector<std::string> infiles){

  std::cout<<"Start copying the essentials."<<std::endl;
  for( std::string infile: infiles){
    std::cout<<infile<<std::endl;

    TFile* rfile = TFile::Open(infile.c_str());
    routfile->cd();
    for(const auto&& k: *rfile->GetListOfKeys()){
      std::string kn(k->GetName());
      std::cout<<kn<<std::endl;
      if (kn == "Events"){
        TTree* t= ((TTree*) rfile->Get(kn.c_str()));
        t->SetBranchStatus("*",0);
        t->SetBranchStatus("run",1);
        otherTrees[kn] = t->CopyTree("1","",1);
      }else if(kn == "MetaData" or kn== "ParameterSets"){
        otherTrees[kn] = ((TTree*) rfile->Get(kn.c_str()))->CopyTree("1");
      }else if(kn == "LuminosityBlocks" or kn == "Runs"){
        otherTrees[kn] = ((TTree*) rfile->Get(kn.c_str()))->CopyTree("1");
      }else if( std::string(k->ClassName()) == "TTree"){
        std::cout<<"Not copying unknown tree kn"<<std::endl;
      }else{
        //otherObjects[kn] = rfile.Get(kn)
      }
    }
    routfile->cd();
    for(auto t : otherTrees){
      t.second->Write();
    }
    rfile->Close();
    delete rfile;
  }
  
  std::cout<<"Finished copying the essentials."<<std::endl;
}

std::unordered_map<CUTS, std::vector<int>*, EnumHash> Analyzer::getArray() {
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> rmap;
  for(auto e: Enum<CUTS>()) {
    rmap[e] = new std::vector<int>();
  }
  return rmap;
}



void Analyzer::create_fillInfo() {

  fillInfo["FillLeadingJet"] = new FillVals(CUTS::eSusyCom, FILLER::Dipart, _Jet, _Jet);
  fillInfo["FillGen"] =        new FillVals(CUTS::eGen, FILLER::Single, _Gen);
  fillInfo["FillTau1"] =       new FillVals(CUTS::eRTau1, FILLER::Single, _Tau);
  fillInfo["FillTau2"] =       new FillVals(CUTS::eRTau2, FILLER::Single, _Tau);
  fillInfo["FillMuon1"] =      new FillVals(CUTS::eRMuon1, FILLER::Single, _Muon);
  fillInfo["FillMuon2"] =      new FillVals(CUTS::eRMuon2, FILLER::Single, _Muon);
  fillInfo["FillElectron1"] =  new FillVals(CUTS::eRElec1, FILLER::Single, _Electron);
  fillInfo["FillElectron2"] =  new FillVals(CUTS::eRElec2, FILLER::Single, _Electron);

  fillInfo["FillJet1"] =       new FillVals(CUTS::eRJet1, FILLER::Single, _Jet);
  fillInfo["FillJet2"] =       new FillVals(CUTS::eRJet2, FILLER::Single, _Jet);
  fillInfo["FillBJet"] =       new FillVals(CUTS::eRBJet, FILLER::Single, _Jet);
  fillInfo["FillCentralJet"] = new FillVals(CUTS::eRCenJet, FILLER::Single, _Jet);
  fillInfo["FillWJet"] =       new FillVals(CUTS::eRWjet, FILLER::Single, _FatJet);

  fillInfo["FillDiElectron"] = new FillVals(CUTS::eDiElec, FILLER::Dipart, _Electron, _Electron);
  fillInfo["FillDiMuon"] =     new FillVals(CUTS::eDiMuon, FILLER::Dipart, _Muon, _Muon);
  fillInfo["FillDiTau"] =      new FillVals(CUTS::eDiTau, FILLER::Dipart, _Tau, _Tau);
  fillInfo["FillMetCuts"] =    new FillVals();
  fillInfo["FillDiJet"] =      new FillVals(CUTS::eDiJet, FILLER::Dipart, _Jet, _Jet);

  fillInfo["FillMuon1Tau1"] =       new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon1Tau2"] =       new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau1"] =       new FillVals(CUTS::eMuon2Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau2"] =       new FillVals(CUTS::eMuon2Tau2, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillElectron1Tau1"] =   new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron1Tau2"] =   new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau1"] =   new FillVals(CUTS::eElec2Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau2"] =   new FillVals(CUTS::eElec2Tau2, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillMuon1Electron1"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon1Electron2"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron1"] =  new FillVals(CUTS::eMuon2Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron2"] =  new FillVals(CUTS::eMuon2Elec2, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillElectron1Jet1"] =   new FillVals(CUTS::eElec1Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron1Jet2"] =   new FillVals(CUTS::eElec1Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron2Jet1"] =   new FillVals(CUTS::eElec2Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron2Jet2"] =   new FillVals(CUTS::eElec2Jet2, FILLER::Dilepjet, _Electron, _Jet);

  //////I hate this solution so much.  Its terrible
  fillInfo["FillElectron1Electron2"] =     new FillVals(CUTS::eDiElec, FILLER::Single, _Electron, _Electron);
  fillInfo["FillMuon1Muon2"] =             new FillVals(CUTS::eDiMuon, FILLER::Single, _Muon, _Muon);
  fillInfo["FillTau1Tau2"] =               new FillVals(CUTS::eDiTau, FILLER::Single, _Tau, _Tau);

  //efficiency plots
  //In principal the efficiency plots should only be used, when also the object is used, but hey nobody knows!
  fillInfo["FillTauEfficiency1"] =       new FillVals(CUTS::eRTau1, FILLER::Single, _Tau);
  fillInfo["FillTauEfficiency2"] =       new FillVals(CUTS::eRTau2, FILLER::Single, _Tau);
  fillInfo["FillMuonEfficiency1"] =      new FillVals(CUTS::eRMuon1, FILLER::Single, _Muon);
  fillInfo["FillMuonEfficiency2"] =      new FillVals(CUTS::eRMuon2, FILLER::Single, _Muon);
  fillInfo["FillElectronEfficiency1"] =  new FillVals(CUTS::eRElec1, FILLER::Single, _Electron);
  fillInfo["FillElectronEfficiency2"] =  new FillVals(CUTS::eRElec2, FILLER::Single, _Electron);
  fillInfo["FillJetEfficiency1"] =       new FillVals(CUTS::eRJet1, FILLER::Single, _Jet);
  fillInfo["FillJetEfficiency2"] =       new FillVals(CUTS::eRJet2, FILLER::Single, _Jet);



  for(auto it: *histo.get_groups()) {
    if(fillInfo[it] == nullptr) fillInfo[it] = new FillVals();
  }

}


////destructor
Analyzer::~Analyzer() {
  clear_values();
  delete BOOM;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
  if(!isData) delete _Gen;

  for(auto fpair: fillInfo) {
    delete fpair.second;
    fpair.second=nullptr;
  }

  for(auto e: Enum<CUTS>()) {
    delete goodParts[e];
    goodParts[e]=nullptr;
  }

}


///resets values so analysis can start
void Analyzer::clear_values() {

  for(auto e: Enum<CUTS>()) {
    goodParts[e]->clear();
  }
  //faster!!
  for(auto &it: syst_parts) {
    if (it.size() == 0) continue;
    for(auto e: Enum<CUTS>()) {
      it[e]->clear();
    }
  }
  if(infoFile!=BOOM->GetFile()){
    std::cout<<"New file!"<<std::endl;
    infoFile=BOOM->GetFile();
  }

  leadIndex=-1;
  maxCut = 0;
}

///Function that does most of the work.  Calculates the number of each particle
bool Analyzer::preprocess(int event) {
  int test= BOOM->GetEntry(event);
  if(test<0){
    std::cout << "Could not read the event from the following file: "<<BOOM->GetFile()->GetNewUrl().Data() << std::endl;
  }
  for(Particle* ipart: allParticles){
    ipart->init();
  }
  _MET->init();

  active_part = &goodParts;

  if(nTruePU > 199) return false;
  pu_weight = (!isData && CalculatePUSystematics) ? hPU[(int)(nTruePU+1)] : 1.0;

  
  // SET NUMBER OF GEN PARTICLES
  if(!isData){
    _Gen->setOrigReco();
    getGoodGen(_Gen->pstats["Gen"]);
    getGoodTauNu();
  }


  //////Triggers and Vertices
  active_part->at(CUTS::eRVertex)->resize(bestVertices);
  //TriggerCuts(*(trigPlace[0]), *(trigName[0]), CUTS::eRTrig1);
  TriggerCuts( CUTS::eRTrig1);
  ////check update met is ok
  for(size_t i=0; i < syst_names.size(); i++) {
     //////Smearing
    smearLepton(*_Electron, CUTS::eGElec, _Electron->pstats["Smear"], distats["Electron_systematics"], i);
    smearLepton(*_Muon, CUTS::eGMuon, _Muon->pstats["Smear"], distats["Muon_systematics"], i);
    smearLepton(*_Tau, CUTS::eGTau, _Tau->pstats["Smear"], distats["Tau_systematics"], i);

    smearJet(*_Jet,CUTS::eGJet,_Jet->pstats["Smear"], i);
    smearJet(*_FatJet,CUTS::eGJet,_FatJet->pstats["Smear"], i);
    updateMet(i);

  }
  for(size_t i=0; i < syst_names.size(); i++) {
    std::string systname = syst_names.at(i);
    for( auto part: allParticles) part->setCurrentP(i);
    _MET->setCurrentP(i);
    getGoodParticles(i);
  }
  active_part = &goodParts;

  if( event < 10 || ( event < 100 && event % 10 == 0 ) ||
    ( event < 1000 && event % 100 == 0 ) ||
    ( event < 10000 && event % 1000 == 0 ) ||
    ( event >= 10000 && event % 10000 == 0 ) ) {
       std::cout << std::setprecision(2)<<event << " Events analyzed "<< static_cast<double>(event)/nentries*100. <<"% done"<<std::endl;
       std::cout << std::setprecision(5);
  }

  return true;
}


void Analyzer::getGoodParticles(int syst){

  std::string systname=syst_names.at(syst);
  if(syst == 0) active_part = &goodParts;
  else active_part=&syst_parts.at(syst);
    //    syst=syst_names[syst];



  // // SET NUMBER OF RECO PARTICLES
  // // MUST BE IN ORDER: Muon/Electron, Tau, Jet
  getGoodRecoLeptons(*_Electron, CUTS::eRElec1, CUTS::eGElec, _Electron->pstats["Elec1"],syst);
  getGoodRecoLeptons(*_Electron, CUTS::eRElec2, CUTS::eGElec, _Electron->pstats["Elec2"],syst);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon1, CUTS::eGMuon, _Muon->pstats["Muon1"],syst);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon2, CUTS::eGMuon, _Muon->pstats["Muon2"],syst);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau1, CUTS::eGTau, _Tau->pstats["Tau1"],syst);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau2, CUTS::eGTau, _Tau->pstats["Tau2"],syst);
  getGoodRecoJets(CUTS::eRBJet, _Jet->pstats["BJet"],syst);
  getGoodRecoJets(CUTS::eRJet1, _Jet->pstats["Jet1"],syst);
  getGoodRecoJets(CUTS::eRJet2, _Jet->pstats["Jet2"],syst);
  getGoodRecoFatJets(CUTS::eRWjet, _FatJet->pstats["WJet"],syst);

}



////Reads cuts from Cuts.in file and see if the event has enough particles
bool Analyzer::fillCuts(bool fillCounter) {
  const std::unordered_map<std::string,std::pair<int,int> >* cut_info = histo.get_cuts();
  const std::vector<std::string>* cut_order = histo.get_cutorder();

  bool prevTrue = true;

  maxCut=0;
  //  std::cout << active_part << std::endl;;

  for(size_t i = 0; i < cut_order->size(); i++) {
    std::string cut = cut_order->at(i);
    if(isData && cut.find("Gen") != std::string::npos){
      maxCut += 1;
      continue;
    }
    int min= cut_info->at(cut).first;
    int max= cut_info->at(cut).second;
    int nparticles = active_part->at(cut_num.at(cut))->size();
    //if(!fillCounter) std::cout << cut << ": " << nparticles << " (" << min << ", " << max << ")" <<std::endl;
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      if(fillCounter && crbins == 1) {
        cuts_per[i]++;
        cuts_cumul[i] += (prevTrue) ? 1 : 0;
        maxCut += (prevTrue) ? 1 : 0;
      }else{
        maxCut += (prevTrue) ? 1 : 0;
      }
    }else {
      //cout<<"here 2  "<<std::endl;
      prevTrue = false;
    }
  }

  return prevTrue;
}



///Prints the number of events that passed each cut per event and cumulatively
//done at the end of the analysis
void Analyzer::printCuts() {
  std::vector<std::string> cut_order;
  if(crbins > 1) cut_order = *(histo.get_folders());
  else cut_order = *(histo.get_cutorder());
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  double run_time_real=elapsed_seconds.count();


  std::cout.setf(std::ios::floatfield,std::ios::fixed);
  std::cout<<std::setprecision(3);
  std::cout << "\n";
  std::cout << "Selection Efficiency " << "\n";
  std::cout << "Total events: " << nentries << "\n";
  std::cout << "\n";
  std::cout << "Run Time (real): " <<run_time_real <<" s\n";
  std::cout << "Time per 1k Events (real): " << run_time_real/(nentries/1000) <<" s\n";
  std::cout << "Events/s: " << static_cast<double>(nentries)/(run_time_real) <<" 1/s (real) \n";
  std::cout << "                        Name                  Indiv.";
  if(crbins == 1) std::cout << "            Cumulative";
  std::cout << std::endl << "---------------------------------------------------------------------------\n";
  for(size_t i = 0; i < cut_order.size(); i++) {
    std::cout << std::setw(28) << cut_order.at(i) << "    ";
    if(isData && cut_order.at(i).find("Gen") != std::string::npos) std::cout << "Skipped" << std::endl;
    else if(crbins != 1 && blinded && i == (size_t)SignalRegion) std::cout << "Blinded Signal Region" << std::endl;
    else {
      std::cout << std::setw(10) << cuts_per.at(i) << "  ( " << std::setw(5) << ((float)cuts_per.at(i)) / nentries << ") ";
      if(crbins == 1) std::cout << std::setw(12) << cuts_cumul.at(i) << "  ( " << std::setw(5) << ((float)cuts_cumul.at(i)) / nentries << ") ";

      std::cout << std::endl;
    }
  }
  std::cout <<std::setprecision(5);
  std::cout << "---------------------------------------------------------------------------\n";

  //write all the histograms
  //attention this is not the fill_histogram method from the Analyser
  histo.fill_histogram(routfile);
  if(doSystematics)
    syst_histo.fill_histogram(routfile);

}

/////////////PRIVATE FUNCTIONS////////////////


///Calculates met from values from each file plus smearing and treating muons as neutrinos
void Analyzer::updateMet(int syst) {
  // _MET->update(distats["Run"], *_Jet,  syst);

  /////MET CUTS

  if(!passCutRange(_MET->pt(), distats["Run"]["MetCut"])) return;
  if(distats["Run"]["DiscrByHT"] && _MET->HT() < distats["Run"]["HtCut"]) return;

  if(syst==0){
    active_part->at(CUTS::eMET)->push_back(1);
  }else{
    syst_parts.at(syst).at(CUTS::eMET)->push_back(1);
  }
}

/////Check if a given branch is not found in the file

void Analyzer::branchException(std::string branch){
  if(BOOM->FindBranch(branch.c_str()) == 0 ){
     throw "Branch not found in the current sample. Check the config files associated to this branch.";
  }
}

/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral() {

  genMaper = {
    {5, new GenFill(2, CUTS::eGJet)},     {6,  new GenFill(2, CUTS::eGTop)},
    {11, new GenFill(1, CUTS::eGElec)},   {13, new GenFill(1, CUTS::eGMuon)},
    {15, new GenFill(2, CUTS::eGTau)},    {23, new GenFill(62, CUTS::eGZ)},
    {24, new GenFill(62, CUTS::eGW)},      {25, new GenFill(2, CUTS::eGHiggs)}
};
  
  isData=true;
  if(BOOM->FindBranch("Pileup_nTrueInt")!=0){
    isData=false;
  }
  if(!isData){
    SetBranch("Pileup_nTrueInt", nTruePU);
    SetBranch("genWeight", gen_weight);
    //SetBranch("rho", rho);
  }else{
    nTruePU=0;
    gen_weight=0;
  }

  SetBranch("PV_npvs", bestVertices);

  read_info(filespace + "Run_info.json");
  read_info(filespace + "Systematics_info.json");
  
  for(std::string trigger : trigNames){
    bool decison=false;
    //std::cout<< "Trigger: "<< trigger<<std::endl;
    
   for( int i=0; i<BOOM->GetListOfBranches()->GetSize(); i++){
	   std::string branch_name(BOOM->GetListOfBranches()->At(i)->GetName());
	   if (branch_name.find("HLT_")!=std::string::npos){
           	std::cout<< branch_name << std::endl;
	   }
   }
   try{
     branchException(trigger.c_str());
   }

   catch (const char* msg){
     std::cout << "ERROR! Trigger " << trigger << ": "  << msg << std::endl;
     std::cout<< "options are:" << std::endl;
     for( int i=0; i<BOOM->GetListOfBranches()->GetSize(); i++){
       std::string branch_name(BOOM->GetListOfBranches()->At(i)->GetName());
       if (branch_name.find("HLT_")!=std::string::npos){
	 std::cout<< branch_name << std::endl;
       }
     }
     continue;
     //std::exit(1);
   }

   SetBranch(trigger.c_str(),decison);       
   trig_decision.push_back(&decison);
    
  }
}



///parsing method that gets info on diparts and basic run info
//put in std::map called "distats"
void Analyzer::read_info(std::string filename) {
  std::ifstream info_file(filename);
  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }

  std::string line;
  std::stringstream ss;
  
  while(getline(info_file, line)) {
    std::size_t found = line.find("//");
    if(found != std::string::npos) {
      line = line.substr(0, found);
    }
    ss << line;
  }

  json jtemp;
  
  ss >> jtemp;

  distats.update(jtemp);
  
  info_file.close();
}


// This code works pretty much (at least in my tests), but dagnabit, its ugly.  They all can't be winners, at least now...
void Analyzer::setCutNeeds() {


  for(auto it: *histo.get_groups()) {
    if(fillInfo[it]->type == FILLER::None) continue;
    neededCuts.loadCuts(fillInfo[it]->ePos);
  }
  for(auto it : *histo.get_cutorder()) {
    try{
      neededCuts.loadCuts(cut_num.at(it));
    }catch(...){
      std::cout<<"The following cut is strange: "<<it<<std::endl;
      exit(2);
    }
  }
  // for(auto it: testVec) {
  //   neededCuts.loadCuts(it->info->ePos);
  // }
  if(!isData and distats["Run"]["ApplyZBoostSF"] and isVSample){
    neededCuts.loadCuts(CUTS::eGen);
    neededCuts.loadCuts(CUTS::eGZ);
    neededCuts.loadCuts(CUTS::eGW);
  }

  neededCuts.loadCuts(_Jet->findExtraCuts());
  if(doSystematics) {
    neededCuts.loadCuts(CUTS::eGen);
  }

  for(auto it: jetCuts) {
    if(!neededCuts.isPresent(it)) continue;
    neededCuts.loadCuts(_Jet->overlapCuts(it));
  }

  if(neededCuts.isPresent(CUTS::eRWjet)) {
    neededCuts.loadCuts(_FatJet->findExtraCuts());
    neededCuts.loadCuts(_FatJet->overlapCuts(CUTS::eRWjet));
  } else {
    std::cout<<"WJets not needed. They will be deactivated!"<<std::endl;
    _FatJet->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRTau1) || neededCuts.isPresent(CUTS::eRTau2) ) {
    neededCuts.loadCuts(_Tau->findExtraCuts());
  } else {
    std::cout<<"Taus not needed. They will be deactivated!"<<std::endl;
    _Tau->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRElec1) || neededCuts.isPresent(CUTS::eRElec2) ) {
    neededCuts.loadCuts(_Electron->findExtraCuts());
  } else {
    std::cout<<"Electrons not needed. They will be deactivated!"<<std::endl;
    _Electron->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRMuon1) || neededCuts.isPresent(CUTS::eRMuon2) ) {
    neededCuts.loadCuts(_Muon->findExtraCuts());
  } else {
    std::cout<<"Muons not needed. They will be deactivated!"<<std::endl;
    _Muon->unBranch();
  }

  if( !neededCuts.isPresent(CUTS::eGen) and !isData) {
    std::cout<<"Gen not needed. They will be deactivated!"<<std::endl;
    _Gen->unBranch();

  }

  std::cout << "Cuts being filled: " << std::endl;
  for(auto cut : neededCuts.getCuts()) {
    std::cout << enumNames.at(static_cast<CUTS>(cut)) << "   ";
  }
  std::cout << std::endl;
}


///Smears lepton only if specified and not a data file.  Otherwise, just filles up lorentz std::vectors
//of the data into the std::vector container smearP with is in each lepton object.
void Analyzer::smearLepton(Lepton& lep, CUTS eGenPos, const json& stats, const json& syst_stats, int syst) {
  if( isData) {
    lep.setOrigReco();
    return;
  }

  std::string systname = syst_names.at(syst);
  if(!lep.needSyst(syst)) return;

  if(systname=="orig" && !stats["SmearTheParticle"]){
    lep.setOrigReco();
  } else {
    systematics.loadScaleRes(stats, syst_stats, systname);
    for(size_t i = 0; i < lep.size(); i++) {
      TLorentzVector lepReco = lep.RecoP4(i);
      TLorentzVector genVec =  matchLeptonToGen(lepReco, lep.pstats["Smear"],eGenPos);
      TLorentzVector newVec = systematics.shiftLepton(lepReco, genVec, _MET->systdeltaMEx[syst], _MET->systdeltaMEy[syst]);
      lep.addP4Syst(newVec, syst);
    }
  }
}

///Same as smearlepton, just jet specific
void Analyzer::smearJet(Particle& jet, const CUTS eGenPos, const json& stats, int syst) {
  //at the moment
  if(isData || jet.type != PType::Jet ){
    //|| !stats["SmearTheJet"]
    jet.setOrigReco();
    return;
  }
  if(!jet.needSyst(syst)){
    return;
  }
  //add energy scale uncertainty

  std::string systname = syst_names.at(syst);

  for(size_t i=0; i< jet.size(); i++) {
    TLorentzVector jetReco = jet.RecoP4(i);
    jet.addP4Syst(jetReco,syst);
    continue;
    // if(JetMatchesLepton(*_Muon, jetReco, stats["MuonMatchingDeltaR"], CUTS::eGMuon) ||
    //    JetMatchesLepton(*_Tau, jetReco, stats["TauMatchingDeltaR"], CUTS::eGTau) ||
    //    JetMatchesLepton(*_Electron, jetReco,stats["ElectronMatchingDeltaR"], CUTS::eGElec)){
    //   jet.addP4Syst(jetReco,syst);
    //   continue;
    // }
    double sf=1.;
    //only apply corrections for jets not for FatJets

    TLorentzVector genJet=matchJetToGen(jetReco, jet.pstats["Smear"],eGenPos);
    if(systname=="orig" && stats["SmearTheJet"]){
      sf=jetScaleRes.GetRes(jetReco,genJet, rho, 0);
    }else if(systname=="Jet_Res_Up"){
      sf=jetScaleRes.GetRes(jetReco,genJet, rho, 1);
    }else if(systname=="Jet_Res_Down"){
      sf=jetScaleRes.GetRes(jetReco,genJet, rho, -1);
    }else if(systname=="Jet_Scale_Up"){
      sf = jetScaleRes.GetScale(jetReco, false, +1.);
    }else if(systname=="Jet_Scale_Down"){
      sf = jetScaleRes.GetScale(jetReco, false, -1) ;
    }
    //    std::cout<<systname<<"  "<<sf<<"  "<<jetReco.Pt()<<"  "<<genJet.Pt()<<std::endl;
    systematics.shiftParticle(jet, jetReco, sf, _MET->systdeltaMEx[syst], _MET->systdeltaMEy[syst], syst);
  }
}


/////checks if jet is close to a lepton and the lepton is a gen particle, then the jet is a lepton object, so
//this jet isn't smeared
bool Analyzer::JetMatchesLepton(const Lepton& lepton, const TLorentzVector& jetV, double partDeltaR, CUTS eGenPos) {
  for(size_t j = 0; j < lepton.size(); j++) {
    if(jetV.DeltaR(lepton.RecoP4(j)) < partDeltaR && matchLeptonToGen(lepton.RecoP4(j), lepton.pstats.at("Smear"), eGenPos) != TLorentzVector(0,0,0,0)) return true;
  }
  return false;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchLeptonToGen(const TLorentzVector& lvec, const json& stats, CUTS ePos) {
  if(ePos == CUTS::eGTau) {
    return matchTauToGen(lvec, stats["GenMatchingDeltaR"]);
  }
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(_Gen->p4(it)) <= stats["GenMatchingDeltaR"]) {
      if(stats["UseMotherID"] && abs(_Gen->pdg_id[_Gen->genPartIdxMother[it]]) != stats["MotherID"]) continue;
      return _Gen->p4(it);
    }
  }
  return TLorentzVector(0,0,0,0);
}


///Tau specific matching function.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically
TLorentzVector Analyzer::matchTauToGen(const TLorentzVector& lvec, double lDeltaR) {
  TLorentzVector genVec(0,0,0,0);
  int i = 0;
  for(vec_iter it=active_part->at(CUTS::eGTau)->begin(); it !=active_part->at(CUTS::eGTau)->end();it++, i++) {
    int nu = active_part->at(CUTS::eNuTau)->at(i);
    if(nu == -1) continue;

    genVec = _Gen->p4(*it) - _Gen->p4(nu);
    if(lvec.DeltaR(genVec) <= lDeltaR) {
      return genVec;
    }
  }
  return genVec;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchJetToGen(const TLorentzVector& lvec, const json& stats, CUTS ePos) {
  //for the future store gen jets
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(_Gen->p4(it)) <= stats["GenMatchingDeltaR"]) {
      //nothing more than b quark or gluon
      if( !( (abs(_Gen->pdg_id[it])<5) || (abs(_Gen->pdg_id[it])==9) ||  (abs(_Gen->pdg_id[it])==21) ) ) continue;
      return _Gen->p4(it);
    }
  }
  return TLorentzVector(0,0,0,0);
}



////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
int Analyzer::matchToGenPdg(const TLorentzVector& lvec, double minDR) {
  double _minDR=minDR;
  int found=-1;
  for(size_t i=0; i< _Gen->size(); i++) {

    if(lvec.DeltaR(_Gen->p4(i)) <=_minDR) {
      //only hard interaction
      if( _Gen->status[i]<10){
        found=i;
        _minDR=lvec.DeltaR(_Gen->p4(i));
      }
    }
  }
  if (found>=0){
    return _Gen->pdg_id[found];
  }
  return 0;
}


////Calculates the number of gen particles.  Based on id number and status of each particle
void Analyzer::getGoodGen(const json& stats) {
  if(! neededCuts.isPresent(CUTS::eGen)) return;
  for(size_t j = 0; j < _Gen->size(); j++) {
    int id = abs(_Gen->pdg_id[j]);
    if(genMaper.find(id) != genMaper.end() && _Gen->status[j] == genMaper.at(id)->status) {
      active_part->at(genMaper.at(id)->ePos)->push_back(j);
    }
    //something special for jet
    if( (id<5 || id==9 ||  id==21) && genMaper.find(id) != genMaper.end() && _Gen->status[j] == genMaper.at(5)->status) {
      active_part->at(genMaper.at(5)->ePos)->push_back(j);
    }
  }
}


////// need to fix
////Tau neutrino specific function used for calculating the number of hadronic taus
void Analyzer::getGoodTauNu() {
  // for(auto it : *active_part->at(CUTS::eGTau)) {
  //   bool leptonDecay = false;
  //   int nu = -1;
  //   for(size_t j = 0; j < _Gen->size(); j++) {
  //     //if(abs(_Gen->BmotherIndex->at(j)) == (it)) {
  //       //if( (abs(_Gen->pdg_id->at(j)) == 16) && (abs(_Gen->motherpdg_id->at(j)) == 15) && (_Gen->status->at(_Gen->BmotherIndex->at(j)) == 2) ) nu = j;
  //       if( (abs(_Gen->pdg_id[j]) == 16) && (abs(_Gen->pdg_id[_Gen->genPartIdxMother[j]]) == 15) ) nu = j;
  //       else if( (abs(_Gen->pdg_id[j]) == 12) || (abs(_Gen->pdg_id[j]) == 14) ) leptonDecay = true;
  //     //}
  //   }
  //   nu = (leptonDecay) ? -1 : nu;
  //   active_part->at(CUTS::eNuTau)->push_back(nu);
  // }

}


////////////DIRTY, YUCK YUCK YUCK
std::vector<std::string> Analyzer::bset(const json& j) {
  std::vector<std::string> output;
// for(auto& val: j) {
  for (auto it = j.cbegin(); it != j.cend(); ++it) {
    if(it.value().is_boolean()) {
      output.push_back(it.key());
    }
  }
   
  return output;
}



///////Only tested for if is Zdecay, can include massptasymmpair later?
/////Tests to see if a light lepton came form a zdecay
bool Analyzer::isZdecay(const TLorentzVector& theObject, const Lepton& lep) {
  bool eventIsZdecay = false;
  const float zMass = 90.1876;
  const float zWidth = 2.4952;
  float zmmPtAsymmetry = -10.;

  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.
  for(std::vector<TLorentzVector>::const_iterator lepit= lep.begin(); lepit != lep.end(); lepit++) {
    if(theObject.DeltaR(*lepit) < 0.3) continue;
    if(theObject == (*lepit)) continue;

    TLorentzVector The_LorentzVect = theObject + (*lepit);
    zmmPtAsymmetry = (theObject.Pt() - lepit->Pt()) / (theObject.Pt() + lepit->Pt());

    if( (abs(The_LorentzVect.M() - zMass) < 3.*zWidth) || (fabs(zmmPtAsymmetry) < 0.20) ) {
      eventIsZdecay = true;
      break;
    }
  }

  return eventIsZdecay;
}


///Function used to find the number of reco leptons that pass the various cuts.
///Divided into if blocks for the different lepton requirements.
void Analyzer::getGoodRecoLeptons(const Lepton& lep, const CUTS ePos, const CUTS eGenPos, const json& stats, const int syst) {
  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!lep.needSyst(syst)) {
    active_part->at(ePos) = goodParts[ePos];
    return;
  }

  int i = 0;
  for(auto lvec: lep) {
    bool passCuts = true;
    if (passCuts) passCuts = fabs(lvec.Eta()) > stats["EtaCut"];
    if (passCuts) passCuts = passCutRange(lvec.Pt(), stats["PtCut"]);

    if((lep.pstats.at("Smear")["MatchToGen"]) && (!isData)) {   /////check
      if(matchLeptonToGen(lvec, lep.pstats.at("Smear") ,eGenPos) == TLorentzVector(0,0,0,0)) continue;
    }

    for( auto cut: bset(stats)) {
      if(!passCuts) break;
      // else if(cut == "DoDiscrByIsolation") {
      //   double firstIso = (! stats["IsoSumPtCutValue"].empty()) ? stats.pmap.at("IsoSumPtCutValue").first : ival(ePos) - ival(CUTS::eRTau1) + 1;
      //   double secondIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").second : stats["FlipIsolationRequirement"];
      //   passCuts = passCuts && lep.get_Iso(i, firstIso, secondIso);
      // }
      else if(cut == "DiscrIfIsZdecay" && lep.type != PType::Tau ) passCuts = isZdecay(lvec, lep);
      else if(cut == "DiscrByMetDphi") passCuts = passCutRange(absnormPhi(lvec.Phi() - _MET->phi()), stats["MetDphiCut"]);
      else if(cut == "DiscrByMetMt") passCuts = passCutRange(calculateLeptonMetMt(lvec), stats["MetMtCut"]);
      /////muon cuts
      else if(lep.type == PType::Muon){
        if(cut == "DoDiscrByTightID") passCuts = _Muon->tight[i];
        else if(cut == "DoDiscrBySoftID") passCuts = _Muon->soft[i];
      }
      ////electron cuts
      else if(lep.type == PType::Electron){
        if(cut == "DoDiscrByHLTID"){
          std::bitset<8> idvariable(_Electron->cutBased_HLTPreSel[i]);
          if(ival(ePos) - ival(CUTS::eRElec1)){
            passCuts = (_Electron->cbHLTIDele1&idvariable).count();
          }else{
            passCuts = (_Electron->cbHLTIDele2&idvariable).count();
          }
        }
        if(cut == "DoDiscrByCBID"){
          std::bitset<8> idvariable(_Electron->cutBased[i]);
          if(ival(ePos) - ival(CUTS::eRElec1)){
            passCuts = (_Electron->cbIDele1&idvariable).count();
          }else{
            passCuts = (_Electron->cbIDele2&idvariable).count();
          }
        }
        else if(cut == "DoDiscrByHEEPID")
         passCuts = _Electron->isPassHEEPId[i];
      }
      else if(lep.type == PType::Tau){
        if(cut == "DoDiscrByCrackCut") passCuts = !isInTheCracks(lvec.Eta());
        /////tau cuts
        else if(cut == "DoDzCut") passCuts = (_Tau->dz[i] <= stats["DzCutThreshold"]);
        else if(cut == "DoDiscrByLeadTrack") passCuts = (_Tau->leadTkPtOverTauPt[i]*_Tau->pt(i) >= stats["LeadTrackThreshold"]);
             // ----Electron and Muon vetos
        else if (cut == "DoDiscrAgainstElectron") passCuts = _Tau->pass_against_Elec(ePos, i);
        else if (cut == "SelectTausThatAreElectrons") passCuts = !_Tau->pass_against_Elec(ePos, i);

        else if (cut == "DoDiscrAgainstMuon") passCuts = _Tau->pass_against_Muon(ePos, i);
        else if (cut == "SelectTausThatAreMuons") passCuts = !_Tau->pass_against_Muon(ePos, i);

        else if(cut == "DiscrByProngType") {
	  std::string prongName = stats["ProngType"];
          passCuts = (prongName.find("hps") == std::string::npos || _Tau->DecayModeNewDMs[i] != 0);
          passCuts = passProng(stats["ProngType"], _Tau->decayMode[i]);
        }
        else if(cut == "decayModeFindingNewDMs") passCuts = _Tau->DecayModeNewDMs[i] != 0;
        else if(cut == "decayModeFinding") passCuts = _Tau->DecayMode[i] != 0;
              // ----anti-overlap requirements
        else if(cut == "RemoveOverlapWithMuon1s") passCuts = !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats["Muon1MatchingDeltaR"]);
        else if(cut == "RemoveOverlapWithMuon2s") passCuts = !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats["Muon2MatchingDeltaR"]);
        else if(cut == "RemoveOverlapWithElectron1s") passCuts = !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats["Electron1MatchingDeltaR"]);
        else if(cut == "RemoveOverlapWithElectron2s") passCuts = !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats["Electron2MatchingDeltaR"]);
      }
      else std::cout << "cut: " << cut << " not listed" << std::endl;
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }

  return;
}

////Jet specific function for finding the number of jets that pass the cuts.
//used to find the nubmer of good jet1, jet2, central jet, 1st and 2nd leading jets and bjet.
void Analyzer::getGoodRecoJets(CUTS ePos, const json& stats, const int syst) {

  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_Jet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;

  for(auto lvec: *_Jet) {
    bool passCuts = true;
    passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats["EtaCut"]);
    passCuts = passCuts && (lvec.Pt() > stats["PtCut"]) ;

    for( auto cut: bset(stats)) {
      if(!passCuts) break;

    /// BJet specific
      else if(cut == "ApplyJetBTagging") passCuts = passCuts && (_Jet->bDiscriminator[i] > stats["JetBTaggingCut"]);
      //else if(cut == "MatchBToGen") passCuts = passCuts && (isData ||  abs(_Jet->partonFlavour->at(i)) == 5);
      else if(cut == "ApplyLooseID") passCuts = passCuts && _Jet->passedLooseJetID(i);

    // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats["Muon1MatchingDeltaR"]);
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats["Muon2MatchingDeltaR"]);
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats["Electron1MatchingDeltaR"]);
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats["Electron2MatchingDeltaR"]);
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats["Tau1MatchingDeltaR"]);
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats["Tau2MatchingDeltaR"]);

      else if(cut == "UseBtagSF") {
        //double bjet_SF = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, lvec.Eta(), lvec.Pt());
        //passCuts = passCuts && (isData || ((double) rand()/(RAND_MAX)) <  bjet_SF);
      }
    }
    if(_Jet->pstats["BJet"]["RemoveBJetsFromJets"] and ePos!=CUTS::eRBJet){
      passCuts = passCuts && find(active_part->at(CUTS::eRBJet)->begin(), active_part->at(CUTS::eRBJet)->end(), i) == active_part->at(CUTS::eRBJet)->end();
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;

  }
}


////FatJet specific function for finding the number of V-jets that pass the cuts.
void Analyzer::getGoodRecoFatJets(CUTS ePos, const json& stats, const int syst) {
  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_FatJet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;

  for(auto lvec: *_FatJet) {
    bool passCuts = true;
    passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats["EtaCut"]);
    passCuts = passCuts && (lvec.Pt() > stats["PtCut"]) ;

    ///if else loop for central jet requirements
    for( auto cut: bset(stats)) {
      if(!passCuts) break;

      else if(cut == "ApplyJetWTagging") passCuts = passCuts && (passCutRange(_FatJet->tau2[i]/_FatJet->tau1[i], stats["JetTau2Tau1Ratio"]) &&
								 passCutRange(_FatJet->PrunedMass[i], stats["JetWmassCut"]));
    // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats["Muon1MatchingDeltaR"]);
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats["Muon2MatchingDeltaR"]);
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats["Electron1MatchingDeltaR"]);
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats["Electron2MatchingDeltaR"]);
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats["Tau1MatchingDeltaR"]);
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats["Tau2MatchingDeltaR"]);

    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }
}

///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(overlapper.p4(it)) < MatchingDeltaR) return true;
  }
  return false;
}

///Tests if tau decays into the specified number of jet prongs.
bool Analyzer::passProng(std::string prong, int value) {
  return ( (prong.find("1") != std::string::npos &&  (value<5)) ||
  (prong.find("2") != std::string::npos &&  (value>=5 && value<10)) ||
  (prong.find("3") != std::string::npos && (value>=10 && value<12)) );
}


////Tests if tau is within the cracks of the detector (the specified eta ranges)
bool Analyzer::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 ||
  (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
  (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
  (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
  (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}


///sees if the event passed one of the two cuts provided
void Analyzer::TriggerCuts(CUTS ePos) {
  if(! neededCuts.isPresent(ePos)) return;
  for(bool* trigger : trig_decision){
   //#std::cout<< "trig_decision: "<< *trigger << std::endl;
    if(*trigger){
      active_part->at(ePos)->push_back(0);
      return;
    }
  }
}


bool Analyzer::passCutRange(double value, const json& cuts) {
  return (value >= cuts.at(0) && value < cuts.at(1));
}

//-----Calculate lepton+met transverse mass
double Analyzer::calculateLeptonMetMt(const TLorentzVector& Tobj) {
  double px = Tobj.Px() + _MET->px();
  double py = Tobj.Py() + _MET->py();
  double et = Tobj.Et() + _MET->energy();
  double mt2 = et*et - (px*px + py*py);
  return (mt2 >= 0) ? sqrt(mt2) : -1;
}





////Grabs a list of the groups of histograms to be filled and asked Fill_folder to fill up the histograms
void Analyzer::fill_histogram() {
  
  if(!isData && distats["Run"]["ApplyGenWeight"] && gen_weight == 0.0) return;

  if(isData && blinded && maxCut == SignalRegion) return;

  const std::vector<std::string>* groups = histo.get_groups();
  wgt = 1.;
  if(!isData){
    if(distats["Run"]["UsePileUpWeight"]) wgt*= pu_weight;
    if(distats["Run"]["ApplyGenWeight"]) wgt *= (gen_weight > 0) ? 1.0 : -1.0;
    //add weight here
  }
  //backup current weight
  backup_wgt=wgt;

  for(size_t i = 0; i < syst_names.size(); i++) {
    for(Particle* ipart: allParticles) ipart->setCurrentP(i);
    _MET->setCurrentP(i);
    active_part =&syst_parts.at(i);
    //////i == 0 is orig or no syst case
    if(i == 0) {
      active_part = &goodParts;
      fillCuts(true);
      for(auto it: *groups) {
        fill_Folder(it, maxCut, histo, false);
      }
      if(!fillCuts(false)) {
        fill_Tree();
      }
    }else{

      wgt=backup_wgt;
      if(syst_names[i].find("weight")!=std::string::npos){
        if(syst_names[i]=="Pileup_weight_Up"){
          if(distats["Run"]["UsePileUpWeight"]) {
            wgt/=   pu_weight;
            wgt *=  hPU_up[(int)(nTruePU+1)];
          }
        }else if(syst_names[i]=="Pileup_weight_Down"){
          if(distats["Run"]["UsePileUpWeight"]) {
            wgt/=   pu_weight;
            wgt *=  hPU_down[(int)(nTruePU+1)];
          }
        }
      }
      //get the non particle conditions:
      for(auto itCut : nonParticleCuts){
        active_part->at(itCut)=goodParts.at(itCut);
      }
      if(!fillCuts(false)) continue;
      for(auto it: *syst_histo.get_groups()) {
        fill_Folder(it, i, syst_histo, true);
      }
      wgt=backup_wgt;
    }
  }
  for(Particle* ipart: allParticles) ipart->setCurrentP(0);
  _MET->setCurrentP(0);
  active_part = &goodParts;
}

///Function that fills up the histograms
void Analyzer::fill_Folder(std::string group, const int max, Histogramer &ihisto, bool issyst) {
  return;

}

void Analyzer::fill_Tree(){

  if(0){
    //do our dirty tree stuff here:
    int p1=-1;
    int p2=-1;
    if(active_part->at(CUTS::eDiTau)->size()==1){
      p1= active_part->at(CUTS::eDiTau)->at(0) / BIG_NUM;
      p2= active_part->at(CUTS::eDiTau)->at(0) % BIG_NUM;
    } else{
      return;
    }
    int j1=-1;
    int j2=-1;
    double mass=0;
    for(auto it : *active_part->at(CUTS::eDiJet)) {
      int j1tmp= (it) / _Jet->size();
      int j2tmp= (it) % _Jet->size();
      if(diParticleMass(_Jet->p4(j1tmp),_Jet->p4(j2tmp),"")>mass){
        j1=j1tmp;
        j2=j2tmp;
        mass=diParticleMass(_Jet->p4(j1tmp),_Jet->p4(j2tmp),"");
      }
    }
    if(p1<0 or p2<0 or j1<0 or j2 <0)
      return;
    zBoostTree["tau1_pt"]   = _Tau->pt(p1);
    zBoostTree["tau1_eta"]  = _Tau->eta(p1);
    zBoostTree["tau1_phi"]  = _Tau->phi(p1);
    zBoostTree["tau2_pt"]   = _Tau->pt(p2);
    zBoostTree["tau2_eta"]  = _Tau->eta(p2);
    zBoostTree["tau2_phi"]  = _Tau->phi(p2);
    zBoostTree["tau_mass"]  = diParticleMass(_Tau->p4(p1),_Tau->p4(p2),"");
    zBoostTree["met"]       = _MET->pt();
    zBoostTree["mt_tau1"]   = calculateLeptonMetMt(_Tau->p4(p1));
    zBoostTree["mt_tau2"]   = calculateLeptonMetMt(_Tau->p4(p2));
    zBoostTree["mt2"]       = _MET->MT2(_Tau->p4(p1),_Tau->p4(p2));
    zBoostTree["cosDphi1"]  = absnormPhi(_Tau->phi(p1) - _MET->phi());
    zBoostTree["cosDphi2"]  = absnormPhi(_Tau->phi(p2) - _MET->phi());
    zBoostTree["jet1_pt"]   = _Jet->pt(j1);
    zBoostTree["jet1_eta"]  = _Jet->eta(j1);
    zBoostTree["jet1_phi"]  = _Jet->phi(j1);
    zBoostTree["jet2_pt"]   = _Jet->pt(j2);
    zBoostTree["jet2_eta"]  = _Jet->eta(j2);
    zBoostTree["jet2_phi"]  = _Jet->phi(j2);
    zBoostTree["jet_mass"]  = mass;
    zBoostTree["weight"]    = wgt;

    //put it accidentally in the tree
    histo.fillTree("TauTauTree");
  }
}


/////need to fix this systematic
void Analyzer::initializePileupInfo(std::string MCHisto, std::string DataHisto, std::string DataHistoName, std::string MCHistoName) {

  TFile *file1 = new TFile((PUSPACE+MCHisto).c_str());
  TH1D* histmc = (TH1D*)file1->FindObjectAny(MCHistoName.c_str());
  if(!histmc) throw std::runtime_error("failed to extract histogram");

  TFile* file2 = new TFile((PUSPACE+DataHisto).c_str());
  TH1D* histdata = (TH1D*)file2->FindObjectAny(DataHistoName.c_str());
  if(!histdata) throw std::runtime_error("failed to extract histogram");
  TH1D* histdata_up = (TH1D*)file2->FindObjectAny((DataHistoName+"Up").c_str());
  TH1D* histdata_down = (TH1D*)file2->FindObjectAny((DataHistoName+"Down").c_str());


  histmc->Scale(1./histmc->Integral());
  histdata->Scale(1./histdata->Integral());
  if(histdata_up){
    histdata_up->Scale(1./histdata_up->Integral());
    histdata_down->Scale(1./histdata_down->Integral());
  }

  //double factor = histmc->Integral() / histdata->Integral();
  double value,valueUp,valueDown;
  int ibin=0;
  for(int bin=0; bin < histmc->GetNbinsX(); bin++) {
    ibin=histdata->FindBin(bin);
    if(histmc->GetBinContent(ibin) == 0){
      value = 1;
      valueUp = 1;
      valueDown = 1;
    }else{
      value = histdata->GetBinContent(ibin) / histmc->GetBinContent(ibin);
      if(histdata_up){
        valueUp = histdata->GetBinContent(ibin) / histmc->GetBinContent(ibin);
        valueDown = histdata->GetBinContent(ibin) / histmc->GetBinContent(ibin);
      }
    }
    hPU[bin]      = value;
    if(histdata_up){
      hPU_up[bin]   = valueUp;
      hPU_down[bin] = valueDown;
    }else{
      hPU_up[bin]   = value;
      hPU_down[bin] = value;
    }

  }

  file1->Close();
  file2->Close();

}


///Normalizes phi to be between -PI and PI
double normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}


///Takes the absolute value of of normPhi (made because constant use)
double absnormPhi(double phi) {
  return abs(normPhi(phi));
}
