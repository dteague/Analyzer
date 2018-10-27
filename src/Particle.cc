#include "Particle.h"
#include <bitset>
#include <signal.h>

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    PARTICLE   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


//particle is a objet that stores multiple versions of the particle candidates
Particle::Particle(TTree* _BOOM, std::string _GenName, std::string filename, std::vector<std::string> _syst_names) : BOOM(_BOOM), GenName(_GenName), syst_names(_syst_names) {
  type = PType::None;
  getPartStats(filename);

  std::regex genName_regex(".*([A-Z][^[:space:]]+)");
  std::regex syst_regex("([A-Za-z]+).+");
  std::smatch mGen, mSyst;
  
  std::regex_match(GenName, mGen, genName_regex);

  for( auto item : syst_names) {
    if(item == "orig") {
      systVec.push_back(new std::vector<TLorentzVector>());
      continue;
    }
    if(!regex_match(item, mSyst, syst_regex)){
      systVec.push_back(nullptr);
      continue;
    }
    if(mGen[1] == mSyst[1]) {
      systVec.push_back(new std::vector<TLorentzVector>());
      std::cout << GenName << ": " << item << std::endl;
    } else {
      systVec.push_back(nullptr);
    }
  }
  m_n=0;

  SetBranch(("n"+GenName).c_str(), m_n);
  SetBranch((GenName+"_pt").c_str(), m_pt);
  SetBranch((GenName+"_eta").c_str(), m_eta);
  SetBranch((GenName+"_phi").c_str(), m_phi);
  SetBranch((GenName+"_mass").c_str(), m_mass);
}

double Particle::pt(uint index)const         {return cur_P->at(index).Pt();}
double Particle::eta(uint index)const        {return cur_P->at(index).Eta();}
double Particle::phi(uint index)const        {return cur_P->at(index).Phi();}
double Particle::energy(uint index)const     {return cur_P->at(index).E();}
double Particle::charge(uint index)const     {return 0;}

uint Particle::size()const                   {return cur_P->size();}
std::vector<TLorentzVector>::iterator Particle::begin(){ return cur_P->begin();}
std::vector<TLorentzVector>::iterator Particle::end(){ return cur_P->end();}
std::vector<TLorentzVector>::const_iterator Particle::begin()const { return cur_P->begin();}
std::vector<TLorentzVector>::const_iterator Particle::end()const { return cur_P->end();}

TLorentzVector Particle::p4(uint index)const {return (cur_P->at(index));}
TLorentzVector& Particle::p4(uint index) {return cur_P->at(index);}
TLorentzVector Particle::RecoP4(uint index)const {return Reco.at(index);}
TLorentzVector& Particle::RecoP4(uint index) {return Reco.at(index);}




void Particle::addPtEtaPhiESyst(double ipt,double ieta, double iphi, double ienergy, int syst){
  TLorentzVector mp4;
  mp4.SetPtEtaPhiE(ipt,ieta,iphi,ienergy);
  systVec.at(syst)->push_back(mp4);
}


void Particle::addP4Syst(TLorentzVector mp4, int syst){
  systVec.at(syst)->push_back(mp4);
}


void Particle::init(){
    //cleanup of the particles
  Reco.clear();
  for(auto it: systVec){
    if(it != nullptr) it->clear();
  }
  TLorentzVector tmp;
  for(uint i=0; i < m_n; i++) {
    tmp.SetPtEtaPhiM(m_pt[i],m_eta[i],m_phi[i],m_mass[i]);
    Reco.push_back(tmp);

  }
  setCurrentP(-1);

}


void Particle::setOrigReco() {
  /////memory loss here if no smear and new std::vector. Only once, so ignore for now...
  systVec.at(0) = &Reco;
}

bool Particle::needSyst(int syst) const {
  return systVec.at(syst) != nullptr;
}


void Particle::setCurrentP(int syst){
  if(syst == -1) {
    cur_P = &Reco;
  } else if(systVec.at(syst) == nullptr || systVec.at(syst)->size() == 0) {
    cur_P = systVec.at(0);  //orig
  } else {
    cur_P = systVec.at(syst);
  }
  //  activeSystematic=syst;
}


void Particle::unBranch() {
  BOOM->SetBranchStatus((GenName+"*").c_str(), 0);
  BOOM->SetBranchStatus(("*"+GenName).c_str(), 0);
}


void Particle::getPartStats(std::string filename) {
  std::ifstream info_file(filename);

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    return;
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
  ss >>pstats;

  info_file.close();
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////    PHOTON  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Photon::Photon(TTree* _BOOM, std::string filename, std::vector<std::string> syst_names) : Particle(_BOOM, "Photon", filename, syst_names) {
  SetBranch("Photon_hoe", hoverE);
  SetBranch("Photon_r9", phoR);
  SetBranch("Photon_sieie", sigmaIEtaIEta);
  SetBranch("Photon_pfRelIso03_all", pfIso_all);
  SetBranch("Photon_pfRelIso03_chg", pfIso_chg);
  SetBranch("Photon_electronVeto", eleVeto);
  SetBranch("Photon_pixelSeed", hasPixelSeed);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    GENERATED   ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Generated::Generated(TTree* _BOOM, std::string filename, std::vector<std::string> syst_names) : Particle(_BOOM, "GenPart", filename, syst_names) {

  SetBranch("GenPart_pdgId", pdg_id);
  SetBranch("GenPart_genPartIdxMother", genPartIdxMother);
  SetBranch("GenPart_status", status);
  SetBranch("GenPart_statusFlags", statusFlags);
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////    JET  ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Jet::Jet(TTree* _BOOM, std::string filename, std::vector<std::string> syst_names) : Particle(_BOOM, "Jet", filename, syst_names) {
  type = PType::Jet;
  
  SetBranch("Jet_jetId", jetId);
  SetBranch("Jet_neHEF", neutralHadEnergyFraction);
  SetBranch("Jet_neEmEF", neutralEmEmEnergyFraction);
  SetBranch("Jet_nConstituents", numberOfConstituents);
  SetBranch("Jet_nMuons", nMuons);
  SetBranch("Jet_chHEF", chargedHadronEnergyFraction);
  SetBranch("Jet_chEmEF", chargedEmEnergyFraction);
  SetBranch("Jet_btagCSVV2", bDiscriminator);
  SetBranch("Jet_puId", puID);
  SetBranch("Jet_area", area);
  if(_BOOM->FindBranch("Jet_partonFlavour")!=0){
    SetBranch("Jet_partonFlavour", partonFlavour);
  }
  
  
}

std::vector<CUTS> Jet::findExtraCuts() {
  std::vector<CUTS> return_vec;
  if(pstats["Smear"]["SmearTheJet"]) {
    return_vec.push_back(CUTS::eGen);
  }
  return return_vec;
}

std::vector<CUTS> Jet::overlapCuts(CUTS ePos) {
  std::vector<CUTS> returnCuts;
  auto& tmpset = pstats[jetNameMap.at(ePos)];
  if(tmpset["RemoveOverlapWithMuon1s"]) returnCuts.push_back(CUTS::eRMuon1);
  if(tmpset["RemoveOverlapWithMuon2s"]) returnCuts.push_back(CUTS::eRMuon2);
  if(tmpset["RemoveOverlapWithElectron1s"]) returnCuts.push_back(CUTS::eRElec1);
  if(tmpset["RemoveOverlapWithElectron2s"]) returnCuts.push_back(CUTS::eRElec2);
  if(tmpset["RemoveOverlapWithTau1s"]) returnCuts.push_back(CUTS::eRTau1);
  if(tmpset["RemoveOverlapWithTau2s"]) returnCuts.push_back(CUTS::eRTau2);

  return returnCuts;
}

bool Jet::passedLooseJetID(int nobj) {
  std::bitset<8> bit_jet(jetId[nobj]);
  return bit_jet[0];
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    FATJET   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


FatJet::FatJet(TTree* _BOOM, std::string filename, std::vector<std::string> syst_names) : Particle(_BOOM, "FatJet", filename, syst_names) {
  type = PType::FatJet;

  SetBranch("FatJet_tau1", tau1);
  SetBranch("FatJet_tau2", tau2);
  SetBranch("FatJet_tau3", tau3);
  SetBranch("FatJet_tau4", tau4);
  SetBranch("FatJet_mass", PrunedMass);
  SetBranch("FatJet_msoftdrop", SoftDropMass);
}

std::vector<CUTS> FatJet::findExtraCuts() {
  std::vector<CUTS> return_vec;
  if(pstats["Smear"]["SmearTheJet"]) {
    return_vec.push_back(CUTS::eGen);
  }
  return return_vec;
}

std::vector<CUTS> FatJet::overlapCuts(CUTS ePos) {
  std::vector<CUTS> returnCuts;
  auto& tmpset = pstats[jetNameMap.at(ePos)];
  if(tmpset["RemoveOverlapWithMuon1s"]) returnCuts.push_back(CUTS::eRMuon1);
  if(tmpset["RemoveOverlapWithMuon2s"]) returnCuts.push_back(CUTS::eRMuon2);
  if(tmpset["RemoveOverlapWithElectron1s"]) returnCuts.push_back(CUTS::eRElec1);
  if(tmpset["RemoveOverlapWithElectron2s"]) returnCuts.push_back(CUTS::eRElec2);
  if(tmpset["RemoveOverlapWithTau1s"]) returnCuts.push_back(CUTS::eRTau1);
  if(tmpset["RemoveOverlapWithTau2s"]) returnCuts.push_back(CUTS::eRTau2);

  return returnCuts;
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////    LEPTON   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Lepton::Lepton(TTree* _BOOM, std::string GenName, std::string EndName, std::vector<std::string> syst_names) : Particle(_BOOM, GenName, EndName, syst_names) {
  SetBranch((GenName+"_charge").c_str(), _charge);
}

std::vector<CUTS> Lepton::findExtraCuts() {
  std::vector<CUTS> return_vec;
  auto& tmpset = pstats["Smear"];
  if(tmpset["SmearTheParticle"]|| tmpset["MatchToGen"]) {
    return_vec.push_back(cutMap.at(type));
  }
  if(tmpset["doEfficiencyPlots"]){
    return_vec.push_back(cutMap.at(type));
  }
  return return_vec;
}

double Lepton::charge(uint index)const     {return _charge[index];}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    ELECTRON   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Electron::Electron(TTree* _BOOM, std::string filename, std::vector<std::string> syst_names) : Lepton(_BOOM, "Electron", filename, syst_names) {
  type = PType::Electron;
  auto& elec1 = pstats["Elec1"];
  auto& elec2 = pstats["Elec2"];
  
  std::bitset<8> tmp(elec1["DiscrByCBID"]);
  cbIDele1=tmp;
  tmp=static_cast<std::bitset<8>>(elec2["DiscrByCBID"]);
  cbIDele2=tmp;
  
  tmp=static_cast<std::bitset<8>>(elec1["DiscrByHLTID"]);
  cbHLTIDele1=tmp;
  tmp=static_cast<std::bitset<8>>(elec2["DiscrByHLTID"]);
  cbHLTIDele2=tmp;

  if(_BOOM->FindBranch("Electron_mvaSpring16GP")!=0){
    std::cout<<"Electron MVA ID: Electron_mvaSpring16"<<std::endl;
  } else{
    std::cout<<"Electron MVA ID: Electron_mvaFall17"<<std::endl;
  }
 
  if((elec1["DoDiscrByIsolation"]|| elec2["DoDiscrByIsolation"]) && _BOOM->FindBranch("Electron_mvaFall17Iso")!=0 ) {
   SetBranch("Electron_miniPFRelIso_all", miniPFRelIso_all);
   SetBranch("Electron_miniPFRelIso_chg", miniPFRelIso_chg);
   SetBranch("Electron_mvaFall17Iso", mvaFall17Iso);
   SetBranch("Electron_mvaFall17noIso", mvaFall17noIso);
   SetBranch("Electron_pfRelIso03_all", pfRelIso03_all);
   SetBranch("Electron_pfRelIso03_chg", pfRelIso03_chg);
  }

  if((elec1["DoDiscrByIsolation"]|| elec2["DoDiscrByIsolation"]) && _BOOM->FindBranch("Electron_mvaSpring16GP")!=0 ) {
   SetBranch("Electron_miniPFRelIso_all", miniPFRelIso_all);
   SetBranch("Electron_miniPFRelIso_chg", miniPFRelIso_chg);
   SetBranch("Electron_mvaSpring16GP", mvaSpring16GP);
   SetBranch("Electron_mvaSpring16HZZ", mvaSpring16HZZ);
   SetBranch("Electron_pfRelIso03_all", pfRelIso03_all);
   SetBranch("Electron_pfRelIso03_chg", pfRelIso03_chg);
  }

  if(elec1["DoDiscrBycutBasedID"]|| elec2["DoDiscrByLooseID"]) {
    SetBranch("Electron_cutBased", cutBased);
    SetBranch("Electron_cutBased_HLTPreSel", cutBased_HLTPreSel);
  }
  
  if((elec1["DoDiscrBymvaID"]|| elec2["DoDiscrByLooseID"]) && _BOOM->FindBranch("Electron_mvaFall17Iso")!=0){
    SetBranch("Electron_mvaFall17Iso_WP90", mvaIso_90);
    SetBranch("Electron_mvaFall17noIso_WP90", mvanoIso_WP90);
    SetBranch("Electron_mvaFall17Iso_WP80", mvaIso_80);
    SetBranch("Electron_mvaFall17noIso_WP80", mvanoIso_WP80);
    SetBranch("Electron_mvaFall17Iso_WPL", mvaIso_WPL);
    SetBranch("Electron_mvaFall17noIso_WPL", mvanoIso_WPL);
  }

  if((elec1["DoDiscrBymvaID"]|| elec2["DoDiscrByLooseID"]) && _BOOM->FindBranch("Electron_mvaSpring16GP_WP90")!=0){
    SetBranch("Electron_mvaSpring16GP_WP90", mvaGP_90);
    SetBranch("Electron_mvaSpring16GP_WP80", mvaGP_80);
    SetBranch("Electron_mvaSpring16HZZ_WPL", mvaHZZ_WPL);
    SetBranch("Electron_mvaTTH", mvaTTH); 
  }

  if(elec1["DoDiscrByHEEPID"]|| elec2["DoDiscrByHEEPID"]) {
    SetBranch("Electron_cutBased_HEEP", isPassHEEPId);
  }
}

bool Electron::get_Iso(int index, double min, double max) const {
  return (pfRelIso03_all[index] >= min && pfRelIso03_all[index] < max);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////   MUON  // ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Muon::Muon(TTree* _BOOM, std::string filename, std::vector<std::string> syst_names) : Lepton(_BOOM, "Muon", filename, syst_names) {
  type = PType::Muon;
  auto& mu1 = pstats["Muon1"];
  auto& mu2 = pstats["Muon2"];

  if(mu1["DoDiscrByTightID"]|| mu2["DoDiscrByTightID"]) {
    SetBranch("Muon_tightId", tight);
     }
  if(mu1["DoDiscrBySoftID"]|| mu2["DoDiscrBySoftID"]) {
    SetBranch("Muon_softId", soft);
  }
  if(mu1["DoDiscrByIsolation"]|| mu2["DoDiscrByIsolation"]) {
    SetBranch("Muon_miniPFRelIso_all", miniPFRelIso_all);
    SetBranch("Muon_miniPFRelIso_chg", miniPFRelIso_chg);
    SetBranch("Muon_pfRelIso03_all"  , pfRelIso03_all);
    SetBranch("Muon_pfRelIso03_chg"  , pfRelIso03_chg);
    SetBranch("Muon_pfRelIso04_all"  , pfRelIso04_all);
  }
}

bool Muon::get_Iso(int index, double min, double max) const {
  return (pfRelIso03_all[index] >= min && pfRelIso03_all[index] < max);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    TAUS    ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


Taus::Taus(TTree* _BOOM, std::string filename, std::vector<std::string> syst_names) : Lepton(_BOOM, "Tau", filename, syst_names) {
  type = PType::Tau;
  
  std::bitset<8> tmp(pstats["Tau1"]["DiscrByMinIsolation"]);
  tau1minIso=tmp;
  tmp=static_cast<std::bitset<8>>(pstats["Tau1"]["DiscrByMaxIsolation"]);
  tau1maxIso=tmp;
  
  tmp=static_cast<std::bitset<8>>(pstats["Tau2"]["DiscrByMinIsolation"]);
  tau2minIso=tmp;
  tmp=static_cast<std::bitset<8>>(pstats["Tau2"]["DiscrByMaxIsolation"]);
  tau2maxIso=tmp;
  
  tmp=static_cast<std::bitset<8>>(pstats["Tau1"]["DiscrAgainstElectron"]);
  tau1ele=tmp;
  tmp=static_cast<std::bitset<8>>( pstats["Tau1"]["DiscrAgainstMuon"]);
  tau1mu= tmp;
  
  tmp=static_cast<std::bitset<8>>(pstats["Tau2"]["DiscrAgainstElectron"]);
  tau2ele=tmp;
  tmp= static_cast<std::bitset<8>>(pstats["Tau2"]["DiscrAgainstMuon"]);
  tau2mu= tmp;

  SetBranch("Tau_idAntiEle", againstElectron);
  SetBranch("Tau_idAntiMu",  againstMuon);
  SetBranch("Tau_idDecayMode",  DecayMode);
  SetBranch("Tau_idDecayModeNewDMs",  DecayModeNewDMs);
  SetBranch("Tau_idMVAoldDM",  MVAoldDM);
  SetBranch("Tau_decayMode", decayMode);
  SetBranch("Tau_leadTkDeltaEta", leadTkDeltaEta);
  SetBranch("Tau_leadTkDeltaPhi", leadTkDeltaPhi);
  SetBranch("Tau_leadTkPtOverTauPt", leadTkPtOverTauPt);
  SetBranch("Tau_dz", dz);
  SetBranch("Tau_dxy", dxy);
  SetBranch("Tau_chargedIso", chargedIsoPtSum);
  SetBranch("Tau_neutralIso", neutralIso);
  SetBranch("Tau_puCorr", puCorr);
  

}
std::vector<CUTS> Taus::findExtraCuts() {
  std::vector<CUTS> return_vec = Lepton::findExtraCuts();

  auto& tau1 = pstats["Tau1"];
  auto& tau2 = pstats["Tau2"];
  if(tau1["RemoveOverlapWithMuon1s"] || tau2["RemoveOverlapWithMuon1s"])
    return_vec.push_back(CUTS::eRMuon1);
  if(tau1["RemoveOverlapWithMuon2s"] || tau2["RemoveOverlapWithMuon2s"])
    return_vec.push_back(CUTS::eRMuon2);
  if(tau1["RemoveOverlapWithElectron1s"] || tau2["RemoveOverlapWithElectron1s"])
    return_vec.push_back(CUTS::eRElec1);
  if(tau1["RemoveOverlapWithElectron2s"] || tau2["RemoveOverlapWithElectron2s"])
    return_vec.push_back(CUTS::eRElec2);

  return return_vec;
}

//onetwo is 1 for the first 0 for the second
bool Taus::get_Iso(int index, double onetwo, double flipisolation) const {
  std::bitset<8> tau_iso(MVAoldDM[index]);
  std::bitset<8> tau_isomin_mask(tau1minIso);
  std::bitset<8> tau_isomax_mask(tau1maxIso);
  
  if(onetwo != 1 ){
    tau_isomin_mask=tau2minIso;
    tau_isomax_mask=tau2maxIso;
  }
  
  if(!flipisolation){
    //cout<<tau_isomax_mask<<"  "<<tau_iso<<"   "<<(tau_isomax_mask& tau_iso)<<"   "<<  (tau_isomax_mask& tau_iso).count()<<std::endl; 
    return (tau_isomax_mask& tau_iso).count();
  }else{
    return(!((tau_isomax_mask&tau_iso).count()) and (tau_isomin_mask&tau_iso).count());
  }
}

bool Taus::pass_against_Elec(CUTS ePos, int index) {
  std::bitset<8> tau_ele(againstElectron[index]);
  std::bitset<8> tau_ele_mask(tau1ele);
  if(ePos == CUTS::eRTau2){
    std::bitset<8> tmp(tau2ele);
    tau_ele_mask=tmp;
  }
  //cout<<tau_ele_mask<<"  "<<tau_ele<<"   "<<(tau_ele_mask&tau_ele)<<"   "<<  (tau_ele_mask&tau_ele).count()<<std::endl; 
  return (tau_ele_mask&tau_ele).count();
}

bool Taus::pass_against_Muon(CUTS ePos, int index) {
  std::bitset<8> tau_mu(againstMuon[index]);
  std::bitset<8> tau_mu_mask(tau1mu);
  if(ePos == CUTS::eRTau2){
    std::bitset<8> tmp(tau2mu);
    tau_mu_mask=tmp;
  }
  //cout<<tau_mu_mask<<"  "<<tau_mu<<"   "<<(tau_mu_mask&tau_mu)<<"   "<<  (tau_mu_mask&tau_mu).count()<<std::endl; 
  return (tau_mu_mask&tau_mu).count();
}
