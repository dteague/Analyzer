/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonCombos(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const json& stats, const int syst) {
  if(! neededCuts.isPresent(ePosFin)) return;
  std::string systname = syst_names.at(syst);

  if(!lep1.needSyst(syst) && !lep2.needSyst(syst)) {
    active_part->at(ePosFin)=goodParts[ePosFin];
    return;
  }

  bool sameParticle = (&lep1 == &lep2);
  TLorentzVector part1, part2;

  for(auto i1 : *active_part->at(ePos1)) {
    part1 = lep1.p4(i1);
    for(auto i2 : *active_part->at(ePos2)) {
      if(sameParticle && i2 <= i1) continue;
      part2 = lep2.p4(i2);
      bool passCuts = true;
      for(auto cut : bset(stats)) {
        if(!passCuts) break;
        else if (cut == "DiscrByDeltaR") passCuts = passCuts && (part1.DeltaR(part2) >= stats["DeltaRCut"]);
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && passCutRange(cos(absnormPhi(part1.Phi() - part2.Phi())), stats["CosDphiCut"]);
        else if(cut == "DiscrByDeltaPt") passCuts = passCuts && passCutRange(part1.Pt() - part2.Pt(), stats["DeltaPtCutValue"]);
        else if(cut == "DiscrByDeltaPtDivSumPt") {
          double ptDiv = (part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt());
          passCuts = passCuts && passCutRange(ptDiv, stats["DeltaPtDivSumPtCutValue"]);
        }
        else if (cut == "DiscrByMassReco") {
          double diMass = diParticleMass(part1,part2, stats["HowCalculateMassReco"]);
          passCuts = passCuts && passCutRange(diMass, stats["MassCut"]);
        }
        else if(cut == "DiscrByCosDphiPtAndMet"){
          double CosDPhi1 = cos(absnormPhi(part1.Phi() - _MET->phi()));
          passCuts = passCuts && passCutRange(CosDPhi1, stats["CosDphiPtAndMetCut"]);
        }


        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }
      if (stats["DiscrByOSLSType"]){
        //   if it is 1 or 0 it will end up in the bool std::map!!
        if(stats["DiscrByOSLSType"] && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
      }else if (! stats["DiscrByOSLSType"].empty()){
        if(lep1.charge(i1) * lep2.charge(i2) > 0) continue;
      }else if (! stats["DiscrByOSLSType"].empty() ){
        if(stats["DiscrByOSLSType"] == "LS" && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
        else if(stats["DiscrByOSLSType"] == "OS" && (lep1.charge(i1) * lep2.charge(i2) >= 0)) continue;
      }

      ///Particles that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2
      if(passCuts)
        active_part->at(ePosFin)->push_back(i1*BIG_NUM + i2);
    }
  }
}

/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonJetCombos(Lepton& lep1, Jet& jet1, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const json& stats, const int syst) {
  if(! neededCuts.isPresent(ePosFin)) return;
  std::string systname = syst_names.at(syst);
  if(!lep1.needSyst(syst) && !jet1.needSyst(syst)) {
    active_part->at(ePosFin)=goodParts[ePosFin];
    return;
  }

  TLorentzVector llep1, ljet1;
  // ----Separation cut between jets (remove overlaps)
  for(auto ij2 : *active_part->at(ePos1)) {
    llep1 = lep1.p4(ij2);
    for(auto ij1 : *active_part->at(ePos2)) {
      ljet1 = _Jet->p4(ij1);

      bool passCuts = true;
      for(auto cut : bset(stats)) {
        if(!passCuts) break;
        else if(cut == "DiscrByDeltaR") passCuts = (ljet1.DeltaR(llep1) >= stats["DeltaRCut"]);
        else if(cut == "DiscrByDeltaEta") passCuts = passCutRange(abs(ljet1.Eta() - llep1.Eta()), stats["DeltaEtaCut"]);
        else if(cut == "DiscrByDeltaPhi") passCuts = passCutRange(absnormPhi(ljet1.Phi() - llep1.Phi()), stats["DeltaPhiCut"]);
        else if(cut == "DiscrByOSEta") passCuts = (ljet1.Eta() * llep1.Eta() < 0);
        else if(cut == "DiscrByMassReco") passCuts = passCutRange((ljet1+llep1).M(), stats["MassCut"]);
        else if(cut == "DiscrByCosDphi") passCuts = passCutRange(cos(absnormPhi(ljet1.Phi() - llep1.Phi())), stats["CosDphiCut"]);
        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2
      if(passCuts) active_part->at(ePosFin)->push_back(ij1*_Jet->size() + ij2);
    }
  }
}


//////////////LOOK INTO DIJET PICKING
///////HOW TO GET RID OF REDUNCENCIES??

/////Same as gooddilepton, just jet specific
void Analyzer::getGoodDiJets(const json& stats, const int syst) {
  if(! neededCuts.isPresent(CUTS::eDiJet)) return;
  std::string systname = syst_names.at(syst);
  if(systname!="orig"){
    //save time to not rerun stuff
    if( systname.find("Jet")==std::string::npos){
      active_part->at(CUTS::eDiJet)=goodParts[CUTS::eDiJet];
      return;
    }
  }
  TLorentzVector jet1, jet2;
  // ----Separation cut between jets (remove overlaps)
  for(auto ij2 : *active_part->at(CUTS::eRJet2)) {
    jet2 = _Jet->p4(ij2);
    for(auto ij1 : *active_part->at(CUTS::eRJet1)) {
      if(ij1 == ij2) continue;
      jet1 = _Jet->p4(ij1);

      bool passCuts = true;
      //cout<<"---------------------"<<std::endl;
      for(auto cut : bset(stats)) {
        //cout<<cut<<"    "<<passCuts<<std::endl;;
        if(!passCuts) break;
        else if(cut == "DiscrByDeltaR") passCuts = passCuts && (jet1.DeltaR(jet2) >= stats["DeltaRCut"]);
        else if(cut == "DiscrByDeltaEta") passCuts = passCuts && passCutRange(abs(jet1.Eta() - jet2.Eta()), stats["DeltaEtaCut"]);
        else if(cut == "DiscrByDeltaPhi") passCuts = passCuts && passCutRange(absnormPhi(jet1.Phi() - jet2.Phi()), stats["DeltaPhiCut"]);
        else if(cut == "DiscrByOSEta") passCuts = passCuts && (jet1.Eta() * jet2.Eta() < 0);
        else if(cut == "DiscrByMassReco") passCuts = passCuts && passCutRange((jet1+jet2).M(), stats["MassCut"]);
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && passCutRange(cos(absnormPhi(jet1.Phi() - jet2.Phi())), stats["CosDphiCut"]);
        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2
      if(passCuts) active_part->at(CUTS::eDiJet)->push_back(ij1*_Jet->size() + ij2);
    }
  }
}


double Analyzer::getZBoostWeight(){
  double boostweigth=1.;
  if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGZ)->size() ==1 || active_part->at(CUTS::eGW)->size() ==1)){
    //cout<<" Z or W " <<std::endl;
    double boostz = 0;
    if(active_part->at(CUTS::eGZ)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGZ)->at(0));
    }
    if(active_part->at(CUTS::eGW)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGW)->at(0));
    }
    if(boostz > 0 && boostz <= 50) {boostweigth = 1.1192;}
    else if (boostz > 50 && boostz <= 100) {boostweigth = 1.1034;}
    else if (boostz > 100 && boostz <= 150) {boostweigth = 1.0675;}
    else if (boostz > 150 && boostz <= 200) {boostweigth = 1.0637;}
    else if (boostz > 200 && boostz <= 300) {boostweigth = 1.0242;}
    else if (boostz > 300 && boostz <= 400) {boostweigth = 0.9453;}
    else if (boostz > 400 && boostz <= 600) {boostweigth = 0.8579;}
    else if (boostz >= 600) {boostweigth = 0.7822;}
    else {boostweigth = 1;}
  }
  return boostweigth;
}


double Analyzer::getWkfactor(){
  double kfactor=1.;
  if(!isWSample)
    return kfactor;
  if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGW)->size() ==1)){
    //this k-factor is not computed for the low mass W!
    double wmass=_Gen->p4(active_part->at(CUTS::eGW)->at(0)).M();
    if(wmass<100){
      return 1.;
    }
    if(active_part->at(CUTS::eGTau)->size()){
      kfactor=k_ele_h->GetBinContent(k_ele_h->FindBin(wmass));
    }
    else if(active_part->at(CUTS::eGMuon)->size()){
      kfactor=k_mu_h->GetBinContent(k_mu_h->FindBin(wmass));
    }
    else if(active_part->at(CUTS::eGElec)->size()){
      kfactor=k_tau_h->GetBinContent(k_tau_h->FindBin(wmass));
    }
  }
  return kfactor;
}

////VBF specific cuts dealing with the leading jets.
void Analyzer::VBFTopologyCut(const json& stats, const int syst) {
  if(! neededCuts.isPresent(CUTS::eSusyCom)) return;
  std::string systname = syst_names.at(syst);


  if(systname!="orig"){
    //only jet stuff is affected
    //save time to not rerun stuff
    if( systname.find("Jet")==std::string::npos){
      active_part->at(CUTS::eSusyCom)=goodParts[CUTS::eSusyCom];
      return;
    }
  }

  TLorentzVector ljet1 = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0));
  TLorentzVector ljet2 = _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));
  TLorentzVector dijet = ljet1 + ljet2;
  double dphi1 = normPhi(ljet1.Phi() - _MET->phi());
  double dphi2 = normPhi(ljet2.Phi() - _MET->phi());

  bool passCuts = true;
  for(auto cut: bset(stats)) {
    if(!passCuts) break;
    else if(cut == "DiscrByMass") passCuts = passCuts && passCutRange(dijet.M(), stats["MassCut"]);
    else if(cut == "DiscrByPt") passCuts = passCuts && passCutRange(dijet.Pt(), stats["PtCut"]);
    else if(cut == "DiscrByDeltaEta") passCuts = passCuts && passCutRange(abs(ljet1.Eta() - ljet2.Eta()), stats["DeltaEtaCut"]);
    else if(cut == "DiscrByDeltaPhi") passCuts = passCuts && passCutRange(absnormPhi(ljet1.Phi() - ljet2.Phi()), stats["DeltaPhiCut"]);
    else if(cut == "DiscrByOSEta") passCuts = passCuts && (ljet1.Eta() * ljet2.Eta() < 0);
    else if(cut == "DiscrByR1") passCuts = passCuts && passCutRange(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0)), stats["R1Cut"]);
    else if(cut == "DiscrByR2") passCuts = passCuts && passCutRange(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), stats["R2Cut"]);
    else if(cut == "DiscrByAlpha") {
      double alpha = (dijet.M() > 0) ? ljet2.Pt() / dijet.M() : -1;
      passCuts = passCuts && passCutRange(alpha, stats["AlphaCut"]);
    }
    else if(cut == "DiscrByDphi1") passCuts = passCuts && passCutRange(abs(dphi1), stats["Dphi1Cut"]);
    else if(cut == "DiscrByDphi2") passCuts = passCuts && passCutRange(abs(dphi2), stats["Dphi2Cut"]);

    else std::cout << "cut: " << cut << " not listed" << std::endl;
  }

  if(passCuts)  active_part->at(CUTS::eSusyCom)->push_back(0);
  return;
}

/////Calculate the diparticle mass based on how to calculate it
///can use Collinear Approximation, which can fail (number failed available in a histogram)
///can use VectorSumOfVisProductAndMet which is sum of particles and met
///Other which is adding without met
double Analyzer::diParticleMass(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, std::string howCalc) {
  bool ratioNotInRange = false;
  TLorentzVector The_LorentzVect;

  if(howCalc == "InvariantMass") {
    return (Tobj1 + Tobj2).M();
  }


  //////check this equation/////
  if(howCalc == "CollinearApprox") {
    double denominator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1 = (Tobj2.Py()*_MET->px() - Tobj2.Px()*_MET->py())/denominator;
    double x2 = (Tobj1.Px()*_MET->py() - Tobj1.Py()*_MET->px())/denominator;
    ratioNotInRange=!((x1 < 0.) && (x2 < 0.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px()*(1 + x1) + Tobj2.Px()*(1+x2)), (Tobj1.Py()*(1+x1) + Tobj2.Py()*(1+x2)), (Tobj1.Pz()*(1+x1) + Tobj2.Pz()*(1+x2)), (Tobj1.Energy()*(1+x1) + Tobj2.Energy()*(1+x2)) );
      return The_LorentzVect.M();
    }
  }

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    return (Tobj1 + Tobj2 + _MET->p4()).M();
  }

  return (Tobj1 + Tobj2).M();
}

////Tests if the CollinearApproximation works for finding the mass of teh particles
bool Analyzer::passDiParticleApprox(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, std::string howCalc) {
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + _MET->px())) - (Tobj2.Px() * (Tobj1.Py() + _MET->py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + _MET->py())) - (Tobj1.Py() * (Tobj2.Px() + _MET->px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    return (x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.);
  } else {
    return true;
  }
}



void Analyzer::initializeWkfactor(std::vector<std::string> infiles) {
  if(infiles[0].find("WJets") != std::string::npos){
    isWSample = true;
  }else{
    isWSample=false;
    return;
  }
  //W-jet k-factor Histograms:
  TFile k_ele("Pileup/k_faktors_ele.root");
  TFile k_mu("Pileup/k_faktors_mu.root");
  TFile k_tau("Pileup/k_faktors_tau.root");

  k_ele_h =dynamic_cast<TH1D*>(k_ele.FindObjectAny("k_fac_m"));
  k_mu_h  =dynamic_cast<TH1D*>(k_mu.FindObjectAny("k_fac_m"));
  k_tau_h =dynamic_cast<TH1D*>(k_tau.FindObjectAny("k_fac_m"));

  k_ele.Close();
  k_mu.Close();
  k_tau.Close();

}
void Analyzer::fill_efficiency() {
  //cut efficiency
  const std::vector<CUTS> goodGenLep={CUTS::eGElec,CUTS::eGMuon,CUTS::eGTau};
  //just the lepton 1 for now
  const std::vector<CUTS> goodRecoLep={CUTS::eRElec1,CUTS::eRMuon1,CUTS::eRTau1};



  for(size_t igen=0;igen<goodGenLep.size();igen++){
    Particle* part =particleCutMap.at(goodGenLep[igen]);
    CUTS cut=goodRecoLep[igen];
    std::smatch mGen;
    std::string tmps=part->getName();
    std::regex_match(tmps, mGen, genName_regex);
    //loop over all gen leptons
    for(int iigen : *active_part->at(goodGenLep[igen])){


      int foundReco=-1;
      for(size_t ireco=0; ireco<part->size(); ireco++){
        if(part->p4(ireco).DeltaR(_Gen->p4(iigen))<0.3){
          foundReco=ireco;
        }
      }
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Pt", _Gen->pt(iigen), foundReco>=0,0);
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Eta",_Gen->eta(iigen),foundReco>=0,0);
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Phi",_Gen->phi(iigen),foundReco>=0,0);
      if(foundReco>=0){
        bool id_particle= (find(active_part->at(cut)->begin(),active_part->at(cut)->end(),foundReco)!=active_part->at(cut)->end());
        histo.addEffiency("eff_"+std::string(mGen[1])+"Pt", _Gen->pt(iigen), id_particle,0);
        histo.addEffiency("eff_"+std::string(mGen[1])+"Eta",_Gen->eta(iigen),id_particle,0);
        histo.addEffiency("eff_"+std::string(mGen[1])+"Phi",_Gen->phi(iigen),id_particle,0);
      }
    }
  }
}

bool Analyzer::select_mc_background(){
  //will return true if Z* mass is smaller than 200GeV
  if(_Gen == nullptr){
    return true;
  }
  if(gen_selection["DY_noMass_gt_200"]){
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      if(abs(_Gen->pdg_id[i])==11 or abs(_Gen->pdg_id[i])==13 or abs(_Gen->pdg_id[i])==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          return (lep1+lep2).M()<200;
        }else{
          lep1= _Gen->p4(i);
        }
      }
    }
  }
  //will return true if Z* mass is smaller than 200GeV
  if(gen_selection["DY_noMass_gt_100"]){
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      if(abs(_Gen->pdg_id[i])==11 or abs(_Gen->pdg_id[i])==13 or abs(_Gen->pdg_id[i])==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          return (lep1+lep2).M()<100;
        }else{
          lep1= _Gen->p4(i);
        }
      }
    }
  }
  //cout<<"Something is rotten in the state of Denmark."<<std::endl;
  //cout<<"could not find gen selection particle"<<std::endl;
  return true;
}
double Analyzer::getTauDataMCScaleFactor(int updown){
  double sf=1.;
  //for(size_t i=0; i<_Tau->size();i++){
  for(auto i : *active_part->at(CUTS::eRTau1)){
    if(matchTauToGen(_Tau->p4(i),0.4)!=TLorentzVector()){

      // if(updown==-1) sf*=  _Tau->pstats["Smear"]["TauSF"] * (1.-(0.35*_Tau->pt(i)/1000.0));
      // else if(updown==0) sf*=  _Tau->pstats["Smear"]["TauSF"];
      // else if(updown==1) sf*=  _Tau->pstats["Smear"]["TauSF"] * (1.+(0.05*_Tau->pt(i)/1000.0));
    }
  }
  return sf;
}

void Analyzer::initializeMCSelection(std::vector<std::string> infiles) {
    // check if we need to make gen level cuts to cross clean the samples:

  isVSample = false;
  if(infiles[0].find("DY") != std::string::npos){
    isVSample = true;
    if(infiles[0].find("DYJetsToLL_M-50_HT-") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
      //gen_selection["DY_noMass_gt_200"]=true;
    //get the DY1Jet DY2Jet ...
    }else if(infiles[0].find("JetsToLL_TuneCUETP8M1_13TeV") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
    }else{
      //set it to false!!
      gen_selection["DY_noMass_gt_100"]=false;
      gen_selection["DY_noMass_gt_200"]=false;
    }

    if(infiles[0].find("DYJetsToLL_M-50_TuneCUETP8M1_13TeV") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
    }else{
      //set it to false!!
      gen_selection["DY_noMass_gt_100"]=false;
      gen_selection["DY_noMass_gt_200"]=false;
    }
  }else{
    //set it to false!!
    gen_selection["DY_noMass_gt_200"]=false;
    gen_selection["DY_noMass_gt_100"]=false;
  }

  if(infiles[0].find("WJets") != std::string::npos){
    isVSample = true;
  }
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
    //   // else if(cut == "DoDiscrByIsolation") {
    //   //   double firstIso = (! stats["IsoSumPtCutValue"].empty()) ? stats.pmap.at("IsoSumPtCutValue").first : ival(ePos) - ival(CUTS::eRTau1) + 1;
    //   //   double secondIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").second : stats["FlipIsolationRequirement"];
    //   //   passCuts = passCuts && lep.get_Iso(i, firstIso, secondIso);
    //   // }
    //   else if(cut == "DiscrIfIsZdecay" && lep.type != PType::Tau ) passCuts = isZdecay(lvec, lep);
    //   else if(cut == "DiscrByMetDphi") passCuts = passCutRange(absnormPhi(lvec.Phi() - _MET->phi()), stats["MetDphiCut"]);
    //   else if(cut == "DiscrByMetMt") passCuts = passCutRange(calculateLeptonMetMt(lvec), stats["MetMtCut"]);
    //   /////muon cuts
    //   else if(lep.type == PType::Muon){
    //     if(cut == "DoDiscrByTightID") passCuts = _Muon->tight[i];
    //     else if(cut == "DoDiscrBySoftID") passCuts = _Muon->soft[i];
    //   }
    //   ////electron cuts
    //   else if(lep.type == PType::Electron){
    //     if(cut == "DoDiscrByHLTID"){
    //       std::bitset<8> idvariable(_Electron->cutBased_HLTPreSel[i]);
    //       if(ival(ePos) - ival(CUTS::eRElec1)){
    //         passCuts = (_Electron->cbHLTIDele1&idvariable).count();
    //       }else{
    //         passCuts = (_Electron->cbHLTIDele2&idvariable).count();
    //       }
    //     }
    //     if(cut == "DoDiscrByCBID"){
    //       std::bitset<8> idvariable(_Electron->cutBased[i]);
    //       if(ival(ePos) - ival(CUTS::eRElec1)){
    //         passCuts = (_Electron->cbIDele1&idvariable).count();
    //       }else{
    //         passCuts = (_Electron->cbIDele2&idvariable).count();
    //       }
    //     }
    //     else if(cut == "DoDiscrByHEEPID")
    //      passCuts = _Electron->isPassHEEPId[i];
    //   }
    //   else if(lep.type == PType::Tau){
    //     if(cut == "DoDiscrByCrackCut") passCuts = !isInTheCracks(lvec.Eta());
    //     /////tau cuts
    //     else if(cut == "DoDzCut") passCuts = (_Tau->dz[i] <= stats["DzCutThreshold"]);
    //     else if(cut == "DoDiscrByLeadTrack") passCuts = (_Tau->leadTkPtOverTauPt[i]*_Tau->pt(i) >= stats["LeadTrackThreshold"]);
    //          // ----Electron and Muon vetos
    //     else if (cut == "DoDiscrAgainstElectron") passCuts = _Tau->pass_against_Elec(ePos, i);
    //     else if (cut == "SelectTausThatAreElectrons") passCuts = !_Tau->pass_against_Elec(ePos, i);

    //     else if (cut == "DoDiscrAgainstMuon") passCuts = _Tau->pass_against_Muon(ePos, i);
    //     else if (cut == "SelectTausThatAreMuons") passCuts = !_Tau->pass_against_Muon(ePos, i);

    //     else if(cut == "DiscrByProngType") {
    // 	  std::string prongName = stats["ProngType"];
    //       passCuts = (prongName.find("hps") == std::string::npos || _Tau->DecayModeNewDMs[i] != 0);
    //       passCuts = passProng(stats["ProngType"], _Tau->decayMode[i]);
    //     }
    //     else if(cut == "decayModeFindingNewDMs") passCuts = _Tau->DecayModeNewDMs[i] != 0;
    //     else if(cut == "decayModeFinding") passCuts = _Tau->DecayMode[i] != 0;
    //           // ----anti-overlap requirements
    //     else if(cut == "RemoveOverlapWithMuon1s") passCuts = !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats["Muon1MatchingDeltaR"]);
    //     else if(cut == "RemoveOverlapWithMuon2s") passCuts = !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats["Muon2MatchingDeltaR"]);
    //     else if(cut == "RemoveOverlapWithElectron1s") passCuts = !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats["Electron1MatchingDeltaR"]);
    //     else if(cut == "RemoveOverlapWithElectron2s") passCuts = !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats["Electron2MatchingDeltaR"]);
    //   }
    //   else std::cout << "cut: " << cut << " not listed" << std::endl;
    // }
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
