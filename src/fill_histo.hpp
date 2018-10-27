  /*be aware in this function
   * the following definition is used:
   * histAddVal(val, name) histo.addVal(val, group, max, name, wgt)
   * so each histogram knows the group, max and weight!
   */
  if(group == "FillRun" && (&ihisto==&histo)) {
    if(crbins != 1) {
      for(int i = 0; i < crbins; i++) {
        ihisto.addVal(false, group, i, "Events", 1);
        if(distats["Run"]["ApplyGenWeight"]) {
          //put the weighted events in bin 3
          ihisto.addVal(2, group,i, "Events", (gen_weight > 0) ? 1.0 : -1.0);
        }
        ihisto.addVal(wgt, group, i, "Weight", 1);
      }
    }
    else{
      ihisto.addVal(false, group,ihisto.get_maxfolder(), "Events", 1);
      if(distats["Run"]["ApplyGenWeight"]) {
        //put the weighted events in bin 3
        ihisto.addVal(2, group,ihisto.get_maxfolder(), "Events", (gen_weight > 0) ? 1.0 : -1.0);
      }
      ihisto.addVal(wgt, group, ihisto.get_maxfolder(), "Weight", 1);
      ihisto.addVal(nTruePU, group, ihisto.get_maxfolder(), "PUWeight", 1);
    }
    histAddVal(true, "Events");
    histAddVal(bestVertices, "NVertices");
  } else if(!isData && group == "FillGen") {

    int nhadtau = 0;
    TLorentzVector genVec;
    int i = 0;
    for(vec_iter it=active_part->at(CUTS::eGTau)->begin(); it!=active_part->at(CUTS::eGTau)->end(); it++, i++) {
      int nu = active_part->at(CUTS::eNuTau)->at(i);
      if(nu != -1) {
        genVec = _Gen->p4(*it) - _Gen->p4(nu);
        histAddVal(genVec.Pt(), "HadTauPt");
        histAddVal(genVec.Eta(), "HadTauEta");
        nhadtau++;
      }
      histAddVal(_Gen->energy(*it), "TauEnergy");
      histAddVal(_Gen->pt(*it), "TauPt");
      histAddVal(_Gen->eta(*it), "TauEta");
      histAddVal(_Gen->phi(*it), "TauPhi");
      for(vec_iter it2=it+1; it2!=active_part->at(CUTS::eGTau)->end(); it2++) {
        histAddVal(diParticleMass(_Gen->p4(*it),_Gen->p4(*it2), "none"), "DiTauMass");
      }
    }
    histAddVal(active_part->at(CUTS::eGTau)->size(), "NTau");
    histAddVal(nhadtau, "NHadTau");

    for(auto it : *active_part->at(CUTS::eGZ)) {
      histAddVal(_Gen->pt(it), "ZPt");
      histAddVal(_Gen->eta(it), "ZEta");
      histAddVal(_Gen->p4(it).M(), "ZMass");
    }
    histAddVal(active_part->at(CUTS::eGZ)->size(), "NZ");

    for(auto it : *active_part->at(CUTS::eGW)) {
      histAddVal(_Gen->pt(it), "WPt");
      histAddVal(_Gen->eta(it), "WEta");
      histAddVal(_Gen->p4(it).M(), "WMass");
    }
    histAddVal(active_part->at(CUTS::eGW)->size(), "NW");


    for(auto it : *active_part->at(CUTS::eGMuon)) {
      histAddVal(_Gen->energy(it), "MuonEnergy");
      histAddVal(_Gen->pt(it), "MuonPt");
      histAddVal(_Gen->eta(it), "MuonEta");
      histAddVal(_Gen->phi(it), "MuonPhi");
    }
    histAddVal(active_part->at(CUTS::eGMuon)->size(), "NMuon");

    double mass=0;
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t igen=0; igen<_Gen->size(); igen++){
      //if a Z boson is explicitly there
      if(abs(_Gen->pdg_id[igen])==11 or abs(_Gen->pdg_id[igen])==13 or abs(_Gen->pdg_id[igen])==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(igen);
          mass=(lep1+lep2).M();
          //cout<<"mass  leptons "<<mass<<std::endl;
          break;
        }else{
          //cout<<_Gen->size()<<"   "<<igen<<std::endl;
          //if(_Gen->size()>_Gen->cur_P.size()){
           //_Gen->init();
          //}
          lep1= _Gen->RecoP4(igen);
        }
      }
    }
    histAddVal(mass, "LeptonMass");
  } else if(fillInfo[group]->type == FILLER::Single) {
    Particle* part = fillInfo[group]->part;
    CUTS ePos = fillInfo[group]->ePos;

    for(auto it : *active_part->at(ePos)) {
      histAddVal(part->p4(it).Energy(), "Energy");
      histAddVal(part->p4(it).Pt(), "Pt");
      histAddVal(part->p4(it).Eta(), "Eta");
      histAddVal(part->p4(it).Phi(), "Phi");
      histAddVal(part->p4(it).DeltaPhi(_MET->p4()), "MetDphi");
      if(part->type == PType::Tau) {
        //if(_Tau->nProngs->at(it) == 1){
          //histAddVal(part->pt(it), "Pt_1prong");
        //}else if(_Tau->nProngs->at(it) == 3){
          //histAddVal(part->pt(it), "Pt_3prong");
        //}
        //histAddVal(_Tau->nProngs->at(it), "NumSignalTracks");
        histAddVal(_Tau->charge(it), "Charge");
        histAddVal(_Tau->againstElectron[it], "againstElectron");
        histAddVal(_Tau->againstMuon[it], "againstMuon");
        histAddVal(_Tau->DecayMode[it], "DecayMode");
        histAddVal(_Tau->DecayModeNewDMs[it], "DecayModeNewDMs");
        histAddVal(_Tau->MVAoldDM[it], "MVAoldDM");
        histAddVal(_Tau->decayMode[it], "decayMode");
        histAddVal(_Tau->leadTkDeltaEta[it], "leadTkDeltaEta");
        histAddVal(_Tau->leadTkDeltaPhi[it], "leadTkDeltaPhi");
        histAddVal(_Tau->leadTkPtOverTauPt[it], "leadTkPtOverTauPt");
        histAddVal(_Tau->dz[it], "dz");
        //histAddVal(_Tau->leadChargedCandPt->at(it), "SeedTrackPt");
        //histAddVal(_Tau->leadChargedCandDz_pv->at(it), "leadChargedCandDz");
      }
      if(part->type != PType::Jet) {
        histAddVal(calculateLeptonMetMt(part->p4(it)), "MetMt");
      }
      if(part->type == PType::FatJet ) {
        histAddVal(_FatJet->PrunedMass[it], "PrunedMass");
        histAddVal(_FatJet->SoftDropMass[it], "SoftDropMass");
        histAddVal(_FatJet->tau1[it], "tau1");
        histAddVal(_FatJet->tau2[it], "tau2");
        histAddVal(_FatJet->tau2[it]/_FatJet->tau1[it], "tau2Overtau1");
      }
    }

    if((part->type != PType::Jet ) && active_part->at(ePos)->size() > 0) {
      std::vector<std::pair<double, int> > ptIndexVector;
      for(auto it : *active_part->at(ePos)) {
        ptIndexVector.push_back(std::make_pair(part->pt(it),it));
      }
      sort(ptIndexVector.begin(),ptIndexVector.end());
      if(ptIndexVector.size()>0){
        histAddVal(part->pt(ptIndexVector.back().second), "FirstLeadingPt");
        histAddVal(part->eta(ptIndexVector.back().second), "FirstLeadingEta");
      }
      if(ptIndexVector.size()>1){
        histAddVal(part->pt(ptIndexVector.at(ptIndexVector.size()-2).second), "SecondLeadingPt");
        histAddVal(part->eta(ptIndexVector.at(ptIndexVector.size()-2).second), "SecondLeadingEta");
      }
    }

    histAddVal(active_part->at(ePos)->size(), "N");


  } else if(group == "FillMetCuts") {
    histAddVal(_MET->MHT(), "MHT");
    histAddVal(_MET->HT(), "HT");
    histAddVal(_MET->HT() + _MET->MHT(), "Meff");
    histAddVal(_MET->pt(), "Met");
    histAddVal(_MET->phi(), "MetPhi");

  } else if(group == "FillLeadingJet" && active_part->at(CUTS::eSusyCom)->size() == 0) {

    if(active_part->at(CUTS::eR1stJet)->size()>0) {
      histAddVal(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)).Pt(), "FirstPt");
      histAddVal(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)).Eta(), "FirstEta");
    }
    if(active_part->at(CUTS::eR2ndJet)->size()>0) {
      histAddVal(_Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0)).Pt(), "SecondPt");
      histAddVal(_Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0)).Eta(), "SecondEta");
    }


  } else if(group == "FillLeadingJet" && active_part->at(CUTS::eSusyCom)->size() != 0) {

    TLorentzVector first = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0));
    TLorentzVector second = _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

    histAddVal(first.Pt(), "FirstPt");
    histAddVal(second.Pt(), "SecondPt");

    histAddVal(first.Eta(), "FirstEta");
    histAddVal(second.Eta(), "SecondEta");

    TLorentzVector LeadDiJet = first + second;

    histAddVal(LeadDiJet.M(), "Mass");
    histAddVal(LeadDiJet.Pt(), "Pt");
    histAddVal(fabs(first.Eta() - second.Eta()), "DeltaEta");
    histAddVal(first.DeltaR(second), "DeltaR");

    double dphiDijets = absnormPhi(first.Phi() - second.Phi());
    double dphi1 = normPhi(first.Phi() - _MET->phi());
    double dphi2 = normPhi(second.Phi() - _MET->phi());
    double alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histAddVal(dphiDijets, "LeadSublDijetDphi");
    histAddVal2(_MET->pt(),dphiDijets, "MetVsDiJetDeltaPhiLeadSubl");
    histAddVal2(fabs(first.Eta()-second.Eta()), dphiDijets, "DeltaEtaVsDeltaPhiLeadSubl");

    histAddVal(absnormPhi(_MET->phi() - LeadDiJet.Phi()), "MetDeltaPhi");


    histAddVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), "R1");
    histAddVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), "R2");
    histAddVal(normPhi(first.Phi() - _MET->MHTphi()), "Dphi1MHT");
    histAddVal(normPhi(second.Phi() - _MET->MHTphi()), "Dphi2MHT");
    histAddVal(dphi1, "Dphi1");
    histAddVal(dphi2, "Dphi2");
    histAddVal2(dphi1,dphi2, "Dphi1VsDphi2");
    histAddVal(alpha, "Alpha");


    //dijet info
  } else if(group == "FillDiJet") {
    double leaddijetmass = 0;
    double leaddijetpt = 0;
    double leaddijetdeltaR = 0;
    double leaddijetdeltaEta = 0;
    double etaproduct = 0;
    for(auto it : *active_part->at(CUTS::eDiJet)) {
      int p1 = (it) / _Jet->size();
      int p2 = (it) % _Jet->size();
      TLorentzVector jet1 = _Jet->p4(p1);
      TLorentzVector jet2 = _Jet->p4(p2);
      TLorentzVector DiJet = jet1 + jet2;

      if(DiJet.M() > leaddijetmass) {
        leaddijetmass = DiJet.M();
        etaproduct = (jet1.Eta() * jet2.Eta() > 0) ? 1 : -1;
      }
      if(DiJet.Pt() > leaddijetpt) leaddijetpt = DiJet.Pt();
      if(fabs(jet1.Eta() - jet2.Eta()) > leaddijetdeltaEta) leaddijetdeltaEta = fabs(jet1.Eta() - jet2.Eta());
      if(jet1.DeltaR(jet2) > leaddijetdeltaR) leaddijetdeltaR = jet1.DeltaR(jet2);

      histAddVal(DiJet.M(), "Mass");
      histAddVal(DiJet.Pt(), "Pt");
      histAddVal(fabs(jet1.Eta() - jet2.Eta()), "DeltaEta");
      histAddVal(absnormPhi(jet1.Phi() - jet2.Phi()), "DeltaPhi");
      histAddVal(jet1.DeltaR(jet2), "DeltaR");
    }


    histAddVal(leaddijetmass, "LargestMass");
    histAddVal(leaddijetpt, "LargestPt");
    histAddVal(leaddijetdeltaEta, "LargestDeltaEta");
    histAddVal(leaddijetdeltaR, "LargestDeltaR");
    histAddVal(etaproduct, "LargestMassEtaProduct");

    for(auto index : *(active_part->at(CUTS::eRTau1)) ) {
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetmass, "mTvsLeadingMass");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetdeltaEta, "mTvsLeadingDeltaEta");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetdeltaR, "mTvsLeadingDeltaR");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetpt, "mTvsLeadingPt");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetmass, "MetDphiVSLeadingMass");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetdeltaEta, "MetDphiVSLeadingDeltaEta");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetdeltaR, "MetDphiVSLeadingDeltaR");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetpt, "MetDphiVSLeadingPt");
    }



    ////diparticle stuff

  } else if(fillInfo[group]->type == FILLER::Dilepjet) {
    Jet* jet = static_cast<Jet*>(fillInfo[group]->part);
    Lepton* lep = static_cast<Lepton*>(fillInfo[group]->part2);
    CUTS ePos = fillInfo[group]->ePos;
    std::string digroup = group;
    digroup.erase(0,4);

    TLorentzVector part1;
    TLorentzVector part2;

    for(auto it : *active_part->at(ePos)) {

      int p1= (it) / _Jet->size();;
      int p2= (it) % _Jet->size();;

      part1 = lep->p4(p1);
      part2 = jet->p4(p2);

      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR");
      if(group.find("Di") != std::string::npos) {
        histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
      } else {
        histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");
      histAddVal(absnormPhi(part1.Phi() - _MET->phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - _MET->phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - _MET->phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - _MET->phi())), "CosDphi_DeltaPtAndMet");

      double diMass = diParticleMass(part1,part2, distats[digroup]["HowCalculateMassReco"]);
      if(passDiParticleApprox(part1,part2, distats[digroup]["HowCalculateMassReco"])) {
        histAddVal(diMass, "ReconstructableMass");
      } else {
        histAddVal(diMass, "NotReconstructableMass");
      }

      if ((active_part->at(CUTS::eR1stJet)->size()>0 && active_part->at(CUTS::eR1stJet)->at(0) != -1) && (active_part->at(CUTS::eR2ndJet)->size()>0 && active_part->at(CUTS::eR2ndJet)->at(0) != -1)) {
        TLorentzVector TheLeadDiJetVect = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)) + _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

        histAddVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), "Part1DiJetDeltaPhi");
        histAddVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), "Part2DiJetDeltaPhi");
        histAddVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), "DiJetReconstructableMass");
      }
    }
  } else if(fillInfo[group]->type == FILLER::Dipart) {
    Lepton* lep1 = static_cast<Lepton*>(fillInfo[group]->part);
    Lepton* lep2 = static_cast<Lepton*>(fillInfo[group]->part2);
    CUTS ePos = fillInfo[group]->ePos;
    std::string digroup = group;
    digroup.erase(0,4);

    TLorentzVector part1;
    TLorentzVector part2;

    for(auto it : *active_part->at(ePos)) {

      int p1= (it) / BIG_NUM;
      int p2= (it) % BIG_NUM;

      part1 = lep1->p4(p1);
      part2 = lep2->p4(p2);

      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR");
      if(group.find("Di") != std::string::npos) {
        histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
      } else {
        histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");

      histAddVal(cos(absnormPhi(part1.Phi() - _MET->phi())), "Part1CosDphiPtandMet");
      histAddVal(cos(absnormPhi(part2.Phi() - _MET->phi())), "Part2CosDphiPtandMet");


      histAddVal(absnormPhi(part1.Phi() - _MET->phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - _MET->phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - _MET->phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - _MET->phi())), "CosDphi_DeltaPtAndMet");

      double diMass = diParticleMass(part1,part2, distats[digroup]["HowCalculateMassReco"]);
      if(passDiParticleApprox(part1,part2, distats[digroup]["HowCalculateMassReco"])) {
        histAddVal(diMass, "ReconstructableMass");
      } else {
        histAddVal(diMass, "NotReconstructableMass");
      }

      double InvMass = diParticleMass(part1,part2, "InvariantMass");
      histAddVal(InvMass, "InvariantMass");

      double ptSum = part1.Pt() + part2.Pt();
      histAddVal(ptSum, "SumOfPt");

      if ((active_part->at(CUTS::eR1stJet)->size()>0 && active_part->at(CUTS::eR1stJet)->at(0) != -1) && (active_part->at(CUTS::eR2ndJet)->size()>0 && active_part->at(CUTS::eR2ndJet)->at(0) != -1)) {
        TLorentzVector TheLeadDiJetVect = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)) + _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

        histAddVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), "Part1DiJetDeltaPhi");
        histAddVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), "Part2DiJetDeltaPhi");
        histAddVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), "DiJetReconstructableMass");
      }

      if(lep1->type != PType::Tau) {
        histAddVal(isZdecay(part1, *lep1), "Part1IsZdecay");
      }
      if(lep2->type != PType::Tau){
        histAddVal(isZdecay(part2, *lep2), "Part2IsZdecay");
      }


      //electron tau stuff:
      if(lep1->type == PType::Electron && lep2->type == PType::Electron){
        //loop over taus to find a match in the unisolated taus:
        int matchedTauInd=-1;
        TLorentzVector matchedEle;
        TLorentzVector unmatchedEle;
        for( size_t itau =0; itau< _Tau->size(); itau++){
          if(part2.DeltaR(_Tau->p4(itau))<0.3){
            //we are sure that part1 passes the tight id
            matchedTauInd=itau;
            matchedEle=part2;
            unmatchedEle=part1;
          }else if(part1.DeltaR(_Tau->p4(itau))<0.3){
            //check if part2 passes the tight id:
            if(find(active_part->at(CUTS::eRElec1)->begin(),active_part->at(CUTS::eRElec1)->end(),p2)!=active_part->at(CUTS::eRElec1)->end()){
              matchedTauInd=itau;
              matchedEle=part1;
              unmatchedEle=part2;
            }
          }
        }
        if(matchedTauInd>=0){
          if(find(active_part->at(CUTS::eRTau1)->begin(),active_part->at(CUTS::eRTau1)->end(),matchedTauInd)!=active_part->at(CUTS::eRTau1)->end()){
            histAddVal(_Tau->p4(matchedTauInd).Pt(), "DiEleGoodTauMatchPt");
            histAddVal(_Tau->p4(matchedTauInd).Pt()-matchedEle.Pt(), "DiEleGoodTauMatchDeltaPt");
            histAddVal((_Tau->p4(matchedTauInd)+unmatchedEle).M(), "DiEleGoodTauMatchMass");
            histAddVal((matchedEle+unmatchedEle).M(), "DiEleEleGoodMatchMass");
            histAddVal(matchedEle.Pt(), "DiEleEleGoodMatchPt");
            //histAddVal(_Tau->leadChargedCandPtError->at(matchedTauInd),"DiEleleadChargedCandPtErrorGoodMatched");
            //histAddVal(_Tau->leadChargedCandValidHits->at(matchedTauInd),"DiEleleadChargedCandValidHitGoodMatched");
            histAddVal2( matchedEle.Pt(),   (_Tau->p4(matchedTauInd).Pt()-matchedEle.Pt())/matchedEle.Pt(), "DiEleTauGoodMatchPt_vs_DeltaPt");
            histAddVal2( matchedEle.Pt(),   matchedEle.Eta(), "DiEleTauGoodMatchPt_vs_eta");
            histAddVal2( _Tau->pt(matchedTauInd),   _Tau->decayMode[matchedTauInd], "DiEleTauGoodMatchPt_vs_Decay");
          }else{
            histAddVal(_Tau->p4(matchedTauInd).Pt(), "DiEleTauMatchPt");
            histAddVal(_Tau->p4(matchedTauInd).Pt()-matchedEle.Pt(), "DiEleTauMatchDeltaPt");
            histAddVal((_Tau->p4(matchedTauInd)+unmatchedEle).M(), "DiEleTauMatchMass");
            histAddVal((matchedEle+unmatchedEle).M(), "DiEleEleMatchMass");
            histAddVal(matchedEle.Pt(), "DiEleEleMatchPt");
            //histAddVal(_Tau->leadChargedCandPtError->at(matchedTauInd),"DiEleleadChargedCandPtErrorMatched");
            //histAddVal(_Tau->leadChargedCandValidHits->at(matchedTauInd),"DiEleleadChargedCandValidHitsMatched");
            histAddVal2( matchedEle.Pt(),   (_Tau->p4(matchedTauInd).Pt()-matchedEle.Pt())/matchedEle.Pt(), "DiEleTauMatchPt_vs_DeltaPt");
            histAddVal2( matchedEle.Pt(),   matchedEle.Eta(), "DiEleTauMatchPt_vs_eta");
            histAddVal2( _Tau->pt(matchedTauInd),   _Tau->decayMode[matchedTauInd], "DiEleTauMatchPt_vs_Decay");
          }
        }else{
          histAddVal((part1+part2).M(), "DiEleEleUnMatchMass");
          histAddVal(part2.Pt(), "DiEleEleUnMatchPt");
          histAddVal2( part2.Pt(),   part2.Eta(), "DiEleUnMatchPt_vs_eta");
          if(!isData){
            histAddVal(part2.Pt(), "DiEleEleUnMatchPt_gen_"+std::to_string(abs(matchToGenPdg(part2,0.3))));
          }
          int found=-1;
          for(size_t i=0; i< _Jet->size(); i++) {
            if(part2.DeltaR(_Jet->p4(i)) <=0.4) {
              found=i;
            }
          }
          if (found>=0){
            //histAddVal(_Jet->chargedMultiplicity[found], "DiEleEleUnMatchJetMultiplicity");
          }else{
            histAddVal(-1, "DiEleEleUnMatchJetMultiplicity");
          }

        }
      }
    }
  }
