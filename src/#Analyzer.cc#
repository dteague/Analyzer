#include "Analyzer.h"
#define ival(x) static_cast<int>(x)
#define BIG_NUM 46340
/*For speed, uncomment these commands and the define statement in the h file.  It will lead to 
  unchecked bounds on values, but will improve speed */
//#define const
//#define at(x) operator[](x)
typedef vector<int>::iterator vec_iter;

//Filespace that has all of the .in files
const string FILESPACE = "PartDet/";
const string PUSPACE = "Pileup/";
//////////PUBLIC FUNCTIONS////////////////////

///Constructor
Analyzer::Analyzer(string infile, string outfile) : hPUmc(new TH1F("hPUmc", "hPUmc", 100, 0, 100)), hPUdata(new TH1F("hPUdata", "hPUdata", 100, 0, 100)), MetCov(2,2) {
  cout << "setup start" << endl;
  f = TFile::Open(infile.c_str());
  f->cd("TNT");
  BOOM = (TTree*)f->Get("TNT/BOOM");
  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "TOTAL EVENTS: " << nentries << std::endl;

  for(int i=0; i < nTrigReq; i++) {
    vector<int>* tmpi = new vector<int>();
    vector<string>* tmps = new vector<string>();
    trigPlace[i] = tmpi;
    trigName[i] = tmps;
  }

  setupGeneral(BOOM,infile);

  isData = distats["Run"].bmap.at("isData");
  CalculatePUSystematics = distats["Run"].bmap.at("CalculatePUSystematics");
  histo = Histogramer(1, FILESPACE+"Hist_entries.in", FILESPACE+"Cuts.in", outfile, isData);


  initializePileupInfo(distats["Run"].smap.at("MCHistos"), distats["Run"].smap.at("DataHistos"));

  //////need to initialize histo and get values for cut arrays

  cuts_per.resize(histo.get_cuts()->size());
  cuts_cumul.resize(histo.get_cuts()->size());

  if(!isData) {
    _Gen = new Generated(BOOM, FILESPACE + "Gen_info.in");
    genStat = _Gen->pstats["Gen"];
    genMap = genStat.dmap;
  }
  _Electron = new Electron(BOOM, FILESPACE + "Electron_info.in");
  _Muon = new Muon(BOOM, FILESPACE + "Muon_info.in");
  _Tau = new Taus(BOOM, FILESPACE + "Tau_info.in");
  _Jet = new Jet(BOOM, FILESPACE + "Jet_info.in");

  std::cout << "setup complete" << std::endl << endl;

}

////destructor
Analyzer::~Analyzer() {
  delete f;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
  if(!isData) delete _Gen;
  
  for(int i=0; i < nTrigReq; i++) {
    delete trigPlace[i];
    delete trigName[i];
  }
}


///resets values so analysis can start
void Analyzer::clear_values() {
  for(int i=0; i < (int)goodParts.size(); i++) {
    goodParts[i].clear();
  }
  deltaMEx=0;
  deltaMEy=0;
  sumpxForMht = 0.0;
  sumpyForMht = 0.0;
  sumptForHt  = 0.0;
  leadIndex=-1;
  maxCut = 0;
}



///Function that does most of the work.  Calculates the number of each particle
void Analyzer::preprocess(int event) {
  BOOM->GetEntry(event);

  //TODO: add in pdf vector(set to 1 for now);
  
  theMETVector.SetPxPyPzE(Met_px, Met_py, Met_pz, sqrt(pow(Met_px,2) + pow(Met_py,2)));
  pu_weight = (!isData && CalculatePUSystematics) ? getPileupWeight(nTruePU) : 1.0;


  //////Setting up vectors for particles.  Put functions after this
  smearLepton(*_Electron, CUTS::eGElec, _Electron->pstats["Smear"]);
  smearLepton(*_Muon, CUTS::eGMuon, _Muon->pstats["Smear"]);
  smearLepton(*_Tau, CUTS::eGTau, _Tau->pstats["Smear"]);
  smearJet(_Jet->pstats["Smear"]);

  //////Triggers and Vertices
  goodParts[ival(CUTS::eRVertex)].resize(bestVertices);
  TriggerCuts(*(trigPlace[0]), *(trigName[0]), CUTS::eRTrig1);
  TriggerCuts(*(trigPlace[1]), *(trigName[1]), CUTS::eRTrig2);

  ////Updates Met and does MET cut
  updateMet();

  ////////////////////put new functions here to fill up the goodParts array for 
  ////// making cuts



  if(event % 50000 == 0) {
    cout << "Event #" << event << endl;
  }
}


////Reads cuts from Cuts.in file and see if the event has enough particles
void Analyzer::fillCuts() {
  unordered_map<string,pair<int,int> >* cut_info = histo.get_cuts();
  vector<string>* cut_order = histo.get_order();

  string cut;
  int min, max;
  bool prevTrue = true;
  int nparticles, i=0;
  maxCut=0;

  

  for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
    if(isData && it->find("Gen") != string::npos) continue;
    cut = *it;
    min= cut_info->at(cut).first;
    max= cut_info->at(cut).second;
    nparticles = goodParts[ival(cut_num[cut])].size();
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      if((cut_num[cut] == CUTS::eR1stJet || cut_num[cut] == CUTS::eR2ndJet) && goodParts[ival(cut_num[cut])].at(0) == -1 ) {
	prevTrue = false;
	continue;  ////dirty dirty hack
      }
      cuts_per[i]++;
      cuts_cumul[i] += (prevTrue) ? 1 : 0;
      maxCut += (prevTrue) ? 1 : 0;
    } else prevTrue = false;
  }

}


///Prints the number of events that passed each cut per event and cumulatively
//done at the end of the analysis
void Analyzer::printCuts() {
  vector<string>* cut_order = histo.get_order();
  int i =0;

  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(3);
  cout << "\n";
  cout << "Selection Efficiency " << "\n";
  cout << "Total events: " << nentries << "\n";
  cout << "               Name                 Indiv.         Cumulative\n";
  cout << "---------------------------------------------------------------------------\n";
  for(vector<string>::iterator it=cut_order->begin(); it != cut_order->end(); it++, i++) {
    cout << setw(28) << *it << "    ";
    if(isData && it->find("Gen") != string::npos) cout << "Skipped" << endl;
    else cout << setw(5) << cuts_per.at(i) << "  ( " << setw(5) << ((float)cuts_per.at(i)) / nentries << ") "
	      << setw(5) << cuts_cumul.at(i) << "  ( " << setw(5) << ((float)cuts_cumul.at(i)) / nentries << ") " << endl;
  }
  cout << "---------------------------------------------------------------------------\n";  
  histo.fill_histogram();
}

/////////////PRIVATE FUNCTIONS////////////////


/// For combinations of particles, use this convention.  Makes your life very easy for filling histograms
      ///Particlesp that lead to good combo are BIG_NUM * part1 + part2
      /// final / BIG_NUM = part1 (make sure is integer)
      /// final % BIG_NUM = part2 




////Grabs a list of the groups of histograms to be filled and asked Fill_folder to fill up the histograms
void Analyzer::fill_histogram() {
  if(distats["Run"].bmap["ApplyGenWeight"] && gen_weight == 0.0) return;
  fillCuts();
  vector<string> groups = *histo.get_groups();
  wgt = pu_weight;
  if(distats["Run"].bmap["ApplyGenWeight"]) wgt *= (gen_weight > 0) ? 1.0 : -1.0;

  for(vector<string>::iterator it = groups.begin(); it!=groups.end(); it++) {
    fill_Folder(*it, maxCut);
  }
}

///Function that fills up the histograms
void Analyzer::fill_Folder(string group, int max) {

  ////example of how to fill histogram.  don't change group, max or wgt

  //histo.addVal(Value, group, max, "Name", wgt);

}




//////-----------------DONT MESS WITH THESE FUNCTIONS-------------------/////

///Smears lepton only if specified and not a data file.  Otherwise, just filles up lorentz vectors
//of the data into the vector container smearP with is in each lepton object.
void Analyzer::smearLepton(Lepton& lepton, CUTS eGenPos, const PartStats& stats) {
  lepton.smearP.clear();

  for(int i = 0; i < (int)lepton.pt->size(); i++) {
    TLorentzVector tmpSmear;
    tmpSmear.SetPtEtaPhiE(lepton.pt->at(i), lepton.eta->at(i), lepton.phi->at(i), lepton.energy->at(i));
    lepton.smearP.push_back(tmpSmear);
  }
}

///Same as smearlepton, just jet specific
void Analyzer::smearJet(const PartStats& stats) {
  _Jet->smearP.clear();
  TLorentzVector jetV;

  for(int i=0; i< (int)_Jet->pt->size(); i++) {
    jetV.SetPtEtaPhiE(_Jet->pt->at(i), _Jet->eta->at(i), _Jet->phi->at(i), _Jet->energy->at(i));
    _Jet->smearP.push_back(jetV);
  }
}

///Calculates met from values from each file plus smearing and treating muons as neutrinos
void Analyzer::updateMet() {
  ////// Neutrino update before calculation
  if(distats["Run"].bmap.at("TreatMuonsAsNeutrinos")) {
    for(vec_iter it=goodParts[ival(CUTS::eRMuon1)].begin(); it!=goodParts[ival(CUTS::eRMuon1)].end(); it++) {
      if(find(goodParts[ival(CUTS::eRMuon2)].begin(), goodParts[ival(CUTS::eRMuon2)].end(), (*it)) != goodParts[ival(CUTS::eRMuon2)].end() ) continue;
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }    
    for(vec_iter it=goodParts[ival(CUTS::eRMuon2)].begin(); it!=goodParts[ival(CUTS::eRMuon2)].end(); it++) {
      deltaMEx += _Muon->smearP.at(*it).Px();
      deltaMEy += _Muon->smearP.at(*it).Py();
    }
  }
  ///---MHT and HT calculations----////
  int i=0;
  for(vector<TLorentzVector>::iterator it=_Jet->smearP.begin(); it!=_Jet->smearP.end(); it++, i++) {
    if( (it->Pt() > distats["Run"].dmap.at("JetPtForMhtAndHt")) && (fabs(it->Eta()) < distats["Run"].dmap.at("JetEtaForMhtAndHt")) ) {
      if(distats["Run"].bmap.at("ApplyJetLooseIDforMhtAndHt") && !passedLooseJetID(i) ) continue;
      
      sumpxForMht -= it->Px();
      sumpyForMht -= it->Py();
      sumptForHt  += it->Pt();
    }
  }
  phiForMht = atan2(sumpyForMht,sumpxForMht);
  //  if(sumpxForMht < 0) phiForMht += (sumpyForMht >= 0) ? TMath::Pi() : -TMath::Pi();

  theMETVector.SetPxPyPzE(theMETVector.Px()+deltaMEx, theMETVector.Py()+deltaMEy, theMETVector.Pz(), 
  			  TMath::Sqrt(pow(theMETVector.Px()+deltaMEx,2) + pow(theMETVector.Py()+deltaMEy,2)));

  /////MET CUTS

  if(distats["Run"].bmap.at("DiscrByMet")) {
    if(theMETVector.Pt() < distats["Run"].pmap.at("RecoMetCut").first) return;
    if(theMETVector.Pt() > distats["Run"].pmap.at("RecoMetCut").second) return;
  }
  if(distats["Run"].bmap.at("DiscrByMHT") && sqrt(pow(sumpxForMht,2) + pow(sumpyForMht,2)) < distats["Run"].dmap.at("MhtCut")) return;

  if(distats["Run"].bmap.at("DiscrByHT") && sumptForHt < distats["Run"].dmap.at("HtCut")) return; 
  
  goodParts[ival(CUTS::eMET)].push_back(1);
}



/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral(TTree* BOOM, string infile) {
  BOOM->SetBranchStatus("Trigger_decision", 1);
  BOOM->SetBranchStatus("Trigger_names", 1);
  BOOM->SetBranchStatus("nTruePUInteractions", 1);
  BOOM->SetBranchStatus("bestVertices", 1);
  BOOM->SetBranchStatus("weightevt", 1);
  BOOM->SetBranchStatus("Met_type1PF_px", 1);
  BOOM->SetBranchStatus("Met_type1PF_py", 1);
  BOOM->SetBranchStatus("Met_type1PF_pz", 1);
  BOOM->SetBranchStatus("Met_type1PF_cov00", 1);
  BOOM->SetBranchStatus("Met_type1PF_cov01", 1);
  BOOM->SetBranchStatus("Met_type1PF_cov10", 1);
  BOOM->SetBranchStatus("Met_type1PF_cov11", 1);

  BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision);
  BOOM->SetBranchAddress("Trigger_names", &Trigger_names);
  BOOM->SetBranchAddress("nTruePUInteractions", &nTruePU);
  BOOM->SetBranchAddress("bestVertices", &bestVertices);
  BOOM->SetBranchAddress("weightevt", &gen_weight);
  BOOM->SetBranchAddress("Met_type1PF_px", &Met_px);
  BOOM->SetBranchAddress("Met_type1PF_py", &Met_py);
  BOOM->SetBranchAddress("Met_type1PF_pz", &Met_pz);
  BOOM->SetBranchAddress("Met_type1PF_cov00", &MetCov[0][0]);
  BOOM->SetBranchAddress("Met_type1PF_cov01", &MetCov[0][1]);
  BOOM->SetBranchAddress("Met_type1PF_cov10", &MetCov[1][0]);
  BOOM->SetBranchAddress("Met_type1PF_cov11", &MetCov[1][1]);

  read_info(FILESPACE + "ElectronTau_info.in");
  read_info(FILESPACE + "MuonTau_info.in");
  read_info(FILESPACE + "MuonElectron_info.in");
  read_info(FILESPACE + "DiParticle_info.in");
  read_info(FILESPACE + "VBFCuts_info.in");
  read_info(FILESPACE + "Run_info.in");
}


///parsing method that gets info on diparts and basic run info
//put in map called "distats"
void Analyzer::read_info(string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }

  vector<string> stemp;
  string group, line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    stemp.clear();
    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      cout << "error in " << filename << "; no groups specified for data" << endl;
      exit(1);
    } else if(stemp.size() == 2) {
      if(stemp.at(0).find("Trigger") != string::npos) {
	int ntrig = (stemp.at(0).find("1") != string::npos) ? 0 : 1;
	trigName[ntrig]->push_back(stemp.at(1));
	trigPlace[ntrig]->push_back(0);
	continue;
      }
	
      char* p;
      strtod(stemp[1].c_str(), &p);
      if(stemp[1] == "1" || stemp[1] == "true") distats[group].bmap[stemp[0]]=true;
      else if(stemp[1] == "0" || stemp[1] == "false") distats[group].bmap[stemp[0]]=false; 
      else if(*p) distats[group].smap[stemp[0]] = stemp[1];
      else  distats[group].dmap[stemp[0]]=stod(stemp[1]);

    } else if(stemp.size() == 3) distats[group].pmap[stemp[0]] = make_pair(stod(stemp[1]), stod(stemp[2]));
    else if(stemp.size() == 4) {
      distats[group].smap["SVHistname"] = stemp[0];
      distats[group].dmap["SVbins"] = stod(stemp[1]);
      distats[group].dmap["SVmin"] = stod(stemp[2]);
      distats[group].dmap["SVmax"] = stod(stemp[3]);
    }
  }
  info_file.close();
}



void Analyzer::initializePileupInfo(string MCHisto, string DataHisto) {
  // Filenames must be c_strings below. Here is the conversion from strings to c_strings
  // As you can see above cstr1 corresponds to MC and cstr2 corresponds to data.

  TFile *file1 = new TFile((PUSPACE+MCHisto).c_str());
  TH1* histmc = static_cast<TH1*>(file1->Get("analyzeHiMassTau/NVertices_0"));
  if(!histmc) {throw std::runtime_error("failed to extract histogram");}
  for(int bin=0; bin<=(histmc->GetXaxis()->GetNbins() + 1); bin++) {
    hPUmc->SetBinContent(bin,histmc->GetBinContent(bin));
  }
  file1->Close();

  TFile* file2 = new TFile((PUSPACE+DataHisto).c_str());
  TH1* histdata = static_cast<TH1*>(file2->Get("analyzeHiMassTau/NVertices_0"));
  if(!histdata) {throw std::runtime_error("failed to extract histogram");}
  for(int bin=0; bin<=(histdata->GetXaxis()->GetNbins() + 1); bin++) {
    hPUdata->SetBinContent(bin,histdata->GetBinContent(bin));
  }
  file2->Close();
}

double Analyzer::getPileupWeight(float ntruePUInt) {
  int bin;
  double MCintegral;
  double MCvalue;
  double Dataintegral;
  double Datavalue;

  // The probability that data (or MC) has N pileup interactions is value / integral
  // The ratio of the data and MC probability density functions gives us our pileup weights

  //std::cout << "Grabbing pileup info. " << std::endl;
  bin = hPUmc->GetBin(ntruePUInt+1);
  MCvalue = hPUmc->GetBinContent(bin);
  MCintegral = hPUmc->Integral();
  Datavalue = hPUdata->GetBinContent(bin);
  Dataintegral = hPUdata->Integral();

  return ((MCvalue * Dataintegral) != 0) ? (Datavalue * MCintegral) / (MCvalue * Dataintegral) : 1.0;
}


///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(vec_iter it=goodParts[ival(ePos)].begin(); it < goodParts[ival(ePos)].end(); it++) {
    if(lvec.DeltaR(overlapper.smearP.at(*it)) < MatchingDeltaR) return true;
  }
  return false;
}

///Tests if tau decays into the specified number of jet prongs.
bool Analyzer::passProng(string prong, int value) {
  return ( (prong.find("1") != string::npos && value == 1) ||
	   (prong.find("2") != string::npos && value == 2) ||
	   (prong.find("3") != string::npos && value == 3) );
}


////Tests if tau is within the cracks of the detector (the specified eta ranges)
bool Analyzer::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 ||
          (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
          (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
          (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
          (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}
 

//Tests if a jet meets a litany of different tests
bool Analyzer::passedLooseJetID(int nobj) {
  if (_Jet->neutralHadEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->neutralEmEmEnergyFraction->at(nobj) >= 0.99) return false;
  if (_Jet->numberOfConstituents->at(nobj) <= 1) return false;
  if (_Jet->muonEnergyFraction->at(nobj) >= 0.80) return false;
  if ( (fabs(_Jet->smearP.at(nobj).Eta()) < 2.4) && 
           ((_Jet->chargedHadronEnergyFraction->at(nobj) <= 0.0) || 
	    (_Jet->chargedMultiplicity->at(nobj) <= 0.0) || 
	    (_Jet->chargedEmEnergyFraction->at(nobj) >= 0.99) )) return false;
  return true;
}


///sees if the event passed one of the two cuts provided
void Analyzer::TriggerCuts(vector<int>& prevTrig, const vector<string>& trigvec, CUTS ePos) {
  if(! need_cut[ePos]) return;
  for(int i = 0; i < (int)trigvec.size(); i++) {
    if(prevTrig[i] >= (int)Trigger_names->size() || trigvec.at(i) != Trigger_names->at(prevTrig.at(i)) ) {
      for(int j = 0; j < (int)Trigger_names->size(); j++) {
	if(Trigger_names->at(j).find(trigvec.at(i)) != string::npos) {
	  prevTrig.at(i) = j;
	  break;
	}
      }
    }
    if(prevTrig.at(i) < (int)Trigger_names->size() && Trigger_decision->at(prevTrig.at(i)) == 1) {
      goodParts[ival(ePos)].push_back(0);
      return;
    }
  }
}


///////Only tested for if is Zdecay, can include massptasymmpair later?
/////Tests to see if a light lepton came form a zdecay
bool Analyzer::isZdecay(const TLorentzVector& theObject, const Lepton& lep) {
  bool eventIsZdecay = false;
  const float zMass = 90.1876;
  const float zWidth = 2.4952;
  float zmmPtAsymmetry = -10.;

  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.
  for(vector<TLorentzVector>::const_iterator lepit= lep.smearP.begin(); lepit != lep.smearP.end(); lepit++) {
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


///Calculates the Pzeta value
double Analyzer::getPZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double px = visPx + theMETVector.Px();
  double py = visPy + theMETVector.Py();
  double pZeta = px*zetaX + py*zetaY;
  return pZeta;
}


///calculates the visible pzeta value
double Analyzer::getPZetaVis(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double pZetaVis = visPx*zetaX + visPy*zetaY;
  return pZetaVis;
}

///Normalizes phi to be between -PI and PI
double Analyzer::normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}

///Takes the absolute value of of normPhi (made because constant use)
double Analyzer::absnormPhi(double phi) {
  return abs(normPhi(phi));
}




