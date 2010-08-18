#include "TFile.h"
#include "TH1F.h"
#include "TVector3.h"
#include "EventSelector.h"
#include <vector>
#include <string>
#include <sstream>

//chain the input trees together
void makeChain(vector<string>& fileList, TChain* chain)
{
  for (vector<string>::const_iterator iFile = fileList.begin(); iFile != fileList.end(); ++iFile) { chain->Add((*iFile).c_str()); }
}

//parse the cfg file
int parseCfgFile(EventSelector& evtProperties, map<string, unsigned int>& haloParams, string& cfgFileName)
{
  //open the cfg file
  ifstream in(cfgFileName.c_str());
  if (!in.is_open()) {
    cerr << "Error opening configuration file " << cfgFileName << ".\n";
    return -1;
  }

  //maps to hold the parameter name/value pairs
  map<string, unsigned int> uintParams;
  map<string, double> doubleParams;
  map<string, string> stringParams;

  //read in the parameter names and their values
  string line = "";
  while (!in.eof()) {
    getline(in, line);
    if (line[0] == '#') continue; //comment

    //extract parameter name, type, and value
    size_t spacePos = line.find(' ');
    if (spacePos == string::npos) {
      cerr << "Error in configuration file " << cfgFileName << ": improper formatting in line \"" << line << "\".\n";
      in.close();
      return -1;
    }
    stringstream paramName;
    paramName << line.substr(0, spacePos); //parameter name is the first word in the line
    string paramNameString = paramName.str();
    size_t spacePos2 = line.find(' ', spacePos + 1);
    if (spacePos2 == string::npos) {
      cerr << "Error in configuration file " << cfgFileName << ": improper formatting in line \"" << line << "\".\n";
      in.close();
      return -1;
    }
    stringstream paramType;
    paramType << line.substr(spacePos + 1, spacePos2 - spacePos - 1); //parameter type is the second word in the line
    string paramTypeString = paramType.str();
    stringstream paramVal;
    if (paramTypeString == "string") {
      if ((line.find('"') == (spacePos2 + 1)) && (line.find('"', spacePos2 + 2) == (line.length() - 1))) {
	paramVal << line.substr(spacePos2 + 2, line.length() - spacePos2 - 3); //to handle properly formatted strings
	string paramValString = paramVal.str();
	if (stringParams.find(paramNameString) == stringParams.end()) stringParams[paramNameString] = paramValString;
	else cerr << "Ignoring duplicate value " << paramValString << " of parameter " << paramNameString << ".\n";
      }
      else {
	cerr << "Error in configuration file " << cfgFileName << ": improper formatting in line \"" << line << "\".\n";
	in.close();
	return -1;
      }
    }
    else paramVal << line.substr(spacePos2 + 1, line.length() - spacePos2 - 1); //parameter value is the third word in the line
    if (paramTypeString == "unsigned_int") {
      unsigned int paramValUint;
      paramVal >> paramValUint;
      if (uintParams.find(paramNameString) == uintParams.end()) uintParams[paramNameString] = paramValUint;
      else cerr << "Ignoring duplicate value " << paramValUint << " of parameter " << paramNameString << ".\n";
    }
    if (paramType.str() == "double") {
      double paramValDouble;
      paramVal >> paramValDouble;
      if (doubleParams.find(paramNameString) == doubleParams.end()) doubleParams[paramNameString] = paramValDouble;
      else cerr << "Ignoring duplicate value " << paramValDouble << " of parameter " << paramNameString << ".\n";
    }
  }
  in.close();

  //create EventSelector object
  if ((uintParams.find("sampleType") != uintParams.end()) && (doubleParams.find("ECALIsoMaxPTMultiplierEB") != doubleParams.end()) && 
      (doubleParams.find("ECALIsoMaxConstantEB") != doubleParams.end()) && (doubleParams.find("ECALIsoMaxPTMultiplierEE") != doubleParams.end()) && 
      (doubleParams.find("ECALIsoMaxConstantEE") != doubleParams.end()) && (doubleParams.find("HCALIsoMaxPTMultiplierEB") != doubleParams.end()) && 
      (doubleParams.find("HCALIsoMaxConstantEB") != doubleParams.end()) && (doubleParams.find("HCALIsoMaxPTMultiplierEE") != doubleParams.end()) && 
      (doubleParams.find("HCALIsoMaxConstantEE") != doubleParams.end()) && (doubleParams.find("HOverEMaxPresel") != doubleParams.end()) && 
      (doubleParams.find("ETMin") != doubleParams.end()) && (uintParams.find("fiducialRegion") != uintParams.end()) && 
      (uintParams.find("useHOverE") != uintParams.end()) && (doubleParams.find("HOverEMax") != doubleParams.end()) && 
      (uintParams.find("useSigmaEtaEta") != uintParams.end()) && (doubleParams.find("sigmaEtaEtaMax") != doubleParams.end()) && 
      (uintParams.find("useTrackIso") != uintParams.end()) && (doubleParams.find("trackIsoMaxPTMultiplier") != doubleParams.end()) && 
      (doubleParams.find("trackIsoMaxConstant") != doubleParams.end()) && (doubleParams.find("trackPTMin") != doubleParams.end()) && 
      (doubleParams.find("eTrackRMin") != doubleParams.end()) && (doubleParams.find("minDRPhotons") != doubleParams.end()) && 
      (uintParams.find("useTimingCut") != uintParams.end()) && (doubleParams.find("maxSeedTime") != doubleParams.end()) && 
      (stringParams.find("debugFileName") != stringParams.end()) && (uintParams.find("debugFlag") != uintParams.end())) {
    unsigned int numReqdCands = 2;
    if (uintParams["sampleType"] == ETRACK) numReqdCands = 1;
    EventSelector theEvtProperties(uintParams["sampleType"], doubleParams["ECALIsoMaxPTMultiplierEB"], doubleParams["ECALIsoMaxConstantEB"], 
				   doubleParams["ECALIsoMaxPTMultiplierEE"], doubleParams["ECALIsoMaxConstantEE"], doubleParams["HCALIsoMaxPTMultiplierEB"], 
				   doubleParams["HCALIsoMaxConstantEB"], doubleParams["HCALIsoMaxPTMultiplierEE"], doubleParams["HCALIsoMaxConstantEE"], 
				   doubleParams["HOverEMaxPresel"], doubleParams["ETMin"], uintParams["fiducialRegion"], uintParams["useHOverE"], 
				   doubleParams["HOverEMax"], uintParams["useSigmaEtaEta"], doubleParams["sigmaEtaEtaMax"], uintParams["useTrackIso"], 
				   doubleParams["trackIsoMaxPTMultiplier"], doubleParams["trackIsoMaxConstant"], doubleParams["trackPTMin"], 
				   doubleParams["eTrackRMin"], doubleParams["minDRPhotons"], uintParams["useTimingCut"], doubleParams["maxSeedTime"], 
				   numReqdCands, 0, 0, 0, stringParams["debugFileName"], uintParams["debugFlag"]);
    evtProperties = theEvtProperties;
  }
  else {
    cerr << "Error in configuration file " << cfgFileName << ": not all necessary parameters supplied.\n";
    return -1;
  }

  //save halo parameters
  if ((uintParams.find("checkHaloCoincidenceWithPassingEBCandsOnly") != uintParams.end()) && (uintParams.find("rejectHalo") != uintParams.end())) {
    haloParams["checkHaloCoincidenceWithPassingEBCandsOnly"] = uintParams["checkHaloCoincidenceWithPassingEBCandsOnly"];
    haloParams["rejectHalo"] = uintParams["rejectHalo"];
  }
  else {
    cerr << "Error in configuration file " << cfgFileName << ": not all necessary parameters supplied.\n";
    return -1;
  }

  //exit
  return 0;
}

//determine if an event passed the HLT selection
bool passedHLT(EvtInfoBranches& evtInfo, vector<int>& HLTBits)
{
  for (vector<int>::const_iterator iHLTBit = HLTBits.begin(); iHLTBit != HLTBits.end(); ++iHLTBit) {
    if (evtInfo.TrgBook[*iHLTBit] == 1) return true;
  }
  return false;
}

//determine if a particular integer is found in a vector of integers
bool foundElemInVec(const vector<int>& vec, const int elem)
{
  vector<int>::const_iterator iVec = vec.begin();
  bool foundElem = false;
  while ((iVec != vec.end()) && (!foundElem)) {
    if (*iVec == elem) foundElem = true;
    ++iVec;
  }
  return foundElem;
}

//determine if an event passed the full selection
bool passedFullSelection(EventSelector& evtProperties, const map<string, unsigned int>& haloParams, const EvtInfoBranches& evtInfo, 
			 const PhoInfoBranches& photonInfo, const TrkInfoBranches& trackInfo, const HEHitInfoBranches& HEInfo, 
			 const CosInfoBranches& cosmicTrackInfo, vector<int>& passingWithoutPixelSeed, vector<int>& passingWithPixelSeed)
{
  //get/set the run, event, and lumi section numbers
  evtProperties.setRun(evtInfo.Run);
  evtProperties.setEvt(evtInfo.Event);
  evtProperties.setLumiSec(evtInfo.LumiBlk);
  evtProperties.printEvtInfo();

  //pass flag
  bool pass = false;

  //decide if the required number of passing objects were found
  bool allCandsFound = false;
  if (evtProperties.foundPhotonCandidates(photonInfo, passingWithoutPixelSeed, passingWithPixelSeed)) {

    //decide if a track was found
    if (evtProperties.getSampleType() == ETRACK) { if (evtProperties.foundTrack(trackInfo, photonInfo, passingWithPixelSeed)) allCandsFound = true; }

    //all samples besides ETRACK have all candidates found
    else allCandsFound = true;
  }

  //only check data quality on events with two passing candidates
  bool halo = false;
  if (allCandsFound && (haloParams.find("rejectHalo")->second == 1)) {

    //create the collection of candidates to check for halo coincidence
    vector<int> passingCands;

    //only check coincidence with EB candidates that pass the selection...
    if (haloParams.find("checkHaloCoincidenceWithPassingEBCandsOnly")->second == 1) {
      for (vector<int>::const_iterator iPhoton = passingWithoutPixelSeed.begin(); iPhoton != passingWithoutPixelSeed.end(); ++iPhoton) {

	//******************REMOVE THIS CHECK ONCE YOU VERIFY THAT IT'S UNNECCESSARY******************//
	if (!foundElemInVec(passingCands, *iPhoton)) passingCands.push_back(*iPhoton);
	else {
	  stringstream infoStream;
	  infoStream << "Already found an element with index " << *iPhoton << " in passingCands.  >1 element in passingWithoutPixelSeed has the same index!";
	  infoStream << endl;
	  evtProperties.printDebug(infoStream.str());
	}
	//********************************************************************************************//

      }
      for (vector<int>::const_iterator iElectron = passingWithPixelSeed.begin(); iElectron != passingWithPixelSeed.end(); ++iElectron) {

	//******************REMOVE THIS CHECK ONCE YOU VERIFY THAT IT'S UNNECCESSARY******************//
	if (!foundElemInVec(passingCands, *iElectron)) passingCands.push_back(*iElectron);
	else {
	  stringstream infoStream;
	  infoStream << "Already found an element with index " << *iElectron << " in passingCands.  An element in passingWithPixelSeed has the ";
	  infoStream << "same index as an element in passingWithoutPixelSeed!\n";
	  evtProperties.printDebug(infoStream.str());
	}
	//********************************************************************************************//

      }
    }

    //...or, check halo coincidence with all EB candidates
    else for (int iPhoton = 0; iPhoton < photonInfo.Size; ++iPhoton) { passingCands.push_back(iPhoton); }

    //calculate the halo tag
    halo = !evtProperties.passesDataQualityCuts(HEInfo, photonInfo, cosmicTrackInfo, passingCands);
  }

  //determine if the event passes
  stringstream infoStream;
  if (allCandsFound && (!halo)) {
    pass = true;
    infoStream << "Event passes -- run " << evtInfo.Run << ", event " << evtInfo.Event << ", lumi section " << evtInfo.LumiBlk;
    infoStream << ".\n---------------------------\n\n";
  }
  else infoStream << "Event fails.\n---------------------------\n\n";
  evtProperties.printDebug(infoStream.str());
  return pass;
}

//run the sample maker
void runSampleMaker(string outputFileName, vector<string>& fileList, string cfgFileName, vector<int>& HLTBits)
{
  //

  //chain the input trees together
  TChain* chain = new TChain("gmsbAna/root");
  //makeChain(fileList, chain);
  for (vector<string>::const_iterator iFile = fileList.begin(); iFile != fileList.end(); ++iFile) { chain->Add((*iFile).c_str()); }

  //open output file
  TFile* out = new TFile(outputFileName.c_str(), "RECREATE");
  if (!out->IsOpen()) {
    cerr << "Error opening output file " << outputFileName << ".\n";
    delete out;
    return;
  }
  out->cd();

  //set branch addresses
  EvtInfoBranches evtInfo;
  evtInfo.Register(chain);
  PhoInfoBranches photonInfo;
  photonInfo.Register(chain);
  TrkInfoBranches trackInfo;
  trackInfo.Register(chain);
  CosInfoBranches cosmicTrackInfo;
  cosmicTrackInfo.Register(chain);
  HEHitInfoBranches HEInfo;
  HEInfo.Register(chain);

  //clone the tree
  out->cd();
  TTree* outputTree = chain->CloneTree(0);

  //parse cfg file
  EventSelector evtProperties;
  map<string, unsigned int> haloParams;
  int success = parseCfgFile(evtProperties, haloParams, cfgFileName);
  if (success != 0) {
    cerr << "Error in parseCfgFile.\n";
    delete chain;
    out->Close();
    delete out;
    return;
  }

  //book photon histograms
  out->cd();
  out->mkdir("photon_histograms");
  out->cd("photon_histograms");
  TH1F* ECALIso = new TH1F("ECALIso", "ECAL isolation of passing e/#gamma candidates;ECAL isolation (GeV);e/#gamma candidates per GeV", 50, 0.0, 50.0);
  TH1F* HCALIso = new TH1F("HCALIso", "HCAL isolation of passing e/#gamma candidates;HCAL isolation (GeV);e/#gamma candidates per GeV", 50, 0.0, 50.0);
  TH1F* ET = new TH1F("ET", "E_{T} of passing e/#gamma candidates;E_{T} (GeV);e/#gamma candidates per GeV", 100, 0.0, 100.0);
  TH1F* HOverE = new TH1F("HOverE", "H/E of passing e/#gamma candidates;H/E;e/#gamma candidates per 0.01", 100, 0.0, 1.0);
  TH1F* fiducialRegionDist = new TH1F("fiducialRegion", "Fiducial region of passing e/#gamma candidates;;e/#gamma candidates per ECAL fiducial region", 2, 
				      0.5, 2.5);
  fiducialRegionDist->GetXaxis()->SetBinLabel(EB, "EB");
  fiducialRegionDist->GetXaxis()->SetBinLabel(EEND, "EE");
  TH1F* etaWidth = new TH1F("etaWidth", "#sigma_{#eta#eta} of passing e/#gamma candidates;#sigma_{#eta#eta};e/#gamma candidates per 0.001", 100, 0.0, 0.1);
  TH1F* trackIso = new TH1F("trackIso", "Track isolation of passing e/#gamma candidates;Track isolation (GeV);e/#gamma candidates per GeV", 50, 0.0, 50.0);
  TH1F* eta = new TH1F("eta", "#eta of passing e/#gamma candidates;#eta;e/#gamma candidates per 0.2", 15, -1.5, 1.5);
  TH1F* phi = new TH1F("phi", "#phi of passing e/#gamma candidates;#phi;e/#gamma candidates per 0.4", 16, -3.2, 3.2);
  TH1F* numCandsPerEvt = new TH1F("numCandsPerEvt", "Number of passing candidates per event;Number per event;Events per 1", 5, -0.5, 4.5);
  TH1F* diEMPT = new TH1F("diEMPT", "di-EM p_{T} of first two passing e/#gamma candidates;di-EM p_{T} (GeV);Events per GeV", 100, 0.0, 100.0);
  TH1F* dPhiCands = new TH1F("dPhiCands", "#Delta#phi between the 2 leading e/#gamma candidates;#Delta#phi (rad);Events/0.4", 8, 0.0, 3.2);
  TH1F* dPhiMETCand1 = new TH1F("dPhiMETCand1", "#Delta#phi between the ME_{T} and the leading e/#gamma candidate;#Delta#phi (rad);Events/0.4", 8, 0.0, 3.2);
  TH1F* dPhiMETCand2 = new TH1F("dPhiMETCand2", "#Delta#phi between the ME_{T} and the trailing e/#gamma candidate;#Delta#phi (rad);Events/0.4", 8, 0.0, 
				3.2);
  TH1F* dPhiMETCand = new TH1F("dPhiMETCand", "#Delta#phi between the ME_{T} and the e/#gamma candidate;#Delta#phi (rad);e/#gamma candidates/0.4", 8, 0.0, 
			       3.2);
  TH1F* dPhiMETDiEMPT = new TH1F("dPhiMETDiEMPT", "#Delta#phi between the ME_{T} and the di-EM p_{T} vector;#Delta#phi (rad);Events/0.4", 8, 0.0, 3.2);

  //book track histograms
  /*out->cd("..");
  out->mkdir("track_histograms");
  out->cd("track_histograms");
  TH1F* trackPT = new TH1F("trackPT", "p_{T} of passing tracks;p_{T} (GeV);Tracks per GeV", 100, 0.0, 100.0);
  TH1F* dRTrackPhoton = new TH1F("dRTrackPhoton", "#DeltaR(track, e/#gamma candidate) of passing tracks-e/#gamma candidate pairs;#DeltaR;", 100, 0.0, 1.0);
  dRTrackPhoton->GetYaxis()->SetTitle("track-e/#gamma candidate pairs per 0.01");*/

  //loop over the events
  unsigned int numPassingHLT = 0;
  unsigned int numPassingAll = 0;
  Long64_t bigNum = chain->GetEntriesFast();
  unsigned int numTot = 0;
  vector<unsigned int> passingRun;
  vector<unsigned int> passingEvt;
  vector<unsigned int> passingLumiSection;
  for (Long64_t iEvt = 0; iEvt < bigNum/*20000*/; ++iEvt) {
    if (chain->LoadTree(iEvt) < 0) {
      cout << "Event #" << (iEvt + 1) << " has no corresponding tree.  There must have been only " << iEvt << " events in the chain.\n";
      numTot = iEvt;
      break;
    }
    chain->GetEntry(iEvt);
    cout << "Processing event #" << (iEvt + 1) << "...\n";

    //did it pass the HLT selection?
    if (passedHLT(evtInfo, HLTBits)) {
      ++numPassingHLT;

      //did it pass the full selection?
      vector<int> passingWithoutPixelSeed;
      vector<int> passingWithPixelSeed;
      if (passedFullSelection(evtProperties, haloParams, evtInfo, photonInfo, trackInfo, HEInfo, cosmicTrackInfo, passingWithoutPixelSeed, 
			      passingWithPixelSeed)) {
	++numPassingAll;
	passingRun.push_back(evtInfo.Run);
	passingEvt.push_back(evtInfo.Event);
	passingLumiSection.push_back(evtInfo.LumiBlk);

	//fill plots of photon quantities by photon
	for (vector<int>::const_iterator iPhoton = passingWithoutPixelSeed.begin(); iPhoton != passingWithoutPixelSeed.end(); ++iPhoton) {
	  ECALIso->Fill(photonInfo.EcalIso[*iPhoton]);
	  HCALIso->Fill(photonInfo.HcalIso[*iPhoton]);
	  ET->Fill(photonInfo.Pt[*iPhoton]);
	  HOverE->Fill(photonInfo.HadOverEM[*iPhoton]);
	  fiducialRegionDist->Fill(evtProperties.ECALFiducialRegion(photonInfo, *iPhoton));
	  etaWidth->Fill(photonInfo.sigmaIetaIeta[*iPhoton]);
	  trackIso->Fill(photonInfo.TrackIsoPtHol[*iPhoton]);
	  eta->Fill(photonInfo.Eta[*iPhoton]);
	  phi->Fill(photonInfo.Phi[*iPhoton]);
	}

	//fill plot of number of passing candidates per event
	numCandsPerEvt->Fill(passingWithoutPixelSeed.size() + passingWithPixelSeed.size());

	//plots specific to gg and ff samples
	if ((evtProperties.getSampleType() == GAMMAGAMMA) || (evtProperties.getSampleType() == FF)) {

	  //identify the two leading candidates
	  int cand1 = 0;
	  float cand1PT = -1.0;
	  for (vector<int>::const_iterator iCand = passingWithoutPixelSeed.begin(); iCand != passingWithoutPixelSeed.end(); ++iCand) {
	    const float pT = photonInfo.Pt[*iCand];
	    if (pT > cand1PT) {
	      cand1 = *iCand;
	      cand1PT = pT;
	    }
	  }
	  int cand2 = 0;
	  float cand2PT = -1.0;
	  for (vector<int>::const_iterator iCand = passingWithoutPixelSeed.begin(); iCand != passingWithoutPixelSeed.end(); ++iCand) {
	    const float pT = photonInfo.Pt[*iCand];
	    if ((pT > cand2PT) && (*iCand != cand1)) {
	      cand2 = *iCand;
	      cand2PT = pT;
	    }
	  }

	  //di-EM pT
	  TVector3 p1(photonInfo.Px[cand1], photonInfo.Py[cand1], photonInfo.Pz[cand1]);
	  TVector3 p2(photonInfo.Px[cand2], photonInfo.Py[cand2], photonInfo.Pz[cand2]);
	  TVector3 diEMPTVec = p1 + p2;
	  diEMPT->Fill(diEMPTVec.Perp());

	  //dPhi
	  dPhiCands->Fill(evtProperties.dPhi(photonInfo.Phi[cand1], photonInfo.Phi[cand2]));
	  float dPhi1 = evtProperties.dPhi(photonInfo.Phi[cand1], evtInfo.METphi);
	  float dPhi2 = evtProperties.dPhi(photonInfo.Phi[cand2], evtInfo.METphi);
	  dPhiMETCand1->Fill(dPhi1);
	  dPhiMETCand2->Fill(dPhi2);
	  dPhiMETCand->Fill(dPhi1);
	  dPhiMETCand->Fill(dPhi2);
	  dPhiMETDiEMPT->Fill(evtProperties.dPhi(evtInfo.METphi, (float)diEMPTVec.Phi()));
	}

	//save the event to the output tree
	out->cd();
	outputTree->Fill();
      }
    }
  }

  //write the output file and exit
  cout << numPassingHLT << "/" << numTot << " events passed the HLT selection.\n";
  cout << numPassingAll << "/" << numTot << " events passed the full selection:\n";
  for (vector<unsigned int>::const_iterator iPassingEvt = passingEvt.begin(); iPassingEvt != passingEvt.end(); ++iPassingEvt) {
    const unsigned int index = iPassingEvt - passingEvt.begin();
    cout << "Run " << passingRun[index] << ", event " << passingEvt[index] << ", lumi section " << passingLumiSection[index] << endl;
  }
  out->cd();
  outputTree->Write();
  out->cd("photon_histograms");
  ECALIso->Write();
  HCALIso->Write();
  ET->Write();
  HOverE->Write();
  fiducialRegionDist->Write();
  etaWidth->Write();
  trackIso->Write();
  eta->Write();
  phi->Write();
  numCandsPerEvt->Write();
  diEMPT->Write();
  dPhiCands->Write();
  dPhiMETCand1->Write();
  dPhiMETCand2->Write();
  dPhiMETCand->Write();
  dPhiMETDiEMPT->Write();
  out->Write();
  out->Close();
  delete out;
  delete chain;
}
