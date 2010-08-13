#include "TFile.h"
#include "EventSelector.h"
#include <vector>
#include <string>

//chain the input trees together
void makeChain(vector<string>& fileList, TChain* chain)
{
  for (vector<string>::const_iterator iFile = fileList.begin(), iFile != fileList.end(); ++iFile) {
    TFile* in = new TFile((*iFile).c_str());
    if (!in->IsOpen()) {
      cerr << "Error opening input file " << *iFile << ".\n";
      delete in;
      return;
    }
    TTree* tree;
    in->GetObject("gmsbAna/root", tree);
    chain->Add((*iFile).c_str());
  }
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
    if (line == "#") continue; //comment

    //extract parameter name, type, and value
    size_t spacePos = line.find(' ');
    if (spacePos == string::npos) {
      cerr << "Error in configuration file " << cfgFileName << ": improper formatting in line \"" << line << "\".\n";
      in.close();
      return -1;
    }
    stringstream paramName;
    paramName << line.substr(0, spacePos); //parameter name is the first word in the line
    size_t spacePos2 = line.find(' ', spacePos + 1);
    if (spacePos2 == string::npos) {
      cerr << "Error in configuration file " << cfgFileName << ": improper formatting in line \"" << line << "\".\n";
      in.close();
      return -1;
    }
    stringstream paramType;
    paramType << line.substr(spacePos + 1, spacePos2 - spacePos - 1); //parameter type is the second word in the line
    stringstream paramVal;
    if (paramType == "string") {
      if ((line.find('"') == (spacePos2 + 1)) && (line.find('"', spacePos2 + 2) == (line.length() - 1))) {
	paramVal << line.substr(spacePos2 + 2, line.length() - spacePos2 - 3); //to handle properly formatted strings
	if (stringParams.find(paramName) == stringParams.end()) stringParams[paramName] = paramVal;
	else cerr << "Ignoring duplicate value " << paramVal << " of parameter " << paramName << ".\n";
      }
      else {
	cerr << "Error in configuration file " << cfgFileName << ": improper formatting in line \"" << line << "\".\n";
	in.close();
	return -1;
      }
    }
    else paramVal << line.substr(spacePos2 + 1, line.length() - spacePos2 - 1); //parameter value is the third word in the line
    if (paramType == "unsigned_int") {
      if (uintParams.find(paramName) == uintParams.end()) uintParams[paramName] = paramVal;
      else cerr << "Ignoring duplicate value " << paramVal << " of parameter " << paramName << ".\n";
    }
    if (paramType == "double") {
      if (doubleParams.find(paramName) == doubleParams.end()) doubleParams[paramName] = paramVal;
      else cerr << "Ignoring duplicate value " << paramVal << " of parameter " << paramName << ".\n";
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
			 const CosInfoBranches& cosmicTrackInfo)
{
  //get/set the run, event, and lumi section numbers
  evtProperties.setRun(evtInfo.Run);
  evtProperties.setEvt(evtInfo.Event);
  evtProperties.setLumiSec(evtInfo.LumiBlk);
  evtProperties.printEvtInfo();

  //pass flag
  bool pass = false;

  //decide if the required number of passing objects were found
  vector<int> passingWithoutPixelSeed;
  vector<int> passingWithPixelSeed;
  bool allCandsFound = false;
  if (evtProperties.foundPhotonCandidates(photonInfo, passingWithoutPixelSeed, passingWithPixelSeed)) {

    //decide if a track was found
    if (evtProperties.getSampleType() == ETRACK) { if (evtProperties.foundTrack(trackInfo, photonInfo, passingWithPixelSeed)) allCandsFound = true; }

    //all samples besides ETRACK have all candidates found
    else allCandsFound = true;
  }

  //only check data quality on events with two passing candidates
  bool halo = false;
  if (allCandsFound && (haloParams["rejectHalo"] == 1)) {

    //create the collection of candidates to check for halo coincidence
    vector<int> passingCands;

    //only check coincidence with EB candidates that pass the selection...
    if (haloParams["checkHaloCoincidenceWithPassingEBCandsOnly"] == 1) {
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
    else for (unsigned int iPhoton = 0; iPhoton < photonInfo.Size; ++iPhoton) { passingCands.push_back(iPhoton); }

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
void runSampleMaker(string& outputFileName, vector<string>& fileList, string& cfgFileName, vector<int>& HLTBits)
{
  //chain the input trees together
  TChain* chain = new TChain("gmsbAna/root");
  makeChain(fileList, chain);

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

  //loop over the events
  unsigned int numPassingHLT = 0;
  unsigned int numPassingAll = 0;
  unsigned int numTot = chain->GetEntriesFast();
  vector<unsigned int> passingRun;
  vector<unsigned int> passingEvt;
  vector<unsigned int> passingLumiSection;
  for (unsigned int iEvt = 0; iEvt < numTot; ++iEvt) {
    cout << "Processing event #" << (iEvt + 1) << "...\n";

    //did it pass the HLT selection?
    if (passedHLT(evtInfo, HLTBits)) {
      ++numPassingHLT;

      //did it pass the full selection?
      if (passedFullSelection(evtProperties, haloParams, evtInfo, photonInfo, trackInfo, HEInfo, cosmicTrackInfo)) {
	++numPassingAll;
	passingRun.push_back(evtInfo.Run);
	passingEvt.push_back(evtInfo.Event);
	passingLumiSection.push_back(evtInfo.LumiBlk);

	//save the event to the output tree
	out->cd();
	outputTree->Fill();
      }
    }
  }

  //write the output file and exit
  out->cd();
  outputTree->Write();
  out->Write();
  out->Close();
  delete out;
  delete chain;
}
