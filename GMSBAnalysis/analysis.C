#include "TChain.h"
#include "TFile.h"
#include "EventSelector.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

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
int parseCfgFile(EventSelector& evtProperties, string& cfgFileName, unsigned int run, unsigned int evt, unsigned int lumiSection)
{
  //initialize parameter variables
  unsigned int sampleType = 0;
  double ECALIsoMaxPTMultiplierEB = 0.0;
  double ECALIsoMaxConstantEB = 0.0;
  double ECALIsoMaxPTMultiplierEE = 0.0;
  double ECALIsoMaxConstantEE = 0.0;
  double HCALIsoMaxPTMultiplierEB = 0.0;
  double HCALIsoMaxConstantEB = 0.0;
  double HCALIsoMaxPTMultiplierEE = 0.0;
  double HCALIsoMaxConstantEE = 0.0;
  double HOverEMaxPresel = 0.0;
  double ETMin = 0.0;
  unsigned int fiducialRegion = 0;
  bool useHOverE = false;
  double HOverEMax = 0.0;
  bool useSigmaEtaEta = false;
  double sigmaEtaEtaMax = 0.0;
  bool useTrackIso = false;
  double trackIsoMaxPTMultiplier = 0.0;
  double trackIsoMaxConstant = 0.0;
  double trackPTMin = 0.0;
  double eTrackRMin = 0.0;
  double minDRPhotons = 0.0;
  bool useTimingCut = false;
  double maxSeedTime = 0.0;
  unsigned int numReqdCands = 2;
  string debugFileName = "";
  bool debugFlag = false;
  bool checkHaloCoincidenceWithPassingEBCandsOnly = false;
  bool rejectHalo = false;

  //read the file and import the parameter values
  ifstream in(cfgFileName.c_str());
  if (!in.is_open()) {
    cerr << "Error opening configuration file " << cfgFileName << ".\n";
    return;
  }
  sring line = "";
  while (!in.eof()) {
    getline(in, line);
    if (line[0] != '#') { //comment line
      size_t spacePos = 0;
      spacePos = line.find(' ');
      if (spacePos != string::npos) {
	stringstream paramName;
	paramName << line.substr(0, spacePos);
	
      }
    }
  }
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

  //clone the tree
  TTree* outputTree = chain->CloneTree(0);

  //set branch addresses


  //parse cfg file
  EventSelector evtProperties;
  int haloPars = parseCfgFile(evtProperties, cfgFileName, 0, 0, 0);
}
