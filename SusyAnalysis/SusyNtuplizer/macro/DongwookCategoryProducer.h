#ifndef SusyAnalysis_SusyNtuplizer_macro_DongwookCategoryProducer_h
#define SusyAnalysis_SusyNtuplizer_macro_DongwookCategoryProducer_h

#include <fstream>
#include "../../../GMSBTools/Filters/interface/Typedefs.h"
#include "TString.h"
#include "TChain.h"
#include "EventAnalyzer.h"
#include "../src/SusyCategory.h"

struct ParameterSet {

  TString photonTag;
  double photon1ETMin;
  double photon2ETMin;
  double photonAbsEtaMax;
  double photonECALIsoMaxPTMultiplier;
  double photonECALIsoMaxConstant;
  double photonECALIsoEffArea;
  double photonHCALIsoMaxPTMultiplier;
  double photonHCALIsoMaxConstant;
  double photonHCALIsoEffArea;
  double photonHOverEMax;
  double photonR9Max;
  double photonTrackIsoMaxPTMultiplier;
  double photonTrackIsoMaxConstant;
  double photonCombinedIsoMax;
  double fakeCombinedIsoMax;
  unsigned int isoConeHLT;
  unsigned int isoConeOffline;
  double photonSigmaIetaIetaMax;
  double photonAbsSeedTimeMax;
  double photonE2OverE9Max;
  double photonDPhiMin;
  double photonDRMin;
  bool pixelVetoOnFake;
  STRING treeName;
  VSTRING input;
  vector<TString> HLT;
  int nEvts;
  STRING JSON;
  STRING outputFile;

};

class DongwookCategoryProducer {

public:

  //constructor
  DongwookCategoryProducer(const ParameterSet&);

  //destructor
  ~DongwookCategoryProducer();

private:

  //input
  TString photonTag_;
  double photon1ETMin_;
  double photon2ETMin_;
  double photonAbsEtaMax_;
  double photonECALIsoMaxPTMultiplier_;
  double photonECALIsoMaxConstant_;
  double photonECALIsoEffArea_;
  double photonHCALIsoMaxPTMultiplier_;
  double photonHCALIsoMaxConstant_;
  double photonHCALIsoEffArea_;
  double photonHOverEMax_;
  double photonR9Max_;
  double photonTrackIsoMaxPTMultiplier_;
  double photonTrackIsoMaxConstant_;
  double photonCombinedIsoMax_;
  double fakeCombinedIsoMax_;
  unsigned int isoConeHLT_;
  unsigned int isoConeOffline_;
  double photonSigmaIetaIetaMax_;
  double photonAbsSeedTimeMax_;
  double photonE2OverE9Max_;
  double photonDPhiMin_;
  double photonDRMin_;
  bool pixelVetoOnFake_;
  STRING treeName_;
  VSTRING input_;
  vector<TString> HLT_;
  int nEvts_;
  STRING JSON_;
  STRING outputFile_;

  //tree access
  TChain* chain_;
  EventAnalyzer* treeReader_;
  TTree* outTree_;
  susy::Category* category_;

  //write selected events to file
  void writeEvents(const RUNEVTLUMIMAP&, ofstream&, const STRING&) const;
};

#endif
