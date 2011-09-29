// adapted from SusyAnalysis/SusyNtuplizer/macro/ana.C

#include "../../../GMSBTools/Filters/interface/Typedefs.h"

void ana_standalone() {

  gSystem->Load("libSusy.so");

  // // Look ../jec/JetMETObjects/README
  gSystem->Load("../jec/lib/libJetMETObjects.so");

  //GMSBTools shared library (contains Categorizer class)
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");

  // Printing utility for ntuple variables
  gROOT->LoadMacro("SusyEventPrinter.cc++");

  // Main analysis code
  gROOT->LoadMacro("EventAnalyzer.cc++");
  gROOT->LoadMacro("DongwookCategoryProducer.cc++");

  //configuration
  ParameterSet pars;
  pars.photonTag = "photons";
  pars.photon1ETMin = 45.0;
  pars.photon2ETMin = 30.0;
  pars.photonAbsEtaMax = 1.4442;
  pars.photonECALIsoMaxPTMultiplier = 0.006;
  pars.photonECALIsoMaxConstant = 4.2;
  // pars.photonECALIsoEffArea = 0.2918;
  pars.photonECALIsoEffArea = 0.0;
  pars.photonHCALIsoMaxPTMultiplier = 0.0025;
  pars.photonHCALIsoMaxConstant = 2.2;
  // pars.photonHCALIsoEffArea = 0.1039;
  pars.photonHCALIsoEffArea = 0.0;
  pars.photonHOverEMax = 0.05;
  pars.photonR9Max = 0.98;
  pars.photonTrackIsoMaxPTMultiplier = 0.001;
  pars.photonTrackIsoMaxPTConstant = 2.0;
  pars.photonSigmaIetaIetaMax = 0.011;
  pars.photonAbsSeedTimeMax = -1.0;
  pars.photonE2OverE9Max = -1.0;
  pars.photonDPhiMin = 0.05;
  pars.photonDRMin = 0.8;
  pars.pixelVetoOnFake = true;
  pars.treeName = "susyTree";
  pars.input = VSTRING(1, "~/RA3/data/Data2011A_ToRun167913_Filter_NoPileupCorr_Photon_NEW.root");
  pars.HLT = vector<TString>();
  pars.HLT.push_back("HLT_Photon32_CaloIdL_Photon26_CaloIdL");
  pars.HLT.push_back("HLT_Photon36_CaloIdL_Photon22_CaloIdL");
  pars.HLT.push_back("HLT_Photon40_CaloIdL_Photon28_CaloIdL");
  pars.nEvts = -1;
  pars.JSON = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON.txt";
  // pars.JSON = "";
  // pars.outputFile = "/Users/rachelyohay/RA3/data/debug.root";
  pars.outputFile = "/Users/rachelyohay/RA3/data/debug.root";

  TStopwatch ts;

  ts.Start();

  DongwookCategoryProducer producer(pars);

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
