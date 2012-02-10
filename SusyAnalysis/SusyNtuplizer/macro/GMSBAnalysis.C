#include "../../../GMSBTools/Filters/interface/Typedefs.h"

void GMSBAnalysis()
{
  //load source code and precompiled shared libraries
  gSystem->Load("libRooFit.so");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("libSusy.so");
  gSystem->Load("../jec/lib/libJetMETObjects.so");
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");
  gSystem->Load("../../../PhysicsTools/lib/libUtilities.so");
  gSystem->Load("../../../PhysicsTools/lib/libTagAndProbe.so");
//   gSystem->Load("SusyEventPrinter_cc.so");
  gROOT->LoadMacro("SusyEventPrinter.cc++");
  gSystem->SetIncludePath("-I../../.. -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms12/include");
  gROOT->LoadMacro("GMSBAnalyzer.C++");

  //instantiate GMSBAnalyzer object
  vector<string> input;
  input.push_back("/data2/yohay/RA3/4684pb-1_categorized_18-Jan-12_skim_v2.root");
  TChain* tree = new TChain("susyTree");
  for (VSTRING_IT iIn = input.begin(); iIn != input.end(); ++iIn) { tree->Add((*iIn).c_str()); }
  GMSBAnalyzer analyzer(tree);
  analyzer.setTag("photons");
  analyzer.setNEvts(10000);
  analyzer.setL1JECFile("/afs/cern.ch/user/y/yohay/data_CMSSW_4_2_8/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L1FastJet.txt");
  analyzer.setL2JECFile("/afs/cern.ch/user/y/yohay/data_CMSSW_4_2_8/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L2Relative.txt");
  analyzer.setL3JECFile("/afs/cern.ch/user/y/yohay/data_CMSSW_4_2_8/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L3Absolute.txt");
  analyzer.setJECErrFile("/afs/cern.ch/user/y/yohay/data_CMSSW_4_2_8/src/SusyAnalysis/SusyNtuplizer/jec/GR_R_42_V19_AK5PF_Uncertainty.txt");
  analyzer.addHLT(TString("HLT_Photon26_IsoVL_Photon18"), 0, 160404, 163261);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_Photon22_CaloIdL"), 0, 161216, 166967);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL"), 0, 166347, 180252);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id"), 3, 166347, 180252);
  analyzer.addHLT(TString("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL"), 3, 166347, 180252);
  analyzer.addHLT(TString("HLT_Photon36_R9Id_Photon22_R9Id"), 3, 166347, 180252);

  //loop over events
  TStopwatch ts;
  ts.Start();
  analyzer.runMETAnalysis("/data2/yohay/RA3/4684pb-1_MET_18-Jan-12_skim_newJetCleaning_relCombIso02_dz1.root");
//   analyzer.runMETAnalysis("/data2/yohay/RA3/4684pb-1_MET_latestEffectiveAreas_MC_EM_enriched_QCD.root");
//   analyzer.runMETAnalysis("/data2/yohay/RA3/4684pb-1_MET_withJESSyst.root");
//   analyzer.runEMFractionAnalysis("/data2/yohay/RA3/4684pb-1_EMF_18-Jan-12.root");
//   analyzer.runMETAnalysis("/data2/yohay/RA3/debug.root");
//   analyzer.runMETAnalysisWithEEBackgroundFit("/data2/yohay/RA3/3558pb-1_MET_eeBkgFit.root");
  // analyzer.runMETAnalysisWithEEBackgroundFit("/Users/rachelyohay/RA3/data/debug.root");
//   analyzer.testFitting("/data2/yohay/RA3/3500pb-1_MET_eeBkgFit.root", "/data2/yohay/RA3/fit.root");
  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

  //clean up
  delete tree;
}
