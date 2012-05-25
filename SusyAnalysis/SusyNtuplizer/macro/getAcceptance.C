#include "../../../GMSBTools/Filters/interface/Typedefs.h"

void getAcceptance_JOB()
{
  //load source code and precompiled shared libraries
  gSystem->Load("libRooFit.so");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("../jec/lib/libJetMETObjects.so");
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");
  gSystem->Load("../../../PhysicsTools/lib/libUtilities.so");
  gSystem->Load("../../../PhysicsTools/lib/libTagAndProbe.so");
  gSystem->Load("libSusy.so");
  gSystem->Load("SusyEventPrinter_cc.so");
//   gROOT->LoadMacro("SusyEventPrinter.cc++");
  gSystem->SetIncludePath("-I../../.. -I/uscmst1/prod/sw/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms12/include");
//   gROOT->LoadMacro("GMSBAnalyzer.C++");
  gSystem->Load("GMSBAnalyzer_C.so");

  //instantiate GMSBAnalyzer object
  vector<string> input;
  input.push_back("FILE");
  TChain* tree = new TChain("susyTree");
  for (VSTRING_IT iIn = input.begin(); iIn != input.end(); ++iIn) { tree->Add((*iIn).c_str()); }
  GMSBAnalyzer analyzer(tree);
  analyzer.setTag("photons");
  analyzer.setNEvts(-1);
  analyzer.setPUFile("160404-180252_PU.root");
  analyzer.initPU("Fall11Obs");
  analyzer.addHLT(TString("HLT_Photon26_IsoVL_Photon18"), 0, 0, -1);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_Photon22_CaloIdL"), 0, 0, -1);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL"), 0, 0, -1);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id"), 3, 0, -1);
  analyzer.addHLT(TString("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL"), 3, 0, -1);
  analyzer.addHLT(TString("HLT_Photon36_R9Id_Photon22_R9Id"), 3, 0, -1);
  analyzer.setPUReweightingFlag(false);
  analyzer.setHLTFlag(false);

  //loop over events
  TStopwatch ts;
  ts.Start();
  analyzer.makeAcceptancePlots("/uscms_data/d2/rpyohay/SCAN/acceptance_LABEL.root");
  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

  //clean up
  delete tree;
}
