#include "../../../GMSBTools/Filters/interface/Typedefs.h"

void GMSBAnalysis()
{
  //load source code and precompiled shared libraries
  gSystem->Load("libSusy.so");
  gSystem->Load("../jec/lib/libJetMETObjects.so");
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");
  gROOT->LoadMacro("SusyEventPrinter.cc+");
  gSystem->SetIncludePath("-I../../..");
  gROOT->LoadMacro("GMSBAnalyzer.C++");

  //instantiate GMSBAnalyzer object
  vector<string> input(1, "~/RA3/data/1140pb-1_categorized.root");
  TChain* tree = new TChain("susyTree");
  for (VSTRING_IT iIn = input.begin(); iIn != input.end(); ++iIn) { tree->Add((*iIn).c_str()); }
  GMSBAnalyzer analyzer(tree);
  analyzer.setTag("photons");
  analyzer.setNEvts(-1);

  //loop over events
  TStopwatch ts;
  ts.Start();
  analyzer.runEEVsFFAnalysis("/Users/rachelyohay/RA3/data/eeVsFFv4.root");
  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

  //clean up
  delete tree;
}
