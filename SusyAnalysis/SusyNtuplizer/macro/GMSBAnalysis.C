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
  gSystem->Load("SusyEventPrinter_cc.so");
//   gROOT->LoadMacro("SusyEventPrinter.cc++");
  gSystem->SetIncludePath("-I../../.. -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms7/include");
  gROOT->LoadMacro("GMSBAnalyzer.C++");

  //instantiate GMSBAnalyzer object
  vector<string> input;
  input.push_back("/data2/yohay/RA3/4684pb-1_categorized_18-Jan-12_skim_v2.root");
  TChain* tree = new TChain("susyTree");
  for (VSTRING_IT iIn = input.begin(); iIn != input.end(); ++iIn) { tree->Add((*iIn).c_str()); }
  GMSBAnalyzer analyzer(tree);
  analyzer.setTag("photons");
  analyzer.setNEvts(-1);
  analyzer.setL1JECFile("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L1FastJet.txt");
  analyzer.setL2JECFile("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L2Relative.txt");
  analyzer.setL3JECFile("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L3Absolute.txt");
  analyzer.setJECErrFile("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/jec/GR_R_42_V19_AK5PF_Uncertainty.txt");
  analyzer.addHLT(TString("HLT_Photon26_IsoVL_Photon18"), 0, 160404, 163261);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_Photon22_CaloIdL"), 0, 161216, 166967);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL"), 0, 166347, 180252);
  analyzer.addHLT(TString("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id"), 3, 166347, 180252);
  analyzer.addHLT(TString("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL"), 3, 166347, 180252);
  analyzer.addHLT(TString("HLT_Photon36_R9Id_Photon22_R9Id"), 3, 166347, 180252);
  analyzer.setDiEMInvMassFitPar(2, "muCB", 91.2); //Z signal Crystal Ball mean
  analyzer.setDiEMInvMassFitPar(2, "sigma", 1.8); //Z signal Crystal Ball width
  analyzer.setDiEMInvMassFitPar(2, "alphaCB", 1.063); //Z signal Crystal Ball power law turn-on
  analyzer.setDiEMInvMassFitPar(2, "n", 143.16); //Z signal Crystal Ball power
  analyzer.setDiEMInvMassFitPar(2, "alphaCMS", 97.0); /*mee background RooCMSShape location of 
							 erf slope*/
  analyzer.setDiEMInvMassFitPar(2, "beta", 0.0922); /*mee background RooCMSShape erf slope 
						       parameter*/
  analyzer.setDiEMInvMassFitPar(2, "gamma", 0.191); //mee background RooCMSShape exponent
  analyzer.setDiEMInvMassFitPar(2, "muCMS", 58.0); /*mee background RooCMSShape location of 
						      exponential*/
  analyzer.setDiEMInvMassFitParLowerBound(2, "m", 71.0); /*mass variable (no fixed value, only 
							    range)*/
  analyzer.setDiEMInvMassFitParLowerBound(2, "muCB", 86.2);
  analyzer.setDiEMInvMassFitParLowerBound(2, "sigma", 1.0);
  analyzer.setDiEMInvMassFitParLowerBound(2, "alphaCB", /*1.063*/0.0);
  analyzer.setDiEMInvMassFitParLowerBound(2, "n", /*143.16*/0.0);
  analyzer.setDiEMInvMassFitParLowerBound(2, "alphaCMS", /*97.0*/0.0);
  analyzer.setDiEMInvMassFitParLowerBound(2, "beta", /*0.0922*/0.0);
  analyzer.setDiEMInvMassFitParLowerBound(2, "gamma", /*0.191*/-0.2);
  analyzer.setDiEMInvMassFitParLowerBound(2, "muCMS", /*58.0*/0.0);
  analyzer.setDiEMInvMassFitParUpperBound(2, "m", 111.0);
  analyzer.setDiEMInvMassFitParUpperBound(2, "muCB", 96.2);
  analyzer.setDiEMInvMassFitParUpperBound(2, "sigma", 5.0);
  analyzer.setDiEMInvMassFitParUpperBound(2, "alphaCB", /*1.063*/5.0);
  analyzer.setDiEMInvMassFitParUpperBound(2, "n", /*143.16*/200.0);
  analyzer.setDiEMInvMassFitParUpperBound(2, "alphaCMS", /*97.0*/200.0);
  analyzer.setDiEMInvMassFitParUpperBound(2, "beta", /*0.0922*/0.2);
  analyzer.setDiEMInvMassFitParUpperBound(2, "gamma", /*0.191*/0.2);
  analyzer.setDiEMInvMassFitParUpperBound(2, "muCMS", /*58.0*/200.0);
  analyzer.setDiEMInvMassFitPar(1, "muCB", 91.2);
  analyzer.setDiEMInvMassFitPar(1, "sigma", 1.8);
  analyzer.setDiEMInvMassFitPar(1, "alphaCB", 1.063);
  analyzer.setDiEMInvMassFitPar(1, "n", 143.16);
  analyzer.setDiEMInvMassFitPar(1, "alphaCMS", 72.02);
  analyzer.setDiEMInvMassFitPar(1, "beta", 0.098);
  analyzer.setDiEMInvMassFitPar(1, "gamma", 0.0375);
  analyzer.setDiEMInvMassFitPar(1, "muCMS", 56.0);
  analyzer.setDiEMInvMassFitParLowerBound(1, "m", 71.0);
  analyzer.setDiEMInvMassFitParLowerBound(1, "muCB", 86.2);
  analyzer.setDiEMInvMassFitParLowerBound(1, "sigma", 1.0);
  analyzer.setDiEMInvMassFitParLowerBound(1, "alphaCB", /*1.063*/0.0);
  analyzer.setDiEMInvMassFitParLowerBound(1, "n", /*143.16*/0.0);
  analyzer.setDiEMInvMassFitParLowerBound(1, "alphaCMS", /*72.02*/0.0);
  analyzer.setDiEMInvMassFitParLowerBound(1, "beta", /*0.098*/0.0);
  analyzer.setDiEMInvMassFitParLowerBound(1, "gamma", /*0.0375*/-0.2);
  analyzer.setDiEMInvMassFitParLowerBound(1, "muCMS", /*56.0*/0.0);
  analyzer.setDiEMInvMassFitParUpperBound(1, "m", 111.0);
  analyzer.setDiEMInvMassFitParUpperBound(1, "muCB", 96.2);
  analyzer.setDiEMInvMassFitParUpperBound(1, "sigma", 5.0);
  analyzer.setDiEMInvMassFitParUpperBound(1, "alphaCB", /*1.063*/5.0);
  analyzer.setDiEMInvMassFitParUpperBound(1, "n", /*143.16*/200.0);
  analyzer.setDiEMInvMassFitParUpperBound(1, "alphaCMS", /*72.02*/200.0);
  analyzer.setDiEMInvMassFitParUpperBound(1, "beta", /*0.098*/0.2);
  analyzer.setDiEMInvMassFitParUpperBound(1, "gamma", /*0.0375*/0.2);
  analyzer.setDiEMInvMassFitParUpperBound(1, "muCMS", /*56.0*/200.0);
  analyzer.setDiEMInvMassFitParUnit("m", "GeV");
  analyzer.setDiEMInvMassFitParUnit("muCB", "GeV");
  analyzer.setDiEMInvMassFitParUnit("sigma", "GeV");
  analyzer.setDiEMInvMassFitParUnit("alphaCB", "");
  analyzer.setDiEMInvMassFitParUnit("n", "");
  analyzer.setDiEMInvMassFitParUnit("alphaCMS", "GeV");
  analyzer.setDiEMInvMassFitParUnit("beta", "");
  analyzer.setDiEMInvMassFitParUnit("gamma", "");
  analyzer.setDiEMInvMassFitParUnit("muCMS", "GeV");

  //loop over events
  TStopwatch ts;
  ts.Start();
  analyzer.
    runMETAnalysis("/data2/yohay/RA3/4684pb-1_MET_18-Jan-12_skim_finalEGMisIDRate_geq1j_nominalMETBins.root");
//   analyzer.runMETAnalysis("/data2/yohay/RA3/4684pb-1_MET_latestEffectiveAreas_MC_EM_enriched_QCD.root");
//   analyzer.runMETAnalysis("/data2/yohay/RA3/4684pb-1_MET_withJESSyst.root");
//   analyzer.runEMFractionAnalysis("/data2/yohay/RA3/4684pb-1_EMF_18-Jan-12_skim.root");
//   analyzer.runMETAnalysis("/data2/yohay/RA3/debug.root");
//   analyzer.runMETAnalysisWithEGBackgroundFit("/data2/yohay/RA3/3558pb-1_MET_eeBkgFit.root");
  // analyzer.runMETAnalysisWithEGBackgroundFit("/Users/rachelyohay/RA3/data/debug.root");
//   analyzer.testFitting("/data2/yohay/RA3/3500pb-1_MET_eeBkgFit.root", "/data2/yohay/RA3/fit.root");
  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

  //clean up
  delete tree;
}
