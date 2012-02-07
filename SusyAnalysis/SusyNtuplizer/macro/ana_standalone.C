// adapted from SusyAnalysis/SusyNtuplizer/macro/ana.C

#include "EventAnalyzer.h"

void ana_standalone() {

  gSystem->Load("libSusy.so");

  // // Look ../jec/JetMETObjects/README
  gSystem->Load("../jec/lib/libJetMETObjects.so");

  //GMSBTools shared library (contains Categorizer class)
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");

  // Printing utility for ntuple variables
  gSystem->Load("SusyEventPrinter_cc.so");
//   gROOT->LoadMacro("SusyEventPrinter.cc++");

  // Main analysis code
  gSystem->SetIncludePath("-I../../..");
//   gROOT->LoadMacro("EventAnalyzer.cc++");
//   gROOT->LoadMacro("DongwookCategoryProducer.cc++");
  gSystem->Load("EventAnalyzer_cc.so");
  gSystem->Load("DongwookCategoryProducer_cc.so");

  /*now, ECAL/HCAL/track isolation cuts are meant to be matched with loose trigger cuts, so don't 
    use them (i.e. track isolation) to distinguish photon from fake!*/

  //configuration
  ParameterSet pars;
  pars.photonTag = "photons";
  pars.photon1ETMin = 40.0;
  pars.photon2ETMin = 25.0;
  pars.photonAbsEtaMax = 1.4442;
  pars.photonECALIsoMaxPTMultiplier = 0.012; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					       EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
					       cone*/
  pars.photonECALIsoMaxConstant = 6.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					 EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonECALIsoEffArea = 0.093; /*https://twiki.cern.ch/twiki/bin/view/CMS/
				      RA3IsolationConePileupCorrections, dR = 0.3 cone, 
				      18-Jan-12, default rho_EtaMax*/
  pars.photonHCALIsoMaxPTMultiplier = 0.005; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					       EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
					       cone*/
  pars.photonHCALIsoMaxConstant = 4.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					 EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonHCALIsoEffArea = 0.0281; /*https://twiki.cern.ch/twiki/bin/view/CMS/
					RA3IsolationConePileupCorrections, dR = 0.3 cone, 
					18-Jan-12, default rho_EtaMax*/
  pars.photonHOverEMax = 0.05;
  pars.photonR9Max = 1.0;
  pars.photonR9Min = 0.8;
  pars.photonTrackIsoMaxPTMultiplier = 0.002; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
						EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
						cone*/
  pars.photonTrackIsoMaxConstant = 4.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					  EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonCombinedIsoMax = 6.0;
  pars.fakeCombinedIsoMax = 20.0;
  pars.isoConeHLT = DR03;
  pars.isoConeOffline = DR03;
  pars.photonSigmaIetaIetaMax = 0.011;
  pars.photonHLTSigmaIetaIetaMax = 0.014;
  pars.photonAbsSeedTimeMax = -1.0;
  pars.photonE2OverE9Max = -1.0;
  pars.photonDPhiMin = 0.05;
  pars.photonDRMin = 0.6;
  pars.pixelVetoOnFake = true;
  pars.treeName = "susyTree";
  pars.input = VSTRING(1, "/data2/yohay/RA3/Data2011_Filter_JsonHLTnVertexTwo40-25GeVBarrelPhosWithR9-HoverE-CaloIdL-R9idOrIsoVL_NoPileupCorr_Photon_AnalysisSkim.root");
  pars.HLT = vector<TString>();
  pars.HLT.push_back("HLT_Photon26_IsoVL_Photon18");
  pars.HLT.push_back("HLT_Photon36_CaloIdL_Photon22_CaloIdL");
  pars.HLT.push_back("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL");
  pars.HLT.push_back("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id");
  pars.HLT.push_back("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL");
  pars.HLT.push_back("HLT_Photon36_R9Id_Photon22_R9Id");
  pars.nEvts = -1;
  pars.JSON = "/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/JSON_toRun180252.txt";
//   pars.JSON = "";
  pars.outputFile = "/data2/yohay/RA3/4684pb-1_categorized_18-Jan-12_skim.root";
  pars.recategorize = false;

  TStopwatch ts;

  ts.Start();

  DongwookCategoryProducer producer(pars);

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
