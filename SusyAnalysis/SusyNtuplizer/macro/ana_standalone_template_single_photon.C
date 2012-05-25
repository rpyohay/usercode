// adapted from SusyAnalysis/SusyNtuplizer/macro/ana.C

#include "/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/EventAnalyzer.h"

void ana_standalone_single_photon_MC() {

  gSystem->Load("/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/libSusy.so");

  // // Look ../jec/JetMETObjects/README
  gSystem->Load("/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/jec/lib/libJetMETObjects.so");

  //GMSBTools shared library (contains Categorizer class)
  gSystem->Load("/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/GMSBTools/lib/libFilters.so");

  // Printing utility for ntuple variables
  gSystem->Load("/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/SusyEventPrinter_cc.so");

  // Main analysis code
  gSystem->SetIncludePath("-I/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src");
  gSystem->Load("/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/EventAnalyzer_cc.so");
  gSystem->Load("/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/DongwookCategoryProducer_cc.so");

  /*now, ECAL/HCAL/track isolation cuts are meant to be matched with loose trigger cuts, so don't 
    use them (i.e. track isolation) to distinguish photon from fake!*/

  //configuration
  ParameterSet pars;
  pars.photonTag = "photons";
  pars.photon1ETMin = 40.0;
  pars.photon2ETMin = 30.0;
  pars.photonAbsEtaMax = 1.4442;
  pars.photonECALIsoMaxPTMultiplier = 0.012; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					       EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
					       cone*/
  pars.photonECALIsoMaxConstant = 6.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					 EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonECALIsoEffArea = 0.1474; /*https://twiki.cern.ch/twiki/bin/view/CMS/
					RA3IsolationConePileupCorrections, dR = 0.3 cone, 
					18-Oct-11*/
  pars.photonHCALIsoMaxPTMultiplier = 0.005; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					       EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
					       cone*/
  pars.photonHCALIsoMaxConstant = 4.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					 EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonHCALIsoEffArea = 0.0467; /*https://twiki.cern.ch/twiki/bin/view/CMS/
					RA3IsolationConePileupCorrections, dR = 0.3 cone, 
					18-Oct-11*/
  pars.photonHOverEMax = 0.05;
  pars.photonR9Max = 0.98;
  pars.photonTrackIsoMaxPTMultiplier = 0.002; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
						EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
						cone*/
  pars.photonTrackIsoMaxConstant = 4.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					  EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonCombinedIsoMax = 6.0;
  pars.fakeCombinedIsoMax = 12.0;
  pars.isoConeHLT = DR03;
  pars.isoConeOffline = DR03;
  pars.photonSigmaIetaIetaMax = 0.011;
  pars.photonAbsSeedTimeMax = -1.0;
  pars.photonE2OverE9Max = -1.0;
  pars.photonDPhiMin = 0.05;
  pars.photonDRMin = 0.8;
  pars.pixelVetoOnFake = true;
  pars.treeName = "susyTree";
//   pars.input = VSTRING(1, "/data2/yohay/RA3/Data2011A_ToRun167913_Filter-JsonHLTtwo43-30GeVPhosWithR9HoverE_NoPileupCorr_Photon_NEW.root");
  pars.input = VSTRING();
FILES
  pars.HLT = vector<TString>();
  pars.HLT.push_back("HLT_Photon70_CaloIdL_HT300");
  pars.HLT.push_back("HLT_Photon70_CaloIdL_HT400");
  pars.HLT.push_back("HLT_Photon135");
  pars.HLT.push_back("HLT_Photon75_CaloIdVL_IsoL");
  pars.HLT.push_back("HLT_Photon125");
  pars.nEvts = -1;
  pars.JSON = "";
//   pars.JSON = "/afs/cern.ch/user/y/yohay/MC_CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/JSON_160431-177878_May10ReReco_Run2011APromptRecov4v6_Aug5ReReco_Run2011BPromptRecov1.txt";
//   pars.outputFile = "/data2/yohay/RA3/Data2011A_ToRun167913_Filter-JsonHLTtwo43-30GeVPhosWithR9HoverE_NoPileupCorr_Photon_NEW_categorized_OR.root";
//   pars.outputFile = "/data2/yohay/RA3/1140pb-1_ff_categorized_new.root";
  pars.outputFile = "DATASET_JSON_HLT_PV_single_photon_skim.root";

  TStopwatch ts;

  ts.Start();

  DongwookCategoryProducer producer(pars);

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
