void countTriggers() {

  // Printing utility for ntuple variables
  gROOT->LoadMacro("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/SusyEventPrinter.cc+");
//   gROOT->LoadMacro("SusyEventPrinter.cc+");

  // Main analysis code
//   gSystem->SetIncludePath("-I../../..");
  gSystem->SetIncludePath("-I/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src");
  gROOT->LoadMacro("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/EventAnalyzer.cc+");
//   gROOT->LoadMacro("EventAnalyzer.cc+");

  //configuration
  TChain chain("susyTree");
  chain.Add("susyEvent_ALL_1.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-May10ReReco-v1/Photon/Runs160442-163869/susyEvent_ALL_1.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-May10ReReco-v1/Photon/Runs160442-163869/susyEvent_ALL_2.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Run166438/susyEvent_ALL.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs165088-166346/susyEvent_ALL_1.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs165088-166346/susyEvent_ALL_2.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166374-166486/susyEvent_ALL.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166502-166530/susyEvent_ALL.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166554-166787/susyEvent_ALL.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166839-166911/susyEvent_ALL.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166921-167078/susyEvent_ALL.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs167098-167284/susyEvent_ALL.root");
//   chain.Add("/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs167551-167913/susyEvent_ALL.root");
  EventAnalyzer treeReader(&chain);
  treeReader.SetPrintInterval(10000);
  treeReader.SetPrintLevel(0);
  treeReader.SetUseTrigger(true);
  treeReader.AddHltName("HLT_Photon32_CaloIdL_Photon26_CaloIdL");
  treeReader.AddHltName("HLT_Photon36_CaloIdL_Photon22_CaloIdL");
  treeReader.AddHltName("HLT_Photon40_CaloIdL_Photon28_CaloIdL");
  treeReader.SetFilter(false);
  treeReader.SetProcessNEvents(-1);
  treeReader.IncludeAJson("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/Cert_160404-172255_7TeV_PromptReco_Collisions11_JSON.txt");
  treeReader.SetPhotonTag("photons");

  //run
  TStopwatch ts;
  ts.Start();
//   treeReader.countTriggers("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/count_May10ReReco_PromptRecov4_testv3");
  treeReader.countTriggers("count_May10ReReco_PromptRecov4_batchTest");
  ts.Stop();
  cout << "Real time : " << ts.RealTime()/60.0 << " minutes" << endl;
  cout << "CPU time  : " << ts.CpuTime()/60.0 << " minutes" << endl;
}
