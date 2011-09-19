void localCountTriggersOUTPUT() {

  // Printing utility for ntuple variables
  gROOT->LoadMacro("SusyEventPrinter.cc+");

  // Main analysis code
  gSystem->SetIncludePath("-I../../..");
  gROOT->LoadMacro("EventAnalyzer.cc+");

  //configuration
  TChain chain("susyTree");
  chain.Add("FILE");
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
  treeReader.countTriggers("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/count_May10ReReco_PromptRecov4OUTPUT");
  ts.Stop();
  cout << "Real time : " << ts.RealTime()/60.0 << " minutes" << endl;
  cout << "CPU time  : " << ts.CpuTime()/60.0 << " minutes" << endl;
}
