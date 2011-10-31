#include "../../../GMSBTools/Filters/interface/Typedefs.h"

void GMSBAnalysis()
{
  //load source code and precompiled shared libraries
  gSystem->Load("libSusy.so");
  gSystem->Load("../jec/lib/libJetMETObjects.so");
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");
  gSystem->Load("../../../PhysicsTools/lib/libUtilities.so");
  gSystem->Load("SusyEventPrinter_cc.so");
//   gROOT->LoadMacro("SusyEventPrinter.cc+");
  gSystem->SetIncludePath("-I../../..");
  gROOT->LoadMacro("GMSBAnalyzer.C++");

  //instantiate GMSBAnalyzer object
  vector<string> input;
  input.push_back("/data2/yohay/RA3/ntuple_DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  TChain* tree = new TChain("susyTree");
  for (VSTRING_IT iIn = input.begin(); iIn != input.end(); ++iIn) { tree->Add((*iIn).c_str()); }
  GMSBAnalyzer analyzer(tree);
  analyzer.setTag("photons");
  analyzer.setNEvts(1);
  analyzer.setIntLumi(1140.0/*pb^-1*/);
  analyzer.setFileMapEntry("QCD_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   8.16E+08/*pb*/, 11000000);
  analyzer.setFileMapEntry("QCD_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   5.31E+07/*pb*/, 6583068);
  analyzer.setFileMapEntry("QCD_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   6360000.0/*pb*/, 6600000);
  analyzer.setFileMapEntry("QCD_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   784000.0/*pb*/, 6589956);
  analyzer.setFileMapEntry("QCD_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   115000.0/*pb*/, 6127528);
  analyzer.setFileMapEntry("QCD_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   24300.0/*pb*/, 6220160);
  analyzer.setFileMapEntry("QCD_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   1170.0/*pb*/, 6432669);
  analyzer.setFileMapEntry("QCD_Pt-470to600_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   70.2/*pb*/, 3990085);
  analyzer.setFileMapEntry("QCD_Pt-600to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   15.6/*pb*/, 4245695);
  analyzer.setFileMapEntry("QCD_Pt-800to1000_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   1.84/*pb*/, 4053888);
  analyzer.setFileMapEntry("QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   0.332/*pb*/, 2093222);
  analyzer.setFileMapEntry("QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   0.0109/*pb*/, 2196200);
  analyzer.setFileMapEntry("QCD_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   3.58E-04/*pb*/, 293139);
  analyzer.setFileMapEntry("G_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   172000.0/*pb*/, 2046119);
  analyzer.setFileMapEntry("G_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   16700.0/*pb*/, 2187260);
  analyzer.setFileMapEntry("G_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   2720.0/*pb*/, 2036704);
  analyzer.setFileMapEntry("G_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
			   447.0/*pb*/, 2046637);
  analyzer.setFileMapEntry("G_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   84.2/*pb*/, 2088216);
  analyzer.setFileMapEntry("G_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
			   22.6/*pb*/, 2069161);
  analyzer.setFileMapEntry("G_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   1.49/*pb*/, 2076880);
  analyzer.setFileMapEntry("G_Pt-470to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   0.132/*pb*/, 2087212);
  analyzer.setFileMapEntry("G_Pt-800to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
			   0.00348/*pb*/, 2131800);
  analyzer.setFileMapEntry("G_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   1.27E-05/*pb*/, 2198160);
  analyzer.setFileMapEntry("G_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			   2.94E-07/*pb*/, 2188301);
  analyzer.setFileMapEntry("WToENu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			   7899.0/*pb*/, 5334220);
  analyzer.setFileMapEntry("WToMuNu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			   7899.0/*pb*/, 5413258);
  analyzer.setFileMapEntry("WToTauNu_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 
			   7899.0/*pb*/, 5500000);
  analyzer.setFileMapEntry("DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			   1300.0/*pb*/, 2262653);
  analyzer.setFileMapEntry("DYToMuMu_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			   1300.0/*pb*/, 2148325);
  analyzer.setFileMapEntry("DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 1300.0/*pb*/, 2032536);
  analyzer.setFileMapEntry("TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 
			   94.0/*pb*/, 1089625);
  analyzer.setFileMapEntry("DiPhotonBox_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   358.2/*pb*/, 528400);
  analyzer.setFileMapEntry("DiPhotonBox_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   12.37/*pb*/, 518288);
  analyzer.setFileMapEntry("DiPhotonBox_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   2.08E-04/*pb*/, 515028);
  analyzer.setFileMapEntry("DiPhotonBorn_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   236.4/*pb*/, 507554);
  analyzer.setFileMapEntry("DiPhotonBorn_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   22.37/*pb*/, 532864);
  analyzer.setFileMapEntry("DiPhotonBorn_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   0.008072/*pb*/, 526240);
  analyzer.setPUFile("/data2/yohay/RA3/1140pb-1_PU.root");
  analyzer.initPU();

//   const unsigned int size = analyzer.getSize();
//   for (unsigned int i = 0; i < size; ++i) {
//     string process(analyzer.getProcess(i));
//     cout << "Dataset: " << process << endl;
//     cout << "Cross section: " << analyzer.getXSec(process) << " pb\n";
//     cout << "No. events: " << analyzer.getNEvtsProcessed(process) << endl;
//     cout << "Weight: " << analyzer.getWeight(process) << endl;
//   }

  //loop over events
  TStopwatch ts;
  ts.Start();
//   analyzer.runMETAnalysis("/data2/yohay/RA3/debug_noPUReweighting.root");
  analyzer.stripBranch("/data2/yohay/RA3/ntuple_DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root", "susyEvent");
  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

  //clean up
  delete tree;
}
