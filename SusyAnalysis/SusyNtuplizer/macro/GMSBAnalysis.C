//note: you may have to change the SusyEvent class to run over Summer11 skims from before updates!

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
  gSystem->SetIncludePath("-I../../.. -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms7/include");
  gROOT->LoadMacro("GMSBAnalyzer.C++");

  //instantiate GMSBAnalyzer object
  vector<string> input;
  //   input.push_back("/data2/yohay/RA3/4684pb-1_categorized_18-Jan-12_skim_v2.root");
  //   input.push_back("/data2/yohay/RA3/DiPhotonJets/ntuple_DiPhotonJets_7TeV-madgraph-Summer11-PU_S4_START42_V11-v1_JSON_HLT_PV_skim.root");
  input.push_back("/data2/yohay/RA3/Zee/ntuple_DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v3.root");
  input.push_back("/data2/yohay/RA3/tt/ntuple_v2_TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v3.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/ntuple_DiPhotonBorn_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v3.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/ntuple_DiPhotonBorn_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v3.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/ntuple_DiPhotonBorn_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v3.root");
  //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_DiPhotonBorn_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_DiPhotonBorn_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_DiPhotonBorn_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_DiPhotonBox_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_DiPhotonBox_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_DiPhotonBox_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-470to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-800to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_G_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S4_START42_V11-v1_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_1.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_3.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_4.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_5.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_6.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_7.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_8.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_9.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_10.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_11.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_12.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_13.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2_14.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-470to600_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-600to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-800to1000_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_QCD_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_WToENu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  // //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_WToTauNu_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  //   input.push_back("/data2/yohay/RA3/MCSkims_latestEffectiveAreas/ntuple_v2_TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim_v2.root");
  TChain* tree = new TChain("susyTree");
  for (VSTRING_IT iIn = input.begin(); iIn != input.end(); ++iIn) { tree->Add((*iIn).c_str()); }
  GMSBAnalyzer analyzer(tree);
  analyzer.setTag("photons");
  analyzer.setNEvts(-1);
  analyzer.setIntLumi(4684.0/*pb^-1*/);
  // //   analyzer.setFileMapEntry("QCD_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   8.16E+08/*pb*/, 11000000);
  // //   analyzer.setFileMapEntry("QCD_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   5.31E+07/*pb*/, 6583068);
  // //   analyzer.setFileMapEntry("QCD_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   6360000.0/*pb*/, 6600000);
  // //   analyzer.setFileMapEntry("QCD_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   784000.0/*pb*/, 6589956);
  // //   analyzer.setFileMapEntry("QCD_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   115000.0/*pb*/, 6127528);
  // //   analyzer.setFileMapEntry("QCD_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   24300.0/*pb*/, 6220160);
  // //   analyzer.setFileMapEntry("QCD_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   1170.0/*pb*/, 6432669);
  // //   analyzer.setFileMapEntry("QCD_Pt-470to600_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   70.2/*pb*/, 3990085);
  // //   analyzer.setFileMapEntry("QCD_Pt-600to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   15.6/*pb*/, 4245695);
  // //   analyzer.setFileMapEntry("QCD_Pt-800to1000_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   1.84/*pb*/, 4053888);
  // //   analyzer.setFileMapEntry("QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   0.332/*pb*/, 2093222);
  // //   analyzer.setFileMapEntry("QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   0.0109/*pb*/, 2196200);
  // //   analyzer.setFileMapEntry("QCD_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   3.58E-04/*pb*/, 293139);
  //   analyzer.setFileMapEntry("GJet_Pt-20_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S4_START42_V11-v1", 501.15/*pb*/, 6757937); //xsec*filterEff
  //   analyzer.setFileMapEntry("QCD_Pt-30to40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 10868/*pb*/, 6094670); //xsec*filterEff
  //   analyzer.setFileMapEntry("QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 43571/*pb*/, 40287002); //xsec*filterEff
  // //   analyzer.setFileMapEntry("G_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   172000.0/*pb*/, 2046119);
  // //   analyzer.setFileMapEntry("G_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   16700.0/*pb*/, 2187260);
  // //   analyzer.setFileMapEntry("G_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   2720.0/*pb*/, 2036704);
  // //   analyzer.setFileMapEntry("G_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
  // // 			   447.0/*pb*/, 2046637);
  // //   analyzer.setFileMapEntry("G_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   84.2/*pb*/, 2088216);
  // //   analyzer.setFileMapEntry("G_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
  // // 			   22.6/*pb*/, 2069161);
  // //   analyzer.setFileMapEntry("G_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   1.49/*pb*/, 2076880);
  // //   analyzer.setFileMapEntry("G_Pt-470to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   0.132/*pb*/, 2087212);
  // //   analyzer.setFileMapEntry("G_Pt-800to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
  // // 			   0.00348/*pb*/, 2131800);
  // //   analyzer.setFileMapEntry("G_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   1.27E-05/*pb*/, 2198160);
  // //   analyzer.setFileMapEntry("G_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
  // // 			   2.94E-07/*pb*/, 2188301);
  //   analyzer.setFileMapEntry("WToENu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
  // 			   7899.0/*pb*/, 5334220);
  //   analyzer.setFileMapEntry("WToMuNu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
  // 			   7899.0/*pb*/, 5413258);
  //   analyzer.setFileMapEntry("WToTauNu_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 
  // 			   7899.0/*pb*/, 5500000);
  analyzer.setFileMapEntry("DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			   1300.0/*pb*/, 2262653);
  //   analyzer.setFileMapEntry("DYToMuMu_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
  // 			   1300.0/*pb*/, 2148325);
  //   analyzer.setFileMapEntry("DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 1300.0/*pb*/, 2032536);
  analyzer.setFileMapEntry("TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 
			   94.0/*pb*/, 1089625);
  //   analyzer.setFileMapEntry("DiPhotonBox_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
  // 			   358.2/*pb*/, 528400);
  //   analyzer.setFileMapEntry("DiPhotonBox_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
  // 			   12.37/*pb*/, 518288);
  //   analyzer.setFileMapEntry("DiPhotonBox_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
  // 			   2.08E-04/*pb*/, 515028);
  analyzer.setFileMapEntry("DiPhotonBorn_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   236.4/*pb*/, 507554);
  analyzer.setFileMapEntry("DiPhotonBorn_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   22.37/*pb*/, 532864);
  analyzer.setFileMapEntry("DiPhotonBorn_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			   0.008072/*pb*/, 526240);
  analyzer.setFileMapEntry("DiPhotonJets_7TeV-madgraph-Summer11-PU_S4_START42_V11-v1", 
			   154.7/*pb*/, 1150800);
  analyzer.setPUFile("/data2/yohay/RA3/160404-180252_PU.root");
  analyzer.initPU();
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
  //   analyzer.
  //     runMETAnalysis("/data2/yohay/RA3/4684pb-1_MET_18-Jan-12_skim_noReweighting.root");
  analyzer.runMETAnalysis("/data2/yohay/RA3/4684pb-1_MET_MC_ggPythia.root");
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
