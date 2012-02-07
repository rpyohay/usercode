// Original Author:  Dongwook Jang
// $Id: ana.C,v 1.8 2011/11/01 22:14:51 dwjang Exp $
//
// Jet energy correction is possible at ntuple level.
// $ cd ../jec/JetMETObjects
// $ make
// This will create a shared library in jec/lib
// which is included below as libJetMETObjects.so
//
// Come back to this directory and do
// $ make
// $ root -b -q -l ana.C
// will produce hist_"physics"_"ds".root

void ana(TString ds="relval", TString physics="ttbar") {

  //precompiled and to-be-compiled libraries
  gSystem->Load("libSusy.so");
  gSystem->Load("../jec/lib/libJetMETObjects.so");
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");
//   gROOT->LoadMacro("SusyEventPrinter.cc++");
  gSystem->Load("SusyEventPrinter_cc.so");
  gSystem->SetIncludePath("-I../../..");
  gROOT->LoadMacro("SusyEventAnalyzer.cc++");

  //ntuples to run over
  vector<string> input;
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_10_1_flU.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_1_1_Dz4.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_2_2_hJ4.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_3_1_RDQ.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_4_2_j8E.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_5_2_Fj5.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_6_1_6kI.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_7_2_JW1.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_8_2_Aam.root");
  input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_9_1_ZtH.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_10_1_Kee.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_1_1_jAk.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_2_1_jxR.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_3_1_pAZ.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_4_1_H56.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_5_1_8XS.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_6_1_bZy.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_7_1_hti.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_8_2_Hyr.root");
  input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_9_2_EHt.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_10_1_itX.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_11_1_kvU.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_1_1_J6F.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_2_1_iun.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_3_1_SPV.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_4_1_p43.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_5_1_KLb.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_6_1_JBN.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_7_1_D5Z.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_8_2_rs7.root");
  input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_9_1_Vvk.root");
  TChain* tree = new TChain("susyTree");
  for (vector<string>::const_iterator iIn = input.begin(); iIn != input.end(); 
       ++iIn) { tree->Add((*iIn).c_str()); }

  //analyzer object
  SusyEventAnalyzer analyzer(tree);
  analyzer.SetProcessNEvents(-1);
  analyzer.setIntLumi(4684.0/*pb^-1*/);
  analyzer.setFileMapEntry("diphoton_10-25", 236.4/*pb*/, 507554);
  analyzer.setFileMapEntry("diphoton_25-250", 22.37/*pb*/, 532864);
  analyzer.setFileMapEntry("diphoton_250", 0.008072/*pb*/, 526240);

  //run
  TStopwatch ts;
  ts.Start();
  analyzer.plot("/data2/yohay/RA3/pThat_Summer11PythiaDiPhotonBorn_noCuts.root");
  ts.Stop();
  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

  //clean up
  delete tree;
}
