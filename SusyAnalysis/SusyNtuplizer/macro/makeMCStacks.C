#include "../../../GMSBTools/Filters/interface/Typedefs.h"

void makeMCStacks()
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
  gSystem->SetIncludePath("-I../../.. -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms7/include");
  gSystem->Load("GMSBAnalyzer_C.so");

  //set up THStacks for all the MC
  THStack mee("mee", "mee");
  THStack nJets("nJets", "nJets");

  //set up vector of pointers to TFile objects
  vector<TFile*> pFileVec;

  //create vector of vectors of MC processes
  vector<vector<string> > MC;
  MC.push_back(vector<string>(1, "/data2/yohay/RA3/ntuple_DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root"));
  MC[1] = vector<string>(1, "/data2/yohay/RA3/ntuple_DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[2] = vector<string>();
  MC[2].push_back("/data2/yohay/RA3/ntuple_DiPhotonBorn_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[2].push_back("/data2/yohay/RA3/ntuple_DiPhotonBorn_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[2].push_back("/data2/yohay/RA3/ntuple_DiPhotonBorn_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[2].push_back("/data2/yohay/RA3/ntuple_DiPhotonBox_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[2].push_back("/data2/yohay/RA3/ntuple_DiPhotonBox_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[2].push_back("/data2/yohay/RA3/ntuple_DiPhotonBox_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[3] = vector<string>();
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-470to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-800to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[3].push_back("/data2/yohay/RA3/ntuple_G_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4] = vector<string>();
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-470to600_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-600to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-800to1000_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[4].push_back("/data2/yohay/RA3/ntuple_QCD_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[5] = vector<string>(1, "/data2/yohay/RA3/ntuple_WToENu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[6] = vector<string>(1, "/data2/yohay/RA3/ntuple_WToTauNu_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");
  MC[7] = vector<string>(1, "/data2/yohay/RA3/ntuple_v2_TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2_JSON_HLT_PV_skim.root");

  //create pointer to GMSBAnalyzer object
  GMSBAnalyzer* pAnalyzer = NULL;

  //loop over vector of MC processes
  for (vector<vector<string> >::const_iterator iMC = MC.begin(); iMC != MC.end(); ++iMC) {

    //set up the TChain
    TChain* tree = new TChain("susyTree");
    for (VSTRING_IT iIn = iMC->begin(); iIn != iMC->end(); ++iIn) {
      tree->Add((*iIn).c_str());
    }

    //instantiate the GMSBAnalyzer object
    pAnalyzer = new GMSBAnalyzer(tree);
    pAnalyzer->setTag("photons");
    pAnalyzer->setNEvts(-1);
    pAnalyzer->setIntLumi(4684.0/*pb^-1*/);
    pAnalyzer->setFileMapEntry("QCD_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       8.16E+08/*pb*/, 11000000);
    pAnalyzer->setFileMapEntry("QCD_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       5.31E+07/*pb*/, 6583068);
    pAnalyzer->setFileMapEntry("QCD_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       6360000.0/*pb*/, 6600000);
    pAnalyzer->setFileMapEntry("QCD_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       784000.0/*pb*/, 6589956);
    pAnalyzer->
      setFileMapEntry("QCD_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      115000.0/*pb*/, 6127528);
    pAnalyzer->
      setFileMapEntry("QCD_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      24300.0/*pb*/, 6220160);
    pAnalyzer->
      setFileMapEntry("QCD_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      1170.0/*pb*/, 6432669);
    pAnalyzer->
      setFileMapEntry("QCD_Pt-470to600_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      70.2/*pb*/, 3990085);
    pAnalyzer->
      setFileMapEntry("QCD_Pt-600to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      15.6/*pb*/, 4245695);
    pAnalyzer->
      setFileMapEntry("QCD_Pt-800to1000_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      1.84/*pb*/, 4053888);
    pAnalyzer->
      setFileMapEntry("QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      0.332/*pb*/, 2093222);
    pAnalyzer->
      setFileMapEntry("QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      0.0109/*pb*/, 2196200);
    pAnalyzer->setFileMapEntry("QCD_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       3.58E-04/*pb*/, 293139);
    pAnalyzer->setFileMapEntry("G_Pt-15to30_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       172000.0/*pb*/, 2046119);
    pAnalyzer->setFileMapEntry("G_Pt-30to50_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       16700.0/*pb*/, 2187260);
    pAnalyzer->setFileMapEntry("G_Pt-50to80_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       2720.0/*pb*/, 2036704);
    pAnalyzer->setFileMapEntry("G_Pt-80to120_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
			       447.0/*pb*/, 2046637);
    pAnalyzer->setFileMapEntry("G_Pt-120to170_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       84.2/*pb*/, 2088216);
    pAnalyzer->setFileMapEntry("G_Pt-170to300_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
			       22.6/*pb*/, 2069161);
    pAnalyzer->setFileMapEntry("G_Pt-300to470_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       1.49/*pb*/, 2076880);
    pAnalyzer->setFileMapEntry("G_Pt-470to800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       0.132/*pb*/, 2087212);
    pAnalyzer->setFileMapEntry("G_Pt-800to1400_TuneZ2_7TeV_pythia6-Summer11-PU_S4_START42_V11-v1", 
			       0.00348/*pb*/, 2131800);
    pAnalyzer->
      setFileMapEntry("G_Pt-1400to1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
		      1.27E-05/*pb*/, 2198160);
    pAnalyzer->setFileMapEntry("G_Pt-1800_TuneZ2_7TeV_pythia6-Summer11-PU_S3_START42_V11-v2", 
			       2.94E-07/*pb*/, 2188301);
    pAnalyzer->setFileMapEntry("WToENu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			       7899.0/*pb*/, 5334220);
    pAnalyzer->setFileMapEntry("WToMuNu_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			       7899.0/*pb*/, 5413258);
    pAnalyzer->
      setFileMapEntry("WToTauNu_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 
		      7899.0/*pb*/, 5500000);
    pAnalyzer->setFileMapEntry("DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			       1300.0/*pb*/, 2262653);
    pAnalyzer->setFileMapEntry("DYToMuMu_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2", 
			       1300.0/*pb*/, 2148325);
    pAnalyzer->
      setFileMapEntry("DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 
		      1300.0/*pb*/, 2032536);
    pAnalyzer->setFileMapEntry("TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2", 
			       94.0/*pb*/, 1089625);
    pAnalyzer->setFileMapEntry("DiPhotonBox_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			       358.2/*pb*/, 528400);
    pAnalyzer->
      setFileMapEntry("DiPhotonBox_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
		      12.37/*pb*/, 518288);
    pAnalyzer->setFileMapEntry("DiPhotonBox_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			       2.08E-04/*pb*/, 515028);
    pAnalyzer->
      setFileMapEntry("DiPhotonBorn_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
		      236.4/*pb*/, 507554);
    pAnalyzer->
      setFileMapEntry("DiPhotonBorn_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
		      22.37/*pb*/, 532864);
    pAnalyzer->setFileMapEntry("DiPhotonBorn_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2", 
			       0.008072/*pb*/, 526240);
    pAnalyzer->setPUFile("/data2/yohay/RA3/160404-180252_PU.root");
    pAnalyzer->initPU();
    pAnalyzer->setL1JECFile("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L1FastJet.txt");
    pAnalyzer->setL2JECFile("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L2Relative.txt");
    pAnalyzer->setL3JECFile("/afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L3Absolute.txt");

    //run analyzer
    const string beginTag = "ntuple_";
    const string endTag = "_7TeV";
    const string name(*(iMC->begin()));
    string extractedName;
    const size_t beginPos = name.find(beginTag);
    const size_t endPos = name.find(endTag);
    if ((beginPos != string::npos) && (endPos != string::npos) && (beginPos <= endPos)) {
      extractedName = 
	name.substr(beginPos + beginTag.length(), endPos - beginPos - beginTag.length());
    }
    else {
      cerr << "Error extracting MC sample name.\n";
      cerr << "     name = " << name << endl;
      cerr << "     beginPos = " << beginPos << endl;
      cerr << "     endPos = " << endPos << endl;
      cerr << "Using name \"dummy\".\n";
      extractedName = "dummy";
    }
    string fileName = "/data2/yohay/RA3/4684pb-1_" + extractedName + ".root";
    
    TStopwatch ts;
    ts.Start();
    pAnalyzer->compareDataToMC(fileName);
    ts.Stop();
    std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
    std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

    //clean up
    delete tree;
    delete pAnalyzer;
    tree = NULL;
    pAnalyzer = NULL;

    //add histograms from this MC dataset to the stack
    pFileVec.push_back(new TFile(fileName.c_str()));
    if ((pFileVec.end() - 1)->IsOpen()) {
      TH1F* pMEE = NULL;
      TH1F* pNJets = NULL;
      (pFileVec.end() - 1)->GetObject("mee", pMEE);
      (pFileVec.end() - 1)->GetObject("nJets", pNJets);
      if ((pMEE != NULL) && (pNJets != NULL)) {
	pMEE->SetFillColor(kRed - iMC - MC.begin());
	pNJets->SetFillColor(kRed - iMC - MC.begin());
	mee.Add(pMEE);
	nJets.Add(pNJets);
      }
      else {
	cerr << "Error: histogram not found.\n";
	cerr << "     pMEE = " << pMEE << endl;
	cerr << "     pNJets = " << pNJets << endl;
	cerr << "     fileName = " << fileName << endl;
	cerr << "Skipping this MC dataset.\n";
	pFileVec.erase((pFileVec.end() - 1), pFileVec.end());
      }
    }
    else cerr << "Error: could not open file " << fileName << ".  Skipping this MC dataset.\n";
  }

  //save stacked histograms
  string outName = "/data2/yohay/RA3/MC_stacks.root";
  TFile out(outName.c_str(), "RECREATE");
  if (out.IsOpen()) {
    out.cd();
    TCanvas meeCanvas("meeCanvas", "meeCanvas", 600, 600);
    meeCanvas.cd();
    mee.Draw();
    TCanvas nJetsCanvas("nJetsCanvas", "nJetsCanvas", 600, 600);
    nJetsCanvas.cd();
    nJets.Draw();
    out.Write();
    out.Close();
  }
  else cerr << "Error: could not open file " << outName << ".\n";

  //deallocate memory for open files
  for (vector<TFile*>::const_iterator i = pFileVec.begin(); i != pFileVec.end(); ++i) {
    delete *i;
    i = NULL;
  }
}
