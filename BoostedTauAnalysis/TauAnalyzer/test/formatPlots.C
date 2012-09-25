{
  //load
  gROOT->Reset();
  gROOT->ProcessLine("#include <utility>");
  string macroPath("/afs/cern.ch/user/y/yohay/CMSSW_5_2_4_patch3_20Sep12/src/BoostedTauAnalysis/");
  macroPath+="TauAnalyzer/test/";
  gSystem->Load((macroPath + "STLDictionary.so").c_str());
  gROOT->LoadMacro((macroPath + "Plot.C++").c_str());
//   gSystem->Load((macroPath + "Plot_C.so").c_str());

  //unit strings
  string unitPTTau("Visible p_{T}^{#tau} (GeV)");
  string unitPTMu("Gen p_{T}^{#mu} (GeV)");
  string unitEtaTau("Visible #eta^{#tau}");
  string unitEtaMu("Gen #eta^{#mu}");
  string unitDR("#DeltaR(visible gen #tau, gen #mu)");
  string noUnit("");

  //map of bin labels for certain efficiency plots
  vector<string> binLabels;
  binLabels.push_back("kNull");
  binLabels.push_back("kOneProng0PiZero");
  binLabels.push_back("kOneProng1PiZero");
  binLabels.push_back("kOneProng2PiZero");
  binLabels.push_back("kOneProng3PiZero");
  binLabels.push_back("kOneProngNPiZero");
  binLabels.push_back("kTwoProng0PiZero");
  binLabels.push_back("kTwoProng1PiZero");
  binLabels.push_back("kTwoProng2PiZero");
  binLabels.push_back("kTwoProng3PiZero");
  binLabels.push_back("kTwoProngNPiZero");
  binLabels.push_back("kThreeProng0PiZero");
  binLabels.push_back("kThreeProng1PiZero");
  binLabels.push_back("kThreeProng2PiZero");
  binLabels.push_back("kThreeProng3PiZero");
  binLabels.push_back("kThreeProngNPiZero");
  binLabels.push_back("kRareDecayMode");
  map<string, vector<string> > binLabelMap;
  binLabelMap["muHadGenDecayMode"] = binLabels;
  binLabelMap["muHadCorrectRecoDecayModeGenDecayMode"] = binLabels;
  binLabelMap["muHadRecoDecayMode"] = binLabels;
  binLabelMap["muHadGen1ProngRecoDecayMode"] = binLabels;
  binLabelMap["muHadGen1Prong1Pi0RecoDecayMode"] = binLabels;
  binLabelMap["muHadGen3ProngRecoDecayMode"] = binLabels;

  //map of inputs to efficiency histograms
  map<string, pair<string, string> > effHistMap;
  effHistMap["numeratorPT"] = make_pair(string("denominatorPT"), unitPTMu);
  map<string, pair<string, string> > effHistMapMu;
  effHistMapMu["numeratorPT"] = make_pair(string("denominatorPT"), unitPTMu);
  effHistMapMu["numeratorEta"] = make_pair(string("denominatorEta"), unitEtaMu);
  map<string, pair<string, string> > effHistMapTau;
  effHistMapTau["numeratorPT"] = make_pair(string("denominatorPT"), unitPTTau);
  effHistMapTau["numeratorEta"] = make_pair(string("denominatorEta"), unitEtaTau);

  //map of inputs to 1D histograms
  map<string, string> hist1DMap;
  hist1DMap["numeratorPT"] = unitPTMu;
  hist1DMap["denominatorPT"] = unitPTMu;
  map<string, string> hist1DMapMu;
  hist1DMapMu["numeratorPT"] = unitPTMu;
  hist1DMapMu["denominatorPT"] = unitPTMu;
  hist1DMapMu["numeratorEta"] = unitEtaMu;
  hist1DMapMu["denominatorEta"] = unitEtaMu;
  map<string, string> hist1DMapTau;
  hist1DMapTau["numeratorPT"] = unitPTTau;
  hist1DMapTau["denominatorPT"] = unitPTTau;
  hist1DMapTau["numeratorEta"] = unitEtaTau;
  hist1DMapTau["denominatorEta"] = unitEtaTau;

//   vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > histMap;
//   histMap.push_back(make_pair(make_pair(), make_pair()));

//   map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, pair<Color_t, Style_t> > > > 
//     canvasMap;
//   canvasMap[make_pair(string("muHadGen1ProngRecoDecayMode"), noUnit)] = ;
//   canvasMap[make_pair(string("muHadGen1Prong1Pi0RecoDecayMode"), noUnit)] = ;
//   canvasMap[make_pair(string("muHadGen3ProngRecoDecayMode"), noUnit)] = ;

  //plot
//   string path("/data1/yohay/efficiency_files/combined/NMSSMHiggs_gg_");
//   string path("/afs/cern.ch/user/y/yohay/CMSSW_5_2_4_patch3_clone/src/NMSSMHiggs_gg_");
//   string path("/data1/yohay/DYMuMu_nonIso_");
//   string path("/data1/yohay/NMSSMHiggs_gg_");
//   string suffixIn("_combined.root");
//   string suffixIn("_JOBNUM.root");
//   string suffixOut("_eff.root");
//   plotNice(path + "muLeg_tauMuX_Mu40eta2p1" + suffixIn, effHistMapMu, binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_tauMuX_Mu40eta2p1" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_tauMuX_IsoMu24eta2p1" + suffixIn, effHistMapMu, binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_tauMuX_IsoMu24eta2p1" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_tauMuX_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixIn, 
// 	   effHistMapMu, binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_tauMuX_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "muLeg_tauMuX_IsoMu18eta2p1LooseIsoPFTau20" + suffixIn, effHistMapMu, 
// 	   binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_tauMuX_IsoMu18eta2p1LooseIsoPFTau20" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_tauMuX_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixIn, 
// 	   effHistMapMu, binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_tauMuX_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "muLeg_tauMuX_IsoMu15eta2p1L1ETM20" + suffixIn, effHistMapMu, binLabelMap, 
// 	   hist1DMapMu, path + "muLeg_tauMuX_IsoMu15eta2p1L1ETM20" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_diTauMu_Mu40eta2p1" + suffixIn, effHistMapMu, binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_diTauMu_Mu40eta2p1" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_diTauMu_IsoMu24eta2p1" + suffixIn, effHistMapMu, binLabelMap, 
// 	   hist1DMapMu, path + "muLeg_diTauMu_IsoMu24eta2p1" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_diTauMu_IsoMu18eta2p1LooseIsoPFTau20" + suffixIn, effHistMapMu, 
// 	   binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_diTauMu_IsoMu18eta2p1LooseIsoPFTau20" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_diTauMu_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixIn, 
// 	   effHistMapMu, binLabelMap, hist1DMapMu, path + 
// 	   "muLeg_diTauMu_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "muLeg_diTauMu_IsoMu15eta2p1L1ETM20" + suffixIn, effHistMapMu, binLabelMap, 
// 	   hist1DMapMu, path + "muLeg_diTauMu_IsoMu15eta2p1L1ETM20" + suffixOut, "noPDF");
//   plotNice(path + "tauLeg_tau1ProngX_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixIn, 
// 	   effHistMapTau, binLabelMap, hist1DMapTau, 
// 	   path + "tauLeg_tau1ProngX_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "tauLeg_tau1ProngX_IsoMu18eta2p1LooseIsoPFTau20" + suffixIn, effHistMapTau, 
// 	   binLabelMap, hist1DMapTau, 
// 	   path + "tauLeg_tau1ProngX_IsoMu18eta2p1LooseIsoPFTau20" + suffixOut, "noPDF");
//   plotNice(path + 
// 	   "tauLeg_tau1ProngX_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixIn, 
// 	   effHistMapTau, binLabelMap, hist1DMapTau, path + 
// 	   "tauLeg_tau1ProngX_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "tauLeg_tau1ProngX_LooseIsoPFTau35Trk20Prong1" + suffixIn, effHistMapTau, 
// 	   binLabelMap, hist1DMapTau, 
// 	   path + "tauLeg_tau1ProngX_LooseIsoPFTau35Trk20Prong1" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_tauMutau1Prong_Mu40eta2p1" + suffixIn, effHistMapMu, binLabelMap, 
// 	   hist1DMapMu, path + "muLeg_tauMutau1Prong_Mu40eta2p1" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_tauMutau1Prong_IsoMu24eta2p1" + suffixIn, effHistMapMu, binLabelMap, 
// 	   hist1DMapMu, path + "muLeg_tauMutau1Prong_IsoMu24eta2p1" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_tauMutau1Prong_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixIn, 
// 	   effHistMapMu, binLabelMap, hist1DMapMu, path + 
// 	   "muLeg_tauMutau1Prong_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "muLeg_tauMutau1Prong_IsoMu18eta2p1LooseIsoPFTau20" + suffixIn, effHistMapMu, 
// 	   binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_tauMutau1Prong_IsoMu18eta2p1LooseIsoPFTau20" + suffixOut, "noPDF");
//   plotNice(path + 
// 	   "muLeg_tauMutau1Prong_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixIn, 
// 	   effHistMapMu, binLabelMap, hist1DMapMu, path + 
// 	   "muLeg_tauMutau1Prong_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "muLeg_tauMutau1Prong_IsoMu15eta2p1L1ETM20" + suffixIn, effHistMapMu, 
// 	   binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_tauMutau1Prong_IsoMu15eta2p1L1ETM20" + suffixOut, "noPDF");
//   plotNice(path + "muLeg_diTauMu_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixIn, 
// 	   effHistMapMu, binLabelMap, hist1DMapMu, 
// 	   path + "muLeg_diTauMu_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "tauLeg_tauXX_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixIn, 
// 	   effHistMapTau, binLabelMap, hist1DMapTau, 
// 	   path + "tauLeg_tauXX_IsoMu18eta2p1MediumIsoPFTau25Trk5eta2p1" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "tauLeg_tauXX_IsoMu18eta2p1LooseIsoPFTau20" + suffixIn, effHistMapTau, 
// 	   binLabelMap, hist1DMapTau, 
// 	   path + "tauLeg_tauXX_IsoMu18eta2p1LooseIsoPFTau20" + suffixOut, "noPDF");
//   plotNice(path + "tauLeg_tauXX_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixIn, 
// 	   effHistMapTau, binLabelMap, hist1DMapTau, 
// 	   path + "tauLeg_tauXX_IsoMu15eta2p1LooseIsoPFTau35Trk20Prong1L1ETM20" + suffixOut, 
// 	   "noPDF");
//   plotNice(path + "tauLeg_tauXX_LooseIsoPFTau35Trk20Prong1" + suffixIn, effHistMapTau, 
// 	   binLabelMap, hist1DMapTau, 
// 	   path + "tauLeg_tauXX_LooseIsoPFTau35Trk20Prong1" + suffixOut, "noPDF");
  //   string beginSignal("/data1/yohay/NMSSMHiggs_gg_trigger_");
  //   string endSignal(".root");
//   string beginDY("/data1/yohay/Summer12_DYToTauTau_");
//   string endDY("_genMuPTGt20GeV_genTauPTGt15GeV.root");
//   plotNice(beginSignal + "analysis" + endSignal, effHistMap, binLabelMap, hist1DMap, 
// 	   beginSignal + "plots" + endSignal, "noPDF");
//   plotNice(beginDY + "analysis" + endDY, effHistMap, binLabelMap, hist1DMap, 
// 	   beginDY + "plots" + endDY, "noPDF");

// plotNiceDifferentFiles(map<pair<string, string>, vector<pair<pair<TFile*, Option_t*>, 
// 			    pair<Color_t, Style_t> > > >& canvasMap, 
// 			    const map<string, vector<string> >& binLabelMap, 
// 			    const string& outputFileName, const string& savePath)


//space-saving constant definitions
  const string analysisFilePath("/data1/yohay/");
  const string objTag("_");
  const string suffix(".root");
  const string effTag("_eff");
  const string controlSample("NMSSMHiggs_WH");
  const string signalSample("NMSSMHiggs_WH");
  const string leg("");
  const string trigger("LooseCombinedIsolationDBSumPtCorr");
  const string signalGenSel("muHad_MCTruthMuonRemoval_");
  const string controlGenSel("muHad_");
  const string outputTag("MCTruthMuonRemoval_" + trigger);

  //make individual efficiency plots for signal and Z-->mumu samples
  vector<string> effInputFiles;
  vector<string> comparisonInputFiles;
  effInputFiles.push_back(analysisFilePath + controlSample + objTag + controlGenSel + trigger);
  effInputFiles.push_back(analysisFilePath + signalSample + objTag + leg + signalGenSel + trigger);
  for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
       ++iFile) {
    comparisonInputFiles.push_back(*iFile + effTag + suffix);
    plotNice(*iFile + suffix, effHistMapTau, binLabelMap, hist1DMapTau, *iFile + effTag + suffix, 
	     "noPDF");
  }

  //compare Z-->mumu and signal muon efficiency
  vector<string> canvasNames;
  canvasNames.push_back("eff_numeratorPT_over_denominatorPT");
  canvasNames.push_back("eff_numeratorEta_over_denominatorEta");
  vector<string> graphNames;
  graphNames.push_back("divide_numeratorPT_by_denominatorPT");
  graphNames.push_back("divide_numeratorEta_by_denominatorEta");
  Color_t colors[2] = {kBlack, kRed};
  Style_t styles[2] = {20, 21};
  drawMultipleEfficiencyGraphsOn1Canvas(analysisFilePath + "effVsPT" + objTag + outputTag + 
					suffix, comparisonInputFiles, canvasNames, graphNames, 
					colors, styles);
}
