{
  //load
  gROOT->Reset();
  gROOT->ProcessLine("#include <utility>");
  string macroPath("/afs/cern.ch/user/y/yohay/CMSSW_5_3_2_patch4/src/BoostedTauAnalysis/");
  macroPath+="TauAnalyzer/test/";
  gSystem->Load((macroPath + "STLDictionary.so").c_str());
  gROOT->LoadMacro((macroPath + "Plot.C++").c_str());
//   gSystem->Load((macroPath + "Plot_C.so").c_str());

  //needed so vector<Color_t> and vector<Style_t> work
  vector<short> dummy;

  //unit strings
  string unitPTTau("Reco #tau p_{T} (GeV)");
  string unitPTMu("Reco #mu p_{T} (GeV)");
  string unitEtaTau("Reco #tau #eta");
  string unitEtaMu("Reco #mu #eta");
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

  //set up canvas and graph names
  vector<string> canvasNames;
  canvasNames.push_back("hadTauAssociatedMuMultiplicityCanvas");
  canvasNames.push_back("muHadMassCanvas");
  canvasNames.push_back("muHadChargeCanvas");
  canvasNames.push_back("METCanvas");
  canvasNames.push_back("WMuMTCanvas");
  canvasNames.push_back("tauMuMTCanvas");
  canvasNames.push_back("dPhiWMuMETCanvas");
  canvasNames.push_back("dPhiTauMuMETCanvas");
  canvasNames.push_back("tauMuTauHadJetHTCanvas");
  canvasNames.push_back("diJetHTCanvas");
  canvasNames.push_back("jetTauJetHTCanvas");
  canvasNames.push_back("tauMuTauHadJetWMuHTCanvas");
  canvasNames.push_back("diJetWMuHTCanvas");
  canvasNames.push_back("jetTauJetWMuHTCanvas");
  canvasNames.push_back("tauMuPTCanvas");
  canvasNames.push_back("tauHadPTCanvas");
  canvasNames.push_back("tauHadIsoCanvas");
  canvasNames.push_back("softMuIsoCandMultiplicityCanvas");
  vector<string> graphNames;
  graphNames.push_back("hadTauAssociatedMuMultiplicity");
  graphNames.push_back("muHadMass");
  graphNames.push_back("muHadCharge");
  graphNames.push_back("MET");
  graphNames.push_back("WMuMT");
  graphNames.push_back("tauMuMT");
  graphNames.push_back("dPhiWMuMET");
  graphNames.push_back("dPhiTauMuMET");
  graphNames.push_back("tauMuTauHadJetHT");
  graphNames.push_back("diJetHT");
  graphNames.push_back("jetTauJetHT");
  graphNames.push_back("tauMuTauHadJetWMuHT");
  graphNames.push_back("diJetWMuHT");
  graphNames.push_back("jetTauJetWMuHT");
  graphNames.push_back("tauMuPT");
  graphNames.push_back("tauHadPT");
  graphNames.push_back("tauHadIso");
  graphNames.push_back("softMuIsoCandMultiplicity");

  //set up plot style options
  vector<string> legendHeaders20InvFb(canvasNames.size(), "Normalized to 20 fb^{-1}");
  vector<string> legendHeaders1(canvasNames.size(), "Normalized to 1");
  vector<Color_t> colors;
  colors.push_back(kBlack);
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kGreen);
  vector<Style_t> styles;
  styles.push_back(20);
  styles.push_back(21);
  styles.push_back(22);
  styles.push_back(23);
  styles.push_back(24);
  vector<string> legendEntriesWNJetsToLNuInd;
  legendEntriesWNJetsToLNuInd.push_back("Wh_{1}");
  legendEntriesWNJetsToLNuInd.push_back("W + 1 jet");
  legendEntriesWNJetsToLNuInd.push_back("W + 2 jets");
  legendEntriesWNJetsToLNuInd.push_back("W + 3 jets");
  legendEntriesWNJetsToLNuInd.push_back("W + 4 jets");
  vector<string> legendEntriesWNJetsToLNu;
  legendEntriesWNJetsToLNu.push_back("Wh_{1}");
  legendEntriesWNJetsToLNu.push_back("W + jets");
  vector<string> legendEntriesSearchVsControl;
  legendEntriesSearchVsControl.push_back("Isolated #tau leptons");
  legendEntriesSearchVsControl.push_back("Non-isolated #tau leptons");
  const bool setLogY = true;
  const bool setLinY = false;
  const bool drawStack = true;
  const bool drawSame = false;

  //weights (sig. figs are probably wrong)
  const float Wh1Weight = 0.0681;
  vector<float> weightsWJetsToLNu;
  weightsWJetsToLNu.push_back(Wh1Weight);
  weightsWJetsToLNu.push_back(10.4);
  vector<float> weights1(5, 0.0);
  vector<float> weightsWNJetsToLNuInd;
  weightsWNJetsToLNuInd.push_back(Wh1Weight);
  weightsWNJetsToLNuInd.push_back(4.57);
  weightsWNJetsToLNuInd.push_back(1.01);
  weightsWNJetsToLNuInd.push_back(0.654);
  weightsWNJetsToLNuInd.push_back(0.313);
  vector<float> weightsWNJetsToLNu;
  weightsWNJetsToLNu.push_back(Wh1Weight);
  weightsWNJetsToLNu.push_back(6524.691358);
  vector<float> WNJetsToLNuRelXSecWeights;
  WNJetsToLNuRelXSecWeights.push_back(0.00070041554);
  WNJetsToLNuRelXSecWeights.push_back(0.00015423878);
  WNJetsToLNuRelXSecWeights.push_back(0.00010017824);
  WNJetsToLNuRelXSecWeights.push_back(0.000047972683);

  //space-saving constant definitions
  const string analysisFilePath("/data1/yohay/");
  const string fileExt(".root");
  const string tag20InvFb("_20fb-1");
  const string tag1("_normalizedTo1");
  const string vTag("_v3");

  //hadd
  string isoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadIsoAnalysis_W");
  string nonIsoPrefix(analysisFilePath + "WNJetsToLNu/analysis/muHadNonIsoAnalysis_W");
  string suffix("JetsToLNu_v3" + fileExt);
  string isoHaddOutputFile(isoPrefix + "N" + suffix);
  string nonIsoHaddOutputFile(nonIsoPrefix + "N" + suffix);
  vector<string> isoHaddInputFiles;
  vector<string> nonIsoHaddInputFiles;
  for (unsigned int iNJets = 1; iNJets <= 4; ++iNJets) {
    stringstream isoName;
    isoName << isoPrefix << iNJets << suffix;
    isoHaddInputFiles.push_back(isoName.str());
    stringstream nonIsoName;
    nonIsoName << nonIsoPrefix << iNJets << suffix;
    nonIsoHaddInputFiles.push_back(nonIsoName.str());
  }
  haddCanvases(isoHaddOutputFile, isoHaddInputFiles, WNJetsToLNuRelXSecWeights, canvasNames, 
	       graphNames);
  haddCanvases(nonIsoHaddOutputFile, nonIsoHaddInputFiles, WNJetsToLNuRelXSecWeights, canvasNames, 
	       graphNames);

  //compare signal to background
  string sigVsBkgOutputFile20InvFb(analysisFilePath + "Wh1VsWNJets_muHadIsoAnalysis" + 
				   tag20InvFb + vTag + fileExt);
  string sigVsBkgOutputFile1(analysisFilePath + "Wh1VsWNJets_muHadIsoAnalysis" + tag1 + vTag + 
			     fileExt);
  vector<string> sigVsBkgInputFiles;
  sigVsBkgInputFiles.push_back(isoHaddOutputFile);
  sigVsBkgInputFiles.push_back(analysisFilePath + 
			       "Wh1_Medium/muHadAnalysisV20_passMediumIso_tauHadPTGt0" + fileExt);
//   colors[1] = kBlack;
//   colors[0] = kRed;
//   styles[1] = 20;
//   styles[0] = 21;
//   legendEntriesWNJetsToLNu[1] = "Wh_{1}";
//   legendEntriesWNJetsToLNu[0] = "W + jets";
//   weightsWNJetsToLNu[1] = Wh1Weight;
//   weightsWNJetsToLNu[0] = 6524.691358;
//   drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile20InvFb, sigVsBkgInputFiles, 
// 					canvasNames, graphNames, legendHeaders20InvFb, colors, 
// 					styles, legendEntriesWNJetsToLNu, weightsWNJetsToLNu, 
// 					setLogY, drawStack);
//   drawMultipleEfficiencyGraphsOn1Canvas(sigVsBkgOutputFile1, sigVsBkgInputFiles, 
// 					canvasNames, graphNames, legendHeaders1, colors, 
// 					styles, legendEntriesWNJetsToLNu, weights1, setLinY, 
// 					drawSame);

  //compare search sample to control sample
  string searchVsControlOutputFile(analysisFilePath + "WNJetsToLNu/analysis/isoVsNonIsoTaus" + 
				   tag1 + vTag + fileExt);
  vector<string> searchVsControlInputFiles;
  searchVsControlInputFiles.push_back(isoHaddOutputFile);
  searchVsControlInputFiles.push_back(nonIsoHaddOutputFile);
  drawMultipleEfficiencyGraphsOn1Canvas(searchVsControlOutputFile, searchVsControlInputFiles, 
					canvasNames, graphNames, legendHeaders1, colors, styles, 
					legendEntriesSearchVsControl, weights1, setLinY, drawSame);

//   //compare dR(W muon, soft muon) for events with mu+had mass > 0 and > 2
//   mergePlotsIn1File("/data1/yohay/Wh1_Medium/muHadAnalysisV8.root", 
// 		    "/data1/yohay/Wh1_Medium/dRWMuSoftMu_muHadMassGt0VsGt2_WMuonExcluded.root");

// //space-saving constant definitions
//   const string analysisFilePath("/data1/yohay/WJets/WJets_tau_analysis");
//   const string fileExt(".root");

//   //make individual efficiency plots for signal and Z-->mumu samples
//   vector<string> effInputFiles;
//   vector<string> comparisonInputFiles;
//   effInputFiles.push_back(analysisFilePath + controlSample + objTag + controlGenSel + trigger);
//   effInputFiles.push_back(analysisFilePath + signalSample + objTag + leg + signalGenSel + trigger);
//   for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
//        ++iFile) {
// //     comparisonInputFiles.push_back(*iFile + effTag + fileExt);
//     comparisonInputFiles.push_back(*iFile + fileExt);
// //     plotNice(*iFile + fileExt, jetAnalysisEffMap, binLabelMap, jetAnalysisMap, 
// // 	     *iFile + effTag + fileExt, "noPDF");
//   }

//   //compare Z-->mumu and signal muon efficiency
//   vector<string> canvasNames;
// //   canvasNames.push_back("eff_numeratorPT_over_denominatorPT");
// //   canvasNames.push_back("eff_numeratorEta_over_denominatorEta");
//   canvasNames.push_back("muEnergyFractionCanvas");
//   vector<string> graphNames;
// //   graphNames.push_back("divide_numeratorPT_by_denominatorPT");
// //   graphNames.push_back("divide_numeratorEta_by_denominatorEta");
//   graphNames.push_back("muEnergyFraction");
//   Color_t colors[2] = {kBlack, kRed};
//   Style_t styles[2] = {20, 21};
// //   drawMultipleEfficiencyGraphsOn1Canvas(analysisFilePath + "eff" + objTag + outputTag + 
// // 					fileExt, comparisonInputFiles, canvasNames, graphNames, 
// // 					colors, styles);
}
