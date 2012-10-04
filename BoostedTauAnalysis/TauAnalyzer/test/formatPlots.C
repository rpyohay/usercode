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

//space-saving constant definitions
  const string analysisFilePath("/data1/yohay/");
  const string objTag("");
  const string suffix(".root");
  const string effTag("_eff");
  const string controlSample("WJetsToLNu_");
  const string signalSample("Wh_");
  const string leg("");
  const string trigger("jet_analysis");
  const string signalGenSel("muHad_");
  const string controlGenSel("");
  const string outputTag("_" + trigger);

  //make individual efficiency plots for signal and Z-->mumu samples
  vector<string> effInputFiles;
  vector<string> comparisonInputFiles;
  effInputFiles.push_back(analysisFilePath + controlSample + objTag + controlGenSel + trigger);
  effInputFiles.push_back(analysisFilePath + signalSample + objTag + leg + signalGenSel + trigger);
  for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
       ++iFile) {
//     comparisonInputFiles.push_back(*iFile + effTag + suffix);
    comparisonInputFiles.push_back(*iFile + suffix);
//     plotNice(*iFile + suffix, jetAnalysisEffMap, binLabelMap, jetAnalysisMap, 
// 	     *iFile + effTag + suffix, "noPDF");
  }

  //compare Z-->mumu and signal muon efficiency
  vector<string> canvasNames;
//   canvasNames.push_back("eff_numeratorPT_over_denominatorPT");
//   canvasNames.push_back("eff_numeratorEta_over_denominatorEta");
  canvasNames.push_back("muEnergyFractionCanvas");
  vector<string> graphNames;
//   graphNames.push_back("divide_numeratorPT_by_denominatorPT");
//   graphNames.push_back("divide_numeratorEta_by_denominatorEta");
  graphNames.push_back("muEnergyFraction");
  Color_t colors[2] = {kBlack, kRed};
  Style_t styles[2] = {20, 21};
  drawMultipleEfficiencyGraphsOn1Canvas(analysisFilePath + "eff" + objTag + outputTag + 
					suffix, comparisonInputFiles, canvasNames, graphNames, 
					colors, styles);
}
