{
  gSystem->Load("libRooFit.so");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("../../../PhysicsTools/lib/libTagAndProbe.so");
  gSystem->SetIncludePath("-I../../..");
  gROOT->LoadMacro("plotCode.C++");
  // compareDataToMC();
  // compareEEToGG();
  // compare2Methods();
  // compare3Methods();
  // compareEEToGG1ToGG2();
  // compare3Samples();
  // vector<string> inFile;
  // inFile.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_UllaJetDefinition_finerGeq1jDijetPTBinning_v2.root");
  // inFile.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_allEE.root");
  // compareEEToFF(inFile, "/Users/rachelyohay/RA3/data/eeVsFF_data_ZeeVsAllEE.root");
  // inFile.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_UllaJetDefinition_finerGeq1jDijetPTBinning_v2.root");
  // inFile.
  //   push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_Type-IMETCorrections.root");
  // compareEEToFF(inFile, 
  // 		"/Users/rachelyohay/RA3/data/eeVsFF_data_uncorrectedVsType-ICorrectedMET.root");
  // inFile.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_rightPFElectronCleaning.root");
  // inFile.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_Yueh-FengDiEMPTBinsAndJetDefinition.root");
  // compareEEToFF(inFile, 
  // 		"/Users/rachelyohay/RA3/data/ggVsFF_data_UllaVsYueh-FengCleaning.root");
  // replotEMF("/Users/rachelyohay/RA3/data/4684pb-1_EMF_18-Jan-12_skim.root", 
  // 	    "/Users/rachelyohay/RA3/data/EMF.root");
  // replotMET("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_UllaJetDefinition_finerGeq1jDijetPTBinning_v2.root", "/Users/rachelyohay/RA3/data/all_ee_plots.root");
  // compareEEToFFForEERegions("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_UllaJetDefinition_finerGeq1jDijetPTBinning_v2.root", "/Users/rachelyohay/RA3/data/eeVsFF_data_ZeeVsSidebands.root");

  // vector<string> inputFiles;
  // inputFiles.
  //   push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_10GeVDiEMPTBins.root");
  // inputFiles.
  //   push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_50GeVDiEMPTBins.root");
  // inputFiles.
  //   push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_wideHighDiEMPTBins.root");
  // inputFiles.
  //   push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_Yueh-FengDiEMPTBins.root");
  // inputFiles.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_rightPFElectronCleaning.root");

  // vector<short> color;
  // color.push_back(kRed);
  // color.push_back(kBlue);
  // color.push_back(kMagenta);
  // color.push_back(kSpring);
  // color.push_back(kBlack);

  // vector<Style_t> markerStyle;
  // markerStyle.push_back(kFullDiamond);
  // markerStyle.push_back(kFullTriangleDown);
  // markerStyle.push_back(kFullTriangleUp);
  // markerStyle.push_back(kFullSquare);
  // markerStyle.push_back(kFullCircle);

  // vector<string> entryLabel;
  // entryLabel.push_back("10 GeV bins");
  // entryLabel.push_back("50 GeV bins");
  // entryLabel.push_back("Wide bin for high di-EM p_{T}");
  // entryLabel.push_back("CMS-PAS-SUS-12-001");
  // entryLabel.push_back("Nominal");

  // vector<Option_t*> drawOption;
  // drawOption.push_back("");
  // for (unsigned int i = 1; i <= 4; ++i) { drawOption.push_back("SAME"); }

  // compareDiEMPTBinSizes(inputFiles, color, markerStyle, entryLabel, drawOption, 
  // 			"/Users/rachelyohay/RA3/data/di-EM_pT_bin_size_comparison.root");

  // vector<unsigned int> binNums;
  // for (unsigned int iBin = 10; iBin <= 14; ++iBin) { binNums.push_back(iBin); }
  // printErrors("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_finalEGMisIDRate_inclusive_limitMETBins.root", "/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_finalEGMisIDRate_inclusive_limitMETBins_errors.txt", binNums);

  // estimateMEEBkgShapeErr();

  makeFinalMETPlots("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_finalEGMisIDRate_geq1j_nominalMETBins.root", "/Users/rachelyohay/RA3/data/final_MET_plots_geq1j.root");

  // calculateElectronPhotonMisIDRate("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_egMisIDRate.root", "/Users/rachelyohay/RA3/data/egMisIDRateFitv18.root");

  // vector<string> vars;
  // vars.push_back("PT");
  // // vars.push_back("Eta");

  // vector<string> formattedVars;
  // formattedVars.push_back("p_{T} (GeV)");
  // // formattedVars.push_back("#eta");

  // plotDependentEGMisIDRate("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_pTDepEGMisIDRateSyst.root", "/Users/rachelyohay/RA3/data/egMisIDRateVsPT.root", vars, formattedVars);

  // replotMET("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_finalEGMisIDRate_inclusive_nominalMETBins.root", "/Users/rachelyohay/RA3/data/egMET.root");

}
