{
  gSystem->Load("libRooFit.so");
  gSystem->Load("libRooFitCore.so");
  gSystem->Load("../../../PhysicsTools/lib/libTagAndProbe.so");
  gSystem->SetIncludePath("-I../../.. -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms7/include");
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

//   makeFinalMETPlots("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_finalEGMisIDRate_geq1j_nominalMETBins.root", "/Users/rachelyohay/RA3/data/final_MET_plots_geq1j.root");

  // calculateElectronPhotonMisIDRate("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_egMisIDRate.root", "/Users/rachelyohay/RA3/data/egMisIDRateFitv18.root");

  // vector<string> vars;
  // vars.push_back("PT");
  // // vars.push_back("Eta");

  // vector<string> formattedVars;
  // formattedVars.push_back("p_{T} (GeV)");
  // // formattedVars.push_back("#eta");

  // plotDependentEGMisIDRate("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_pTDepEGMisIDRateSyst.root", "/Users/rachelyohay/RA3/data/egMisIDRateVsPT.root", vars, formattedVars);

  // replotMET("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_finalEGMisIDRate_inclusive_nominalMETBins.root", "/Users/rachelyohay/RA3/data/egMET.root");

  vector<string> inputFileNames;
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_0.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_1.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_2.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_3.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_4.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_5.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_6.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_7.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_8.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_9.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_10.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_11.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_12.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_13.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_14.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_15.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_16.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_17.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_18.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_19.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_20.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_21.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_22.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_23.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_24.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_25.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_26.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_27.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_28.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_29.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_30.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_31.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_32.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_33.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_34.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_35.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_36.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_37.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_38.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_39.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_40.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_41.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_42.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_43.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_44.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_45.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_46.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_47.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_48.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_49.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_50.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_51.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_52.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_53.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_54.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_55.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_56.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_57.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_58.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_59.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_60.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_61.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_62.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_63.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_64.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_65.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_66.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_67.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_68.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_69.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_70.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_71.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_72.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_73.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_74.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_75.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_76.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_77.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_78.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_79.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_80.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_81.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_82.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_83.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_84.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_85.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_86.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_87.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_88.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_89.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_90.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_91.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_92.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_93.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_94.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_95.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_96.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_97.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_98.root");
  inputFileNames.push_back("/data2/yohay/RA3/JESUncertainty/4684pb-1_MET_18-Jan-12_skim_JESUncertainty_geq1j_nominalBins_99.root");

  vector<float> METBins;
  METBins.push_back(0.0);
  METBins.push_back(5.0);
  METBins.push_back(10.0);
  METBins.push_back(15.0);
  METBins.push_back(20.0);
  METBins.push_back(25.0);
  METBins.push_back(30.0);
  METBins.push_back(35.0);
  METBins.push_back(40.0);
  METBins.push_back(50.0);
  METBins.push_back(70.0);
  METBins.push_back(100.0);
  METBins.push_back(150.0);
  METBins.push_back(500.0);
  estimateJESUncertaintyOnBkg(inputFileNames, METBins, "/data2/yohay/RA3/JESUncertainty/nominalBins_geq1j.root");
}
