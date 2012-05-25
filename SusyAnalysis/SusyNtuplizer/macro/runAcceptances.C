{
  //compile
  gROOT->Reset();
  gROOT->LoadMacro("Acceptance.C++");

  //gsq_B scan binning
  const pair<Int_t, Int_t> squarkGluinoNBins(17, 17);
  const pair<Double_t, Double_t> squarkGluinoMinBin(350.0, 370.0);
  const pair<Double_t, Double_t> squarkGluinoMaxBin(2050.0, 2070.0);
  const Double_t mSquarkStepSize = (squarkGluinoMaxBin.first - squarkGluinoMinBin.first)/(Double_t)squarkGluinoNBins.first;
  const Double_t mSquarkStepSizeOver2 = mSquarkStepSize/2.0;
  const Double_t mGluino1StepSize = (squarkGluinoMaxBin.second - squarkGluinoMinBin.second)/(Double_t)squarkGluinoNBins.second;
  const Double_t mGluino1StepSizeOver2 = mGluino1StepSize/2.0;

  //gB scan binning
  const pair<Int_t, Int_t> binoGluinoNBins(20, 15);
  const pair<Double_t, Double_t> binoGluinoMinBin(-20.0, 290.0);
  const pair<Double_t, Double_t> binoGluinoMaxBin(980.0, 1040.0);
  const Double_t mBino1StepSize = (binoGluinoMaxBin.first - binoGluinoMinBin.first)/(Double_t)binoGluinoNBins.first;
  const Double_t mBino1StepSizeOver2 = mBino1StepSize/2.0;
  const Double_t mGluino2StepSize = (binoGluinoMaxBin.second - binoGluinoMinBin.second)/(Double_t)binoGluinoNBins.second;
  const Double_t mGluino2StepSizeOver2 = mGluino2StepSize/2.0;

  //WB scan binning
  const pair<Int_t, Int_t> winoBinoNBins(17, 21);
  const pair<Double_t, Double_t> winoBinoMinBin(102.5, -7.5);
  const pair<Double_t, Double_t> winoBinoMaxBin(527.5, 517.5);
  const Double_t mWinoStepSize = (winoBinoMaxBin.first - winoBinoMinBin.first)/(Double_t)winoBinoNBins.first;
  const Double_t mWinoStepSizeOver2 = mWinoStepSize/2.0;
  const Double_t mBino2StepSize = (winoBinoMaxBin.second - winoBinoMinBin.second)/(Double_t)winoBinoNBins.second;
  const Double_t mBino2StepSizeOver2 = mBino2StepSize/2.0;

  //binolikegrid2 scan binning
  const pair<Int_t, Int_t> squarkGluinoBinoOldNBins(21, 21);
  const pair<Double_t, Double_t> squarkGluinoBinoOldMinBin(360.0, 360.0);
  const pair<Double_t, Double_t> squarkGluinoBinoOldMaxBin(2040.0, 2040.0);
  const Double_t mSquark2StepSize = 
    (squarkGluinoBinoOldMaxBin.first - squarkGluinoBinoOldMinBin.first)/(Double_t)squarkGluinoBinoOldNBins.first;
  const Double_t mSquark2StepSizeOver2 = mSquark2StepSize/2.0;
  const Double_t mGluino3StepSize = 
    (squarkGluinoBinoOldMaxBin.second - squarkGluinoBinoOldMinBin.second)/(Double_t)squarkGluinoBinoOldNBins.second;
  const Double_t mGluino3StepSizeOver2 = mGluino3StepSize/2.0;

  //winolikegrid2 scan binning
  const pair<Int_t, Int_t> squarkGluinoWinoOldNBins(21, 21);
  const pair<Double_t, Double_t> squarkGluinoWinoOldMinBin(360.0, 360.0);
  const pair<Double_t, Double_t> squarkGluinoWinoOldMaxBin(2040.0, 2040.0);
  const Double_t mSquark3StepSize = 
    (squarkGluinoWinoOldMaxBin.first - squarkGluinoWinoOldMinBin.first)/(Double_t)squarkGluinoWinoOldNBins.first;
  const Double_t mSquark3StepSizeOver2 = mSquark3StepSize/2.0;
  const Double_t mGluino4StepSize = 
    (squarkGluinoWinoOldMaxBin.second - squarkGluinoWinoOldMinBin.second)/(Double_t)squarkGluinoWinoOldNBins.second;
  const Double_t mGluino4StepSizeOver2 = mGluino4StepSize/2.0;

  //binochigrids scan binning
  const pair<Int_t, Int_t> binoGluinoOldNBins(10, 24);
  const pair<Double_t, Double_t> binoGluinoOldMinBin(100, 120);
  const pair<Double_t, Double_t> binoGluinoOldMaxBin(1100, 2040);
  const Double_t mBino3StepSize = (binoGluinoOldMaxBin.first - binoGluinoOldMinBin.first)/(Double_t)binoGluinoOldNBins.first;
  const Double_t mBino3StepSizeOver2 = mBino3StepSize/2.0;
  const Double_t mGluino5StepSize = (binoGluinoOldMaxBin.second - binoGluinoOldMinBin.second)/(Double_t)binoGluinoOldNBins.second;
  const Double_t mGluino5StepSizeOver2 = mGluino5StepSize/2.0;

  //gsq_B acceptance files
  vector<string> squarkGluinoAccFiles;
//   string squarkGluinoPath("/uscms_data/d2/rpyohay/Spectra_gsq_B/");
//   for (UInt_t mSquark = (UInt_t)(squarkGluinoMinBin.first + mSquarkStepSizeOver2); 
//        mSquark <= (UInt_t)(squarkGluinoMaxBin.first - mSquarkStepSizeOver2); mSquark+=(UInt_t)mSquarkStepSize) {
//     for (UInt_t mGluino = (UInt_t)(squarkGluinoMinBin.second + mGluino1StepSizeOver2); 
// 	 mGluino <= (UInt_t)(squarkGluinoMaxBin.second - mGluino1StepSizeOver2); mGluino+=(UInt_t)mGluino1StepSize) {
//       stringstream file;
//       file << squarkGluinoPath << "acceptance_" << mSquark << "_" << mGluino << "_375.root";
//       squarkGluinoAccFiles.push_back(file.str().c_str());
//     }
//   }

  //gB acceptance files
  vector<string> binoGluinoAccFiles;
//   string binoGluinoPath("/uscms_data/d2/rpyohay/Spectra_gB/");
//   for (UInt_t mBino = (UInt_t)(binoGluinoMinBin.first + mBino1StepSizeOver2); 
//        mBino <= (UInt_t)(binoGluinoMaxBin.first - mBino1StepSizeOver2); mBino+=(UInt_t)mBino1StepSize) {
//     for (UInt_t mGluino = (UInt_t)(binoGluinoMinBin.second + mGluino2StepSizeOver2); 
// 	 mGluino <= (UInt_t)(binoGluinoMaxBin.second - mGluino2StepSizeOver2); mGluino+=(UInt_t)mGluino2StepSize) {
//       stringstream file;
//       file << binoGluinoPath << "acceptance_2000_" << mGluino << "_" << mBino << ".root";
//       binoGluinoAccFiles.push_back(file.str().c_str());
//     }
//   }

  //WB acceptance files
  vector<string> winoBinoAccFiles;
//   string winoBinoPath("/uscms_data/d2/rpyohay/Spectra_WB/");
//   for (UInt_t mWino = (UInt_t)(winoBinoMinBin.first + mWinoStepSizeOver2); 
//        mWino <= (UInt_t)(winoBinoMaxBin.first - mWinoStepSizeOver2); mWino+=(UInt_t)mWinoStepSize) {
//     for (UInt_t mBino = (UInt_t)(winoBinoMinBin.second + mBino2StepSizeOver2); 
// 	 mBino <= (UInt_t)(winoBinoMaxBin.second - mBino2StepSizeOver2); mBino+=(UInt_t)mBino2StepSize) {
//       stringstream file;
//       file << winoBinoPath << "acceptance_" << mWino << "_" << mBino << ".root";
//       winoBinoAccFiles.push_back(file.str().c_str());
//     }
//   }

  //binolikegrid2 acceptance files
  vector<string> squarkGluinoBinoOldAccFiles;
  string squarkGluinoBinoOldPath("/uscms_data/d2/rpyohay/binolikegrid2/");
  for (UInt_t mSquark = (UInt_t)(squarkGluinoBinoOldMinBin.first + mSquark2StepSizeOver2); 
       mSquark <= (UInt_t)(squarkGluinoBinoOldMaxBin.first - mSquark2StepSizeOver2); mSquark+=(UInt_t)mSquark2StepSize) {
    for (UInt_t mGluino = (UInt_t)(squarkGluinoBinoOldMinBin.second + mGluino3StepSizeOver2); 
	 mGluino <= (UInt_t)(squarkGluinoBinoOldMaxBin.second - mGluino3StepSizeOver2); mGluino+=(UInt_t)mGluino3StepSize) {
      stringstream file;
      file << squarkGluinoBinoOldPath << "acceptance_" << mSquark << "_" << mGluino << "_375.root";
      squarkGluinoBinoOldAccFiles.push_back(file.str().c_str());
    }
  }

  //winolikegrid2 acceptance files
  vector<string> squarkGluinoWinoOldAccFiles;
  string squarkGluinoWinoOldPath("/uscms_data/d2/rpyohay/winolikegrid2/");
  for (UInt_t mSquark = (UInt_t)(squarkGluinoWinoOldMinBin.first + mSquark3StepSizeOver2); 
       mSquark <= (UInt_t)(squarkGluinoWinoOldMaxBin.first - mSquark3StepSizeOver2); mSquark+=(UInt_t)mSquark3StepSize) {
    for (UInt_t mGluino = (UInt_t)(squarkGluinoWinoOldMinBin.second + mGluino4StepSizeOver2); 
	 mGluino <= (UInt_t)(squarkGluinoWinoOldMaxBin.second - mGluino4StepSizeOver2); mGluino+=(UInt_t)mGluino4StepSize) {
      stringstream file;
      file << squarkGluinoWinoOldPath << "acceptance_" << mSquark << "_" << mGluino << "_375.root";
      squarkGluinoWinoOldAccFiles.push_back(file.str().c_str());
    }
  }

  //binochigrids acceptance files
  vector<string> binoGluinoOldAccFiles;
  string binoGluinoOldPath("/uscms_data/d2/rpyohay/binochigrids/");
  for (UInt_t mBino = (UInt_t)(binoGluinoOldMinBin.first + mBino3StepSizeOver2); 
       mBino <= (UInt_t)(binoGluinoOldMaxBin.first - mBino3StepSizeOver2); mBino+=(UInt_t)mBino3StepSize) {
    for (UInt_t mGluino = (UInt_t)(binoGluinoOldMinBin.second + mGluino5StepSizeOver2); 
	 mGluino <= (UInt_t)(binoGluinoOldMaxBin.second - mGluino5StepSizeOver2); mGluino+=(UInt_t)mGluino5StepSize) {
      stringstream file;
      file << binoGluinoOldPath << "acceptance_2500_" << mGluino << "_" << mBino << ".root";
      binoGluinoOldAccFiles.push_back(file.str().c_str());
    }
  }

  //run acceptances
  writeAcceptanceFile("/uscms_data/d2/rpyohay/acc_old.root", squarkGluinoNBins, binoGluinoNBins, winoBinoNBins, squarkGluinoBinoOldNBins, 
		      squarkGluinoWinoOldNBins, binoGluinoOldNBins, squarkGluinoMinBin, binoGluinoMinBin, winoBinoMinBin, 
		      squarkGluinoBinoOldMinBin, squarkGluinoWinoOldMinBin, binoGluinoOldMinBin, squarkGluinoMaxBin, binoGluinoMaxBin, 
		      winoBinoMaxBin, squarkGluinoBinoOldMaxBin, squarkGluinoWinoOldMaxBin, binoGluinoOldMaxBin, squarkGluinoAccFiles, 
		      binoGluinoAccFiles, winoBinoAccFiles, squarkGluinoBinoOldAccFiles, squarkGluinoWinoOldAccFiles, 
		      binoGluinoOldAccFiles);

  //MET bins
  vector<float> METBins;
  METBins.push_back(50.0);
  METBins.push_back(60.0);
  METBins.push_back(70.0);
  METBins.push_back(80.0);
  METBins.push_back(100.0);

//   //print ranges of nGG rel. stat. error
//   sortNGGRelStatErr(squarkGluinoBinoOldAccFiles, squarkGluinoWinoOldAccFiles, binoGluinoOldAccFiles, METBins);

  //signal files for limit setting
  vector<string> outputROOTFiles;
  outputROOTFiles.
    push_back("/uscms/home/rpyohay/UserCode/LPCPJM/LimitSetting_Jan2012/inputHists.20120520/signal_contamination_bino_chi0375.root");
  outputROOTFiles.
    push_back("/uscms/home/rpyohay/UserCode/LPCPJM/LimitSetting_Jan2012/inputHists.20120520/signal_contamination_wino_chi0375.root");
  outputROOTFiles.
    push_back("/uscms/home/rpyohay/UserCode/LPCPJM/LimitSetting_Jan2012/inputHists.20120520/signal_contamination_bino_chi0.root");

  //no. generated events files for limit setting
  vector<string> outputTxtFiles;
  outputTxtFiles.
    push_back("/uscms/home/rpyohay/UserCode/LPCPJM/LimitSetting_Jan2012/xsecdat/mcAccMap_bino_mN375_20120520.dat");
  outputTxtFiles.
    push_back("/uscms/home/rpyohay/UserCode/LPCPJM/LimitSetting_Jan2012/xsecdat/mcAccMap_wino_mN375_20120520.dat");
  outputTxtFiles.
    push_back("/uscms/home/rpyohay/UserCode/LPCPJM/LimitSetting_Jan2012/xsecdat/mcAccMap_bino_mNScan_20120520.dat");

  //vector of vectors of acceptance files names (1 vector for each grid)
  vector<vector<string> > accFiles;
  accFiles.push_back(squarkGluinoBinoOldAccFiles);
  accFiles.push_back(squarkGluinoWinoOldAccFiles);
  accFiles.push_back(binoGluinoOldAccFiles);

  //scans
  vector<string> scans;
  scans.push_back("gsq_B");
  scans.push_back("gsq_W");
  scans.push_back("gB");

//   makeFilesForLimitSetting(outputROOTFiles, accFiles, scans, "nojet");

  makeNGenFiles(outputTxtFiles, accFiles, scans);
}
