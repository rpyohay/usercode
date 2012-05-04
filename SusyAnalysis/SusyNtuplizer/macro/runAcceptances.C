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

  //gsq_B acceptance files
  vector<string> squarkGluinoAccFiles;
  string squarkGluinoPath("/uscms_data/d2/rpyohay/Spectra_gsq_B/");
  for (UInt_t mSquark = (UInt_t)(squarkGluinoMinBin.first + mSquarkStepSizeOver2); 
       mSquark <= (UInt_t)(squarkGluinoMaxBin.first - mSquarkStepSizeOver2); mSquark+=(UInt_t)mSquarkStepSize) {
    for (UInt_t mGluino = (UInt_t)(squarkGluinoMinBin.second + mGluino1StepSizeOver2); 
	 mGluino <= (UInt_t)(squarkGluinoMaxBin.second - mGluino1StepSizeOver2); mGluino+=(UInt_t)mGluino1StepSize) {
      stringstream file;
      file << squarkGluinoPath << "acceptance_" << mSquark << "_" << mGluino << "_375.root";
      squarkGluinoAccFiles.push_back(file.str().c_str());
    }
  }

  //gB acceptance files
  vector<string> binoGluinoAccFiles;
  string binoGluinoPath("/uscms_data/d2/rpyohay/Spectra_gB/");
  for (UInt_t mBino = (UInt_t)(binoGluinoMinBin.first + mBino1StepSizeOver2); 
       mBino <= (UInt_t)(binoGluinoMaxBin.first - mBino1StepSizeOver2); mBino+=(UInt_t)mBino1StepSize) {
    for (UInt_t mGluino = (UInt_t)(binoGluinoMinBin.second + mGluino2StepSizeOver2); 
	 mGluino <= (UInt_t)(binoGluinoMaxBin.second - mGluino2StepSizeOver2); mGluino+=(UInt_t)mGluino2StepSize) {
      stringstream file;
      file << binoGluinoPath << "acceptance_2000_" << mGluino << "_" << mBino << ".root";
      binoGluinoAccFiles.push_back(file.str().c_str());
    }
  }

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

  //run acceptances
  writeAcceptanceFile("/uscms_data/d2/rpyohay/acc_squarkGluino_binoGluino_METGeq50GeV.root", squarkGluinoNBins, binoGluinoNBins, 
		      winoBinoNBins, squarkGluinoMinBin, binoGluinoMinBin, winoBinoMinBin, squarkGluinoMaxBin, binoGluinoMaxBin, 
		      winoBinoMaxBin, squarkGluinoAccFiles, binoGluinoAccFiles, winoBinoAccFiles);
}
