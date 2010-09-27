{
  gROOT->Reset();
  gROOT->LoadMacro("/afs/cern.ch/user/y/yohay/UserCode/yohay/GMSBAnalysis/EventSelector.cc++");
  gROOT->LoadMacro("/afs/cern.ch/user/y/yohay/UserCode/yohay/GMSBAnalysis/analysis.C++");

  //input files holding trees
  vector<string> fileList;
  //fileList.push_back("rfio:/castor/cern.ch/user/h/heyburn/EGSD_June14_JSON_vSept13.root");
  fileList.push_back("rfio:/castor/cern.ch/user/h/heyburn/June14ReReco_JSON_vSept13.root");
  //fileList.push_back("rfio:/castor/cern.ch/user/h/heyburn/gmsbNtuple_EGpromptv4Aug18_JSON_vSept13.root");
  //fileList.push_back("rfio:/castor/cern.ch/user/h/heyburn/July16ReReco_JSON_vSept13.root");
  fileList.push_back("/data/yohay/gmsbNtuple_EGpromptv4Sept2_JSON_vSept13.root");

  //desired HLT selection
  vector<int> HLTBits;
  HLTBits.push_back(12);
  HLTBits.push_back(13);

  //run!
  runSampleMaker("/data/yohay/etrack_test.root", fileList, "etrack_cfg.txt", HLTBits, 0, 99999);
  /*runSampleMaker("/data/yohay/etrack_2.root", fileList, "etrack_cfg.txt", HLTBits, 100000, 199999);
  runSampleMaker("/data/yohay/etrack_3.root", fileList, "etrack_cfg.txt", HLTBits, 200000, 299999);
  runSampleMaker("/data/yohay/etrack_4.root", fileList, "etrack_cfg.txt", HLTBits, 300000, 399999);
  runSampleMaker("/data/yohay/etrack_5.root", fileList, "etrack_cfg.txt", HLTBits, 400000, 499999);
  runSampleMaker("/data/yohay/etrack_6.root", fileList, "etrack_cfg.txt", HLTBits, 500000, 599999);
  runSampleMaker("/data/yohay/etrack_7.root", fileList, "etrack_cfg.txt", HLTBits, 600000, 699999);
  runSampleMaker("/data/yohay/etrack_8.root", fileList, "etrack_cfg.txt", HLTBits, 700000, 799999);
  runSampleMaker("/data/yohay/etrack_9.root", fileList, "etrack_cfg.txt", HLTBits, 800000, 899999);
  runSampleMaker("/data/yohay/etrack_10.root", fileList, "etrack_cfg.txt", HLTBits, 900000, 999999);*/
}
