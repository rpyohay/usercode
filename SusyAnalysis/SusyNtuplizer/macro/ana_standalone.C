// adapted from SusyAnalysis/SusyNtuplizer/macro/ana.C

#include "EventAnalyzer.h"

void ana_standalone() {

  gSystem->Load("libSusy.so");

  // // Look ../jec/JetMETObjects/README
  gSystem->Load("../jec/lib/libJetMETObjects.so");

  //GMSBTools shared library (contains Categorizer class)
  gSystem->Load("../../../GMSBTools/lib/libFilters.so");

  // Printing utility for ntuple variables
  gSystem->Load("SusyEventPrinter_cc.so");
//   gROOT->LoadMacro("SusyEventPrinter.cc++");

  // Main analysis code
  gSystem->SetIncludePath("-I../../..");
//   gROOT->LoadMacro("EventAnalyzer.cc++");
//   gROOT->LoadMacro("DongwookCategoryProducer.cc++");
  gSystem->Load("EventAnalyzer_cc.so");
  gSystem->Load("DongwookCategoryProducer_cc.so");

  /*now, ECAL/HCAL/track isolation cuts are meant to be matched with loose trigger cuts, so don't 
    use them (i.e. track isolation) to distinguish photon from fake!*/

  //configuration
  ParameterSet pars;
  pars.photonTag = "photons";
  pars.photon1ETMin = 40.0;
  pars.photon2ETMin = 25.0;
  pars.photonAbsEtaMax = 1.4442;
  pars.photonECALIsoMaxPTMultiplier = 0.012; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					       EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
					       cone*/
  pars.photonECALIsoMaxConstant = 6.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					 EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonECALIsoEffArea = 0.093; /*https://twiki.cern.ch/twiki/bin/view/CMS/
				      RA3IsolationConePileupCorrections, dR = 0.3 cone, 
				      18-Jan-12, default rho_EtaMax*/
  pars.photonHCALIsoMaxPTMultiplier = 0.005; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					       EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
					       cone*/
  pars.photonHCALIsoMaxConstant = 4.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					 EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonHCALIsoEffArea = 0.0281; /*https://twiki.cern.ch/twiki/bin/view/CMS/
					RA3IsolationConePileupCorrections, dR = 0.3 cone, 
					18-Jan-12, default rho_EtaMax*/
  pars.photonHOverEMax = 0.05;
  pars.photonR9Max = 1.0;
  pars.photonR9Min = 0.8;
  pars.photonTrackIsoMaxPTMultiplier = 0.002; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
						EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 
						cone*/
  pars.photonTrackIsoMaxConstant = 4.0; /*https://twiki.cern.ch/twiki/bin/viewauth/CMS/
					  EgammaWorkingPointsv3, IsoVL, assuming dR = 0.3 cone*/
  pars.photonCombinedIsoMax = 6.0;
  pars.fakeCombinedIsoMax = 20.0;
  pars.isoConeHLT = DR03;
  pars.isoConeOffline = DR03;
  pars.photonSigmaIetaIetaMax = 0.011;
  pars.photonHLTSigmaIetaIetaMax = 0.014;
  pars.photonAbsSeedTimeMax = -1.0;
  pars.photonE2OverE9Max = -1.0;
  pars.photonDPhiMin = 0.05;
  pars.photonDRMin = 0.6;
  pars.pixelVetoOnFake = true;
  pars.treeName = "susyTree";
  pars.input = vector<string>();
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_10_1_flU.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_1_1_Dz4.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_2_2_hJ4.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_3_1_RDQ.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_4_2_j8E.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_5_2_Fj5.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_6_1_6kI.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_7_2_JW1.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_8_2_Aam.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_10-25/susyEvent_9_1_ZtH.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_10_1_Kee.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_1_1_jAk.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_2_1_jxR.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_3_1_pAZ.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_4_1_H56.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_5_1_8XS.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_6_1_bZy.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_7_1_hti.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_8_2_Hyr.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_25-250/susyEvent_9_2_EHt.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_10_1_itX.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_11_1_kvU.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_1_1_J6F.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_2_1_iun.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_3_1_SPV.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_4_1_p43.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_5_1_KLb.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_6_1_JBN.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_7_1_D5Z.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_8_2_rs7.root");
//   pars.input.push_back("/data2/yohay/RA3/diphoton_250/susyEvent_9_1_Vvk.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_10_1_hb1.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_11_1_psH.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_12_1_EAt.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_13_1_ky3.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_14_1_fRr.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_15_1_FVF.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_16_1_TiA.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_17_1_kHJ.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_18_1_S3f.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_19_1_YT4.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_1_1_BlL.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_20_1_pio.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_21_1_qfx.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_22_1_MAR.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_23_1_wxS.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_24_1_wC0.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_25_1_9Cl.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_26_1_vGX.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_27_1_Kpa.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_28_1_9Mq.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_29_1_yO2.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_2_2_ZaE.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_30_1_Meq.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_31_1_QZi.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_32_1_Vfh.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_33_1_NIM.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_34_1_Bhi.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_35_1_uQ5.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_36_1_zN9.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_37_1_MnR.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_38_1_CS3.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_39_1_xm0.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_3_1_2EO.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_40_1_RcU.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_41_1_QAC.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_42_1_aFP.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_4_1_xkH.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_5_1_crh.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_6_1_Vuo.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_7_1_RBc.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_8_1_ZHW.root");
//   pars.input.push_back("/data2/yohay/RA3/Zee/susyEvent_9_1_cPm.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_10_1_Fai.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_11_1_4Xo.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_12_1_bAK.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_13_1_ANA.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_14_1_xXq.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_15_1_fay.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_16_1_r5s.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_17_1_naR.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_18_1_UvQ.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_19_1_MvE.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_1_1_lYG.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_2_2_ld4.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_3_1_LXc.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_4_3_QyT.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_5_1_jZV.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_6_1_Hr6.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_7_1_jdd.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_8_1_4mZ.root");
//   pars.input.push_back("/data2/yohay/RA3/tt/susyEvent_9_2_qkF.root");
  pars.HLT = vector<TString>();
  pars.HLT.push_back("HLT_Photon26_IsoVL_Photon18");
  pars.HLT.push_back("HLT_Photon36_CaloIdL_Photon22_CaloIdL");
  pars.HLT.push_back("HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL");
  pars.HLT.push_back("HLT_Photon36_CaloIdL_IsoVL_Photon22_R9Id");
  pars.HLT.push_back("HLT_Photon36_R9Id_Photon22_CaloIdL_IsoVL");
  pars.HLT.push_back("HLT_Photon36_R9Id_Photon22_R9Id");
  pars.nEvts = -1;
  pars.JSON = "";
  pars.outputFile = "/data2/yohay/RA3/diphoton_10-25/skim.root";
  pars.recategorize = false;

  TStopwatch ts;

  ts.Start();

  DongwookCategoryProducer producer(pars);

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
