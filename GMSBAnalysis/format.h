#include "TChain.h"

#define W_MASS             80.403
#define Z_MASS             91.1876

#define MAX_PHOTONS        100
#define MAX_JETS           100
#define MAX_TRACKS         1000
#define MAX_MUONS          100
#define MAX_COSMICS        100
#define MAX_RECHITS        30000
#define MAX_HEHITS         10000
#define MAX_VERTEXS        100
#define MAX_GENPARTICLES   10000
#define MAX_GENDAUGHTERS   10

#define N_TRIGGER_BOOKINGS 70  
char TriggerBooking[N_TRIGGER_BOOKINGS][50] = {
//photon triggers
"HLT_DoublePhoton5_L1R",                       //bit 0
"HLT_Photon10_L1R",                           // bit 1
"HLT_Photon10_Cleaned_L1R",                   // bit 2
"HLT_DoublePhoton10_L1R",                     // bit 3
"HLT_Photon15_L1R",                           // bit 
"HLT_Photon15_Cleaned_L1R",                   // bit 5
"HLT_Photon15_LooseEcalIso_Cleaned_L1R",      // bit 6
"HLT_Photon15_TrackIso_Cleaned_L1R",          // bit 7
"HLT_DoublePhoton15_L1R",                     // bit 8
"HLT_Photon20_L1R",                           // bit 9
"HLT_DoublePhoton20_L1R ",                    // bit 10
"HLT_Photon20_Cleaned_L1R",                   // bit 11
"HLT_Photon30_L1R",                           // bit 12
"HLT_Photon30_Cleaned_L1R",                   // bit 13
"HLT_Photon50_L1R",                           // bit 14
"HLT_Photon50_Cleaned_L1R",                   // bit 15
"lots",                                       // bit 16
"more",                                       // bit 17
"triggers",                                   // bit 18
"lots",                                       // bit 19
"more",                                       // bit 20
"triggers",                                   // bit 21
"lots",                                       // bit 22
"more",                                       // bit 23
"triggers",                                   // bit 24
"lots",                                       // bit 25
"more",                                       // bit 26
"triggers",                                   // bit 27
"lots",                                       // bit 28
"more",                                       // bit 29
"triggers",                                   // bit 30
"lots",                                       // bit 31
"more",                                       // bit 32
"triggers",                                   // bit 33
"lots",                                       // bit 34
"more",                                       // bit 35
"triggers",                                   // bit 36
"lots",                                       // bit 37
"more",                                       // bit 38
"triggers",                                   // bit 39
"lots",                                       // bit 40
"more",                                       // bit 41
"triggers",                                   // bit 42
"lots",                                       // bit 43
//other
"HLT_DoublePhoton4_eeRes_L1R",                // bit 44
"HLT_DoublePhoton5_CEP_L1R",                  // bit 45
//beam halo trigger
"HLT_CSCBeamHalo",                            // bit 46
//muon triggers
"HLT_L1Mu",                                   // bit 47
"HLT_L1MuOpen",                               // bit 48
"HLT_L2Mu0",                                  // bit 49
"HLT_L2Mu3",                                  // bit 50
"HLT_L2Mu5",                                  // bit 51
"HLT_L1MuOpen_DT",                            // bit 52
"HLT_L1Mu14_L1ETM30",                         // bit 53
"HLT_L1Mu14_L1SingleJet6U",                   // bit 54
"HLT_L1Mu14_L1SingleEG10",                    // bit 55
"HLT_L1Mu20",                                 // bit 56
"HLT_DoubleMu3",                              // bit 57
"HLT_Mu3",                                    // bit 58
"HLT_Mu5",                                    // bit 59
"HLT_Mu9",                                    // bit 60
"HLT_IsoMu3",                                 // bit 61
"HLT_L2Mu9",                                  // bit 62
"HLT_L2Mu11",                                 // bit 63
"HLT_Mu7",                                    // bit 64
//some speciall HLT for Summer10 and Spring10 MC
"HLT_Photon25_L1R",			      // bit 65
"HLT_Photon30_L1R_8E29",		      // bit 66
"HLT_Photon30_L1R_1E31",		      // bit 67
"HLT_Photon70_L1R",			      // bit 68
"HLT_DoublePhoton15_L1R"		      // bit 69
};

class EvtInfoBranches {
public:
  int   Run;
  int   Event;
  int   LumiBlk;
  int   BX;    	//bunchCrossing
  int   BO;   	//orbitNumber
  int   isMC;   //McFlag
  char  Proc[20];  //name of the process
  float weight;	//event weight for MC
  float MET;	//CaloMET
  float METx;
  float METy;
  float METphi;
  float pfMET;
  float pfMETphi;
  float tcMET;
  float tcMETphi;
  float SumEt;
  float METSignificance;
  float tcSumEt;
  float tcMETSignificance;
  float pfSumEt;
  float pfMETSignificance;
  int   TrgCount;	// No. of fired booking bits
  int   TrgBook[N_TRIGGER_BOOKINGS];	// Trigger bits

  void Register(TTree *root) {
	root->SetBranchAddress("EvtInfo_Run"	    , &Run);
	root->SetBranchAddress("EvtInfo_Event"	    , &Event);
        root->SetBranchAddress("EvtInfo_LumiBlk"      , &LumiBlk);
        root->SetBranchAddress("EvtInfo_BX"           , &BX);
        root->SetBranchAddress("EvtInfo_BO"           , &BO);
        root->SetBranchAddress("EvtInfo_isMC"         , &isMC);
        root->SetBranchAddress("EvtInfo_Proc"         , &Proc[0]);
        root->SetBranchAddress("EvtInfo_weight"       , &weight);
	root->SetBranchAddress("EvtInfo_MET"	      , &MET);
        root->SetBranchAddress("EvtInfo_METx"         , &METx);
        root->SetBranchAddress("EvtInfo_METy"         , &METy);
	root->SetBranchAddress("EvtInfo_METphi"	      , &METphi);
        root->SetBranchAddress("EvtInfo_pfMET"        , &pfMET);
        root->SetBranchAddress("EvtInfo_pfMETphi"     , &pfMETphi);
        root->SetBranchAddress("EvtInfo_tcMET"        , &tcMET);
        root->SetBranchAddress("EvtInfo_tcMETphi"     , &tcMETphi);
	root->SetBranchAddress("EvtInfo_SumEt"	    , &SumEt);
	root->SetBranchAddress("EvtInfo_METSignificance", &METSignificance);
	root->SetBranchAddress("EvtInfo_tcSumEt"	    , &tcSumEt);
	root->SetBranchAddress("EvtInfo_tcMETSignificance", &tcMETSignificance);
	root->SetBranchAddress("EvtInfo_pfSumEt"	    , &pfSumEt);
	root->SetBranchAddress("EvtInfo_pfMETSignificance", &pfMETSignificance);
	root->SetBranchAddress("EvtInfo_TrgCount"     , &TrgCount);
	root->SetBranchAddress("EvtInfo_TrgBook"      , &TrgBook[0]);
  }										    
};

class PhoInfoBranches {
public:
  int	Size; 
  int	Index[MAX_PHOTONS];
  float Px[MAX_PHOTONS];
  float Py[MAX_PHOTONS];
  float Pz[MAX_PHOTONS];
  float Pt[MAX_PHOTONS];
  float E[MAX_PHOTONS];
  float Eta[MAX_PHOTONS];
  float Phi[MAX_PHOTONS];
  float ScE[MAX_PHOTONS];
  float ScRawEnergy[MAX_PHOTONS];
  float ScEta[MAX_PHOTONS];
  float ScPhi[MAX_PHOTONS];
  float ScX[MAX_PHOTONS];
  float ScY[MAX_PHOTONS];
  float ScZ[MAX_PHOTONS];
  float HadOverEM[MAX_PHOTONS];
  float EcalIso[MAX_PHOTONS];
  float HcalIso[MAX_PHOTONS];
  float TrackIsoPtHol[MAX_PHOTONS];
  float TrackIsoPtSol[MAX_PHOTONS];
  int   nTrackHol[MAX_PHOTONS];
  int   nTrackSol[MAX_PHOTONS];
  float sigmaIetaIeta[MAX_PHOTONS];
  float r9[MAX_PHOTONS];
  int   hasPixelSeed[MAX_PHOTONS];
  int   isEBGap[MAX_PHOTONS];
  int   isEB[MAX_PHOTONS];
  int   isEBEEGap[MAX_PHOTONS];
  int   isConv[MAX_PHOTONS];
  float VtxX[MAX_PHOTONS];
  float VtxY[MAX_PHOTONS];
  float VtxZ[MAX_PHOTONS];
  float phiwid[MAX_PHOTONS];
  float etaphiwid[MAX_PHOTONS];
  float drminjet[MAX_PHOTONS];
  float roundness[MAX_PHOTONS];
  float angle[MAX_PHOTONS];
  int   ncry[MAX_PHOTONS];
  int   seedIeta[MAX_PHOTONS];
  int   seedIphi[MAX_PHOTONS];
  float seedE[MAX_PHOTONS];
  float seedTime[MAX_PHOTONS];

  void Register(TTree *root) {
	root->SetBranchAddress("PhoInfo_Size"               , &Size);
	root->SetBranchAddress("PhoInfo_Index"              , &Index[0]);
	root->SetBranchAddress("PhoInfo_Px"                 , &Px[0]);
        root->SetBranchAddress("PhoInfo_Py"                 , &Py[0]);
        root->SetBranchAddress("PhoInfo_Pz"                 , &Pz[0]);
        root->SetBranchAddress("PhoInfo_Pt"                 , &Pt[0]);
        root->SetBranchAddress("PhoInfo_E"                  , &E[0]);
	root->SetBranchAddress("PhoInfo_Eta"                , &Eta[0]);
	root->SetBranchAddress("PhoInfo_Phi"                , &Phi[0]);
        root->SetBranchAddress("PhoInfo_ScE"                , &ScE[0]);
        root->SetBranchAddress("PhoInfo_ScRawEnergy"        , &ScRawEnergy[0]);
        root->SetBranchAddress("PhoInfo_ScEta"              , &ScEta[0]);
        root->SetBranchAddress("PhoInfo_ScPhi"              , &ScPhi[0]);
        root->SetBranchAddress("PhoInfo_ScX"                , &ScX[0]);
        root->SetBranchAddress("PhoInfo_ScY"                , &ScY[0]);
        root->SetBranchAddress("PhoInfo_ScZ"                , &ScZ[0]);
        root->SetBranchAddress("PhoInfo_HadOverEM"          , &HadOverEM[0]);
        root->SetBranchAddress("PhoInfo_EcalIso"            , &EcalIso[0]);
        root->SetBranchAddress("PhoInfo_HcalIso"            , &HcalIso[0]);
        root->SetBranchAddress("PhoInfo_TrackIsoPtHol"      , &TrackIsoPtHol[0]);
        root->SetBranchAddress("PhoInfo_TrackIsoPtSol"      , &TrackIsoPtSol[0]);
        root->SetBranchAddress("PhoInfo_nTrackHol"          , &nTrackHol[0]);
        root->SetBranchAddress("PhoInfo_nTrackSol"          , &nTrackSol[0]);
        root->SetBranchAddress("PhoInfo_sigmaIetaIeta"      , &sigmaIetaIeta[0]);
        root->SetBranchAddress("PhoInfo_r9"                 , &r9[0]);
        root->SetBranchAddress("PhoInfo_hasPixelSeed"       , &hasPixelSeed[0]);
        root->SetBranchAddress("PhoInfo_isEBGap"            , &isEBGap[0]);
        root->SetBranchAddress("PhoInfo_isEB"               , &isEB[0]);
        root->SetBranchAddress("PhoInfo_isEBEEGap"          , &isEBEEGap[0]);
        root->SetBranchAddress("PhoInfo_isConv"             , &isConv[0]);
        root->SetBranchAddress("PhoInfo_VtxX"               , &VtxX[0]);
        root->SetBranchAddress("PhoInfo_VtxY"               , &VtxY[0]);
        root->SetBranchAddress("PhoInfo_VtxZ"               , &VtxZ[0]);
        root->SetBranchAddress("PhoInfo_phiwid"             , &phiwid[0]);
        root->SetBranchAddress("PhoInfo_etaphiwid"          , &etaphiwid[0]);
        root->SetBranchAddress("PhoInfo_drminjet"           , &drminjet[0]);
        root->SetBranchAddress("PhoInfo_roundness"          , &roundness[0]);
        root->SetBranchAddress("PhoInfo_angle"              , &angle[0]);
        root->SetBranchAddress("PhoInfo_ncry"               , &ncry[0]);
        root->SetBranchAddress("PhoInfo_seedIeta"           , &seedIeta[0]);
        root->SetBranchAddress("PhoInfo_seedIphi"           , &seedIphi[0]);
        root->SetBranchAddress("PhoInfo_seedE"              , &seedE[0]);
        root->SetBranchAddress("PhoInfo_seedTime"           , &seedTime[0]);
  }  
};

class JetInfoBranches {
public:
  int	Size; 
  int   Index[MAX_JETS];
  float Et[MAX_JETS];
  float Pt[MAX_JETS];
  float Eta[MAX_JETS];
  float Phi[MAX_JETS];
  int   nCaloTowers[MAX_JETS];
  float EMEnergyFraction[MAX_JETS];
  int   N90[MAX_JETS];
  int   N60[MAX_JETS];
  float fHPD[MAX_JETS];
  float fRBX[MAX_JETS];
  float rHF[MAX_JETS];
 
  void Register(TTree *root) {
	root->SetBranchAddress("JetInfo_Size"		   , &Size);
	root->SetBranchAddress("JetInfo_Index"		   , &Index[0]);
	root->SetBranchAddress("JetInfo_Et"		   , &Et[0]);
	root->SetBranchAddress("JetInfo_Pt"		   , &Pt[0]);
	root->SetBranchAddress("JetInfo_Eta"		   , &Eta[0]);
	root->SetBranchAddress("JetInfo_Phi"		   , &Phi[0]);
	root->SetBranchAddress("JetInfo_nCaloTowers"	   , &nCaloTowers[0]);
        root->SetBranchAddress("JetInfo_EMEnergyFraction"    , &EMEnergyFraction[0]);
        root->SetBranchAddress("JetInfo_N90"                 , &N90[0]);
        root->SetBranchAddress("JetInfo_N60"                 , &N60[0]);
	root->SetBranchAddress("JetInfo_fHPD"                , &fHPD[0]);
	root->SetBranchAddress("JetInfo_fRBX"                , &fRBX[0]);
	root->SetBranchAddress("JetInfo_rHF"                 , &rHF[0]);
  }  
};

class PATJetInfoBranches {
public:
  int   Size;
  int   Index[MAX_JETS];

  //uncorrected jet 4-vector
  float  uncorrJetET[MAX_JETS];
  float  uncorrJetPX[MAX_JETS];
  float  uncorrJetPY[MAX_JETS];
  float  uncorrJetPZ[MAX_JETS];
  float  uncorrJetEta[MAX_JETS];
  float  uncorrJetPhi[MAX_JETS];

  //L1 (L1 corrects for pileup, noise, etc.) corrected jet 4-vector
  float  L1CorrJetET[MAX_JETS];
  float  L1CorrJetPX[MAX_JETS];
  float  L1CorrJetPY[MAX_JETS];
  float  L1CorrJetPZ[MAX_JETS];
  float  L1CorrJetEta[MAX_JETS];
  float  L1CorrJetPhi[MAX_JETS];

  //L1 + L2 (L2 corrects for nonuniform jet response as a function of eta) corrected jet 4-vector
  float  L1L2CorrJetET[MAX_JETS];
  float  L1L2CorrJetPX[MAX_JETS];
  float  L1L2CorrJetPY[MAX_JETS];
  float  L1L2CorrJetPZ[MAX_JETS];
  float  L1L2CorrJetEta[MAX_JETS];
  float  L1L2CorrJetPhi[MAX_JETS];

  //L1 + L2 + L3 (L3 corrects for nonuniform jet response as a function of pT) corrected jet 4-vector
  //L1+L2+L3 is the default full correction
  float  L1L2L3CorrJetET[MAX_JETS];
  float  L1L2L3CorrJetPX[MAX_JETS];
  float  L1L2L3CorrJetPY[MAX_JETS];
  float  L1L2L3CorrJetPZ[MAX_JETS];
  float  L1L2L3CorrJetEta[MAX_JETS];
  float  L1L2L3CorrJetPhi[MAX_JETS];

  //L1+L2+L3 corrected jet information
  int    L1L2L3CorrJetN90[MAX_JETS]; //# calo towers carrying 90% of the jet's total energy
  int    L1L2L3CorrJetN60[MAX_JETS]; //# calo towers carrying 60% of the jet's total energy
  float  L1L2L3CorrJetEMF[MAX_JETS]; //EM energy fraction
  float  L1L2L3CorrJetHadEnergyInHB[MAX_JETS];
  float  L1L2L3CorrJetHadEnergyInHE[MAX_JETS];
  float  L1L2L3CorrJetHadEnergyInHO[MAX_JETS];
  float  L1L2L3CorrJetHadEnergyInHF[MAX_JETS];
  float  L1L2L3CorrJetEMEnergyInEB[MAX_JETS];
  float  L1L2L3CorrJetEMEnergyInEE[MAX_JETS];

  //MC matching
  float  L1L2L3CorrJetMatchedGenJetET[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetPX[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetPY[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetPZ[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetEta[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetPhi[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetEMEnergy[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetHadEnergy[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetInvisibleEnergy[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenJetOtherEnergy[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenPartonET[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenPartonPX[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenPartonPY[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenPartonPZ[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenPartonEta[MAX_JETS];
  float  L1L2L3CorrJetMatchedGenPartonPhi[MAX_JETS];
  int    L1L2L3CorrJetMatchedGenPartonPDGID[MAX_JETS];

  //jet energy corrections (JEC)
  float  JECUncorrToL1L2L3[MAX_JETS]; //factor applied to correct an uncorrected jet to an L1+L2+L3 corrected jet
  float  JECL1ToL1L2L3[MAX_JETS]; //factor applied to correct an L1 corrected jet to an L1+L2+L3 corrected jet
  float  JECL1L2ToL1L2L3[MAX_JETS]; //factor applied to correct an L1+L2 corrected jet to an L1+L2+L3 corrected jet

  //jet ID (computed for the UNCORRECTED jet)
  float  uncorrJetFHPD[MAX_JETS]; //fraction of jet energy carried by the hottest HPD
  float  uncorrJetFRBX[MAX_JETS]; //fraction of jet energy carried by the hottest RBX
  int    uncorrJetN90Hits[MAX_JETS]; //min # RecHits containing 90% of the jet energy
  int    uncorrJetNHCALTowers[MAX_JETS]; //# HCAL towers in the jet
  int    uncorrJetNECALTowers[MAX_JETS]; //# ECAL towers in the jet
  float  uncorrJetRHF[MAX_JETS]; //calculated in GMSBTuple.cc; something to do with HF noise rejection
  int    passesLOOSEPURE09JetID[MAX_JETS]; //1 if jet passed LOOSE PURE09 jet ID, 0 otherwise (https://twiki.cern.ch/twiki/bin/view/CMS/JetID)
  int    passesTIGHTPURE09JetID[MAX_JETS]; //1 if jet passed TIGHT PURE09 jet ID, 0 otherwise (https://twiki.cern.ch/twiki/bin/view/CMS/JetID)

  //cross-cleaning
  int    noPhotonOverlap09[MAX_JETS]; //true if the jet does not overlap with a photon within dR = 0.9 (M. Balazs's prescription)
  int    noPhotonOverlap05[MAX_JETS]; //true if the jet does not overlap with a photon within dR = 0.5 (R. Nandi's prescription)

  void Register(TTree *root, bool includeDetailedUncorrJetInfo, bool includeDetailedJetIDInfo, bool MC) {
	root->SetBranchAddress("PATJetInfo_Size"		      , &Size);
	root->SetBranchAddress("PATJetInfo_Index"		      , &Index[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetET"           , &L1L2L3CorrJetET[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetPX"           , &L1L2L3CorrJetPX[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetPY"           , &L1L2L3CorrJetPY[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetPZ"           , &L1L2L3CorrJetPZ[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetEta"          , &L1L2L3CorrJetEta[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetPhi"          , &L1L2L3CorrJetPhi[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetN90"          , &L1L2L3CorrJetN90[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetN60"          , &L1L2L3CorrJetN60[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetEMF"          , &L1L2L3CorrJetEMF[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetHadEnergyInHB", &L1L2L3CorrJetHadEnergyInHB[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetHadEnergyInHE", &L1L2L3CorrJetHadEnergyInHE[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetHadEnergyInHO", &L1L2L3CorrJetHadEnergyInHO[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetHadEnergyInHF", &L1L2L3CorrJetHadEnergyInHF[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetEMEnergyInEB" , &L1L2L3CorrJetEMEnergyInEB[0]);
	root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetEMEnergyInEE" , &L1L2L3CorrJetEMEnergyInEE[0]);
	root->SetBranchAddress("PATJetInfo_passesLOOSEPURE09JetID"    , &passesLOOSEPURE09JetID[0]);
	root->SetBranchAddress("PATJetInfo_passesTIGHTPURE09JetID"    , &passesTIGHTPURE09JetID[0]);
	root->SetBranchAddress("PATJetInfo_noPhotonOverlap09"         , &noPhotonOverlap09[0]);
	root->SetBranchAddress("PATJetInfo_noPhotonOverlap05"         , &noPhotonOverlap05[0]);

	if (includeDetailedUncorrJetInfo) {
	  root->SetBranchAddress("PATJetInfo_uncorrJetET"      , &uncorrJetET[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetPX"      , &uncorrJetPX[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetPY"      , &uncorrJetPY[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetPZ"      , &uncorrJetPZ[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetEta"     , &uncorrJetEta[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetPhi"     , &uncorrJetPhi[0]);
	  root->SetBranchAddress("PATJetInfo_L1CorrJetET"      , &L1CorrJetET[0]);
	  root->SetBranchAddress("PATJetInfo_L1CorrJetPX"      , &L1CorrJetPX[0]);
	  root->SetBranchAddress("PATJetInfo_L1CorrJetPY"      , &L1CorrJetPY[0]);
	  root->SetBranchAddress("PATJetInfo_L1CorrJetPZ"      , &L1CorrJetPZ[0]);
	  root->SetBranchAddress("PATJetInfo_L1CorrJetEta"     , &L1CorrJetEta[0]);
	  root->SetBranchAddress("PATJetInfo_L1CorrJetPhi"     , &L1CorrJetPhi[0]);
	  root->SetBranchAddress("PATJetInfo_L1L2CorrJetET"    , &L1L2CorrJetET[0]);
	  root->SetBranchAddress("PATJetInfo_L1L2CorrJetPX"    , &L1L2CorrJetPX[0]);
	  root->SetBranchAddress("PATJetInfo_L1L2CorrJetPY"    , &L1L2CorrJetPY[0]);
	  root->SetBranchAddress("PATJetInfo_L1L2CorrJetPZ"    , &L1L2CorrJetPZ[0]);
	  root->SetBranchAddress("PATJetInfo_L1L2CorrJetEta"   , &L1L2CorrJetEta[0]);
	  root->SetBranchAddress("PATJetInfo_L1L2CorrJetPhi"   , &L1L2CorrJetPhi[0]);
	  root->SetBranchAddress("PATJetInfo_JECUncorrToL1L2L3", &JECUncorrToL1L2L3[0]);
	  root->SetBranchAddress("PATJetInfo_JECL1ToL1L2L3"    , &JECL1ToL1L2L3[0]);
	  root->SetBranchAddress("PATJetInfo_JECL1L2ToL1L2L3"  , &JECL1L2ToL1L2L3[0]);
	}

	if (includeDetailedJetIDInfo) {
	  root->SetBranchAddress("PATJetInfo_uncorrJetFHPD"       , &uncorrJetFHPD[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetFRBX"       , &uncorrJetFRBX[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetN90Hits"    , &uncorrJetN90Hits[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetNHCALTowers", &uncorrJetNHCALTowers[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetNECALTowers", &uncorrJetNECALTowers[0]);
	  root->SetBranchAddress("PATJetInfo_uncorrJetRHF"        , &uncorrJetRHF[0]);
	}

	if (MC) {
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetET"             , &L1L2L3CorrJetMatchedGenJetET[0]);
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetPX"             , &L1L2L3CorrJetMatchedGenJetPX[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetPY"             , &L1L2L3CorrJetMatchedGenJetPY[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetPZ"             , &L1L2L3CorrJetMatchedGenJetPZ[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetEta"            , &L1L2L3CorrJetMatchedGenJetEta[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetPhi"            , &L1L2L3CorrJetMatchedGenJetPhi[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetEMEnergy"       , &L1L2L3CorrJetMatchedGenJetEMEnergy[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetHadEnergy"      , &L1L2L3CorrJetMatchedGenJetHadEnergy[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetInvisibleEnergy", &L1L2L3CorrJetMatchedGenJetInvisibleEnergy[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenJetOtherEnergy"    , &L1L2L3CorrJetMatchedGenJetOtherEnergy[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenPartonET"          , &L1L2L3CorrJetMatchedGenPartonET[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenPartonPX"          , &L1L2L3CorrJetMatchedGenPartonPX[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenPartonPY"          , &L1L2L3CorrJetMatchedGenPartonPY[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenPartonPZ"          , &L1L2L3CorrJetMatchedGenPartonPZ[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenPartonEta"         , &L1L2L3CorrJetMatchedGenPartonEta[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenPartonPhi"         , &L1L2L3CorrJetMatchedGenPartonPhi[0]); 
	  root->SetBranchAddress("PATJetInfo_L1L2L3CorrJetMatchedGenPartonPDGID"       , &L1L2L3CorrJetMatchedGenPartonPDGID[0]); 
	}
  }
};

class TrkInfoBranches {
public:
  int   Size;
  int   Index[MAX_TRACKS];
  int   Charge[MAX_TRACKS];
  float Pt[MAX_TRACKS];
  float Phi[MAX_TRACKS];
  float Eta[MAX_TRACKS];
  int   Nhit[MAX_TRACKS];
  float Chi2[MAX_TRACKS];
  float PtErr[MAX_TRACKS];
  float dxy[MAX_TRACKS];
  float dz[MAX_TRACKS];
  float outerPt[MAX_TRACKS];
  float outerPhi[MAX_TRACKS];
  float outerEta[MAX_TRACKS];
  float outerX[MAX_TRACKS];
  float outerY[MAX_TRACKS];
  float outerZ[MAX_TRACKS];

  void Register(TTree *root) {
        root->SetBranchAddress("TrkInfo_Size"           	, &Size);
        root->SetBranchAddress("TrkInfo_Index"          	, &Index[0]);
        root->SetBranchAddress("TrkInfo_Charge"         	, &Charge[0]);
        root->SetBranchAddress("TrkInfo_Pt"              	, &Pt[0]);
        root->SetBranchAddress("TrkInfo_Phi"              	, &Phi[0]);
        root->SetBranchAddress("TrkInfo_Eta"              	, &Eta[0]);
        root->SetBranchAddress("TrkInfo_Nhit"            	, &Nhit[0]);
        root->SetBranchAddress("TrkInfo_Chi2"            	, &Chi2[0]);
        root->SetBranchAddress("TrkInfo_PtErr"           	, &PtErr[0]);
        root->SetBranchAddress("TrkInfo_dxy"              	, &dxy[0]);
        root->SetBranchAddress("TrkInfo_dz"              	, &dz[0]);
        root->SetBranchAddress("TrkInfo_outerPt"         	, &outerPt[0]);
        root->SetBranchAddress("TrkInfo_outerPhi"         	, &outerPhi[0]);
        root->SetBranchAddress("TrkInfo_outerEta"         	, &outerEta[0]);
        root->SetBranchAddress("TrkInfo_outerX"           	, &outerX[0]);
        root->SetBranchAddress("TrkInfo_outerY"           	, &outerY[0]);
        root->SetBranchAddress("TrkInfo_outerZ"           	, &outerZ[0]);
  }
};

class MuonInfoBranches {
public:
  int   Size;
  int   Index[MAX_MUONS];
  int   Charge[MAX_MUONS];
  float Pt[MAX_MUONS];
  float Phi[MAX_MUONS];
  float Eta[MAX_MUONS];
  int   isStandAlone[MAX_MUONS];
  int   isGlobal[MAX_MUONS];
  int   Stations[MAX_MUONS];

  void Register(TTree *root) {
        root->SetBranchAddress("MuonInfo_Size"          	, &Size);
        root->SetBranchAddress("MuonInfo_Index"         	, &Index[0]);
        root->SetBranchAddress("MuonInfo_Charge"         	, &Charge[0]);
        root->SetBranchAddress("MuonInfo_Pt"              	, &Pt[0]);
        root->SetBranchAddress("MuonInfo_Phi"             	, &Phi[0]);
        root->SetBranchAddress("MuonInfo_Eta"             	, &Eta[0]);
        root->SetBranchAddress("MuonInfo_isStandAlone"    	, &isStandAlone[0]);
        root->SetBranchAddress("MuonInfo_isGlobal"        	, &isGlobal[0]);
        root->SetBranchAddress("MuonInfo_Stations"        	, &Stations[0]);
  }
};

class CosInfoBranches {
public:
  int   Size;
  int   Index[MAX_COSMICS];
  float Pt[MAX_COSMICS];
  float Eta[MAX_COSMICS];
  float Phi[MAX_COSMICS];
  int   Nhit[MAX_COSMICS];
  float Chi2[MAX_COSMICS];
  float PtErr[MAX_COSMICS];
  float X0[MAX_COSMICS];
  float Y0[MAX_COSMICS];
  float Z0[MAX_COSMICS];
  float Phi0[MAX_COSMICS];
  float Eta0[MAX_COSMICS];
  float Pt0[MAX_COSMICS];
  float Pz0[MAX_COSMICS];
  float Py0[MAX_COSMICS];
  float Px0[MAX_COSMICS];
  float XF[MAX_COSMICS];
  float YF[MAX_COSMICS];
  float ZF[MAX_COSMICS];
  float PhiF[MAX_COSMICS];
  float EtaF[MAX_COSMICS];
  float PtF[MAX_COSMICS];
  float PzF[MAX_COSMICS];
  float PyF[MAX_COSMICS];
  float PxF[MAX_COSMICS];

  void Register(TTree *root) {
        root->SetBranchAddress("CosInfo_Size"            	, &Size);
        root->SetBranchAddress("CosInfo_Index"           	, &Index[0]);
        root->SetBranchAddress("CosInfo_Pt"            		, &Pt[0]);
        root->SetBranchAddress("CosInfo_Eta"            	, &Eta[0]);
        root->SetBranchAddress("CosInfo_Phi"            	, &Phi[0]);
        root->SetBranchAddress("CosInfo_Nhit"            	, &Nhit[0]);
        root->SetBranchAddress("CosInfo_Chi2"            	, &Chi2[0]);
        root->SetBranchAddress("CosInfo_PtErr"           	, &PtErr[0]);
        root->SetBranchAddress("CosInfo_X0"            		, &X0[0]);
        root->SetBranchAddress("CosInfo_Y0"            		, &Y0[0]);
        root->SetBranchAddress("CosInfo_Z0"            		, &Z0[0]);
        root->SetBranchAddress("CosInfo_Phi0"            	, &Phi0[0]);
        root->SetBranchAddress("CosInfo_Eta0"            	, &Eta0[0]);
        root->SetBranchAddress("CosInfo_Pt0"            	, &Pt0[0]);
        root->SetBranchAddress("CosInfo_Pz0"            	, &Pz0[0]);
        root->SetBranchAddress("CosInfo_Py0"            	, &Py0[0]);
        root->SetBranchAddress("CosInfo_Px0"            	, &Px0[0]);
        root->SetBranchAddress("CosInfo_XF"            		, &XF[0]);
        root->SetBranchAddress("CosInfo_YF"            		, &YF[0]);
        root->SetBranchAddress("CosInfo_ZF"            		, &ZF[0]);
        root->SetBranchAddress("CosInfo_PhiF"            	, &PhiF[0]);
        root->SetBranchAddress("CosInfo_EtaF"            	, &EtaF[0]);
        root->SetBranchAddress("CosInfo_PtF"            	, &PtF[0]);
        root->SetBranchAddress("CosInfo_PzF"            	, &PzF[0]);
        root->SetBranchAddress("CosInfo_PyF"            	, &PyF[0]);
        root->SetBranchAddress("CosInfo_PxF"            	, &PxF[0]);

  }
};

class HitInfoBranches {
public:
  int   Size;
  int   Index[MAX_RECHITS];
  float Eta[MAX_RECHITS];
  float Phi[MAX_RECHITS];
  int   IEta[MAX_RECHITS];
  int   IPhi[MAX_RECHITS];
  float E[MAX_RECHITS];
  float Et[MAX_RECHITS];
  int   Clustered[MAX_RECHITS];
  int   Flag[MAX_RECHITS];
  float Time[MAX_RECHITS];

  void Register(TTree *root) {
        root->SetBranchAddress("HitInfo_Size"            	, &Size);
        root->SetBranchAddress("HitInfo_Index"           	, &Index[0]);
        root->SetBranchAddress("HitInfo_Eta"            	, &Eta[0]);
        root->SetBranchAddress("HitInfo_Phi"            	, &Phi[0]);
        root->SetBranchAddress("HitInfo_IEta"            	, &IEta[0]);
        root->SetBranchAddress("HitInfo_IPhi"            	, &IPhi[0]);
        root->SetBranchAddress("HitInfo_E"            		, &E[0]);
        root->SetBranchAddress("HitInfo_Et"            		, &Et[0]);
        root->SetBranchAddress("HitInfo_Clustered"      	, &Clustered[0]);
        root->SetBranchAddress("HitInfo_Flag"            	, &Flag[0]);
        root->SetBranchAddress("HitInfo_Time"            	, &Time[0]);

  }
};

class HEHitInfoBranches {
public:
  int   Size;
  int   Index[MAX_HEHITS];
  float Energy[MAX_HEHITS];
  float Et[MAX_HEHITS];
  float Time[MAX_HEHITS];
  int   Ieta[MAX_HEHITS];
  int   Iphi[MAX_HEHITS];
  int   Depth[MAX_HEHITS];
  float Phi[MAX_HEHITS];
  float Eta[MAX_HEHITS];
  float X[MAX_HEHITS];
  float Y[MAX_HEHITS];
  float Z[MAX_HEHITS];

  void Register(TTree *root) {
        root->SetBranchAddress("HEHitInfo_Size"          	, &Size);
        root->SetBranchAddress("HEHitInfo_Index"         	, &Index[0]);
        root->SetBranchAddress("HEHitInfo_Energy"       	, &Energy[0]);
        root->SetBranchAddress("HEHitInfo_Et"            	, &Et[0]);
        root->SetBranchAddress("HEHitInfo_Time"          	, &Time[0]);
        root->SetBranchAddress("HEHitInfo_Ieta"         	, &Ieta[0]);
        root->SetBranchAddress("HEHitInfo_Iphi"         	, &Iphi[0]);
        root->SetBranchAddress("HEHitInfo_Depth"        	, &Depth[0]);
        root->SetBranchAddress("HEHitInfo_Phi"            	, &Phi[0]);
        root->SetBranchAddress("HEHitInfo_Eta"            	, &Eta[0]);
        root->SetBranchAddress("HEHitInfo_X"            	, &X[0]);
        root->SetBranchAddress("HEHitInfo_Y"            	, &Y[0]);
        root->SetBranchAddress("HEHitInfo_Z"            	, &Z[0]);
  }
};

class VtxInfoBranches {
public:
  int   Size;
  int   Index[MAX_VERTEXS];
  float X[MAX_VERTEXS];
  float Y[MAX_VERTEXS];
  float Z[MAX_VERTEXS];
  float Chi2[MAX_VERTEXS];
  int   nTrack[MAX_VERTEXS];
  float TrkPtSum[MAX_VERTEXS];

  void Register(TTree *root) {
        root->SetBranchAddress("VtxInfo_Size"        		, &Size);
        root->SetBranchAddress("VtxInfo_Index"            	, &Index[0]);
        root->SetBranchAddress("VtxInfo_X"              	, &X[0]);
        root->SetBranchAddress("VtxInfo_Y"              	, &Y[0]);
        root->SetBranchAddress("VtxInfo_Z"              	, &Z[0]);
        root->SetBranchAddress("VtxInfo_Chi2"             	, &Chi2[0]);
        root->SetBranchAddress("VtxInfo_nTrack"           	, &nTrack[0]);
        root->SetBranchAddress("VtxInfo_TrkPtSum"         	, &TrkPtSum[0]);
  }
};

class GenInfoBranches {
public:
  int   Size;
  int   Index[MAX_GENPARTICLES];
  int   Charge[MAX_GENPARTICLES];
  int   Status[MAX_GENPARTICLES];
  int   pdgId[MAX_GENPARTICLES];
  float Pt[MAX_GENPARTICLES];
  float Pz[MAX_GENPARTICLES];
  float Mass[MAX_GENPARTICLES];
  float Phi[MAX_GENPARTICLES];
  float Eta[MAX_GENPARTICLES];
  float Vx[MAX_GENPARTICLES];
  float Vy[MAX_GENPARTICLES];
  float Vz[MAX_GENPARTICLES];
  int   nMothers[MAX_GENPARTICLES];
  int   nDaughters[MAX_GENPARTICLES];
  int   momIndex[MAX_GENPARTICLES];
  int   dauIndex[MAX_GENPARTICLES][MAX_GENDAUGHTERS];

  void Register(TTree *root) {
        root->SetBranchAddress("GenInfo_Size"             , &Size);
        root->SetBranchAddress("GenInfo_Index"            , &Index[0]);
        root->SetBranchAddress("GenInfo_Charge"           , &Charge[0]);
        root->SetBranchAddress("GenInfo_Status"           , &Status[0]);
        root->SetBranchAddress("GenInfo_pdgId"            , &pdgId[0]);
        root->SetBranchAddress("GenInfo_Pt"               , &Pt[0]);
        root->SetBranchAddress("GenInfo_Pz"               , &Pz[0]);
        root->SetBranchAddress("GenInfo_Mass"             , &Mass[0]);
        root->SetBranchAddress("GenInfo_Phi"              , &Phi[0]);
        root->SetBranchAddress("GenInfo_Eta"              , &Eta[0]);
        root->SetBranchAddress("GenInfo_Vx"               , &Vx[0]);
        root->SetBranchAddress("GenInfo_Vy"               , &Vy[0]);
        root->SetBranchAddress("GenInfo_Vz"               , &Vz[0]);
        root->SetBranchAddress("GenInfo_nMothers"         , &nMothers[0]);
        root->SetBranchAddress("GenInfo_nDaughters"       , &nDaughters[0]);
        root->SetBranchAddress("GenInfo_momIndex"         , &momIndex[0]);
        root->SetBranchAddress("GenInfo_dauIndex"         , &dauIndex[0][0]);
  }
};

