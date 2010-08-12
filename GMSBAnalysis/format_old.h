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

#define N_TRIGGER_BOOKINGS 32  
char TriggerBooking[N_TRIGGER_BOOKINGS][38] = {
//single photon triggers
"HLT_Photon10_L1R",                           // bit 0
"HLT_Photon10_Cleaned_L1R",                   // bit 1
"HLT_Photon15_L1R",                           // bit 2
"HLT_Photon15_Cleaned_L1R",                   // bit 3
"HLT_Photon15_LooseEcalIso_Cleaned_L1R",      // bit 4
"HLT_Photon15_TrackIso_Cleaned_L1R",          // bit 5
"HLT_Photon20_L1R",                           // bit 6
"HLT_Photon20_Cleaned_L1R",                   // bit 7
"HLT_Photon30_L1R",                           // bit 8
//double photon triggers
"HLT_DoublePhoton4_eeRes_L1R",                // bit 9
"HLT_DoublePhoton5_L1R",                      // bit 10
"HLT_DoublePhoton5_CEP_L1R",                  // bit 11
"HLT_DoublePhoton10_L1R",                     // bit 12
//beam halo trigger
"HLT_CSCBeamHalo",                            // bit 13
//muon triggers
"HLT_L1Mu",                                   // bit 14
"HLT_L1MuOpen",                               // bit 15
"HLT_L2Mu0",                                  // bit 16
"HLT_L2Mu3",                                  // bit 17
"HLT_L2Mu5",                                  // bit 18
"HLT_L1MuOpen_DT",                            // bit 19
"HLT_L1Mu14_L1ETM30",                         // bit 20
"HLT_L1Mu14_L1SingleJet6U",                   // bit 21
"HLT_L1Mu14_L1SingleEG10",                    // bit 22
"HLT_L1Mu20",                                 // bit 23
"HLT_DoubleMu3",                              // bit 24
"HLT_Mu3",                                    // bit 25
"HLT_Mu5",                                    // bit 26
"HLT_Mu9",                                    // bit 27
"HLT_IsoMu3",                                 // bit 28
"HLT_L2Mu9",                                  // bit 29
"HLT_L2Mu11",                                 // bit 30
"HLT_Mu7"                                     // bit 31
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

  template<typename T>
    void Register(T *root) {
	root->SetBranchAddress("EvtInfo.Run"	    , &Run);
	root->SetBranchAddress("EvtInfo.Event"	    , &Event);
        root->SetBranchAddress("EvtInfo.LumiBlk"      , &LumiBlk);
        root->SetBranchAddress("EvtInfo.BX"           , &BX);
        root->SetBranchAddress("EvtInfo.BO"           , &BO);
        root->SetBranchAddress("EvtInfo.isMC"         , &isMC);
        //root->SetBranchAddress("EvtInfo.Proc"         , &Proc[0]);
        root->SetBranchAddress("EvtInfo.weight"       , &weight);
	root->SetBranchAddress("EvtInfo.MET"	      , &MET);
        root->SetBranchAddress("EvtInfo.METx"         , &METx);
        root->SetBranchAddress("EvtInfo.METy"         , &METy);
	root->SetBranchAddress("EvtInfo.METphi"	      , &METphi);
        root->SetBranchAddress("EvtInfo.pfMET"        , &pfMET);
        root->SetBranchAddress("EvtInfo.pfMETphi"     , &pfMETphi);
        root->SetBranchAddress("EvtInfo.tcMET"        , &tcMET);
        root->SetBranchAddress("EvtInfo.tcMETphi"     , &tcMETphi);
	root->SetBranchAddress("EvtInfo.SumEt"	    , &SumEt);
	root->SetBranchAddress("EvtInfo.METSignificance", &METSignificance);
	root->SetBranchAddress("EvtInfo.tcSumEt"	    , &tcSumEt);
	root->SetBranchAddress("EvtInfo.tcMETSignificance", &tcMETSignificance);
	root->SetBranchAddress("EvtInfo.pfSumEt"	    , &pfSumEt);
	root->SetBranchAddress("EvtInfo.pfMETSignificance", &pfMETSignificance);
	root->SetBranchAddress("EvtInfo.TrgCount"     , &TrgCount);
	root->SetBranchAddress("EvtInfo.TrgBook"      , &TrgBook[0]);
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

  template<typename T>
  void Register(T *root) {
	root->SetBranchAddress("PhoInfo.Size"               , &Size);
	root->SetBranchAddress("PhoInfo.Index"              , &Index[0]);
	root->SetBranchAddress("PhoInfo.Px"                 , &Px[0]);
        root->SetBranchAddress("PhoInfo.Py"                 , &Py[0]);
        root->SetBranchAddress("PhoInfo.Pz"                 , &Pz[0]);
        root->SetBranchAddress("PhoInfo.Pt"                 , &Pt[0]);
        root->SetBranchAddress("PhoInfo.E"                  , &E[0]);
	root->SetBranchAddress("PhoInfo.Eta"                , &Eta[0]);
	root->SetBranchAddress("PhoInfo.Phi"                , &Phi[0]);
        root->SetBranchAddress("PhoInfo.ScE"                , &ScE[0]);
        root->SetBranchAddress("PhoInfo.ScRawEnergy"        , &ScRawEnergy[0]);
        root->SetBranchAddress("PhoInfo.ScEta"              , &ScEta[0]);
        root->SetBranchAddress("PhoInfo.ScPhi"              , &ScPhi[0]);
        root->SetBranchAddress("PhoInfo.ScX"                , &ScX[0]);
        root->SetBranchAddress("PhoInfo.ScY"                , &ScY[0]);
        root->SetBranchAddress("PhoInfo.ScZ"                , &ScZ[0]);
        root->SetBranchAddress("PhoInfo.HadOverEM"          , &HadOverEM[0]);
        root->SetBranchAddress("PhoInfo.EcalIso"            , &EcalIso[0]);
        root->SetBranchAddress("PhoInfo.HcalIso"            , &HcalIso[0]);
        root->SetBranchAddress("PhoInfo.TrackIsoPtHol"      , &TrackIsoPtHol[0]);
        root->SetBranchAddress("PhoInfo.TrackIsoPtSol"      , &TrackIsoPtSol[0]);
        root->SetBranchAddress("PhoInfo.nTrackHol"          , &nTrackHol[0]);
        root->SetBranchAddress("PhoInfo.nTrackSol"          , &nTrackSol[0]);
        root->SetBranchAddress("PhoInfo.sigmaIetaIeta"      , &sigmaIetaIeta[0]);
        root->SetBranchAddress("PhoInfo.r9"                 , &r9[0]);
        root->SetBranchAddress("PhoInfo.hasPixelSeed"       , &hasPixelSeed[0]);
        root->SetBranchAddress("PhoInfo.isEBGap"            , &isEBGap[0]);
        root->SetBranchAddress("PhoInfo.isEB"               , &isEB[0]);
        root->SetBranchAddress("PhoInfo.isEBEEGap"          , &isEBEEGap[0]);
        root->SetBranchAddress("PhoInfo.isConv"             , &isConv[0]);
        root->SetBranchAddress("PhoInfo.VtxX"               , &VtxX[0]);
        root->SetBranchAddress("PhoInfo.VtxY"               , &VtxY[0]);
        root->SetBranchAddress("PhoInfo.VtxZ"               , &VtxZ[0]);
        root->SetBranchAddress("PhoInfo.phiwid"             , &phiwid[0]);
        root->SetBranchAddress("PhoInfo.etaphiwid"          , &etaphiwid[0]);
        root->SetBranchAddress("PhoInfo.drminjet"           , &drminjet[0]);
        root->SetBranchAddress("PhoInfo.roundness"          , &roundness[0]);
        root->SetBranchAddress("PhoInfo.angle"              , &angle[0]);
        root->SetBranchAddress("PhoInfo.ncry"               , &ncry[0]);
        root->SetBranchAddress("PhoInfo.seedIeta"           , &seedIeta[0]);
        root->SetBranchAddress("PhoInfo.seedIphi"           , &seedIphi[0]);
        root->SetBranchAddress("PhoInfo.seedE"              , &seedE[0]);
        //root->SetBranchAddress("PhoInfo.seedTime"           , &seedTime[0]);
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
  
  template<typename T>
  void Register(T *root) {
	root->SetBranchAddress("JetInfo.Size"		   , &Size);
	root->SetBranchAddress("JetInfo.Index"		   , &Index[0]);
	root->SetBranchAddress("JetInfo.Et"		   , &Et[0]);
	root->SetBranchAddress("JetInfo.Pt"		   , &Pt[0]);
	root->SetBranchAddress("JetInfo.Eta"		   , &Eta[0]);
	root->SetBranchAddress("JetInfo.Phi"		   , &Phi[0]);
	root->SetBranchAddress("JetInfo.nCaloTowers"	   , &nCaloTowers[0]);
        root->SetBranchAddress("JetInfo.EMEnergyFraction"    , &EMEnergyFraction[0]);
        root->SetBranchAddress("JetInfo.N90"                 , &N90[0]);
        root->SetBranchAddress("JetInfo.N60"                 , &N60[0]);
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

  template<typename T>
  void Register(T *root) {
        root->SetBranchAddress("TrkInfo.Size"           	, &Size);
        root->SetBranchAddress("TrkInfo.Index"          	, &Index[0]);
        root->SetBranchAddress("TrkInfo.Charge"         	, &Charge[0]);
        root->SetBranchAddress("TrkInfo.Pt"              	, &Pt[0]);
        root->SetBranchAddress("TrkInfo.Phi"              	, &Phi[0]);
        root->SetBranchAddress("TrkInfo.Eta"              	, &Eta[0]);
        root->SetBranchAddress("TrkInfo.Nhit"            	, &Nhit[0]);
        root->SetBranchAddress("TrkInfo.Chi2"            	, &Chi2[0]);
        root->SetBranchAddress("TrkInfo.PtErr"           	, &PtErr[0]);
        root->SetBranchAddress("TrkInfo.dxy"              	, &dxy[0]);
        root->SetBranchAddress("TrkInfo.dz"              	, &dz[0]);
        root->SetBranchAddress("TrkInfo.outerPt"         	, &outerPt[0]);
        root->SetBranchAddress("TrkInfo.outerPhi"         	, &outerPhi[0]);
        root->SetBranchAddress("TrkInfo.outerEta"         	, &outerEta[0]);
        root->SetBranchAddress("TrkInfo.outerX"           	, &outerX[0]);
        root->SetBranchAddress("TrkInfo.outerY"           	, &outerY[0]);
        root->SetBranchAddress("TrkInfo.outerZ"           	, &outerZ[0]);
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

  template<typename T>
  void Register(T *root) {
        root->SetBranchAddress("MuonInfo.Size"          	, &Size);
        root->SetBranchAddress("MuonInfo.Index"         	, &Index[0]);
        root->SetBranchAddress("MuonInfo.Charge"         	, &Charge[0]);
        root->SetBranchAddress("MuonInfo.Pt"              	, &Pt[0]);
        root->SetBranchAddress("MuonInfo.Phi"             	, &Phi[0]);
        root->SetBranchAddress("MuonInfo.Eta"             	, &Eta[0]);
        root->SetBranchAddress("MuonInfo.isStandAlone"    	, &isStandAlone[0]);
        root->SetBranchAddress("MuonInfo.isGlobal"        	, &isGlobal[0]);
        root->SetBranchAddress("MuonInfo.Stations"        	, &Stations[0]);
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

  template<typename T>
  void Register(T *root) {
        root->SetBranchAddress("CosInfo.Size"            	, &Size);
        root->SetBranchAddress("CosInfo.Index"           	, &Index[0]);
        root->SetBranchAddress("CosInfo.Pt"            		, &Pt[0]);
        root->SetBranchAddress("CosInfo.Eta"            	, &Eta[0]);
        root->SetBranchAddress("CosInfo.Phi"            	, &Phi[0]);
        root->SetBranchAddress("CosInfo.Nhit"            	, &Nhit[0]);
        root->SetBranchAddress("CosInfo.Chi2"            	, &Chi2[0]);
        root->SetBranchAddress("CosInfo.PtErr"           	, &PtErr[0]);
        root->SetBranchAddress("CosInfo.X0"            		, &X0[0]);
        root->SetBranchAddress("CosInfo.Y0"            		, &Y0[0]);
        root->SetBranchAddress("CosInfo.Z0"            		, &Z0[0]);
        root->SetBranchAddress("CosInfo.Phi0"            	, &Phi0[0]);
        root->SetBranchAddress("CosInfo.Eta0"            	, &Eta0[0]);
        root->SetBranchAddress("CosInfo.Pt0"            	, &Pt0[0]);
        root->SetBranchAddress("CosInfo.Pz0"            	, &Pz0[0]);
        root->SetBranchAddress("CosInfo.Py0"            	, &Py0[0]);
        root->SetBranchAddress("CosInfo.Px0"            	, &Px0[0]);
        root->SetBranchAddress("CosInfo.XF"            		, &XF[0]);
        root->SetBranchAddress("CosInfo.YF"            		, &YF[0]);
        root->SetBranchAddress("CosInfo.ZF"            		, &ZF[0]);
        root->SetBranchAddress("CosInfo.PhiF"            	, &PhiF[0]);
        root->SetBranchAddress("CosInfo.EtaF"            	, &EtaF[0]);
        root->SetBranchAddress("CosInfo.PtF"            	, &PtF[0]);
        root->SetBranchAddress("CosInfo.PzF"            	, &PzF[0]);
        root->SetBranchAddress("CosInfo.PyF"            	, &PyF[0]);
        root->SetBranchAddress("CosInfo.PxF"            	, &PxF[0]);

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

  template<typename T>
  void Register(T *root) {
        root->SetBranchAddress("HitInfo.Size"            	, &Size);
        root->SetBranchAddress("HitInfo.Index"           	, &Index[0]);
        root->SetBranchAddress("HitInfo.Eta"            	, &Eta[0]);
        root->SetBranchAddress("HitInfo.Phi"            	, &Phi[0]);
        root->SetBranchAddress("HitInfo.IEta"            	, &IEta[0]);
        root->SetBranchAddress("HitInfo.IPhi"            	, &IPhi[0]);
        root->SetBranchAddress("HitInfo.E"            		, &E[0]);
        root->SetBranchAddress("HitInfo.Et"            		, &Et[0]);
        root->SetBranchAddress("HitInfo.Clustered"      	, &Clustered[0]);
        root->SetBranchAddress("HitInfo.Flag"            	, &Flag[0]);
        root->SetBranchAddress("HitInfo.Time"            	, &Time[0]);

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

  template<typename T>
  void Register(T *root) {
        root->SetBranchAddress("HEHitInfo.Size"          	, &Size);
        root->SetBranchAddress("HEHitInfo.Index"         	, &Index[0]);
        root->SetBranchAddress("HEHitInfo.Energy"       	, &Energy[0]);
        root->SetBranchAddress("HEHitInfo.Et"            	, &Et[0]);
        root->SetBranchAddress("HEHitInfo.Time"          	, &Time[0]);
        root->SetBranchAddress("HEHitInfo.Ieta"         	, &Ieta[0]);
        root->SetBranchAddress("HEHitInfo.Iphi"         	, &Iphi[0]);
        root->SetBranchAddress("HEHitInfo.Depth"        	, &Depth[0]);
        root->SetBranchAddress("HEHitInfo.Phi"            	, &Phi[0]);
        root->SetBranchAddress("HEHitInfo.Eta"            	, &Eta[0]);
        root->SetBranchAddress("HEHitInfo.X"            	, &X[0]);
        root->SetBranchAddress("HEHitInfo.Y"            	, &Y[0]);
        root->SetBranchAddress("HEHitInfo.Z"            	, &Z[0]);
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

  template<typename T>
  void Register(T *root) {
        root->SetBranchAddress("VtxInfo.Size"        		, &Size);
        root->SetBranchAddress("VtxInfo.Index"            	, &Index[0]);
        root->SetBranchAddress("VtxInfo.X"              	, &X[0]);
        root->SetBranchAddress("VtxInfo.Y"              	, &Y[0]);
        root->SetBranchAddress("VtxInfo.Z"              	, &Z[0]);
        root->SetBranchAddress("VtxInfo.Chi2"             	, &Chi2[0]);
        root->SetBranchAddress("VtxInfo.nTrack"           	, &nTrack[0]);
        root->SetBranchAddress("VtxInfo.TrkPtSum"         	, &TrkPtSum[0]);
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

  template<typename T>
  void Register(T *root) {
        root->SetBranchAddress("GenInfo.Size"             , &Size);
        root->SetBranchAddress("GenInfo.Index"            , &Index[0]);
        root->SetBranchAddress("GenInfo.Charge"           , &Charge[0]);
        root->SetBranchAddress("GenInfo.Status"           , &Status[0]);
        root->SetBranchAddress("GenInfo.pdgId"            , &pdgId[0]);
        root->SetBranchAddress("GenInfo.Pt"               , &Pt[0]);
        root->SetBranchAddress("GenInfo.Pz"               , &Pz[0]);
        root->SetBranchAddress("GenInfo.Mass"             , &Mass[0]);
        root->SetBranchAddress("GenInfo.Phi"              , &Phi[0]);
        root->SetBranchAddress("GenInfo.Eta"              , &Eta[0]);
        root->SetBranchAddress("GenInfo.Vx"               , &Vx[0]);
        root->SetBranchAddress("GenInfo.Vy"               , &Vy[0]);
        root->SetBranchAddress("GenInfo.Vz"               , &Vz[0]);
        root->SetBranchAddress("GenInfo.nMothers"         , &nMothers[0]);
        root->SetBranchAddress("GenInfo.nDaughters"       , &nDaughters[0]);
        root->SetBranchAddress("GenInfo.momIndex"         , &momIndex[0]);
        root->SetBranchAddress("GenInfo.dauIndex"         , &dauIndex[0][0]);
  }
};

