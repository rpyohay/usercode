import FWCore.ParameterSet.Config as cms

CutAnalyzer = cms.EDAnalyzer('CutAnalyzer',

                             #optional output file name (default is out.root)
                             outFileName = cms.untracked.string("/data/yohay/GGM/tan_beta_10/RECO/skim/out_GGM_HLT15_pT20.root"),
                             #outFileName = cms.untracked.string("out_photon_jet_pT15_HLT15_pT20.root"),

                             #1: gammagamma (candidate) sample
                             #2: egamma sample
                             #3: ee sample
                             #4: ff sample
                             #5: e + track sample
                             sampleType = cms.uint32(4),

                             #object (gamma, e, f) selection cut values
                             
                             ECALIsoMaxPTMultiplierEB = cms.double(0.0),
                             ECALIsoMaxConstantEB = cms.double(999999999.0), #GeV
                             ECALIsoMaxPTMultiplierEE = cms.double(0.0),
                             ECALIsoMaxConstantEE = cms.double(999999999.0), #GeV
                             HCALIsoMaxEB = cms.double(999999999.0), #GeV
                             HCALIsoMaxEE = cms.double(999999999.0), #GeV
                             HOverEMaxPresel = cms.double(1.0),
                             ETMin = cms.double(20.0), #GeV

                             #1: EB
                             #2: EE
                             #3: all ECAL
                             fiducialRegion = cms.uint32(1),

                             #unused at the moment
                             HOverEMax = cms.double(0.0),

                             sigmaEtaEtaMax = cms.double(0.0),
                             trackIsoMax = cms.double(0.0), #GeV

                             #optional track selection cuts (default is -1.0 for both, meaning the track collection is not of interest)
                             #trackPTMin = cms.untracked.double(15.0), #GeV
                             #eTrackRMin = cms.untracked.double(0.8),

                             #reco::Photon tag
                             photonTag = cms.InputTag("photons"),

                             #reco::Track tag
                             trackTag = cms.InputTag("generalTracks", "", "RECO"),

                             #optional debug file name (file name is debug.txt by default)
                             debugFileName = cms.untracked.string("CutAnalyzer_debug_GGM_HLT15_pT20.txt"),

                             #optional debug flag (debugging off by default)
                             debugFlag = cms.untracked.bool(True)                             
)
