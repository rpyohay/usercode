import FWCore.ParameterSet.Config as cms

CutAnalyzer = cms.EDAnalyzer('CutAnalyzer',

                             #optional output file name (default is out.root)
                             outFileName = cms.untracked.string("out_ff_keepHalo.root"),

                             #1: gammagamma (candidate) sample
                             #2: egamma sample
                             #3: ee sample
                             #4: ff sample
                             #5: e + track sample
                             sampleType = cms.uint32(4),

                             #object (gamma, e, f) selection cut values

                             ECALIsoMaxPTMultiplierEB = cms.double(0.004),
                             ECALIsoMaxConstantEB = cms.double(4.2), #GeV
                             #ECALIsoMaxPTMultiplierEB = cms.double(0.0),
                             #ECALIsoMaxConstantEB = cms.double(999999999.0), #GeV
                             ECALIsoMaxPTMultiplierEE = cms.double(0.0),
                             ECALIsoMaxConstantEE = cms.double(999999999.0), #GeV
                             HCALIsoMaxPTMultiplierEB = cms.double(0.001),
                             HCALIsoMaxConstantEB = cms.double(2.2), #GeV
                             #HCALIsoMaxPTMultiplierEB = cms.double(0.0),
                             #HCALIsoMaxConstantEB = cms.double(999999999.0), #GeV
                             HCALIsoMaxPTMultiplierEE = cms.double(0.0),
                             HCALIsoMaxConstantEE = cms.double(999999999.0), #GeV
                             HOverEMaxPresel = cms.double(0.05),
                             #HOverEMaxPresel = cms.double(1.0),
                             ETMin = cms.double(10.0), #GeV

                             #1: EB
                             #2: EE
                             #3: all ECAL
                             fiducialRegion = cms.uint32(1),

                             #photon ID cuts
                             #useHOverE = cms.untracked.bool(True),
                             HOverEMax = cms.double(0.0),
                             #useSigmaEtaEta = cms.untracked.bool(True),
                             sigmaEtaEtaMax = cms.double(0.0),
                             useTrackIso = cms.untracked.bool(True),
                             trackIsoMaxPTMultiplier = cms.double(0.001), #GeV
                             trackIsoMaxConstant = cms.double(2.0), #GeV
                             #trackIsoMaxPTMultiplier = cms.double(0.0), #GeV
                             #trackIsoMaxConstant = cms.double(0.0), #GeV
                             minDRPhotons = cms.double(0.8),
                             maxSeedTime = cms.double(3.0), #ns

                             #option flag to turn off halo rejection (on by default)
                             rejectHalo = cms.untracked.bool(False),

                             #optional track selection cuts (default is -1.0 for both, meaning the track collection is not of interest)
                             #trackPTMin = cms.untracked.double(15.0), #GeV
                             #eTrackRMin = cms.untracked.double(0.8),

                             #reco::Photon tag
                             photonTag = cms.InputTag("photons", "", "RECOCleaned"),

                             #reco::Track tag
                             trackTag = cms.InputTag("generalTracks", "", "RECO"),

                             #reco::HBHERecHit tag
                             HBHERecHitTag = cms.InputTag("hbhereco"),
                             
                             #reco::Track tag for cosmic tracks
                             cosmicTrackTag = cms.InputTag("cosmicMuons"),

                             #EB RecHit tag
                             EBRecHitTag = cms.InputTag("ecalRecHit", "EcalRecHitsEB", "RECO"),

                             #EE RecHit tag
                             EERecHitTag = cms.InputTag("ecalRecHit", "EcalRecHitsEE", "RECO"),

                             #optional debug file name (file name is debug.txt by default)
                             debugFileName = cms.untracked.string("CutAnalyzer_debug_ff_keepHalo.txt"),

                             #optional debug flag (debugging off by default)
                             debugFlag = cms.untracked.bool(True)
                             )
