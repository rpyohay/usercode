import FWCore.ParameterSet.Config as cms

SampleMaker = cms.EDFilter("SampleMaker",

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

                           trackPTMin = cms.double(15.0), #GeV
                           eTrackRMin = cms.double(0.8),
                           minDRPhotons = cms.double(0.8),
                           maxSeedTime = cms.double(3.0), #ns

                           #reco::Photon tag
                           photonTag = cms.InputTag("photons", "", "RECOCleaned"),

                           #reco::Track tag for general tracks
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
                           debugFileName = cms.untracked.string("SampleMaker_debug_ff_newCode_keepHalo_sameMuonEndcapReq.txt"),

                           #optional debug flag (debugging off by default)
                           debugFlag = cms.untracked.bool(True),

                           #optional flag to tell the code how to check halo coincidence
                           #checks against passing EB objects only by default; setting this flag false causes it to check against all EB reco::Photon objects
                           #checkHaloCoincidenceWithPassingEBCandsOnly = cms.untracked.bool(False)

                           #option flag to turn off halo rejection (on by default)
                           rejectHalo = cms.untracked.bool(False)
                           )
