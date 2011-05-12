import FWCore.ParameterSet.Config as cms

EDMCategoryProducer = cms.EDProducer('EDMCategoryProducer'#,

                                     #input
                                     #photonTag = cms.untracked.InputTag("photons", "", "RECOCleaned"),
                                     #recHitTagEB = cms.untracked.InputTag("reducedEcalRecHitsEB", "", "RECOCleaned"),
                                     #recHitTagEE = cms.untracked.InputTag("reducedEcalRecHitsEE", "", "RECOCleaned"),

                                     #cuts
                                     #photonETMin = cms.untracked.double(30.0), #GeV
                                     #photonAbsEtaMax = cms.untracked.double(1.379),
                                     #photonECALIsoMaxPTMultiplier = cms.untracked.double(0.006),
                                     #photonECALIsoMaxConstant = cms.untracked.double(4.2), #GeV
                                     #photonHCALIsoMaxPTMultiplier = cms.untracked.double(0.0025),
                                     #photonHCALIsoMaxConstant = cms.untracked.double(2.2), #GeV
                                     #photonHOverEMax = cms.untracked.double(0.05),
                                     #photonTrackIsoMaxPTMultiplier = cms.untracked.double(0.001),
                                     #photonTrackIsoMaxPTConstant = cms.untracked.double(2.0), #GeV
                                     #photonSigmaIetaIetaMax = cms.untracked.double(0.013),
                                     #photonAbsSeedTimeMax = cms.untracked.double(3.0), #ns
                                     #photonE2OverE9Max = cms.untracked.double(0.95),
                                     #photonDPhiMin = cms.untracked.double(0.05),
                                     #channelStatuses = cms.untracked.vuint32(0)
)
