import FWCore.ParameterSet.Config as cms

SkimAnalyzer = cms.EDAnalyzer('SkimAnalyzer',

                              #input
#                              outputFile = cms.untracked.string("skim_analysis.root"),
#                              photonTag = cms.untracked.InputTag(
#    "photons", "", "RECOCleaned"
#    ),
#                              categoryTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "eventCategory", "CATEGORY"
#    ),
#                              diEMETTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "evtDiEMET", "CATEGORY"
#    ),
##                               photonSeedTimeTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "photonSeedTime", "OWNPARTICLES"
##     ),
##                               photonE2OverE9Tag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "photonE2OverE9", "OWNPARTICLES"
##     ),
#                              trgResultsTag = cms.untracked.InputTag(
#    "TriggerResults", "", "HLT"
#    ),
#                              tcMETTag = cms.untracked.InputTag(
#    "tcMet", "", "RECO"
#    ),
#                              ETMinScanDecisionTag = cms.untracked.InputTag(
#    "CutScanProducer", "passETMin", "OWNPARTICLES"
#    ),
#                              ECALIsoMaxScanDecisionTag = cms.untracked.InputTag(
#    "CutScanProducer", "passECALIsoMax", "OWNPARTICLES"
#    ),
#                              HCALIsoMaxScanDecisionTag = cms.untracked.InputTag(
#    "CutScanProducer", "passHCALIsoMax", "OWNPARTICLES"
#    ),
#                              HOverEMaxScanDecisionTag = cms.untracked.InputTag(
#    "CutScanProducer", "passHOverEMax", "OWNPARTICLES"
#    ),
#			      trackIsoMaxScanDecisionTag = cms.untracked.InputTag(
#    "CutScanProducer", "passTrackIsoMax", "OWNPARTICLES"
#    ),
#                              sigmaIetaIetaMaxScanDecisionTag = cms.untracked.InputTag(
#    "CutScanProducer", "passSigmaIetaIetaMax", "OWNPARTICLES"
#    ),
#                              passETMinTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "passETMin", "OWNPARTICLES"
#    ),
#                              passAbsEtaMaxTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "passAbsEtaMax", "OWNPARTICLES"
#    ),
#                              passECALIsoMaxTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "passECALIsoMax", "OWNPARTICLES"
#    ),
#                              passHCALIsoMaxTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "passHCALIsoMax", "OWNPARTICLES"
#    ),
#                              passHOverEMaxTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "passHOverEMax", "OWNPARTICLES"
#    ),
#                              passTrackIsoMaxTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "passTrackIsoMax", "OWNPARTICLES"
#    ),
#                              passSigmaIetaIetaMaxTag = cms.untracked.InputTag(
#    "EDMCategoryProducer", "passSigmaIetaIetaMax", "OWNPARTICLES"
#    ),
##                              passAbsSeedTimeMaxTag = cms.untracked.InputTag(
##    "EDMCategoryProducer", "passAbsSeedTimeMax", "OWNPARTICLES"
##    ),
##                              passE2OverE9MaxTag = cms.untracked.InputTag(
##    "EDMCategoryProducer", "passE2OverE9Max", "OWNPARTICLES"
##    ),

                              #scan points
                              ETMinScan = cms.vdouble(35.0, 43.0, 50.0), #GeV
                              ECALIsoMaxPTMultiplierScan = cms.vdouble(0.012, 0.012, 0.012, 0.006),
                              ECALIsoMaxConstantScan = cms.vdouble(6.0, 5.5, 5.0, 4.2), #GeV
                              HCALIsoMaxPTMultiplierScan = cms.vdouble(0.005, 0.005, 0.005,
                                                                       0.0025),
                              HCALIsoMaxConstantScan = cms.vdouble(4.0, 3.5, 3.0, 2.2), #GeV
                              HOverEMaxScan = cms.vdouble(0.1, 0.05, 0.03, 0.05),
                              trackIsoMaxPTMultiplierScan = cms.vdouble(0.002, 0.002, 0.002,
                                                                        0.001),
                              trackIsoMaxConstantScan = cms.vdouble(4.0, 3.5, 3.0, 2.0), #GeV
                              sigmaIetaIetaMaxScan = cms.vdouble(0.024, 0.013, 0.011, 0.009)#,

                              #triggers
#                              HLTProcessName = cms.untracked.string("HLT"),
#                              HLTPathGG = cms.untracked.vstring(
#    "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1"
#    ),
#                              HLTPathFF = cms.untracked.vstring(
#    "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1"
#    )
)
