import FWCore.ParameterSet.Config as cms

GMSBAnalyzer = cms.EDAnalyzer('GMSBAnalyzer',

##                               #input
##                               outputFile = cms.untracked.string("skim_analysis.root"),
##                               photonTag = cms.untracked.InputTag("photons", "", "RECOCleaned"),
##                               categoryTag = cms.untracked.InputTag("EDMCategoryProducer",
##                                                                    "eventCategory", "CATEGORY"
##                                                                    ),
##                               diEMETTag = cms.untracked.InputTag("EDMCategoryProducer",
##                                                                  "evtDiEMET", "CATEGORY"
##                                                                  ),
##                               photonSeedTimeTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "photonSeedTime", "OWNPARTICLES"
##     ),
##                               photonE2OverE9Tag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "photonE2OverE9", "OWNPARTICLES"
##     ),
##                               passAllTag = cms.untracked.InputTag("EDMCategoryProducer",
##                                                                   "passingPhotons", "CATEGORY"
##                                                                   ),
##                               trgResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
##                               tcMETTag = cms.untracked.InputTag("tcMet", "", "RECO"),
##                               passETMinTag = cms.untracked.InputTag("EDMCategoryProducer",
##                                                                     "passETMin", "OWNPARTICLES"
##                                                                     ),
##                               passAbsEtaMaxTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "passAbsEtaMax", "OWNPARTICLES"
##     ),
##                               passECALIsoMaxTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "passECALIsoMax", "OWNPARTICLES"
##     ),
##                               passHCALIsoMaxTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "passHCALIsoMax", "OWNPARTICLES"
##     ),
##                               passHOverEMaxTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "passHOverEMax", "OWNPARTICLES"
##     ),
##                               passTrackIsoMaxTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "passTrackIsoMax", "OWNPARTICLES"
##     ),
##                               passSigmaIetaIetaMaxTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "passSigmaIetaIetaMax", "OWNPARTICLES"
##     ),
##                               passAbsSeedTimeMaxTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "passAbsSeedTimeMax", "OWNPARTICLES"
##     ),
##                               passE2OverE9MaxTag = cms.untracked.InputTag(
##     "EDMCategoryProducer", "passE2OverE9Max", "OWNPARTICLES"
##     ),

##                               #triggers
##                               HLTProcessName = cms.untracked.string("HLT"),
##                               HLTPathGG = cms.untracked.vstring(
##     "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1"
##     ),
##                               HLTPathFF = cms.untracked.vstring(
##     "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1"
##     ),
##                               HLTPathAllPhotons = cms.untracked.vstring(
##     "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1"
##     )
)
