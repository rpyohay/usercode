import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#desired INFO messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('Error')
process.MessageLogger.categories.append('Statistics')
process.MessageLogger.categories.append('SkimAnalyzer')
#process.MessageLogger.categories.append('EDMCategoryProducer')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

#events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#input
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/y/yohay/421/skim_studies/asymDoublePhotonSkim.root'
    ),
                            skipEvents = cms.untracked.uint32(0)
)

#global tag
process.load('Geometry.CaloEventSetup.CaloTopology_cfi')
process.load('Geometry.CaloEventSetup.CaloGeometry_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('GR_P_V17::All')

#stuff you need for some unknown reason
process.load('Configuration/StandardSequences/EndOfProcess_cff')

#di-photon filter
process.load('EGamma/EGammaSkims/diPhotonFilter_cfi')
process.diPhotonFilter.photonFirstEtCut = cms.double(30.0) #GeV
process.diPhotonFilter.photonSecondEtCut = cms.double(30.0) #GeV

#trigger analyzer
process.load('TriggerAnalysis/TriggerAnalyzer/triggeranalyzer_cfi')
process.TriggerAnalyzer.outputFile = cms.untracked.string(
    "/data2/yohay/skim_studies/CS_test_trigger_analysis_35GeV-30GeV.root"
    )
process.TriggerAnalyzer.unprescaledHLTPaths = cms.untracked.vstring(
    'HLT_DoublePhoton33', 'HLT_DoublePhoton5_IsoVL_CEP', 'HLT_Photon125_NoSpikeFilter',
    'HLT_Photon20_CaloIdVL_IsoL', 'HLT_Photon20_R9Id_Photon18_R9Id',
    'HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL',
    'HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id', 'HLT_Photon26_CaloIdL_IsoVL_Photon18',
    'HLT_Photon26_IsoVL_Photon18_IsoVL', 'HLT_Photon26_IsoVL_Photon18', 'HLT_Photon26_Photon18',
    'HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL', 'HLT_Photon30_CaloIdVL_IsoL',
    'HLT_Photon30_CaloIdVL', 'HLT_Photon32_CaloIdL_Photon26_CaloIdL',
    'HLT_Photon36_CaloIdL_Photon22_CaloIdL', 'HLT_Photon50_CaloIdVL_IsoL',
    'HLT_Photon75_CaloIdVL_IsoL', 'HLT_Photon75_CaloIdVL'
    )

#cut scan producer
process.load('Skims/CutScanProducer/cutscanproducer_cfi')
process.CutScanProducer.photonTag = cms.untracked.InputTag("photons", "", "RECO")
process.CutScanProducer.extETMin = cms.untracked.double(35.0) #GeV
process.CutScanProducer.HOverEMaxScan = cms.vdouble(0.1, 0.07, 0.05, 0.03)
process.CutScanProducer.trackIsoMaxPTMultiplierScan = cms.vdouble(0.005, 0.002, 0.002, 0.001)
process.CutScanProducer.trackIsoMaxConstantScan = cms.vdouble(17.5, 3.5, 3.0, 2.0) #GeV

#category producer
process.load('GMSBTools/Producers/edmcategoryproducer_cfi')
process.EDMCategoryProducer.photonTag = cms.untracked.InputTag("photons", "", "RECO")
process.EDMCategoryProducer.recHitTagEB = cms.untracked.InputTag(
    "reducedEcalRecHitsEB", "", "RECO"
    )
process.EDMCategoryProducer.recHitTagEE = cms.untracked.InputTag(
    "reducedEcalRecHitsEE", "", "RECO"
    )

#skim analyzer module
process.load('Skims/SkimAnalyzer/skimanalyzer_cfi')
process.SkimAnalyzer.outputFile = cms.untracked.string(
    '/data2/yohay/skim_studies/CS_test_analysis_35GeV-30GeV.root'
    )
process.SkimAnalyzer.photonTag = cms.untracked.InputTag("photons", "", "RECO")
process.SkimAnalyzer.categoryTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "eventCategory", "ANALYSIS"
    )
process.SkimAnalyzer.diEMETTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "evtDiEMET", "ANALYSIS"
    )
process.SkimAnalyzer.photonSeedTimeTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "photonSeedTime", "ANALYSIS"
    )
process.SkimAnalyzer.photonE2OverE9Tag = cms.untracked.InputTag(
    "EDMCategoryProducer", "photonE2OverE9", "ANALYSIS"
    )
process.SkimAnalyzer.ETMinScanDecisionTag = cms.untracked.InputTag(
    "CutScanProducer", "passETMin", "ANALYSIS"
    )
process.SkimAnalyzer.ECALIsoMaxScanDecisionTag = cms.untracked.InputTag(
    "CutScanProducer", "passECALIsoMax", "ANALYSIS"
    )
process.SkimAnalyzer.HCALIsoMaxScanDecisionTag = cms.untracked.InputTag(
    "CutScanProducer", "passHCALIsoMax", "ANALYSIS"
    )
process.SkimAnalyzer.HOverEMaxScanDecisionTag = cms.untracked.InputTag(
    "CutScanProducer", "passHOverEMax", "ANALYSIS"
    )
process.SkimAnalyzer.trackIsoMaxScanDecisionTag = cms.untracked.InputTag(
    "CutScanProducer", "passTrackIsoMax", "ANALYSIS"
    )
process.SkimAnalyzer.sigmaIetaIetaMaxScanDecisionTag = cms.untracked.InputTag(
    "CutScanProducer", "passSigmaIetaIetaMax", "ANALYSIS"
    )
process.SkimAnalyzer.passETMinTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passETMin", "ANALYSIS"
    )
process.SkimAnalyzer.passAbsEtaMaxTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passAbsEtaMax", "ANALYSIS"
    )
process.SkimAnalyzer.passECALIsoMaxTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passECALIsoMax", "ANALYSIS"
    )
process.SkimAnalyzer.passHCALIsoMaxTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passHCALIsoMax", "ANALYSIS"
    )
process.SkimAnalyzer.passHOverEMaxTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passHOverEMax", "ANALYSIS"
    )
process.SkimAnalyzer.passTrackIsoMaxTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passTrackIsoMax", "ANALYSIS"
    )
process.SkimAnalyzer.passSigmaIetaIetaMaxTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passSigmaIetaIetaMax", "ANALYSIS"
    )
process.SkimAnalyzer.HOverEMaxScan = cms.vdouble(0.1, 0.07, 0.05, 0.03)
process.SkimAnalyzer.trackIsoMaxPTMultiplierScan = cms.vdouble(0.005, 0.002, 0.002, 0.001)
process.SkimAnalyzer.trackIsoMaxConstantScan = cms.vdouble(17.5, 3.5, 3.0, 2.0) #GeV
process.SkimAnalyzer.HLTPathGG = cms.untracked.vstring(
    "HLT_Photon32_CaloIdL_Photon26_CaloIdL"
    )
process.SkimAnalyzer.HLTPathFF = cms.untracked.vstring(
    "HLT_Photon32_CaloIdL_Photon26_CaloIdL"
    )

#event content to save
process.eventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'drop *',
    
    #trigger
    'keep *_hltL1GtObjectMap_*_*',
    'keep *_TriggerResults_*_*',
    'keep *_hltTriggerSummaryAOD_*_*',
    'keep *_gtDigis_*_*',
    'keep *_scalersRawToDigi_*_*',
    'keep *_l1extraParticles_*_*',

    #photons (and PF)
    'keep *_selectDigi_*_*',
    'keep *_reducedEcalRecHitsEB_*_*',
    'keep *_reducedEcalRecHitsEE_*_*',
    'keep *_PhotonIDProd_*_*',
    'keep *_hybridSuperClusters_*_*',
    'keep *_multi5x5BasicClusters_*_*',
    'keep *_particleFlow_*_*',
    'keep *_particleFlowClusterECAL_*_*',
    'keep *_particleFlowClusterHCAL_*_*',
    'keep *_particleFlowClusterHFEM_*_*',
    'keep *_particleFlowClusterHFHAD_*_*',
    'keep *_particleFlowClusterPS_*_*',
    'keep *_photons_*_*',
    'keep *_photonCore_*_*',
    'keep *_pfElectronTranslator_*_*',
    'keep *_multi5x5SuperClustersWithPreshower_*_*',
    'keep *_multi5x5PreshowerClusterShape_*_*',
    'keep *_correctedHybridSuperClusters_*_*',
    'keep *_correctedMulti5x5SuperClustersWithPreshower_*_*',

    #jets (and PF)
    #use of JPT jets requires running JPT corrections to the AK5 calo jets
    'keep *_ak5CaloJets_*_*',
    'keep *_ak5PFJets_*_*',
    'keep *_ak5JetTracksAssociatorAtVertex_*_*',
    'keep *_ak5JetExtender_*_*',
    'keep *_towerMaker_*_*',
    'keep *_ak5JetID_*_*',
    'keep *_ak5PFJetsRecoTauPiZeros_*_*',
    'keep *_JetPlusTrackZSPCorJetAntiKt5_*_*',
    'keep *_TCTauJetPlusTrackZSPCorJetAntiKt5_*_*',
    'keep *_trackRefsForJets_*_*',

    #cosmics, halo, beamspot, and primary vertices
    'keep *_cosmicsVeto_*_*',
    'keep *_BeamHaloSummary_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_GlobalHaloData_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep *_offlinePrimaryVerticesWithBS_*_*',

    #MET (and PF)
    'keep *_corMetGlobalMuons_*_*',
    'keep *_met_*_*',
    'keep *_tcMet_*_*',
    'keep *_tcMetWithPFclusters_*_*',
    'keep *_pfMet_*_*',

    #cut scan value maps
    #'keep *_CutScanProducer_*_*',

    #category information
    'keep *_EDMCategoryProducer_*_*'
    ),

    splitLevel = cms.untracked.int32(0)
    )

#output
process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string(
    '/data2/yohay/skim_studies/CS_test_events_35GeV-30GeV.root'
    ),
                                  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('skim')
    ),
                                  outputCommands = process.eventContent.outputCommands
                                  )
#debug information
process.options = cms.untracked.PSet(wantSummary  = cms.untracked.bool(True))

#execution path
process.skim = cms.Path(
    process.CutScanProducer*process.EDMCategoryProducer*process.SkimAnalyzer*
    process.TriggerAnalyzer
    )

#output path
process.e = cms.EndPath(process.output)
