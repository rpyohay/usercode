import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

#append CutScanProducer to get info messages
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.categories.append('CutScanProducer')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/data2/yohay/diPhoton_skim_Run2011A_30GeV-30GeV_9_1_WNE.root'
    )
                            )

process.load('Skims/CutScanProducer/cutscanproducer_cfi')
process.CutScanProducer.photonTag = cms.untracked.InputTag("photons", "", "RECO")

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
    'keep *_CutScanProducer_*_*'
    ),

    splitLevel = cms.untracked.int32(0)
    )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(
    '/data2/yohay/diPhoton_skim_Run2011A_30GeV-30GeV_cutScan.root'
    ),
                               outputCommands = process.eventContent.outputCommands
                               )

process.p = cms.Path(process.CutScanProducer)

process.e = cms.EndPath(process.out)
