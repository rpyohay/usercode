import FWCore.ParameterSet.Config as cms

## process = cms.Process("SKIM")
process = cms.Process("HADTAUSKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.GlobalTag.globaltag = "START52_V9::All"

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
##     "file:/data1/yohay/044C78B0-BD95-E111-973C-0026189438D2.root"
##     'file:/data1/yohay/DYToTauTauSkim_latestTauID.root'
##     'file:/data1/yohay/DYToTauTauSkim_latestTauID_someClustersDropped.root'
    'file:/data1/yohay/DYToTauTauSkim_latestTauID_nearlyAllClustersDropped.root'
    )
    )

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

process.HadronicTauDecayFinder = cms.EDFilter('HadronicTauDecayFinder',
                                              genParticleTag = cms.InputTag('genParticles'),
                                              momPDGID = cms.uint32(23)
                                              )

process.output = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring(
    "keep *",
    "drop *_*ak7*_*_*",
    "drop *_*kt4*_*_*",
    "drop *_kt6GenJets_*_*",
    "drop *_kt6CaloJets*_*_*",
    "drop *_kt6PFJetsCentral*_*_*",
    "drop *_kt6PFJets_sigma*_*",
    "drop *_kt6JetID_*_*",
    "drop *_kt6PFJets_rhos_*",
    "drop recoPFJets_kt6PFJets_*_*",
    "drop *_fixedGridRho*_*_*",
    "drop *_hfEMClusters_*_*",
    "drop *_tevMuons_*_*",
    "drop *_eid*_*_*",
    "drop *_muonMETValueMapProducer_muCorrData_*",
    "drop *_muons_muonShowerInformation_*",
    "drop *_muons_combined_*",
    "drop *_muons_csc_*",
    "drop *_muons_dt_*",
    "drop l1extraL1HFRingss_*_*_*",
    "drop *_muonsFromCosmics*_*_*",
    "drop recoCaloClusters_*_*_*",
    "drop recoPreshowerCluster*_*_*_*",
    "drop *_hfRecoEcalCandidate_*_*",
    "drop *_regionalCosmicTracks_*_*",
    "drop *_generalV0Candidates_*_*",
    "drop *_selectDigi_*_*",
    "drop *_*BJetTags*_*_*",
    "drop *_castorreco_*_*",
    "drop *_reduced*RecHits*_*_*",
    "drop *_PhotonIDProd_*_*",
    "drop *_*_*photons*_*",
    "drop *_dedx*_*_*",
    "drop *_*_cosmicsVeto_*",
    "drop *_muonTCMETValueMapProducer_*_*",
    "drop *_BeamHaloSummary_*_*",
    "drop *_caloRecoTau*_*_*",
    "drop *_GlobalHaloData_*_*",
    "drop *_hpsTancTau*_*_*",
    "drop *_shrinkingConePFTau*_*_*",
    "drop *_ak5CaloJets_*_*",
    "drop *_ak5TrackJets_*_*",
    "drop *_*_uncleanOnly*_*",
    "drop recoCaloMETs_*_*_*",
    "drop recoConversions_*_*_*",
    "drop *_CastorTowerReco_*_*",
    "drop *_gsfElectron*_*_*",
    "drop *_uncleanedOnlyGsfElectron*_*_*",
    "drop recoJPTJets_*_*_*",
    "drop recoMETs_*_*_*",
    "drop *_photons_*_*",
    "drop *_photonCore_*_*",
    "drop *_*osmicMuons*_*_*",
    "drop *_uncleanedOnly*_*_*",
    "drop *_TriggerResults_*_RECO",
    "drop *_ak5PFJetsRecoTauPiZeros_*_RECO",
    "drop *_hpsPFTauDiscriminationByDecayModeFinding_*_RECO",
    "drop *_hpsPFTauDiscriminationByLooseChargedIsolation_*_RECO",
    "drop *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr_*_*",
    "drop *_hpsPFTauDiscriminationByLooseElectronRejection_*_*",
    "drop *_hpsPFTauDiscriminationByLooseIsolation_*_*",
    "drop *_hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr_*_*",
    "drop *_hpsPFTauDiscriminationByLooseMuonRejection_*_RECO",
    "drop *_hpsPFTauDiscriminationByMVAElectronRejection_*_RECO",
    "drop *_hpsPFTauDiscriminationByMediumChargedIsolation_*_RECO",
    "drop *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr_*_*",
    "drop *_hpsPFTauDiscriminationByMediumElectronRejection_*_*",
    "drop *_hpsPFTauDiscriminationByMediumIsolation_*_*",
    "drop *_hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr_*_*",
    "drop *_hpsPFTauDiscriminationByMediumMuonRejection_*_RECO",
    "drop *_hpsPFTauDiscriminationByRawChargedIsolationDBSumPtCorr_*_RECO",
    "drop *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr_*_RECO",
    "drop *_hpsPFTauDiscriminationByRawGammaIsolationDBSumPtCorr_*_RECO",
    "drop *_hpsPFTauDiscriminationByTightChargedIsolation_*_RECO",
    "drop *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr_*_*",
    "drop *_hpsPFTauDiscriminationByTightElectronRejection_*_*",
    "drop *_hpsPFTauDiscriminationByTightIsolation_*_*",
    "drop *_hpsPFTauDiscriminationByTightIsolationDBSumPtCorr_*_*",
    "drop *_hpsPFTauDiscriminationByTightMuonRejection_*_RECO",
    "drop *_hpsPFTauDiscriminationByVLooseChargedIsolation_*_RECO",
    "drop *_hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr_*_*",
    "drop *_hpsPFTauDiscriminationByVLooseIsolation_*_*",
    "drop *_hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr_*_*",
    "drop *_hpsPFTauProducer_*_RECO",
    "drop *_recoTauAK5PFJets08Region_*_SKIM",
    "drop *_ak5PFJetTracksAssociatorAtVertex_*_SKIM",
    "drop *_ak5PFJetsLegacyHPSPiZeros_*_SKIM",
    "drop *_combinatoricRecoTausDiscriminationByLeadingPionPtCut_*_SKIM",
    "drop *_combinatoricRecoTausHPSSelector_*_SKIM",
    "drop *_hpsSelectionDiscriminator_*_SKIM",
    "drop *_combinatoricRecoTaus_*_SKIM",
    "drop *_hpsPFTauProducerSansRefs_*_SKIM",
    "drop *_pfRecoTauTagInfoProducer_*_SKIM",
    "drop *_recoTauPileUpVertices_*_SKIM",
    "drop *_towerMaker_*_*", #HCAL towers
    "drop *_particleFlowClusterHF*_*_*",
    "drop *_pfPhotonTranslator_*_*",
    "drop *_pfElectronTranslator_*_*",
    "drop *_correctedHybridSuperClusters_*_*",
    "drop *_correctedMulti5x5SuperClustersWithPreshower_*_*",
    ),
    fileName = cms.untracked.string(
    '/data1/yohay/DYToTauTauSkim_latestTauID_fastHadTauFinder.root'
    )
    )

## process.p = cms.Path(process.PFTau*process.HadronicTauDecayFinder)
process.p = cms.Path(process.HadronicTauDecayFinder)

process.end = cms.EndPath(process.output)