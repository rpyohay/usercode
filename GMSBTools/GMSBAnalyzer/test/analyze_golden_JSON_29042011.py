import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS2")

#desired INFO messages
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('Error')
process.MessageLogger.categories.append('GMSBAnalyzer')
process.MessageLogger.categories.append('GMSBAnalyzerNullPointer')
#process.MessageLogger.categories.append('GMSBAnalyzerEfficiencyCalculation')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

#events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#input
process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_10_1_9d2.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_11_1_f3h.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_12_1_Qgs.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_13_1_g0R.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_14_1_LFc.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_15_1_6tw.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_16_1_HOL.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_17_1_3y1.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_18_1_tg4.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_19_1_1i8.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_1_1_aWQ.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_2_1_Fxr.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_3_1_aGE.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_4_1_B7x.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_5_1_VLW.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_6_1_H1v.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_7_1_SBx.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_8_1_1B2.root',
    'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v1/events_9_1_JGV.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_10_1_GEk.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_11_1_WQs.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_12_2_3zk.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_13_1_T8W.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_14_1_oaL.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_15_1_hHO.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_16_1_EbN.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_17_1_XHZ.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_18_1_ADe.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_19_1_dY4.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_1_1_kyY.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_20_1_o70.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_21_1_AWM.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_22_1_8y5.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_23_1_axI.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_24_1_ZgV.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_25_1_x0f.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_26_1_WBr.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_27_1_iDW.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_28_1_ubX.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_29_1_F7Z.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_2_1_dFg.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_30_1_F61.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_31_1_WwG.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_32_1_I7S.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_33_1_2IZ.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_34_1_neb.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_35_1_uVB.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_36_1_qnn.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_37_1_wyz.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_38_2_VdC.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_39_1_8qg.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_3_1_Hgy.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_40_1_MzL.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_41_1_u1L.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_42_1_AnF.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_43_1_eqp.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_44_1_cwV.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_45_1_CGV.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_46_1_Db9.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_47_1_2MQ.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_48_1_sBx.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_49_1_LVb.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_4_1_J31.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_5_1_wFt.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_6_3_0O8.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_7_1_9Iz.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_8_1_zHT.root',
     'rfio:/castor/cern.ch/user/y/yohay/421/analysis_golden_JSON_29042011_Photon-Run2011A-PromptReco-v2/events_9_1_1rU.root'
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

#trigger analyzer
process.load('TriggerAnalysis/TriggerAnalyzer/triggeranalyzer_cfi')
#process.TriggerAnalyzer.outputFile = cms.untracked.string("/data2/yohay/trigger_analysis_10.root")
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

#category producer
process.load('GMSBTools/Producers/edmcategoryproducer_cfi')
process.EDMCategoryProducer.photonTag = cms.untracked.InputTag("photons", "", "RECO")
process.EDMCategoryProducer.recHitTagEB = cms.untracked.InputTag("reducedEcalRecHitsEB", "",
                                                                 "RECO"
                                                                 )
process.EDMCategoryProducer.recHitTagEE = cms.untracked.InputTag("reducedEcalRecHitsEE", "",
                                                                 "RECO"
                                                                 )
process.EDMCategoryProducer.photonETMin = cms.untracked.double(35.0) #GeV

#remove timing cut due to time reco bias in 2011 prompt reco before May technical stop
#see https://hepex.phys.virginia.edu:8080/rpyohay/2 for details
process.EDMCategoryProducer.photonAbsSeedTimeMax = cms.untracked.double(-1.0) #ns

#GMSB analyzer
process.load('GMSBTools/GMSBAnalyzer/gmsbanalyzer_cfi')
process.GMSBAnalyzer.outputFile = cms.untracked.string("/data2/yohay/analysis_11052011.root")
process.GMSBAnalyzer.photonTag = cms.untracked.InputTag("photons", "", "RECO")
process.GMSBAnalyzer.categoryTag = cms.untracked.InputTag("EDMCategoryProducer", "eventCategory",
                                                          "ANALYSIS2")
process.GMSBAnalyzer.diEMETTag = cms.untracked.InputTag("EDMCategoryProducer", "evtDiEMET",
                                                        "ANALYSIS2")
process.GMSBAnalyzer.photonSeedTimeTag = cms.untracked.InputTag("EDMCategoryProducer",
                                                                "photonSeedTime", "ANALYSIS2")
process.GMSBAnalyzer.photonE2OverE9Tag = cms.untracked.InputTag("EDMCategoryProducer",
                                                                "photonE2OverE9", "ANALYSIS2")
process.GMSBAnalyzer.passAllTag = cms.untracked.InputTag("EDMCategoryProducer", "passingPhotons",
                                                         "ANALYSIS2")
process.GMSBAnalyzer.passETMinTag = cms.untracked.InputTag("EDMCategoryProducer", "passETMin",
                                                           "ANALYSIS2")
process.GMSBAnalyzer.passAbsEtaMaxTag = cms.untracked.InputTag("EDMCategoryProducer",
                                                               "passAbsEtaMax", "ANALYSIS2")
process.GMSBAnalyzer.passECALIsoMaxTag = cms.untracked.InputTag("EDMCategoryProducer",
                                                                "passECALIsoMax", "ANALYSIS2")
process.GMSBAnalyzer.passHCALIsoMaxTag = cms.untracked.InputTag("EDMCategoryProducer",
                                                                "passHCALIsoMax", "ANALYSIS2")
process.GMSBAnalyzer.passHOverEMaxTag = cms.untracked.InputTag("EDMCategoryProducer",
                                                               "passHOverEMax", "ANALYSIS2")
process.GMSBAnalyzer.passTrackIsoMaxTag = cms.untracked.InputTag("EDMCategoryProducer",
                                                                 "passTrackIsoMax", "ANALYSIS2")
process.GMSBAnalyzer.passSigmaIetaIetaMaxTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passSigmaIetaIetaMax", "ANALYSIS2"
    )
process.GMSBAnalyzer.passAbsSeedTimeMaxTag = cms.untracked.InputTag(
    "EDMCategoryProducer", "passAbsSeedTimeMax", "ANALYSIS2"
    )
process.GMSBAnalyzer.passE2OverE9MaxTag = cms.untracked.InputTag("EDMCategoryProducer",
                                                                 "passE2OverE9Max", "ANALYSIS2")
process.GMSBAnalyzer.HLTPathGG = cms.untracked.vstring("HLT_Photon32_CaloIdL_Photon26_CaloIdL")
process.GMSBAnalyzer.HLTPathFF = cms.untracked.vstring("HLT_Photon32_CaloIdL_Photon26_CaloIdL")
process.GMSBAnalyzer.HLTPathAllPhotons = cms.untracked.vstring(
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

    #category information
    'keep *_EDMCategoryProducer_*_*'
    ),

    splitLevel = cms.untracked.int32(0)
    )

#output
process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string('events.root'),
                                  outputCommands = process.eventContent.outputCommands
                                  )
#debug information
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#execution path
#process.p = cms.Path(process.EDMCategoryProducer*process.GMSBAnalyzer*process.TriggerAnalyzer)
process.p = cms.Path(process.EDMCategoryProducer*process.GMSBAnalyzer)

#output path
#process.e = cms.EndPath(process.output)
