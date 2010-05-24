import FWCore.ParameterSet.Config as cms

process = cms.Process("CutAnalyzer")

#keep the logging output to a nice level
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet( reportEvery = cms.untracked.int32(100) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_0.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_1.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_10.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_11.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_12.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_13.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_14.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_15.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_16.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_17.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_18.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_19.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_2.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_20.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_21.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_22.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_23.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_24.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_25.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_26.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_27.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_28.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_29.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_3.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_30.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_31.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_32.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_33.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_34.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_35.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_36.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_37.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_38.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_39.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_4.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_40.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_41.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_42.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_43.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_44.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_45.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_46.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_47.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_48.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_49.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_5.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_50.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_51.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_52.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_53.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_54.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_55.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_56.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_57.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_58.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_59.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_6.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_60.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_61.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_62.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_63.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_64.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_65.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_66.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_67.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_68.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_69.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_7.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_70.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_71.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_72.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_73.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_74.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_75.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_76.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_77.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_78.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_79.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_8.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_80.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_81.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_82.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_83.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_84.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_85.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_86.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_87.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_88.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_89.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_9.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_90.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_91.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_92.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_93.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_94.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_95.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_96.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_97.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_98.root',
    'file:/data/yohay/GGM/tan_beta_10/RECO/GGM_tan_beta_10_RECO_99.root'
    #'file:/data/yohay/GGM/tan_beta_10/RECO/skim/ff_test_GGM_HLT15_pT20.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("/data/yohay/GGM/tan_beta_10/RECO/skim/ff_test_GGM_HLT15_pT20.root"),
                               #fileName = cms.untracked.string("ff_test_photon_jet_pT15_HLT15_pT20.root"),
                               SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring( 'p' ) )
                               )

#load the HLT filter
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT8E29")
process.hltHighLevel.HLTPaths = cms.vstring("HLT_Photon15_L1R")

#load the sample maker
process.load("GMSBTools.SampleMaker.samplemaker_cfi")

#load the sample analyzer
process.load("GMSBTools.CutAnalyzer.cutanalyzer_cfi")

#show TrigReport
#process.options = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True)
    #)

#make the ff sample and check the distributions
process.p = cms.Path(process.hltHighLevel*process.SampleMaker*process.CutAnalyzer)
#process.p = cms.Path(process.CutAnalyzer)
process.e = cms.EndPath(process.out)
