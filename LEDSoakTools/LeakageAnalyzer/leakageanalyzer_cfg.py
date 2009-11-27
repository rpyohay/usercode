import FWCore.ParameterSet.Config as cms

process = cms.Process("LeakageAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_1.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_10001-15000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_15001-20000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_20001-25000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_25001-30000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_30001-35000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_35001-40000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_40001-45000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_45001-50000.root', broken
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_50001-55000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_55001-60000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_60001-65000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_65001-70000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_70001-75000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_75001-80000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_80001-85000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_85001-90000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_90001-95000.root',
        #'rfio:/castor/cern.ch/user/x/xquan/digis_121620_95001-100000.root',
    #'rfio:/castor/cern.ch/user/y/yohay/calo-triggered_ECAL_digis_121620_70000.root'
    #'file:/afs/cern.ch/user/y/yohay/scratch0/calo-triggered_ECAL_digis_121620_test2.root'

    #121620 randoms filtered by BX
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_1-5000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_10001-15000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_15001-20000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_20001-25000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_25001-30000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_30001-35000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_35001-40000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_40001-45000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_45001-50000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_50001-55000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_5001-10000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_55001-60000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_60001-65000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_65001-70000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_70001-75000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_75001-80000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_80001-85000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_85001-90000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_90001-95000.root',
    'rfio:/castor/cern.ch/user/x/xquan/digis_filtered_121620_95001-100000.root'
    )
)

process.load("L1TriggerConfig.L1GtConfigProducers.L1GtConfig_cff")
process.load("LEDSoakTools.LeakageAnalyzer.leakageanalyzer_cfi")

process.p = cms.Path(process.LeakageAnalyzer)
