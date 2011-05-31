import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIGGERANALYZER")

process.load("FWCore.MessageService.MessageLogger_cfi")

#stuff you need for the trigger for some unknown reason
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')

#global tag: needed for trigger stuff
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('GR_P_V14::All')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/data2/yohay/diPhoton_skim_Run2011A_30GeV-30GeV_prelimAnalysis.root'
    )
                            )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load('TriggerAnalysis/TriggerAnalyzer/triggeranalyzer_cfi')
process.TriggerAnalyzer.outputFile = cms.untracked.string(
    '/data2/yohay/diPhoton_skim_Run2011A_30GeV-30GeV_trgAnalysis.root'
    )
process.TriggerAnalyzer.unprescaledHLTPaths = cms.untracked.vstring(
    'HLT_Photon32_CaloIdL_Photon26_CaloIdL'
    )

process.p = cms.Path(process.TriggerAnalyzer)
