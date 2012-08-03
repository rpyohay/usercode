import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
##     'file:/data1/yohay/DYToTauTauSkim_latestTauID_nearlyAllClustersDropped.root'
    #'file:/data1/yohay/044C78B0-BD95-E111-973C-0026189438D2.root'
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_1.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_10.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_2.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_3.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_4.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_5.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_6.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_7.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_8.root',
    'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125/Summer12_NMSSMHiggs_9.root'
##     'file:Summer12_NMSSMHiggs_1.root'
    ),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

#process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.GeometryExtended_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load('Configuration.EventContent.EventContent_cff')
#process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

#process.GlobalTag.globaltag = "START52_V9::All"

process.load("BoostedTauAnalysis.TauAnalyzer.tauanalyzer_cfi")

process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('/data1/yohay/test.root')
    )

process.p = cms.Path(process.TauAnalyzer)

#process.end = cms.EndPath(process.output)
