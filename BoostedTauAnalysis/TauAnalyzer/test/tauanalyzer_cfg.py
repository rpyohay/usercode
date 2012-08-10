import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/data1/yohay/NMSSMHiggs_HaaOnly_skim_v4.root'
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_1.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_10.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_2.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_3.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_4.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_5.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_6.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_7.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_8.root',
##     'root://eoscms//eos/cms/store/user/yohay/NMSSM_Higgs_a9_H1115_H2125_H3500/Summer12_NMSSMHiggs_9.root'
##     'file:/data1/yohay/FCDE987D-859B-E111-B445-0025B3E05DDA.root'
    ),
##                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.load("BoostedTauAnalysis.TauAnalyzer.tauanalyzer_cfi")

process.p = cms.Path(process.TauAnalyzer)
