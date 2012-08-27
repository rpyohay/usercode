import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data1/yohay/Summer12_DYToTauTau_skim.root'
    )
)

process.load("BoostedTauAnalysis.TauAnalyzer.tauanalyzer_cfi")
process.TauAnalyzer.outFileName = cms.string(
    '/data1/yohay/Summer12_DYToTauTau_analysis_genMuPTGt20GeV_genTauPTGt15GeV.root'
    )
process.TauAnalyzer.momPDGID = cms.int32(23)
process.TauAnalyzer.genMuPTMin = cms.double(20.0)
process.TauAnalyzer.effVsEtaPTMin = cms.double(15.0)

process.p = cms.Path(process.TauAnalyzer)
