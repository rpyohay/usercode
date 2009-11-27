import FWCore.ParameterSet.Config as cms

LeakageAnalyzer = cms.EDAnalyzer('LeakageAnalyzer',
                                 EEDigiCollection = cms.InputTag("ecalEBunpacker","eeDigis"),
                                 outName = cms.string("/afs/cern.ch/user/y/yohay/scratch0/leakage_121620_random-triggered_filtered.root")
)
