import FWCore.ParameterSet.Config as cms

CalibrationSequenceAnalyzer = cms.EDAnalyzer('CalibrationSequenceAnalyzer',
                                             EEDigiCollection = cms.InputTag("ecalEBunpacker","eeDigis"),
                                             unpacker = cms.InputTag("ecalEBunpacker"),
                                             outFile = cms.string("/afs/cern.ch/user/y/yohay/scratch0/tree_RUN-ped_1st_sample-events_MIN_EVT-MAX_EVT.root"),
                                             alpha = cms.double(1.2812),
                                             tRise = cms.double(2.0862)
)
