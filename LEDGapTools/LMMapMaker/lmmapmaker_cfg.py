import FWCore.ParameterSet.Config as cms

process = cms.Process("LMMapMaker")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("EmptySource"
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(0) )

#keep the logging output to a nice level
process.MessageLogger = cms.Service("MessageLogger")

#make the map
process.load("LEDGapTools.LMMapMaker.lmmapmaker_cfi")

process.p = cms.Path(process.LMMapMaker)
