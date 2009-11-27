import FWCore.ParameterSet.Config as cms

LeakageFilter = cms.EDFilter('LeakageFilter',
                             unpacker = cms.InputTag("ecalEBunpacker")
)
