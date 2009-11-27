import FWCore.ParameterSet.Config as cms

TEFilter = cms.EDFilter('TEFilter',
                        unpacker = cms.InputTag("ecalEBunpacker"),
                        eventType = cms.string("LED"),
                        partition1 = cms.uint32(1), #1 - EE+, 2 - EB+, 3 - EB-, 4 - EE-, 0 - none
                        partition2 = cms.uint32(4),
                        partition3 = cms.uint32(0),
                        partition4 = cms.uint32(0)
)
