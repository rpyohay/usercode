import FWCore.ParameterSet.Config as cms

OrbitNumberAnalyzer = cms.EDAnalyzer('OrbitNumberAnalyzer',
                                     #0 for the calibration orbit, 1 for the next-to-calibration orbit, 2 for the next-to-next-to-
                                     #calibration orbit, etc. up to 111 for the before-calibration orbit
                                     orbitType = cms.uint32(0)
)
