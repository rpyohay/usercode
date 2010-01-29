import FWCore.ParameterSet.Config as cms

PATVerifier = cms.EDAnalyzer('PATVerifier',
                             outName = cms.string("PATVerifier_test_15-20_GeV.root"),
                             photonSrc = cms.InputTag("cleanLayer1Photons"),
                             # Some extra cuts you might wish to make
                             #  before histograms/TTrees are filled.
                             # Minimum Et
                             #minPhotonEt     = cms.double(10.0),
                             minPhotonEt     = cms.double(0.0),
                             # Minimum and max abs(eta)
                             minPhotonAbsEta = cms.double(0.0),
                             maxPhotonAbsEta = cms.double(3.0),
                             # Minimum R9 = E(3x3) / E(SuperCluster)
                             #minPhotonR9     = cms.double(0.3),
                             minPhotonR9     = cms.double(0.0),
                             # Maximum HCAL / ECAL Energy
                             #maxPhotonHoverE = cms.double(0.2)
                             maxPhotonHoverE = cms.double(1.0)
)
