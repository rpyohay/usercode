import FWCore.ParameterSet.Config as cms

process = cms.Process("PATVerifier")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PATTools.PATVerifier.patverifier_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/y/yohay/PhotonJet_Pt15to20-Summer09-MC_31X_V9_7TeV-v1-GEN-SIM-RECO/photonjet_1.root',
    'rfio:/castor/cern.ch/user/y/yohay/PhotonJet_Pt15to20-Summer09-MC_31X_V9_7TeV-v1-GEN-SIM-RECO/photonjet_2.root'
    )
)

process.p = cms.Path(process.PATVerifier)
