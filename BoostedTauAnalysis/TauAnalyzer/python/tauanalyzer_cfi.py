import FWCore.ParameterSet.Config as cms

TauAnalyzer = cms.EDAnalyzer(
    'TauAnalyzer',
    outFileName = cms.string(
    "/data1/yohay/DYToTauTau_analysis_NMSSMHiggs_a9_H1115_H2125_H3500.root"
##     "/data1/yohay/DYToTauTau_analysis_pT.root"
    ),
    genParticleTag = cms.InputTag("genParticles"),
    tauTag = cms.InputTag("hpsPFTauProducer"),
    HPSDiscriminatorTags = cms.VInputTag(
    cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
    cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA")
    ),
##     momPDGID = cms.uint32(23), #Z for Drell-Yan
    momPDGID = cms.uint32(36) #a for signal
    )
