import FWCore.ParameterSet.Config as cms

TauAnalyzer = cms.EDAnalyzer(
    'TauAnalyzer',
    outFileName = cms.string(
    '/data1/yohay/NMSSMHiggs_HaaOnly_analysis_newCodeTest.root'
##     '/data1/yohay/NMSSMHiggs_analysis.root'
##     '/data1/yohay/DY_analysis.root'
    ),
    genParticleTag = cms.InputTag("genParticles"),
    tauTag = cms.InputTag("hpsPFTauProducer"),
    muonTag = cms.InputTag("muons"),
    vtxTag = cms.InputTag("offlinePrimaryVertices"),
    HPSDiscriminatorTags = cms.VInputTag(
    cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
    cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA")
    ),
##     momPDGID = cms.uint32(23), #Z for Drell-Yan
    momPDGID = cms.uint32(36), #a for signal
    genMuTauPTMin = cms.double(0.0), #GeV
    genMuPTMin = cms.double(20.0), #GeV
    effVsEtaPTMin = cms.double(10.0) #GeV
    )
