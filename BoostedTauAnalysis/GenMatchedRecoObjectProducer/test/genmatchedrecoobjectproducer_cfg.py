import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

#pseudoscalar a and Z PDG IDs
A_PDGID = 36
Z_PDGID = 23

#tau decay types
TAU_HAD = 0
TAU_MU = 1
TAU_E = 2

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data1/yohay/NMSSMHiggs_gg_skim.root'
    )
    )

#add a collection of tight muons
#(cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId) to the event
process.tightMuonSelector = cms.EDFilter('CustomMuonSelector',
                                         muonTag = cms.InputTag('muons'),
                                         vtxTag = cms.InputTag('offlinePrimaryVertices')
                                         )

#add a collection of tight muons matched to gen muonic tau decays
process.genMuTauMatchedMuonSelector = cms.EDFilter('GenMatchedMuonProducer',
                                                   genParticleTag = cms.InputTag('genParticles'),
                                                   recoObjTag = cms.InputTag('tightMuonSelector'),
                                                   momPDGID = cms.int32(A_PDGID),
                                                   primaryTauDecayType = cms.uint32(TAU_MU),
                                                   sisterTauDecayType = cms.uint32(TAU_HAD),
                                                   chargedHadronPTMin = cms.double(0.0),
                                                   neutralHadronPTMin = cms.double(0.0),
                                                   chargedLeptonPTMin = cms.double(0.0),
                                                   totalPTMin = cms.double(0.0),
                                                   applyPTCuts = cms.bool(False),
                                                   countKShort = cms.bool(True)
                                                   )

#add a collection of tight muons matched to gen muonic tau decays matched to trigger objects
process.triggerMatchedMuonSelector = cms.EDProducer(
    'trgMatchedMuonProducer',
    InputProducer = cms.InputTag('genMuTauMatchedMuonSelector'),
    hltTags = cms.VInputTag(cms.InputTag('HLT_Mu17_Mu8_v16', '', 'HLT')),
    

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('/data1/yohay/debug.root'),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p'))
                               )

  
process.p = cms.Path(process.tightMuonSelector*process.genMuTauMatchedMuonSelector)

process.e = cms.EndPath(process.out)
