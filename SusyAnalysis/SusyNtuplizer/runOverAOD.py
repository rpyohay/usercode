import FWCore.ParameterSet.Config as cms

# change this to 0 if you run on MC files
realData = 1

process = cms.Process("RA3")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring( )
                            )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections/Configuration/JetCorrectionServices_cff")
process.load("JetMETCorrections/Type1MET/caloMETCorrections_cff")
process.metAnalysisSequence = cms.Sequence(process.producePFMETCorrections*
                                           process.produceCaloMETCorrections)

# calculate rho & jetArea
#process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.kt6PFJets.doRhoFastjet = True
#process.ak5PFJetsWithFastjetArea = process.ak5PFJets.clone()
#process.ak5PFJetsWithFastjetArea.doAreaFastjet = True

# SusyNtuplizer options
process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
process.susyNtuplizer.debugLevel = cms.int32(0)
#process.susyNtuplizer.outputFileName = cms.string(
#    "/data2/yohay/Summer11_TT_debug/Type-I_MET_debug_ntuple.root"
#    )

if realData:
    process.source.fileNames = cms.untracked.vstring(
        'file:/data2/yohay/Photon_Aug5ReReco.root'
        )
    process.GlobalTag.globaltag = 'GR_R_42_V19::All'
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3Residual")
    # JEC for data
    process.jet = cms.Sequence(
#        process.kt6PFJets * process.ak5PFJetsWithFastjetArea *
        # CaloJets
        process.ak5CaloJetsL2L3Residual * process.ak5CaloJetsL1L2L3Residual *
        # PFJets
        process.ak5PFJetsL2L3Residual * process.ak5PFJetsL1FastL2L3Residual
    )
else:
    process.source.fileNames = cms.untracked.vstring(
        'file:/data2/yohay/Summer11_TT_debug/26F2AB20-348C-E011-BCF1-0017A4770400.root'
        )
    process.GlobalTag.globaltag = 'START42_V13::All'
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3")
    # JEC for MC
    process.jet = cms.Sequence(
#        process.kt6PFJets * process.ak5PFJetsWithFastjetArea *
        # CaloJets
        process.ak5CaloJetsL2L3 * process.ak5CaloJetsL1L2L3 *
        # PFJets
        process.ak5PFJetsL2L3 * process.ak5PFJetsL1FastL2L3
    )

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string(
#    "/data2/yohay/Summer11_TT_debug/Type-I_MET_debug_AOD.root"
#    )
#                               )

process.p = cms.Path( process.metAnalysisSequence * process.jet * process.susyNtuplizer )
#process.e = cms.EndPath(process.out)
