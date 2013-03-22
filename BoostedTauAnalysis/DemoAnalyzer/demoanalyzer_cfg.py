import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/data1/yohay/NMSSMHiggs_WH_files251-500_24Sep12.root'
    )
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root") )

process.demo = cms.EDAnalyzer('DemoAnalyzer',
#                              pdgid = cms.int32(15)
)


process.p = cms.Path(process.demo)
