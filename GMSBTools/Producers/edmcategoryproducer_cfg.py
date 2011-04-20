import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/data2/yohay/diPhoton_skim_Run2011A_30GeV-30GeV_cutScan.root'
    )
)

process.load('GMSBTools/Producers/edmcategoryproducer_cfi')
process.CutScanProducer.photonTag = cms.untracked.InputTag("photons", "", "OWNPARTICLES")
process.CutScanProducer.recHitTagEB = cms.untracked.InputTag("reducedEcalRecHitsEB", "", "OWNPARTICLES"),
process.CutScanProducer.recHitTagEE = cms.untracked.InputTag("reducedEcalRecHitsEE", "", "OWNPARTICLES"),
                                     
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('/data2/yohay/diPhoton_skim_Run2011A_30GeV-30GeV_categorized.root')
)

  
process.p = cms.Path(process.EDMCategoryProducer)

process.e = cms.EndPath(process.out)
