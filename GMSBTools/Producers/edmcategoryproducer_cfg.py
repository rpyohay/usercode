import FWCore.ParameterSet.Config as cms

process = cms.Process("CATEGORY")

#append EDMCategoryProducer to get info messages
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.categories.append('EDMCategoryProducer')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/data2/yohay/diPhoton_skim_Run2011A_30GeV-30GeV_cutScan.root'
    )
)

#global tag
process.load('Geometry.CaloEventSetup.CaloTopology_cfi')
process.load('Geometry.CaloEventSetup.CaloGeometry_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('GR_P_V14::All')

process.load('GMSBTools/Producers/edmcategoryproducer_cfi')
process.EDMCategoryProducer.photonTag = cms.untracked.InputTag("photons", "", "RECO")
process.EDMCategoryProducer.recHitTagEB = cms.untracked.InputTag("reducedEcalRecHitsEB", "", "RECO")
process.EDMCategoryProducer.recHitTagEE = cms.untracked.InputTag("reducedEcalRecHitsEE", "", "RECO")
                                     
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('/data2/yohay/diPhoton_skim_Run2011A_30GeV-30GeV_categorized.root')
)

  
process.p = cms.Path(process.EDMCategoryProducer)

process.e = cms.EndPath(process.out)
