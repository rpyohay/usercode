import FWCore.ParameterSet.Config as cms

process = cms.Process("CalibrationSequenceAnalyzer")

#calo geometry service model
#process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
#process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

#ECAL mapping
#process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
#process.load("Geometry.EcalMapping.EcalMapping_cfi")
#process.eegeom = cms.ESSource("EmptyESSource",
    #recordName = cms.string('EcalMappingRcd'),
    #iovIsRunNotTime = cms.bool(True),
    #firstValid = cms.vuint32(1)
#)

#database stuff
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR09_E_V3::All'

#keep the logging output to a nice level
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAX_EVTS) )

#for large input files
#process.AdaptorConfig = cms.Service("AdaptorConfig")

process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32(NUM_SKIPPED_EVTS),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/y/yohay/LED_digis_RUN.root'
    )
)

#run the analysis
process.load("prompt_feedback.CalibrationSequenceAnalyzer.calibrationsequenceanalyzer_cfi_RUN_MIN_EVT-MAX_EVT")

process.p = cms.Path(process.CalibrationSequenceAnalyzer)
