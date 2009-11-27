import FWCore.ParameterSet.Config as cms

process = cms.Process("TEFilter")

process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
    #118358
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/358/F639FC38-7FC2-DE11-AA85-0019B9F70468.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/358/F044B73C-7FC2-DE11-960F-001D09F24D4E.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/358/EEB24FA5-80C2-DE11-ACAC-0019B9F709A4.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/358/E03B78F3-7FC2-DE11-AB43-001D09F2924F.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/358/D2EB3888-7EC2-DE11-A73C-000423D9997E.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#keep the logging output to a nice level
process.MessageLogger = cms.Service("MessageLogger")

#load all the unpacker stuff
process.load("CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi")
process.load("EventFilter.EcalRawToDigiDev.EcalUnpackerMapping_cfi")
process.load("EventFilter.EcalRawToDigiDev.EcalUnpackerData_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")

#do the filtering
process.load("LEDGapTools.TEFilter.tefilter_cfi")

#create a file to which to write the EE LED events
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('pedestal_digis_118358.root'),
                               SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
                               #outputCommands = cms.untracked.vstring('drop *',
                                                                       #"keep *_*_eeDigiSkim_*")
                               )

#run the unpacker and skim off the LED events
process.p = cms.Path(process.ecalEBunpacker*process.TEFilter)

#write the resulting collections to the output file
process.e = cms.EndPath(process.out)
