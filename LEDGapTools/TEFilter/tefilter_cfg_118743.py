import FWCore.ParameterSet.Config as cms

process = cms.Process("TEFilter")

process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
    #118743
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/FCF224C7-F8C3-DE11-9269-001D09F25438.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/F07F8F8F-F4C3-DE11-A5FC-001617C3B6C6.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/E6D100F3-F5C3-DE11-90EA-001D09F252E9.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/E4B1537D-F9C3-DE11-97B8-0019DB29C5FC.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/DC1332A7-F6C3-DE11-9627-001D09F24399.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/C88C8F78-F9C3-DE11-817E-001D09F24024.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/BE7A994B-FCC3-DE11-BADF-001D09F2438A.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/B87D4E28-FFC3-DE11-A356-001617DBD230.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/B01D4035-FAC3-DE11-9486-001D09F2AD4D.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/AE5ACD14-F8C3-DE11-9452-001D09F241B9.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/9AF9A250-FCC3-DE11-8290-001D09F25393.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/84209643-00C4-DE11-BE6F-001D09F23944.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/7EEDB6B3-FDC3-DE11-A837-000423D999CA.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/7443F36B-FEC3-DE11-874F-000423D8F63C.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/724BD760-F7C3-DE11-96A8-001D09F29619.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/70DA151C-F8C3-DE11-B394-001D09F25456.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/5A9CD5AC-FBC3-DE11-8210-000423D98EA8.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/58C5D5A0-F4C3-DE11-8CD5-000423D99B3E.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/2C7CE1A8-F6C3-DE11-8F68-0019B9F72CE5.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/2889D840-F5C3-DE11-BCF6-000423D99B3E.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/1C2175C1-FBC3-DE11-A099-000423D992A4.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/1068E26B-FEC3-DE11-AA86-000423D987E0.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/743/066EECF5-FAC3-DE11-B3F2-001617DBCF6A.root'
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
                               fileName = cms.untracked.string('LED_digis_118743.root'),
                               SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
                               #outputCommands = cms.untracked.vstring('drop *',
                                                                       #"keep *_*_eeDigiSkim_*")
                               )

#run the unpacker and skim off the LED events
process.p = cms.Path(process.ecalEBunpacker*process.TEFilter)

#write the resulting collections to the output file
process.e = cms.EndPath(process.out)
