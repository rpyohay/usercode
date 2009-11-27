import FWCore.ParameterSet.Config as cms

process = cms.Process("TEFilter")

process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
    #110397
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/FCB74D03-5584-DE11-A2CB-001D09F290BF.root', #639
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/F6C8D2D1-4B84-DE11-B1B8-001D09F2438A.root', #614
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/F6254404-5584-DE11-8B6D-001D09F2503C.root', #614
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/F4F165E5-4684-DE11-A875-000423D94494.root', #628
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/F40506E6-4684-DE11-A6AA-000423D6006E.root', #625
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/F2BC58D7-5784-DE11-BE0D-001D09F295FB.root', #612
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/EEE3853C-5284-DE11-9C68-000423D944F8.root', #616
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/EE459739-4684-DE11-BDEF-001D09F24EE3.root' #601
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/E0AA9433-5284-DE11-AE9F-000423D6CA42.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/DE4ED638-4684-DE11-B6D4-001D09F23A02.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/DCB515AC-4284-DE11-83C6-001D09F24D67.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/DA291311-5084-DE11-B133-001D09F2932B.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/CCBBC79A-4784-DE11-B0DA-000423D94494.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/CC932BAF-4284-DE11-A052-001D09F253FC.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/CA159E13-5084-DE11-9EBF-001D09F2527B.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/C6C708B9-5584-DE11-ADD5-001D09F23A84.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/C694A204-4E84-DE11-B0FD-000423D996B4.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/BEFC175B-4F84-DE11-A7B3-000423D6B2D8.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/BCBF6BBA-4984-DE11-A5ED-000423D985B0.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/B8773E02-4984-DE11-A60D-000423D99160.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/B45501AD-4284-DE11-8D11-001D09F23A3E.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/AED7A147-5484-DE11-81F3-000423D9A212.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/AC25EFDB-3F84-DE11-8BE7-001D09F2546F.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/AAA908F8-4184-DE11-BB97-001D09F24D67.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/A26F0535-5984-DE11-8CBE-0019B9F71A6B.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/98C21635-5984-DE11-B531-001D09F242EA.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/98319A28-4B84-DE11-A6EB-000423D98A44.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/90EC8DFC-5284-DE11-AAD3-000423D98BE8.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/9010A69C-5384-DE11-86D7-001D09F24493.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/82DA1225-5784-DE11-8458-000423DD2F34.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/80BB5424-4B84-DE11-A2A4-000423D94534.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/7AB67453-4884-DE11-944F-001D09F25401.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/70CF448E-4C84-DE11-97DA-000423D99160.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/68217147-5484-DE11-BCA3-001D09F2A465.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/6806AB3E-4D84-DE11-AEEB-001D09F23A3E.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/60791574-5684-DE11-B28F-000423D99896.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/583BB56C-4A84-DE11-9579-001D09F28EC1.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/508C45AF-4284-DE11-8936-001D09F24F65.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/4A7750C0-5084-DE11-A1F2-001D09F23A34.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/48A2BC2E-5284-DE11-8441-000423D9517C.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/46C41A36-5984-DE11-ACC7-001D09F2910A.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/469BAB8C-4584-DE11-ACE5-000423D6B2D8.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/46634E89-5184-DE11-BB12-000423D6CA42.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/40807E01-7484-DE11-8AB7-001D09F2906A.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/36FFEDCE-4484-DE11-B4B5-001D09F2527B.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/30680F33-5784-DE11-9FEA-000423D9870C.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/2EE60EB9-5584-DE11-ADF5-0019B9F71A6B.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/2EABA0D7-3F84-DE11-AB75-001D09F24DDF.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/2CADCF9B-4784-DE11-A539-000423D6006E.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/2C3127F1-4D84-DE11-9040-000423D9A212.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/28828A43-4184-DE11-AA34-001D09F295A1.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/26C0D301-4984-DE11-9EB9-000423D99264.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/225A878C-4C84-DE11-BBBD-000423D98834.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/1EC1378F-4084-DE11-B53D-001D09F24DDF.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/1C6FC17A-5184-DE11-9E47-000423D9A212.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/1C5DC2DE-4484-DE11-8699-001D09F295FB.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/124E4818-4484-DE11-9199-000423D98834.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/0E352BAC-4E84-DE11-A96C-000423D985B0.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/0C0D7F6D-4A84-DE11-835A-001D09F28F25.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/0AAFFFD5-5784-DE11-B1AA-001D09F242EA.root',
#    '/store/data/CRAFT09/TestEnables/RAW/v1/000/110/397/0677DF43-4184-DE11-A142-000423D951D4.root'

    #116815
    #'/store/data/Commissioning09/TestEnables/RAW/v3/000/116/815/FAB76AAB-DEB4-DE11-9BBF-000423D98634.root',
    #'/store/data/Commissioning09/TestEnables/RAW/v3/000/116/815/F6559C87-DCB4-DE11-9613-001D09F25401.root'

    #117951 (Fri. Oct. 23, ~7-8 AM)
    #DQM implies that EE+ some EE+ events are out of time, and some are no more than pedestal (varies by sector)
    #one sector appears to have no recorded data, not even pedestal
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/117/951/ECDE0338-9FBF-DE11-92BA-0016177CA7A0.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/117/951/E8E874A8-94BF-DE11-ABA4-001D09F24353.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/117/951/C4048E0E-96BF-DE11-81CF-000423D6AF24.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/117/951/B25195DE-97BF-DE11-A760-001D09F2A690.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/117/951/685F2646-9ABF-DE11-9475-001D09F242EA.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/117/951/1C0FA1BB-9BBF-DE11-8519-003048D2BF1C.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/117/951/00AA92DF-98BF-DE11-AD46-001D09F28D4A.root'

    #118346
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/FCFF0033-78C2-DE11-A10E-001D09F2546F.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/F4EC66A5-79C2-DE11-B404-001D09F29146.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/F445B23E-73C2-DE11-88B5-000423D94A20.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/EE68ED7A-70C2-DE11-BD0D-001617DBCF6A.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/EA2C07AC-74C2-DE11-ABF9-001D09F2B30B.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/E87CD04D-7AC2-DE11-8AAC-000423D6C8EE.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/DAD09023-71C2-DE11-8DCB-001617E30CD4.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/D2431DD6-71C2-DE11-BC85-001617E30CD4.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/C0831E11-76C2-DE11-8A48-001D09F25041.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/B2CA7104-6FC2-DE11-81B9-000423D98EC8.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/AA7CC734-6CC2-DE11-BFF8-001D09F2910A.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/8486FCBD-6FC2-DE11-A482-000423D986C4.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/7AC2A219-82C2-DE11-AFD5-001D09F250AF.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/7813137E-77C2-DE11-B553-001D09F24D8A.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/6EA07E06-6FC2-DE11-B057-000423D98750.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/54EDFFCC-76C2-DE11-AF44-0030487A18A4.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/54D653CC-76C2-DE11-AF0E-0030487D0D3A.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/524A024E-7AC2-DE11-A2CC-000423D98B28.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/50031330-6CC2-DE11-9AE6-001D09F28EA3.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/48A12A16-74C2-DE11-8462-001617DBCF6A.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/400919A1-72C2-DE11-9696-001D09F252F3.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/3E4C3A9E-6DC2-DE11-9C89-0030487D0D3A.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/2A18929D-6DC2-DE11-B570-000423D6C8EE.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/1AC915A1-6DC2-DE11-884E-000423D98750.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/1852EBE3-78C2-DE11-99AB-001D09F2423B.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/0EC54961-75C2-DE11-85E0-001D09F2545B.root',
#    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/346/06D51A74-70C2-DE11-8A2A-001617C3B6DC.root'

    #118333
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/EC540605-61C2-DE11-9310-0030487C5CFA.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/E8F74111-5CC2-DE11-A1DD-001D09F29597.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/B4603AA7-61C2-DE11-B3F3-001617E30D52.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/AE70E76C-5DC2-DE11-8B46-000423D94A20.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/9093EB7A-64C2-DE11-ACE8-001617C3B6DC.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/8A251C5E-60C2-DE11-9878-001617E30E28.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/84CAC287-5FC2-DE11-A6C7-000423D94E70.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/6AF7095D-62C2-DE11-8504-0019DB29C5FC.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/6A5515C0-63C2-DE11-9B6B-000423D99BF2.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/46CAC0C3-63C2-DE11-93F2-001D09F254CE.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/46C785A9-61C2-DE11-B2FB-000423D99658.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/42BDC0D4-5EC2-DE11-85CA-001617E30CD4.root',
    '/store/data/Commissioning09/TestEnables/RAW/v3/000/118/333/3C956B12-63C2-DE11-B23B-001617E30D12.root'
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
                               fileName = cms.untracked.string('LED_digis_118333.root'),
                               SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
                               #outputCommands = cms.untracked.vstring('drop *',
                                                                       #"keep *_*_eeDigiSkim_*")
                               )

#run the unpacker and skim off the LED events
process.p = cms.Path(process.ecalEBunpacker*process.TEFilter)

#write the resulting collections to the output file
process.e = cms.EndPath(process.out)
