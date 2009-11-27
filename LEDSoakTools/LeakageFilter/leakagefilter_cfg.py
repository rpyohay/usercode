import FWCore.ParameterSet.Config as cms

process = cms.Process("LeakageFilter")

process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
    #121620
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/FECAC37A-D9D4-DE11-8010-001D09F291D2.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/FE35F1B7-D6D4-DE11-8CB7-001D09F28F0C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/FCFC25D5-CED4-DE11-B18F-001D09F291D7.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/FC280951-CBD4-DE11-AAA7-000423D98B6C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/F88C5DCB-D3D4-DE11-87ED-001D09F23944.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/F68BD38C-E0D4-DE11-ACF8-001D09F253FC.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/F604ED34-C9D4-DE11-81AE-001D09F2A49C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/F4B05172-CDD4-DE11-97D4-0030487A1990.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/F2DCDA8D-E0D4-DE11-8427-001617C3B706.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/F0E42AD6-CED4-DE11-8C24-000423D985E4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/F030D353-C6D4-DE11-87B5-001D09F2441B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/F019C0E6-DAD4-DE11-A45E-001D09F24353.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/EEE87097-CAD4-DE11-B36E-00304879FBB2.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/EED72803-CCD4-DE11-9C4A-000423D9890C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/EEC7E0B9-CCD4-DE11-83EA-0030487A18A4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/ECCCCB12-C7D4-DE11-9420-000423D98868.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/EC27D0B7-C7D4-DE11-BED1-0030487A1FEC.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/EA877D35-DAD4-DE11-94FB-001D09F24600.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/EA7E8B53-C6D4-DE11-875C-001D09F232B9.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/EA77FB88-C8D4-DE11-8C21-001D09F28F0C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/E8341241-D5D4-DE11-9DAB-001617C3B706.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/E6F41BC1-D8D4-DE11-BF95-000423D33970.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/E6ADF527-DFD4-DE11-A87A-001D09F34488.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/E4E4F8AA-E2D4-DE11-A0AF-001D09F2514F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/E2F3C853-CBD4-DE11-A2E0-001D09F29533.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/E2A4A477-E8D4-DE11-87A7-0019B9F581C9.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/E28D5AB8-CCD4-DE11-B0FB-0030487A3232.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/E0BF2253-C6D4-DE11-8D2F-0019B9F72BAA.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/DE3D4873-D7D4-DE11-87AA-001D09F2841C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/DCB26371-E1D4-DE11-9CAB-001D09F2905B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/DC67D733-C9D4-DE11-BE84-001D09F251B8.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/DC19D458-D7D4-DE11-9C23-0019B9F707D8.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D85B78E8-C9D4-DE11-92BB-003048D2C020.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D81DE8E6-DAD4-DE11-980D-001D09F28D4A.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D81166E8-DAD4-DE11-BDDA-0019B9F7312C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D4FB9A49-E1D4-DE11-8D71-001D09F28F11.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D4CBD329-DFD4-DE11-8B30-001D09F2910A.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D48B0F17-D8D4-DE11-812E-003048D2BF1C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D4673C87-C8D4-DE11-B4BB-001D09F291D2.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D2EA9138-E6D4-DE11-BDC1-000423D951D4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D2BD33F7-D0D4-DE11-A31F-003048D2C020.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/D09666F7-D0D4-DE11-BD53-001D09F2441B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/CCBABC7A-C8D4-DE11-ADA4-001D09F2423B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/CC286D80-E5D4-DE11-92C4-000423D99CEE.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/CA72E585-DED4-DE11-B384-001D09F252F3.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/C6981FD5-CED4-DE11-8691-0019B9F730D2.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/C47707AE-D1D4-DE11-BAFD-001D09F29849.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/C219ED08-DDD4-DE11-8862-001D09F2423B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/C0D6B942-D5D4-DE11-86FF-000423D985B0.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/C047FE88-D4D4-DE11-BFB6-001D09F29321.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/BEA6DCB9-CCD4-DE11-852F-00304879FA4A.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/BEA5555A-DCD4-DE11-B4B0-001D09F2960F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/BE78ADB0-E2D4-DE11-A1A9-001D09F2841C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/BC24068B-CED4-DE11-A484-0019B9F7312C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/B853AED0-DDD4-DE11-AE44-001D09F23D1D.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/B815C400-CCD4-DE11-ABBB-003048D2BED6.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/B4EC70B0-E2D4-DE11-B61F-001D09F2462D.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/B4C5F089-CED4-DE11-9024-0019B9F72BFF.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/B284188C-CED4-DE11-98B4-001D09F23944.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/AEE46B44-D0D4-DE11-9752-000423D9517C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/AE8C39D8-DFD4-DE11-8440-001D09F241F0.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/AAFE8283-C8D4-DE11-8B4F-000423D986A8.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/AABF9F18-D8D4-DE11-A0DB-0019B9F704D6.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/A8A4ED11-C7D4-DE11-8851-000423DD2F34.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/A849FCED-E6D4-DE11-9367-001D09F2915A.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/A818B450-CBD4-DE11-AEFC-001D09F2545B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/A2C6B1EB-D5D4-DE11-97EB-001D09F2924F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/9EE6D04F-DCD4-DE11-83E3-001D09F2B30B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/9EE57B99-CAD4-DE11-B2F5-000423D6A6F4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/9E6FB845-D5D4-DE11-B25E-001617C3B6E2.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/9A674778-D9D4-DE11-9959-001D09F2437B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/96DFD868-E3D4-DE11-A440-001617E30F48.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/946CAF00-CCD4-DE11-8B4E-003048D2BF1C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/92FF1D1D-D3D4-DE11-8FE0-001D09F251D1.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/929DFF3E-E6D4-DE11-A5C1-000423D98E6C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/924A78ED-D5D4-DE11-A541-001D09F2447F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/9013921E-D3D4-DE11-B661-001D09F25208.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/8C6F08FB-E1D4-DE11-AB4C-001D09F23F2A.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/8C62BF4C-DCD4-DE11-B7B1-0019B9F72D71.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/8C438A09-DDD4-DE11-B5FB-001D09F291D7.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/8A2EC0D8-DFD4-DE11-8A93-001D09F2512C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/86E9A5AE-D1D4-DE11-A840-001D09F29169.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/86D8CF34-C9D4-DE11-82DF-001D09F24353.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/864E64F0-E4D4-DE11-9C54-001D09F241B9.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/864945B4-D6D4-DE11-B617-0019B9F709A4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/862A06E8-C9D4-DE11-8DBA-0030487A3C9A.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/8072F4F1-E4D4-DE11-A27A-001D09F2437B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/803D02EE-C4D4-DE11-A8C6-0019B9F709A4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/7EDBCC63-E3D4-DE11-8E03-000423D9863C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/7CAD9766-D7D4-DE11-8FCF-001D09F28F1B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/7A4D59C9-E7D4-DE11-9BD4-001D09F23944.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/7A27727B-D9D4-DE11-9FE5-001D09F2447F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/76DC0D6D-CDD4-DE11-AFD5-0030487C5CFA.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/76CCAA3C-DAD4-DE11-AB34-001D09F2A690.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/7647D07B-C8D4-DE11-A898-003048D374F2.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/761769A7-C5D4-DE11-8EC2-000423D992DC.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/74903955-E9D4-DE11-9D1A-003048D2C108.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/70BA0243-D0D4-DE11-83F0-0030486733D8.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/70AF5B3E-E6D4-DE11-B4FF-000423D99896.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/7082128F-CFD4-DE11-B623-000423D998BA.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/706E3050-DCD4-DE11-862D-001D09F24FEC.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/7016200F-E5D4-DE11-9485-000423D98634.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/6E6166EC-E6D4-DE11-8B62-001D09F28F1B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/6AC8900B-E5D4-DE11-B19A-000423D99996.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/66E47B73-DED4-DE11-A0EC-0019B9F730D2.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/62DE8096-CAD4-DE11-ADA7-0030486730C6.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/62D395F4-E4D4-DE11-BD8C-001D09F2AD84.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/5EB6E783-DED4-DE11-A620-001D09F29321.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/5E525332-C9D4-DE11-9330-001D09F2305C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/5A8A22D7-DFD4-DE11-B5F6-001D09F25109.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/56B6F564-E3D4-DE11-B206-000423D991D4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/56A47940-E1D4-DE11-ADBD-001D09F25438.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/5614AC31-DAD4-DE11-A836-001D09F2AD7F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/54FA3F81-E5D4-DE11-AEF9-000423D98BC4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/54A4E37B-D9D4-DE11-BDB5-001D09F251FE.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/521E2FCD-E7D4-DE11-9CB4-001D09F276CF.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/4C12C183-D4D4-DE11-B220-001D09F2841C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/4877D8F8-E1D4-DE11-8839-001D09F2A49C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/48064A18-E5D4-DE11-AC60-000423D94534.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/46F683B6-DDD4-DE11-82D8-001D09F253C0.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/46E7E8CE-D3D4-DE11-9C70-001617DBD224.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/46A4CB50-CBD4-DE11-B1FA-0019B9F707D8.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/44FB5DAE-D1D4-DE11-B21D-001D09F24493.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/442FA617-D3D4-DE11-85B4-001D09F251B8.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3EE58671-DED4-DE11-AC6D-001D09F24E39.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3E2D3870-CDD4-DE11-A413-0030487A322E.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3CFF3868-D2D4-DE11-A307-001D09F24934.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3AFF2C5A-C6D4-DE11-9451-003048D2BE12.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3AB76BF7-E1D4-DE11-A872-001D09F282F5.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3A777032-C4D4-DE11-9540-000423D9853C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3827CF50-CBD4-DE11-B8FD-003048D2C0F0.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3812FC21-D3D4-DE11-B331-001D09F2B2CF.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/36525778-D9D4-DE11-AB69-001D09F2514F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/36413EEC-D5D4-DE11-8147-001D09F24303.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3620684D-CBD4-DE11-AB6F-003048D2BF1C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/3082B96B-D2D4-DE11-BFE8-001D09F251B8.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/303F98BA-DDD4-DE11-8FEB-001D09F26C5C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/2C87FDED-C4D4-DE11-B1BA-003048D374F2.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/2C37D516-D8D4-DE11-8C4B-001D09F250AF.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/2C256FE6-DAD4-DE11-AF67-001D09F24DA8.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/2A5CA2B4-D6D4-DE11-AA87-001D09F2546F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/261106EE-C4D4-DE11-9ECC-001D09F29619.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/260A258A-D4D4-DE11-96F0-001D09F28EC1.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/224F72A6-C5D4-DE11-8D19-003048678098.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/204B0A79-E8D4-DE11-B5A8-001D09F2960F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/1EA853CB-E7D4-DE11-9C2E-001617C3B76A.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/1AF67208-DDD4-DE11-8A58-001D09F23A20.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/18B3FB98-DBD4-DE11-BD89-001D09F24F1F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/1854FEEB-D5D4-DE11-B3A0-001D09F2532F.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/1839C649-D0D4-DE11-A156-003048D2C1C4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/182420A6-C5D4-DE11-A519-0030487C6090.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/16EEEC93-CFD4-DE11-A853-001D09F24682.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/16745B6F-CDD4-DE11-A3B2-0030487D1BCC.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/141BDC97-D4D4-DE11-B046-001D09F34488.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/12F25713-C7D4-DE11-96F2-000423D985E4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/1226D268-D2D4-DE11-80F8-001D09F29619.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/107155A4-D9D4-DE11-9150-001D09F24259.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/104BF7ED-E6D4-DE11-8FEB-001D09F24FBA.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/1047E7E7-C9D4-DE11-B7E7-003048D2C0F4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/0E669C49-E1D4-DE11-AA03-001D09F2441B.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/0E47ECFB-D0D4-DE11-8653-003048D37580.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/0AE20D26-DFD4-DE11-BC90-001D09F29619.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/0A42418A-CED4-DE11-8E1C-001D09F24600.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/0A140E7F-E5D4-DE11-A313-0019DB29C5FC.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/08E27D02-D1D4-DE11-9A2E-003048D2BED6.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/082E6AF7-E4D4-DE11-A85D-001D09F27067.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/0808B376-E8D4-DE11-AD0A-0019B9F7312C.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/0655E990-CFD4-DE11-BAD6-0016177CA7A0.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/04771EED-E6D4-DE11-8F78-001617C3B5E4.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/04332996-DBD4-DE11-8FED-001D09F24493.root',
    'rfio:/castor/cern.ch/cms/store/data/BeamCommissioning09/Calo/RAW/v1/000/121/620/042AA98D-E0D4-DE11-A8ED-000423D98C20.root'
    )
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) #2805835 events total

#keep the logging output to a nice level
process.MessageLogger = cms.Service("MessageLogger")

#load the ECAL unpacker
process.load("CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi")
process.load("EventFilter.EcalRawToDigiDev.EcalUnpackerMapping_cfi")
process.load("EventFilter.EcalRawToDigiDev.EcalUnpackerData_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")

#StandardSequences configurations
#from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CCEcal/CRUZET2/CaloOnlineTools/EcalTools/python/ecalCosmicsHists_cfg.py?revision=1.14&view=markup
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
#process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
#process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")

#global tag: see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Global_Run_data
#from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CCEcal/CRUZET2/CaloOnlineTools/EcalTools/python/ecalCosmicsHists_cfg.py?revision=1.14&view=markup
#process.GlobalTag.globaltag = 'GR09_31X_V4P::All'

#unpack the GT raw data
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtBoardMapsConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtConfig_cff")
import EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi
process.gtDigis = EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi.l1GtUnpack.clone()
process.gtDigis.DaqGtInputTag = 'source'

#load the L1 filter
process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.l1Filter.algorithms = cms.vstring('L1_SingleEG2', 'L1_SingleEG5', 'L1_SingleEG8', 'L1_SingleEG10', 'L1_SingleEG12', 'L1_SingleEG15', 'L1_SingleEG20', 'L1_SingleEG25', 'L1_DoubleNoIsoEGBTBtight', 'L1_DoubleNoIsoEGBTBloose', 'L1_DoubleNoIsoEGTopBottom', 'L1_DoubleNoIsoEGTopBottomCen', 'L1_DoubleNoIsoEGTopBottomCen2', 'L1_DoubleNoIsoEGTopBottomCenVert')

#load the leakage filter
process.load("LEDSoakTools.LeakageFilter.leakagefilter_cfi")

#create a file to which to write the events
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('ECAL-triggered_121620.root'),
                               SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
                               #outputCommands = cms.untracked.vstring('drop *',
                                                                       #"keep *_*_eeDigiSkim_*")
                               )

#run the unpackers, filter on L1, and filter on BX
#process.p = cms.Path(process.ecalDigis*process.LeakageFilter)
process.p = cms.Path(process.ecalEBunpacker*process.gtDigis*process.l1Filter*process.LeakageFilter)

#write the resulting collections to the output file
process.e = cms.EndPath(process.out)
