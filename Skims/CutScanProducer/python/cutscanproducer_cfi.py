import FWCore.ParameterSet.Config as cms

#the first 3 ECAL, HCAL, and track isolation values correspond to photon isolation working points
#defined for HLT studies (https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaWorkingPointsv3)
#from left to right, they go very loose (VL) --> loose (L) --> tight (T)
#the last point for ECAL and HCAL isolation and H/E is the RA3 2010 selection cut
#the 2nd point for sigmaIetaIeta is the RA3 2010 selection cut
#values of sigmaIetaIeta < 0.013 are tighter than the 2010 RA3 selection cut
#values of sigmaIetaIeta > 0.014 are possibly looser than the unprescaled triggers

CutScanProducer = cms.EDProducer('CutScanProducer',
                                 #photonTag = cms.untracked.InputTag("photons", "", "RECOCleaned"),
                                 ETMinScan = cms.vdouble(35.0, 43.0, 50.0), #GeV
                                 ECALIsoMaxPTMultiplierScan = cms.vdouble(0.012, 0.012, 0.012,
                                                                          0.006),
                                 ECALIsoMaxConstantScan = cms.vdouble(6.0, 5.5, 5.0, 4.2), #GeV
                                 HCALIsoMaxPTMultiplierScan = cms.vdouble(0.005, 0.005, 0.005,
                                                                          0.0025),
                                 HCALIsoMaxConstantScan = cms.vdouble(4.0, 3.5, 3.0, 2.2), #GeV
                                 HOverEMaxScan = cms.vdouble(0.1, 0.05),
                                 trackIsoMaxPTMultiplierScan = cms.vdouble(0.002, 0.002, 0.002),
                                 trackIsoMaxConstantScan = cms.vdouble(4.0, 3.5, 3.0), #GeV
                                 sigmaIetaIetaMaxScan = cms.vdouble(0.024, 0.013, 0.011)#,
                                 #doETMinScan = cms.untracked.bool(),
                                 #doECALIsoMaxScan = cms.untracked.bool(),
                                 #doHCALIsoMaxScan = cms.untracked.bool(),
                                 #doHOverEMaxScan = cms.untracked.bool(),
                                 #doTrackIsoMaxScan = cms.untracked.bool(),
                                 #doSigmaIetaIetaScan = cms.untracked.bool()
                                 )
