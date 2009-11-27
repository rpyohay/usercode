import FWCore.ParameterSet.Config as cms

LMMapMaker = cms.EDAnalyzer('LMMapMaker',
                            mapFileName = cms.string("/afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/prompt_feedback/diffusing_sphere_index.txt"),
                            dee1FileName = cms.string("/afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/prompt_feedback/LM_regions_dee1_formatted.txt"),
                            dee2FileName = cms.string("/afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/prompt_feedback/LM_regions_dee2_formatted.txt"),
                            dee3FileName = cms.string("/afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/prompt_feedback/LM_regions_dee3_formatted.txt"),
                            dee4FileName = cms.string("/afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/prompt_feedback/LM_regions_dee4_formatted.txt")
)
