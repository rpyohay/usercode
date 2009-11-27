#!/bin/tcsh

cd /afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/prompt_feedback/CalibrationSequenceAnalyzer/cfg
eval `scramv1 runtime -csh`
cd -
cmsRun /afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/prompt_feedback/CalibrationSequenceAnalyzer/cfg/calibrationsequenceanalyzer_cfg_RUN_MIN_EVT-MAX_EVT.py
rfcp /afs/cern.ch/user/y/yohay/scratch0/tree_RUN-ped_1st_sample-events_MIN_EVT-MAX_EVT.root /castor/cern.ch/user/y/yohay/run_RUN_trees
