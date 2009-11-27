#!/bin/tcsh

#submit 1 job of 5000 events

#go to the working directory
cd ~/LED/CMSSW_3_2_1/src/prompt_feedback/CalibrationSequenceAnalyzer/batch_jobs/submission

#generate the cfi file
if ( ! -e ../../python/calibrationsequenceanalyzer_cfi_THIS_RUN_THIS_MIN_EVT-THIS_MAX_EVT.py ) then
    sed -e "s/RUN/THIS_RUN/g" -e "s/MIN_EVT/THIS_MIN_EVT/g" -e "s/MAX_EVT/THIS_MAX_EVT/g" ../../python/calibrationsequenceanalyzer_cfi.py > ../../python/calibrationsequenceanalyzer_cfi_THIS_RUN_THIS_MIN_EVT-THIS_MAX_EVT.py
endif

#generate the cfg file
if ( ! -e ../../cfg/calibrationsequenceanalyzer_cfg_THIS_RUN_THIS_MIN_EVT-THIS_MAX_EVT.py ) then
    sed -e "s/MAX_EVTS/THIS_MAX_EVTS/g" -e "s/NUM_SKIPPED_EVTS/THIS_NUM_SKIPPED_EVTS/g" -e "s/RUN/THIS_RUN/g" -e "s/MIN_EVT/THIS_MIN_EVT/g" -e "s/MAX_EVT/THIS_MAX_EVT/g" ../../cfg/calibrationsequenceanalyzer_cfg.py > ../../cfg/calibrationsequenceanalyzer_cfg_THIS_RUN_THIS_MIN_EVT-THIS_MAX_EVT.py
endif

#generate the script
if ( ! -e ../../sh/calibrationsequenceanalyzer_cfg_THIS_RUN_THIS_MIN_EVT-THIS_MAX_EVT.sh ) then
    sed -e "s/RUN/THIS_RUN/g" -e "s/MIN_EVT/THIS_MIN_EVT/g" -e "s/MAX_EVT/THIS_MAX_EVT/g" ../../sh/calibrationsequenceanalyzer_cfg.sh > ../../sh/calibrationsequenceanalyzer_cfg_THIS_RUN_THIS_MIN_EVT-THIS_MAX_EVT.sh
endif

#make the script executable
chmod a+x ../../sh/calibrationsequenceanalyzer_cfg_THIS_RUN_THIS_MIN_EVT-THIS_MAX_EVT.sh

#submit the job
cd ../../sh
bsub -q 1nh -J tree_THIS_RUN-alpha_greater_than_1-events_THIS_MIN_EVT-THIS_MAX_EVT < calibrationsequenceanalyzer_cfg_THIS_RUN_THIS_MIN_EVT-THIS_MAX_EVT.sh

#exit
exit 0
