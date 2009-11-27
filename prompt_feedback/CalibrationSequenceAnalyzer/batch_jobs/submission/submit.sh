#!/bin/tcsh

#submit all jobs for a particular run

#job parameters that don't change
set MAX_EVTS = 5000
set NUM_BLOCKS = 20
set RUN = 118358

#go to the working directory
cd ~/LED/CMSSW_3_2_1/src/prompt_feedback/CalibrationSequenceAnalyzer/batch_jobs/submission
cmsenv

#make the CASTOR directory in which to store the trees
set LINE = `rfdir | grep "run_${RUN}_trees"`
if ( "$LINE" == "" ) then
    rfmkdir /castor/cern.ch/user/y/yohay/run_${RUN}_trees
endif

#loop over all blocks of MAX_EVTS events
set SEQ = `seq $NUM_BLOCKS`
foreach block ($SEQ)

    #job parameters that change for each job
    @ MIN_EVT = ( $block - 1 ) * $MAX_EVTS + 1
    @ MAX_EVT = $MIN_EVT + $MAX_EVTS - 1
    @ NUM_SKIPPED_EVTS = ( $block - 1) * $MAX_EVTS

    #generate the submission script
    if ( ! -e submit_1_job_$RUN-alpha_greater_than_1-events_$MIN_EVT-$MAX_EVT.sh ) then
	sed -e "s/THIS_MAX_EVTS/$MAX_EVTS/g" -e "s/THIS_RUN/$RUN/g" -e "s/THIS_MIN_EVT/$MIN_EVT/g" -e "s/THIS_MAX_EVT/$MAX_EVT/g" -e "s/THIS_NUM_SKIPPED_EVTS/$NUM_SKIPPED_EVTS/g" submit_1_job.sh > submit_1_job_$RUN-alpha_greater_than_1-events_$MIN_EVT-$MAX_EVT.sh
    endif

    #make the script executable
    chmod a+x submit_1_job_$RUN-alpha_greater_than_1-events_$MIN_EVT-$MAX_EVT.sh

    #run the script
    ./submit_1_job_$RUN-alpha_greater_than_1-events_$MIN_EVT-$MAX_EVT.sh

end

#exit
exit 0
