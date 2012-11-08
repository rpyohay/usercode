#!/bin/bash

jobDir="/afs/cern.ch/user/y/yohay/CMSSW_5_2_4_patch3_20Sep12/src/BoostedTauAnalysis/TauAnalyzer/test/DIR"
fileNamePrefix="analyzeSelectionTemplate_TAUISOWORKINGPOINT_JOB"
EDMFile="DIR_selected_events_JOB.root"

cd $jobDir
eval `scramv1 runtime -sh`
cd -
cp $jobDir/$fileNamePrefix.py .
cmsRun $fileNamePrefix.py >& $fileNamePrefix.txt
cp $fileNamePrefix.txt $jobDir
cmsStage -f $EDMFile /store/user/yohay/
rm $fileNamePrefix.* NMSSMSignal_MuProperties_JOB.root $EDMFile

exit 0
