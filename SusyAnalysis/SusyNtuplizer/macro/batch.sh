#!/bin/bash

cd /afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src
eval `scramv1 runtime -sh`
cd -
echo "scping..."
scp -B -o GSSAPIAuthentication=yes -o StrictHostKeyChecking=no ndpc3:/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-May10ReReco-v1/Photon/Runs160442-163869/susyEvent_ALL_1.root .
echo "done scping, running..."
root -l -q /afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/countTriggers.C
scp count_May10ReReco_PromptRecov4_batchTest.txt count_May10ReReco_PromptRecov4_batchTest.root /afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer
echo "done running, rming..."
rm susyEvent_ALL_1.root count_May10ReReco_PromptRecov4_batchTest.txt count_May10ReReco_PromptRecov4_batchTest.root
echo "done"

exit 0
