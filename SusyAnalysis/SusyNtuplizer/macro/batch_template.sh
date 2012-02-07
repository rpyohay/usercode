#!/bin/bash

cd /afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src
eval `scramv1 runtime -sh`
cd -
num_files_head=NUM_FILES_HEAD
for file in `rfdir $CASTOR_HOME/424p2/DIR | awk '{ print $9 }' | head -n $num_files_head | tail -n 50`
  do
  rfcp $CASTOR_HOME/424p2/DIR/$file .
  echo "$file"
done
file_list=`ls | grep susyEvent | sed -e "s%\(susyEvent_[0-9_A-Za-z]*\.root\)%pars.input.push_back(\"\1\");%g"`
file_list_2=`echo $file_list | sed -e "s%\.%\\\.%"`
sed -e "s%FILES%$file_list_2%" -e "s%DATASET%DIR%" -e "s%NUM%JOB%" /afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/ana_standalone_template.C > ana_standalone_MC_JOB.C
root -b -q ana_standalone_MC_JOB.C
rfcp DIR_JSON_HLT_PV_skim_v3_JOB.txt $CASTOR_HOME
rfcp DIR_JSON_HLT_PV_skim_v3_JOB.root $CASTOR_HOME
rm DIR_JSON_HLT_PV_skim_v3_JOB.txt DIR_JSON_HLT_PV_skim_v3_JOB.root susyEvent* ana_standalone_MC_JOB.C*

exit 0
