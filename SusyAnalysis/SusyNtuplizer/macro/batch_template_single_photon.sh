#!/bin/bash

cd /afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src
eval `scramv1 runtime -sh`
cd -
for file in `rfdir $CASTOR_HOME/424p2/DIR | awk '{ print $9 }'`
  do
  rfcp $CASTOR_HOME/424p2/DIR/$file .
done
file_list=`ls | grep susyEvent | sed -e "s%\(susyEvent_[0-9_A-Za-z]*\.root\)%pars.input.push_back(\"\1\");%g"`
file_list_2=`echo $file_list | sed -e "s%\.%\\\.%"`
sed -e "s%FILES%$file_list_2%" -e "s%DATASET%DIR%" /afs/cern.ch/user/y/yohay/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro/ana_standalone_template_single_photon.C > ana_standalone_single_photon_MC.C
root -b -q ana_standalone_single_photon_MC.C
rfcp DIR_JSON_HLT_PV_single_photon_skim.txt $CASTOR_HOME
rfcp DIR_JSON_HLT_PV_single_photon_skim.root $CASTOR_HOME
rm DIR_JSON_HLT_PV_single_photon_skim.txt DIR_JSON_HLT_PV_single_photon_skim.root susyEvent* ana_standalone_single_photon_MC.C*

exit 0
