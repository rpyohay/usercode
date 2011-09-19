#!/bin/bash

cd ~/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro
eval `scramv1 runtime -sh`

file_list="/data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-May10ReReco-v1/Photon/Runs160442-163869/susyEvent_ALL_1.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-May10ReReco-v1/Photon/Runs160442-163869/susyEvent_ALL_2.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Run166438/susyEvent_ALL.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs165088-166346/susyEvent_ALL_1.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs165088-166346/susyEvent_ALL_2.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166374-166486/susyEvent_ALL.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166502-166530/susyEvent_ALL.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166554-166787/susyEvent_ALL.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166839-166911/susyEvent_ALL.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs166921-167078/susyEvent_ALL.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs167098-167284/susyEvent_ALL.root /data/ndpc3/c/dmorse/RA3/SusyNtuples/cms423v5_v1/Run2011A-PromptReco-v4/Photon/Runs167551-167913/susyEvent_ALL.root"

count=0
for file in $file_list
  do
  count=`expr $count + 1`
  sed -e "s%FILE%$file%g" -e "s%OUTPUT%_$count%g" localCountTriggers_template.C > localCountTriggers_${count}.C
  root -l -q localCountTriggers_${count}.C
done

exit 0
