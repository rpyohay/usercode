#!/bin/bash

#setup CMSSW environment
cd ~/scratch0/CMSSW_4_2_4_patch2/src
eval `scramv1 runtime -sh`
cd ~/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro

#get total size of each dataset
for dir in `rfdir $CASTOR_HOME/424p2 | awk '{ print $9 }'`
  do
  sum=0
  for file_size in `rfdir $CASTOR_HOME/424p2/$dir | awk '{ print $5 }'`
    do
    sum=`expr $sum + $file_size`
  done
  echo "Size of $dir: $sum"

#  if [ "$dir" != "ntuple_DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2" ]
#      then

#      #rfcp files locally
#      mkdir /data2/yohay/RA3/tmp/$dir
#      for file in `rfdir $CASTOR_HOME/424p2/$dir | awk '{ print $9 }'`
#	do
#	rfcp $CASTOR_HOME/424p2/$dir/$file /data2/yohay/RA3/tmp/$dir
#      done
#
#      #process dataset
#      file_list=`ls /data2/yohay/RA3/tmp/$dir | sed -e "s%\(susyEvent_[0-9_A-Za-z]*\.root\)%pars.input.push_back(\"/data2/yohay/RA3/tmp/$dir/\1\");%g"`
#      file_list_2=`echo $file_list | sed -e "s%\.%\\\.%"`
#      sed -e "s%FILES%$file_list_2%" -e "s%DATASET%$dir%" ana_standalone_template.C > ana_standalone_MC.C
#      root -l -q ana_standalone_MC.C
#
#      #delete files
##      rm -rf /data2/yohay/RA3/tmp/$dir
#  else
#      echo "Skipping $dir"
#  fi
done

exit 0
