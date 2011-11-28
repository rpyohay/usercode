#!/bin/bash

cd ~/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro
eval `scramv1 runtime -sh`

total=0
for dir in `rfdir $CASTOR_HOME | grep Summer11 | grep root | grep single_photon | awk '{ print $9 }'`
#`rfdir $CASTOR_HOME/424p2 | awk '{ print $9 }'`
  do
  sum=0
  for file_size in `rfdir $CASTOR_HOME/$dir | awk '{ print $5 }'`
#`rfdir $CASTOR_HOME/424p2/$dir | awk '{ print $5 }'`
    do
    sum=`expr $sum + $file_size`
  done
#  echo "Size of $dir: $sum"
#  sum=`expr $sum + 20000000000`
  sum_MB=`expr $sum / 1000000`
  sum_GB=`expr $sum_MB / 1000`
  sum_GB_rounded_up=`expr $sum_GB + 1`
  total=`expr $total + $sum_GB_rounded_up`
  echo "Approximate size of $dir: $sum_GB_rounded_up GB"
#  if [ $sum_MB -gt 100000 ]
#      then
#      sum_MB=100000
#  fi
#  sum_MB=150000
#  if [  ]
#  if [ "$dir" = "ntuple_QCD_Pt-40_doubleEMEnriched_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2" ]
#      then
#      sed -e "s%DIR%$dir%g" batch_template.sh > batch_$dir.sh
#      chmod a+x batch_$dir.sh
#      echo "Submitting job ana_standalone_$dir..."
#      queue="2nd"
#      echo "bsub -R \"rusage[pool=$sum_MB]\" -q $queue -J ana_standalone_$dir < batch_$dir.sh"
#      bsub -R "rusage[pool=$sum_MB]" -q $queue -J ana_standalone_$dir < batch_$dir.sh
#      sed -e "s%DIR%$dir%g" batch_template_single_photon.sh > batch_single_photon_$dir.sh
#      chmod a+x batch_single_photon_$dir.sh
#      echo "Submitting job ana_standalone_single_photon_$dir..."
#      queue="2nd"
#      echo "bsub -R \"rusage[pool=$sum_MB]\" -q $queue -J ana_standalone_single_photon_$dir < batch_single_photon_$dir.sh"
#      bsub -R "rusage[pool=$sum_MB]" -q $queue -J ana_standalone_single_photon_$dir < batch_single_photon_$dir.sh
#  else
#      echo "Skipping $dir"
#  fi
done
echo "Approximate total size: $total GB"

exit 0
