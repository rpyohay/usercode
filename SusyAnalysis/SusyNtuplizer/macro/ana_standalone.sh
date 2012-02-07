#!/bin/bash

cd ~/scratch0/CMSSW_4_2_4_patch2/src/SusyAnalysis/SusyNtuplizer/macro
eval `scramv1 runtime -sh`

total=0
queue="8nh"
for dir in `rfdir $CASTOR_HOME/424p2 | awk '{ print $9 }'`
  do
  if [ "$dir" = "ntuple_DYToEE_M-20_TuneZ2_7TeV-pythia6-Summer11-PU_S3_START42_V11-v2" ] || [ "$dir" = "ntuple_DiPhotonBorn_Pt-10To25_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2" ] || [ "$dir" = "ntuple_DiPhotonBorn_Pt-250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2" ] || [ "$dir" = "ntuple_DiPhotonBorn_Pt-25To250_7TeV-pythia6-Summer11-PU_S4_START42_V11-v2" ] || [ "$dir" = "ntuple_v2_TT_TuneZ2_7TeV-pythia6-tauola-Summer11-PU_S3_START42_V11-v2" ]
      then
      sum=0
      for file_size in `rfdir $CASTOR_HOME/424p2/$dir | awk '{ print $5 }'`
	do
	sum=`expr $sum + $file_size`
      done
      sum_MB=`expr $sum / 1000000`
      sum_GB=`expr $sum_MB / 1000`
      sum_GB_rounded_up=`expr $sum_GB + 1`
      total=`expr $total + $sum_GB_rounded_up`
      echo "Approximate size of $dir: $sum_GB_rounded_up GB"
      queue="1nd"
      sum_MB=`expr $sum_MB + 20000`
      if [ $sum_MB -gt 80000 ]
	  then
	  sum_MB=150000
	  queue="1nd"
      elif [ $sum_MB -gt 120000 ]
	  then
	  sum_MB=200000
	  queue="2nd"
      fi
      queue="8nh"
      sum_MB=100000
      num_files=`rfdir $CASTOR_HOME/424p2/$dir | awk '{ print $9 }' | wc -l`
      num_jobs=`expr $num_files / 50`
      num_jobs=`expr $num_jobs + 1`
      for job in `seq 1 $num_jobs`
	do
	let num_files_head=$job*50
	sed -e "s%DIR%$dir%g" -e "s%NUM_FILES_HEAD%$num_files_head%" -e "s%JOB%$job%g" batch_template.sh > batch_${dir}_${job}.sh
	chmod a+x batch_${dir}_${job}.sh
	echo "Submitting job ana_standalone_${dir}_${job}..."
	echo "bsub -R \"rusage[pool=$sum_MB]\" -q $queue -J ana_standalone_${dir}_${job} < batch_${dir}_${job}.sh"
	bsub -R "rusage[pool=$sum_MB]" -q $queue -J ana_standalone_${dir}_${job} < batch_${dir}_${job}.sh
	rm batch_${dir}_${job}.sh*
      done
  else
      echo "Skipping $dir"
  fi
done
echo "Approximate total size: $total GB"

exit 0
