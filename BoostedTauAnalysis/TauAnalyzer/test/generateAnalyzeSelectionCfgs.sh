#!/bin/bash

#arguments
discriminator=$1
sample=$2

#macros
dir="${sample}_${discriminator}"

#make working directory
mkdir -p $dir
cd $dir

#for background case
if [ "$sample" = "WJets" ]
    then

    #set up CMSSW environment so cmsLs can be used
    eval `scramv1 runtime -sh`

    #macros
    nFilesPerJob=88
    first_file=1
    last_file=1000
    prefix="'root://eoscms//eos/cms"
    EOSDir="/store/user/friccita/WToMuNu_skim/"
    file_list=`cmsLs $EOSDir | grep root | awk '{ print $5 }'`
    identifier="${EOSDir}Summer12_WJetsToLNu_WMuNuSkim"

    #loop over number of input files
    input_file_list=""
    nJobs=1
    for iFile in `seq ${first_file} ${last_file}`
      do

      #if proposed file exists...
      file=`echo $file_list | sed -e "s%.*\(${identifier}_${iFile}_[0-9]*_[0-9A-Za-z]*\.root\).*%\1%g"`
      file_list_mod=`echo $file_list | sed -e "s%\(.*\)%\1%g"`
      if [ "$file" != "$file_list_mod" ]
	  then
	  file_name="${prefix}${file}'"
	  file_name=`echo $file_name | sed -e "s%/%\/%g" | sed -e "s%\.%\.%g"`

          #format file name
	  if [ `expr $iFile % $nFilesPerJob` != "0" ] && [ "$iFile" != "$last_file" ]
	      then

              #append comma
	      file_name="${file_name},\n    "
	      input_file_list="${input_file_list}${file_name}"

          #create cfg and sh files if you're at the max number of files per job
	  else
	      input_file_list="${input_file_list}${file_name}"
	      sed -e "s%FILES%$input_file_list%" -e "s%JOB%$nJobs%" -e "s%XXX%$discriminator%" ../analyzeSelectionTemplate.py > analyzeSelectionTemplate_${discriminator}_${nJobs}.py
	      sed -e "s%DIR%$dir%" -e "s%TAUISOWORKINGPOINT%$discriminator%" -e "s%JOB%$nJobs%g" ../analyzeSelectionTemplate.sh > analyzeSelectionTemplate_${discriminator}_${nJobs}.sh
	      chmod a+x analyzeSelectionTemplate_${discriminator}_${nJobs}.sh
	      bsub -q 1nh -J analyzeSelectionTemplate_${discriminator}_${nJobs} < analyzeSelectionTemplate_${discriminator}_${nJobs}.sh
	      input_file_list=""
	      nJobs=`expr $nJobs + 1`
	  fi
      fi
    done
fi

#for signal case
if [ "$sample" = "Wh1" ]
    then
    input_file_list="'file:/data1/yohay/NMSSMHiggs_WH_files1-250_24Sep12.root',\n    'file:/data1/yohay/NMSSMHiggs_WH_files251-500_24Sep12.root',\n    'file:/data1/yohay/NMSSMHiggs_WH_files501-750_24Sep12.root',\n    'file:/data1/yohay/NMSSMHiggs_WH_files751-1000_24Sep12.root'"
    input_file_list=`echo $input_file_list | sed -e "s%/%\/%g" | sed -e "s%\.%\.%g"`
    nJobs=0
    sed -e "s%FILES%$input_file_list%" -e "s%JOB%$nJobs%" -e "s%XXX%$discriminator%" ../analyzeSelectionTemplate.py > analyzeSelectionTemplate_${discriminator}_${nJobs}.py
fi

exit 0
