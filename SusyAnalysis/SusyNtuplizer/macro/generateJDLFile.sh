#!/bin/bash

#arguments
scan=$1
script=$2
output_path=$3
path_before_scan=$4
path_after_scan=$5

#output directory
output_dir="$output_path$scan"

#CONDOR output directory
condor_output_dir="$output_dir/CONDOROutput"

#what to grep for to find the right input file?
found=`echo $script | grep getAcceptance`
if [ "$found" != "" ]
    then
    what_to_grep="skim"
else
    what_to_grep="tree"
fi

#no. jobs
num_jobs=`ls $path_before_scan$scan$path_after_scan | grep $what_to_grep | grep root | wc -l`

#make output directory
if [ ! -d "$output_dir" ]
    then
    mkdir $output_dir
fi

#make CONDOR output directory
if [ ! -d "$condor_output_dir" ]
    then
    mkdir $condor_output_dir
fi

#generate JDL file
cat <<EOF > $script.jdl
universe = vanilla
Executable = $script.sh

Requirements   = Memory >= 199 && OpSys == "LINUX" && Arch != "DUMMY" && TARGET.FileSystemDomain == "fnal.gov"

Output = $condor_output_dir/batch_\$(cluster)_\$(process).stdout
Error  = $condor_output_dir/batch_\$(cluster)_\$(process).stderr
Log    = $condor_output_dir/batch_\$(cluster)_\$(process).condor

notification = Never

Arguments = \$(process)
Queue $num_jobs
EOF

exit 0
