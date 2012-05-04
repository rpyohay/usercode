#!/bin/bash

#arguments
scan=$1
script=$2
path_before_scan=$3
path_after_scan=$4

#what to grep for to find the right input file?
found=`echo $script | grep getAcceptance`
if [ "$found" != "" ]
    then
    what_to_grep="skim"
else
    what_to_grep="root"
fi

#generate CONDOR executable
cat <<EOF > $script.sh
#!/bin/bash
source /uscmst1/prod/sw/cms/bashrc prod
cd ~/CMSSW_4_2_8/src/SusyAnalysis/SusyNtuplizer/macro
eval \`scramv1 runtime -sh\`
job=\$1
dir="$path_before_scan$scan$path_after_scan"
file_list=\`ls -alh \$dir | grep $what_to_grep | grep root | awk '{ print \$9 }'\`
file_array=( \$file_list )
file=\${file_array[\$job]}
label=\`echo \$file | sed -e "s%[a-zA-Z]*_\(.*\)\.root%\1%"\`
sed -e "s%JOB%\$job%g" -e "s%SCAN%$scan%g" -e "s%FILE%\$dir\$file%g" -e "s%LABEL%\$label%g" $script.C > ${script}_\$job.C
root -l -b -q ${script}_\$job.C
rm ${script}_\$job.C
exit 0
EOF
chmod a+x $script.sh

exit 0
