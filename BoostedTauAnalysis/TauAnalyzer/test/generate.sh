#!/bin/bash

#parse arguments
if [ $# -ne 3 ]
    then
    echo "Usage: ./generate.sh script_name dir_name num_jobs"
    exit 0
fi
script_name=$1
dir_name=$2
num_jobs=$3

#make directory on EOS
cmsMkdir /store/user/yohay/${dir_name}

#generate cfg and sh files for each job
for i in `seq 1 $num_jobs`
  do
  sed -e "s%JOBNUM%$i%g" ${script_name}.py > ${script_name}_${i}.py
  sed -e "s%JOBNUM%$i%g" -e "s%DIRNAME%$dir_name%g" -e "s%SCRIPTNAME%$script_name%g" ${script_name}.sh > ${script_name}_${i}.sh
  chmod a+x ${script_name}_${i}.sh
  bsub -q 8nh -J ${script_name}_${i} < ${script_name}_${i}.sh
done

exit 0
