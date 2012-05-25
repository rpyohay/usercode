#!/bin/bash

#set up CMSSW environment
cd ~/scratch0/CMSSW_4_2_4_patch2/src
eval `scramv1 runtime -sh`

#generate the random shifts
cd SusyAnalysis/SusyNtuplizer/macro
echo "Generating random shifts..."
root -l -b -q 'generateRandomShiftsForJESUncertainty.C+(100, "randomShiftsForJESUncertainty.txt")'

#create an analyzer job for each toy
n_toys=`cat randomShiftsForJESUncertainty.txt | wc -l`
max_toy=`expr $n_toys - 1`
for i_toy in `seq 0 $max_toy`
  do
  line=`expr $i_toy + 1`
  shift=`cat randomShiftsForJESUncertainty.txt | head -n $line | tail -n 1`
  sed -e "s%ITOY%$i_toy%g" -e "s%ISHIFT%$shift%g" generateToysForJESUncertainty.sh > batch/generateToysForJESUncertainty_$i_toy.sh
  chmod a+x batch/generateToysForJESUncertainty_$i_toy.sh
  echo "Submitting job generateToysForJESUncertainty_$i_toy..."
  bsub -q 1nd -R "rusage[pool=25000]" -J generateToysForJESUncertainty_$i_toy < batch/generateToysForJESUncertainty_$i_toy.sh
done

exit 0
