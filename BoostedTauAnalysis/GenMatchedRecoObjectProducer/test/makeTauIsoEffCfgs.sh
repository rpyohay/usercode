#!/bin/bash

#input files
inputFiles="NMSSMHiggs_gg_skim_1000Files Summer12DrellYan"

#sequences
sequences=( signalSequence ZTauTauSequence )
medIsoSequences=( medIsoSignalSequence medIsoZTauTauSequence )

#generate a cfg file for each input file
i=0
for iFile in $inputFiles
  do
  tauRecoSeq=""
  process="SKIM"
  if [ -n "`echo ${sequences[${i}]} | grep ZTauTau`" ]
      then
      #tauRecoSeq="process.PFTau*"
      #process="OWNPARTICLES"
      process="RECO"
  fi
  #sed -e "s%FILE_NAME%$iFile%" -e "s%PROCESS%${process}%g" -e "s%TAURECOSEQUENCE%${tauRecoSeq}%" -e "s%SEQUENCE%${sequences[${i}]}%" tauIsoEff.py > tauIsoEff_${sequences[${i}]}.py
  sed -e "s%FILE_NAME%$iFile%" -e "s%PROCESS%${process}%g" -e "s%TAURECOSEQUENCE%${tauRecoSeq}%" -e "s%SEQUENCE%${medIsoSequences[${i}]}%" tauIsoEff.py > tauIsoEff_${medIsoSequences[${i}]}.py
  i=`expr $i + 1`
done

exit 0

