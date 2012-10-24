#!/bin/bash

#arguments
discriminator=$1
sample=$2
minJob=$3
maxJob=$4

#macros
filePrefix="${sample}_${discriminator}/analyzeSelectionTemplate_${discriminator}_"

#totals
nWMuNuTot=0
nHLTTot=0
nWMuonPTTot=0
nWMuonIsoTot=0
nJetTot=0
nTauMuonPTTot=0
nTauMuonSoftTot=0
nMuHadIsoTot=0

#loop over job outputs
for iJob in `seq $minJob $maxJob`
  do
  echo "Job #$iJob"

  #get the numbers of events passing each cut
  nWMuNu=`grep genWMuNuSelector ${filePrefix}${iJob}.txt | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nHLT=`grep IsoMu24eta2p1Selector ${filePrefix}${iJob}.txt | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nWMuonPT=`grep WMuonPTSelector ${filePrefix}${iJob}.txt | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nWMuonIso=`grep WIsoMuonSelector ${filePrefix}${iJob}.txt | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nJet=`grep jetSelector ${filePrefix}${iJob}.txt | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nTauMuonPT=`grep tauMuonPTSelector ${filePrefix}${iJob}.txt | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nTauMuonSoft=`grep tauMuonSelector ${filePrefix}${iJob}.txt | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`
  nMuHadIso=`grep muHadIsoTauSelector ${filePrefix}${iJob}.txt | head -n 1 | sed -e "s%TrigReport[ ]*1[ ]*0[ ]*[0-9]*[ ]*\([0-9]*\).*%\1%"`

  #increment totals
  nWMuNuTot=`expr $nWMuNuTot + $nWMuNu`
  nHLTTot=`expr $nHLTTot + $nHLT`
  nWMuonPTTot=`expr $nWMuonPTTot + $nWMuonPT`
  nWMuonIsoTot=`expr $nWMuonIsoTot + $nWMuonIso`
  nJetTot=`expr $nJetTot + $nJet`
  nTauMuonPTTot=`expr $nTauMuonPTTot + $nTauMuonPT`
  nTauMuonSoftTot=`expr $nTauMuonSoftTot + $nTauMuonSoft`
  nMuHadIsoTot=`expr $nMuHadIsoTot + $nMuHadIso`
done

#print totals
echo "nWMuNuTot = $nWMuNuTot"
echo "nHLTTot = $nHLTTot"
echo "nWMuonPTTot = $nWMuonPTTot"
echo "nWMuonIsoTot = $nWMuonIsoTot"
echo "nJetTot = $nJetTot"
echo "nTauMuonPTTot = $nTauMuonPTTot"
echo "nTauMuonSoftTot = $nTauMuonSoftTot"
echo "nMuHadIsoTot = $nMuHadIsoTot"

exit 0