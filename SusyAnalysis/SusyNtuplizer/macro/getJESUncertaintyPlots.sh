#!/bin/bash

if [ ! -d /data2/yohay/RA3/JESUncertainty ]
    then
    mkdir /data2/yohay/RA3/JESUncertainty
fi

for file in `rfdir $CASTOR_HOME/JESUncertainty/ | grep 4684pb-1_MET_18-Jan-12_skim_JESUncertainty_nominalBins_ | grep root | awk '{ print $9 }'`
  do
  rfcp $CASTOR_HOME/JESUncertainty/$file /data2/yohay/RA3/JESUncertainty
done

exit 0
