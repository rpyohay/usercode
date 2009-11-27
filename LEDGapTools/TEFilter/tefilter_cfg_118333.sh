#!/bin/csh

cd /afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/LEDGapTools/TEFilter
eval `scramv1 runtime -csh`
cd -
cmsRun /afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/LEDGapTools/TEFilter/tefilter_cfg_118333.py
rfcp LED_digis_118333.root /castor/cern.ch/user/y/yohay/
