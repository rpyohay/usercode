#!/bin/csh

cd /afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/LEDGapTools/TEFilter
eval `scramv1 runtime -csh`
cd -
cmsRun /afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/LEDGapTools/TEFilter/tefilter_cfg_118358.py
rfcp pedestal_digis_118358.root /castor/cern.ch/user/y/yohay/
