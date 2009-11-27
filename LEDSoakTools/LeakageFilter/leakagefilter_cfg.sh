#!/bin/csh

cd /afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/LEDSoakTools/LeakageFilter
eval `scramv1 runtime -csh`
cd -
cmsRun /afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/LEDSoakTools/LeakageFilter/leakagefilter_cfg.py
rfcp ECAL-triggered_121620.root /castor/cern.ch/user/y/yohay/
