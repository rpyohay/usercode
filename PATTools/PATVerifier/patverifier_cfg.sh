#!/bin/csh

cd /afs/cern.ch/user/y/yohay/scratch0/CMSSW_3_3_6/src/PATTools/PATVerifier
eval `scramv1 runtime -csh`
cd -
cmsRun /afs/cern.ch/user/y/yohay/scratch0/CMSSW_3_3_6/src/PATTools/PATVerifier/patverifier_cfg.py
rfcp PATVerifier_test_15-20_GeV.root /castor/cern.ch/user/y/yohay/PhotonJet_Pt15to20-Summer09-MC_31X_V9_7TeV-v1-GEN-SIM-RECO
