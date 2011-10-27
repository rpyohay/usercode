#!/bin/tcsh

cd ../../PhysicsTools/Utilities
make
cd ../../GMSBTools/Filters
make
cd ../../SusyAnalysis/SusyNtuplizer/jec/JetMETObjects
make
cd ../../macro
make

exit 0
