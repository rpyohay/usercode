#!/bin/tcsh

cd ../../PhysicsTools/Utilities
#make clean
make
cd ../TagAndProbe
#make clean
make
cd ../../GMSBTools/Filters
#make clean
make
cd ../../SusyAnalysis/SusyNtuplizer/jec/JetMETObjects
#make clean
make
cd ../../macro
#make clean
make

exit 0
