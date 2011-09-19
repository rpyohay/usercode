#!/bin/bash

count=0
#for line in `jot - 3 6446` #gg
#for line in `jot - 6449 10700` #eg
#for line in `jot - 95100 99380` #ff
#for line in `jot - 1 4280` #ff
#for line in `jot - 1 4251` #eg
#for line in `jot - 1 6444` #gg
#for line in `jot - 1 84395` #ee
for line in `jot - 1 77315` #ee
do
#  run_evt=`head -n $line /Users/rachelyohay/RA3/data/1140pb-1_categorized.txt | tail -n 1 | sed -e "s%\([0-9]* [0-9]*\).*%\1%"`
#    run_evt_bin=`head -n $line /Users/rachelyohay/RA3/data/ee.txt | tail -n 1`
    run_evt=`head -n $line /Users/rachelyohay/RA3/data/ee_Dave_usingDetectorP4.txt | tail -n 1`
    found=`grep "$run_evt" /Users/rachelyohay/RA3/data/ee.txt`
    bin=`echo $found | sed -e "s%[0-9]* [0-9]* \([0-9]\)%\1%"`
#    found=`grep "$run_evt" /Users/rachelyohay/RA3/data/ff_Dave_usingDetectorP4.txt`
#    run_evt=`echo $run_evt_bin | sed -e "s%\([0-9]* [0-9]*\).*%\1%"`
#    bin=`echo $run_evt_bin | sed -e "s%[0-9]* [0-9]* \(.*\)%\1%"`
    if [ "$bin" != "4" ]
#	then
#	found=`grep "$run_evt" /Users/rachelyohay/RA3/data/ee_Dave_usingDetectorP4.txt`
#  run_evt=`head -n $line /Users/rachelyohay/RA3/data/gg_Dave.txt | tail -n 1 | sed -e "s%\([0-9]* [0-9]*\).*%\1%"`
#  found=`head -n 6446 /Users/rachelyohay/RA3/data/1140pb-1_categorized.txt | tail -n 6444 | grep "$run_evt"`
#	if [ "$found" == "" ]
	then
	    count=`expr $count + 1`
	    echo $run_evt
#	fi
    fi
done
echo "$count events"

#evts[0]="Run: 167102  Event: 41367892"
#evts[1]="Run: 167284  Event: 1044948847"
#evts[2]="Run: 167098  Event: 226514608"
#evts[3]="Run: 167282  Event: 364681666"
#evts[4]="Run: 167807  Event: 163352644"
#evts[5]="Run: 167807  Event: 1071372467"
#evts[6]="Run: 167807  Event: 1711604905"
#evts[7]="Run: 167807  Event: 1920901465"
#evts[8]="Run: 167807  Event: 2002886211"
#evts[9]="Run: 167830  Event: 633065492"
#evts[10]="Run: 167830  Event: 872265467"
#evts[11]="Run: 167830  Event: 1074870406"
#evts[12]="Run: 167898  Event: 498944544"
#evts[13]="Run: 167898  Event: 841183462"
#evts[14]="Run: 167898  Event: 1939001688"
#evts[15]="Run: 167674  Event: 381879303"
#evts[16]="Run: 167676  Event: 78292789"
#evts[17]="Run: 167754  Event: 88484623"
#for iEvt in `seq 0 17`
#  do
#  evt=${evts[$iEvt]}
#  gg=`grep "$evt" ggEvents_Data2011A_ToRun167913_NoPileupCorr_Filter_Photon.txt`
#  if [ "$gg" != "" ]
#      then
#      echo "$evt is gg"
#  fi
#  eg=`grep "$evt" egEvents_Data2011A_ToRun167913_NoPileupCorr_Filter_Photon.txt`
#  if [ "$eg" != "" ]
#      then
#      echo "$evt is eg"
#  fi
#  ee=`grep "$evt" eeEvents_Data2011A_ToRun167913_NoPileupCorr_Filter_Photon.txt`
#  if [ "$ee" != "" ]
#      then
#      echo "$evt is ee"
#  fi
#  ff=`grep "$evt" ffEvents_Data2011A_ToRun167913_NoPileupCorr_Filter_Photon.txt`
#  if [ "$ff" != "" ]
#      then
#      echo "$evt is ff"
#  fi
#done

exit 0
