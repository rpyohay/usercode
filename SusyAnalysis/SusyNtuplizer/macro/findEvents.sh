#!/bin/bash

count=0
#for line in `jot - 3 6446` #gg
#for line in `jot - 6449 10700` #eg
#for line in `jot - 95100 99380` #ff
for line in `jot - 1 2196` #ff
#for line in `jot - 1 2485` #eg
#for line in `jot - 1 3614` #gg
#for line in `jot - 1 84395` #ee
#for line in `jot - 1 77315` #ee
do
#  run_evt_lumi=`head -n $line /Users/rachelyohay/RA3/data/1140pb-1_categorized.txt | tail -n 1 | sed -e "s%\([0-9]* [0-9]* [0-9]*\).*%\1%"`
#    run_evt_bin=`head -n $line /Users/rachelyohay/RA3/data/ee.txt | tail -n 1`
    run_evt_lumi=`head -n $line /Users/rachelyohay/RA3/data/ff_Brian.txt | tail -n 1`
#    found=`grep "$run_evt" /Users/rachelyohay/RA3/data/ee.txt`
#    bin=`echo $found | sed -e "s%[0-9]* [0-9]* \([0-9]\)%\1%"`
#    found=`grep "$run_evt_lumi" /Users/rachelyohay/RA3/data/ff_Brian.txt`
#    run_evt=`echo $run_evt_bin | sed -e "s%\([0-9]* [0-9]*\).*%\1%"`
#    bin=`echo $run_evt_bin | sed -e "s%[0-9]* [0-9]* \(.*\)%\1%"`
#    if [ "$bin" == "4" ]
#	then
#	found=`grep "$run_evt" /Users/rachelyohay/RA3/data/eeEvents_Data2011A_ToRun167913_Filter_NoPileupCorr_Photon.txt`
#  run_evt=`head -n $line /Users/rachelyohay/RA3/data/gg_Dave.txt | tail -n 1 | sed -e "s%\([0-9]* [0-9]*\).*%\1%"`
  found=`head -n 99380 /Users/rachelyohay/RA3/data/1140pb-1_categorized.txt | tail -n 4281 | grep "$run_evt_lumi"`
#    if [ "$found" != "$run_evt_lumi" ]
    if [ "$found" != "$run_evt_lumi" ] && [ "$found" != "" ]
#    if [ "$found" == "" ]
    then
	count=`expr $count + 1`
	echo $run_evt_lumi
    fi
#    fi
done
echo "$count events found in Brian's ff sample but not Rachel's"

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

#evts[0]="166841 483584716"
#evts[1]="166864 394801663"
#evts[2]="166864 416227374"
#evts[3]="163817 54364321"
#evts[4]="163589 54614720"
#evts[5]="163596 3521103"
#evts[6]="163660 18512069"
#evts[7]="163664 11099963"
#evts[8]="163758 340434605"
#evts[9]="163817 39521312"
#evts[10]="165567 660879884"
#evts[11]="165993 78420443"
#evts[12]="166380 380536201"
#evts[13]="166380 1455639992"
#evts[14]="166408 645967293"
#evts[15]="166554 681773881"
#evts[16]="166784 136897277"
#evts[17]="166787 9069334"
#evts[18]="166565 127439496"
#evts[19]="166701 69918333"
#evts[20]="167039 146433156"
#evts[21]="167041 120702599"
#evts[22]="166946 200822952"
#evts[23]="166950 587885196"
#evts[24]="167284 429291704"
#evts[25]="167282 22617452"
#evts[26]="167898 485562364"
#evts[27]="167675 622183354"
#evts[28]="167675 646543455"
#evts[29]="167746 90136125"
#evts[30]="167786 45391607"
#evts[31]="167786 50156901"
#evts[32]="167786 54739364"
#for iEvt in `jot - 0 32`
#do
#    evt=${evts[$iEvt]}
#    bin=`grep "$evt" /Users/rachelyohay/RA3/data/ee.txt`
#    echo $bin
#done

exit 0
