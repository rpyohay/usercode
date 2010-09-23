#!/bin/bash

#get cfg file name
if [ -z $1 ]
    then
    echo "Usage: ./makeAnalysisMacro.sh <cfg_file_name> <events_per_call> <number_of_calls>"
    exit 0
fi
cfg=$1

#get events per call
if [ -z $2 ]
    then
    echo "Usage: ./makeAnalysisMacro.sh <cfg_file_name> <events_per_call> <number_of_calls>"
    exit 0
fi
evts_per_call=$2

#get number of calls
if [ -z $3 ]
    then
    echo "Usage: ./makeAnalysisMacro.sh <cfg_file_name> <events_per_call> <number_of_calls>"
    exit 0
fi
number_of_calls=$3

#parse cfg file
parse_cfg_file () {
    return_value=`cat $cfg | grep $1 | grep -v \# | sed -e "s%$1 \(.*\)%\1%"`
}
parse_cfg_file "macro"
macro_name=$return_value
parse_cfg_file "ntuple"
ntuples=$return_value
parse_cfg_file "HLT_bit"
HLT_bits=$return_value
parse_cfg_file "output"
output_name=$return_value
parse_cfg_file "cfg"
macro_cfg=$return_value

#start writing macro file
if [ -e $macro_name ]
    then
    echo "Macro $macro_name already exists."
    exit 0
fi
cat >> $macro_name <<EOF
{
  gROOT->Reset();
  gROOT->LoadMacro("/afs/cern.ch/user/y/yohay/UserCode/yohay/GMSBAnalysis/EventSelector.cc++");
  gROOT->LoadMacro("/afs/cern.ch/user/y/yohay/UserCode/yohay/GMSBAnalysis/analysis.C++");

  //input files holding trees
  vector<string> fileList;
EOF

#write file names to macro file
for ntuple_name in $ntuples
  do
  if [ -z `echo $ntuple_name | grep castor` ]
      then
      stored_on_CASTOR=0
      file_list=$ntuple_name
  else
      stored_on_CASTOR=1
      file_list=`rfdir $ntuple_name | awk '{ print $NF }'`
  fi
  for file in $file_list
    do
    if [ $stored_on_CASTOR -eq 1 ]
	then
	if [ -z `echo $ntuple_name | grep root` ]
	    then
	    formatted_file_name_with_path="rfio:$ntuple_name/$file"
	else
	    formatted_file_name_with_path="rfio:$file"
	fi
    else
	formatted_file_name_with_path=$file
    fi
    cat >> $macro_name <<EOF
  fileList.push_back("$formatted_file_name_with_path");
EOF
  done
done

#write rest of macro file
cat >> $macro_name <<EOF

  //desired HLT selection
  vector<int> HLTBits;
EOF
for HLT_bit in $HLT_bits
  do
  cat >> $macro_name <<EOF
  HLTBits.push_back($HLT_bit);
EOF
done
cat >> $macro_name <<EOF

  //run!
EOF
for call in `seq $number_of_calls`
  do
  call_minus_1=`expr $call - 1`
  let min_evt=$(($evts_per_call*$call_minus_1))
  max_evt=`expr $min_evt + $evts_per_call - 1`
  cat >> $macro_name <<EOF
  runSampleMaker("$output_name", fileList, "$macro_cfg", HLTBits, $min_evt, $max_evt);
EOF
done
cat >> $macro_name <<EOF
}
EOF

exit 0