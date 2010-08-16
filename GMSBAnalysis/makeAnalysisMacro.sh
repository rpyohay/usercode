#!/bin/bash

#get cfg file name
if [ -z $1 ]
    then
    echo "Usage: ./makeAnalysisMacro.sh <cfg_file_name>"
    exit 0
fi
cfg=$1

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
  for file in `rfdir $ntuple_name | awk '{ print $NF }'`
    do
    formatted_file_name_with_path="rfio:$ntuple_name/$file"
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
  runSampleMaker("$output_name", fileList, "$macro_cfg", HLTBits);
}
EOF

exit 0