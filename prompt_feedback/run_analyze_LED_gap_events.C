{
#include <vector>
#include <string>
#include <ctime>

gROOT->Reset();
gROOT->LoadMacro("/afs/cern.ch/user/y/yohay/LED/CMSSW_3_2_1/src/prompt_feedback/analyze_LED_gap_events.C++");

vector<string> sources;
sources.push_back("rfio:/castor/cern.ch/user/y/yohay/run_118358_trees/tree_118358-ped_1st_sample-events_1-5000.root");
sources.push_back("rfio:/castor/cern.ch/user/y/yohay/run_118358_trees/tree_118358-ped_1st_sample-events_5001-10000.root");
sources.push_back("rfio:/castor/cern.ch/user/y/yohay/run_118358_trees/tree_118358-ped_1st_sample-events_10001-15000.root");

struct tm start;
start.tm_isdst = 0; //not Daylight Savings Time
start.tm_yday = 298; //Oct. 26
start.tm_wday = 1; //Mon.
start.tm_year = 109; //2009
start.tm_mon = 9; //Oct.
start.tm_mday = 26; //26
start.tm_hour = 23; //23:
start.tm_min = 04; //04:
start.tm_sec = 36; //36
time_t rawStart = mktime(&start);

analyzeLEDGapEvents(sources, "/afs/cern.ch/user/y/yohay/scratch0/analysis_118358-ped_1st_sample.root", rawStart);
}
