// -*- C++ -*-
//
// Package:    GMSBTools
// Class:      EventAnalyzer
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
//

#ifndef EventAnalyzer_h
#define EventAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <map>

#include "../src/SusyEvent.h"
#include "../../../GMSBTools/Filters/interface/Categorizer.h"
#include "../src/SusyCategory.h"

//isolation cone sizes
#define DR03 0
#define DR04 1

class EventAnalyzer {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  susy::Event     *event;

  // List of branches
  TBranch        *b_Event;

  //events passing particular selections
  RUNEVTLUMIMAP ggEvts_;
  RUNEVTLUMIMAP egEvts_;
  RUNEVTLUMIMAP eeEvts_;
  RUNEVTLUMIMAP ffEvts_;

  //PU subtraction effective areas
  double photonECALIsoEffArea_;
  double photonHCALIsoEffArea_;

  //isolation cone size
  unsigned int isoConeHLT_;
  unsigned int isoConeOffline_;

  EventAnalyzer(TTree *tree=0);
  virtual ~EventAnalyzer();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop(TTree*, Categorizer, susy::Category*);  //event loop for main analysis
  virtual void     countTriggers(const std::string&);           //event loop for counting

  // utility functions
  bool isSameObject(TLorentzVector& p1, TLorentzVector& p2);
  float d0correction(TVector3& beamSpot, susy::Track& track) const;
  void IncludeAJson(std::string jsonfile);  // Call to pull in a json file 
  bool isInJson(Int_t run,Int_t lumi);      // JSON based good run list cut...
  bool PassTrigger(TString v); // return true if path v is fired
  bool PassTriggers(); // return true if any of names in hltNames are fired
  std::string PassLowestUnprescaledTrigger() const; /*return name of the lowest ET threshold 
						      (second leg) unprescaled trigger in hltNames 
						      fired by this event (empty string if none 
						      exists)*/
  std::string PassUnprescaledTrigger() const; //hack until ntuple prescales are fixed
  
  // parameter configuration functions
  void Initialize();         // global variables needed to be initialized just once
  void InitializePerEvent(); // global variables needed to be initialized per event
  void SetPrintInterval(int v) {          printInterval = v; }
  void SetPrintLevel(int v) {             printLevel = v; }
  void SetProcessNEvents(int v) {         processNEvents = v; }
  void SetUseTrigger(bool v) {            useTrigger = v; }
  void AddHltName(TString v) {            hltNames.push_back(v); }
  void SetFilter(bool v) {                enableFilter = v; }
  void SetFilteredFileName(TString v) {   filtered_file_name = v; }
  void SetPhotonTag(const TString& tag) { tag_ = tag; }
  void SetPhotonECALIsoEffArea(const double photonECALIsoEffArea)
  {
    photonECALIsoEffArea_ = photonECALIsoEffArea;
  }
  void SetPhotonHCALIsoEffArea(const double photonHCALIsoEffArea)
  {
    photonHCALIsoEffArea_ = photonHCALIsoEffArea;
  }
  void SetIsoConeHLT(const unsigned int isoConeHLT) { isoConeHLT_ = isoConeHLT; }
  void SetIsoConeOffline(const unsigned int isoConeOffline) { isoConeOffline_ = isoConeOffline; }
  void debugPrint(const Categorizer&, susy::Event*) const;

 private:

  // printLevel
  // 0 : default - no printout
  // 1 : print functional step in every event
  // 2 : print values in collections
  int printLevel;           // print frequency

  int printInterval;        // print level for event content: defined in Event.h
  int processNEvents;       // number of events to be processed
  bool useTrigger;          // flag for using trigger bit selection.
  std::vector<TString> hltNames;          // HLT trigger path names
  bool enableFilter;        // filter events of interest
  TString filtered_file_name; // filtered output file name
  TString tag_;             // name of photon collection tag

  typedef std::map<int,std::map<int,bool> > RunLumiFlagHolder;  //define map that holds json list
  RunLumiFlagHolder goodrunlumilist;  // instantiate it

};

#endif

#ifdef EventAnalyzer_cxx
EventAnalyzer::EventAnalyzer(TTree *tree)
{
  if (tree == 0) {
    std::cout << "Error!!! There is no file containing a tree." << std::endl;
  }
  Init(tree);
  Initialize();
}

EventAnalyzer::~EventAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t EventAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t EventAnalyzer::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
  }
  return centry;
}

void EventAnalyzer::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  //   fChain->SetMakeClass(1);

  event = new susy::Event;

  fChain->SetBranchAddress("susyEvent", &event, &b_Event);
}

void EventAnalyzer::Initialize() {

  printLevel = 0;
  printInterval = 1000;
  processNEvents = -1;
  useTrigger = false;
  enableFilter = false;
  filtered_file_name = "filtered.root";
  isoConeHLT_ = DR03;
  isoConeOffline_ = DR03;
}

void EventAnalyzer::IncludeAJson(std::string jsonfile) {


// Fairly primitive brute force json parser -- opens the json file named in the argument
// and adds that to the goodrunlumilist map.  Overlapping jsons are merged inclusively.

 char thing;

 ifstream jsonInput;

 std::cout << "Sucking in Json file: " << jsonfile << " which includes: " << std::endl;

 jsonInput.open(jsonfile.c_str());
 
 if (!jsonInput.good()) {
   std::cout << "Problem reading Json file...  Didn't suck anything in... " << std::endl;
   return;
 }
 
 jsonInput.get(thing);
 
 while (jsonInput.good()) {
   if (thing=='{') {  // start of list
     while (thing != '}') {
       int runnum;
       if (thing == '"') {
         std::string srunnum;
         jsonInput.get(thing); // get stuff inside ""

         while (thing != '"') {
           srunnum+=thing; // get stuff inside ""
           jsonInput.get(thing);

	   }
         sscanf(srunnum.c_str(),"%i",&runnum);
         std::cout << " runnum: " << runnum << std::endl;
         //bool newrun=true; //never gets used
         
       } // inside ""
       if (thing == '[') {
          jsonInput.get(thing); // get stuff inside []
	 while (thing != ']') {
           if (thing == '[') {
             jsonInput.get(thing); // get stuff inside series []

             std::string lumiseries;
             int firstlumi,lastlumi;
             while (thing !=']') {
               lumiseries+=thing;
                jsonInput.get(thing); // get stuff inside series []
             }
             sscanf(lumiseries.c_str(),"%i,%i",&firstlumi,&lastlumi);
             std::cout << "  lumis  " << firstlumi << " to " << lastlumi << std::endl;

	     // At this point have runnum, first lumi, last lumi -- so can fill map here...
	     for (int l=firstlumi;l<=lastlumi;l++) {
               goodrunlumilist[runnum][l]=true;
	     }

           } // inside actual series []
             jsonInput.get(thing); // get stuff inside []
         }
       } // inside []
         jsonInput.get(thing); // get another char looking for "

     } 
   } // inside {}
    jsonInput.get(thing); // get another char looking for {

 } // EOF 

 jsonInput.close();

}


bool EventAnalyzer::isInJson(Int_t run,Int_t lumi) {

  //exit true if no JSON was applied
  if (goodrunlumilist.size() == 0) return true;

  //use the JSON to determine the passing condition
  RunLumiFlagHolder::const_iterator iRunEvtPair = goodrunlumilist.find(run);
  if (iRunEvtPair != goodrunlumilist.end()) {
    std::map<int,bool> evtFlagMap = iRunEvtPair->second;
    std::map<int,bool>::const_iterator iRunEvtFlagTriplet = evtFlagMap.find(lumi);
    if (iRunEvtFlagTriplet != evtFlagMap.end()) {
      if (iRunEvtFlagTriplet->second) return true;
    }
  }

  //if the event is not in the JSON, exit false
  return false;

}

#endif // #ifdef EventAnalyzer_cxx
