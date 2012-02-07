// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.h
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.h,v 1.4 2011/06/08 16:28:40 dmason Exp $
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <map>

#include "../src/SusyEvent.h"


class SusyEventAnalyzer {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  susy::Event     *event;

  // List of branches
  TBranch        *b_Event;

  SusyEventAnalyzer(TTree *tree=0);
  virtual ~SusyEventAnalyzer();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();                          // event loop for main analysis

  // utility functions
  bool isSameObject(TLorentzVector& p1, TLorentzVector& p2);
  float d0correction(TVector3& beamSpot, susy::Track& track) const;
  void IncludeAJson(std::string jsonfile);  // Call to pull in a json file 
  bool isInJson(Int_t run,Int_t lumi);      // JSON based good run list cut...
  bool PassTrigger(TString v); // return true if path v is fired
  bool PassTriggers(); // return true if any of names in hltNames are fired

  // parameter configuration functions
  void Initialize();         // global variables needed to be initialized just once
  void InitializePerEvent(); // global variables needed to be initialized per event
  void SetDataset(TString& v) {          ds = v; }
  void SetPrintInterval(int v) {         printInterval = v; }
  void SetPrintLevel(int v) {            printLevel = v; }
  void SetProcessNEvents(int v) {        processNEvents = v; }
  void SetUseTrigger(bool v) {           useTrigger = v; }
  void AddHltName(TString v) {           hltNames.push_back(v); }
  void SetFilter(bool v) {               enableFilter = v; }
  void SetFilteredFileName(TString v) {  filtered_file_name = v; }

  //functions added by Rachel
  template<typename T>
  void setHistogramOptions(T& hist, const string& xAxisTitle, const string& yAxisTitle, 
			   const string& zAxisTitle, const bool sumW2 = false) const
  {
    hist.GetXaxis()->SetTitle(xAxisTitle.c_str());
    hist.GetYaxis()->SetTitle(yAxisTitle.c_str());
    hist.GetZaxis()->SetTitle(zAxisTitle.c_str());
    if (sumW2) hist.Sumw2();
  }
  template<typename T, typename U>
  bool indexOutOfRange(const T& iProcess, const U* vec) const
  {
    const unsigned int vecSize = vec->size();
    if (iProcess->second >= vecSize) {
      cerr << "Error: vector index " << iProcess->second << " corresponding to process ";
      cerr << iProcess->first << " is larger than max vector index " << (vecSize - 1) << ".\n";
      return true;
    }
    else return false;
  }
  template<typename T>
  T getValue(const string& processName, const vector<T>& val) const
  {
    T value = 0;
    map<string, unsigned int>::const_iterator iProcess = fileMap_.find(processName);
    if (iProcess != fileMap_.end()) {
      if (!indexOutOfRange(iProcess, &val)) value = val[iProcess->second];
    }
    return value;
  }
  void addEntry(const string&, const float, const unsigned int, const unsigned int);
  void setFileMapEntry(const string&, const float, const unsigned int);
  void setFileMap(const map<string, unsigned int>&);
  void setXSec(const vector<float>&);
  void setNEvtsProcessed(const vector<unsigned int>&);
  void setWeight(const vector<float>&);
  void setIntLumi(const float);
  string getProcess(const unsigned int) const;
  float getXSec(const string&) const;
  unsigned int getNEvtsProcessed(const string&) const;
  float getWeight(const string&) const;
  unsigned int getSize() const;
  float getIntLumi() const;
  void reset();
  void plot(const string&);

 private:

  TString ds;               // dataset name to be used for output histfile name

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

  typedef std::map<int,std::map<int,bool> > RunLumiFlagHolder;  //define map that holds json list
  RunLumiFlagHolder goodrunlumilist;  // instantiate it

  //added by Rachel
  map<string, unsigned int> fileMap_;     /*key is the MC process name (which should be found 
					    within the absolute file name in the TChain), mapped 
					    value is the index of the file in vectors of cross 
					    sections, number of events processed, and 
					    corresponding event weights*/
  vector<float> xSec_;                    //MC process cross section, 1 element per process
  vector<unsigned int> nEvtsProcessed_;   /*no. events processed per MC dataset, 1 element per 
					    dataset*/
  vector<float> weight_;                  //MC process weight, 1 element per process
  float intLumi_;                         //integrated luminosity to normalize events to
  unsigned int size_;                     //size of map/vectors
};

#endif

#ifdef SusyEventAnalyzer_cxx
SusyEventAnalyzer::SusyEventAnalyzer(TTree *tree)
{
  if (tree == 0) {
    std::cout << "Error!!! There is no file containing a tree." << std::endl;
  }
  Init(tree);
  Initialize();
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();

  //added by Rachel
  intLumi_ = 0.0;
  reset();
}

Int_t SusyEventAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t SusyEventAnalyzer::LoadTree(Long64_t entry)
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

void SusyEventAnalyzer::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  //   fChain->SetMakeClass(1);

  event = new susy::Event;

  fChain->SetBranchAddress("susyEvent", &event, &b_Event);

  //added by Rachel
  intLumi_ = 0.0;
  reset();
}

void SusyEventAnalyzer::Initialize() {

  ds = "test";
  printLevel = 0;
  printInterval = 1000;
  processNEvents = -1;
  useTrigger = false;
  enableFilter = false;
  filtered_file_name = "filtered.root";

}

void SusyEventAnalyzer::IncludeAJson(std::string jsonfile) {


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


bool SusyEventAnalyzer::isInJson(Int_t run,Int_t lumi) {

//#ifdef MC
//  return 1;
//#endif

if (goodrunlumilist[run][lumi]) return true;

return false;

}

void SusyEventAnalyzer::addEntry(const string& processName, const float xSec, 
				 const unsigned int nEvts, const unsigned int index)
{
  fileMap_[processName] = index;
  xSec_.push_back(xSec);
  nEvtsProcessed_.push_back(nEvts);
  weight_.push_back(xSec*intLumi_/nEvts);
  ++size_;
}

void SusyEventAnalyzer::setFileMapEntry(const string& processName, const float xSec, 
					const unsigned int nEvts)
{
  //check that the file map and all weight-related vectors are of equal size
  const unsigned int fileMapSize = fileMap_.size();
  const unsigned int xSecSize = xSec_.size();
  const unsigned int nEvtsSize = nEvtsProcessed_.size();
  const unsigned int weightSize = weight_.size();
  if ((fileMapSize != xSecSize) || (xSecSize != nEvtsSize) || (nEvtsSize != weightSize) || 
      (weightSize != size_)) {
    cerr << "Error: size mismatch.\n";
    cerr << "fileMapSize = " << fileMapSize << endl;
    cerr << "xSecSize = " << xSecSize << endl;
    cerr << "nEvtsSize = " << nEvtsSize << endl;
    cerr << "weightSize = " << weightSize << endl;
    cerr << "size_ = " << size_ << endl;
    cerr << "Clearing everything and starting over.\n";
    reset();
  }

  //check that nEvts != 0
  if (nEvts == 0) {
    cerr << "Error: nEvts = 0.  Event weight is meaningless.  Exiting.\n";
    return;
  }

  //if the process name exists in the map, change its attributes
  map<string, unsigned int>::iterator iProcess = fileMap_.find(processName);
  if (iProcess != fileMap_.end()) {

    //check that the mapped value (the vector index) doesn't refer to a nonexistent vector element
    if (indexOutOfRange(iProcess, &xSec_)) {
      cerr << "Clearing everything and starting over.\n";
      reset();
      addEntry(processName, xSec, nEvts, 0);
    }
    else {

      //update the map
      xSec_[iProcess->second] = xSec;
      nEvtsProcessed_[iProcess->second] = nEvts;
      weight_[iProcess->second] = xSec*intLumi_/nEvts;
    }
  }
  else {

    //add this process to the map
    addEntry(processName, xSec, nEvts, xSecSize);
  }
}

void SusyEventAnalyzer::setFileMap(const map<string, 
				   unsigned int>& fileMap) { fileMap_ = fileMap; }

void SusyEventAnalyzer::setXSec(const vector<float>& xSec) { xSec_ = xSec; }

void SusyEventAnalyzer::setNEvtsProcessed(const vector<unsigned int>& nEvts)
{
  nEvtsProcessed_ = nEvts;
}

void SusyEventAnalyzer::setWeight(const vector<float>& weight) { weight_ = weight; }

void SusyEventAnalyzer::setIntLumi(const float intLumi) { intLumi_ = intLumi; }

string SusyEventAnalyzer::getProcess(const unsigned int index) const
{
  bool found = false;
  string process = "";
  map<string, unsigned int>::const_iterator iMap = fileMap_.begin();
  while ((iMap != fileMap_.end()) && !found) {
    if (iMap->second == index) {
      found = true;
      process = iMap->first;
    }
    ++iMap;
  }
  return process;
}

float SusyEventAnalyzer::getXSec(const string& processName) const
{
  return getValue<float>(processName, xSec_);
}

unsigned int SusyEventAnalyzer::getNEvtsProcessed(const string& processName) const
{
  return getValue<unsigned int>(processName, nEvtsProcessed_);
}

float SusyEventAnalyzer::getWeight(const string& processName) const
{
  return getValue<float>(processName, weight_);
}

unsigned int SusyEventAnalyzer::getSize() const { return size_; }

float SusyEventAnalyzer::getIntLumi() const { return intLumi_; }

void SusyEventAnalyzer::reset()
{
  fileMap_.clear();
  xSec_.clear();
  nEvtsProcessed_.clear();
  weight_.clear();
  size_ = 0;
}

#endif // #ifdef SusyEventAnalyzer_cxx
