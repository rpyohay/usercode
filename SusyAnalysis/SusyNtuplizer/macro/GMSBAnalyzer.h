//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 12 17:13:41 2011 by ROOT version 5.27/06b
// from TTree susyTree/SUSY Event
// found on file: /data2/yohay/categorized.root
//////////////////////////////////////////////////////////

#ifndef GMSBAnalyzer_h
#define GMSBAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//the tree has 2 branches: an object of class susy::Event and an object of class susy::Category
#include "../src/SusyEvent.h"
#include "../src/SusyCategory.h"

class GMSBAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   susy::Event     *susyEvent;
   susy::Category  *susyCategory;

   // List of branches
   TBranch        *b_susyEvent_;   //!
   TBranch        *b_susyCategory_;   //!

   GMSBAnalyzer(TTree *tree=0);
   virtual ~GMSBAnalyzer();
   virtual Int_t    Cut(/*Long64_t entry*/);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const std::string&);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void setTag(const TString&);
   void setNEvts(const int);
   void addEntry(const string&, const float, const unsigned int, const unsigned int);
   void setFileMapEntry(const string&, const float, const unsigned int);
   void setFileMap(const map<string, unsigned int>&);
   void setXSec(const vector<float>&);
   void setNEvts(const vector<unsigned int>&);
   void setWeight(const vector<float>&);
   void setIntLumi(const float);

   TString getTag() const;
   int getNEvts() const;
   map<string, unsigned int>* getFileMap() const;
   float getXSec(const string&) const;
   vector<float>* getXSec() const;
   unsigned int getNEvts(const string&) const;
   vector<unsigned int>* getNEvts() const;
   float getWeight(const string&) const;
   vector<float>* getWeight() const;
   float getIntLumi() const;

   bool matchedToGenParticle(const susy::Photon&, const VINT&, const UChar_t) const;
   bool passesDenominatorSelection(const unsigned int, const susy::Photon&, const VINT&) const;
   bool passesNumeratorSelection(const unsigned int) const;
   unsigned int numGoodVertices() const;
   void countEE(string&);
   template<typename T>
     void setHistogramOptions(T& hist, const string& xAxisTitle, const string& yAxisTitle, 
			      const bool sumW2 = false) const
     {
       hist.GetXaxis()->SetTitle(xAxisTitle.c_str());
       hist.GetYaxis()->SetTitle(yAxisTitle.c_str());
       if (sumW2) hist.Sumw2();
     }
   void setCanvasOptions(TCanvas&, const string&, const string&, const float, const float, 
			 const float, const float, const bool setGrid = false) const;
   float fillWeightsHistograms(TH2F&, TH2F&, TH1F&, TH1F&) const;
   float normAndErrorSquared(const TH2F&, const TH1F&, const TH1D*, const unsigned int, 
			     float&) const;
   void generateToys(vector<TH1F*>&, vector<TH1F*>&, const TH1F&, const unsigned int, 
		     const string&, const unsigned int, const Double_t*, const unsigned int, 
		     const Double_t*, const VFLOAT&, const VFLOAT&) const;
   void fillToyDistributions(vector<TH1F*>&, const TH1F&, TCanvas&, vector<TH1F*>&, 
			     const unsigned int, const string&, const unsigned int) const;
   void fillDifferenceHistograms(const TH1F&, TH1F&, const vector<TH1F*>&) const;
   void reweight(const VFLOAT&, const VFLOAT&, TH1F&, const TH1F&) const;
   void setMETErrorBars(TH1F&, const vector<TH1F*>&, const vector<TH1F*>&, const vector<TH1F*>&, 
			const float, const float) const;
   void makeFinalCanvas(TH1*, const Color_t, const Width_t, const Style_t, const Color_t, 
			const Size_t, const string&) const;
   void deallocateMemory(VTH1F&) const;
   void runMETAnalysis(const string&);
   void runEEVsFFAnalysis(const std::string&);
   void runCutFlowAnalysis();
   void skim(const string&, const int);
   void debugPrint(const unsigned int) const;
   void reset();
   template<typename T, typename U>
     bool GMSBAnalyzer::indexOutOfRange(const T& iProcess, const U* vec) const
   {
     const unsigned int vecSize = vec.size();
     if (index >= vecSize) {
       cerr << "Error: vector index " << index << " corresponding to process " << iProcess->first;
       cerr << " is larger than max vector index " << (vecSize - 1) << ".\n";
     }
   }
   template<typename T>
     T GMSBAnalyzer::getValue(const string& processName, const vector<T>& val) const
     {
       T value = 0;
       map<string, unsigned int>::const_iterator iProcess = fileMap_.find(processName);
       if (iProcess != fileMap_.end()) {
	 if (!indexOutOfRange(iProcess, &val_)) value = val_[iProcess->second];
       }
       return value;
     }
 private:
   TString tag_;                       //photon collection tag
   int     nEvts_;                     //number of events to process
   map<string, unsigned int> fileMap_; /*key is the MC process name (which should be found within 
					 the absolute file name in the TChain), mapped value is 
					 the index of the file in vectors of cross sections, 
					 number of events processed, and corresponding event 
					 weights*/
   vector<float> xSec_;                //MC process cross section, 1 element per process
   vector<unsigned int> nEvts_;        //no. events processed per MC dataset, 1 element per dataset
   vector<float> weight_;              //MC process weight, 1 element per process
   float intLumi_;                     //integrated luminosity to normalize events to
};

#endif

#ifdef GMSBAnalyzer_cxx
GMSBAnalyzer::GMSBAnalyzer(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data2/yohay/categorized.root");
      if (!f) {
         f = new TFile("/data2/yohay/categorized.root");
      }
      tree = (TTree*)gDirectory->Get("susyTree");

   }
   Init(tree);

   //initialize private data members
   tag_ = "";
   nEvts_ = 0;
   intLumi_ = 0.0;
   reset();
}

GMSBAnalyzer::~GMSBAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();

   //clear private data members
   tag_ = "";
   nEvts_ = 0;
   intLumi_ = 0.0;
   reset();
}

Int_t GMSBAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GMSBAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GMSBAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;

   susyEvent = new susy::Event;
   fChain->SetBranchAddress("susyEvent", &susyEvent, &b_susyEvent_);
   susyCategory = new susy::Category;
   fChain->SetBranchAddress("susyCategory", &susyCategory, &b_susyCategory_);
   Notify();
}

Bool_t GMSBAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GMSBAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GMSBAnalyzer::Cut(/*Long64_t entry*/)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void GMSBAnalyzer::setTag(const TString& tag) { tag_ = tag; }
void GMSBAnalyzer::setNEvts(const int nEvts) { nEvts_ = nEvts; }
void GMSBAnalyzer::addEntry(const string& processName, const float xSec, const unsigned int nEvts, 
			    const unsigned int index)
{
  fileMap_[processName] = index;
  xSec_.push_back(xSec);
  nEvts_.push_back(nEvts);
  weight_.push_back(xSec*intLumi/nEvts);
}
void GMSBAnalyzer::setFileMapEntry(const string& processName, const float xSec, 
				   const unsigned int nEvts)
{
  //check that the file map and all weight-related vectors are of equal size
  const unsigned int fileMapSize = fileMap_.size();
  const unsigned int xSecSize = xSec_.size();
  const unsigned int nEvtsSize = nEvts_.size();
  const unsigned int weightSize = weight_.size();
  if ((fileMapSize != xSecSize) || (xSecSize != nEvtsSize) || (nEvtsSize != weightSize)) {
    cerr << "Error: size mismatch.\n";
    cerr << "fileMapSize = " << fileMapSize << endl;
    cerr << "xSecSize = " << xSecSize << endl;
    cerr << "nEvtsSize = " << nEvtsSize << endl;
    cerr << "weightSize = " << weightSize << endl;
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
      xSec_[index] = xSec;
      nEvts_[index] = nEvts;
      weight_[index] = xSec*intLumi/nEvts;
    }
  }
  else {

    //add this process to the map
    addEntry(processName, xSec, nEvts, xSecSize);
  }
}
void GMSBAnalyzer::setFileMap(const map<string, unsigned int>& fileMap) { fileMap_ = fileMap; }
void GMSBAnalyzer::setXSec(const vector<float>& xSec) { xSec_ = xSec; }
void GMSBAnalyzer::setNEvts(const vector<unsigned int>& nEvts) { nEvts_ = nEvts; }
void GMSBAnalyzer::setWeight(const vector<float>& weight) { weight_ = weight; }
void GMSBAnalyzer::setIntLumi(const float intLumi) { intLumi_ = intLumi; }

TString GMSBAnalyzer::getTag() const { return tag_; }
int GMSBAnalyzer::getNEvts() const { return nEvts_; }
map<string, unsigned int>* GMSBAnalyzer::getFileMap() const { return &fileMap_; }
float GMSBAnalyzer::getXSec(const string& processName) const
{
  return getValue<float>(processName, xSec_);
}
vector<float>* GMSBAnalyzer::getXSec() const { return &xSec_; }
unsigned int GMSBAnalyzer::getNEvts(const string&) const
{
  return getValue<unsigned int>(processName, nEvts_);
}
vector<unsigned int>* GMSBAnalyzer::getNEvts() const { return &nEvts_; }
float GMSBAnalyzer::getWeight(const string&) const
{
  return getValue<float>(processName, weight_);
}
vector<float>* GMSBAnalyzer::getWeight() const { return &weight_; }
float GMSBAnalyzer::getIntLumi() const { return intLumi_; }

void GMSBAnalyzer::reset()
{
  fileMap_.clear();
  xSec_.clear();
  nEvts_.clear();
  weight_.clear();
}
#endif // #ifdef GMSBAnalyzer_cxx
