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
#include <iostream>

//the tree has 2 branches: an object of class susy::Event and an object of class susy::Category
#include "../src/SusyEvent.h"
#include "../src/SusyCategory.h"
#include "../../../PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"

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
   void setNEvtsProcessed(const vector<unsigned int>&);
   void setWeight(const vector<float>&);
   void setIntLumi(const float);
   void setMCPU(const Double_t*, const unsigned int);
   void setDataPU(TFile&);
   void setPUFile(const string&);
   void setL1JECFile(const string&);
   void setL2JECFile(const string&);
   void setL3JECFile(const string&);

   TString getTag() const;
   int getNEvts() const;
   string getProcess(const unsigned int) const;
   float getXSec(const string&) const;
   unsigned int getNEvtsProcessed(const string&) const;
   float getWeight(const string&) const;
   unsigned int getSize() const;
   float getIntLumi() const;
   double getMCPU(const unsigned int) const;
   double getDataPU(const unsigned int) const;
   string getPUFile() const;
   string getL1JECFile() const;
   string getL2JECFile() const;
   string getL3JECFile() const;

   bool matchedToGenParticle(const susy::Photon&, const VINT&, const UChar_t) const;
   bool passesDenominatorSelection(const unsigned int, const susy::Photon&, const VINT&) const;
   bool passesNumeratorSelection(const unsigned int) const;
   unsigned int numGoodVertices() const;
   void countEE(string&);
   template<typename T>
     void setHistogramOptions(T& hist, const string& xAxisTitle, const string& yAxisTitle, 
			      const string& zAxisTitle, const bool sumW2 = false) const
     {
       hist.GetXaxis()->SetTitle(xAxisTitle.c_str());
       hist.GetYaxis()->SetTitle(yAxisTitle.c_str());
       hist.GetZaxis()->SetTitle(zAxisTitle.c_str());
       if (sumW2) hist.Sumw2();
     }
   void setCanvasOptions(TCanvas&, const string&, const string&, const float, const float, 
			 const float, const float, const bool setGrid = false) const;
   float fillWeightsHistograms(const TH1D*, TH1D*, TH1F&, TH1F&) const;
   string histName(const string&, const string&, const unsigned int) const;
   void generateBackgroundSubtractedSpectra(TH3F&, TH2F&) const;
   float normAndErrorSquared(const TH3F&, const TH1F&, const TH1D*, const unsigned int, 
			     float&) const;
   void bookToyDiEMETWeightsHistograms(const vector<TH1F*>, const string&, vector<TH1F*>&, 
				       const Int_t nMETBins = 1) const;
   void makeToyDiEMETWeightsHistograms(TRandom3&, TH1F&, const TH1F&, vector<TH1F*>&, 
				       const unsigned int iMETBin = 1) const;
   void reweightDefault(const VFLOAT&, const VFLOAT&, const TH1F&, TH1F*) const;
   void reweightBinned(const TH2F&, const TH1F&, TH1F*, const unsigned int) const;
   void generateToys(vector<TH1F*>&, vector<TH1F*>&, const vector<TH1F*>, const unsigned int, 
		     const string&, const Double_t*, const unsigned int, const Double_t*, 
		     const VFLOAT&, const VFLOAT&) const;
   void generateToys(vector<TH1F*>&, vector<TH1F*>&, const vector<TH1F*>, const unsigned int, 
		     const string&, const Double_t*, const unsigned int, const Double_t*, 
		     const TH2F&) const;
   void fillToyDistributions(vector<TH1F*>&, const TH1F&, TCanvas&, vector<TH1F*>&, 
			     const unsigned int, const string&, const unsigned int) const;
   void setMETErrorBars(TH1F&, const vector<TH1F*>&, const vector<TH1F*>&, const vector<TH1F*>&, 
			const float, const float, const bool doDefault = true) const;
   void makeFinalCanvas(TH1*, const Color_t, const Width_t, const Style_t, const Color_t, 
			const Size_t, const string&) const;
   void deallocateMemory(VTH1F&) const;
   void runMETAnalysis(const string&);
   void runMETAnalysisWithEEBackgroundFit(const std::string&);
   void testFitting(const string&, const string&) const;
   void runEEVsFFAnalysis(const std::string&);
   void runCutFlowAnalysis();
   void skim(const string&, const int);
   void stripBranch(const string&, const string&);
   void debugPrint(const unsigned int) const;
   void reset();
   void initPU();
   void clearPU();
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
 private:
   TString tag_;                           //photon collection tag
   int     nEvts_;                         //number of events to process
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
   vector<float> MCPU_;                    //MC PU distribution
   vector<float> dataPU_;                  //data PU distribution
   unsigned int PUSize_;                   //number of PU bins
   string PUFile_;                         //name of data PU file
   reweight::LumiReWeighting lumiWeights_; //PU reweighting object
   string L1JECFile_;                      //file name for L1 JEC
   string L2JECFile_;                      //file name for L2 JEC
   string L3JECFile_;                      //file name for L3 JEC
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
   nEvts_ = 0;
   intLumi_ = 0.0;
   reset();
   clearPU();
}

GMSBAnalyzer::~GMSBAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();

   //clear private data members
   tag_ = "";
   nEvts_ = 0;
   intLumi_ = 0.0;
   PUFile_ = "";
   L1JECFile_ = "";
   L2JECFile_ = "";
   L3JECFile_ = "";
   reset();
   clearPU();
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
  nEvtsProcessed_.push_back(nEvts);
  weight_.push_back(xSec*intLumi_/nEvts);
  ++size_;
}
void GMSBAnalyzer::setFileMapEntry(const string& processName, const float xSec, 
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
void GMSBAnalyzer::setFileMap(const map<string, unsigned int>& fileMap) { fileMap_ = fileMap; }
void GMSBAnalyzer::setXSec(const vector<float>& xSec) { xSec_ = xSec; }
void GMSBAnalyzer::setNEvtsProcessed(const vector<unsigned int>& nEvts)
{
  nEvtsProcessed_ = nEvts;
}
void GMSBAnalyzer::setWeight(const vector<float>& weight) { weight_ = weight; }
void GMSBAnalyzer::setIntLumi(const float intLumi) { intLumi_ = intLumi; }
void GMSBAnalyzer::setMCPU(const Double_t* PUArray, const unsigned int arraySize)
{
  for (unsigned int iPV = 0; iPV < arraySize; ++iPV) { MCPU_.push_back(PUArray[iPV]); }
}
void GMSBAnalyzer::setDataPU(TFile& PUFile)
{
  TH1D* PUDist = (TH1D*)PUFile.Get("pileup");
  for (Int_t iPV = 1; iPV <= PUDist->GetNbinsX(); ++iPV) {
    dataPU_.push_back(PUDist->GetBinContent(iPV));
  }
}
void GMSBAnalyzer::setPUFile(const string& PUFile) { PUFile_ = PUFile; }
void GMSBAnalyzer::setL1JECFile(const string& L1JECFile) { L1JECFile_ = L1JECFile; }
void GMSBAnalyzer::setL2JECFile(const string& L2JECFile) { L2JECFile_ = L2JECFile; }
void GMSBAnalyzer::setL3JECFile(const string& L3JECFile) { L3JECFile_ = L3JECFile; }

TString GMSBAnalyzer::getTag() const { return tag_; }
int GMSBAnalyzer::getNEvts() const { return nEvts_; }
string GMSBAnalyzer::getProcess(const unsigned int index) const
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
float GMSBAnalyzer::getXSec(const string& processName) const
{
  return getValue<float>(processName, xSec_);
}
unsigned int GMSBAnalyzer::getNEvtsProcessed(const string& processName) const
{
  return getValue<unsigned int>(processName, nEvtsProcessed_);
}
float GMSBAnalyzer::getWeight(const string& processName) const
{
  return getValue<float>(processName, weight_);
}
unsigned int GMSBAnalyzer::getSize() const { return size_; }
float GMSBAnalyzer::getIntLumi() const { return intLumi_; }
double GMSBAnalyzer::getMCPU(const unsigned int iPV) const
{
  if (iPV < PUSize_) return MCPU_[iPV];
  else return -1.0;
}
double GMSBAnalyzer::getDataPU(const unsigned int iPV) const
{
  if (iPV < PUSize_) return dataPU_[iPV];
  else return -1.0;
}
string GMSBAnalyzer::getPUFile() const { return PUFile_; }
string GMSBAnalyzer::getL1JECFile() const { return L1JECFile_; }
string GMSBAnalyzer::getL2JECFile() const { return L2JECFile_; }
string GMSBAnalyzer::getL3JECFile() const { return L3JECFile_; }

void GMSBAnalyzer::reset()
{
  fileMap_.clear();
  xSec_.clear();
  nEvtsProcessed_.clear();
  weight_.clear();
  size_ = 0;
}
void GMSBAnalyzer::initPU()
{
  const Double_t PoissonOneXDist[25] = {
     0.14551,
     0.0644453,
     0.0696412,
     0.0700311,
     0.0694257,
     0.0685655,
     0.0670929,
     0.0646049,
     0.0609383,
     0.0564597,
     0.0508014,
     0.0445226,
     0.0378796,
     0.0314746,
     0.0254139,
     0.0200091,
     0.0154191,
     0.0116242,
     0.00846857,
     0.00614328,
     0.00426355,
     0.00300632,
     0.00203485,
     0.00133045,
     0.000893794
  };
  setMCPU(PoissonOneXDist, 25);
  TFile PUFile(PUFile_.c_str());
  setDataPU(PUFile);
  PUFile.Close();
  const unsigned int MCPUSize = MCPU_.size();
  const unsigned int dataPUSize = dataPU_.size();
  cout << "MCPU_.size() = " << MCPU_.size() << endl;
  cout << "dataPU_.size() = " << dataPU_.size() << endl;
  if (MCPUSize < dataPUSize) dataPU_.erase(dataPU_.begin() + MCPUSize, dataPU_.end());
  else if (MCPUSize > dataPUSize) MCPU_.erase(MCPU_.begin() + dataPUSize, MCPU_.end());
  cout << "MCPU_.size() = " << MCPU_.size() << endl;
  cout << "dataPU_.size() = " << dataPU_.size() << endl;
  PUSize_ = dataPU_.size();
  lumiWeights_ = reweight::LumiReWeighting(MCPU_, dataPU_);
}
void GMSBAnalyzer::clearPU()
{
  MCPU_.clear();
  dataPU_.clear();
  PUSize_ = 0;
  lumiWeights_ = reweight::LumiReWeighting();
}
#endif // #ifdef GMSBAnalyzer_cxx
