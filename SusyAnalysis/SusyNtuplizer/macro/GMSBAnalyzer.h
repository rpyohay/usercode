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
#include "../../../GMSBTools/Filters/interface/Categorizer.h"
#include "../../../PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"
#include "../../../DataFormats/Math/interface/deltaR.h"
#include "RooFitResult.h"
	 	 
#define LOWER_BOUND 0
#define UPPER_BOUND 1
#define VALUE 2
#define UNIT 3
	 	 
class minDR {
 public:
  void setPhoton(const susy::Photon* photon) { photon_ = const_cast<susy::Photon*>(photon); }
  void deletePhoton() { photon_ = NULL; }
  susy::Photon* getPhoton() const { return photon_; }
  bool operator()(const susy::PFParticle* PFCandidate1, const susy::PFParticle* PFCandidate2)
  {
    const float photonEta = photon_->momentum.Eta();
    const float photonPhi = photon_->momentum.Phi();
    return (deltaR(photonEta,
		   photonPhi,
		   PFCandidate1->momentum.Eta(),
		   PFCandidate1->momentum.Phi()) <
	    deltaR(photonEta,
		   photonPhi,
		   PFCandidate2->momentum.Eta(),
		   PFCandidate2->momentum.Phi()));
  }
 private:
  susy::Photon* photon_;
};

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
  void setJECErrFile(const string&);
  void addHLT(const TString&, const unsigned int, const unsigned int, const unsigned int);
  void setDiEMInvMassFitParLowerBound(const unsigned int, const string, const float,
				      const bool overwrite = false);
  void setDiEMInvMassFitParUpperBound(const unsigned int, const string, const float,
				      const bool overwrite = false);
  void setDiEMInvMassFitPar(const unsigned int, const string, const float,
			    const bool overwrite = false);
  void setDiEMInvMassFitParUnit(const string, const string, const bool overwrite = false);

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
  string getJECErrFile() const;
  const map<pair<TString, unsigned int>, pair<unsigned int, unsigned int> >* getHLT() const;
  unsigned int getMinRun(const TString&) const;
  unsigned int getMaxRun(const TString&) const;
  float getDiEMInvMassFitParLowerBound(const unsigned int, const string) const;
  float getDiEMInvMassFitParUpperBound(const unsigned int, const string) const;
  float getDiEMInvMassFitPar(const unsigned int, const string) const;
  string getDiEMInvMassFitParUnit(const string) const;

  bool matchedToGenParticle(const susy::Photon&, const VINT&, const UChar_t) const;
  bool passesDenominatorSelection(const unsigned int, const susy::Photon&, const VINT&) const;
  bool passesNumeratorSelection(const unsigned int) const;
  unsigned int numGoodVertices() const;
  void countEE(string&);
  void plotTemp(string&);
  template<typename T>
    void setHistogramOptions(T& hist, const string& xAxisTitle, const string& yAxisTitle, 
			     const string& zAxisTitle, const bool sumW2 = false) const
    {
      hist.GetXaxis()->SetTitle(xAxisTitle.c_str());
      hist.GetYaxis()->SetTitle(yAxisTitle.c_str());
      hist.GetZaxis()->SetTitle(zAxisTitle.c_str());
      if (sumW2) hist.Sumw2();
    }
  template<typename T>
    TLorentzVector correctedJet4Momentum(T& iJet, const string& corrLabel1, float& JES, 
					 const string& corrLabel2 = "") const
    {
      map<TString, Float_t>::const_iterator iCorr1 = iJet->jecScaleFactors.find(corrLabel1);
      map<TString, Float_t>::const_iterator iCorr2 = iJet->jecScaleFactors.find(corrLabel2);
      TLorentzVector corrP4;
      if (iCorr1 != iJet->jecScaleFactors.end()) {
	if (iCorr2 != iJet->jecScaleFactors.end()) {
	  corrP4 = (iCorr1->second/iCorr2->second)*iJet->momentum;
	  JES = iCorr1->second/iCorr2->second;
	}
	else {
	  corrP4 = iCorr1->second*iJet->momentum;
	  JES = iCorr1->second;
	}
      }
      return corrP4;
    }
  void setCanvasOptions(TCanvas&, const string&, const string&, const float, const float, 
			const float, const float, const bool setGrid = false) const;
  float fillWeightsHistograms(const TH1D*, TH1D*, TH1F&, TH1F&, const float) const;
  string histName(const string&, const string&, const unsigned int) const;
  void format(const char*, const char*, const char*, const char*, const char*) const;
  void format(const char*, const float, const float, const float) const;
  void generateBackgroundSubtractedSpectra(TH3F&, TH2F&) const;
  float normAndErrorSquared(const TH3F&, const TH1F&, const TH1D*, const vector<TH1F*>&,
			    const vector<TH1F*>&, const vector<TH1F*>&, const float, const float,
			    const unsigned int, float&, float&, float&, float&,
			    const bool doDefault = true) const;
  float razorNormAndErrorSquared(const TH2F&, const TH2F&, const TH2F&, const unsigned int, 
				 const unsigned int, float&) const;
  void bookToyDiEMETWeightsHistograms(const vector<TH1F*>, const string&, vector<TH1F*>&, 
				      const Int_t nMETBins = 1) const;
  void makeToyDiEMETWeightsHistograms(TRandom3&, TH1F&, const TH1F&, vector<TH1F*>&, 
				      const unsigned int iMETBin = 1) const;
  void reweightDefault(const VFLOAT&, const VFLOAT&, const TH1F&, TH1F*, const VFLOAT&, 
		       const VFLOAT&) const;
  void reweightDefault(const VFLOAT&, const VFLOAT&, const VFLOAT&, const VFLOAT&, const TH1F&, 
		       TH1F*, TH2F*, const VFLOAT&, const VFLOAT&) const;
  void reweightDefault(const VFLOAT&, const VFLOAT&, const VFLOAT&, const VFLOAT&, const VFLOAT&, 
		       const unsigned int, const Double_t*, const vector<vector<TH1F*> >&, TH1F*, 
		       TH2F*, const VFLOAT&, const VFLOAT&, const unsigned int iToy = 0) const;
  void reweightBinned(const TH2F&, const TH1F&, TH1F*, const unsigned int) const;
  void generateToys(vector<TH1F*>&, vector<TH1F*>&, const vector<TH1F*>, const unsigned int, 
		    const string&, const Double_t*, const unsigned int, const Double_t*, 
		    const VFLOAT&, const VFLOAT&, const VFLOAT&, const VFLOAT&) const;
  void generateToys(vector<TH1F*>&, vector<TH2F*>&, vector<TH1F*>&, const vector<TH1F*>, 
		    const unsigned int, const string&, const Double_t*, 
		    const unsigned int, const Double_t*, const unsigned int, const Double_t*, 
		    const unsigned int, const Double_t*, const VFLOAT&, const VFLOAT&, 
		    const VFLOAT&, const VFLOAT&, const VFLOAT&, const VFLOAT&) const;
  void generateToys(vector<TH1F*>&, vector<TH1F*>&, const vector<TH1F*>, const unsigned int, 
		    const string&, const Double_t*) const;
  void generateToys(vector<TH1F*>&, vector<TH1F*>&, const vector<TH1F*>, const unsigned int, 
		    const string&, const Double_t*, const unsigned int, const Double_t*, 
		    const TH2F&) const;
  void makeToyMETDists(const unsigned int, const string&, const unsigned int, const Double_t*, 
		       const unsigned int, const Double_t*, const unsigned int, const Double_t*, 
		       const VFLOAT&, const VFLOAT&, const VFLOAT&, const VFLOAT&, const VFLOAT&, 
		       const unsigned int, const Double_t*, const vector<vector<TH1F*> >&, 
		       vector<TH1F*>&, vector<TH2F*>&, const VFLOAT&, const VFLOAT&) const;
  void generateToyDijets(const unsigned int, TRandom&, const vector<float>&, 
			 const vector<float>&, const vector<float>&, const vector<float>&, 
			 const vector<TLorentzVector>&, const vector<TLorentzVector>&, 
			 TH1F&) const;
  void estimateJESError(const unsigned int, vector<TH1F*>&, const vector<TH1F*>&, 
			const unsigned int, const Double_t*, const vector<float>&, 
			const vector<float>&, const vector<float>&, const vector<float>&, 
			const vector<TLorentzVector>&, const vector<TLorentzVector>&, 
			const vector<float>&, const vector<float>&, const vector<float>&, 
			const vector<float>&, const vector<TLorentzVector>&, 
			const vector<TLorentzVector>&) const;
  void fillToyDistributions(vector<TH1F*>&, const TH1F&, TCanvas&, vector<TH1F*>&, 
			    const unsigned int, const string&, const unsigned int) const;
  void fillToyDistributions(vector<TH1F*>&, vector<TH1F*>&, const TH1F&, const TH2F&, TCanvas&, 
			    vector<TH1F*>&, vector<TH2F*>&, const unsigned int, const string&, 
			    const unsigned int, const unsigned int, const unsigned int) const;
  void setMETErrorBars(TH1F&, const TH1F&, const TH1F&, const vector<TH1F*>&,
		       const vector<TH1F*>&, const vector<TH1F*>&, const float, const float,
		       const float, const float, const float, ofstream&,
		       const bool doDefault = true) const;
  void setRazorErrorBars(TH2F&, const vector<TH1F*>&, const vector<TH1F*>&, const vector<TH1F*>&, 
			 const float, const float, const bool doDefault = true) const;
  void makeFinalCanvas(TH1*, const Color_t, const Width_t, const Style_t, const Color_t, 
		       const Size_t, const Color_t, const string&) const;
  template<typename T>
    void deallocateMemory(vector<T>& vec) const
    {
      for (typename vector<T>::iterator i = vec.begin(); i != vec.end(); ++i) {
	delete *i;
	*i = NULL;
      }
    }
  template<typename T>
    void deallocateMemory(vector<vector<T> >& vec) const
    {
      for (typename vector<vector<T> >::iterator i = vec.begin(); i != vec.end(); ++i) {
	for (typename vector<T>::iterator j = i->begin(); j != i->end(); ++j) {
	  delete *j;
	  *j = NULL;
	}
      }
    }
  void printDiObjectErrorMessage(const unsigned int, const string&, const int, 
				 const unsigned int) const;
  string eventFileName(string, const string&) const;
  float HT() const;
  float MHT() const;
  unsigned int numJets(const float) const;
  float leadingJetET() const;
  unsigned int isJet(const susy::PFJet&, const float, const float ETMin = 30.0) const;
  susy::PFParticle* nearestPFCandidate(const susy::Photon*);
  void fillDRHistogram(TH2F&) const;
  float fitDiEMInvariantMass(vector<TH3F*>&, const unsigned int, const int,
			     const unsigned int, const int, float&, const unsigned int) const;
  float electronPhotonMisIDRate(vector<TH3F*>&, vector<TH3F*>&, const unsigned int,
				const int, const unsigned int, const int,
				float&) const;
  void runMETAnalysis(const string);
  void runMETAnalysisWithEEBackgroundFit(const std::string&);
  void testFitting(const string&, const string&) const;
  void runEEVsFFAnalysis(const std::string&);
  bool passJetFiducialCuts(const TLorentzVector&, const Double_t, const Double_t) const;
  bool passPFJetID(susy::PFJetCollection::const_iterator&) const;
  bool passCaloJetID(susy::CaloJetCollection::const_iterator&) const;
  unsigned int numPhotonOverlaps(const susy::PhotonCollection&, const TLorentzVector&, int&, 
				 const float) const;
  template<typename T>
    T* rightJetHistogram(const string& photonType, const string& identifier, 
			 const string& jetType, map<string, T*>& histograms) const
    {
      stringstream name;
      name << photonType << identifier << "EMFraction" << jetType;
      typename map<string, T*>::const_iterator iHist = histograms.find(name.str());
      if (iHist != histograms.end()) return iHist->second;
      else return NULL;
    }
  template<typename T>
    T* EMFractionHistogram(const int photon1Index, const int photon2Index, 
			   const susy::PhotonCollection& photons, const unsigned int numOverlaps, 
			   map<string, T*>& histograms, const string& identifier, 
			   const string& jetType) const
    {
      T* pHist = NULL;
      if ((photon1Index == -1) && (numOverlaps > 0)) {
	cerr << "Error: numOverlaps = " << numOverlaps;
	cerr << " but photon1Index = -1.  Skipping this jet.\n";
	return pHist;
      }
      int photon1Type = FAIL;
      if (photon1Index != -1) photon1Type = susyCategory->getPhotonType(tag_, photon1Index);
      string photonType;
      if ((photon1Index != -1) && (photon2Index != -1)) {
	if (photons[photon1Index].momentum.Et() > 
	    photons[photon2Index].momentum.Et()) photonType = "Leading";
	else photonType = "Trailing";
      }
      switch (numOverlaps) {
      case 1:
	switch (photon1Type) {
	case G:
	  pHist = rightJetHistogram("g" + photonType, identifier, jetType, histograms);
	  break;
	case E:
	  pHist = rightJetHistogram("e" + photonType, identifier, jetType, histograms);
	  break;
	case F:
	  pHist = rightJetHistogram("f" + photonType, identifier, jetType, histograms);
	  break;
	case FAIL:
	  pHist = rightJetHistogram("fail", identifier, jetType, histograms);
	  break;
	default:
	  cerr << "Error: photon type " << photon1Type << " invalid for photon ";
	  cerr << photon1Index << ".  Skipping this photon.\n";
	  break;
	}
	break;
      case 0:
	pHist = rightJetHistogram("other", identifier, jetType, histograms);
	break;
      default:
	//count number of jets per event with >1 EM overlap
	break;
      }
      return pHist;
    }
  vector<unsigned int> getDecidingPhotonIndices(const unsigned int) const;
  void runEMFractionAnalysis(const string&);
  void compareDataToMC(const string&);
  void runCutFlowAnalysis();
  void skim(const string&, const int);
  void stripBranch(const string&, const string&);
  void debugPrint(const unsigned int) const;
  void reset();
  void initPU();
  void clearPU();
  string printCategory(const int) const;
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
  bool passUserHLT() const; /*if trigger fires in specified active run range supplied by user, 
			      return true*/

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
  string JECErrFile_;                     //file name for JEC errors
  map<pair<TString, unsigned int>, pair<unsigned int, unsigned int> > HLT_; /*map of HLT paths 
									      and corresponding 
									      run ranges to 
									      consider*/
  map<string, float> eeInvMassFitParLowerBounds_;
  map<string, float> eeInvMassFitParUpperBounds_;
  map<string, float> eeInvMassFitPars_;
  map<string, float> egInvMassFitParLowerBounds_;
  map<string, float> egInvMassFitParUpperBounds_;
  map<string, float> egInvMassFitPars_;
  map<string, string> diEMInvMassFitParUnits_;
	 	 
  const map<string, float>* getDiEMInvMassFitParMapFloat(const unsigned int,
							 const unsigned int) const;
	 	 
  map<string, float>* getDiEMInvMassFitParMapFloat(const unsigned int, const unsigned int);
	 	 
  const map<string, string>* getDiEMInvMassFitParMapString(const unsigned int) const;
	 	 
  map<string, string>* getDiEMInvMassFitParMapString(const unsigned int);
	 	 
  void setDiEMInvMassFitParPropertyFloat(const unsigned int, const unsigned int, const string,
					 const float, const bool);
	 	 
  void setDiEMInvMassFitParPropertyString(const unsigned int, const string, const string,
					  const bool);
	 	 
  float getDiEMInvMassFitParPropertyFloat(const unsigned int, const unsigned int,
					  const string) const;
	 	 
  string getDiEMInvMassFitParPropertyString(const unsigned int, const string) const;
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
}

GMSBAnalyzer::~GMSBAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  fChain = NULL;

  //clear private data members
  tag_ = "";
  nEvts_ = 0;
  intLumi_ = 0.0;
  PUFile_ = "";
  L1JECFile_ = "";
  L2JECFile_ = "";
  L3JECFile_ = "";
  JECErrFile_ = "";
  HLT_.clear();
  eeInvMassFitParLowerBounds_.clear();
  eeInvMassFitParUpperBounds_.clear();
  eeInvMassFitPars_.clear();
  egInvMassFitParLowerBounds_.clear();
  egInvMassFitParUpperBounds_.clear();
  egInvMassFitPars_.clear();
  reset();
  clearPU();
  delete susyEvent;
  delete susyCategory;
  susyEvent = NULL;
  susyCategory = NULL;
  b_susyEvent_ = NULL;
  b_susyCategory_ = NULL;
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

  //initialize private data members
  nEvts_ = 0;
  intLumi_ = 0.0;
  reset();
  clearPU();
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
void GMSBAnalyzer::setJECErrFile(const string& JECErrFile) { JECErrFile_ = JECErrFile; }
void GMSBAnalyzer::addHLT(const TString& HLTPath, const unsigned int allowAbove, 
			  const unsigned int minRun, const unsigned int maxRun)
{
  map<pair<TString, unsigned int>, pair<unsigned int, unsigned int> >::const_iterator iHLT = 
    HLT_.find(pair<TString, unsigned int>(HLTPath, allowAbove));
  if (iHLT != HLT_.end()) {
    cout << "Replacing min. run " << iHLT->second.first << " with " << minRun << " and max. run ";
    cout << iHLT->second.second << " with " << maxRun << " for " << HLTPath;
    cout << ", allowed above category threshold " << allowAbove << endl;
  }
  HLT_[pair<TString, unsigned int>(HLTPath, allowAbove)] = 
    pair<unsigned int, unsigned int>(minRun, maxRun);
}
void GMSBAnalyzer::setDiEMInvMassFitParLowerBound(const unsigned int sample, const string name,
						  const float val, const bool overwrite)
{
  setDiEMInvMassFitParPropertyFloat(sample, LOWER_BOUND, name, val, overwrite);
}
void GMSBAnalyzer::setDiEMInvMassFitParUpperBound(const unsigned int sample, const string name,
						  const float val, const bool overwrite)
{
  setDiEMInvMassFitParPropertyFloat(sample, UPPER_BOUND, name, val, overwrite);
}
void GMSBAnalyzer::setDiEMInvMassFitPar(const unsigned int sample, const string name,
					const float val, const bool overwrite)
{
  setDiEMInvMassFitParPropertyFloat(sample, VALUE, name, val, overwrite);
}
void GMSBAnalyzer::setDiEMInvMassFitParUnit(const string name, const string val,
					    const bool overwrite)
{
  setDiEMInvMassFitParPropertyString(UNIT, name, val, overwrite);
}

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
string GMSBAnalyzer::getJECErrFile() const { return JECErrFile_; }
const map<pair<TString, unsigned int>, pair<unsigned int, unsigned int> >* 
  GMSBAnalyzer::getHLT() const { return &HLT_; }
unsigned int GMSBAnalyzer::getMinRun(const TString& HLT) const
{
  unsigned int minRun = 0;
  map<pair<TString, unsigned int>, pair<unsigned int, unsigned int> >::const_iterator iHLT = 
    HLT_.begin();
  while ((iHLT != HLT_.end()) && (minRun == 0)) {
    if (iHLT->first.first == HLT) minRun = iHLT->second.first;
    else ++iHLT;
  }
  return minRun;
}
unsigned int GMSBAnalyzer::getMaxRun(const TString& HLT) const
{
  unsigned int maxRun = 0;
  map<pair<TString, unsigned int>, pair<unsigned int, unsigned int> >::const_iterator iHLT = 
    HLT_.begin();
  while ((iHLT != HLT_.end()) && (maxRun == 0)) {
    if (iHLT->first.first == HLT) maxRun = iHLT->second.second;
    else ++iHLT;
  }
  return maxRun;
}
float GMSBAnalyzer::getDiEMInvMassFitParLowerBound(const unsigned int sample,
						   const string name) const
{
  float lowerBound = 0.0;
  try { lowerBound = getDiEMInvMassFitParPropertyFloat(sample, LOWER_BOUND, name); }
  catch (string badName) {}
  return lowerBound;
}
float GMSBAnalyzer::getDiEMInvMassFitParUpperBound(const unsigned int sample,
						   const string name) const
{
  float upperBound = 0.0;
  try { upperBound = getDiEMInvMassFitParPropertyFloat(sample, UPPER_BOUND, name); }
  catch (string badName) {}
  return upperBound;
}
float GMSBAnalyzer::getDiEMInvMassFitPar(const unsigned int sample,
					 const string name) const
{
  float par = 0.0;
  try { par = getDiEMInvMassFitParPropertyFloat(sample, VALUE, name); }
  catch (string badName) {}
  return par;
}
string GMSBAnalyzer::getDiEMInvMassFitParUnit(const string name) const
{
  string par;
  try { par = getDiEMInvMassFitParPropertyString(UNIT, name); }
  catch (string badName) {}
  return par;
}

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
string GMSBAnalyzer::printCategory(const int categoryInt) const
{
  string categoryString;
  switch (categoryInt) {
  case GG:
    categoryString = "GG";
    break;
  case EG:
    categoryString = "EG";
    break;
  case EE:
    categoryString = "EE";
    break;
  case FF:
    categoryString = "FF";
    break;
  case FAIL:
    categoryString = "FAIL";
    break;
  default:
    break;
  }
  return categoryString;
}
bool GMSBAnalyzer::passUserHLT() const
{
  bool pass = false;

  //loop over supplied HLT paths
  map<pair<TString, unsigned int>, pair<unsigned int, unsigned int> >::const_iterator iHLT = 
    HLT_.begin();
  while ((iHLT != HLT_.end()) && !pass) {

    //get applicable run range and category threshold for this path
    const unsigned int minRun = iHLT->second.first;
    const unsigned int maxRun = iHLT->second.second;
    const unsigned int allowAbove = iHLT->first.second;

    //loop over trigger information for this event
    susy::TriggerMap::const_iterator iInfo = susyEvent->hltMap.begin();
    while ((iInfo != susyEvent->hltMap.end()) && !pass) {

      /*if the event fired the path in question AND it is in the correct run range for the path 
	AND it is the correct category, return true*/
      if (iInfo->first.Contains(iHLT->first.first) && (int(iInfo->second.second)) && 
	  (susyCategory->getEventCategory(tag_) >= (int)allowAbove) && 
	  (susyEvent->runNumber >= (int)minRun) && (susyEvent->runNumber <= (int)maxRun)) {
	pass = true;
      }
      else ++iInfo;
    }

    //increment supplied HLT path counter
    ++iHLT;
  }

  return pass;
}

const map<string, float>*
  GMSBAnalyzer::getDiEMInvMassFitParMapFloat(const unsigned int sample,
					     const unsigned int property) const
{
  const map<string, float>* diEMInvMassFitParMap = NULL;
  switch (property) {
  case LOWER_BOUND:
    switch (sample) {
    case EE:
      diEMInvMassFitParMap = &eeInvMassFitParLowerBounds_;
      break;
    case EG:
      diEMInvMassFitParMap = &egInvMassFitParLowerBounds_;
      break;
    default:
      cerr << "Invalid sample " << sample << ".  No action taken.\n";
      break;
    }
    break;
  case UPPER_BOUND:
    switch (sample) {
    case EE:
      diEMInvMassFitParMap = &eeInvMassFitParUpperBounds_;
      break;
    case EG:
      diEMInvMassFitParMap = &egInvMassFitParUpperBounds_;
      break;
    default:
      cerr << "Invalid sample " << sample << ".  No action taken.\n";
      break;
    }
    break;
  case VALUE:
    switch (sample) {
    case EE:
      diEMInvMassFitParMap = &eeInvMassFitPars_;
      break;
    case EG:
      diEMInvMassFitParMap = &egInvMassFitPars_;
      break;
    default:
      cerr << "Invalid sample " << sample << ".  No action taken.\n";
      break;
    }
    break;
  default:
    cerr << "Invalid di-EM invariant mass fit parameter property " << property;
    cerr << ".  No action taken.\n";
    break;
  }
  return diEMInvMassFitParMap;
}
	 	 
map<string, float>* GMSBAnalyzer::getDiEMInvMassFitParMapFloat(const unsigned int sample,
							       const unsigned int property)
{
  map<string, float>* diEMInvMassFitParMap = NULL;
  switch (property) {
  case LOWER_BOUND:
    switch (sample) {
    case EE:
      diEMInvMassFitParMap = &eeInvMassFitParLowerBounds_;
      break;
    case EG:
      diEMInvMassFitParMap = &egInvMassFitParLowerBounds_;
      break;
    default:
      cerr << "Invalid sample " << sample << ".  No action taken.\n";
      break;
    }
    break;
  case UPPER_BOUND:
    switch (sample) {
    case EE:
      diEMInvMassFitParMap = &eeInvMassFitParUpperBounds_;
      break;
    case EG:
      diEMInvMassFitParMap = &egInvMassFitParUpperBounds_;
      break;
    default:
      cerr << "Invalid sample " << sample << ".  No action taken.\n";
      break;
    }
    break;
  case VALUE:
    switch (sample) {
    case EE:
      diEMInvMassFitParMap = &eeInvMassFitPars_;
      break;
    case EG:
      diEMInvMassFitParMap = &egInvMassFitPars_;
      break;
    default:
      cerr << "Invalid sample " << sample << ".  No action taken.\n";
      break;
    }
    break;
  default:
    cerr << "Invalid di-EM invariant mass fit parameter property " << property;
    cerr << ".  No action taken.\n";
    break;
  }
  return diEMInvMassFitParMap;
}
	 	 
const map<string, string>*
  GMSBAnalyzer::getDiEMInvMassFitParMapString(const unsigned int property) const
{
  const map<string, string>* diEMInvMassFitParMap = NULL;
  if (property == UNIT) return &diEMInvMassFitParUnits_;
  else {
    cerr << "Invalid di-EM invariant mass fit parameter property " << property;
    cerr << ".  No action taken.\n";
  }
  return diEMInvMassFitParMap;
}
	 	 
map<string, string>* GMSBAnalyzer::getDiEMInvMassFitParMapString(const unsigned int property)
{
  map<string, string>* diEMInvMassFitParMap = NULL;
  if (property == UNIT) diEMInvMassFitParMap = &diEMInvMassFitParUnits_;
  else {
    cerr << "Invalid di-EM invariant mass fit parameter property " << property;
    cerr << ".  No action taken.\n";
  }
  return diEMInvMassFitParMap;
}
	 	 
void GMSBAnalyzer::setDiEMInvMassFitParPropertyFloat(const unsigned int sample,
						     const unsigned int property,
						     const string name, const float val,
						     const bool overwrite)
{
  map<string, float>* diEMInvMassFitParMap = getDiEMInvMassFitParMapFloat(sample, property);
  if (diEMInvMassFitParMap == NULL) return;
  map<string, float>::iterator iMap = diEMInvMassFitParMap->find(name);
  if ((iMap == diEMInvMassFitParMap->end()) || overwrite) {
    diEMInvMassFitParMap->insert(pair<string, float>(name, val));
  }
  else {
    cerr << "Property " << property << " for parameter " << name << " already set to ";
    cerr << iMap->second << " and overwrites not allowed.\n";
  }
}
	 	 
void GMSBAnalyzer::setDiEMInvMassFitParPropertyString(const unsigned int property,
						      const string name, const string val,
						      const bool overwrite)
{
  map<string, string>* diEMInvMassFitParMap = getDiEMInvMassFitParMapString(property);
  if (diEMInvMassFitParMap == NULL) return;
  map<string, string>::iterator iMap = diEMInvMassFitParMap->find(name);
  if ((iMap == diEMInvMassFitParMap->end()) || overwrite) {
    diEMInvMassFitParMap->insert(pair<string, string>(name, val));
  }
  else {
    cerr << "Property " << property << " for parameter " << name << " already set to ";
    cerr << iMap->second << " and overwrites not allowed.\n";
  }
}
	 	 
float GMSBAnalyzer::getDiEMInvMassFitParPropertyFloat(const unsigned int sample,
						      const unsigned int property,
						      const string name) const
{
  const map<string, float>* diEMInvMassFitParMap = getDiEMInvMassFitParMapFloat(sample, property);
  if (diEMInvMassFitParMap == NULL) throw property;
  map<string, float>::const_iterator iMap = diEMInvMassFitParMap->find(name);
  if (iMap != diEMInvMassFitParMap->end()) {
    return iMap->second;
  }
  else {
    cerr << "Parameter " << name << " not found.\n";
    throw name;
  }
}
	 	 
string GMSBAnalyzer::getDiEMInvMassFitParPropertyString(const unsigned int property,
							const string name) const
{
  const map<string, string>* diEMInvMassFitParMap = getDiEMInvMassFitParMapString(property);
  if (diEMInvMassFitParMap == NULL) throw property;
  map<string, string>::const_iterator iMap = diEMInvMassFitParMap->find(name);
  if (iMap != diEMInvMassFitParMap->end()) {
    return iMap->second;
  }
  else {
    cerr << "Parameter " << name << " not found.\n";
    throw name;
  }
}

#endif // #ifdef GMSBAnalyzer_cxx
