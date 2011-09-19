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

   TString getTag() const;
   int getNEvts() const;

   bool matchedToGenParticle(const susy::Photon&, const VINT&, const UChar_t) const;
   bool passesDenominatorSelection(const unsigned int, const susy::Photon&, const VINT&) const;
   bool passesNumeratorSelection(const unsigned int) const;
   unsigned int numGoodVertices() const;
   void countEE(string&);
   void debugPrint(const unsigned int) const;
 private:
   TString tag_;   //photon collection tag
   int     nEvts_; //number of events to process
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

TString GMSBAnalyzer::getTag() const { return tag_; }
int GMSBAnalyzer::getNEvts() const { return nEvts_; }
#endif // #ifdef GMSBAnalyzer_cxx
