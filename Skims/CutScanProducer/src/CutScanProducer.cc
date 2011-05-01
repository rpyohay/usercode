// -*- C++ -*-
//
// Package:    CutScanProducer
// Class:      CutScanProducer
// 
/**\class CutScanProducer CutScanProducer.cc Skims/CutScanProducer/src/CutScanProducer.cc

 Description: produce several value maps encoding which supplied cuts were passed by the photons

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Mar 30 18:03:27 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "Skims/SkimTypes/interface/CutScanKey.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// class declaration
//

class CutScanProducer : public edm::EDProducer {
   public:
      explicit CutScanProducer(const edm::ParameterSet&);
      ~CutScanProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  //retrieve collection from the event
  template<typename T>
  const bool getCollection_(T& pCollection, const edm::InputTag& tag, const edm::Event& iEvent)
  {
    bool collectionFound = false;
    try { collectionFound = iEvent.getByLabel(tag, pCollection); }
    catch (cms::Exception& ex) {}
    if (!collectionFound) {
      std::stringstream err;
      err << "No collection of type " << tag << " found in run " << iEvent.run() << ", event ";
      err << iEvent.id().event() << ", lumi section ";
      err << iEvent.getLuminosityBlock().luminosityBlock() << ".\n";
      edm::LogInfo("Error") << err.str();
    }
    return collectionFound;
  }

  //check that the two given vectors match in size.  return true if they do, false if they don't
  template <typename T>
  bool checkVectorSizeMatch(const T& vec1, const T& vec2, const std::string& vec1Name, 
			    const std::string& vec2Name) const
  {
    bool doScan = true;
    if (vec1.size() != vec2.size()) {
      std::stringstream err;
      err << "Error: " << vec1Name << " has " << vec1.size();
      err << " elements, " << vec2Name << " has " << vec2.size();
      err << " elements.  Scan will not be performed.\n";
      edm::LogError("BadInput") << err;
      doScan = false;
    }
    return doScan;
  }

  //implementation for a standard less/greater than cut
  template <typename T, typename U>
  unsigned int decisionWord(const T var, const U& scanCuts, const bool lessThan, 
			    const bool equalTo, const double ET) const
  {
    unsigned int theDecisionWord = 0x0;
    for (typename U::const_iterator iScanCut = scanCuts.begin(); iScanCut != scanCuts.end(); 
	 ++iScanCut) {
      std::stringstream mess;
      mess << "var = " << var << std::endl;
      mess << "*iScanCut = " << *iScanCut << std::endl;
      mess << "lessThan = " << lessThan << std::endl;
      mess << "equalTo = " << equalTo << std::endl;
      mess << "theDecisionWord = " << theDecisionWord << std::endl;
      edm::LogInfo("Debug") << mess.str();
      if (((lessThan && equalTo && (var <= *iScanCut)) || 
	   (lessThan && !equalTo && (var < *iScanCut)) || 
	   (!lessThan && equalTo && (var >= *iScanCut)) || 
	   (!lessThan && !equalTo && (var > *iScanCut))) && (ET > extETMin_)) {
	theDecisionWord = theDecisionWord | (0x1 << (iScanCut - scanCuts.begin()));
      }
    }
    return theDecisionWord;
  }

  //implementation for an isolation cut that depends on ET
  template <typename T, typename U>
  unsigned int decisionWord(const T var, const T ET, const U& scanCuts1, const U& scanCuts2) const
  {
    unsigned int theDecisionWord = 0x0;
    if (checkVectorSizeMatch(scanCuts1, scanCuts2, "scanCuts1", "scanCuts2")) {
      for (typename U::const_iterator iScanCut = scanCuts1.begin(); iScanCut != scanCuts1.end(); 
	   ++iScanCut) {
	std::stringstream mess;
	mess << "var = " << var << std::endl;
	mess << "ET = " << ET << std::endl;
	mess << "*iScanCut = " << *iScanCut << std::endl;
	mess << "scanCuts2[iScanCut - scanCuts1.begin()] = ";
	mess << scanCuts2[iScanCut - scanCuts1.begin()] << std::endl;
	mess << "theDecisionWord = " << theDecisionWord << std::endl;
	edm::LogInfo("Debug") << mess.str();
	if ((var < ((*iScanCut)*ET + scanCuts2[iScanCut - scanCuts1.begin()])) && 
	    (ET > extETMin_)) {
	  theDecisionWord = theDecisionWord | (0x1 << (iScanCut - scanCuts1.begin()));
	}
      }
    }
    else throw cms::Exception("BadInput") << "Error: vector size mismatch.\n";
    return theDecisionWord;
  }

  //initialize counters
  void initializeCounters(const bool, const std::vector<double>&, std::vector<unsigned int>&);

  //increment counters
  void incrementCounters(const std::vector<unsigned int>&, std::vector<unsigned int>&, const bool);

  //print info on percentage of events passing the scan points
  void printStats(std::vector<unsigned int>&, const bool, const std::string&) const;

  //put the product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, 
			   const std::vector<unsigned int>&, edm::Event&, const std::string&, 
			   const bool) const;
      
      // ----------member data ---------------------------

  //input
  edm::InputTag photonTag_;
  double extETMin_;
  std::vector<double> ETMinScan_;
  std::vector<double> ECALIsoMaxPTMultiplierScan_;
  std::vector<double> ECALIsoMaxConstantScan_;
  std::vector<double> HCALIsoMaxPTMultiplierScan_;
  std::vector<double> HCALIsoMaxConstantScan_;
  std::vector<double> HOverEMaxScan_;
  std::vector<double> trackIsoMaxPTMultiplierScan_;
  std::vector<double> trackIsoMaxConstantScan_;
  std::vector<double> sigmaIetaIetaMaxScan_;
  bool doETMinScan_;
  bool doECALIsoMaxScan_;
  bool doHCALIsoMaxScan_;
  bool doHOverEMaxScan_;
  bool doTrackIsoMaxScan_;
  bool doSigmaIetaIetaMaxScan_;

  //counters
  unsigned int numTot_;
  std::vector<unsigned int> numPassingETMinScan_;
  std::vector<unsigned int> numPassingECALIsoMaxScan_;
  std::vector<unsigned int> numPassingHCALIsoMaxScan_;
  std::vector<unsigned int> numPassingHOverEMaxScan_;
  std::vector<unsigned int> numPassingTrackIsoMaxScan_;
  std::vector<unsigned int> numPassingSigmaIetaIetaMaxScan_;

  //cross-counters
  std::vector<unsigned int> numPassingECALHCALIsoMaxScan_;
  std::vector<unsigned int> numPassingECALHCALIsoHOverEMaxScan_;
  std::vector<unsigned int> numPassingECALTrackIsoMaxScan_;
  std::vector<unsigned int> numPassingECALIsoSigmaIetaIetaMaxScan_;
  std::vector<unsigned int> numPassingHCALTrackIsoMaxScan_;
  std::vector<unsigned int> numPassingHCALIsoSigmaIetaIetaMaxScan_;
  std::vector<unsigned int> numPassingHOverETrackIsoMaxScan_;
  std::vector<unsigned int> numPassingHOverESigmaIetaIetaMaxScan_;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
CutScanProducer::CutScanProducer(const edm::ParameterSet& iConfig) :

  //input
  photonTag_(iConfig.getUntrackedParameter<edm::InputTag>("photonTag", 
							  edm::InputTag("photons", "", 
									"RECOCleaned"))),
  extETMin_(iConfig.getUntrackedParameter<double>("extETMin", 30.0/*GeV*/)),
  ETMinScan_(iConfig.getParameter<std::vector<double> >("ETMinScan")),
  ECALIsoMaxPTMultiplierScan_(iConfig.getParameter<std::vector<double> >
			      ("ECALIsoMaxPTMultiplierScan")),
  ECALIsoMaxConstantScan_(iConfig.getParameter<std::vector<double> >("ECALIsoMaxConstantScan")),
  HCALIsoMaxPTMultiplierScan_(iConfig.getParameter<std::vector<double> >
			      ("HCALIsoMaxPTMultiplierScan")),
  HCALIsoMaxConstantScan_(iConfig.getParameter<std::vector<double> >("HCALIsoMaxConstantScan")),
  HOverEMaxScan_(iConfig.getParameter<std::vector<double> >("HOverEMaxScan")),
  trackIsoMaxPTMultiplierScan_(iConfig.getParameter<std::vector<double> >
			       ("trackIsoMaxPTMultiplierScan")),
  trackIsoMaxConstantScan_(iConfig.getParameter<std::vector<double> >("trackIsoMaxConstantScan")),
  sigmaIetaIetaMaxScan_(iConfig.getParameter<std::vector<double> >("sigmaIetaIetaMaxScan")),
  doETMinScan_(iConfig.getUntrackedParameter<bool>("doETMinScan", true)),
  doECALIsoMaxScan_(iConfig.getUntrackedParameter<bool>("doECALIsoMaxScan", true)),
  doHCALIsoMaxScan_(iConfig.getUntrackedParameter<bool>("doHCALIsoMaxScan", true)),
  doHOverEMaxScan_(iConfig.getUntrackedParameter<bool>("doHOverEMaxScan", true)),
  doTrackIsoMaxScan_(iConfig.getUntrackedParameter<bool>("doTrackIsoMaxScan", true)),
  doSigmaIetaIetaMaxScan_(iConfig.getUntrackedParameter<bool>("doSigmaIetaIetaMaxScan", true))
{
  //sanity checks
  doECALIsoMaxScan_ = 
    checkVectorSizeMatch(ECALIsoMaxPTMultiplierScan_, ECALIsoMaxConstantScan_, 
			 "ECALIsoMaxPTMultiplierScan", "ECALIsoMaxConstantScan");
  doHCALIsoMaxScan_ = 
    checkVectorSizeMatch(HCALIsoMaxPTMultiplierScan_, HCALIsoMaxConstantScan_, 
			 "HCALIsoMaxPTMultiplierScan", "HCALIsoMaxConstantScan");
  doTrackIsoMaxScan_ = 
    checkVectorSizeMatch(trackIsoMaxPTMultiplierScan_, trackIsoMaxConstantScan_, 
			 "trackIsoMaxPTMultiplierScan", "trackIsoMaxConstantScan");

  //register products
  /*the value map variable is a short integer, with each bit holding the decision for a 
    particular 
    value of the cut (i.e. a value map of integers for the ECAL isolation scan, where the LSB 
    holds the decision for the 1st ECAL isolation scan point, the next bit holds the decision for 
    the 2nd scan point, etc.)
    probably only 3ish scan points per cut variable, so a char (8 bits) would be plenty
    this scheme reduces the number of value maps to be stored in the event by a factor ~ the 
    number of scan points per cut variable
    just need a key (another product) to get the cut values for each bit -- 
    cut_var bit value1 value2
    string  char double double
    ECALIso 0   0.X    3.Y
    ECALIso 1   0.Z    3.A
    H/E     0   0.B    -1.0
    store this as...a special datatype, a vector of them per event (1 for each cut variable)
   */
  if (doETMinScan_) produces<edm::ValueMap<unsigned int> >("passETMin");
  if (doECALIsoMaxScan_) produces<edm::ValueMap<unsigned int> >("passECALIsoMax");
  if (doHCALIsoMaxScan_) produces<edm::ValueMap<unsigned int> >("passHCALIsoMax");
  if (doHOverEMaxScan_) produces<edm::ValueMap<unsigned int> >("passHOverEMax");
  if (doTrackIsoMaxScan_) produces<edm::ValueMap<unsigned int> >("passTrackIsoMax");
  if (doSigmaIetaIetaMaxScan_) produces<edm::ValueMap<unsigned int> >("passSigmaIetaIetaMax");
  produces<CutScanKeyCollection>();
}


CutScanProducer::~CutScanProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CutScanProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /*write CutScanKey object to event (really only needed once per job, but simplest to do once 
    per event right now)*/
  std::auto_ptr<CutScanKeyCollection> out(new CutScanKeyCollection());
  if (doETMinScan_) out->push_back(CutScanKey("ETMin", ETMinScan_, std::vector<double>()));
  if (doECALIsoMaxScan_) {
    out->push_back(CutScanKey("ECALIsoMax", ECALIsoMaxPTMultiplierScan_, 
			      ECALIsoMaxConstantScan_));
  }
  if (doHCALIsoMaxScan_) {
    out->push_back(CutScanKey("HCALIsoMax", HCALIsoMaxPTMultiplierScan_, 
			      HCALIsoMaxConstantScan_));
  }
  if (doHOverEMaxScan_) out->push_back(CutScanKey("HOverEMax", HOverEMaxScan_, 
						  std::vector<double>()));
  if (doTrackIsoMaxScan_) {
    out->push_back(CutScanKey("trackIsoMax", trackIsoMaxPTMultiplierScan_, 
			      trackIsoMaxConstantScan_));
  }
  if (doSigmaIetaIetaMaxScan_) {
    out->push_back(CutScanKey("sigmaIetaIetaMax", sigmaIetaIetaMaxScan_, std::vector<double>()));
  }
  iEvent.put(out);

  //vectors to store the decision words for each photon, each cut variable
  std::vector<unsigned int> ETMinScanDecision;
  std::vector<unsigned int> ECALIsoMaxScanDecision;
  std::vector<unsigned int> HCALIsoMaxScanDecision;
  std::vector<unsigned int> HOverEMaxScanDecision;
  std::vector<unsigned int> trackIsoMaxScanDecision;
  std::vector<unsigned int> sigmaIetaIetaMaxScanDecision;

  //loop over photons
  edm::Handle<reco::PhotonCollection> pPhotons;
  if (getCollection_(pPhotons, photonTag_, iEvent)) {
    for (reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin(); 
	 iPhoton != pPhotons->end(); ++iPhoton) {
      std::stringstream mess;
      mess << "Photon index: " << iPhoton - pPhotons->begin() << std::endl;
      edm::LogInfo("Debug") << mess.str();

      //plot time per photon, time per event?

      //chars to hold the decision words for each cut variable
      unsigned int ETMinDecisionWord = decisionWord(iPhoton->et(), ETMinScan_, false, false, 
						    iPhoton->et());
      mess.str("");
      mess << "ETMinDecisionWord = " << ETMinDecisionWord << std::endl;
      edm::LogInfo("Debug") << mess.str();
      unsigned int ECALIsoMaxDecisionWord = 
	decisionWord((double)iPhoton->ecalRecHitSumEtConeDR04(), iPhoton->et(), 
		     ECALIsoMaxPTMultiplierScan_, ECALIsoMaxConstantScan_);
      mess.str("");
      mess << "ECALIsoMaxDecisionWord = " << ECALIsoMaxDecisionWord << std::endl;
      edm::LogInfo("Debug") << mess.str();
      unsigned int HCALIsoMaxDecisionWord = 
	decisionWord((double)iPhoton->hcalTowerSumEtConeDR04(), iPhoton->et(), 
		     HCALIsoMaxPTMultiplierScan_, HCALIsoMaxConstantScan_);
      mess.str("");
      mess << "HCALIsoMaxDecisionWord = " << HCALIsoMaxDecisionWord << std::endl;
      edm::LogInfo("Debug") << mess.str();
      unsigned int HOverEMaxDecisionWord = 
	decisionWord(iPhoton->hadronicOverEm(), HOverEMaxScan_, true, false, iPhoton->et());
      mess.str("");
      mess << "HOverEMaxDecisionWord = " << HOverEMaxDecisionWord << std::endl;
      edm::LogInfo("Debug") << mess.str();
      unsigned int trackIsoMaxDecisionWord = 
	decisionWord((double)iPhoton->trkSumPtHollowConeDR04(), iPhoton->et(), 
		     trackIsoMaxPTMultiplierScan_, trackIsoMaxConstantScan_);
      mess.str("");
      mess << "trackIsoMaxDecisionWord = " << trackIsoMaxDecisionWord << std::endl;
      edm::LogInfo("Debug") << mess.str();
      unsigned int sigmaIetaIetaMaxDecisionWord = 
	decisionWord(iPhoton->sigmaIetaIeta(), sigmaIetaIetaMaxScan_, true, false, iPhoton->et());
      mess.str("");
      mess << "sigmaIetaIetaMaxDecisionWord = " << sigmaIetaIetaMaxDecisionWord;
      mess << std::endl;
      edm::LogInfo("Debug") << mess.str();

      //save the decision words
      ETMinScanDecision.push_back(ETMinDecisionWord);
      ECALIsoMaxScanDecision.push_back(ECALIsoMaxDecisionWord);
      HCALIsoMaxScanDecision.push_back(HCALIsoMaxDecisionWord);
      HOverEMaxScanDecision.push_back(HOverEMaxDecisionWord);
      trackIsoMaxScanDecision.push_back(trackIsoMaxDecisionWord);
      sigmaIetaIetaMaxScanDecision.push_back(sigmaIetaIetaMaxDecisionWord);

    }/*for (reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin(); 
       iPhoton != pPhotons->end(); ++iPhoton)*/
  }//if (getCollection_(pPhotons, photonTag_, iEvent))

  //increment the pass counters
  ++numTot_;
  incrementCounters(ETMinScanDecision, numPassingETMinScan_, doETMinScan_);
  incrementCounters(ECALIsoMaxScanDecision, numPassingECALIsoMaxScan_, doECALIsoMaxScan_);
  incrementCounters(HCALIsoMaxScanDecision, numPassingHCALIsoMaxScan_, doHCALIsoMaxScan_);
  incrementCounters(HOverEMaxScanDecision, numPassingHOverEMaxScan_, doHOverEMaxScan_);
  incrementCounters(trackIsoMaxScanDecision, numPassingTrackIsoMaxScan_, doTrackIsoMaxScan_);
  incrementCounters(sigmaIetaIetaMaxScanDecision, numPassingSigmaIetaIetaMaxScan_, 
		    doSigmaIetaIetaMaxScan_);
  std::vector<unsigned int> ECALHCALIsoMaxScanDecision;
  for (std::vector<unsigned int>::const_iterator i = ECALIsoMaxScanDecision.begin(); 
       i != ECALIsoMaxScanDecision.end(); ++i) {
    ECALHCALIsoMaxScanDecision.push_back(*i & 
					 HCALIsoMaxScanDecision[i - 
								ECALIsoMaxScanDecision.begin()]);
  }
  incrementCounters(ECALHCALIsoMaxScanDecision, numPassingECALHCALIsoMaxScan_, 
		    doECALIsoMaxScan_ && doHCALIsoMaxScan_);
  std::vector<unsigned int> ECALHCALIsoHOverEMaxScanDecision;
  for (std::vector<unsigned int>::const_iterator i = ECALIsoMaxScanDecision.begin(); 
       i != ECALIsoMaxScanDecision.end(); ++i) {
    ECALHCALIsoHOverEMaxScanDecision.push_back(*i & 
					       HCALIsoMaxScanDecision
					       [i - ECALIsoMaxScanDecision.begin()] && 
					       HOverEMaxScanDecision
					       [i - ECALIsoMaxScanDecision.begin()]);
  }
  incrementCounters(ECALHCALIsoHOverEMaxScanDecision, numPassingECALHCALIsoHOverEMaxScan_, 
		     doECALIsoMaxScan_ && doHCALIsoMaxScan_ && doHOverEMaxScan_);
  std::vector<unsigned int> ECALTrackIsoMaxScanDecision;
  for (std::vector<unsigned int>::const_iterator i = ECALIsoMaxScanDecision.begin(); 
       i != ECALIsoMaxScanDecision.end(); ++i) {
    ECALTrackIsoMaxScanDecision.push_back(*i & 
					  trackIsoMaxScanDecision
					  [i - ECALIsoMaxScanDecision.begin()]);
  }
  incrementCounters(ECALTrackIsoMaxScanDecision, numPassingECALTrackIsoMaxScan_, 
		     doECALIsoMaxScan_ && doTrackIsoMaxScan_);
  std::vector<unsigned int> ECALIsoSigmaIetaIetaMaxScanDecision;
  for (std::vector<unsigned int>::const_iterator i = ECALIsoMaxScanDecision.begin(); 
       i != ECALIsoMaxScanDecision.end(); ++i) {
    ECALIsoSigmaIetaIetaMaxScanDecision.push_back(*i & 
						  sigmaIetaIetaMaxScanDecision
						  [i - ECALIsoMaxScanDecision.begin()]);
  }
  incrementCounters(ECALIsoSigmaIetaIetaMaxScanDecision, numPassingECALIsoSigmaIetaIetaMaxScan_, 
		     doECALIsoMaxScan_ && doSigmaIetaIetaMaxScan_);
  std::vector<unsigned int> HCALTrackIsoMaxScanDecision;
  for (std::vector<unsigned int>::const_iterator i = HCALIsoMaxScanDecision.begin(); 
       i != HCALIsoMaxScanDecision.end(); ++i) {
    HCALTrackIsoMaxScanDecision.push_back(*i & 
					  trackIsoMaxScanDecision
					  [i - HCALIsoMaxScanDecision.begin()]);
  }
  incrementCounters(HCALTrackIsoMaxScanDecision, numPassingHCALTrackIsoMaxScan_, 
		     doHCALIsoMaxScan_ && doTrackIsoMaxScan_);
  std::vector<unsigned int> HCALIsoSigmaIetaIetaMaxScanDecision;
  for (std::vector<unsigned int>::const_iterator i = HCALIsoMaxScanDecision.begin(); 
       i != HCALIsoMaxScanDecision.end(); ++i) {
    HCALIsoSigmaIetaIetaMaxScanDecision.push_back(*i & 
						  sigmaIetaIetaMaxScanDecision
						  [i - HCALIsoMaxScanDecision.begin()]);
  }
  incrementCounters(HCALIsoSigmaIetaIetaMaxScanDecision, numPassingHCALIsoSigmaIetaIetaMaxScan_, 
		     doHCALIsoMaxScan_ && doSigmaIetaIetaMaxScan_);
  std::vector<unsigned int> HOverETrackIsoMaxScanDecision;
  for (std::vector<unsigned int>::const_iterator i = HOverEMaxScanDecision.begin(); 
       i != HOverEMaxScanDecision.end(); ++i) {
    HOverETrackIsoMaxScanDecision.push_back(*i & 
					    trackIsoMaxScanDecision
					    [i - HOverEMaxScanDecision.begin()]);
  }
  incrementCounters(HOverETrackIsoMaxScanDecision, numPassingHOverETrackIsoMaxScan_, 
		     doHOverEMaxScan_ && doTrackIsoMaxScan_);
  std::vector<unsigned int> HOverESigmaIetaIetaMaxScanDecision;
  for (std::vector<unsigned int>::const_iterator i = HOverEMaxScanDecision.begin(); 
       i != HOverEMaxScanDecision.end(); ++i) {
    HOverESigmaIetaIetaMaxScanDecision.push_back(*i & 
						 sigmaIetaIetaMaxScanDecision
						 [i - HOverEMaxScanDecision.begin()]);
  }
  incrementCounters(HOverESigmaIetaIetaMaxScanDecision, numPassingHOverESigmaIetaIetaMaxScan_, 
		     doHOverEMaxScan_ && doSigmaIetaIetaMaxScan_);

  //put the collections into the event
  putProductIntoEvent(pPhotons, ETMinScanDecision, iEvent, "passETMin", doETMinScan_);
  putProductIntoEvent(pPhotons, ECALIsoMaxScanDecision, iEvent, "passECALIsoMax", 
		      doECALIsoMaxScan_);
  putProductIntoEvent(pPhotons, HCALIsoMaxScanDecision, iEvent, "passHCALIsoMax", 
		      doHCALIsoMaxScan_);
  putProductIntoEvent(pPhotons, HOverEMaxScanDecision, iEvent, "passHOverEMax", doHOverEMaxScan_);
  putProductIntoEvent(pPhotons, trackIsoMaxScanDecision, iEvent, "passTrackIsoMax", 
		      doTrackIsoMaxScan_);
  putProductIntoEvent(pPhotons, sigmaIetaIetaMaxScanDecision, iEvent, "passSigmaIetaIetaMax", 
		      doSigmaIetaIetaMaxScan_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CutScanProducer::beginJob()
{
  //initialize counters
  numTot_ = 0;
  initializeCounters(doETMinScan_, ETMinScan_, numPassingETMinScan_);
  initializeCounters(doECALIsoMaxScan_, ECALIsoMaxPTMultiplierScan_, numPassingECALIsoMaxScan_);
  initializeCounters(doHCALIsoMaxScan_, HCALIsoMaxPTMultiplierScan_, numPassingHCALIsoMaxScan_);
  initializeCounters(doHOverEMaxScan_, HOverEMaxScan_, numPassingHOverEMaxScan_);
  initializeCounters(doTrackIsoMaxScan_, trackIsoMaxPTMultiplierScan_, 
		     numPassingTrackIsoMaxScan_);
  initializeCounters(doSigmaIetaIetaMaxScan_, sigmaIetaIetaMaxScan_, 
		     numPassingSigmaIetaIetaMaxScan_);
  initializeCounters(doECALIsoMaxScan_ && doHCALIsoMaxScan_, ECALIsoMaxPTMultiplierScan_, 
		     numPassingECALHCALIsoMaxScan_);
  initializeCounters(doECALIsoMaxScan_ && doHCALIsoMaxScan_ && doHOverEMaxScan_, 
		     ECALIsoMaxPTMultiplierScan_, numPassingECALHCALIsoHOverEMaxScan_);
  initializeCounters(doECALIsoMaxScan_ && doTrackIsoMaxScan_, ECALIsoMaxPTMultiplierScan_, 
		     numPassingECALTrackIsoMaxScan_);
  initializeCounters(doECALIsoMaxScan_ && doSigmaIetaIetaMaxScan_, ECALIsoMaxPTMultiplierScan_, 
		     numPassingECALIsoSigmaIetaIetaMaxScan_);
  initializeCounters(doHCALIsoMaxScan_ && doTrackIsoMaxScan_, HCALIsoMaxPTMultiplierScan_, 
		     numPassingHCALTrackIsoMaxScan_);
  initializeCounters(doHCALIsoMaxScan_ && doSigmaIetaIetaMaxScan_, HCALIsoMaxPTMultiplierScan_, 
		     numPassingHCALIsoSigmaIetaIetaMaxScan_);
  initializeCounters(doHOverEMaxScan_ && doTrackIsoMaxScan_, HOverEMaxScan_, 
		     numPassingHOverETrackIsoMaxScan_);
  initializeCounters(doHOverEMaxScan_ && doSigmaIetaIetaMaxScan_, HOverEMaxScan_, 
		     numPassingHOverESigmaIetaIetaMaxScan_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CutScanProducer::endJob() {

  //print statistics (i.e. data volume reduction for each scan point)
  printStats(numPassingETMinScan_, doETMinScan_, "ETMin");
  printStats(numPassingECALIsoMaxScan_, doECALIsoMaxScan_, "ECALIsoMax");
  printStats(numPassingHCALIsoMaxScan_, doHCALIsoMaxScan_, "HCALIsoMax");
  printStats(numPassingHOverEMaxScan_, doHOverEMaxScan_, "HOverEMax");
  printStats(numPassingTrackIsoMaxScan_, doTrackIsoMaxScan_, "trackIsoMax");
  printStats(numPassingSigmaIetaIetaMaxScan_, doSigmaIetaIetaMaxScan_, "sigmaIetaIetaMax");
  printStats(numPassingECALHCALIsoMaxScan_, doECALIsoMaxScan_ && doHCALIsoMaxScan_, 
	     "ECALHCALIsoMax");
  printStats(numPassingECALHCALIsoHOverEMaxScan_, 
	     doECALIsoMaxScan_ && doHCALIsoMaxScan_ && doHOverEMaxScan_, "ECALHCALIsoHOverEMax");
  printStats(numPassingECALTrackIsoMaxScan_, doECALIsoMaxScan_ && doTrackIsoMaxScan_, 
	     "ECALTrackIsoMax");
  printStats(numPassingECALIsoSigmaIetaIetaMaxScan_, 
	     doECALIsoMaxScan_ && doSigmaIetaIetaMaxScan_, "ECALIsoSigmaIetaIetaMax");
  printStats(numPassingHCALTrackIsoMaxScan_, doHCALIsoMaxScan_ && doTrackIsoMaxScan_, 
	     "HCALTrackIsoMax");
  printStats(numPassingHCALIsoSigmaIetaIetaMaxScan_, 
	     doHCALIsoMaxScan_ && doSigmaIetaIetaMaxScan_, "HCALIsoSigmaIetaIetaMax");
  printStats(numPassingHOverETrackIsoMaxScan_, doHOverEMaxScan_ && doTrackIsoMaxScan_, 
	     "HOverETrackIsoMax");
  printStats(numPassingHOverESigmaIetaIetaMaxScan_, doHOverEMaxScan_ && doSigmaIetaIetaMaxScan_, 
	     "HOverESigmaIetaIetaMax");
}

void CutScanProducer::initializeCounters(const bool doScan, const std::vector<double>& scanVec, 
					 std::vector<unsigned int>& counterVec)
{
  if (doScan) {
    for (std::vector<double>::const_iterator i = scanVec.begin(); i != scanVec.end(); ++i) {
      counterVec.push_back(0);
    }
  }
}

void CutScanProducer::incrementCounters(const std::vector<unsigned int>& decisionVec, 
					std::vector<unsigned int>& counterVec, const bool doScan)
{
  if (doScan) {
    std::vector<unsigned int> numPhotonsPassing(counterVec.size(), 0);
    for (std::vector<unsigned int>::const_iterator i = decisionVec.begin(); 
	 i != decisionVec.end(); ++i) {
      for (unsigned int j = 0; j < counterVec.size(); ++j) {
	if (((*i >> j) & 0x1) == 1) ++numPhotonsPassing[j];
      }
    }
    for (std::vector<unsigned int>::const_iterator iPass = numPhotonsPassing.begin(); 
	 iPass != numPhotonsPassing.end(); ++iPass) {
      if (*iPass >= 2) ++counterVec[iPass - numPhotonsPassing.begin()];
    }
  }
}

void CutScanProducer::printStats(std::vector<unsigned int>& counterVec, const bool doScan, 
				 const std::string& label) const
{
  if (doScan) {
    std::stringstream stats;
    for (std::vector<unsigned int>::const_iterator i = counterVec.begin(); i != counterVec.end(); 
	 ++i) {
      stats << "Number passing " << label << " working point " << (i - counterVec.begin());
      stats << ": " << *i << " (";
      if (numTot_ > 0) stats << ((float)(*i)*100.0/(float)numTot_) << "%)\n";
      else stats << "N/A)\n";
    }
    stats << "-------------------\n";
    edm::LogInfo("Statistics") << stats.str();
  }
}

void CutScanProducer::putProductIntoEvent(const edm::Handle<reco::PhotonCollection>& pPhotons, 
					  const std::vector<unsigned int>& decision, 
					  edm::Event& iEvent, const std::string& label, 
					  const bool doScan) const
{
  if (doScan) {
    std::auto_ptr<edm::ValueMap<unsigned int> > out(new edm::ValueMap<unsigned int>());
    edm::ValueMap<unsigned int>::Filler filler(*out);
    filler.insert(pPhotons, decision.begin(), decision.end());
    filler.fill();
    iEvent.put(out, label);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(CutScanProducer);
