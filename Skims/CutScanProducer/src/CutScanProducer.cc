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
      edm::LogInfo("CutScanProducer") << err.str();
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
			    const bool equalTo) const
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
      edm::LogInfo("CutScanProducer") << mess.str();
      if ((lessThan && equalTo && (var <= *iScanCut)) || 
	  (lessThan && !equalTo && (var < *iScanCut)) || 
	  (!lessThan && equalTo && (var >= *iScanCut)) || 
	  (!lessThan && !equalTo && (var > *iScanCut))) {
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
	edm::LogInfo("CutScanProducer") << mess.str();
	if (var < ((*iScanCut)*ET + scanCuts2[iScanCut - scanCuts1.begin()])) {
	  theDecisionWord = theDecisionWord | (0x1 << (iScanCut - scanCuts1.begin()));
	}
      }
    }
    else throw cms::Exception("BadInput") << "Error: vector size mismatch.\n";
    return theDecisionWord;
  }

  //put the product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, 
			   const std::vector<unsigned int>&, edm::Event&, const std::string&, 
			   const bool) const;
      
      // ----------member data ---------------------------

  //input
  edm::InputTag photonTag_;
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
  ETMinScan_(iConfig.getParameter<std::vector<double> >("ETMinScan")),
  ECALIsoMaxPTMultiplierScan_(iConfig.getParameter<std::vector<double> >("ECALIsoMaxPTMultiplierScan")),
  ECALIsoMaxConstantScan_(iConfig.getParameter<std::vector<double> >("ECALIsoMaxConstantScan")),
  HCALIsoMaxPTMultiplierScan_(iConfig.getParameter<std::vector<double> >("HCALIsoMaxPTMultiplierScan")),
  HCALIsoMaxConstantScan_(iConfig.getParameter<std::vector<double> >("HCALIsoMaxConstantScan")),
  HOverEMaxScan_(iConfig.getParameter<std::vector<double> >("HOverEMaxScan")),
  trackIsoMaxPTMultiplierScan_(iConfig.getParameter<std::vector<double> >("trackIsoMaxPTMultiplierScan")),
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
      edm::LogInfo("CutScanProducer") << mess.str();

      //plot time per photon, time per event?

      //chars to hold the decision words for each cut variable
      unsigned int ETMinDecisionWord = decisionWord(iPhoton->et(), ETMinScan_, false, false);
      mess.str("");
      mess << "ETMinDecisionWord = " << ETMinDecisionWord << std::endl;
      edm::LogInfo("CutScanProducer") << mess.str();
      unsigned int ECALIsoMaxDecisionWord = 
	decisionWord((double)iPhoton->ecalRecHitSumEtConeDR04(), iPhoton->et(), 
		     ECALIsoMaxPTMultiplierScan_, ECALIsoMaxConstantScan_);
      mess.str("");
      mess << "ECALIsoMaxDecisionWord = " << ECALIsoMaxDecisionWord << std::endl;
      edm::LogInfo("CutScanProducer") << mess.str();
      unsigned int HCALIsoMaxDecisionWord = 
	decisionWord((double)iPhoton->hcalTowerSumEtConeDR04(), iPhoton->et(), 
		     HCALIsoMaxPTMultiplierScan_, HCALIsoMaxConstantScan_);
      mess.str("");
      mess << "HCALIsoMaxDecisionWord = " << HCALIsoMaxDecisionWord << std::endl;
      edm::LogInfo("CutScanProducer") << mess.str();
      unsigned int HOverEMaxDecisionWord = 
	decisionWord(iPhoton->hadronicOverEm(), HOverEMaxScan_, true, false);
      mess.str("");
      mess << "HOverEMaxDecisionWord = " << HOverEMaxDecisionWord << std::endl;
      edm::LogInfo("CutScanProducer") << mess.str();
      unsigned int trackIsoMaxDecisionWord = 
	decisionWord((double)iPhoton->trkSumPtHollowConeDR04(), iPhoton->et(), 
		     trackIsoMaxPTMultiplierScan_, trackIsoMaxConstantScan_);
      mess.str("");
      mess << "trackIsoMaxDecisionWord = " << trackIsoMaxDecisionWord << std::endl;
      edm::LogInfo("CutScanProducer") << mess.str();
      unsigned int sigmaIetaIetaMaxDecisionWord = 
	decisionWord(iPhoton->sigmaIetaIeta(), sigmaIetaIetaMaxScan_, true, false);
      mess.str("");
      mess << "sigmaIetaIetaMaxDecisionWord = " << sigmaIetaIetaMaxDecisionWord;
      mess << std::endl;
      edm::LogInfo("CutScanProducer") << mess.str();

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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CutScanProducer::endJob() {
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
