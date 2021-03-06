// -*- C++ -*-
//
// Package:    SampleMaker
// Class:      SampleMaker
// 
/**\class SampleMaker SampleMaker.cc GMSBTools/SampleMaker/src/SampleMaker.cc

 Description: create a sample of interest to the GMSB analysis community

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu May  6 15:10:56 CEST 2010
// $Id: SampleMaker.cc,v 1.3 2010/06/23 11:10:39 yohay Exp $
//
//

// user include files
#include "FWCore/Framework/interface/EDFilter.h"

//other relevant header files
#include "GMSBTools/EventSelection/interface/EventSelector.h"

//relevant namespaces
using namespace std;
using namespace reco;
using namespace edm;
using namespace cms;

//
// class declaration
//

class SampleMaker : public edm::EDFilter {
   public:
      explicit SampleMaker(const edm::ParameterSet&);
      ~SampleMaker();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  unsigned int sampleType_;
  double ECALIsoMaxPTMultiplierEB_;
  double ECALIsoMaxConstantEB_;
  double ECALIsoMaxPTMultiplierEE_;
  double ECALIsoMaxConstantEE_;
  double HCALIsoMaxPTMultiplierEB_;
  double HCALIsoMaxConstantEB_;
  double HCALIsoMaxPTMultiplierEE_;
  double HCALIsoMaxConstantEE_;
  double HOverEMaxPresel_;
  double ETMin_;
  unsigned int fiducialRegion_;
  bool useHOverE_;
  double HOverEMax_;
  bool useSigmaEtaEta_;
  double sigmaEtaEtaMax_;
  bool useTrackIso_;
  double trackIsoMaxPTMultiplier_;
  double trackIsoMaxConstant_;
  double trackPTMin_;
  double eTrackRMin_;
  double minDRPhotons_;
  double maxSeedTime_;
  InputTag photonTag_;
  InputTag trackTag_;
  InputTag HBHERecHitTag_;
  InputTag cosmicTrackTag_;
  InputTag EBRecHitTag_;
  InputTag EERecHitTag_;
  string debugFileName_;
  bool debugFlag_;
  bool checkHaloCoincidenceWithPassingEBCandsOnly_;
  bool rejectHalo_;
  unsigned int numReqdCands_;
  unsigned int numTot_;
  unsigned int numPassing_;
  ESHandle<CaloGeometry> caloGeometryHandle_;
  EventSelector evtProperties_;
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
SampleMaker::SampleMaker(const edm::ParameterSet& iConfig) :
  sampleType_(iConfig.getParameter<unsigned int>("sampleType")),
  ECALIsoMaxPTMultiplierEB_(iConfig.getParameter<double>("ECALIsoMaxPTMultiplierEB")),
  ECALIsoMaxConstantEB_(iConfig.getParameter<double>("ECALIsoMaxConstantEB")),
  ECALIsoMaxPTMultiplierEE_(iConfig.getParameter<double>("ECALIsoMaxPTMultiplierEE")),
  ECALIsoMaxConstantEE_(iConfig.getParameter<double>("ECALIsoMaxConstantEE")),
  HCALIsoMaxPTMultiplierEB_(iConfig.getParameter<double>("HCALIsoMaxPTMultiplierEB")),
  HCALIsoMaxConstantEB_(iConfig.getParameter<double>("HCALIsoMaxConstantEB")),
  HCALIsoMaxPTMultiplierEE_(iConfig.getParameter<double>("HCALIsoMaxPTMultiplierEE")),
  HCALIsoMaxConstantEE_(iConfig.getParameter<double>("HCALIsoMaxConstantEE")),
  HOverEMaxPresel_(iConfig.getParameter<double>("HOverEMaxPresel")),
  ETMin_(iConfig.getParameter<double>("ETMin")),
  fiducialRegion_(iConfig.getParameter<unsigned int>("fiducialRegion")),
  useHOverE_(iConfig.getUntrackedParameter<bool>("useHOverE", false)),
  HOverEMax_(iConfig.getParameter<double>("HOverEMax")),
  useSigmaEtaEta_(iConfig.getUntrackedParameter<bool>("useSigmaEtaEta", false)),
  sigmaEtaEtaMax_(iConfig.getParameter<double>("sigmaEtaEtaMax")),
  useTrackIso_(iConfig.getUntrackedParameter<bool>("useTrackIso", false)),
  trackIsoMaxPTMultiplier_(iConfig.getParameter<double>("trackIsoMaxPTMultiplier")),
  trackIsoMaxConstant_(iConfig.getParameter<double>("trackIsoMaxConstant")),
  trackPTMin_(iConfig.getParameter<double>("trackPTMin")),
  eTrackRMin_(iConfig.getParameter<double>("eTrackRMin")),
  minDRPhotons_(iConfig.getParameter<double>("minDRPhotons")),
  maxSeedTime_(iConfig.getParameter<double>("maxSeedTime")),
  photonTag_(iConfig.getParameter<InputTag>("photonTag")),
  trackTag_(iConfig.getParameter<InputTag>("trackTag")),
  HBHERecHitTag_(iConfig.getParameter<InputTag>("HBHERecHitTag")),
  cosmicTrackTag_(iConfig.getParameter<InputTag>("cosmicTrackTag")),
  EBRecHitTag_(iConfig.getParameter<InputTag>("EBRecHitTag")),
  EERecHitTag_(iConfig.getParameter<InputTag>("EERecHitTag")),
  debugFileName_(iConfig.getUntrackedParameter<string>("debugFileName", "debug.txt")),
  debugFlag_(iConfig.getUntrackedParameter<bool>("debugFlag", false)),
  checkHaloCoincidenceWithPassingEBCandsOnly_(iConfig.getUntrackedParameter<bool>("checkHaloCoincidenceWithPassingEBCandsOnly", true)),
  rejectHalo_(iConfig.getUntrackedParameter<bool>("rejectHalo", true))
{
   //now do what ever initialization is needed
  numReqdCands_ = 2;
  if (sampleType_ == ETRACK) numReqdCands_ = 1;
  EventSelector evtProperties(sampleType_, ECALIsoMaxPTMultiplierEB_, ECALIsoMaxConstantEB_, ECALIsoMaxPTMultiplierEE_, ECALIsoMaxConstantEE_, 
			      HCALIsoMaxPTMultiplierEB_, HCALIsoMaxConstantEB_, HCALIsoMaxPTMultiplierEE_, HCALIsoMaxConstantEE_, HOverEMaxPresel_, 
			      ETMin_, fiducialRegion_, useHOverE_, HOverEMax_, useSigmaEtaEta_, sigmaEtaEtaMax_, useTrackIso_, trackIsoMaxPTMultiplier_, 
			      trackIsoMaxConstant_, trackPTMin_, eTrackRMin_, minDRPhotons_, maxSeedTime_, numReqdCands_, 0, 0, 0, debugFileName_, 
			      debugFlag_);
  evtProperties_ = evtProperties;
}


SampleMaker::~SampleMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
SampleMaker::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get the run, event, and lumi section numbers
  const unsigned int runNum = iEvent.run();
  const unsigned int evtNum = iEvent.id().event();
  const unsigned int lumiNum = iEvent.getLuminosityBlock().luminosityBlock();
  evtProperties_.setRun(runNum);
  evtProperties_.setEvt(evtNum);
  evtProperties_.setLumiSec(lumiNum);
  evtProperties_.printEvtInfo();

  //increment total event counter
  ++numTot_;

  //pass flag
  bool pass = false;

  //get the ECAL RecHit collections
  vector<EcalRecHit*> ECALRecHits;
  Handle<EBRecHitCollection> pEBRecHits;
  bool foundEBRecHits = false;
  if ((fiducialRegion_ == EB) || (fiducialRegion_ == ECAL)) {
    try { foundEBRecHits = (iEvent.getByLabel(EBRecHitTag_, pEBRecHits)) && (pEBRecHits->size() > 0); }
    catch (cms::Exception& ex) {}
    if (!foundEBRecHits) {
      stringstream infoStream;
      infoStream << "No EB RecHit collection found in run" << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n";
      evtProperties_.printDebug(infoStream.str());
    }
    else {
      for (EBRecHitCollection::const_iterator iEBRecHit = pEBRecHits->begin(); iEBRecHit != pEBRecHits->end(); ++iEBRecHit) {
	ECALRecHits.push_back(const_cast<EcalRecHit*>(&*iEBRecHit));
      }
    }
  }
  if ((fiducialRegion_ == EEND) || (fiducialRegion_ == ECAL)) {
    Handle<EERecHitCollection> pEERecHits;
    bool foundEERecHits = false;
    try { foundEERecHits = (iEvent.getByLabel(EERecHitTag_, pEERecHits)) && (pEERecHits->size() > 0); }
    catch (cms::Exception& ex) {}
    if (!foundEERecHits) {
      stringstream infoStream;
      infoStream << "No EE RecHit collection found in run" << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n";
      evtProperties_.printDebug(infoStream.str());
    }
    else {
      for (EERecHitCollection::const_iterator iEERecHit = pEERecHits->begin(); iEERecHit != pEERecHits->end(); ++iEERecHit) {
	ECALRecHits.push_back(const_cast<EcalRecHit*>(&*iEERecHit));
      }
    }
  }

  //get the reco::Photon collection
  Handle<PhotonCollection> pPhotons;
  bool foundPhotons = false;
  bool allCandsFound = false;
  map<unsigned int, Photon*> passingWithoutPixelSeed;
  map<unsigned int, Photon*> passingWithPixelSeed;
  try { foundPhotons = (iEvent.getByLabel(photonTag_, pPhotons)) && (pPhotons->size() > 0); }
  catch (cms::Exception& ex) {}
  if (!foundPhotons) {
    stringstream infoStream;
    infoStream << "No reco::Photon collection found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n";
    evtProperties_.printDebug(infoStream.str());
  }
  else {

    //decide if the required number of passing objects were found
    if (evtProperties_.foundPhotonCandidates(pPhotons, ECALRecHits, passingWithoutPixelSeed, passingWithPixelSeed)) {

      //decide if a track was found
      if (sampleType_ == ETRACK) {
	Handle<TrackCollection> pTracks;
	bool foundTracks = false;
	try { foundTracks = (iEvent.getByLabel(trackTag_, pTracks)) && (pTracks->size() > 0); }
	catch (cms::Exception& ex) {}
	if (!foundTracks) {
	  stringstream infoStream;
	  infoStream << "No reco::Track collection found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n";
	  evtProperties_.printDebug(infoStream.str());
	}
	else {
	  if (evtProperties_.foundTrack(pTracks, passingWithPixelSeed)) allCandsFound = true;
	}
      }

      //all samples besides ETRACK have all candidates found
      else allCandsFound = true;
    }
  }

  //only check data quality on events with two passing candidates
  bool halo = false; //if required collections aren't available in the event, the halo tag stays false
  if (allCandsFound && rejectHalo_) {

    //get the HB/HE RecHits for calculating beam halo tags
    Handle<HBHERecHitCollection> pHBHERecHits;
    bool foundHBHERecHits = false;
    try { foundHBHERecHits = (iEvent.getByLabel(HBHERecHitTag_, pHBHERecHits)) && (pHBHERecHits->size() > 0); }
    catch (cms::Exception& ex) {}
    if (!foundHBHERecHits) {
      stringstream infoStream;
      infoStream << "No reco::HBHERecHit collection found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n";
      evtProperties_.printDebug(infoStream.str());
    }
    else {

      //get the HCAL geometry for calculating the rho and phi of the HB/HE RecHit (there's got to be a better way to do this!)
      iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle_);
      const CaloGeometry* pGeometry = caloGeometryHandle_.product();

      //create the collection of candidates to check for halo coincidence
      map<unsigned int, Photon*> passingCands;

      //only check coincidence with EB candidates that pass the selection
      if (checkHaloCoincidenceWithPassingEBCandsOnly_) {
	for (map<unsigned int, Photon*>::const_iterator iPhoton = passingWithoutPixelSeed.begin(); iPhoton != passingWithoutPixelSeed.end(); ++iPhoton) {
	  if (passingCands.find((*iPhoton).first) == passingCands.end()) passingCands[(*iPhoton).first] = (*iPhoton).second;
	  else {
	    stringstream infoStream;
	    infoStream << "Already found an element with index " << (*iPhoton).first << " in passingCands.  >1 element in passingWithoutPixelSeed has the ";
	    infoStream << "same index!\n";
	    evtProperties_.printDebug(infoStream.str());
	  }
	}
	for (map<unsigned int, Photon*>::const_iterator iElectron = passingWithPixelSeed.begin(); iElectron != passingWithPixelSeed.end(); ++iElectron) {
	  if (passingCands.find((*iElectron).first) == passingCands.end()) passingCands[(*iElectron).first] = (*iElectron).second;
	  else {
	    stringstream infoStream;
	    infoStream << "Already found an element with index " << (*iElectron).first << " in passingCands.  An element in passingWithPixelSeed has the ";
	    infoStream << "same index as an element in passingWithoutPixelSeed!\n";
	    evtProperties_.printDebug(infoStream.str());
	  }
	}
      }

      //check halo coincidence with all EB candidates
      else {
	for (PhotonCollection::const_iterator iPhoton = pPhotons->begin(); iPhoton != pPhotons->end(); ++iPhoton) {
	  passingCands[iPhoton - pPhotons->begin()] = const_cast<Photon*>(&*iPhoton);
	}
      }

      //only calculate the muon tag if the HE tag fails
      if (!evtProperties_.passesHEBeamHaloTag(pHBHERecHits, passingCands, pGeometry)) {

	//get the cosmic tracks for calculating beam halo tags
	Handle<TrackCollection> pCosmicTracks;
	bool foundCosmicTracks = false;
	try { foundCosmicTracks = (iEvent.getByLabel(cosmicTrackTag_, pCosmicTracks)) && (pCosmicTracks->size() > 0); }
	catch (cms::Exception& ex) {}
	if (!foundCosmicTracks) {
	  stringstream infoStream;
	  infoStream << "No reco::Track collection for cosmic tracks found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum;
	  evtProperties_.printDebug(infoStream.str());
	}
	else {
	  if (evtProperties_.passesMuonBeamHaloTag(pCosmicTracks, passingCands)) halo = true;
	}
      }

      //passed HE halo tag
      else halo = true;
    }//else
  }//if (allCandsFound)

  stringstream infoStream;
  if (allCandsFound && (!halo)) {
    pass = true;
    infoStream << "Event passes -- run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n---------------------------\n\n";
  }
  else infoStream << "Event fails.\n---------------------------\n\n";
  evtProperties_.printDebug(infoStream.str());

  //filter the event
  if (pass) ++numPassing_;
  return pass;
}
  
// ------------ method called once each job just before starting event loop  ------------
void 
SampleMaker::beginJob()
{
  //initialize event counters
  numTot_ = 0;
  numPassing_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SampleMaker::endJob() {
  cout << "Number of events in " << evtProperties_.sampleTypeString() << " sample: " << numPassing_ << "/" << numTot_ << " = ";
  cout << ((float)numPassing_/(float)numTot_);
  cout << endl;
  if (debugFlag_) {
    stringstream infoStream;
    infoStream << "Number of events in " << evtProperties_.sampleTypeString() << " sample: " << numPassing_ << "/" << numTot_ << " = ";
    infoStream << ((float)numPassing_/(float)numTot_) << endl;
    evtProperties_.printDebug(infoStream.str());
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SampleMaker);
