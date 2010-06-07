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
// $Id$
//
//

// user include files
#include "FWCore/Framework/interface/EDFilter.h"

//other relevant header files
#include "GMSBTools/SampleMaker/interface/SampleMaker.h"

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
  const unsigned int ECALFiducialRegion(const Photon*) const;
  const string sampleTypeString() const;
  const bool passesPreselection(const Photon*) const;

  //calculate dPhi between two objects
  template <typename T>
  const T dPhi(T phi1, T phi2) const
  {
    T dPhi = fabs(phi1 - phi2);
    if (dPhi > TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
    return dPhi;
  }

  //calculate dEta between two objects
  template <typename T>
  const T dEta(T eta1, T eta2) const { return fabs(eta1 - eta2); }

  //calculate dR between two objects
  template <typename T>
  const T dR(T eta1, T eta2, T phi1, T phi2) const { return sqrt(dEta(eta1, eta2)*dEta(eta1, eta2) + dPhi(phi1, phi2)*dPhi(phi1, phi2)); }
      
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
  InputTag photonTag_;
  InputTag trackTag_;
  InputTag HBHERecHitTag_;
  InputTag cosmicTrackTag_;
  string debugFileName_;
  ofstream debug_;
  bool debugFlag_;
  unsigned int numReqdCands_;
  unsigned int numTot_;
  unsigned int numPassing_;
  ESHandle<CaloGeometry> caloGeometryHandle_;
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
  photonTag_(iConfig.getParameter<InputTag>("photonTag")),
  trackTag_(iConfig.getParameter<InputTag>("trackTag")),
  HBHERecHitTag_(iConfig.getParameter<InputTag>("HBHERecHitTag")),
  cosmicTrackTag_(iConfig.getParameter<InputTag>("cosmicTrackTag")),
  debugFileName_(iConfig.getUntrackedParameter<string>("debugFileName", "debug.txt")),
  debugFlag_(iConfig.getUntrackedParameter<bool>("debugFlag", false))
{
   //now do what ever initialization is needed
  if (debugFlag_) {
    debug_.open(debugFileName_.c_str());
    if (!debug_.is_open()) {
      cerr << "Error opening file " << debugFileName_ << ".  Debugging will be turned off for this job.\n";
      debugFlag_ = false;
    }
  }
  if (debugFlag_) {
    debug_ << "Sample type: " << sampleType_ << " (" << sampleTypeString() << ")\n";
    debug_ << "Maximum ECAL isolation in EB: " << ECALIsoMaxPTMultiplierEB_ << "*pT + " << ECALIsoMaxConstantEB_ << " GeV\n";
    debug_ << "Maximum ECAL isolation in EE: " << ECALIsoMaxPTMultiplierEE_ << "*pT + " << ECALIsoMaxConstantEE_ << " GeV\n";
    debug_ << "Maximum HCAL isolation in EB: " << HCALIsoMaxPTMultiplierEB_ << "*pT + " << HCALIsoMaxConstantEB_ << " GeV\n";
    debug_ << "Maximum HCAL isolation in EE: " << HCALIsoMaxPTMultiplierEE_ << "*pT + " << HCALIsoMaxConstantEE_ << " GeV\n";
    debug_ << "Minimum photon ET: " << ETMin_ << " GeV\n";
    debug_ << "Fiducial region: " << fiducialRegion_ << " (";
    switch (fiducialRegion_) {
    case EB:
      debug_ << "EB";
      break;
    case EEND:
      debug_ << "EE";
      break;
    case ECAL:
      debug_ << "ECAL";
      break;
    default:
      debug_ << "invalid fiducial region";
      break;
    }
    debug_ << ")\n";
    if (useHOverE_) debug_ << "Maximum H/E: " << HOverEMax_ << endl;
    if (useSigmaEtaEta_) debug_ << "Maximum photon supercluster eta width: " << sigmaEtaEtaMax_ << endl;
    if (useTrackIso_) {
      debug_ << "Maximum photon track isolation: " << trackIsoMaxPTMultiplier_ << "*pT + " << trackIsoMaxConstant_ << " GeV\n";
    }
    debug_ << "Minimum track pT: " << trackPTMin_ << endl;
    debug_ << "Minimum dR(photon, track): " << eTrackRMin_ << endl << endl;
  }
  if ((sampleType_ != GAMMAGAMMA) && (sampleType_ != EGAMMA) && (sampleType_ != EE) && (sampleType_ != FF) && (sampleType_ != ETRACK)) {
    if (debugFlag_) debug_ << "Invalid sample type chosen.  Defaulting to \"FF\".\n\n";
    sampleType_ = FF;
  }
  if ((fiducialRegion_ != EB) && (fiducialRegion_ != EEND) && (fiducialRegion_ != ECAL)) {
    if (debugFlag_) debug_ << "Invalid fiducial region chosen.  Defaulting to \"EB\"\n\n";
    fiducialRegion_ = EB;
  }
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

  //increment total event counter
  ++numTot_;

  //pass flag
  bool pass = false;

  //get the reco::Photon collection
  Handle<PhotonCollection> pPhotons;
  bool foundPhotons = false;
  try { foundPhotons = (iEvent.getByLabel(photonTag_, pPhotons)) && (pPhotons->size() > 0); }
  catch (cms::Exception& ex) {}
  if (!foundPhotons) {
    if (debugFlag_) {
      debug_ << "No reco::Photon collection found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n\n";
    }
  }
  else {

    //loop over the photons looking for objects passing specified criteria
    PhotonCollection::const_iterator iPhoton = pPhotons->begin();
    unsigned int numCands = 0;
    unsigned int numPhotonsProcessed = 0;
    bool foundPhoton = false;
    bool foundElectron = false;
    Photon* e = NULL;
    vector<Photon*> passingCands;
    while ((iPhoton != pPhotons->end()) && (numCands < numReqdCands_)) {

      /*all samples require at least one object satisfying:
	ECAL isolation < ECALIsoMaxPTMultiplierEB_*pT + ECALIsoMaxConstantEB_ GeV (EB), 
	ECALIsoMaxPTMultiplierEE_*pT + ECALIsoMaxConstantEE_ GeV (EE)
	HCAL isolation < HCALIsoMaxPTMultiplierEB_*pT + HCALIsoMaxConstantEB_ GeV (EB), 
	HCALIsoMaxPTMultiplierEE_*pT + HCALIsoMaxConstantEE_ GeV (EE)
	H/E < HOverEMaxPresel_
	ET > ETMin_
	supercluster in fiducial region
      */
      const double pT = iPhoton->pt();
      const double ECALIso = (double)iPhoton->ecalRecHitSumEtConeDR04();
      double ECALIsoMax = ECALIsoMaxPTMultiplierEB_*pT + ECALIsoMaxConstantEB_;
      const double HCALIso = (double)iPhoton->hcalTowerSumEtConeDR04();
      double HCALIsoMax = HCALIsoMaxPTMultiplierEB_*pT + HCALIsoMaxConstantEB_;
      const unsigned int fiducialRegion = ECALFiducialRegion(const_cast<const Photon*>(&*iPhoton));
      if (fiducialRegion == EEND) {
	ECALIsoMax = ECALIsoMaxPTMultiplierEE_*pT + ECALIsoMaxConstantEE_;
	HCALIsoMax = HCALIsoMaxPTMultiplierEE_*pT + HCALIsoMaxConstantEE_;
      }
      const double HOverE = (double)iPhoton->hadronicOverEm();
      const double ET = iPhoton->et();
      const unsigned int index = iPhoton - pPhotons->begin();
      if ((ECALIso < ECALIsoMax) && (HCALIso < HCALIsoMax) && (HOverE < HOverEMaxPresel_) && (ET > ETMin_) && 
	  ((fiducialRegion_ == ECAL) || (fiducialRegion_ == fiducialRegion))) {
	if (debugFlag_) {
	  debug_ << "Run: " << runNum << endl;
	  debug_ << "Event: " << evtNum << endl;
	  debug_ << "Lumi section: " << lumiNum << endl;
	  debug_ << "Photon index: " << index << endl;
	  debug_ << "Photon pT: " << pT << " GeV\n";
	  debug_ << "Photon ECAL isolation: " << ECALIso << " GeV\n";
	  debug_ << "Photon HCAL isolation: " << HCALIso << " GeV\n";
	  debug_ << "Photon H/E: " << HOverE << " GeV\n";
	  debug_ << "Photon ET: " << ET << " GeV\n";
	  debug_ << "Photon fiducial region: " << fiducialRegion << " (";
	  if (fiducialRegion == EB) debug_ << "EB";
	  else debug_ << "EE";
	  debug_ << ")\n\n";
	}

	//decide if the candidate passes based on the photon ID cuts and flags
	const double sigmaEtaEta = (double)iPhoton->sigmaEtaEta();
	const double trackIso = (double)iPhoton->trkSumPtHollowConeDR04();
	const double trackIsoMax = trackIsoMaxPTMultiplier_*pT + trackIsoMaxConstant_;
	bool passedPhotonID = false;
	if ((useHOverE_ && useSigmaEtaEta_ && useTrackIso_) && 
	    ((sigmaEtaEta < sigmaEtaEtaMax_) && (HOverE < HOverEMax_) && (trackIso < trackIsoMax))) passedPhotonID = true;
	else if ((useHOverE_ && useSigmaEtaEta_ && (!useTrackIso_)) && 
		 ((sigmaEtaEta < sigmaEtaEtaMax_) && (HOverE < HOverEMax_))) passedPhotonID = true;
	else if ((useHOverE_ && (!useSigmaEtaEta_) && useTrackIso_) && 
		 ((HOverE < HOverEMax_) && (trackIso < trackIsoMax))) passedPhotonID = true;
	else if (((!useHOverE_) && useSigmaEtaEta_ && useTrackIso_) && 
		 ((sigmaEtaEta < sigmaEtaEtaMax_) && (trackIso < trackIsoMax))) passedPhotonID = true;
	else if ((useHOverE_ && (!useSigmaEtaEta_) && (!useTrackIso_)) && (HOverE < HOverEMax_)) passedPhotonID = true;
	else if (((!useHOverE_) && useSigmaEtaEta_ && (!useTrackIso_)) && (sigmaEtaEta < sigmaEtaEtaMax_)) passedPhotonID = true;
	else if (((!useHOverE_) && (!useSigmaEtaEta_) && useTrackIso_) && (trackIso < trackIsoMax)) passedPhotonID = true;
	else if ((!useHOverE_) && (!useSigmaEtaEta_) && (!useTrackIso_)) passedPhotonID = true;

	/*gammagamma, egamma, ee, and etrack samples require:
	 sigmaEtaEta < sigmaEtaEtaMax_
	 H/E < HOverEMax_
	 track isolation < trackIsoMaxPTMultiplier_*pT + trackIsoMaxConstant_ GeV
	*/
	if (((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA) || (sampleType_ == EE) || (sampleType_ == ETRACK)) && passedPhotonID) {
	  if (debugFlag_) {
	    debug_ << "Photon supercluster eta width: " << sigmaEtaEta << endl;
	    debug_ << "Photon track isolation: " << trackIso << endl;
	  }

	  //egamma, ee, and etrack samples require the candidate to have at least one pixel seed
	  if (iPhoton->hasPixelSeed()) {

	    /*ee and etrack samples require only electrons
	      for the etrack sample, we also need to save the candidate*/
	    if ((sampleType_ == EE) || (sampleType_ == ETRACK)) {
	      ++numCands;
	      if (debugFlag_) debug_ << "Number of found EE/ETRACK candidates: " << numCands << endl;
	      e = const_cast<Photon*>(&*iPhoton);
	      passingCands.push_back(const_cast<Photon*>(&*iPhoton));
	    }

	    //egamma sample counts this electron only if another one wasn't previously counted
	    if ((sampleType_ == EGAMMA) && (!foundElectron)) {
	      foundElectron = true;
	      ++numCands;
	      passingCands.push_back(const_cast<Photon*>(&*iPhoton));
	      if (debugFlag_) {
		debug_ << "Found electron.\n";
		debug_ << "Number of found EGAMMA candidates: " << numCands << endl;
	      }
	    }
	  }

	  //gammagamma and egamma samples require no pixel seed
	  else {

	    //gammagamma sample requires only photons
	    if (sampleType_ == GAMMAGAMMA) {
	      ++numCands;
	      passingCands.push_back(const_cast<Photon*>(&*iPhoton));
	      if (debugFlag_) debug_ << "Number of found GAMMAGAMMA candidates: " << numCands << endl;
	    }

	    //egamma sample counts this photons only if another one wasn't previously counted
	    if ((sampleType_ == EGAMMA) && (!foundPhoton)) {
	      foundPhoton = true;
	      ++numCands;
	      passingCands.push_back(const_cast<Photon*>(&*iPhoton));
	      if (debugFlag_) {
		debug_ << "Found photon.\n";
		debug_ << "Number of found EGAMMA candidates: " << numCands << endl;
	      }
	    }
	  }
	}

	/*ff sample requires:
	 sigmaEtaEta >= sigmaEtaEtaMax_ or
	 H/E >= HOverEMax_ or
	 track isolation >= trackIsoMaxPTMultiplier_*pT + trackIsoMaxConstant_ GeV

	 NB. No requirement on number of pixel seeds in the fake definition!
	*/
	//HACK 2-Jun-10: require pixel seed until we figure out if we care about it or not
	else if ((sampleType_ == FF) && (!passedPhotonID) && (iPhoton->hasPixelSeed())) {
	  ++numCands;
	  passingCands.push_back(const_cast<Photon*>(&*iPhoton));
	  if (debugFlag_) {
	    debug_ << "Number of found FF candidates: " << numCands << endl;
	    debug_ << "Photon supercluster eta width: " << sigmaEtaEta << endl;
	    debug_ << "Photon track isolation: " << trackIso << endl;
	  }
	}
      }

      //advance to the next photon in the collection
      ++iPhoton;
      ++numPhotonsProcessed;
      if (debugFlag_) debug_ << endl;
    }
    if (debugFlag_) debug_ << "Number of photons processed: " << numPhotonsProcessed << "/" << pPhotons->size() << endl << endl;

    //no matter the sample, no event can pass if the minimum number of required candidates wasn't found
    if (numCands == numReqdCands_) {
      if (debugFlag_) debug_ << "Total number of found candidates: " << numCands << endl;

      //get the HB/HE RecHits for calculating beam halo tags
      Handle<HBHERecHitCollection> pHBHERecHits;
      bool foundHBHERecHits = false;
      bool HEHalo = false; //if no HB/HE RecHits are available in the event, the HE halo tag stays false
      bool muonHalo = false; //if no cosmic muon tracks are available in the event, the muon halo tag stays false
      try { foundHBHERecHits = (iEvent.getByLabel(HBHERecHitTag_, pHBHERecHits)) && (pHBHERecHits->size() > 0); }
      catch (cms::Exception& ex) {}
      if (!foundHBHERecHits) {
	if (debugFlag_) {
	  debug_ << "No reco::HBHERecHit collection found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n\n";
	}
      }
      else {

	//get the HCAL geometry for calculating the rho and phi of the HB/HE RecHit (there's got to be a better way to do this!)
	iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle_);
	const CaloGeometry* pGeometry = caloGeometryHandle_.product();

	//calculate the HE tag
	vector<Photon*>::const_iterator iPassingCand = passingCands.begin();
	while ((iPassingCand != passingCands.end()) && (!HEHalo)) {
	  HBHERecHitCollection::const_iterator iHBHERecHit = pHBHERecHits->begin();
	  while ((iHBHERecHit != pHBHERecHits->end()) && (!HEHalo)) {
	    const DetId HBHERecHitID = iHBHERecHit->detid();
	    const float rho = pGeometry->getPosition(HBHERecHitID).perp();
	    if (debugFlag_) debug_ << "HB/HE RecHit rho: " << rho << " cm\n";
	    const double phi = (double)pGeometry->getPosition(HBHERecHitID).phi();
	    const double E = iHBHERecHit->energy();
	    const double dPhiPhotonHBHERecHit = dPhi((*iPassingCand)->phi(), phi);
	    if ((iHBHERecHit->id().subdet() == HcalEndcap) && (E > 1/*GeV*/) && (rho < 130/*cm*/) && (rho > 115/*cm*/) && (dPhiPhotonHBHERecHit <= 0.2)) {
	      if (debugFlag_) {
		debug_ << "Found HE halo candidate:\n";
		debug_ << "DetID: " << HBHERecHitID.rawId() << endl;
		debug_ << "Rho: " << rho << " cm\n";
		debug_ << "Energy: " << E << " GeV\n";
		debug_ << "dPhi: " << dPhiPhotonHBHERecHit << endl << endl;
	      }
	      HEHalo = true;
	    }
	    ++iHBHERecHit;
	  }
	  ++iPassingCand;
	}
      }

      //only calculate the muon halo tag if the HB/HE tag is false to save time
      if (!HEHalo) {

	//get the cosmic tracks for calculating beam halo tags
	Handle<TrackCollection> pCosmicTracks;
	bool foundCosmicTracks = false;
	try { foundCosmicTracks = (iEvent.getByLabel(cosmicTrackTag_, pCosmicTracks)) && (pCosmicTracks->size() > 0); }
	catch (cms::Exception& ex) {}
	if (!foundCosmicTracks) {
	  if (debugFlag_) {
	    debug_ << "No reco::Track collection for cosmic tracks found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum;
	    debug_ << ".\n\n";
	  }
	}
	else {

	  //calculate the muon tag, using only the innermost hit position of the track as per Andrew's e-mail on 6-Jun-10
	  vector<Photon*>::const_iterator iPassingCand = passingCands.begin();
	  while ((iPassingCand != passingCands.end()) && (!muonHalo)) {
	    TrackCollection::const_iterator iCosmicTrack = pCosmicTracks->begin();
	    while ((iCosmicTrack != pCosmicTracks->end()) && (!muonHalo)) {
	      const double rhoInner = sqrt(iCosmicTrack->innerPosition().perp2());
	      const double dPhiInner = dPhi((*iPassingCand)->phi(), iCosmicTrack->innerPosition().phi());
	      if (debugFlag_) {
		debug_ << "Cosmic track innermost hit rho: " << rhoInner << endl;
		debug_ << "Cosmic track innermost hit phi: " << (iCosmicTrack->innerPosition().phi()) << endl;
	      }
	      bool isCSC = true;
	      CSCDetId CSCID;
	      try { CSCID = CSCDetId(iCosmicTrack->innerDetId()); }
	      catch (cms::Exception& ex) { isCSC = false; }
	      bool isRPCEndcap = true;
	      if (!isCSC) {
		RPCDetId RPCID;
		try { RPCID = RPCDetId(iCosmicTrack->innerDetId()); }
		catch (cms::Exception& ex) { isRPCEndcap = false; }
		if ((isRPCEndcap) && (RPCID.region() == 0)) isRPCEndcap = false; //region 0 is the barrel
	      }
	      else isRPCEndcap = false;
	      if ((isCSC || isRPCEndcap) && (rhoInner < 170/*cm*/) && (rhoInner > 115/*cm*/) && (dPhiInner <= 0.2)) {
		if (debugFlag_) {
		  debug_ << "Found muon halo candidate:\n";
		  debug_ << "DetID: " << iCosmicTrack->innerDetId() << endl;
		  debug_ << "Rho: " << rhoInner << " cm\n";
		  debug_ << "dPhi: " << dPhiInner << endl << endl;
		}
		muonHalo = true;
	      }
	      ++iCosmicTrack;
	    }
	    ++iPassingCand;
	  }
	}
      }

      //no matter the sample, no event can pass if it passes the halo tags
      if ((!HEHalo) && (!muonHalo)) {

	//gammagamma, egamma, ee, and ff samples require separation between the candidates of 0.8
	if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA) || (sampleType_ == EE) || (sampleType_ == FF)) {
	  Photon* photon1 = *(passingCands.begin());
	  Photon* photon2 = *(passingCands.begin() + 1);
	  const double dRPhoton1Photon2 = dR(photon1->eta(), photon2->eta(), photon1->phi(), photon2->phi());
	  if (debugFlag_) debug_ << "dR between candidates: " << dRPhoton1Photon2 << endl;
	  if (dRPhoton1Photon2 > 0.8) {
	    pass = true;
	    if (debugFlag_) {
	      debug_ << "Event passes.\n";
	      debug_ << "Run: " << runNum << ", event: " << evtNum << ", lumi section: " << lumiNum << endl;
	    }
	  }
	  else if (debugFlag_) debug_ << "Event fails the dR requirement.\n";
	}

	//etrack sample requires a track
	else {

	  //get the GeneralTrack collection
	  Handle<TrackCollection> pTracks;
	  bool tracksFound = false;
	  try { tracksFound = (iEvent.getByLabel(trackTag_, pTracks)) && (pTracks->size() > 0); }
	  catch (cms::Exception& ex) {}
	  if (!tracksFound) {
	    if (debugFlag_) {
	      debug_ << "No reco::Track collection for general tracks found in run " << runNum << ", event " << evtNum << ".\n";
	    }
	  }
	  else {
	    TrackCollection::const_iterator iTrack = pTracks->begin();
	    bool trackFound = false;

	    /*etrack sample requires:
	      track pT > trackPTMin_
	      dR(track, electron) > eTrackRMin_
	    */
	    while ((iTrack != pTracks->end()) && (!trackFound)) {
	      double dRETrack = dR(e->eta(), iTrack->eta(), e->phi(), iTrack->phi());
	      double pT = iTrack->pt();
	      if ((pT > trackPTMin_) && (dRETrack > eTrackRMin_)) {
		trackFound = true;
		if (debugFlag_) {
		  debug_ << "Track pT: " << pT << " GeV\n";
		  debug_ << "dR(track, photon): " << dRETrack << endl;
		}
	      }
	      else {
		++iTrack;
		if (debugFlag_) debug_ << "This track doesn't meet the criteria; advancing to the next track.\n";
	      }
	    }

	    //if the track was found, pass the event; otherwise, fail the event
	    pass = trackFound;
	    if (debugFlag_) {
	      if (trackFound) {
		debug_ << "Event passes.\n";
		debug_ << "Run: " << runNum << ", event: " << evtNum << ", lumi section: " << lumiNum << endl;
	      }
	      else debug_ << "Event fails.\n";
	    }
	  }
	}
      }
    }
  }

  //filter the event
  if (debugFlag_) debug_ << endl;
  if (pass) ++numPassing_;
  return pass;
}
  
// ------------ method called once each job just before starting event loop  ------------
void 
SampleMaker::beginJob()
{
  //set the minimum number of candidate objects required to be found for each event
  numReqdCands_ = 2;
  if (sampleType_ == ETRACK) numReqdCands_ = 1;
  if (debugFlag_) debug_ << "Number of required candidates: " << numReqdCands_ << endl;

  //initialize event counters
  numTot_ = 0;
  numPassing_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SampleMaker::endJob() {
  cout << "Number of events in " << sampleTypeString() << " sample: " << numPassing_ << "/" << numTot_ << " = ";
  cout << ((float)numPassing_/(float)numTot_);
  cout << endl;
  if (debugFlag_) {
    debug_ << "Number of events in " << sampleTypeString() << " sample: " << numPassing_ << "/" << numTot_ << " = ";
    debug_ << ((float)numPassing_/(float)numTot_) << endl;
    debug_.close();
  }
}

//get the ECAL fiducial region in units this code understands
const unsigned int SampleMaker::ECALFiducialRegion(const Photon* photon) const
{
  unsigned int reg = -1;
  if (!photon->isEBEEGap()) {
    if (photon->isEB()) reg = EB;
    if (photon->isEE()) reg = EEND;
  }
  return reg;
}

//get the sample type as a string
const string SampleMaker::sampleTypeString() const
{
  string sampleType = "";
  switch (sampleType_) {
  case GAMMAGAMMA:
    sampleType = "GAMMAGAMMA";
    break;
  case EGAMMA:
    sampleType = "EGAMMA";
    break;
  case EE:
    sampleType = "EE";
    break;
  case FF:
    sampleType = "FF";
    break;
  case ETRACK:
    sampleType = "ETRACK";
    break;
  default:
    sampleType = "invalid sample type";
    break;
  }
  return sampleType;
}

//determine if the EM object passes or fails preselection criteria
/*const bool SampleMaker::passesPreselection(const Photon* pPhoton) const
  {*/
  /*all samples require at least one object satisfying:
    ECAL isolation < ECALIsoMaxPTMultiplierEB_*pT + ECALIsoMaxConstantEB_ GeV (EB), 
    ECALIsoMaxPTMultiplierEE_*pT + ECALIsoMaxConstantEE_ GeV (EE)
    HCAL isolation < HCALIsoMaxPTMultiplierEB_*pT + HCALIsoMaxConstantEB_ GeV (EB), 
    HCALIsoMaxPTMultiplierEE_*pT + HCALIsoMaxConstantEE_ GeV (EE)
    H/E < HOverEMaxPresel_
    ET > ETMin_
    supercluster in fiducial region
  */
  /*bool pass = false;
  const double pT = pPhoton->pt();
  const double ECALIso = (double)pPhoton->ecalRecHitSumEtConeDR04();
  double ECALIsoMax = ECALIsoMaxPTMultiplierEB_*pT + ECALIsoMaxConstantEB_;
  const double HCALIso = (double)pPhoton->hcalTowerSumEtConeDR04();
  double HCALIsoMax = HCALIsoMaxPTMultiplierEB_*pT + HCALIsoMaxConstantEB_;
  const unsigned int fiducialRegion = ECALFiducialRegion(pPhoton);
  if (fiducialRegion == EEND) {
    ECALIsoMax = ECALIsoMaxPTMultiplierEE_*pT + ECALIsoMaxConstantEE_;
    HCALIsoMax = HCALIsoMaxPTMultiplierEE_*pT + HCALIsoMaxConstantEE_;
  }
  const double HOverE = (double)pPhoton->hadronicOverEm();
  const double ET = pPhoton->et();
  const unsigned int index = pPhoton - pPhotons->begin();
  if ((ECALIso < ECALIsoMax) && (HCALIso < HCALIsoMax) && (HOverE < HOverEMaxPresel_) && (ET > ETMin_) && 
      ((fiducialRegion_ == ECAL) || (fiducialRegion_ == fiducialRegion))) pass = true;
  return pass;
  }*/

//define this as a plug-in
DEFINE_FWK_MODULE(SampleMaker);
