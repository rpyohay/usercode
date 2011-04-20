// -*- C++ -*-
//
// Package:    Producers
// Class:      EDMCategoryProducer
// 
/**\class EDMCategoryProducer EDMCategoryProducer.cc GMSBTools/Producers/src/EDMCategoryProducer.cc

 Description: produce information about the event category and store it in the EDM event

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Mon Apr 18 18:28:36 CEST 2011
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
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "GMSBTools/Filters/interface/Categorizer.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"


//
// class declaration
//

class EDMCategoryProducer : public edm::EDProducer {
   public:
      explicit EDMCategoryProducer(const edm::ParameterSet&);
      ~EDMCategoryProducer();

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

  //return the DetId of the seed crystal
  DetId seedID(const reco::Photon*) const;

  //return the right RecHit collection based on the given DetId
  const edm::Handle<EcalRecHitCollection>& 
  recHitCollectionForHit(const DetId&, const edm::Handle<EBRecHitCollection>&, 
			 const edm::Handle<EERecHitCollection>&) const;

  //return the time of the RecHit with given DetId
  float recHitTime(const DetId&, const edm::Handle<EcalRecHitCollection>&) const;

  //return E2/E9 of the RecHit with the given DetId
  float recHitE2OverE9(const DetId&, CaloNavigator<DetId>&, 
		       const edm::Handle<EcalRecHitCollection>&) const;

  //return the seed time of the given photon
  float seedTime(const reco::Photon*, const edm::Handle<EBRecHitCollection>&, 
		 const edm::Handle<EERecHitCollection>&) const;

  //return E2/E9 of the given photon
  float e2OverE9(const reco::Photon*, edm::ESHandle<CaloTopology>&, const edm::Handle<EBRecHitCollection>&, 
  const edm::Handle<EERecHitCollection>&) const;

  //put the edm::ValueMap<bool> product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, const VBOOL&, edm::Event&, 
			   const STRING&) const;

  //put the edm::ValueMap<int> product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, const VINT&, edm::Event&, 
			   const STRING&) const;

  //put int product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, const int, edm::Event&, 
			   const STRING&) const;

  //put bool product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, const bool, edm::Event&, 
			   const STRING&) const;
      
      // ----------member data ---------------------------

  //input
  edm::InputTag photonTag_;
  edm::InputTag recHitTagEB_;
  edm::InputTag recHitTagEE_;

  //cuts
  double photonETMin_;
  double photonAbsEtaMax_;
  double photonECALIsoMaxPTMultiplier_;
  double photonECALIsoMaxConstant_;
  double photonHCALIsoMaxPTMultiplier_;
  double photonHCALIsoMaxConstant_;
  double photonHOverEMax_;
  double photonTrackIsoMaxPTMultiplier_;
  double photonTrackIsoMaxPTConstant_;
  double photonSigmaIetaIetaMax_;
  double photonAbsSeedTimeMax_;
  double photonE2OverE9Max_;
  double photonDPhiMin_;
  VUINT channelStatuses_;
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
EDMCategoryProducer::EDMCategoryProducer(const edm::ParameterSet& iConfig) :

  //input
  photonTag_(iConfig.getUntrackedParameter<edm::InputTag>("photonTag", 
							  edm::InputTag("photons", "", 
									"RECOCleaned"))),
  recHitTagEB_(iConfig.getUntrackedParameter<edm::InputTag>("recHitTagEB", 
							    edm::InputTag("reducedEcalRecHitsEB", 
									  "", "RECOCleaned"))),
  recHitTagEE_(iConfig.getUntrackedParameter<edm::InputTag>("recHitTagEE", 
							    edm::InputTag("reducedEcalRecHitsEE", 
									  "", "RECOCleaned"))),

  //cuts
  photonETMin_(iConfig.getUntrackedParameter<double>("photonETMin", 30.0/*GeV*/)),
  photonAbsEtaMax_(iConfig.getUntrackedParameter<double>("photonAbsEtaMax", 1.379)),
  photonECALIsoMaxPTMultiplier_(iConfig.getUntrackedParameter<double>
				("photonECALIsoMaxPTMultiplier", 0.006)),
  photonECALIsoMaxConstant_(iConfig.getUntrackedParameter<double>
			    ("photonECALIsoMaxConstant", 4.2/*GeV*/)),
  photonHCALIsoMaxPTMultiplier_(iConfig.getUntrackedParameter<double>
				("photonHCALIsoMaxPTMultiplier", 0.0025)),
  photonHCALIsoMaxConstant_(iConfig.getUntrackedParameter<double>
			    ("photonHCALIsoMaxConstant", 2.2/*GeV*/)),
  photonHOverEMax_(iConfig.getUntrackedParameter<double>("photonHOverEMax", 0.05)),
  photonTrackIsoMaxPTMultiplier_(iConfig.getUntrackedParameter<double>
				 ("photonTrackIsoMaxPTMultiplier", 0.001)),
  photonTrackIsoMaxPTConstant_(iConfig.getUntrackedParameter<double>
			       ("photonTrackIsoMaxPTConstant", 2.0/*GeV*/)),
  photonSigmaIetaIetaMax_(iConfig.getUntrackedParameter<double>("photonSigmaIetaIetaMax", 0.013)),
  photonAbsSeedTimeMax_(iConfig.getUntrackedParameter<double>("photonAbsSeedTimeMax", 3.0/*ns*/)),
  photonE2OverE9Max_(iConfig.getUntrackedParameter<double>("photonE2OverE9Max", 0.95)),
  photonDPhiMin_(iConfig.getUntrackedParameter<double>("photonDPhiMin", 0.05)),
  channelStatuses_(iConfig.getUntrackedParameter<VUINT>("channelStatuses", VUINT(EcalRecHit::kGood, 1)))
{
   //register your products
  produces<edm::ValueMap<unsigned int> >("photonType");
  produces<unsigned int>("eventCategory");
  produces<bool>("passDPhiMin");
  produces<edm::ValueMap<bool> >("passETMin");
  produces<edm::ValueMap<bool> >("passAbsEtaMax");
  produces<edm::ValueMap<bool> >("passECALIsoMax");
  produces<edm::ValueMap<bool> >("passHCALIsoMax");
  produces<edm::ValueMap<bool> >("passHOverEMax");
  produces<edm::ValueMap<bool> >("passTrackIsoMax");
  produces<edm::ValueMap<bool> >("passSigmaIetaIetaMax");
  produces<edm::ValueMap<bool> >("passAbsSeedTimeMax");
  produces<edm::ValueMap<bool> >("passE2OverE9Max");
  produces<edm::ValueMap<bool> >("hasPixelSeed");

   //now do what ever other initialization is needed
  
}


EDMCategoryProducer::~EDMCategoryProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EDMCategoryProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //photon variables to pass to the Categorizer object
  VDOUBLE photonET;
  VDOUBLE photonEta;
  VDOUBLE photonECALIso;
  VDOUBLE photonHCALIso;
  VDOUBLE photonHOverE;
  VDOUBLE photonTrackIso;
  VDOUBLE photonSigmaIetaIeta;
  VDOUBLE photonSeedTime;
  VDOUBLE photonE2OverE9;
  VDOUBLE photonPhi;
  VBOOL photonHasPixelSeed;
  Categorizer categorizer(photonET, photonEta, photonECALIso, photonHCALIso, photonHOverE, 
			  photonTrackIso, photonSigmaIetaIeta, photonSeedTime, photonE2OverE9, 
			  photonPhi, photonHasPixelSeed, photonETMin_, photonAbsEtaMax_, 
			  photonECALIsoMaxPTMultiplier_, photonECALIsoMaxConstant_, 
			  photonHCALIsoMaxPTMultiplier_, photonHCALIsoMaxConstant_, 
			  photonHOverEMax_, photonTrackIsoMaxPTMultiplier_,
			  photonTrackIsoMaxPTConstant_, photonSigmaIetaIetaMax_, 
			  photonAbsSeedTimeMax_, photonE2OverE9Max_, photonDPhiMin_);
  //loop over photons
  edm::Handle<reco::PhotonCollection> pPhotons;
  edm::Handle<EBRecHitCollection> pRecHitsEB;
  edm::Handle<EERecHitCollection> pRecHitsEE;
  edm::ESHandle<CaloTopology> caloTopology;
  iSetup.get<CaloTopologyRecord>().get(caloTopology);
  if (getCollection_(pPhotons, photonTag_, iEvent) && 
      getCollection_(pRecHitsEB, recHitTagEB_, iEvent) && 
      getCollection_(pRecHitsEE, recHitTagEE_, iEvent)) {
    for (reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin(); 
	 iPhoton != pPhotons->end(); ++iPhoton) {

      //fill photon variable vectors
      photonET.push_back(iPhoton->et());
      photonEta.push_back(iPhoton->eta());
      photonECALIso.push_back(iPhoton->ecalRecHitSumEtConeDR04());
      photonHCALIso.push_back(iPhoton->hcalTowerSumEtConeDR04());
      photonHOverE.push_back(iPhoton->hadronicOverEm());
      photonTrackIso.push_back(iPhoton->trkSumPtHollowConeDR04());
      photonSigmaIetaIeta.push_back(iPhoton->sigmaIetaIeta());
      photonSeedTime.push_back(seedTime(const_cast<const reco::Photon*>(&*iPhoton), pRecHitsEB, 
					pRecHitsEE));
      photonE2OverE9.push_back(e2OverE9(const_cast<const reco::Photon*>(&*iPhoton), caloTopology, pRecHitsEB, pRecHitsEE));
      photonPhi.push_back(iPhoton->phi());
      photonHasPixelSeed.push_back(iPhoton->hasPixelSeed());
    }

    //categorize the event
    categorizer.setPhotonET(photonET);
    categorizer.setPhotonEta(photonEta);
    categorizer.setPhotonECALIso(photonECALIso);
    categorizer.setPhotonHCALIso(photonHCALIso);
    categorizer.setPhotonHOverE(photonHOverE);
    categorizer.setPhotonTrackIso(photonTrackIso);
    categorizer.setPhotonSigmaIetaIeta(photonSigmaIetaIeta);
    /*categorizer.setPhotonSeedTime(photonSeedTime);
      categorizer.setPhotonE2OverE9(photonE2OverE9);*/
    categorizer.setPhotonPhi(photonPhi);
    categorizer.setPhotonHasPixelSeed(photonHasPixelSeed);
    try {
      categorizer.decideAll();
      categorizer.findPassingPhotons();
      categorizer.classify();
    }
    catch (STRING& badInput) { throw cms::Exception("EDMCategoryProducer") << badInput; }
  }

  //write the category information to the event
  try {
    putProductIntoEvent(pPhotons, categorizer.getPassingPhotonType(), iEvent, "photonType");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassETMin(), iEvent, "passETMin");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassAbsEtaMax(), iEvent, "passAbsEtaMax");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassECALIsoMax(), iEvent, 
			"passECALIsoMax");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassHCALIsoMax(), iEvent, 
			"passHCALIsoMax");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassHOverEMax(), iEvent, "passHOverEMax");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassTrackIsoMax(), iEvent, 
			"passTrackIsoMax");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassSigmaIetaIetaMax(), iEvent, 
			"passSigmaIetaIetaMax");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassAbsSeedTimeMax(), iEvent, 
			"passAbsSeedTimeMax");
    putProductIntoEvent(pPhotons, categorizer.getPhotonPassE2OverE9Max(), iEvent, 
			"passE2OverE9Max");
    putProductIntoEvent(pPhotons, categorizer.getPhotonHasPixelSeed(), iEvent, "hasPixelSeed");
    putProductIntoEvent(pPhotons, categorizer.getCategory(), iEvent, "eventCategory");
    putProductIntoEvent(pPhotons, categorizer.getEvtPassDPhiMin(), iEvent, "passDPhiMin");
  }
  catch (STRING& badInput) { throw cms::Exception("EDMCategoryProducer") << badInput; }
}

// ------------ method called once each job just before starting event loop  ------------
void 
EDMCategoryProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EDMCategoryProducer::endJob() {
}

DetId EDMCategoryProducer::seedID(const reco::Photon* photon) const
{
  return photon->superCluster()->seed()->seed();
}

const edm::Handle<EcalRecHitCollection>& 
EDMCategoryProducer::recHitCollectionForHit(const DetId& ID, 
					    const edm::Handle<EBRecHitCollection>& pRecHitsEB, 
					    const edm::Handle<EERecHitCollection>& pRecHitsEE) 
  const
{
  if (ID.det() == DetId::Ecal) {
    if (ID.subdetId() == EcalBarrel) return pRecHitsEB;
    else if (ID.subdetId() == EcalEndcap) return pRecHitsEE;
    else throw "Error: subdetector ID not EcalBarrel or EcalEndcap.\n";
  }
  else throw "Error: Detector ID not DetId::Ecal.\n";
}

float EDMCategoryProducer::recHitTime(const DetId& ID, 
				      const edm::Handle<EcalRecHitCollection>& pRecHits) const
{
  float recHitTime = -1.0;
  EcalRecHitCollection::const_iterator hit = pRecHits->find(ID);
  if ((hit != pRecHits->end()) && hit->isTimeValid()) recHitTime = fabs(hit->time());
  return recHitTime;
}

float EDMCategoryProducer::recHitE2OverE9(const DetId& ID, CaloNavigator<DetId>& cursor, 
					  const edm::Handle<EcalRecHitCollection>& pRecHits) const
{
  float highestE = 0.0;
  float secondHighestE = 0.0;
  float E3x3 = 0.0;
  float e2OverE9 = -1.0;
  bool missingHit = false;
  int i = -1;
  int j = -1;
  while ((i <= 1) && !missingHit) {
    while ((j <= 1) && !missingHit) {
      EcalRecHitCollection::const_iterator hit = pRecHits->find(cursor.offsetBy(i, j));
      if (hit == pRecHits->end()) missingHit = true;
      else {
	uint32_t recoFlag = hit->recoFlag();
	bool useChannel = false;
	std::vector<unsigned int>::const_iterator iFlag = channelStatuses_.begin();
	while ((iFlag != channelStatuses_.end()) && !useChannel) {
	  useChannel = useChannel || (recoFlag == *iFlag);
	  ++iFlag;
	}
	if (useChannel) {
	  float hitE = 0.0;
	  if (recoFlag == EcalRecHit::kOutOfTime) hitE = hit->outOfTimeEnergy();
	  else hitE = hit->energy();
	  if (hitE > highestE) {
	    secondHighestE = highestE;
	    highestE = hitE;
	  }
	  else if (hitE > secondHighestE) secondHighestE = hitE;
	  E3x3+=hitE;
	}
	else missingHit = true;
      }
      ++j;
    }
    ++i;
  }
  if (!missingHit) e2OverE9 = (highestE + secondHighestE)/E3x3;
  return e2OverE9;
}

float EDMCategoryProducer::seedTime(const reco::Photon* photon, 
				    const edm::Handle<EBRecHitCollection>& pRecHitsEB, 
				    const edm::Handle<EERecHitCollection>& pRecHitsEE) const
{
  float theSeedTime = -1.0;
  const DetId theSeedID = seedID(photon);
  try {
    theSeedTime = recHitTime(theSeedID, recHitCollectionForHit(theSeedID, pRecHitsEB, 
							       pRecHitsEE));
  }
  catch (std::string& notECAL) { std::cerr << notECAL; }
  return theSeedTime;
}

float EDMCategoryProducer::e2OverE9(const reco::Photon* photon, edm::ESHandle<CaloTopology>& caloTopology, 
				    const edm::Handle<EBRecHitCollection>& pRecHitsEB, 
				    const edm::Handle<EERecHitCollection>& pRecHitsEE) const
{
  float theE2OverE9 = -1.0;
  const DetId theSeedID = seedID(photon);
  CaloNavigator<DetId> cursor = 
    CaloNavigator<DetId>(theSeedID, caloTopology->getSubdetectorTopology(theSeedID));
  try {
    theE2OverE9 = recHitE2OverE9(theSeedID, cursor, 
				 recHitCollectionForHit(theSeedID, pRecHitsEB, pRecHitsEE));
  }
  catch (std::string& notECAL) { std::cerr << notECAL; }
  return theE2OverE9;
}

void 
EDMCategoryProducer::putProductIntoEvent(const edm::Handle<reco::PhotonCollection>& pPhotons, 
					 const VBOOL& vals, edm::Event& iEvent, 
					 const STRING& label) const
{
  std::auto_ptr<edm::ValueMap<bool> > out(new edm::ValueMap<bool>());
  edm::ValueMap<bool>::Filler filler(*out);
  filler.insert(pPhotons, vals.begin(), vals.end());
  filler.fill();
  iEvent.put(out, label);
}

void 
EDMCategoryProducer::putProductIntoEvent(const edm::Handle<reco::PhotonCollection>& pPhotons, 
					 const VINT& vals, edm::Event& iEvent, 
					 const STRING& label) const
{
  std::auto_ptr<edm::ValueMap<int> > out(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler(*out);
  filler.insert(pPhotons, vals.begin(), vals.end());
  filler.fill();
  iEvent.put(out, label);
}

void 
EDMCategoryProducer::putProductIntoEvent(const edm::Handle<reco::PhotonCollection>& pPhotons, 
					 const int vals, edm::Event& iEvent, 
					 const STRING& label) const
{
  std::auto_ptr<int> out(new int());
  *out = vals;
  iEvent.put(out, label);
}

void 
EDMCategoryProducer::putProductIntoEvent(const edm::Handle<reco::PhotonCollection>& pPhotons, 
					 const bool vals, edm::Event& iEvent, 
					 const STRING& label) const
{
  std::auto_ptr<bool> out(new bool());
  *out = vals;
  iEvent.put(out, label);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EDMCategoryProducer);
