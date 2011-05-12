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
#include "DataFormats/EcalDetId/interface/EBDetId.h"


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
      edm::LogInfo("Error") << err.str();
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

  //return the seed ieta of the given photon if it is in EB
  int seedIeta(const reco::Photon*) const;

  //return E2/E9 of the given photon
  float e2OverE9(const reco::Photon*, edm::ESHandle<CaloTopology>&, 
		 const edm::Handle<EBRecHitCollection>&, 
		 const edm::Handle<EERecHitCollection>&) const;

  //put the edm::ValueMap<bool> product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, const VBOOL&, edm::Event&, 
			   const STRING&) const;

  //put the edm::ValueMap<int> product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, const VINT&, edm::Event&, 
			   const STRING&) const;

  //put the edm::ValueMap<double> product into the event
  void putProductIntoEvent(const edm::Handle<reco::PhotonCollection>&, const VDOUBLE&, 
			   edm::Event&, const STRING&) const;

  //put int product into the event
  void putProductIntoEvent(const int, edm::Event&, const STRING&) const;

  //put bool product into the event
  void putProductIntoEvent(const bool, edm::Event&, const STRING&) const;

  //put double product into the event
  void putProductIntoEvent(const double, edm::Event&, const STRING&) const;

  //print debug info
  void debugPrint(const Categorizer&, edm::Event&) const;
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
  channelStatuses_(iConfig.getUntrackedParameter<VUINT>("channelStatuses", 
							VUINT(EcalRecHit::kGood, 1)))
{
   //register your products
  produces<edm::ValueMap<int> >("photonType");
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
  produces<edm::ValueMap<bool> >("passingPhotons");
  produces<int>("eventCategory");
  produces<bool>("passDPhiMin");
  produces<double>("evtDiEMET");
  produces<double>("evtInvMass");
  produces<edm::ValueMap<double> >("photonSeedTime");
  produces<edm::ValueMap<double> >("photonE2OverE9");
  produces<edm::ValueMap<int> >("photonSeedIeta");

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
  VINT photonSeedIeta;
  VDOUBLE photonPhi;
  VBOOL photonHasPixelSeed;
  VBOOL photonPasses;
  math::XYZTLorentzVector diEMP4;
  unsigned int passingPhotonCount = 0;
  double evtInvMass = -1.0;
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
      photonE2OverE9.push_back(e2OverE9(const_cast<const reco::Photon*>(&*iPhoton), caloTopology, 
					pRecHitsEB, pRecHitsEE));
      photonSeedIeta.push_back(seedIeta(const_cast<const reco::Photon*>(&*iPhoton)));
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
    categorizer.setPhotonSeedTime(photonSeedTime);
    categorizer.setPhotonE2OverE9(photonE2OverE9);
    categorizer.setPhotonPhi(photonPhi);
    categorizer.setPhotonHasPixelSeed(photonHasPixelSeed);
    try {
      categorizer.decideAll();
      categorizer.findPassingPhotons();
      categorizer.classify();

      //find the 2 photons that decided the event and the di-EM Lorentz vector
      VINT passingPhotons = categorizer.getPassingPhotons();
      for (reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin(); 
	   iPhoton != pPhotons->end(); ++iPhoton) {
	if (((iPhoton - pPhotons->begin()) == passingPhotons[0]) || 
	    ((iPhoton - pPhotons->begin()) == passingPhotons[1])) {
	  photonPasses.push_back(true);
	  diEMP4+=iPhoton->p4();
	  ++passingPhotonCount;
	}
	else photonPasses.push_back(false);
      }
    }
    catch (STRING& badInput) { throw cms::Exception("EDMCategoryProducer") << badInput; }
  }

  //write the category information to the event
  try {
    putProductIntoEvent(pPhotons, categorizer.getPhotonType(), iEvent, "photonType");
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
    putProductIntoEvent(pPhotons, photonPasses, iEvent, "passingPhotons");
    putProductIntoEvent(categorizer.getCategory(), iEvent, "eventCategory");
    putProductIntoEvent(categorizer.getEvtPassDPhiMin(), iEvent, "passDPhiMin");
    putProductIntoEvent(categorizer.getEvtDiEMET(), iEvent, "evtDiEMET");
    if (passingPhotonCount > 2) {
      STRINGSTREAM err;
      err << "Error: " << passingPhotonCount << " passing photons.\n";
      throw cms::Exception("EDMCategoryProducer") << err.str();
    }
    if (passingPhotonCount == 2) evtInvMass = diEMP4.M();
    putProductIntoEvent(evtInvMass, iEvent, "evtInvMass");
    putProductIntoEvent(pPhotons, photonSeedTime, iEvent, "photonSeedTime");
    putProductIntoEvent(pPhotons, photonE2OverE9, iEvent, "photonE2OverE9");
    putProductIntoEvent(pPhotons, photonSeedIeta, iEvent, "photonSeedIeta");
  }
  catch (STRING& badInput) { throw cms::Exception("EDMCategoryProducer") << badInput; }
  catch (cms::Exception& badInput) { throw cms::Exception("EDMCategoryProducer") << badInput; }

  //print debug info
  if (categorizer.getCategory() == GG) debugPrint(categorizer, iEvent);
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
  if ((hit != pRecHits->end()) && hit->isTimeValid()) recHitTime = hit->time();
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

int EDMCategoryProducer::seedIeta(const reco::Photon* photon) const
{
  unsigned int theSeedIeta = 86;
  const DetId theSeedID = seedID(photon);
  try {
    theSeedIeta = EBDetId(theSeedID).ieta();
  }
  catch (cms::Exception& ex) { edm::LogInfo("InvalidEBDetId") << ex.what(); }
  return theSeedIeta;
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
  try { filler.insert(pPhotons, vals.begin(), vals.end()); }
  catch (cms::Exception& ex) {
    STRINGSTREAM err;
    err << "putProductIntoEvent/" << label << ": " << ex.what();
    throw cms::Exception("EDMCategoryProducer") << err.str();
  }
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
  try { filler.insert(pPhotons, vals.begin(), vals.end()); }
  catch (cms::Exception& ex) {
    STRINGSTREAM err;
    err << "putProductIntoEvent/" << label << ": " << ex.what();
    throw cms::Exception("EDMCategoryProducer") << err.str();
  }
  filler.fill();
  iEvent.put(out, label);
}

void 
EDMCategoryProducer::putProductIntoEvent(const edm::Handle<reco::PhotonCollection>& pPhotons, 
					 const VDOUBLE& vals, edm::Event& iEvent, 
					 const STRING& label) const
{
  std::auto_ptr<edm::ValueMap<double> > out(new edm::ValueMap<double>());
  edm::ValueMap<double>::Filler filler(*out);
  try { filler.insert(pPhotons, vals.begin(), vals.end()); }
  catch (cms::Exception& ex) {
    STRINGSTREAM err;
    err << "putProductIntoEvent/" << label << ": " << ex.what();
    throw cms::Exception("EDMCategoryProducer") << err.str();
  }
  filler.fill();
  iEvent.put(out, label);
}

void EDMCategoryProducer::putProductIntoEvent(const int vals, edm::Event& iEvent, 
					      const STRING& label) const
{
  std::auto_ptr<int> out(new int());
  *out = vals;
  iEvent.put(out, label);
}

void EDMCategoryProducer::putProductIntoEvent(const bool vals, edm::Event& iEvent, 
					      const STRING& label) const
{
  std::auto_ptr<bool> out(new bool());
  *out = vals;
  iEvent.put(out, label);
}

void EDMCategoryProducer::putProductIntoEvent(const double vals, edm::Event& iEvent, 
					      const STRING& label) const
{
  std::auto_ptr<double> out(new double());
  *out = vals;
  iEvent.put(out, label);
}

void EDMCategoryProducer::debugPrint(const Categorizer& categorizer, edm::Event& iEvent) const
{
  //don't print unless event processing is completely finished
  if (!categorizer.done()) {
    edm::LogInfo("EDMCategoryProducer") << "Error: finish event processing.\n";
    return;
  }

  //get all event information
  VDOUBLE photonET = categorizer.getPhotonET();
  VDOUBLE photonEta = categorizer.getPhotonEta();
  VDOUBLE photonECALIso = categorizer.getPhotonECALIso();
  VDOUBLE photonHCALIso = categorizer.getPhotonHCALIso();
  VDOUBLE photonHOverE = categorizer.getPhotonHOverE();
  VDOUBLE photonTrackIso = categorizer.getPhotonTrackIso();
  VDOUBLE photonSigmaIetaIeta = categorizer.getPhotonSigmaIetaIeta();
  VDOUBLE photonSeedTime = categorizer.getPhotonSeedTime();
  VDOUBLE photonE2OverE9 = categorizer.getPhotonE2OverE9();
  VDOUBLE photonPhi = categorizer.getPhotonPhi();
  VBOOL photonHasPixelSeed = categorizer.getPhotonHasPixelSeed();
  double photonETMin = categorizer.getPhotonETMin();
  double photonAbsEtaMax = categorizer.getPhotonAbsEtaMax();
  double photonECALIsoMaxPTMultiplier = categorizer.getPhotonECALIsoMaxPTMultiplier();
  double photonECALIsoMaxConstant = categorizer.getPhotonECALIsoMaxConstant();
  double photonHCALIsoMaxPTMultiplier = categorizer.getPhotonHCALIsoMaxPTMultiplier();
  double photonHCALIsoMaxConstant = categorizer.getPhotonHCALIsoMaxConstant();
  double photonHOverEMax = categorizer.getPhotonHOverEMax();
  double photonTrackIsoMaxPTMultiplier = categorizer.getPhotonTrackIsoMaxPTMultiplier();
  double photonTrackIsoMaxConstant = categorizer.getPhotonTrackIsoMaxConstant();
  double photonSigmaIetaIetaMax = categorizer.getPhotonSigmaIetaIetaMax();
  double photonAbsSeedTimeMax = categorizer.getPhotonAbsSeedTimeMax();
  double photonE2OverE9Max = categorizer.getPhotonE2OverE9Max();
  double photonDPhiMin = categorizer.getPhotonDPhiMin();
  VBOOL photonPassETMin = categorizer.getPhotonPassETMin();
  VBOOL photonPassAbsEtaMax = categorizer.getPhotonPassAbsEtaMax();
  VBOOL photonPassECALIsoMax = categorizer.getPhotonPassECALIsoMax();
  VBOOL photonPassHCALIsoMax = categorizer.getPhotonPassHCALIsoMax();
  VBOOL photonPassHOverEMax = categorizer.getPhotonPassHOverEMax();
  VBOOL photonPassTrackIsoMax = categorizer.getPhotonPassTrackIsoMax();
  VBOOL photonPassSigmaIetaIetaMax = categorizer.getPhotonPassSigmaIetaIetaMax();
  VBOOL photonPassAbsSeedTimeMax = categorizer.getPhotonPassAbsSeedTimeMax();
  VBOOL photonPassE2OverE9Max = categorizer.getPhotonPassE2OverE9Max();
  VBOOL photonPassPreselection = categorizer.getPhotonPassPreselection();
  bool evtPassDPhiMin = categorizer.getEvtPassDPhiMin();
  VINT passingPhotons = categorizer.getPassingPhotons();
  VINT photonType = categorizer.getPhotonType();
  int category = categorizer.getCategory();

  //loop over photons
  STRINGSTREAM debug;
  debug << "*****************\n";
  debug << "Run " << iEvent.run() << ", event " << iEvent.id().event() << ", lumi section ";
  debug << iEvent.getLuminosityBlock().luminosityBlock() << std::endl;
  for (VDOUBLE_IT iET = photonET.begin(); iET != photonET.end(); ++iET) {

    //print photon quantity, cut value, and pass flag
    debug << "%%%%%%%%%%%%%%%%%\n";
    const unsigned int i = iET - photonET.begin();
    debug << "Photon index: " << i << std::endl;
    debug << "-----------------\n";
    debug << "Photon ET: " << *iET << " GeV\n";
    debug << "Cut: photon ET > " << photonETMin << " GeV\n";
    debug << "Result: ";
    if (photonPassETMin[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon eta: " << photonEta[i] << std::endl;
    debug << "Cut: photon |eta| < " << photonAbsEtaMax << std::endl;
    debug << "Result: ";
    if (photonPassAbsEtaMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon ECAL isolation: " << photonECALIso[i] << " GeV\n";
    debug << "Cut: photon ECAL isolation < (" << photonECALIsoMaxPTMultiplier << "ET + ";
    debug << photonECALIsoMaxConstant << ") GeV = ";
    debug << (photonECALIsoMaxPTMultiplier*(*iET) + photonECALIsoMaxConstant) << " GeV\n";
    debug << "Result: ";
    if (photonPassECALIsoMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon HCAL isolation: " << photonHCALIso[i] << " GeV\n";
    debug << "Cut: photon HCAL isolation < (" << photonHCALIsoMaxPTMultiplier << "ET + ";
    debug << photonHCALIsoMaxConstant << ") GeV = ";
    debug << (photonHCALIsoMaxPTMultiplier*(*iET) + photonHCALIsoMaxConstant) << " GeV\n";
    debug << "Result: ";
    if (photonPassHCALIsoMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon H/E: " << photonHOverE[i] << std::endl;
    debug << "Cut: photon H/E < " << photonHOverEMax << std::endl;
    debug << "Result: ";
    if (photonPassHOverEMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon track isolation: " << photonTrackIso[i] << " GeV\n";
    debug << "Cut: photon track isolation < (" << photonTrackIsoMaxPTMultiplier << "ET + ";
    debug << photonTrackIsoMaxConstant << ") GeV = ";
    debug << (photonTrackIsoMaxPTMultiplier*(*iET) + photonTrackIsoMaxConstant) << " GeV\n";
    debug << "Result: ";
    if (photonPassTrackIsoMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon sigmaIetaIeta: " << photonSigmaIetaIeta[i] << std::endl;
    debug << "Cut: photon sigmaIetaIeta < " << photonSigmaIetaIetaMax << std::endl;
    debug << "Result: ";
    if (photonPassSigmaIetaIetaMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon seed time: " << photonSeedTime[i] << " ns\n";
    debug << "Cut: photon |seed time| < " << photonAbsSeedTimeMax << " ns\n";
    debug << "Result: ";
    if (photonPassAbsSeedTimeMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon E2/E9: " << photonE2OverE9[i] << std::endl;
    debug << "Cut: photon E2/E9 < " << photonE2OverE9Max << std::endl;
    debug << "Result: ";
    if (photonPassE2OverE9Max[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Result: ";
    if (photonPassPreselection[i]) debug << "photon passed preselection\n";
    else debug << "photon failed preselection\n";
    debug << "-----------------\n";
    debug << "Photon has pixel seed: ";
    if (photonHasPixelSeed[i]) debug << "yes\n";
    else debug << "no\n";
    debug << "-----------------\n";
    debug << "Photon type: ";
    switch (photonType[i]) {
    case FAIL:
      debug << "fail\n";
      break;
    case G:
      debug << "photon\n";
      break;
    case E:
      debug << "electron\n";
      break;
    case F:
      debug << "fake\n";
      break;
    default:
      debug << "unknown\n";
      break;
    }
    debug << "%%%%%%%%%%%%%%%%%\n";
  }

  //print event quantities
  debug << "Passing photon indices: " << passingPhotons[0] << ", " << passingPhotons[1];
  debug << std::endl;
  debug << "#################\n";
  debug << "Passing photon phi: ";
  if (passingPhotons[0] >= 0) debug << photonPhi[passingPhotons[0]] << ", ";
  else debug << "N/A" << ", ";
  if (passingPhotons[1] >= 0) debug << photonPhi[passingPhotons[1]] << std::endl;
  else debug << "N/A" << std::endl;
  debug << "Cut: dPhi(photon, photon) > " << photonDPhiMin << std::endl;
  debug << "Result: ";
  if (evtPassDPhiMin) debug << "pass\n";
  else debug << "fail\n";
  debug << "#################\n";
  debug << "Event category: ";
  switch (category) {
  case FAIL:
    debug << "fail\n";
    break;
  case GG:
    debug << "candidate\n";
    break;
  case EE:
    debug << "ee\n";
    break;
  case FF:
    debug << "ff\n";
    break;
  case EG:
    debug << "eg\n";
    break;
  default:
    debug << "unknown\n";
    break;
  }
  debug << "*****************\n";

  //log
  edm::LogInfo("EDMCategoryProducer") << debug.str();
}

//define this as a plug-in
DEFINE_FWK_MODULE(EDMCategoryProducer);
