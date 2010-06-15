#include "GMSBTools/EventSelection/interface/EventSelector.h"
#include <fstream>

//relevant namespaces
using namespace std;
using namespace reco;
using namespace edm;
using namespace cms;

//default constructor with fake-fake values as of 7-Jun-10
EventSelector::EventSelector() :
  sampleType_(FF),
  ECALIsoMaxPTMultiplierEB_(0.004),
  ECALIsoMaxConstantEB_(4.2),
  ECALIsoMaxPTMultiplierEE_(0.0),
  ECALIsoMaxConstantEE_(999999999.0),
  HCALIsoMaxPTMultiplierEB_(0.001),
  HCALIsoMaxConstantEB_(2.2),
  HCALIsoMaxPTMultiplierEE_(0.0),
  HCALIsoMaxConstantEE_(999999999.0),
  HOverEMaxPresel_(0.05),
  ETMin_(10.0),
  fiducialRegion_(EB),
  useHOverE_(false),
  HOverEMax_(0.0),
  useSigmaEtaEta_(false),
  sigmaEtaEtaMax_(0.0),
  useTrackIso_(true),
  trackIsoMaxPTMultiplier_(0.001),
  trackIsoMaxConstant_(2.0),
  trackPTMin_(15.0),
  eTrackRMin_(0.8),
  minDRPhotons_(0.8),
  numReqdCands_(2),
  run_(0),
  evt_(0),
  lumiSec_(0),
  debugFileName_("debug.txt"),
  debugFlag_(false) { initializeDebugging(); }

//constructor from a list of cut values
EventSelector::EventSelector(const unsigned int sampleType, const double ECALIsoMaxPTMultiplierEB, const double ECALIsoMaxConstantEB, 
			     const double ECALIsoMaxPTMultiplierEE, const double ECALIsoMaxConstantEE, const double HCALIsoMaxPTMultiplierEB, 
			     const double HCALIsoMaxConstantEB, const double HCALIsoMaxPTMultiplierEE, const double HCALIsoMaxConstantEE, 
			     const double HOverEMaxPresel, const double ETMin, const unsigned int fiducialRegion, const bool useHOverE, 
			     const double HOverEMax, const bool useSigmaEtaEta, const double sigmaEtaEtaMax, const bool useTrackIso, 
			     const double trackIsoMaxPTMultiplier, const double trackIsoMaxConstant, const double trackPTMin, const double eTrackRMin, 
			     const double minDRPhotons, const unsigned int numReqdCands, const unsigned int run, const unsigned int evt, 
			     const unsigned int lumiSec, const string debugFileName, const bool debugFlag) :
  sampleType_(sampleType),
  ECALIsoMaxPTMultiplierEB_(ECALIsoMaxPTMultiplierEB),
  ECALIsoMaxConstantEB_(ECALIsoMaxConstantEB),
  ECALIsoMaxPTMultiplierEE_(ECALIsoMaxPTMultiplierEE),
  ECALIsoMaxConstantEE_(ECALIsoMaxConstantEE),
  HCALIsoMaxPTMultiplierEB_(HCALIsoMaxPTMultiplierEB),
  HCALIsoMaxConstantEB_(HCALIsoMaxConstantEB),
  HCALIsoMaxPTMultiplierEE_(HCALIsoMaxPTMultiplierEE),
  HCALIsoMaxConstantEE_(HCALIsoMaxConstantEE),
  HOverEMaxPresel_(HOverEMaxPresel),
  ETMin_(ETMin),
  fiducialRegion_(fiducialRegion),
  useHOverE_(useHOverE),
  HOverEMax_(HOverEMax),
  useSigmaEtaEta_(useSigmaEtaEta),
  sigmaEtaEtaMax_(sigmaEtaEtaMax),
  useTrackIso_(useTrackIso),
  trackIsoMaxPTMultiplier_(trackIsoMaxPTMultiplier),
  trackIsoMaxConstant_(trackIsoMaxConstant),
  trackPTMin_(trackPTMin),
  eTrackRMin_(eTrackRMin),
  minDRPhotons_(minDRPhotons),
  numReqdCands_(numReqdCands),
  run_(run),
  evt_(evt),
  lumiSec_(lumiSec),
  debugFileName_(debugFileName),
  debugFlag_(debugFlag) { initializeDebugging(); }

//copy constructor
EventSelector::EventSelector(/*const */EventSelector& other) :
  sampleType_(other.getSampleType()),
  ECALIsoMaxPTMultiplierEB_(other.getECALIsoMaxPTMultiplierEB()),
  ECALIsoMaxConstantEB_(other.getECALIsoMaxConstantEB()),
  ECALIsoMaxPTMultiplierEE_(other.getECALIsoMaxPTMultiplierEE()),
  ECALIsoMaxConstantEE_(other.getECALIsoMaxConstantEE()),
  HCALIsoMaxPTMultiplierEB_(other.getHCALIsoMaxPTMultiplierEB()),
  HCALIsoMaxConstantEB_(other.getHCALIsoMaxConstantEB()),
  HCALIsoMaxPTMultiplierEE_(other.getHCALIsoMaxPTMultiplierEE()),
  HCALIsoMaxConstantEE_(other.getHCALIsoMaxConstantEE()),
  HOverEMaxPresel_(other.getHOverEMaxPresel()),
  ETMin_(other.getETMin()),
  fiducialRegion_(other.getFiducialRegion()),
  useHOverE_(other.getUseHOverE()),
  HOverEMax_(other.getHOverEMax()),
  useSigmaEtaEta_(other.getUseSigmaEtaEta()),
  sigmaEtaEtaMax_(other.getSigmaEtaEtaMax()),
  useTrackIso_(other.getUseTrackIso()),
  trackIsoMaxPTMultiplier_(other.getTrackIsoMaxPTMultiplier()),
  trackIsoMaxConstant_(other.getTrackIsoMaxConstant()),
  trackPTMin_(other.getTrackPTMin()),
  eTrackRMin_(other.getETrackRMin()),
  minDRPhotons_(other.getMinDRPhotons()),
  numReqdCands_(other.getNumReqdCands()),
  run_(other.getRun()),
  evt_(other.getEvt()),
  lumiSec_(other.getLumiSec()),
  debugFileName_(other.getDebugFileName()),
  debugFlag_(other.getDebugFlag())
{
  if (other.debugFileOpen()) other.closeDebugFile();
  initializeDebugging();
}

//destructor
EventSelector::~EventSelector() { closeDebugFile(); }

//assignment operator
EventSelector& EventSelector::operator=(/*const */EventSelector& other)
{
  if (this != &other) {
    sampleType_ = other.getSampleType();
    ECALIsoMaxPTMultiplierEB_ = other.getECALIsoMaxPTMultiplierEB();
    ECALIsoMaxConstantEB_ = other.getECALIsoMaxConstantEB();
    ECALIsoMaxPTMultiplierEE_ = other.getECALIsoMaxPTMultiplierEE();
    ECALIsoMaxConstantEE_ = other.getECALIsoMaxConstantEE();
    HCALIsoMaxPTMultiplierEB_ = other.getHCALIsoMaxPTMultiplierEB();
    HCALIsoMaxConstantEB_ = other.getHCALIsoMaxConstantEB();
    HCALIsoMaxPTMultiplierEE_ = other.getHCALIsoMaxPTMultiplierEE();
    HCALIsoMaxConstantEE_ = other.getHCALIsoMaxConstantEE();
    HOverEMaxPresel_ = other.getHOverEMaxPresel();
    ETMin_ = other.getETMin();
    fiducialRegion_ = other.getFiducialRegion();
    useHOverE_ = other.getUseHOverE();
    HOverEMax_ = other.getHOverEMax();
    useSigmaEtaEta_ = other.getUseSigmaEtaEta();
    sigmaEtaEtaMax_ = other.getSigmaEtaEtaMax();
    useTrackIso_ = other.getUseTrackIso();
    trackIsoMaxPTMultiplier_ = other.getTrackIsoMaxPTMultiplier();
    trackIsoMaxConstant_ = other.getTrackIsoMaxConstant();
    trackPTMin_ = other.getTrackPTMin();
    eTrackRMin_ = other.getETrackRMin();
    minDRPhotons_ = other.getMinDRPhotons();
    numReqdCands_ = other.getNumReqdCands();
    run_ = other.getRun();
    evt_ = other.getEvt();
    lumiSec_ = other.getLumiSec();
    debugFileName_ = other.getDebugFileName();
    debugFlag_ = other.getDebugFlag();
    if (other.debugFileOpen()) other.closeDebugFile();
    initializeDebugging();
  }
  return *this;
}

//getters
const unsigned int EventSelector::getSampleType() const { return sampleType_; }
const double EventSelector::getECALIsoMaxPTMultiplierEB() const { return ECALIsoMaxPTMultiplierEB_; }
const double EventSelector::getECALIsoMaxConstantEB() const { return ECALIsoMaxConstantEB_; }
const double EventSelector::getECALIsoMaxPTMultiplierEE() const { return ECALIsoMaxPTMultiplierEE_; }
const double EventSelector::getECALIsoMaxConstantEE() const { return ECALIsoMaxConstantEE_; }
const double EventSelector::getHCALIsoMaxPTMultiplierEB() const { return HCALIsoMaxPTMultiplierEB_; }
const double EventSelector::getHCALIsoMaxConstantEB() const { return HCALIsoMaxConstantEB_; }
const double EventSelector::getHCALIsoMaxPTMultiplierEE() const { return HCALIsoMaxPTMultiplierEE_; }
const double EventSelector::getHCALIsoMaxConstantEE() const { return HCALIsoMaxConstantEE_; }
const double EventSelector::getHOverEMaxPresel() const { return HOverEMaxPresel_; }
const double EventSelector::getETMin() const { return ETMin_; }
const unsigned int EventSelector::getFiducialRegion() const { return fiducialRegion_; }
const bool EventSelector::getUseHOverE() const { return useHOverE_; }
const double EventSelector::getHOverEMax() const { return HOverEMax_; }
const bool EventSelector::getUseSigmaEtaEta() const { return useSigmaEtaEta_; }
const double EventSelector::getSigmaEtaEtaMax() const { return sigmaEtaEtaMax_; }
const bool EventSelector::getUseTrackIso() const { return useTrackIso_; }
const double EventSelector::getTrackIsoMaxPTMultiplier() const { return trackIsoMaxPTMultiplier_; }
const double EventSelector::getTrackIsoMaxConstant() const { return trackIsoMaxConstant_; }
const double EventSelector::getTrackPTMin() const { return trackPTMin_; }
const double EventSelector::getETrackRMin() const { return eTrackRMin_; }
const double EventSelector::getMinDRPhotons() const { return minDRPhotons_; }
const unsigned int EventSelector::getNumReqdCands() const { return numReqdCands_; }
const unsigned int EventSelector::getRun() const { return run_; }
const unsigned int EventSelector::getEvt() const { return evt_; }
const unsigned int EventSelector::getLumiSec() const { return lumiSec_; }
const string EventSelector::getDebugFileName() const { return debugFileName_; }
const bool EventSelector::getDebugFlag() const { return debugFlag_; }

//setters
void EventSelector::setRun(const unsigned int run) { run_ = run; }
void EventSelector::setEvt(const unsigned int evt) { evt_ = evt; }
void EventSelector::setLumiSec(const unsigned int lumiSec) {lumiSec_ = lumiSec; }

//open the debug file and write some stuff to it
void EventSelector::initializeDebugging()
{
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
    debug_ << "Minimum dR(photon, track): " << eTrackRMin_ << endl;
    debug_ << "Minimum dR(photon, photon): " << minDRPhotons_ << endl;
    debug_ << "Number of required candidates: " << numReqdCands_ << endl;
  }
  if ((sampleType_ != GAMMAGAMMA) && (sampleType_ != EGAMMA) && (sampleType_ != EE) && (sampleType_ != FF) && (sampleType_ != ETRACK)) {
    if (debugFlag_) debug_ << "Invalid sample type chosen.  Defaulting to \"FF\".\n";
    sampleType_ = FF;
  }
  if ((fiducialRegion_ != EB) && (fiducialRegion_ != EEND) && (fiducialRegion_ != ECAL)) {
    if (debugFlag_) debug_ << "Invalid fiducial region chosen.  Defaulting to \"EB\"\n";
    fiducialRegion_ = EB;
  }
  if (debugFlag_) debug_ << "\n----------------------\n";
}

//print event information
void EventSelector::printEvtInfo() { if (debugFlag_) debug_ << "Run: " << run_ << ", event: " << evt_ << ", lumi section: " << lumiSec_ << endl; }

//print debugging information
void EventSelector::printDebug(const string info) { if (debugFlag_) debug_ << info << endl; }

//determine if the debug file is open
const bool EventSelector::debugFileOpen() const { return debug_.is_open(); }

//close the debug file
void EventSelector::closeDebugFile() { debug_.close(); }

//get the sample type as a string
const string EventSelector::sampleTypeString() const
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

//get the ECAL fiducial region in units this code understands
const unsigned int EventSelector::ECALFiducialRegion(const Photon* photon) const
{
  unsigned int reg = EEND;
  if ((photon->isEB()) && (fabs(fabs(photon->superCluster()->position().eta()) - 1.479) >= 0.1)) reg = EB;
  return reg;
}

//calculate the data quality flag for the event
const bool EventSelector::passesDataQualityCuts(const Handle<HBHERecHitCollection>& pHBHERecHits, const vector<Photon*>& passingCands, 
						const CaloGeometry* pGeometry, const Handle<TrackCollection>& pCosmicTracks)
{
  return !(passesHEBeamHaloTag(pHBHERecHits, passingCands, pGeometry) || passesMuonBeamHaloTag(pCosmicTracks, passingCands));
}

//calculate the HE tag for the event
const bool EventSelector::passesHEBeamHaloTag(const Handle<HBHERecHitCollection>& pHBHERecHits, const vector<Photon*>& passingCands, 
					      const CaloGeometry* pGeometry)
{
  bool HEHalo = false;
  vector<Photon*>::const_iterator iPassingCand = passingCands.begin();
  while ((iPassingCand != passingCands.end()) && (!HEHalo)) {
    HBHERecHitCollection::const_iterator iHBHERecHit = pHBHERecHits->begin();
    while ((iHBHERecHit != pHBHERecHits->end()) && (!HEHalo)) {
      const DetId HBHERecHitID = iHBHERecHit->detid();
      const float rho = pGeometry->getPosition(HBHERecHitID).perp();
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
  return HEHalo;
}

//calculate the muon tag for the event, using only the innermost hit position of the track as per Andrew's e-mail on 6-Jun-10
const bool EventSelector::passesMuonBeamHaloTag(const Handle<TrackCollection>& pCosmicTracks, const vector<Photon*>& passingCands)
{
  bool muonHalo = false;
  vector<Photon*>::const_iterator iPassingCand = passingCands.begin();
  while ((iPassingCand != passingCands.end()) && (!muonHalo)) {
    TrackCollection::const_iterator iCosmicTrack = pCosmicTracks->begin();
    while ((iCosmicTrack != pCosmicTracks->end()) && (!muonHalo)) {
      const double rhoInner = sqrt(iCosmicTrack->innerPosition().perp2());
      const double dPhiInner = dPhi((*iPassingCand)->phi(), iCosmicTrack->innerPosition().phi());
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
  return muonHalo;
}

//calculate the preselection flag for an individual reco::Photon
const bool EventSelector::passesPreselection(const Photon* photon, const unsigned int index)
{
  bool pass = false;
  const double pT = photon->pt();
  const double ECALIso = (double)photon->ecalRecHitSumEtConeDR04();
  double ECALIsoMax = ECALIsoMaxPTMultiplierEB_*pT + ECALIsoMaxConstantEB_;
  const double HCALIso = (double)photon->hcalTowerSumEtConeDR04();
  double HCALIsoMax = HCALIsoMaxPTMultiplierEB_*pT + HCALIsoMaxConstantEB_;
  const unsigned int fiducialRegion = ECALFiducialRegion(photon);
  if (fiducialRegion == EEND) {
    ECALIsoMax = ECALIsoMaxPTMultiplierEE_*pT + ECALIsoMaxConstantEE_;
    HCALIsoMax = HCALIsoMaxPTMultiplierEE_*pT + HCALIsoMaxConstantEE_;
  }
  const double HOverE = (double)photon->hadronicOverEm();
  const double ET = photon->et();
  if (debugFlag_) {
    debug_ << "Photon index: " << index << endl;
    debug_ << "Photon pT: " << pT << " GeV\n";
    debug_ << "Photon ECAL isolation: " << ECALIso << " GeV\n";
    debug_ << "Photon HCAL isolation: " << HCALIso << " GeV\n";
    debug_ << "Photon H/E: " << HOverE << " GeV\n";
    debug_ << "Photon ET: " << ET << " GeV\n";
    debug_ << "Photon fiducial region: " << fiducialRegion << " (";
    if (fiducialRegion == EB) debug_ << "EB";
    else debug_ << "EE";
    debug_ << ")\n";
    debug_ << "Photon eta: " << photon->eta() << endl;
    debug_ << "Photon phi: " << photon->phi() << endl;
  }
  if ((ECALIso < ECALIsoMax) && (HCALIso < HCALIsoMax) && (HOverE < HOverEMaxPresel_) && (ET > ETMin_) && 
      ((fiducialRegion_ == ECAL) || (fiducialRegion_ == fiducialRegion))) {
    if (debugFlag_) debug_ << "Photon object passes preselection.\n\n";
    pass = true;
  }
  else if (debugFlag_) debug_ << "Photon object fails preselection.\n\n";
  return pass;
}

//calculate the candidate ID flag for an individual reco::Photon
const bool EventSelector::passesCandidateID(const Photon* photon, bool& foundElectron, bool& foundPhoton)
{
  //pass flag
  bool pass = false;

  //decide if the candidate passes based on the photon ID cuts and flags
  const double pT = photon->pt();
  const double HOverE = (double)photon->hadronicOverEm();
  const double sigmaEtaEta = (double)photon->sigmaEtaEta();
  const double trackIso = (double)photon->trkSumPtHollowConeDR04();
  const double trackIsoMax = trackIsoMaxPTMultiplier_*pT + trackIsoMaxConstant_;
  if (debugFlag_) {
    debug_ << "Photon supercluster eta width: " << sigmaEtaEta << endl;
    debug_ << "Photon track isolation: " << trackIso << endl;
  }
  bool passedPhotonID = false;
  if ((useHOverE_ && useSigmaEtaEta_ && useTrackIso_) && ((sigmaEtaEta < sigmaEtaEtaMax_) && (HOverE < HOverEMax_) && (trackIso < trackIsoMax))) {
    passedPhotonID = true;
  }
  else if ((useHOverE_ && useSigmaEtaEta_ && (!useTrackIso_)) && ((sigmaEtaEta < sigmaEtaEtaMax_) && (HOverE < HOverEMax_))) passedPhotonID = true;
  else if ((useHOverE_ && (!useSigmaEtaEta_) && useTrackIso_) && ((HOverE < HOverEMax_) && (trackIso < trackIsoMax))) passedPhotonID = true;
  else if (((!useHOverE_) && useSigmaEtaEta_ && useTrackIso_) && ((sigmaEtaEta < sigmaEtaEtaMax_) && (trackIso < trackIsoMax))) passedPhotonID = true;
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

    //egamma, ee, and etrack samples require the candidate to have at least one pixel seed
    if (photon->hasPixelSeed()) {

      /*ee and etrack samples require only electrons
	for the etrack sample, we also need to save the candidate*/
      if ((sampleType_ == EE) || (sampleType_ == ETRACK)) {
	pass = true;
	if (debugFlag_) debug_ << "Found an EE/ETRACK candidate.\n\n";
      }
      
      //egamma sample counts this electron only if another one wasn't previously counted
      if ((sampleType_ == EGAMMA) && (!foundElectron)) {
	pass = true;
	foundElectron = true;
	if (debugFlag_) debug_ << "Found an EGAMMA candidate electron.\n\n";
      }
    }

    //gammagamma and egamma samples require no pixel seed
    else {

      //gammagamma sample requires only photons
      if (sampleType_ == GAMMAGAMMA) {
	pass = true;
	if (debugFlag_) debug_ << "Found a GAMMAGAMMA candidate.\n\n";
      }

      //egamma sample counts this photons only if another one wasn't previously counted
      if ((sampleType_ == EGAMMA) && (!foundPhoton)) {
	pass = true;
	foundPhoton = true;
	if (debugFlag_) debug_ << "Found an EGAMMA candidate photon.\n\n";
      }
    }
  }

  /*ff sample requires:
    sigmaEtaEta >= sigmaEtaEtaMax_ or
    H/E >= HOverEMax_ or
    track isolation >= trackIsoMaxPTMultiplier_*pT + trackIsoMaxConstant_ GeV
    
    NB. No requirement on number of pixel seeds in the fake definition!
  */
  else if ((sampleType_ == FF) && (!passedPhotonID)) {
    pass = true;
    if (debugFlag_) debug_ << "Found an FF candidate.\n\n";
  }

  //exit
  return pass;
}

//determine whether the event has the right number of candidate reco::Photon objects or not
const bool EventSelector::foundPhotonCandidates(const Handle<PhotonCollection>& pPhotons, vector<Photon*>& passingCands)
{
  //pass flag
  bool pass = false;

  //loop over the photons looking for objects passing specified criteria
  PhotonCollection::const_iterator iPhoton = pPhotons->begin();
  unsigned int numCands = 0;
  unsigned int numPhotonsProcessed = 0;
  bool foundPhoton = false;
  bool foundElectron = false;
  while ((iPhoton != pPhotons->end()) && (numCands < numReqdCands_)) {
    const unsigned int index = iPhoton - pPhotons->begin();
    const Photon* photon = const_cast<const Photon*>(&*iPhoton);

    //does the candidate pass the preselection and ID criteria?
    if ((passesPreselection(photon, index)) && (passesCandidateID(photon, foundElectron, foundPhoton))) {
      passingCands.push_back(const_cast<Photon*>(&*iPhoton));
      ++numCands;
      if (debugFlag_) debug_ << "Number of found candidates: " << numCands << endl << endl;
    }

    //advance to the next photon in the collection
    ++iPhoton;
    ++numPhotonsProcessed;
  }
  if (debugFlag_) {
    debug_ << "Number of photons processed: " << numPhotonsProcessed << "/" << pPhotons->size() << endl;
    debug_ << "Total number of found candidates: " << numCands << endl << endl;
  }

  //exit
  if (numCands == numReqdCands_) {
    pass = true;
    if (debugFlag_) debug_ << "Required number of candidates found.\n\n";
  }
  else if (debugFlag_) debug_ << "Required number of candidates not found.\n\n";
  return pass;
}

//determine if a passing track was found for the ETRACK sample
const bool EventSelector::foundTrack(const Handle<TrackCollection>& pTracks, const vector<Photon*>& passingCands)
{
  //check that passingCands only has one element
  const unsigned int numElectrons = passingCands.size();
  if (numElectrons != 1) throw numElectrons;
  Photon* e = passingCands[0];

  //pass flag
  bool trackFound = false;

  /*etrack sample requires:
    track pT > trackPTMin_
    dR(track, electron) > eTrackRMin_
  */
  TrackCollection::const_iterator iTrack = pTracks->begin();
  while ((iTrack != pTracks->end()) && (!trackFound)) {
    double dRETrack = dR(e->eta(), iTrack->eta(), e->phi(), iTrack->phi());
    double pT = iTrack->pt();
    if (debugFlag_) {
      debug_ << "Track pT: " << pT << " GeV\n";
      debug_ << "dR(track, photon): " << dRETrack << endl;
    }
    if ((pT > trackPTMin_) && (dRETrack > eTrackRMin_)) {
      trackFound = true;
      if (debugFlag_) debug_ << "Found a candidate track.\n\n";
    }
    else {
      ++iTrack;
      if (debugFlag_) debug_ << "This track doesn't meet the criteria; advancing to the next track.\n\n";
    }
  }

  //exit
  return trackFound;
}

//determine if the candidates are separated by the dR cut value
const bool EventSelector::passDRCut(const vector<Photon*>& passingCands)
{
  //check that passingCands has two or greater elements
  const unsigned int numCands = passingCands.size();
  if (numCands < 2) throw numCands;

  bool pass = false;
  if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA) || (sampleType_ == EE) || (sampleType_ == FF)) {
    Photon* photon1 = *(passingCands.begin());
    Photon* photon2 = *(passingCands.begin() + 1);
    const double dRPhoton1Photon2 = dR(photon1->eta(), photon2->eta(), photon1->phi(), photon2->phi());
    if (debugFlag_) debug_ << "dR between candidates: " << dRPhoton1Photon2 << endl;
    if (dRPhoton1Photon2 > minDRPhotons_) {
      pass = true;
      if (debugFlag_) debug_ << "reco::Photon objects pass the dR cut.\n";
    }
    else if (debugFlag_) debug_ << "reco::Photon objects fail the dR cut.\n";
  }
  else pass = true; //no dR requirement on ETRACK sample
  return pass;
}

//determine if all required candidates were found
const bool EventSelector::foundAllCandidates(const Handle<PhotonCollection>& pPhotons, const Handle<TrackCollection>& pTracks, vector<Photon*>& passingCands)
{
  bool pass = true;
  if (foundPhotonCandidates(pPhotons, passingCands)) {
    if (sampleType_ == ETRACK) {
      try { if (foundTrack(pTracks, passingCands)) pass = true; }
      catch (const unsigned int badSize) {
	if (debugFlag_) {
	  debug_ << "Number of electron candidates for the ETRACK sample is " << badSize << ", not one.  This event will fail.  Check your code.\n";
	}
      }
    }
    else pass = true;
  }
  return pass;
}
