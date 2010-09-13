#include "GMSBTools/EventSelection/interface/EventSelector.h"

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
  useTimingCut_(true),
  maxSeedTime_(3.0),
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
			     const double minDRPhotons, const bool useTimingCut, const double maxSeedTime, const unsigned int numReqdCands, 
			     const unsigned int run, const unsigned int evt, const unsigned int lumiSec, const string debugFileName, const bool debugFlag) :
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
  useTimingCut_(useTimingCut),
  maxSeedTime_(maxSeedTime),
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
  useTimingCut_(other.getUseTimingCut()),
  maxSeedTime_(other.getMaxSeedTime()),
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
    useTimingCut_ = other.getUseTimingCut();
    maxSeedTime_ = other.getMaxSeedTime();
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
const bool EventSelector::getUseTimingCut() const { return useTimingCut_; }
const double EventSelector::getMaxSeedTime() const { return maxSeedTime_; }
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
    if (useTimingCut_) debug_ << "Maximum seed time: " << maxSeedTime_ << " ns\n";
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
const bool EventSelector::passesDataQualityCuts(const Handle<HBHERecHitCollection>& pHBHERecHits, const map<unsigned int, Photon*>& passingCands, 
						const CaloGeometry* pGeometry, const Handle<TrackCollection>& pCosmicTracks)
{
  return !(passesHEBeamHaloTag(pHBHERecHits, passingCands, pGeometry) || passesMuonBeamHaloTag(pCosmicTracks, passingCands));
}

//calculate the HE tag for the event
const bool EventSelector::passesHEBeamHaloTag(const Handle<HBHERecHitCollection>& pHBHERecHits, const map<unsigned int, Photon*>& passingCands, 
					      const CaloGeometry* pGeometry)
{
  bool HEHalo = false;
  map<unsigned int, Photon*>::const_iterator iPassingCand = passingCands.begin();
  while ((iPassingCand != passingCands.end()) && (!HEHalo)) {
    if (photonIsHEHalo(pHBHERecHits, const_cast<const Photon*>((*iPassingCand).second), pGeometry)) HEHalo = true;
    ++iPassingCand;
  }
  return HEHalo;
}

//determine if a particular photon is due to HE-tagged halo bremsstrahlung
const bool EventSelector::photonIsHEHalo(const Handle<HBHERecHitCollection>& pHBHERecHits, const Photon* photon, const CaloGeometry* pGeometry)
{
  bool isHEHalo = false;
  HBHERecHitCollection::const_iterator iHBHERecHit = pHBHERecHits->begin();
  while ((iHBHERecHit != pHBHERecHits->end()) && (!isHEHalo)) {
    const DetId HBHERecHitID = iHBHERecHit->detid();
    const float rho = pGeometry->getPosition(HBHERecHitID).perp();
    const double phi = (double)pGeometry->getPosition(HBHERecHitID).phi();
    const double E = iHBHERecHit->energy();
    const double dPhiPhotonHBHERecHit = dPhi(photon->phi(), phi);
    if ((iHBHERecHit->id().subdet() == HcalEndcap) && (E > 1/*GeV*/) && (rho < 130/*cm*/) && (rho > 115/*cm*/) && (dPhiPhotonHBHERecHit <= 0.2)) {
      if (debugFlag_) {
	debug_ << "Found HE halo candidate:\n";
	debug_ << "DetID: " << HBHERecHitID.rawId() << endl;
	debug_ << "Rho: " << rho << " cm\n";
	debug_ << "Energy: " << E << " GeV\n";
	debug_ << "dPhi: " << dPhiPhotonHBHERecHit << endl << endl;
      }
      isHEHalo = true;
    }
    ++iHBHERecHit;
  }
  return isHEHalo;
}

//calculate the muon tag for the event, using only the innermost hit position of the track as per Andrew's e-mail on 6-Jun-10
const bool EventSelector::passesMuonBeamHaloTag(const Handle<TrackCollection>& pCosmicTracks, const map<unsigned int, Photon*>& passingCands)
{
  bool muonHalo = false;
  map<unsigned int, Photon*>::const_iterator iPassingCand = passingCands.begin();
  while ((iPassingCand != passingCands.end()) && (!muonHalo)) {
    if (photonIsMuonHalo(pCosmicTracks, const_cast<const Photon*>((*iPassingCand).second))) muonHalo = true;
    ++iPassingCand;
  }
  return muonHalo;
}

//determine if a particular photon is due to muon-tagged halo bremsstrahlung
const bool EventSelector::photonIsMuonHalo(const Handle<TrackCollection>& pCosmicTracks, const Photon* photon)
{
  bool isMuonHalo = false;
  TrackCollection::const_iterator iCosmicTrack = pCosmicTracks->begin();
  while ((iCosmicTrack != pCosmicTracks->end()) && (!isMuonHalo)) {
    const double rhoInner = sqrt(iCosmicTrack->innerPosition().perp2());
    const double dPhiInner = dPhi(photon->phi(), iCosmicTrack->innerPosition().phi());
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
      isMuonHalo = true;
    }
    ++iCosmicTrack;
  }
  return isMuonHalo;
}

//calculate the preselection flag for an individual reco::Photon
const bool EventSelector::passesPreselection(const Photon* photon, const vector<EcalRecHit*>& ECALRecHits, const unsigned int index)
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
  DetId seedCrystalID = photon->superCluster()->seed()->seed();
  vector<EcalRecHit*>::const_iterator iECALRecHit = ECALRecHits.begin();
  bool foundSeed = false;
  float seedTime = 0.0;
  if (((fiducialRegion_ == ECAL) || (fiducialRegion_ == fiducialRegion)) && (useTimingCut_)) {
    while ((iECALRecHit != ECALRecHits.end()) && (!foundSeed)) {
      if ((*iECALRecHit)->detid() == seedCrystalID) {
	foundSeed = true;
	seedTime = (*iECALRecHit)->time();
      }
      ++iECALRecHit;
    }
    if ((!foundSeed) && (debugFlag_)) {
      debug_ << "Error: seed crystal was not found in supplied ECAL RecHit collection.  Timing cut will not be applied for this event.\n";
    }
  }
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
    if (foundSeed && useTimingCut_) debug_ << "Photon seed crystal time: " << seedTime << " ns\n";
  }
  if ((ECALIso < ECALIsoMax) && (HCALIso < HCALIsoMax) && (HOverE < HOverEMaxPresel_) && (ET > ETMin_) && 
      ((fiducialRegion_ == ECAL) || (fiducialRegion_ == fiducialRegion)) && ((seedTime < maxSeedTime_) || (!useTimingCut_))) {
    if (debugFlag_) debug_ << "Photon object passes preselection.\n\n";
    pass = true;
  }
  else if (debugFlag_) debug_ << "Photon object fails preselection.\n\n";
  return pass;
}

//calculate the candidate ID flag for an individual reco::Photon
const bool EventSelector::passesCandidateID(const Photon* photon, unsigned int& type)
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
      if ((sampleType_ == EE) || (sampleType_ == ETRACK) || (sampleType_ == EGAMMA)) {
	pass = true;
	if (debugFlag_) debug_ << "Found an EE/ETRACK/EGAMMA candidate.\n\n";
	type = HAS_PIXEL_SEED;
      }
    }

    //gammagamma and egamma samples require no pixel seed
    else {

      //gammagamma sample requires only photons
      if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA)) {
	pass = true;
	if (debugFlag_) debug_ << "Found a GAMMAGAMMA/EGAMMA candidate.\n\n";
	type = LACKS_PIXEL_SEED;
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
    type = LACKS_PIXEL_SEED; //actually, it can have a pixel seed, but isn't _required_ to have a pixel seed
  }

  //exit
  return pass;
}

//determine whether the event has the right number of candidate reco::Photon objects or not
const bool EventSelector::foundPhotonCandidates(const Handle<PhotonCollection>& pPhotons, const vector<EcalRecHit*>& ECALRecHits, 
						map<unsigned int, Photon*>& passingWithoutPixelSeed, map<unsigned int, Photon*>& passingWithPixelSeed)
{
  //pass flag
  bool pass = false;

  //loop over the photons looking for objects passing specified criteria
  PhotonCollection::const_iterator iPhoton = pPhotons->begin();
  unsigned int numCands = 0;
  unsigned int numPhotonsProcessed = 0;
  while (iPhoton != pPhotons->end()) {
    const unsigned int index = iPhoton - pPhotons->begin();
    const Photon* constPhoton = const_cast<const Photon*>(&*iPhoton);

    //does the candidate pass the preselection and ID criteria?
    unsigned int type = 0;
    if ((passesPreselection(constPhoton, ECALRecHits, index)) && (passesCandidateID(constPhoton, type))) {
      Photon* photon = const_cast<Photon*>(constPhoton);
      if (type == LACKS_PIXEL_SEED) passingWithoutPixelSeed[index] = photon;
      if (type == HAS_PIXEL_SEED) passingWithPixelSeed[index] = photon;
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
  if (numCands /*==*/>= numReqdCands_) {
    if (passDRCut(passingWithoutPixelSeed, passingWithPixelSeed)) pass = true;
  }
  else if (debugFlag_) debug_ << "Required number of candidates not found.\n\n";
  return pass;
}

//determine if a passing track was found for the ETRACK sample
const bool EventSelector::foundTrack(const Handle<TrackCollection>& pTracks, const map<unsigned int, Photon*>& passingWithPixelSeed)
{
  //pass flag
  bool trackFound = false;

  /*etrack sample requires:
    track pT > trackPTMin_
    dR(track, electron) > eTrackRMin_
  */
  map<unsigned int, Photon*>::const_iterator iElectron = passingWithPixelSeed.begin();
  while ((iElectron != passingWithPixelSeed.end()) && (!trackFound)) {
    if (electronHasPassingTrack(pTracks, const_cast<const Photon*>((*iElectron).second), (*iElectron).first)) trackFound = true;
    if ((!trackFound) && (debugFlag_)) {
      debug_ << "This electron doesn't meet the criteria when paired with any tracks; advancing to the next electron.\n\n";
    }
    ++iElectron;
  }

  //exit
  return trackFound;
}

//determine if there is an energetic track separated by a minimum dR value from a particular electron
const bool EventSelector::electronHasPassingTrack(const Handle<TrackCollection>& pTracks, const Photon* electron, const unsigned int index)
{
  bool trackFound = false;
  TrackCollection::const_iterator iTrack = pTracks->begin();
  while ((iTrack != pTracks->end()) && (!trackFound)) {
    double dRETrack = dR(electron->eta(), iTrack->eta(), electron->phi(), iTrack->phi());
    double pT = iTrack->pt();
    if (debugFlag_) {
      debug_ << "Track pT: " << pT << " GeV\n";
      debug_ << "dR(track " << (iTrack - pTracks->begin()) << ", photon " << index << "): " << dRETrack << endl;
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
  return trackFound;
}

//determine if the candidates are separated by the dR cut value
const bool EventSelector::passDRCut(const map<unsigned int, Photon*>& passingWithoutPixelSeed, const map<unsigned int, Photon*>& passingWithPixelSeed)
{
  //pass flag
  unsigned int numNonoverlapping = 0;

  //copy the correct input map for the calculation of dR between two objects in a single vector
  map<unsigned int, Photon*> passingCandsMap;
  if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == FF)) {
    passingCandsMap = passingWithoutPixelSeed;
  }
  else if (sampleType_ == EE) passingCandsMap = passingWithPixelSeed;

  //calculate dR between objects in one vector
  if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == EE) || (sampleType_ == FF)) {
    map<unsigned int, Photon*>::const_iterator iCand = passingCandsMap.begin();
    map<unsigned int, Photon*> passingCandsMapReduced = passingCandsMap;
    passingCandsMapReduced.erase(passingCandsMapReduced.begin());
    while ((iCand != passingCandsMap.end()) && (numNonoverlapping < 2)) {
      const unsigned int index = (*iCand).first;
      if (photonIsNonoverlapping(const_cast<const Photon*>((*iCand).second), index, passingCandsMapReduced)) {
	++numNonoverlapping;
	if (debugFlag_) debug_ << "Photon " << index << " is nonoverlapping.\n";
      }
      ++iCand;
      if (passingCandsMapReduced.size() > 0) passingCandsMapReduced.erase(passingCandsMapReduced.begin());
    }
    if ((numNonoverlapping < 2) && debugFlag_) debug_ << passingCandsMap.size() << " candidates were found, but none pass the dR cut.\n\n";
  }

  //calculate dR between objects in two vectors (i.e. separate the photons and electrons so we don't calculate dR between two photons or two electrons)
  else if (sampleType_ == EGAMMA) {
    map<unsigned int, Photon*>::const_iterator iPhoton = passingWithoutPixelSeed.begin();
    map<unsigned int, Photon*>::const_iterator iElectron = passingWithPixelSeed.begin();
    while ((iPhoton != passingWithoutPixelSeed.end()) && (numNonoverlapping < 2)) {
      const unsigned int index = (*iPhoton).first;
      if (photonIsNonoverlapping(const_cast<const Photon*>((*iPhoton).second), index, passingWithPixelSeed)) {
	++numNonoverlapping;
	if (debugFlag_) debug_ << "Photon " << index << " is nonoverlapping.\n";
      }
      ++iPhoton;
    }
    if ((numNonoverlapping < 2) && debugFlag_) {
      const unsigned int numPhotons = passingWithoutPixelSeed.size();
      const unsigned int numElectrons = passingWithPixelSeed.size();
      debug_ << numPhotons << " photon and " << numElectrons << " electron candidate";
      if (numElectrons != 1) debug_ << "s";
      debug_ << " were found, but none pass the dR cut.\n\n";
    }
  }

  //no dR requirement on ETRACK sample
  else numNonoverlapping = 2;

  //exit
  if (numNonoverlapping >= 2) return true;
  else return false;
}

//determine if there are any other candidates within a minimum dR value of a particular candidate
const bool EventSelector::photonIsNonoverlapping(const Photon* photon, const unsigned int index, const map<unsigned int, Photon*>& otherPassingCands)
{
  bool isOverlapping = false;
  map<unsigned int, Photon*>::const_iterator iOtherPassingCand = otherPassingCands.begin();
  while ((iOtherPassingCand != otherPassingCands.end()) && (!isOverlapping)) {
    Photon* cand = (*iOtherPassingCand).second;
    const double dRPhotons = dR(photon->eta(), cand->eta(), photon->phi(), cand->phi());
    if (debugFlag_) {
      debug_ << "Photons " << index << " and " << (*iOtherPassingCand).first << " are separated by dR = " << dRPhotons << ".\n";
    }
    if (dRPhotons <= minDRPhotons_) isOverlapping = true;
    ++iOtherPassingCand;
  }
  return (!isOverlapping);
}

//determine if all required candidates were found
const bool EventSelector::foundAllCandidates(const Handle<PhotonCollection>& pPhotons, const vector<EcalRecHit*>& ECALRecHits, 
					     const Handle<TrackCollection>& pTracks, map<unsigned int, Photon*>& passingWithoutPixelSeed, 
					     map<unsigned int, Photon*>& passingWithPixelSeed)
{
  bool pass = true;
  if (foundPhotonCandidates(pPhotons, ECALRecHits, passingWithoutPixelSeed, passingWithPixelSeed)) {
    if (sampleType_ == ETRACK) {
      try { if (foundTrack(pTracks, passingWithPixelSeed)) pass = true; }
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
