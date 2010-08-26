// -*- C++ -*-
//
// Package:    CutAnalyzer
// Class:      CutAnalyzer
// 
/**\class CutAnalyzer CutAnalyzer.cc GMSBTools/CutAnalyzer/src/CutAnalyzer.cc

 Description: plot distributions of analysis cuts

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Tue May 18 13:42:51 CEST 2010
// $Id$
//
//

//all the relevant header files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "GMSBTools/EventSelection/interface/EventSelector.h"
#include "DataFormats/PatCandidates/interface/Flags.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <numeric>
#include "DataFormats/JetReco/interface/CaloJet.h"

//ROOT include files
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TTree.h"

#define MAX_NUM_JETS 50

//
// constants, enums and typedefs
//
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzV;

//
// class declaration
//

class CutAnalyzer : public edm::EDAnalyzer {
   public:
      explicit CutAnalyzer(const edm::ParameterSet&);
      ~CutAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  //for alphaT calculation
  /*std::vector<double> deltaSumPt_permutations(const std::vector<LorentzV>& p4s);
    double alphaT(const std::vector<LorentzV>& p4s);*/
  double alphaT(const vector<LorentzV>&, vector<const Photon*>&, const unsigned int, const unsigned int) const;
  void minDHTRecursive(double&, vector<double>&, const vector<double>&, const unsigned int) const;

  //retrieve collection from the event
  template <typename T>
  const bool getCollection_(T& pCollection, const InputTag tag, const Event& iEvent)
  {
    bool collectionFound = false;
    try { collectionFound = (iEvent.getByLabel(tag, pCollection)) && (pCollection->size() > 0); }
    catch (cms::Exception& ex) {}
    if (!collectionFound) cerr << "No collection of type " << tag << " found in event " << iEvent.id().event() << ".\n";
    return collectionFound;
  }

      // ----------member data ---------------------------

  //output
  string outFileName_;
  TFile* out_;
  string debugFileName_;
  ofstream debug_;
  bool debugFlag_;

  //cut values
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
  bool useTimingCut_;
  double maxSeedTime_;
  double minJetPT_;
  double maxJetAbsEta_;
  double minJetAbsEta_;
  bool rejectHalo_;
  unsigned int numReqdCands_;
  EventSelector evtProperties_;

  //input
  InputTag photonTag_;
  InputTag trackTag_;
  InputTag HBHERecHitTag_;
  InputTag cosmicTrackTag_;
  InputTag EBRecHitTag_;
  InputTag EERecHitTag_;
  InputTag caloMETTag_;
  InputTag tcMETTag_;
  ESHandle<CaloGeometry> caloGeometryHandle_;
  InputTag PATAK5CaloJetTag_;
  InputTag RECOAK5CaloJetTag_;
  InputTag AK5JetIDTag_;

  //histograms
  TH1F* ECALIso_;
  TH1F* HCALIso_;
  TH1D* ET_;
  TH1F* HOverE_;
  TH1F* fiducialRegionDist_;
  TH1F* etaWidth_;
  TH1F* trackIso_;
  TH1D* trackPT_;
  TH1D* dRTrackPhoton_;
  TH1F* numCandsPerEvt_;
  TH1F* numTracksPerEvt_;
  TH1D* eta_;
  TH1D* phi_;
  TH1D* diEMPT_;
  TH1D* openingAngle_;
  TH1D* dPhiCands_;
  TH1D* dPhiMETCand1_;
  TH1D* dPhiMETCand2_;
  TH1D* dPhiMETCand_;
  TH1D* openingAngleMETDiEMPT_;
  TH1D* dPhiMETDiEMPT_;
  TH1F* numJets_;
  TH1F* numJets30_;
  TH1D* jetPTJetIDApplied_;
  TH1D* jetPTJetIDNotApplied_;
  TH1D* jetEtaJetIDApplied_;
  TH1D* jetEtaPlus20JetIDApplied_;
  TH1D* jetEtaPlus40JetIDApplied_;
  TH1D* jetEtaJetIDNotApplied_;
  TH1D* jetEtaPlus20JetIDNotApplied_;
  TH1D* jetEtaPlus40JetIDNotApplied_;
  TH1D* alphaT2Photons0Jets_;
  TH1D* alphaT2Photons1Jets_;
  TH1D* alphaT2Photons2Jets_;
  TH1D* alphaT2Photons3Jets_;
  TH1D* alphaT3Photons0Jets_;
  TH1D* alphaT3Photons1Jets_;
  TH1D* alphaT3Photons2Jets_;
  TH1D* alphaTPhotons_;
  TH2D* alphaTVsNJetsSelected_;
  TH1D* evtAlphaT_;
  TH1D* evtAlphaT50_;
  TH1F* hadEnergyInHBJetIDApplied_;
  TH1F* hadEnergyInHEJetIDApplied_;
  TH1F* hadEnergyInHBJetIDNotApplied_;
  TH1F* hadEnergyInHEJetIDNotApplied_;
  TH2F* hadEnergyInHBVsEtaJetIDApplied_;
  TH2F* hadEnergyInHEVsEtaJetIDApplied_;
  TH2F* hadEnergyInHBVsEtaJetIDNotApplied_;
  TH2F* hadEnergyInHEVsEtaJetIDNotApplied_;
  TH2D* jetCorrVsPTJetIDApplied_;
  TH2D* jetCorrVsPTJetIDNotApplied_;
  TH1F* jetCorrJetIDApplied_;
  TH1F* jetCorrJetIDNotApplied_;
  TH1D* jetPTRawJetIDApplied_;
  TH1D* jetPTRawJetIDNotApplied_;
  TH1D* jetPTRECOJetIDApplied_;
  TH1D* jetPTRECOJetIDNotApplied_;

  //graphs
  TGraphErrors* alphaTFraction_;
  TCanvas* alphaTFractionCanvas_;

  //tree
  UInt_t run_;
  UInt_t evt_;
  UInt_t lumiSection_;
  Double_t diEMPTVal_;
  Double_t tcMET_;
  TTree* evtTree_;
};

//
// static data member definitions
//

//
// constructors and destructor
//
CutAnalyzer::CutAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getUntrackedParameter<string>("outFileName", "out.root")),
  debugFileName_(iConfig.getUntrackedParameter<string>("debugFileName", "debug.txt")),
  debugFlag_(iConfig.getUntrackedParameter<bool>("debugFlag", false)),
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
  trackPTMin_(iConfig.getUntrackedParameter<double>("trackPTMin", -1.0)),
  eTrackRMin_(iConfig.getUntrackedParameter<double>("eTrackRMin", -1.0)),
  minDRPhotons_(iConfig.getParameter<double>("minDRPhotons")),
  useTimingCut_(iConfig.getUntrackedParameter<bool>("useTimingCut", true)),
  maxSeedTime_(iConfig.getParameter<double>("maxSeedTime")),
  minJetPT_(iConfig.getParameter<double>("minJetPT")),
  maxJetAbsEta_(iConfig.getParameter<double>("maxJetAbsEta")),
  minJetAbsEta_(iConfig.getParameter<double>("minJetAbsEta")),
  rejectHalo_(iConfig.getUntrackedParameter<bool>("rejectHalo", true)),
  photonTag_(iConfig.getParameter<InputTag>("photonTag")),
  trackTag_(iConfig.getParameter<InputTag>("trackTag")),
  HBHERecHitTag_(iConfig.getParameter<InputTag>("HBHERecHitTag")),
  cosmicTrackTag_(iConfig.getParameter<InputTag>("cosmicTrackTag")),
  EBRecHitTag_(iConfig.getParameter<InputTag>("EBRecHitTag")),
  EERecHitTag_(iConfig.getParameter<InputTag>("EERecHitTag")),
  caloMETTag_(iConfig.getParameter<InputTag>("caloMETTag")),
  tcMETTag_(iConfig.getParameter<InputTag>("tcMETTag")),
  PATAK5CaloJetTag_(iConfig.getParameter<InputTag>("PATAK5CaloJetTag")),
  RECOAK5CaloJetTag_(iConfig.getParameter<InputTag>("RECOAK5CaloJetTag")),
  AK5JetIDTag_(iConfig.getParameter<InputTag>("AK5JetIDTag"))

{
   //now do what ever initialization is needed
  numReqdCands_ = 2;
  if (sampleType_ == ETRACK) numReqdCands_ = 1;
  EventSelector evtProperties(sampleType_, ECALIsoMaxPTMultiplierEB_, ECALIsoMaxConstantEB_, ECALIsoMaxPTMultiplierEE_, ECALIsoMaxConstantEE_, 
			      HCALIsoMaxPTMultiplierEB_, HCALIsoMaxConstantEB_, HCALIsoMaxPTMultiplierEE_, HCALIsoMaxConstantEE_, HOverEMaxPresel_, 
			      ETMin_, fiducialRegion_, useHOverE_, HOverEMax_, useSigmaEtaEta_, sigmaEtaEtaMax_, useTrackIso_, trackIsoMaxPTMultiplier_, 
			      trackIsoMaxConstant_, trackPTMin_, eTrackRMin_, minDRPhotons_, useTimingCut_, maxSeedTime_, numReqdCands_, 0, 0, 0, 
			      debugFileName_, debugFlag_);
  evtProperties_ = evtProperties;
}


CutAnalyzer::~CutAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CutAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get the run, event, and lumi section numbers
  const unsigned int runNum = iEvent.run();
  const unsigned int evtNum = iEvent.id().event();
  const unsigned int lumiNum = iEvent.getLuminosityBlock().luminosityBlock();
  run_ = (UInt_t)runNum;
  evt_ = (UInt_t)evtNum;
  lumiSection_ = (UInt_t)lumiNum;
  evtProperties_.setRun(runNum);
  evtProperties_.setEvt(evtNum);
  evtProperties_.setLumiSec(lumiNum);
  evtProperties_.printEvtInfo();

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
  map<unsigned int, Photon*> passingWithoutPixelSeed;
  map<unsigned int, Photon*> passingWithPixelSeed;
  unsigned int numCands = 0;
  vector<LeafCandidate::Vector> photonMomentum;
  vector<const Photon*> passingPhotons;
  Handle<PhotonCollection> pPhotons;
  bool foundPhotons = false;
  try { foundPhotons = (iEvent.getByLabel(photonTag_, pPhotons)) && (pPhotons->size() > 0); }
  catch (cms::Exception& ex) {}
  if (!foundPhotons) {
    stringstream infoStream;
    infoStream << "No reco::Photon collection found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n";
    evtProperties_.printDebug(infoStream.str());
  }
  else {

    //find candidates passing preselection, photon ID, and dR(photon, photon)
    if (evtProperties_.foundPhotonCandidates(pPhotons, ECALRecHits, passingWithoutPixelSeed, passingWithPixelSeed)) {

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

      //get the HCAL geometry for calculating the rho and phi of the HB/HE RecHit (there's got to be a better way to do this!)
      iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle_);
      const CaloGeometry* pGeometry = caloGeometryHandle_.product();

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

      //loop over candidate list
      for (map<unsigned int, Photon*>::const_iterator iPassingWithoutPixelSeed = passingWithoutPixelSeed.begin(); 
	   iPassingWithoutPixelSeed != passingWithoutPixelSeed.end(); ++iPassingWithoutPixelSeed) {

	//find candidates failing halo tags
	bool passHalo = false;
	const Photon* photon = const_cast<const Photon*>((*iPassingWithoutPixelSeed).second);
	if (foundHBHERecHits) {
	  if (!evtProperties_.photonIsHEHalo(pHBHERecHits, photon, pGeometry)) {
	    if (foundCosmicTracks) {
	      if (!evtProperties_.photonIsMuonHalo(pCosmicTracks, photon)) passHalo = true;
	    }
	    else passHalo = true;
	  }
	}
	else if (foundCosmicTracks) {
	  if (!evtProperties_.photonIsMuonHalo(pCosmicTracks, photon)) passHalo = true;
	}
	else passHalo = true;

	//find candidates passing all requirements
	bool passAll = false;
	if (passHalo || (!rejectHalo_)) {

	  //find ETRACK candidates passing dR(photon, track)
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
	      if (evtProperties_.electronHasPassingTrack(pTracks, photon, (*iPassingWithoutPixelSeed).first)) {
		passAll = true;
		++numCands;
		//TODO: plot pT of non-overlapping tracks that pass the pT requirement
		//TODO: plot dR(photon, track) of non-overlapping tracks that pass the pT requirement
	      }
	    }
	  }

	  //other samples pass
	  else {
	    passAll = true;
	    ++numCands;
	  }
	}

	//fill histograms
	if (passAll) {
	  if (out_->IsOpen()) {
	    out_->cd();
	    ECALIso_->Fill(photon->ecalRecHitSumEtConeDR04());
	    HCALIso_->Fill(photon->hcalTowerSumEtConeDR04());
	    ET_->Fill(photon->et());
	    HOverE_->Fill(photon->hadronicOverEm());
	    fiducialRegionDist_->Fill(evtProperties_.ECALFiducialRegion(photon));
	    etaWidth_->Fill(photon->sigmaEtaEta());
	    trackIso_->Fill(photon->trkSumPtHollowConeDR04());
	    eta_->Fill(photon->eta());
	    phi_->Fill(photon->phi());
	    photonMomentum.push_back(photon->momentum());
	    passingPhotons.push_back(photon);
	  }
	}
      }//end loop over candidates
    }//end if found photon candidates
  }//end else

  //get PAT AK5 calo jets
  JetIDSelectionFunctor jetIDMinimal(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE);
  pat::strbitset ret = jetIDMinimal.getBitTemplate();
  Handle<vector<pat::Jet> > pPATAK5Jets;
  unsigned int numJets = 0;
  unsigned int numJets30 = 0;
  std::vector<LorentzV> jetLorentzVectors;
  std::vector<LorentzV> jetLorentzVectors30;
  if (getCollection_(pPATAK5Jets, PATAK5CaloJetTag_, iEvent)) {
    for (vector<pat::Jet>::const_iterator iJet = pPATAK5Jets->begin(); iJet != pPATAK5Jets->end(); ++iJet) {
      ret.set(false);
      bool passJetID = jetIDMinimal(*iJet, ret);
      bool noPhotonOverlap = pat::Flags::test(*iJet, pat::Flags::Overlap::Photons);

      //select jets with no photon overlap for analysis (https://twiki.cern.ch/twiki/bin/view/CMS/JetID)
      if (noPhotonOverlap) {

	//determine if jet passes pT cuts
	bool passPTCut = false;
	if (iJet->pt() > minJetPT_) passPTCut = true;
	bool passPTPlus20Cut = false;
	if (iJet->pt() > (minJetPT_ + 20/*GeV*/)) passPTPlus20Cut = true;
	bool passPTPlus40Cut = false;
	if (iJet->pt() > (minJetPT_ + 40/*GeV*/)) passPTPlus40Cut = true;

	//determine if jet passes eta cut
	//only consider jets entirely contained within HB or HE
	bool passEtaCut = false;
	bool HB = (fabs(iJet->eta()) >= 0.0) && (fabs(iJet->eta()) <= 1.3) && (iJet->hadEnergyInHE() == 0.0);
	bool HE = (fabs(iJet->eta()) > 1.3) && (fabs(iJet->eta()) <= 3.0) && (iJet->hadEnergyInHB() == 0.0) && (iJet->hadEnergyInHO() == 0.0);
	if ((fabs(iJet->eta()) <= maxJetAbsEta_) && (fabs(iJet->eta()) >= minJetAbsEta_) && (HB || HE)) passEtaCut = true;

	//jets passing LOOSE PURE09 ID
	if (passJetID) {

	  //fill N - 1 plots of jet pT and eta
	  if (out_->IsOpen()) {
	    out_->cd();
	    if (passEtaCut) {
	      jetPTJetIDApplied_->Fill(iJet->pt());
	      jetCorrVsPTJetIDApplied_->Fill(iJet->pt(), iJet->corrFactor("raw"));
	      jetCorrJetIDApplied_->Fill(iJet->corrFactor("raw"));
	      jetPTRawJetIDApplied_->Fill(iJet->correctedJet("raw").pt());
	    }
	    if (passPTCut) jetEtaJetIDApplied_->Fill(fabs(iJet->eta()));
	    if (passPTPlus20Cut) jetEtaPlus20JetIDApplied_->Fill(fabs(iJet->eta()));
	    if (passPTPlus40Cut) jetEtaPlus40JetIDApplied_->Fill(fabs(iJet->eta()));
	  }

	  //count jets
	  if (passPTCut && passEtaCut) {
	    ++numJets;
	    jetLorentzVectors.push_back(iJet->p4());
	  }
	  if (iJet->pt() > 30.0/*GeV*/)	{
	    ++numJets30;
	    jetLorentzVectors30.push_back(iJet->p4());
	  }

	}

	//jets failing LOOSE PURE09 ID
	else {

	  //fill N - 1 plots of jet pT and eta
	  if (out_->IsOpen()) {
	    out_->cd();
	    if (passEtaCut) {
	      jetPTJetIDNotApplied_->Fill(iJet->pt());
	      jetCorrVsPTJetIDNotApplied_->Fill(iJet->pt(), iJet->corrFactor("raw"));
	      jetCorrJetIDNotApplied_->Fill(iJet->corrFactor("raw"));
	      jetPTRawJetIDNotApplied_->Fill(iJet->correctedJet("raw").pt());
	    }
	    if (passPTCut) jetEtaJetIDNotApplied_->Fill(fabs(iJet->eta()));
	    if (passPTPlus20Cut) jetEtaPlus20JetIDNotApplied_->Fill(fabs(iJet->eta()));
	    if (passPTPlus40Cut) jetEtaPlus40JetIDNotApplied_->Fill(fabs(iJet->eta()));
	  }
	}
      }
    }
    if (out_->IsOpen()) {
      out_->cd();
      numJets_->Fill(numJets);
      numJets30_->Fill(numJets30);
      //evtAlphaT_->Fill(alphaT(jetLorentzVectors));
      //evtAlphaT50_->Fill(alphaT(jetLorentzVectors50));
    }
  }

  //get the jet ID variables
  JetIDSelectionFunctor jetIDFunctor(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE);
  pat::strbitset retRECO = jetIDFunctor.getBitTemplate();
  Handle<JetIDValueMap> pJetIDMap;
  if (getCollection_(pJetIDMap, AK5JetIDTag_, iEvent)) {

    //get RECO AK5 calo jets
    Handle<View<CaloJet> > pRECOAK5Jets;
    if (getCollection_(pRECOAK5Jets, RECOAK5CaloJetTag_, iEvent)) {
      for (View<CaloJet>::const_iterator iJet = pRECOAK5Jets->begin(); iJet != pRECOAK5Jets->end(); ++iJet) {

	//determine if the jet passes LOOSE PURE09 jet ID
	RefToBase<CaloJet> jetRef = pRECOAK5Jets->refAt(iJet - pRECOAK5Jets->begin());
	JetID const& jetID = (*pJetIDMap)[jetRef];
	retRECO.set(false);
	bool passed = jetIDFunctor(*iJet, jetID, retRECO);
	CaloJet const& jet = *jetRef;

	//determine if jet passes pT cuts
	bool passPTCut = false;
	if (jet.pt() > minJetPT_) passPTCut = true;
	bool passPTPlus20Cut = false;
	if (jet.pt() > (minJetPT_ + 20/*GeV*/)) passPTPlus20Cut = true;
	bool passPTPlus40Cut = false;
	if (jet.pt() > (minJetPT_ + 40/*GeV*/)) passPTPlus40Cut = true;

	//determine if jet passes eta cut
	//only consider jets entirely contained within HB or HE
	bool passEtaCut = false;
	bool HB = (fabs(jet.eta()) >= 0.0) && (fabs(jet.eta()) <= 1.3) && (jet.hadEnergyInHE() == 0.0);
	bool HE = (fabs(jet.eta()) > 1.3) && (fabs(jet.eta()) <= 3.0) && (jet.hadEnergyInHB() == 0.0) && (iJet->hadEnergyInHO() == 0.0);
	if ((fabs(jet.eta()) <= maxJetAbsEta_) && (fabs(jet.eta()) >= minJetAbsEta_) && (HB || HE)) passEtaCut = true;

	//determine if jet overlaps with a photon
	bool noPhotonOverlap = true;
	PhotonCollection::const_iterator iPhoton = pPhotons->begin();
	while ((iPhoton != pPhotons->end()) && (noPhotonOverlap)) {
	  if (evtProperties_.dR(jet.eta(), iPhoton->eta(), jet.phi(), iPhoton->phi()) < 0.5) noPhotonOverlap = false;
	  ++iPhoton;
	}
	if (noPhotonOverlap) {

	  //jets passing LOOSE PURE09 jet ID
	  if (passed) {

	    //fill N - 1 plots of jet pT and eta
	    if (out_->IsOpen()) {
	      out_->cd();
	      if (passEtaCut) {
		//jetPTJetIDApplied_->Fill(iJet->pt());
		//jetCorrVsPTJetIDApplied_->Fill(iJet->pt(), iJet->corrFactor("raw"));
		//jetCorrJetIDApplied_->Fill(iJet->corrFactor("raw"));
		jetPTRECOJetIDApplied_->Fill(jet.pt());
	      }
	      //if (passPTCut) jetEtaJetIDApplied_->Fill(fabs(iJet->eta()));
	      //if (passPTPlus20Cut) jetEtaPlus20JetIDApplied_->Fill(fabs(iJet->eta()));
	      //if (passPTPlus40Cut) jetEtaPlus40JetIDApplied_->Fill(fabs(iJet->eta()));
	    }
	  }

	  //jets failing LOOSE PURE09 jet ID
	  else {

	    //fill N - 1 plots of jet pT and eta
	    if (out_->IsOpen()) {
	      out_->cd();
	      if (passEtaCut) {
		//jetPTJetIDApplied_->Fill(iJet->pt());
		//jetCorrVsPTJetIDApplied_->Fill(iJet->pt(), iJet->corrFactor("raw"));
		//jetCorrJetIDApplied_->Fill(iJet->corrFactor("raw"));
		jetPTRECOJetIDNotApplied_->Fill(jet.pt());
	      }
	      //if (passPTCut) jetEtaJetIDApplied_->Fill(fabs(iJet->eta()));
	      //if (passPTPlus20Cut) jetEtaPlus20JetIDApplied_->Fill(fabs(iJet->eta()));
	      //if (passPTPlus40Cut) jetEtaPlus40JetIDApplied_->Fill(fabs(iJet->eta()));
	    }
	  }
	}
      }
    }
  }

  //calculate event quantities
  if (numCands >= 2) {

    //di-EM pT
    LeafCandidate::Vector diEMPTVec = photonMomentum[0];
    diEMPTVec+=photonMomentum[1];
    const double diEMPT = sqrt(diEMPTVec.Perp2());
    if (out_->IsOpen()) {
      out_->cd();
      diEMPT_->Fill(diEMPT);
    }
    diEMPTVal_ = diEMPT;

    //idiot check
    const unsigned int numPassing = passingPhotons.size();
    if (numPassing < 2) cerr << "Error: passingPhotons.size() = " << numPassing << ".\n";
    else {

      //angles between candidates
      const double dPhiPhotons = evtProperties_.dPhi(passingPhotons[0]->phi(), passingPhotons[1]->phi());
      if (out_->IsOpen()) {
	out_->cd();
	openingAngle_->Fill(acos(photonMomentum[0].Dot(photonMomentum[1])/(sqrt(photonMomentum[0].Mag2())*sqrt(photonMomentum[1].Mag2()))));
	dPhiCands_->Fill(dPhiPhotons);
      }

      //get the leading and trailing photons
      double cand1ET = -1.0;
      unsigned int cand1Index = 0;
      for (vector<const Photon*>::const_iterator iCand = passingPhotons.begin(); iCand != passingPhotons.end(); ++iCand) {
	const double ET = (*iCand)->et();
	if (ET > cand1ET) {
	  cand1ET = ET;
	  cand1Index = iCand - passingPhotons.begin();
	}
      }
      double cand2ET = -1.0;
      unsigned int cand2Index = 0;
      for (vector<const Photon*>::const_iterator iCand = passingPhotons.begin(); iCand != passingPhotons.end(); ++iCand) {
	const double ET = (*iCand)->et();
	const unsigned int index = iCand - passingPhotons.begin();
	if ((ET > cand2ET) && (index != cand1Index)) {
	  cand2ET = ET;
	  cand2Index = index;
	}
      }

      //angles between caloMET and candidates
      Handle<CaloMETCollection> pCaloMET;
      bool caloMETFound = false;
      try { caloMETFound = (iEvent.getByLabel(caloMETTag_, pCaloMET)) && (pCaloMET->size() > 0); }
      catch (cms::Exception& ex) {}
      if (!caloMETFound) {
	stringstream infoStream;
	infoStream << "No reco::CaloMET collection found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n";
	evtProperties_.printDebug(infoStream.str());
      }
      else {
	const double caloMETPhi = pCaloMET->begin()->phi();
	for (vector<const Photon*>::const_iterator iCand = passingPhotons.begin(); iCand != passingPhotons.end(); ++iCand) {
	  if (out_->IsOpen()) {
	    out_->cd();
	    dPhiMETCand_->Fill(evtProperties_.dPhi(caloMETPhi, (*iCand)->phi()));
	  }
	}

	LeafCandidate::Vector caloMETMomentum = pCaloMET->begin()->momentum();
	if (out_->IsOpen()) {
	  out_->cd();
	  dPhiMETCand1_->Fill(evtProperties_.dPhi(caloMETPhi, passingPhotons[cand1Index]->phi()));
	  dPhiMETCand2_->Fill(evtProperties_.dPhi(caloMETPhi, passingPhotons[cand2Index]->phi()));

	  //angles between MET and di-EM pT
	  openingAngleMETDiEMPT_->Fill(acos(caloMETMomentum.Dot(diEMPTVec)/(sqrt(caloMETMomentum.Mag2())*sqrt(diEMPTVec.Mag2()))));
	  dPhiMETDiEMPT_->Fill(evtProperties_.dPhi(caloMETPhi, diEMPTVec.Phi()));
	}
      }

      //tcMET
      Handle<METCollection> pTCMET;
      bool tcMETFound = false;
      try { tcMETFound = (iEvent.getByLabel(tcMETTag_, pTCMET)) && (pTCMET->size() > 0); }
      catch (cms::Exception& ex) {}
      if (!tcMETFound) {
	stringstream infoStream;
	infoStream << "No tcMET collection found in run " << runNum << ", event " << evtNum << ", lumi section " << lumiNum << ".\n";
	evtProperties_.printDebug(infoStream.str());
      }
      else tcMET_ = pTCMET->begin()->et();

      //alphaT bins:
      //2 photons 0 jets
      //2 photons 1 jet
      //2 photons 2 jets
      //2 photons 3 jets
      //3 photons 0 jets
      //3 photons 1 jet
      //3 photons 2 jets

      //then, add MHT ratio requirement

      //alphaT
      map<string, double> alphaTMap;
      alphaTMap["2 photons 0 jets"] = alphaT(jetLorentzVectors, passingPhotons, 0, 2);
      alphaTMap["2 photons 1 jets"] = alphaT(jetLorentzVectors, passingPhotons, 1, 2);
      alphaTMap["2 photons 2 jets"] = alphaT(jetLorentzVectors, passingPhotons, 2, 2);
      alphaTMap["2 photons 3 jets"] = alphaT(jetLorentzVectors, passingPhotons, 3, 2);
      alphaTMap["3 photons 0 jets"] = alphaT(jetLorentzVectors, passingPhotons, 0, 3);
      alphaTMap["3 photons 1 jets"] = alphaT(jetLorentzVectors, passingPhotons, 1, 3);
      alphaTMap["3 photons 2 jets"] = alphaT(jetLorentzVectors, passingPhotons, 2, 3);
      if (out_->IsOpen()) {
	out_->cd();
	out_->cd("alphaT");
	alphaT2Photons0Jets_->Fill(alphaTMap["2 photons 0 jets"]);
	alphaT2Photons1Jets_->Fill(alphaTMap["2 photons 1 jets"]);
	alphaT2Photons2Jets_->Fill(alphaTMap["2 photons 2 jets"]);
	alphaT2Photons3Jets_->Fill(alphaTMap["2 photons 3 jets"]);
	alphaT3Photons0Jets_->Fill(alphaTMap["3 photons 0 jets"]);
	alphaT3Photons1Jets_->Fill(alphaTMap["3 photons 1 jets"]);
	alphaT3Photons2Jets_->Fill(alphaTMap["3 photons 2 jets"]);
	std::vector<LorentzV> dummy;
	vector<const Photon*> photons;
	photons.push_back(passingPhotons[0]);
	photons.push_back(passingPhotons[1]);
	const double alphaTPhotons = alphaT(dummy, photons, 0, 2);
	alphaTPhotons_->Fill(alphaTPhotons);
	out_->cd("..");
	alphaTVsNJetsSelected_->Fill(numJets, alphaTPhotons);
      }
    }
  }

  //2 cands not found
  else {
    diEMPTVal_ = -1.0;
    tcMET_ = -1.0;
  }

  numCandsPerEvt_->Fill(numCands);
  //TODO: fill numTracksPerEvt_

  //fill event tree
  if (out_->IsOpen()) {
    out_->cd();
    evtTree_->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
CutAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
  if (!out_->IsOpen() && debugFlag_) debug_ << "Error opening file " << outFileName_ << ".\n";
  else {
    out_->cd();
    out_->mkdir("alphaT");
  }

  //book histograms
  ECALIso_ = new TH1F("ECALIso_", "ECAL isolation of passing e/#gamma candidates;ECAL isolation (GeV);e/#gamma candidates per GeV", 50, 
		      0.0, 50.0);
  HCALIso_ = new TH1F("HCALIso_", "HCAL isolation of passing e/#gamma candidates;HCAL isolation (GeV);e/#gamma candidates per GeV", 50, 
		      0.0, 50.0);
  ET_ = new TH1D("ET_", "E_{T} of passing e/#gamma candidates;E_{T} (GeV);e/#gamma candidates per GeV", 100, 0.0, 100.0);
  HOverE_ = new TH1F("HOverE_", "H/E of passing e/#gamma candidates;H/E;e/#gamma candidates per 0.01", 100, 0.0, 1.0);
  fiducialRegionDist_ = new TH1F("fiducialRegion_", 
				 "Fiducial region of passing e/#gamma candidates;;e/#gamma candidates per ECAL fiducial region", 2, 0.5, 
				 2.5);
  fiducialRegionDist_->GetXaxis()->SetBinLabel(EB, "EB");
  fiducialRegionDist_->GetXaxis()->SetBinLabel(EEND, "EE");
  etaWidth_ = new TH1F("etaWidth_", "#sigma_{#eta#eta} of passing e/#gamma candidates;#sigma_{#eta#eta};e/#gamma candidates per 0.001", 
		       100, 0.0, 0.1);
  trackIso_ = new TH1F("trackIso_", "Track isolation of passing e/#gamma candidates;Track isolation (GeV);e/#gamma candidates per GeV", 
		       50, 0.0, 50.0);
  numCandsPerEvt_ = new TH1F("numCandsPerEvt_", "Number of passing candidates per event;Number per event;Events per 1", 5, -0.5, 4.5);
  if (sampleType_ == ETRACK) {
    trackPT_ = new TH1D("trackPT_", "p_{T} of passing tracks;p_{T} (GeV);Tracks per GeV", 100, 0.0, 100.0);
    dRTrackPhoton_ = new TH1D("dRTrackPhoton_", "#DeltaR(track, e/#gamma candidate) of passing tracks-e/#gamma candidate pairs;#DeltaR;", 100, 0.0, 1.0);
    dRTrackPhoton_->GetYaxis()->SetTitle("track-e/#gamma candidate pairs per 0.01");
  }
  eta_ = new TH1D("eta_", "#eta of passing e/#gamma candidates;#eta;e/#gamma candidates per 0.2", 15, -1.5, 1.5);
  phi_ = new TH1D("phi_", "#phi of passing e/#gamma candidates;#phi;e/#gamma candidates per 0.4", 16, -3.2, 3.2);
  diEMPT_ = new TH1D("diEMPT_", "di-EM p_{T} of first two passing e/#gamma candidates;di-EM p_{T} (GeV);Events per GeV", 100, 0.0, 100.0);
  openingAngle_ = new TH1D("openingAngle_", "Angle between 2 leading e/#gamma candidates;Angle (rad);Events/0.4", 8, 0.0, 3.2);
  dPhiCands_ = new TH1D("dPhiCands_", "#Delta#phi between the 2 leading e/#gamma candidates;#Delta#phi (rad);Events/0.4", 8, 0.0, 3.2);
  dPhiMETCand1_ = new TH1D("dPhiMETCand1_", "#Delta#phi between the ME_{T} and the leading e/#gamma candidate;#Delta#phi (rad);Events/0.4", 8, 0.0, 3.2);
  dPhiMETCand2_ = new TH1D("dPhiMETCand2_", "#Delta#phi between the ME_{T} and the trailing e/#gamma candidate;#Delta#phi (rad);Events/0.4", 8, 0.0, 3.2);
  dPhiMETCand_ = new TH1D("dPhiMETCand_", "#Delta#phi between the ME_{T} and the e/#gamma candidate;#Delta#phi (rad);e/#gamma candidates/0.4", 8, 0.0, 3.2);
  openingAngleMETDiEMPT_ = new TH1D("openingAngleMETDiEMPT_", "Angle between the ME_{T} and di-EM p_{T} vector;Angle (rad);Events/0.4", 8, 0.0, 3.2);
  dPhiMETDiEMPT_ = new TH1D("dPhiMETDiEMPT_", "#Delta#phi between the ME_{T} and the di-EM p_{T} vector;#Delta#phi (rad);Events/0.4", 8, 0.0, 3.2);
  numJets_ = new TH1F("numJets_", "numJets_", MAX_NUM_JETS, -0.5, -0.5 + MAX_NUM_JETS);
  numJets30_ = new TH1F("numJets30_", "numJets30_", MAX_NUM_JETS, -0.5, -0.5 + MAX_NUM_JETS);
  jetPTJetIDApplied_ = new TH1D("jetPTJetIDApplied_", "jetPTJetIDApplied_;p_{T} (GeV);Jets/5 GeV", 40, 0.0, 200.0);
  jetPTJetIDNotApplied_ = new TH1D("jetPTJetIDNotApplied_", "jetPTJetIDNotApplied_;p_{T} (GeV);Jets/5 GeV", 40, 0.0, 200.0);
  jetEtaJetIDApplied_ = new TH1D("jetEtaJetIDApplied_", "jetEtaJetIDApplied_;#eta;Jets/0.1", 30, 0.0, 3.0);
  jetEtaPlus20JetIDApplied_ = new TH1D("jetEtaPlus20JetIDApplied_", "jetEtaPlus20JetIDApplied_;#eta;Jets/0.1", 30, 0.0, 3.0);
  jetEtaPlus40JetIDApplied_ = new TH1D("jetEtaPlus40JetIDApplied_", "jetEtaPlus40JetIDApplied_;#eta;Jets/0.1", 30, 0.0, 3.0);
  jetEtaJetIDNotApplied_ = new TH1D("jetEtaJetIDNotApplied_", "jetEtaJetIDNotApplied_;#eta;Jets/0.1", 30, 0.0, 3.0);
  jetEtaPlus20JetIDNotApplied_ = new TH1D("jetEtaPlus20JetIDNotApplied_", "jetEtaPlus20JetIDNotApplied_;#eta;Jets/0.1", 30, 0.0, 3.0);
  jetEtaPlus40JetIDNotApplied_ = new TH1D("jetEtaPlus40JetIDNotApplied_", "jetEtaPlus40JetIDNotApplied_;#eta;Jets/0.1", 30, 0.0, 3.0);
  alphaT2Photons0Jets_ = new TH1D("alphaT2Photons0Jets_", "alphaT2Photons0Jets_", 20, 0.0, 1.0);
  alphaT2Photons1Jets_ = new TH1D("alphaT2Photons1Jets_", "alphaT2Photons1Jets_", 20, 0.0, 1.0);
  alphaT2Photons2Jets_ = new TH1D("alphaT2Photons2Jets_", "alphaT2Photons2Jets_", 20, 0.0, 1.0);
  alphaT2Photons3Jets_ = new TH1D("alphaT2Photons3Jets_", "alphaT2Photons3Jets_", 20, 0.0, 1.0);
  alphaT3Photons0Jets_ = new TH1D("alphaT3Photons0Jets_", "alphaT3Photons0Jets_", 20, 0.0, 1.0);
  alphaT3Photons1Jets_ = new TH1D("alphaT3Photons1Jets_", "alphaT3Photons1Jets_", 20, 0.0, 1.0);
  alphaT3Photons2Jets_ = new TH1D("alphaT3Photons2Jets_", "alphaT3Photons2Jets_", 20, 0.0, 1.0);
  alphaTPhotons_ = new TH1D("alphaTPhotons_", "#alpha_{T} between the 2 photons;#alpha_{T};Events/0.05", 20, 0.0, 1.0);
  stringstream histTitle;
  histTitle.str("");
  histTitle << "#alpha_{T} between the 2 photons vs. number of LOOSE PURE09 AK5 calo jets per event, excluding photons, jet p_{T} > " << minJetPT_;
  histTitle << ", jet |#eta| < " << maxJetAbsEta_ << ";Number per event;#alpha_{T}";
  alphaTVsNJetsSelected_ = new TH2D("alphaTVsNJetsSelected_", histTitle.str().c_str(), MAX_NUM_JETS, -0.5, -0.5 + MAX_NUM_JETS, 40, 0.0, 2.0);
  histTitle.str("");
  histTitle << "#Event alpha_{T} calculated from LOOSE PURE09 AK5 calo jets per event, excluding photons, jet p_{T} > " << minJetPT_;
  histTitle << ", jet |#eta| < " << maxJetAbsEta_ << ";Number per event;#alpha_{T}";
  evtAlphaT_ = new TH1D("evtAlphaT_", histTitle.str().c_str(), 40, 0.0, 2.0);
  histTitle.str("");
  histTitle << "#Event alpha_{T} calculated from LOOSE PURE09 AK5 calo jets per event, excluding photons, jet p_{T} > 50 GeV, jet |#eta| < ";
  histTitle << maxJetAbsEta_ << ";Number per event;#alpha_{T}";
  evtAlphaT50_ = new TH1D("evtAlphaT50_", histTitle.str().c_str(), 40, 0.0, 2.0);
  hadEnergyInHBJetIDApplied_ = new TH1F("hadEnergyInHBJetIDApplied_", "hadEnergyInHBJetIDApplied_", 40, 0.0, 200.0);
  hadEnergyInHEJetIDApplied_ = new TH1F("hadEnergyInHEJetIDApplied_", "hadEnergyInHEJetIDApplied_", 40, 0.0, 200.0);
  hadEnergyInHBJetIDNotApplied_ = new TH1F("hadEnergyInHBJetIDNotApplied_", "hadEnergyInHBJetIDNotApplied_", 40, 0.0, 200.0);
  hadEnergyInHEJetIDNotApplied_ = new TH1F("hadEnergyInHEJetIDNotApplied_", "hadEnergyInHEJetIDNotApplied_", 40, 0.0, 200.0);
  hadEnergyInHBVsEtaJetIDApplied_ = new TH2F("hadEnergyInHBVsEtaJetIDApplied_", "hadEnergyInHBVsEtaJetIDApplied_", 30, 0.0, 3.0, 40, 0.0, 200.0);
  hadEnergyInHEVsEtaJetIDApplied_ = new TH2F("hadEnergyInHEVsEtaJetIDApplied_", "hadEnergyInHEVsEtaJetIDApplied_", 30, 0.0, 3.0, 40, 0.0, 200.0);
  hadEnergyInHBVsEtaJetIDNotApplied_ = new TH2F("hadEnergyInHBVsEtaJetIDNotApplied_", "hadEnergyInHBVsEtaJetIDNotApplied_", 30, 0.0, 3.0, 40, 0.0, 200.0);
  hadEnergyInHEVsEtaJetIDNotApplied_ = new TH2F("hadEnergyInHEVsEtaJetIDNotApplied_", "hadEnergyInHEVsEtaJetIDNotApplied_", 30, 0.0, 3.0, 40, 0.0, 200.0);
  jetCorrVsPTJetIDApplied_ = new TH2D("jetCorrVsPTJetIDApplied_", "jetCorrVsPTJetIDApplied_;p_{T} (GeV);Correction factor", 200, 0.0, 200.0, 200, 0.0, 2.0);
  jetCorrVsPTJetIDNotApplied_ = new TH2D("jetCorrVsPTJetIDNotApplied_", "jetCorrVsPTJetIDNotApplied_;p_{T} (GeV);Correction factor", 200, 0.0, 200.0, 200, 
					 0.0, 2.0);
  jetCorrJetIDApplied_ = new TH1F("jetCorrJetIDApplied_", "jetCorrJetIDApplied_", 200, 0.0, 2.0);
  jetCorrJetIDNotApplied_ = new TH1F("jetCorrJetIDNotApplied_", "jetCorrJetIDNotApplied_", 200, 0.0, 2.0);
  jetPTRawJetIDApplied_ = new TH1D("jetPTRawJetIDApplied_", "jetPTRawJetIDApplied_", 40, 0.0, 200.0);
  jetPTRawJetIDNotApplied_ = new TH1D("jetPTRawJetIDNotApplied_", "jetPTRawJetIDNotApplied_", 40, 0.0, 200.0);
  jetPTRECOJetIDApplied_ = new TH1D("jetPTRECOJetIDApplied_", "jetPTRECOJetIDApplied_;p_{T} (GeV);", 40, 0.0, 200.0);
  jetPTRECOJetIDNotApplied_ = new TH1D("jetPTRECOJetIDNotApplied_", "jetPTRECOJetIDNotApplied_;p_{T} (GeV);", 40, 0.0, 200.0);
  
  //instantiate tree
  run_ = 0;
  evt_ = 0;
  lumiSection_ = 0;
  diEMPTVal_ = 0.0;
  tcMET_ = 0.0;
  evtTree_ = new TTree("evtTree_", "evtTree_");
  evtTree_->Branch("run", &run_, "run/i");
  evtTree_->Branch("evt", &evt_, "evt/i");
  evtTree_->Branch("lumiSection", &lumiSection_, "lumiSection/i");
  evtTree_->Branch("diEMPTVal", &diEMPTVal_, "diEMPTVal/D");
  evtTree_->Branch("tcMET", &tcMET_, "tcMET/D");
  cout << "Run: " << run_ << ", event: " << evt_ << ", lumi section: " << lumiSection_ << ", di-EM pT: " << diEMPTVal_ << " GeV, tcMET: " << tcMET_;
  cout << " GeV\n";
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CutAnalyzer::endJob() {

  //create alphaT fractions graph
  //implementation depends on alphaTVsNJetsSelected having no bins that straddle alphaT = 0.55
  double x[MAX_NUM_JETS];
  double xErr[MAX_NUM_JETS];
  double y[MAX_NUM_JETS];
  double yErr[MAX_NUM_JETS];
  for (unsigned int i = 0; i < MAX_NUM_JETS; ++i) {
    x[i] = i;
    xErr[i] = 0.0;
    unsigned int numAlphaTGreaterThan055 = 0;
    double numAlphaTGreaterThan055Err2 = 0.0;
    unsigned int numAlphaTLessThan055 = 0;
    double numAlphaTLessThan055Err2 = 0.0;
    for (unsigned int j = 1; j <= (unsigned int)alphaTVsNJetsSelected_->GetNbinsY(); ++j) {
      const unsigned int binContent = alphaTVsNJetsSelected_->GetBinContent(i + 1, j);
      const double binErr = alphaTVsNJetsSelected_->GetBinError(i + 1, j);
      if ((alphaTVsNJetsSelected_->GetYaxis()->GetBinLowEdge(j) + alphaTVsNJetsSelected_->GetYaxis()->GetBinWidth(j)) <= 0.55) {
	numAlphaTLessThan055+=binContent;
	numAlphaTLessThan055Err2+=binErr*binErr;
      }
      else {
	numAlphaTGreaterThan055+=binContent;
	numAlphaTGreaterThan055Err2+=binErr*binErr;
      }
    }
    if (numAlphaTLessThan055 > 0) {
      y[i] = (double)(numAlphaTGreaterThan055/numAlphaTLessThan055);
      yErr[i] = y[i]*sqrt((numAlphaTGreaterThan055Err2/(numAlphaTGreaterThan055*numAlphaTGreaterThan055)) + 
			  (numAlphaTLessThan055Err2/(numAlphaTLessThan055*numAlphaTLessThan055)));
    }
    else {
      y[i] = 0.0;
      yErr[i] = 0.0;
    }
  }
  alphaTFraction_ = new TGraphErrors(MAX_NUM_JETS, x, y, xErr, yErr);
  alphaTFractionCanvas_ = new TCanvas("alphaTFractionCanvas_", "alphaTFractionCanvas_", 600, 600);
  alphaTFractionCanvas_->SetWindowSize(1200 - alphaTFractionCanvas_->GetWw(), 1200 - alphaTFractionCanvas_->GetWh());
  alphaTFractionCanvas_->cd();
  alphaTFraction_->Draw("A*");

  //write histograms to file
  if (out_->IsOpen()) {
    out_->cd();
    ECALIso_->Write();
    HCALIso_->Write();
    ET_->Write();
    HOverE_->Write();
    fiducialRegionDist_->Write();
    etaWidth_->Write();
    trackIso_->Write();
    eta_->Write();
    phi_->Write();
    if (sampleType_ == ETRACK) {
      trackPT_->Write();
      dRTrackPhoton_->Write();
    }
    numCandsPerEvt_->Write();
    diEMPT_->Write();
    openingAngle_->Write();
    dPhiCands_->Write();
    dPhiMETCand1_->Write();
    dPhiMETCand2_->Write();
    dPhiMETCand_->Write();
    openingAngleMETDiEMPT_->Write();
    dPhiMETDiEMPT_->Write();
    numJets_->Write();
    numJets30_->Write();
    jetPTJetIDApplied_->Write();
    jetPTJetIDNotApplied_->Write();
    jetEtaJetIDApplied_->Write();
    jetEtaPlus20JetIDApplied_->Write();
    jetEtaPlus40JetIDApplied_->Write();
    jetEtaJetIDNotApplied_->Write();
    jetEtaPlus20JetIDNotApplied_->Write();
    jetEtaPlus40JetIDNotApplied_->Write();
    alphaT2Photons0Jets_->Write();
    alphaT2Photons1Jets_->Write();
    alphaT2Photons2Jets_->Write();
    alphaT2Photons3Jets_->Write();
    alphaT3Photons0Jets_->Write();
    alphaT3Photons1Jets_->Write();
    alphaT3Photons2Jets_->Write();
    alphaTPhotons_->Write();
    alphaTVsNJetsSelected_->Write();
    evtAlphaT_->Write();
    evtAlphaT50_->Write();
    alphaTFractionCanvas_->Write();
    hadEnergyInHBJetIDApplied_->Write();
    hadEnergyInHEJetIDApplied_->Write();
    hadEnergyInHBJetIDNotApplied_->Write();
    hadEnergyInHEJetIDNotApplied_->Write();
    hadEnergyInHBVsEtaJetIDApplied_->Write();
    hadEnergyInHEVsEtaJetIDApplied_->Write();
    hadEnergyInHBVsEtaJetIDNotApplied_->Write();
    hadEnergyInHEVsEtaJetIDNotApplied_->Write();
    jetCorrVsPTJetIDApplied_->Write();
    jetCorrVsPTJetIDNotApplied_->Write();
    jetCorrJetIDApplied_->Write();
    jetCorrJetIDNotApplied_->Write();
    alphaT2Photons0Jets_->Write();
    jetPTRawJetIDApplied_->Write();
    jetPTRawJetIDNotApplied_->Write();
    jetPTRECOJetIDApplied_->Write();
    jetPTRECOJetIDNotApplied_->Write();
    evtTree_->Write();
    out_->Write();
    out_->Close();
  }
  delete out_;
}

double CutAnalyzer::alphaT(const vector<LorentzV>& jetLorentzVectors, vector<const Photon*>& passingPhotons, const unsigned int numJets, 
			   const unsigned int numPhotons) const
{
  //return if wrong number of objects
  if ((jetLorentzVectors.size() != numJets) || (passingPhotons.size() != numPhotons)) return -1.0;

  //make vectors of object pT and calculate x and y components of the multi-object system
  vector<double> pT;
  double sumPX = 0.0;
  double sumPY = 0.0;
  for (vector<LorentzV>::const_iterator iJet = jetLorentzVectors.begin(); iJet != jetLorentzVectors.end(); ++iJet) {
    pT.push_back((*iJet).pt());
    cout << "pT jet " << (iJet - jetLorentzVectors.begin()) << ": " << (*iJet).pt() << endl;
    sumPX+=(*iJet).px();
    cout << "px jet " << (iJet - jetLorentzVectors.begin()) << ": " << (*iJet).px() << endl;
    sumPY+=(*iJet).py();
    cout << "py jet " << (iJet - jetLorentzVectors.begin()) << ": " << (*iJet).py() << endl;
  }
  for (vector<const Photon*>::const_iterator iPhoton = passingPhotons.begin(); iPhoton != passingPhotons.end(); ++iPhoton) {
    pT.push_back((*iPhoton)->pt());
    cout << "pT photon " << (iPhoton - passingPhotons.begin()) << ": " << (*iPhoton)->pt() << endl;
    sumPX+=(*iPhoton)->px();
    cout << "px photon " << (iPhoton - passingPhotons.begin()) << ": " << (*iPhoton)->px() << endl;
    sumPY+=(*iPhoton)->py();
    cout << "py photon " << (iPhoton - passingPhotons.begin()) << ": " << (*iPhoton)->py() << endl;
  }

  //calculate minimum dHT
  double minDHT = 999.0;
  vector<double> pseudojet1;
  for (unsigned int i = 0; i < pT.size(); ++i) { pseudojet1.push_back(-1.0); }
  minDHTRecursive(minDHT, pseudojet1, pT, 0);
  if (minDHT == 999.0) return -1.0;
  cout << "minDHT: " << minDHT << endl;

  //calculate HT
  const double HT = accumulate(pT.begin(), pT.end(), 0.0);
  cout << "HT: " << HT << endl;

  //calculate MHT squared
  const double MHT2 = sumPX*sumPX + sumPY*sumPY;
  cout << "MHT2: " << MHT2 << endl;

  //return alphaT
  cout << "alphaT: " << 0.5*((HT - minDHT)/sqrt(HT*HT - MHT2)) << endl;
  return 0.5*((HT - minDHT)/sqrt(HT*HT - MHT2));
}

//calculate minimum dHT for an arbitrary number of objects
//the answer is in minDHT
void CutAnalyzer::minDHTRecursive(double& minDHT, vector<double>& pseudojet1, const vector<double>& pT, const unsigned int iteration) const
{
  if (pseudojet1.size() != pT.size()) return;
  const unsigned int numObjects = pT.size();
  for (unsigned int i = iteration; i < numObjects; ++i) {
    pseudojet1[i] = pT[i];
    double pseudojet1HT = 0.0;
    double pseudojet2HT = 0.0;
    for (vector<double>::const_iterator iJet = pseudojet1.begin(); iJet != pseudojet1.end(); ++iJet) {
      if (*iJet == -1.0) {
	pseudojet1HT+=0.0;
	pseudojet2HT+=pT[iJet - pseudojet1.begin()];
      }
      else {
	pseudojet1HT+=*iJet;
	pseudojet2HT+=0.0;
      }
    }
    const double dHT = fabs(pseudojet1HT - pseudojet2HT);
    if (dHT < minDHT) minDHT = dHT;
    minDHTRecursive(minDHT, pseudojet1, pT, i + 1);
    pseudojet1[i] = -1.0;
  }
}

//for alphaT calculation
/*std::vector<double> CutAnalyzer::deltaSumPt_permutations(const std::vector<LorentzV>& p4s) {
   std::vector<std::vector<double> > ht(1 << (p4s.size() - 1), std::vector<double>(2, 0.));
   for (unsigned i = 0; i < ht.size(); i++) {
      for (unsigned j = 0; j < p4s.size(); j++) {
         ht[i][(i / (1 << j)) % 2] += p4s[j].pt();
      }
   }
   std::vector<double> deltaHT;
   for (unsigned i = 0; i < ht.size(); i++)
      deltaHT.push_back(fabs(ht[i][0] - ht[i][1]));
   return deltaHT;
}

//for alphaT calculation
double CutAnalyzer::alphaT(const std::vector<LorentzV>& p4s) {
   if (p4s.size() < 2)
      return 0;

   std::vector<double> pTs;
   for (unsigned i = 0; i < p4s.size(); i++)
      pTs.push_back(p4s[i].pt());
   for (unsigned i = 0; i < p4s.size(); i++)
      pTs.push_back(p4s[i].pt());
   const std::vector<double> DHTper(deltaSumPt_permutations(p4s));

   const double mDHT = *(std::min_element(DHTper.begin(), DHTper.end()));
   const double sumPT = accumulate(pTs.begin(), pTs.end(), double(0));
   const LorentzV sumP4 = accumulate(p4s.begin(), p4s.end(), LorentzV());

   return 0.5 * (sumPT - mDHT) / sqrt(sumPT * sumPT - sumP4.perp2());
   }*/

//define this as a plug-in
DEFINE_FWK_MODULE(CutAnalyzer);
