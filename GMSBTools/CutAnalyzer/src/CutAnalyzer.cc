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
#include "GMSBTools/EventSelection/interface/EventSelector.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

//ROOT include files
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"

//relevant namespaces
using namespace edm;
using namespace cms;
using namespace std;
using namespace reco;

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
  double maxSeedTime_;
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
  ESHandle<CaloGeometry> caloGeometryHandle_;

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
  maxSeedTime_(iConfig.getParameter<double>("maxSeedTime")),
  rejectHalo_(iConfig.getUntrackedParameter<bool>("rejectHalo", true)),
  photonTag_(iConfig.getParameter<InputTag>("photonTag")),
  trackTag_(iConfig.getParameter<InputTag>("trackTag")),
  HBHERecHitTag_(iConfig.getParameter<InputTag>("HBHERecHitTag")),
  cosmicTrackTag_(iConfig.getParameter<InputTag>("cosmicTrackTag")),
  EBRecHitTag_(iConfig.getParameter<InputTag>("EBRecHitTag")),
  EERecHitTag_(iConfig.getParameter<InputTag>("EERecHitTag"))

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
	  ECALIso_->Fill(photon->ecalRecHitSumEtConeDR04());
	  HCALIso_->Fill(photon->hcalTowerSumEtConeDR04());
	  ET_->Fill(photon->et());
	  HOverE_->Fill(photon->hadronicOverEm());
	  fiducialRegionDist_->Fill(evtProperties_.ECALFiducialRegion(photon));
	  etaWidth_->Fill(photon->sigmaEtaEta());
	  trackIso_->Fill(photon->trkSumPtHollowConeDR04());
	  eta_->Fill(photon->eta());
	  phi_->Fill(photon->phi());
	}
      }//end loop over candidates
    }//end if found photon candidates
  }//end else

  //fill histograms
  numCandsPerEvt_->Fill(numCands);
  //TODO: fill numTracksPerEvt_
}

// ------------ method called once each job just before starting event loop  ------------
void 
CutAnalyzer::beginJob()
{
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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CutAnalyzer::endJob() {

  //write histograms to file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
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
    out_->Close();
  }
  else if (debugFlag_) {
    debug_ << "Error opening file " << outFileName_ << ".\n";
    debug_.close();
  }
  delete out_;
}

//define this as a plug-in
DEFINE_FWK_MODULE(CutAnalyzer);
