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
#include "GMSBTools/SampleMaker/interface/SampleMaker.h"
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
  const unsigned int ECALFiducialRegion(const Photon*) const;
  const string sampleTypeString() const;

      // ----------member data ---------------------------

  //output
  string outFileName_;
  TFile* out_;

  //cut values
  unsigned int sampleType_;
  double ECALIsoMaxPTMultiplierEB_;
  double ECALIsoMaxConstantEB_;
  double ECALIsoMaxPTMultiplierEE_;
  double ECALIsoMaxConstantEE_;
  double HCALIsoMaxEB_;
  double HCALIsoMaxEE_;
  double HOverEMaxPresel_;
  double ETMin_;
  unsigned int fiducialRegion_;
  double HOverEMax_;
  double sigmaEtaEtaMax_;
  double trackIsoMax_;
  double trackPTMin_;
  double eTrackRMin_;
  InputTag photonTag_;
  InputTag trackTag_;
  string debugFileName_;
  ofstream debug_;
  bool debugFlag_;

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
  TH1F* numPhotonsPerEvt_;
  TH1F* numElectronsPerEvt_;
  TH1F* numFakesPerEvt_;
  TH1F* numTracksPerEvt_;
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
  sampleType_(iConfig.getParameter<unsigned int>("sampleType")),
  ECALIsoMaxPTMultiplierEB_(iConfig.getParameter<double>("ECALIsoMaxPTMultiplierEB")),
  ECALIsoMaxConstantEB_(iConfig.getParameter<double>("ECALIsoMaxConstantEB")),
  ECALIsoMaxPTMultiplierEE_(iConfig.getParameter<double>("ECALIsoMaxPTMultiplierEE")),
  ECALIsoMaxConstantEE_(iConfig.getParameter<double>("ECALIsoMaxConstantEE")),
  HCALIsoMaxEB_(iConfig.getParameter<double>("HCALIsoMaxEB")),
  HCALIsoMaxEE_(iConfig.getParameter<double>("HCALIsoMaxEE")),
  HOverEMaxPresel_(iConfig.getParameter<double>("HOverEMaxPresel")),
  ETMin_(iConfig.getParameter<double>("ETMin")),
  fiducialRegion_(iConfig.getParameter<unsigned int>("fiducialRegion")),
  HOverEMax_(iConfig.getParameter<double>("HOverEMax")),
  sigmaEtaEtaMax_(iConfig.getParameter<double>("sigmaEtaEtaMax")),
  trackIsoMax_(iConfig.getParameter<double>("trackIsoMax")),
  trackPTMin_(iConfig.getUntrackedParameter<double>("trackPTMin", -1.0)),
  eTrackRMin_(iConfig.getUntrackedParameter<double>("eTrackRMin", -1.0)),
  photonTag_(iConfig.getParameter<InputTag>("photonTag")),
  trackTag_(iConfig.getParameter<InputTag>("trackTag")),
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
    debug_ << "Sample type: " << sampleType_ << "(" << sampleTypeString() << ")\n";
    debug_ << "Maximum ECAL isolation in EB: " << ECALIsoMaxPTMultiplierEB_ << "*pT + " << ECALIsoMaxConstantEB_ << " GeV\n";
    debug_ << "Maximum ECAL isolation in EE: " << ECALIsoMaxPTMultiplierEE_ << "*pT + " << ECALIsoMaxConstantEE_ << " GeV\n";
    debug_ << "Maximum HCAL isolation in EB: " << HCALIsoMaxEB_ << " GeV\n";
    debug_ << "Maximum HCAL isolation in EE: " << HCALIsoMaxEE_ << " GeV\n";
    debug_ << "Minimum photon ET: " << ETMin_ << " GeV\n";
    debug_ << "Fiducial region: " << fiducialRegion_ << "(";
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
    debug_ << "Maximum H/E: " << HOverEMax_ << endl;
    debug_ << "Maximum photon supercluster eta width: " << sigmaEtaEtaMax_ << endl;
    debug_ << "Maximum photon track isolation: " << trackIsoMax_ << endl;
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


CutAnalyzer::~CutAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete out_;
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
    unsigned int numPhotons = 0;
    unsigned int numElectrons = 0;
    unsigned int numFakes = 0;
    vector<Photon*> electrons;
    unsigned int numPhotonsProcessed = 0;
    for (PhotonCollection::const_iterator iPhoton = pPhotons->begin(); iPhoton != pPhotons->end(); ++iPhoton) {
      bool foundPassingCand = false;

      /*all samples require at least one object satisfying:
	ECAL isolation < ECALIsoMaxPTMultiplierEB_*pT + ECALIsoMaxConstantEB_ GeV (EB), 
	ECALIsoMaxPTMultiplierEE_*pT + ECALIsoMaxConstantEE_ GeV (EE)
	HCAL isolation < HCALIsoMaxEB_ GeV (EB), HCALIsoMaxEE_ GeV (EE)
	H/E < HOverEMaxPresel_
	ET > ETMin_
	supercluster in fiducial region
      */
      const double pT = iPhoton->pt();
      const double ECALIso = (double)iPhoton->ecalRecHitSumEtConeDR04();
      double ECALIsoMax = ECALIsoMaxPTMultiplierEB_*pT + ECALIsoMaxConstantEB_;
      const double HCALIso = (double)iPhoton->hcalTowerSumEtConeDR04();
      double HCALIsoMax = HCALIsoMaxEB_;
      const unsigned int fiducialRegion = ECALFiducialRegion(const_cast<const Photon*>(&*iPhoton));
      if (fiducialRegion == EEND) {
	ECALIsoMax = ECALIsoMaxPTMultiplierEE_*pT + ECALIsoMaxConstantEE_;
	HCALIsoMax = HCALIsoMaxEE_;
      }
      const double HOverE = (double)iPhoton->hadronicOverEm();
      const double ET = iPhoton->et();
      const double sigmaEtaEta = (double)iPhoton->sigmaEtaEta();
      const double trackIso = (double)iPhoton->trkSumPtHollowConeDR04();
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

	/*gammagamma, egamma, ee, and etrack samples require:
	 sigmaEtaEta < sigmaEtaEtaMax_
	 H/E < HOverEMax_
	 track isolation < trackIsoMax_
	*/
	if (((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA) || (sampleType_ == EE) || (sampleType_ == ETRACK)) && 
	    ((sigmaEtaEta < sigmaEtaEtaMax_) && (HOverE < HOverEMax_) && (trackIso < trackIsoMax_))) {
	  if (debugFlag_) {
	    debug_ << "Photon supercluster eta width: " << sigmaEtaEta << endl;
	    debug_ << "Photon track isolation: " << trackIso << endl;
	  }

	  //egamma, ee, and etrack samples require the candidate to have at least one pixel seed
	  if (iPhoton->hasPixelSeed()) {

	    /*ee and etrack samples require only electrons
	      for the etrack sample, we also need to save the candidate*/
	    if ((sampleType_ == EE) || (sampleType_ == ETRACK) || (sampleType_ == EGAMMA)) {
	      ++numElectrons;
	      if (debugFlag_) debug_ << "Number of found electrons: " << numElectrons << endl;
	      electrons.push_back(const_cast<Photon*>(&*iPhoton));
	      foundPassingCand = true;
	    }
	  }

	  //gammagamma and egamma samples require no pixel seed
	  else {

	    //gammagamma sample requires only photons
	    if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA)) {
	      ++numPhotons;
	      if (debugFlag_) debug_ << "Number of found photons: " << numPhotons << endl;
	      foundPassingCand = true;
	    }
	  }
	}

	/*ff sample requires:
	 sigmaEtaEta > sigmaEtaEtaMax_ or
	 H/E > HOverEMax_ or
	 track isolation > trackIsoMax_

	 NB. No requirement on number of pixel seeds in the fake definition!
	*/
	else if ((sampleType_ == FF) && ((sigmaEtaEta >= sigmaEtaEtaMax_) || (HOverE >=  HOverEMax_) || (trackIso >= trackIsoMax_))) {
	  ++numFakes;
	  foundPassingCand = true;
	  if (debugFlag_) {
	    debug_ << "Number of found fakes: " << numFakes << endl;
	    debug_ << "Photon supercluster eta width: " << sigmaEtaEta << endl;
	    debug_ << "Photon track isolation: " << trackIso << endl;
	  }
	}
      }

      //fill histograms
      if (foundPassingCand) {
	ECALIso_->Fill(ECALIso);
	HCALIso_->Fill(HCALIso);
	ET_->Fill(ET);
	HOverE_->Fill(HOverE);
	fiducialRegionDist_->Fill(fiducialRegion);
	etaWidth_->Fill(sigmaEtaEta);
	trackIso_->Fill(trackIso);
      }

      //advance to the next photon in the collection
      ++numPhotonsProcessed;
      if (debugFlag_) debug_ << endl;
    }
    if (debugFlag_) debug_ << "Number of photons processed: " << numPhotonsProcessed << "/" << pPhotons->size() << endl << endl;

    //fill histograms
    if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA)) numPhotonsPerEvt_->Fill(numPhotons);
    if ((sampleType_ == EGAMMA) || (sampleType_ == EE) || (sampleType_ == ETRACK)) numElectronsPerEvt_->Fill(numElectrons);
    if (sampleType_ == FF) numFakesPerEvt_->Fill(numFakes);

    //etrack sample requires a track
    if (sampleType_ == ETRACK) {

      //get the GeneralTrack collection
      Handle<TrackCollection> pTracks;
      bool tracksFound = false;
      try { tracksFound = (iEvent.getByLabel(trackTag_, pTracks)) && (pTracks->size() > 0); }
      catch (cms::Exception& ex) {}
      if (!tracksFound) {
	if (debugFlag_) debug_ << "No track collection found in run " << iEvent.run() << ", event " << iEvent.id().event() << ".\n";
      }
      else {

	/*etrack sample requires:
	  track pT > trackPTMin_
	  dR(track, electron) > eTrackRMin_
	*/
	unsigned int numTracks = 0;
	for (TrackCollection::const_iterator iTrack = pTracks->begin(); iTrack != pTracks->end(); ++iTrack) {
	  double pT = iTrack->pt();
	  if (pT > trackPTMin_) {
	    for (vector<Photon*>::const_iterator iElectron = electrons.begin(); iElectron != electrons.end(); ++iElectron) {
	      double dEtaETrack = (*iElectron)->eta() - iTrack->eta();
	      double dPhiETrack = (*iElectron)->phi() - iTrack->phi();
	      double dRETrack = sqrt((dEtaETrack)*(dEtaETrack) + (dPhiETrack)*(dPhiETrack));
	      if (dRETrack > eTrackRMin_) {
		++numTracks;
		if (debugFlag_) {
		  debug_ << "Track pT: " << pT << " GeV\n";
		  debug_ << "dR(track, photon): " << dRETrack << endl;
		}
	      }
	      else if (debugFlag_) debug_ << "This track/electron pair doesn't meet the criteria; advancing to the next electron.\n";
	    }
	  }
	  else if (debugFlag_) debug_ << "This track doesn't meet the pT criteria; advancing to the next track.\n";
	}

	//fill histogram
	numTracksPerEvt_->Fill(numTracks);
      }
    }
  }
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
  if (sampleType_ == ETRACK) {
    trackPT_ = new TH1D("trackPT_", "p_{T} of passing tracks;p_{T} (GeV);Tracks per GeV", 100, 0.0, 100.0);
    dRTrackPhoton_ = new TH1D("dRTrackPhoton_", "#DeltaR(track, e/#gamma candidate) of passing tracks-e/#gamma candidate pairs;#DeltaR;", 
			      100, 0.0, 1.0);
    dRTrackPhoton_->GetYaxis()->SetTitle("track-e/#gamma candidate pairs per 0.01");
    numTracksPerEvt_ = new TH1F("numTracksPerEvt_", "Number of passing tracks per event;Number per event;Events per 1", 5, -0.5, 4.5);
  }
  if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA)) {
    numPhotonsPerEvt_ = new TH1F("numPhotonsPerEvt_", 
				 "Number of passing #gamma candidates per event;Number per event;Events per 1", 5, -0.5, 4.5);
  }
  if ((sampleType_ == EE) || (sampleType_ == ETRACK) || (sampleType_ == EGAMMA)) {
    numElectronsPerEvt_ = new TH1F("numElectronsPerEvt_", "Number of passing e candidates per event;Number per event;Events per 1", 5, 
				   -0.5, 4.5);
  }
  if (sampleType_ == FF) {
    numFakesPerEvt_ = new TH1F("numFakesPerEvt_", "Number of passing fake candidates per event;Number per event;Events per 1", 5, -0.5, 
			       4.5);
  }
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
    if (sampleType_ == ETRACK) {
      trackPT_->Write();
      dRTrackPhoton_->Write();
      numTracksPerEvt_->Write();
    }
    if ((sampleType_ == GAMMAGAMMA) || (sampleType_ == EGAMMA)) numPhotonsPerEvt_->Write();
    if ((sampleType_ == EE) || (sampleType_ == ETRACK) || (sampleType_ == EGAMMA)) numElectronsPerEvt_->Write();
    if (sampleType_ == FF) numFakesPerEvt_->Write();
  }
  else if (debugFlag_) {
    debug_ << "Error opening file " << outFileName_ << ".\n";
    debug_.close();
  }
}

//get the ECAL fiducial region in units this code understands
const unsigned int CutAnalyzer::ECALFiducialRegion(const Photon* photon) const
{
  unsigned int reg = EB;
  if (photon->isEE()) reg = EEND;
  return reg;
}

//get the sample type as a string
const string CutAnalyzer::sampleTypeString() const
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

//define this as a plug-in
DEFINE_FWK_MODULE(CutAnalyzer);
