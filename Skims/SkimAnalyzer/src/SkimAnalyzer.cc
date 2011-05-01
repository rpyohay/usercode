//normalize to no cuts?  prescaled trigger with no iso?

//put workflow in master script?  elog?

//Need to:
//  - will want to do this for the Dongwook format too -- can we make it portable?
//  - make it portable by having a class that takes primitive data types to do the classifying
//  - then have a converter to/from edm formats, converter to/from Dongwook formats from primitive type
//  - have an EDMCategoryProducer and DongwookCategoryProducer
//- later, combine cuts and repeat previous step

//get rid of CutScanKey
//make scan points a parameter set to pass to CutScanProducer and SkimAnalyzer

// -*- C++ -*-
//
// Package:    Skims
// Class:      SkimAnalyzer
// 
/*
 Description: analysis of data volume reduction with different skim definitions

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay
//         Created:  Wed Mar 23 17:16:33 CET 2011
// $Id: SkimAnalyzer.cc,v 1.1 2011/03/23 17:16:33 yohay Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "GMSBTools/Filters/interface/Typedefs.h"
#include "GMSBTools/Filters/interface/Categorizer.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TFile.h"
#include "TH2F.h"


//
// class declaration
//

class SkimAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SkimAnalyzer(const edm::ParameterSet&);
      ~SkimAnalyzer();

   private:
      virtual void beginRun(const edm::Run&, edm::EventSetup const&);
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  //return a properly formatted histogram name
  STRING name(const STRING&, const STRING&, const STRING&, const unsigned int, 
	      const STRING&) const;

  //book any scan point histogram
  void bookScanPointHistograms(VTH1F&, const VDOUBLE&, const STRING&, const STRING&, 
			       const STRING&, const STRING&, const unsigned int, const float, 
			       const float);

  //book all the gg MET histograms
  void bookGGMETHistograms(const STRING&);

  //book all the ff MET histograms
  void bookFFMETHistograms(const STRING&);

  //book all the eg MET histograms
  void bookEGMETHistograms(const STRING&);

  //book all the ee MET histograms
  void bookEEMETHistograms(const STRING&);

  //book all the gg diEMET histograms
  void bookGGDiEMETHistograms(const STRING&);

  //book all the ff diEMET histograms
  void bookFFDiEMETHistograms(const STRING&);

  //book all the eg diEMET histograms
  void bookEGDiEMETHistograms(const STRING&);

  //book all the ee diEMET histograms
  void bookEEDiEMETHistograms(const STRING&);

  //book all ET N-1 histograms
  void bookETNMinus1Histograms(const STRING&);

  //book all eta N-1 histograms
  void bookEtaNMinus1Histograms(const STRING&);

  //book all phi N-1 histograms
  void bookPhiNMinus1Histograms(const STRING&);

  //book all ECAL isolation N-1 histograms
  void bookECALIsoNMinus1Histograms(const STRING&);

  //book all HCAL isolation N-1 histograms
  void bookHCALIsoNMinus1Histograms(const STRING&);

  //book all H/E N-1 histograms
  void bookHOverENMinus1Histograms(const STRING&);

  //book all track isolation N-1 histograms
  void bookTrackIsoNMinus1Histograms(const STRING&);

  //book all sigmaIetaIeta N-1 histograms
  void bookSigmaIetaIetaNMinus1Histograms(const STRING&);

  /*increment a vector of counters representing the number of photons passing a particular scan 
    point, 1 vector for each scan cut type (i.e. 1 for ECAL isolation, 1 for H/E, etc.)*/
  void incrementCounters(const VDOUBLE&, VUINT&, const unsigned int) const;

  //extract the scan point as an int from the histogram name
  int scanPoint(const STRING&, const STRING&) const;

  /*loop over the vector of histograms of a single physical quantity corresponding to the 
    given scan cut type -- each histogram corresponds to a particular scan point and trigger 
    combination -- and fill the histograms properly*/
  void fillIndHistogram(VTH1F*, const float, const bool, const STRING&, 
			const STRING&, const unsigned int);

  //fill histograms of event and photon quantities properly
  void fillHistograms(const VUINT&, std::vector<VTH1F*>&, const VFLOAT&, 
		      const VBOOL&, const STRING&, const STRING&);

  //return true if the photon should go into the N-1 plot, false otherwise
  bool passNMinus1(std::vector<edm::Handle<edm::ValueMap<bool> >*>&, 
		   edm::Ref<reco::PhotonCollection>&, edm::Handle<edm::ValueMap<bool> >*) const;

  //extract just the file name from the full path name
  STRING extractFileFromFullPath(const STRING&) const;

  //make ROOT directories in the output file corresponding to the different scan cut types
  void makeOutputDirectories(const unsigned int, const STRING&);

  //fill ROOT directories in the output file with corresponding histograms
  void fillOutputDirectories(VTH1F&, const STRING&);

  //retrieve collection from the event
  template<typename T>
  const bool getCollection_(T& pCollection, const edm::InputTag& tag, const edm::Event& iEvent)
  {
    bool collectionFound = false;
    try { collectionFound = iEvent.getByLabel(tag, pCollection); }
    catch (cms::Exception& ex) {}
    if (!collectionFound) {
      STRINGSTREAM err;
      err << "No collection of type " << tag << " found in run " << iEvent.run();
      err << ", event " << iEvent.id().event() << ", lumi section ";
      err << iEvent.getLuminosityBlock().luminosityBlock() << ".\n";
      edm::LogInfo("Error") << err.str();
    }
    return collectionFound;
  }
      
      // ----------member data ---------------------------

  //input
  STRING outputFile_;
  edm::InputTag photonTag_;
  edm::InputTag categoryTag_;
  edm::InputTag diEMETTag_;
  edm::InputTag trgResultsTag_;
  edm::InputTag tcMETTag_;
  edm::InputTag ETMinScanDecisionTag_;
  edm::InputTag ECALIsoMaxScanDecisionTag_;
  edm::InputTag HCALIsoMaxScanDecisionTag_;
  edm::InputTag HOverEMaxScanDecisionTag_;
  edm::InputTag trackIsoMaxScanDecisionTag_;
  edm::InputTag sigmaIetaIetaMaxScanDecisionTag_;
  edm::InputTag passETMinTag_;
  edm::InputTag passAbsEtaMaxTag_;
  edm::InputTag passECALIsoMaxTag_;
  edm::InputTag passHCALIsoMaxTag_;
  edm::InputTag passHOverEMaxTag_;
  edm::InputTag passTrackIsoMaxTag_;
  edm::InputTag passSigmaIetaIetaMaxTag_;
  /*edm::InputTag passAbsSeedTimeMaxTag_;
    edm::InputTag passE2OverE9MaxTag_;*/

  //scan points
  VDOUBLE ETMinScan_;
  VDOUBLE ECALIsoMaxPTMultiplierScan_;
  VDOUBLE ECALIsoMaxConstantScan_;
  VDOUBLE HCALIsoMaxPTMultiplierScan_;
  VDOUBLE HCALIsoMaxConstantScan_;
  VDOUBLE HOverEMaxScan_;
  VDOUBLE trackIsoMaxPTMultiplierScan_;
  VDOUBLE trackIsoMaxConstantScan_;
  VDOUBLE sigmaIetaIetaMaxScan_;

  //triggers
  STRING HLTProcessName_;
  HLTConfigProvider HLTCfg_;
  VSTRING HLTPathGG_;
  VSTRING HLTPathFF_;

  //MET plots
  VTH1F ggMETETMin_;
  VTH1F ggMETECALIsoMax_;
  VTH1F ggMETHCALIsoMax_;
  VTH1F ggMETHOverEMax_;
  VTH1F ggMETTrackIsoMax_;
  VTH1F ggMETSigmaIetaIetaMax_;
  VTH1F ffMETETMin_;
  VTH1F ffMETECALIsoMax_;
  VTH1F ffMETHCALIsoMax_;
  VTH1F ffMETHOverEMax_;
  VTH1F ffMETTrackIsoMax_;
  VTH1F ffMETSigmaIetaIetaMax_;
  VTH1F eeMETETMin_;
  VTH1F eeMETECALIsoMax_;
  VTH1F eeMETHCALIsoMax_;
  VTH1F eeMETHOverEMax_;
  VTH1F eeMETTrackIsoMax_;
  VTH1F eeMETSigmaIetaIetaMax_;
  VTH1F egMETETMin_;
  VTH1F egMETECALIsoMax_;
  VTH1F egMETHCALIsoMax_;
  VTH1F egMETHOverEMax_;
  VTH1F egMETTrackIsoMax_;
  VTH1F egMETSigmaIetaIetaMax_;

  //di-EM ET plots
  VTH1F ggDiEMETETMin_;
  VTH1F ggDiEMETECALIsoMax_;
  VTH1F ggDiEMETHCALIsoMax_;
  VTH1F ggDiEMETHOverEMax_;
  VTH1F ggDiEMETTrackIsoMax_;
  VTH1F ggDiEMETSigmaIetaIetaMax_;
  VTH1F ffDiEMETETMin_;
  VTH1F ffDiEMETECALIsoMax_;
  VTH1F ffDiEMETHCALIsoMax_;
  VTH1F ffDiEMETHOverEMax_;
  VTH1F ffDiEMETTrackIsoMax_;
  VTH1F ffDiEMETSigmaIetaIetaMax_;
  VTH1F eeDiEMETETMin_;
  VTH1F eeDiEMETECALIsoMax_;
  VTH1F eeDiEMETHCALIsoMax_;
  VTH1F eeDiEMETHOverEMax_;
  VTH1F eeDiEMETTrackIsoMax_;
  VTH1F eeDiEMETSigmaIetaIetaMax_;
  VTH1F egDiEMETETMin_;
  VTH1F egDiEMETECALIsoMax_;
  VTH1F egDiEMETHCALIsoMax_;
  VTH1F egDiEMETHOverEMax_;
  VTH1F egDiEMETTrackIsoMax_;
  VTH1F egDiEMETSigmaIetaIetaMax_;

  /*1 N-1 plot per scan point and trigger combination (don't bin by event category because that 
    introduces unwanted cuts)*/

  //ET N-1 plots
  VTH1F ETNMinus1ETMin_;
  VTH1F ETNMinus1ECALIsoMax_;
  VTH1F ETNMinus1HCALIsoMax_;
  VTH1F ETNMinus1HOverEMax_;
  VTH1F ETNMinus1TrackIsoMax_;
  VTH1F ETNMinus1SigmaIetaIetaMax_;

  //eta N-1 plots
  VTH1F etaNMinus1ETMin_;
  VTH1F etaNMinus1ECALIsoMax_;
  VTH1F etaNMinus1HCALIsoMax_;
  VTH1F etaNMinus1HOverEMax_;
  VTH1F etaNMinus1TrackIsoMax_;
  VTH1F etaNMinus1SigmaIetaIetaMax_;

  //phi N-1 plots
  VTH1F phiNMinus1ETMin_;
  VTH1F phiNMinus1ECALIsoMax_;
  VTH1F phiNMinus1HCALIsoMax_;
  VTH1F phiNMinus1HOverEMax_;
  VTH1F phiNMinus1TrackIsoMax_;
  VTH1F phiNMinus1SigmaIetaIetaMax_;

  //ECAL isolation N-1 plots
  VTH1F ECALIsoNMinus1ETMin_;
  VTH1F ECALIsoNMinus1ECALIsoMax_;
  VTH1F ECALIsoNMinus1HCALIsoMax_;
  VTH1F ECALIsoNMinus1HOverEMax_;
  VTH1F ECALIsoNMinus1TrackIsoMax_;
  VTH1F ECALIsoNMinus1SigmaIetaIetaMax_;

  //HCAL isolation N-1 plots
  VTH1F HCALIsoNMinus1ETMin_;
  VTH1F HCALIsoNMinus1ECALIsoMax_;
  VTH1F HCALIsoNMinus1HCALIsoMax_;
  VTH1F HCALIsoNMinus1HOverEMax_;
  VTH1F HCALIsoNMinus1TrackIsoMax_;
  VTH1F HCALIsoNMinus1SigmaIetaIetaMax_;

  //H/E N-1 plots
  VTH1F HOverENMinus1ETMin_;
  VTH1F HOverENMinus1ECALIsoMax_;
  VTH1F HOverENMinus1HCALIsoMax_;
  VTH1F HOverENMinus1HOverEMax_;
  VTH1F HOverENMinus1TrackIsoMax_;
  VTH1F HOverENMinus1SigmaIetaIetaMax_;

  //track isolation N-1 plots
  VTH1F trackIsoNMinus1ETMin_;
  VTH1F trackIsoNMinus1ECALIsoMax_;
  VTH1F trackIsoNMinus1HCALIsoMax_;
  VTH1F trackIsoNMinus1HOverEMax_;
  VTH1F trackIsoNMinus1TrackIsoMax_;
  VTH1F trackIsoNMinus1SigmaIetaIetaMax_;

  //sigmaIetaIeta N-1 plots
  VTH1F sigmaIetaIetaNMinus1ETMin_;
  VTH1F sigmaIetaIetaNMinus1ECALIsoMax_;
  VTH1F sigmaIetaIetaNMinus1HCALIsoMax_;
  VTH1F sigmaIetaIetaNMinus1HOverEMax_;
  VTH1F sigmaIetaIetaNMinus1TrackIsoMax_;
  VTH1F sigmaIetaIetaNMinus1SigmaIetaIetaMax_;

  /*//seedTime N-1 plots
  VTH1F seedTimeNMinus1ETMin_;
  VTH1F seedTimeNMinus1ECALIsoMax_;
  VTH1F seedTimeNMinus1HCALIsoMax_;
  VTH1F seedTimeNMinus1HOverEMax_;
  VTH1F seedTimeNMinus1TrackIsoMax_;
  VTH1F seedTimeNMinus1SigmaIetaIetaMax_;

  //sigmaIetaIeta N-1 plots
  VTH1F e2OverE9NMinus1ETMin_;
  VTH1F e2OverE9NMinus1ECALIsoMax_;
  VTH1F e2OverE9NMinus1HCALIsoMax_;
  VTH1F e2OverE9NMinus1HOverEMax_;
  VTH1F e2OverE9NMinus1TrackIsoMax_;
  VTH1F e2OverE9NMinus1SigmaIetaIetaMax_;*/

  //hack
  TH2F* HOverEVsAbsEta_;

  //output
  TFile* out_;
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
SkimAnalyzer::SkimAnalyzer(const edm::ParameterSet& iConfig) :

  //input
  outputFile_(iConfig.getUntrackedParameter<STRING>("outputFile", "skim_analysis.root")),
  photonTag_(iConfig.getUntrackedParameter<edm::InputTag>
	     ("photonTag", edm::InputTag("photons", "", "RECOCleaned"))),
  categoryTag_(iConfig.getUntrackedParameter<edm::InputTag>
	       ("categoryTag", edm::InputTag("EDMCategoryProducer", "eventCategory", 
					     "CATEGORY"))),
  diEMETTag_(iConfig.getUntrackedParameter<edm::InputTag>
	     ("diEMETTag", edm::InputTag("EDMCategoryProducer", "evtDiEMET", "CATEGORY"))),
  trgResultsTag_(iConfig.getUntrackedParameter<edm::InputTag>
		 ("trgResultsTag", edm::InputTag("TriggerResults", "", "HLT"))),
  tcMETTag_(iConfig.getUntrackedParameter<edm::InputTag>
	    ("tcMETTag", edm::InputTag("tcMet", "", "RECO"))),
  ETMinScanDecisionTag_(iConfig.getUntrackedParameter<edm::InputTag>
			("ETMinScanDecisionTag", 
			 edm::InputTag("CutScanProducer", "passETMin", "OWNPARTICLES"))),
  ECALIsoMaxScanDecisionTag_(iConfig.getUntrackedParameter<edm::InputTag>
			     ("ECALIsoMaxScanDecisionTag", 
			      edm::InputTag("CutScanProducer", "passECALIsoMax", "OWNPARTICLES"))),
  HCALIsoMaxScanDecisionTag_(iConfig.getUntrackedParameter<edm::InputTag>
			     ("HCALIsoMaxScanDecisionTag", 
			      edm::InputTag("CutScanProducer", "passHCALIsoMax", "OWNPARTICLES"))),
  HOverEMaxScanDecisionTag_(iConfig.getUntrackedParameter<edm::InputTag>
			    ("HOverEMaxScanDecisionTag", 
			     edm::InputTag("CutScanProducer", "passHOverEMax", "OWNPARTICLES"))),
  trackIsoMaxScanDecisionTag_(iConfig.getUntrackedParameter<edm::InputTag>
			      ("trackIsoMaxScanDecisionTag", 
			       edm::InputTag("CutScanProducer", "passTrackIsoMax", 
					     "OWNPARTICLES"))),
  sigmaIetaIetaMaxScanDecisionTag_(iConfig.getUntrackedParameter<edm::InputTag>
				   ("sigmaIetaIetaMaxScanDecisionTag", 
				    edm::InputTag("CutScanProducer", "passSigmaIetaIetaMax", 
						  "OWNPARTICLES"))),
  passETMinTag_(iConfig.getUntrackedParameter<edm::InputTag>
		("passETMinTag", edm::InputTag("EDMCategoryProducer", 
					       "passETMin", "OWNPARTICLES"))),
  passAbsEtaMaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
		    ("passAbsEtaMaxTag", edm::InputTag("EDMCategoryProducer", 
						       "passAbsEtaMax", "OWNPARTICLES"))),
  passECALIsoMaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
		     ("passECALIsoMaxTag", edm::InputTag("EDMCategoryProducer", 
							 "passECALIsoMax", "OWNPARTICLES"))),
  passHCALIsoMaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
		     ("passHCALIsoMaxTag", edm::InputTag("EDMCategoryProducer", 
							 "passHCALIsoMax", "OWNPARTICLES"))),
  passHOverEMaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
		    ("passHOverEMaxTag", edm::InputTag("EDMCategoryProducer", 
						       "passHOverEMax", "OWNPARTICLES"))),
  passTrackIsoMaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
		      ("passTrackIsoMaxTag", edm::InputTag("EDMCategoryProducer", 
							   "passTrackIsoMax", "OWNPARTICLES"))),
  passSigmaIetaIetaMaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
			   ("passSigmaIetaIetaMaxTag", edm::InputTag("EDMCategoryProducer", 
								     "passSigmaIetaIetaMax", 
								     "OWNPARTICLES"))),
  /*passAbsSeedTimeMaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
			 ("passAbsSeedTimeMaxTag", edm::InputTag("EDMCategoryProducer", 
								 "passAbsSeedTimeMax", 
								 "OWNPARTICLES"))),
  passE2OverE9MaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
		      ("passE2OverE9MaxTag", edm::InputTag("EDMCategoryProducer", 
		      "passE2OverE9Max", "OWNPARTICLES"))),*/

  //scan points
  ETMinScan_(iConfig.getParameter<VDOUBLE>("ETMinScan")),
  ECALIsoMaxPTMultiplierScan_(iConfig.getParameter<VDOUBLE>("ECALIsoMaxPTMultiplierScan")),
  ECALIsoMaxConstantScan_(iConfig.getParameter<VDOUBLE>("ECALIsoMaxConstantScan")),
  HCALIsoMaxPTMultiplierScan_(iConfig.getParameter<VDOUBLE>("HCALIsoMaxPTMultiplierScan")),
  HCALIsoMaxConstantScan_(iConfig.getParameter<VDOUBLE>("HCALIsoMaxConstantScan")),
  HOverEMaxScan_(iConfig.getParameter<VDOUBLE>("HOverEMaxScan")),
  trackIsoMaxPTMultiplierScan_(iConfig.getParameter<VDOUBLE>("trackIsoMaxPTMultiplierScan")),
  trackIsoMaxConstantScan_(iConfig.getParameter<VDOUBLE>("trackIsoMaxConstantScan")),
  sigmaIetaIetaMaxScan_(iConfig.getParameter<VDOUBLE>("sigmaIetaIetaMaxScan")),

  //triggers
  HLTProcessName_(iConfig.getUntrackedParameter<STRING>("HLTProcessName", "HLT")),
  HLTPathGG_(iConfig.getUntrackedParameter<VSTRING>
	     ("HLTPathGG", VSTRING(1, "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1"))),
  HLTPathFF_(iConfig.getUntrackedParameter<VSTRING>
	     ("HLTPathFF", VSTRING(1, "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1")))
{
}

SkimAnalyzer::~SkimAnalyzer()
{}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
SkimAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get the event category, trigger information, tcMET, cut scan results, and photon collection
  edm::Handle<reco::PhotonCollection> pPhotons;
  edm::Handle<int> pCategory;
  edm::Handle<double> pDiEMET;
  edm::Handle<edm::TriggerResults> pTrgResults;
  edm::Handle<reco::METCollection> pTCMET;
  edm::Handle<edm::ValueMap<bool> > pPassETMin;
  edm::Handle<edm::ValueMap<bool> > pPassAbsEtaMax;
  edm::Handle<edm::ValueMap<bool> > pPassECALIsoMax;
  edm::Handle<edm::ValueMap<bool> > pPassHCALIsoMax;
  edm::Handle<edm::ValueMap<bool> > pPassHOverEMax;
  edm::Handle<edm::ValueMap<bool> > pPassTrackIsoMax;
  edm::Handle<edm::ValueMap<bool> > pPassSigmaIetaIetaMax;
  /*edm::Handle<edm::ValueMap<bool> > pPassAbsSeedTimeMax;
    edm::Handle<edm::ValueMap<bool> > pPassE2OverE9Max;*/
  if (getCollection_(pPhotons, photonTag_, iEvent) && 
      getCollection_(pCategory, categoryTag_, iEvent) && 
      getCollection_(pDiEMET, diEMETTag_, iEvent) && 
      getCollection_(pTrgResults, trgResultsTag_, iEvent) && 
      getCollection_(pTCMET, tcMETTag_, iEvent) &&  
      getCollection_(pPassETMin, passETMinTag_, iEvent) && 
      getCollection_(pPassAbsEtaMax, passAbsEtaMaxTag_, iEvent) && 
      getCollection_(pPassECALIsoMax, passECALIsoMaxTag_, iEvent) && 
      getCollection_(pPassHCALIsoMax, passHCALIsoMaxTag_, iEvent) && 
      getCollection_(pPassHOverEMax, passHOverEMaxTag_, iEvent) && 
      getCollection_(pPassTrackIsoMax, passTrackIsoMaxTag_, iEvent) && 
      getCollection_(pPassSigmaIetaIetaMax, passSigmaIetaIetaMaxTag_, iEvent)/* && 
      getCollection_(pPassAbsSeedTimeMax, passAbsSeedTimeMaxTag_, iEvent) && 
      getCollection_(pPassE2OverE9Max, passE2OverE9MaxTag_, iEvent)*/) {

    /*for each scan point, loop over value maps to determine how many photons passed the scan 
      point criteria*/
    edm::Handle<edm::ValueMap<unsigned int> > pETMinScanDecision;
    edm::Handle<edm::ValueMap<unsigned int> > pECALIsoMaxScanDecision;
    edm::Handle<edm::ValueMap<unsigned int> > pHCALIsoMaxScanDecision;
    edm::Handle<edm::ValueMap<unsigned int> > pHOverEMaxScanDecision;
    edm::Handle<edm::ValueMap<unsigned int> > pTrackIsoMaxScanDecision;
    edm::Handle<edm::ValueMap<unsigned int> > pSigmaIetaIetaMaxScanDecision;
    bool gotETMinScanDecision = getCollection_(pETMinScanDecision, ETMinScanDecisionTag_, iEvent);
    bool gotECALIsoMaxScanDecision = getCollection_(pECALIsoMaxScanDecision, 
						    ECALIsoMaxScanDecisionTag_, iEvent);
    bool gotHCALIsoMaxScanDecision = getCollection_(pHCALIsoMaxScanDecision, 
						    HCALIsoMaxScanDecisionTag_, iEvent);
    bool gotHOverEMaxScanDecision = getCollection_(pHOverEMaxScanDecision, 
						   HOverEMaxScanDecisionTag_, iEvent);
    bool gotTrackIsoMaxScanDecision = getCollection_(pTrackIsoMaxScanDecision, 
						     trackIsoMaxScanDecisionTag_, iEvent);
    bool gotSigmaIetaIetaMaxScanDecision = getCollection_(pSigmaIetaIetaMaxScanDecision, 
							  sigmaIetaIetaMaxScanDecisionTag_, 
							  iEvent);
    VUINT numPassingETMin(ETMinScan_.size(), 0);
    VUINT numPassingECALIsoMax(ECALIsoMaxPTMultiplierScan_.size(), 0);
    VUINT numPassingHCALIsoMax(HCALIsoMaxPTMultiplierScan_.size(), 0);
    VUINT numPassingHOverEMax(HOverEMaxScan_.size(), 0);
    VUINT numPassingTrackIsoMax(trackIsoMaxPTMultiplierScan_.size(), 0);
    VUINT numPassingSigmaIetaIetaMax(sigmaIetaIetaMaxScan_.size(), 0);
    for (reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin(); 
	 iPhoton != pPhotons->end(); ++iPhoton) {
      edm::Ref<reco::PhotonCollection> ref(pPhotons, iPhoton - pPhotons->begin());
      unsigned int ETMinScanDecisionWord = 0x0;
      unsigned int ECALIsoMaxScanDecisionWord = 0x0;
      unsigned int HCALIsoMaxScanDecisionWord = 0x0;
      unsigned int HOverEMaxScanDecisionWord = 0x0;
      unsigned int trackIsoMaxScanDecisionWord = 0x0;
      unsigned int sigmaIetaIetaMaxScanDecisionWord = 0x0;
      if (gotETMinScanDecision) ETMinScanDecisionWord = (*(pETMinScanDecision.product()))[ref];
      if (gotECALIsoMaxScanDecision) {
	ECALIsoMaxScanDecisionWord = (*(pECALIsoMaxScanDecision.product()))[ref];
      }
      if (gotHCALIsoMaxScanDecision) {
	HCALIsoMaxScanDecisionWord = (*(pHCALIsoMaxScanDecision.product()))[ref];
      }
      if (gotHOverEMaxScanDecision) {
	HOverEMaxScanDecisionWord = (*(pHOverEMaxScanDecision.product()))[ref];
      }
      if (gotTrackIsoMaxScanDecision) {
	trackIsoMaxScanDecisionWord = (*(pTrackIsoMaxScanDecision.product()))[ref];
      }
      if (gotSigmaIetaIetaMaxScanDecision) {
	sigmaIetaIetaMaxScanDecisionWord = (*(pSigmaIetaIetaMaxScanDecision.product()))[ref];
      }
      incrementCounters(ETMinScan_, numPassingETMin, ETMinScanDecisionWord);
      incrementCounters(ECALIsoMaxPTMultiplierScan_, numPassingECALIsoMax, 
			ECALIsoMaxScanDecisionWord);
      incrementCounters(HCALIsoMaxPTMultiplierScan_, numPassingHCALIsoMax, 
			HCALIsoMaxScanDecisionWord);
      incrementCounters(HOverEMaxScan_, numPassingHOverEMax, HOverEMaxScanDecisionWord);
      incrementCounters(trackIsoMaxPTMultiplierScan_, numPassingTrackIsoMax, 
			trackIsoMaxScanDecisionWord);
      incrementCounters(sigmaIetaIetaMaxScan_, numPassingSigmaIetaIetaMax, 
			sigmaIetaIetaMaxScanDecisionWord);
    }

    //pointers to triggers and plots relevant to this event
    VSTRING* HLTPaths = NULL;
    VTH1F* METETMin = NULL;
    VTH1F* METECALIsoMax = NULL;
    VTH1F* METHCALIsoMax = NULL;
    VTH1F* METHOverEMax = NULL;
    VTH1F* METTrackIsoMax = NULL;
    VTH1F* METSigmaIetaIetaMax = NULL;
    VTH1F* diEMETETMin = NULL;
    VTH1F* diEMETECALIsoMax = NULL;
    VTH1F* diEMETHCALIsoMax = NULL;
    VTH1F* diEMETHOverEMax = NULL;
    VTH1F* diEMETTrackIsoMax = NULL;
    VTH1F* diEMETSigmaIetaIetaMax = NULL;
    std::vector<VTH1F*> evtHistPtrVecETMin;
    std::vector<VTH1F*> evtHistPtrVecECALIsoMax;
    std::vector<VTH1F*> evtHistPtrVecHCALIsoMax;
    std::vector<VTH1F*> evtHistPtrVecHOverEMax;
    std::vector<VTH1F*> evtHistPtrVecTrackIsoMax;
    std::vector<VTH1F*> evtHistPtrVecSigmaIetaIetaMax;
    std::vector<VTH1F*> photonHistPtrVecETMin;
    std::vector<VTH1F*> photonHistPtrVecECALIsoMax;
    std::vector<VTH1F*> photonHistPtrVecHCALIsoMax;
    std::vector<VTH1F*> photonHistPtrVecHOverEMax;
    std::vector<VTH1F*> photonHistPtrVecTrackIsoMax;
    std::vector<VTH1F*> photonHistPtrVecSigmaIetaIetaMax;

    //need a vector of pointers to the value maps for the N-1 plot decision
    std::vector<edm::Handle<edm::ValueMap<bool> >*> passMapVec;
    passMapVec.push_back(&pPassETMin);
    passMapVec.push_back(&pPassAbsEtaMax);
    passMapVec.push_back(&pPassECALIsoMax);
    passMapVec.push_back(&pPassHCALIsoMax);
    passMapVec.push_back(&pPassHOverEMax);
    passMapVec.push_back(&pPassTrackIsoMax);
    passMapVec.push_back(&pPassSigmaIetaIetaMax);
    /*passMapVec.push_back(&pPassAbsSeedTimeMax);
    passMapVec.push_back(&pPassE2OverE9Max);*/

    //get the correct plots and triggers for event category
    switch (*(pCategory.product())) {
    case FAIL:
      HLTPaths = &HLTPathGG_; /*to do: address which triggers to require for N-1 plots of photons 
				in failing events*/
      break;
    case GG:
      HLTPaths = &HLTPathGG_;
      METETMin = &ggMETETMin_;
      METECALIsoMax = &ggMETECALIsoMax_;
      METHCALIsoMax = &ggMETHCALIsoMax_;
      METHOverEMax = &ggMETHOverEMax_;
      METTrackIsoMax = &ggMETTrackIsoMax_;
      METSigmaIetaIetaMax = &ggMETSigmaIetaIetaMax_;
      diEMETETMin = &ggDiEMETETMin_;
      diEMETECALIsoMax = &ggDiEMETECALIsoMax_;
      diEMETHCALIsoMax = &ggDiEMETHCALIsoMax_;
      diEMETHOverEMax = &ggDiEMETHOverEMax_;
      diEMETTrackIsoMax = &ggDiEMETTrackIsoMax_;
      diEMETSigmaIetaIetaMax = &ggDiEMETSigmaIetaIetaMax_;
      break;
    case FF:
      HLTPaths = &HLTPathFF_;
      METETMin = &ffMETETMin_;
      METECALIsoMax = &ffMETECALIsoMax_;
      METHCALIsoMax = &ffMETHCALIsoMax_;
      METHOverEMax = &ffMETHOverEMax_;
      METTrackIsoMax = &ffMETTrackIsoMax_;
      METSigmaIetaIetaMax = &ffMETSigmaIetaIetaMax_;
      diEMETETMin = &ffDiEMETETMin_;
      diEMETECALIsoMax = &ffDiEMETECALIsoMax_;
      diEMETHCALIsoMax = &ffDiEMETHCALIsoMax_;
      diEMETHOverEMax = &ffDiEMETHOverEMax_;
      diEMETTrackIsoMax = &ffDiEMETTrackIsoMax_;
      diEMETSigmaIetaIetaMax = &ffDiEMETSigmaIetaIetaMax_;
      break;
    case EG:
      HLTPaths = &HLTPathGG_;
      METETMin = &egMETETMin_;
      METECALIsoMax = &egMETECALIsoMax_;
      METHCALIsoMax = &egMETHCALIsoMax_;
      METHOverEMax = &egMETHOverEMax_;
      METTrackIsoMax = &egMETTrackIsoMax_;
      METSigmaIetaIetaMax = &egMETSigmaIetaIetaMax_;
      diEMETETMin = &egDiEMETETMin_;
      diEMETECALIsoMax = &egDiEMETECALIsoMax_;
      diEMETHCALIsoMax = &egDiEMETHCALIsoMax_;
      diEMETHOverEMax = &egDiEMETHOverEMax_;
      diEMETTrackIsoMax = &egDiEMETTrackIsoMax_;
      diEMETSigmaIetaIetaMax = &egDiEMETSigmaIetaIetaMax_;
      break;
    case EE:
      HLTPaths = &HLTPathGG_;
      METETMin = &eeMETETMin_;
      METECALIsoMax = &eeMETECALIsoMax_;
      METHCALIsoMax = &eeMETHCALIsoMax_;
      METHOverEMax = &eeMETHOverEMax_;
      METTrackIsoMax = &eeMETTrackIsoMax_;
      METSigmaIetaIetaMax = &eeMETSigmaIetaIetaMax_;
      diEMETETMin = &eeDiEMETETMin_;
      diEMETECALIsoMax = &eeDiEMETECALIsoMax_;
      diEMETHCALIsoMax = &eeDiEMETHCALIsoMax_;
      diEMETHOverEMax = &eeDiEMETHOverEMax_;
      diEMETTrackIsoMax = &eeDiEMETTrackIsoMax_;
      diEMETSigmaIetaIetaMax = &eeDiEMETSigmaIetaIetaMax_;
      break;
    default:
      break;
    }

    //get the real trigger names from the approximate inputs
    const VSTRING HLTPathNames = HLTCfg_.triggerNames();
    VSTRING realHLT;
    if (HLTPaths != NULL) {
      for (VSTRING_IT iHLT = HLTPaths->begin(); iHLT != HLTPaths->end(); ++iHLT) {
	VSTRING_IT iName = HLTPathNames.begin();
	bool foundMatch = false;
	while ((iName != HLTPathNames.end()) && (!foundMatch)) {
	  if ((*iName).find(*iHLT) != STRING::npos) {
	    realHLT.push_back(*iName);
	    foundMatch = true;
	  }
	  ++iName;
	}
      }
    }
    else {
      STRINGSTREAM mess;
      mess << "In analyze, HLTPaths is null and category is " << *(pCategory.product()) << ".\n";
      throw cms::Exception("SkimAnalyzer") << mess.str();
    }

    //fill the vectors of pointers to vectors of histograms of event quantities
    evtHistPtrVecETMin.push_back(METETMin);
    evtHistPtrVecETMin.push_back(diEMETETMin);
    evtHistPtrVecECALIsoMax.push_back(METECALIsoMax);
    evtHistPtrVecECALIsoMax.push_back(diEMETECALIsoMax);
    evtHistPtrVecHCALIsoMax.push_back(METHCALIsoMax);
    evtHistPtrVecHCALIsoMax.push_back(diEMETHCALIsoMax);
    evtHistPtrVecHOverEMax.push_back(METHOverEMax);
    evtHistPtrVecHOverEMax.push_back(diEMETHOverEMax);
    evtHistPtrVecTrackIsoMax.push_back(METTrackIsoMax);
    evtHistPtrVecTrackIsoMax.push_back(diEMETTrackIsoMax);
    evtHistPtrVecSigmaIetaIetaMax.push_back(METSigmaIetaIetaMax);
    evtHistPtrVecSigmaIetaIetaMax.push_back(diEMETSigmaIetaIetaMax);

    //fill the vector of event quantities
    VFLOAT evtQuantityVec;
    evtQuantityVec.push_back((float)(pTCMET->begin()->et()));
    evtQuantityVec.push_back((float)(*(pDiEMET.product())));

    //fill the vectors of pointers to vectors of histograms of photon quantities
    photonHistPtrVecETMin.push_back(&ETNMinus1ETMin_);
    photonHistPtrVecETMin.push_back(&etaNMinus1ETMin_);
    photonHistPtrVecETMin.push_back(&phiNMinus1ETMin_);
    photonHistPtrVecETMin.push_back(&ECALIsoNMinus1ETMin_);
    photonHistPtrVecETMin.push_back(&HCALIsoNMinus1ETMin_);
    photonHistPtrVecETMin.push_back(&HOverENMinus1ETMin_);
    photonHistPtrVecETMin.push_back(&trackIsoNMinus1ETMin_);
    photonHistPtrVecETMin.push_back(&sigmaIetaIetaNMinus1ETMin_);
    /*photonHistPtrVecETMin.push_back(&seedTimeNMinus1ETMin_);
      photonHistPtrVecETMin.push_back(&e2OverE9NMinus1ETMin_);*/
    photonHistPtrVecECALIsoMax.push_back(&ETNMinus1ECALIsoMax_);
    photonHistPtrVecECALIsoMax.push_back(&etaNMinus1ECALIsoMax_);
    photonHistPtrVecECALIsoMax.push_back(&phiNMinus1ECALIsoMax_);
    photonHistPtrVecECALIsoMax.push_back(&ECALIsoNMinus1ECALIsoMax_);
    photonHistPtrVecECALIsoMax.push_back(&HCALIsoNMinus1ECALIsoMax_);
    photonHistPtrVecECALIsoMax.push_back(&HOverENMinus1ECALIsoMax_);
    photonHistPtrVecECALIsoMax.push_back(&trackIsoNMinus1ECALIsoMax_);
    photonHistPtrVecECALIsoMax.push_back(&sigmaIetaIetaNMinus1ECALIsoMax_);
    /*photonHistPtrVecETMin.push_back(&seedTimeNMinus1ECALIsoMax_);
      photonHistPtrVecETMin.push_back(&e2OverE9NMinus1ECALIsoMax_);*/
    photonHistPtrVecHCALIsoMax.push_back(&ETNMinus1HCALIsoMax_);
    photonHistPtrVecHCALIsoMax.push_back(&etaNMinus1HCALIsoMax_);
    photonHistPtrVecHCALIsoMax.push_back(&phiNMinus1HCALIsoMax_);
    photonHistPtrVecHCALIsoMax.push_back(&ECALIsoNMinus1HCALIsoMax_);
    photonHistPtrVecHCALIsoMax.push_back(&HCALIsoNMinus1HCALIsoMax_);
    photonHistPtrVecHCALIsoMax.push_back(&HOverENMinus1HCALIsoMax_);
    photonHistPtrVecHCALIsoMax.push_back(&trackIsoNMinus1HCALIsoMax_);
    photonHistPtrVecHCALIsoMax.push_back(&sigmaIetaIetaNMinus1HCALIsoMax_);
    /*photonHistPtrVecETMin.push_back(&seedTimeNMinus1HCALIsoMax_);
      photonHistPtrVecETMin.push_back(&e2OverE9NMinus1HCALIsoMax_);*/
    photonHistPtrVecHOverEMax.push_back(&ETNMinus1HOverEMax_);
    photonHistPtrVecHOverEMax.push_back(&etaNMinus1HOverEMax_);
    photonHistPtrVecHOverEMax.push_back(&phiNMinus1HOverEMax_);
    photonHistPtrVecHOverEMax.push_back(&ECALIsoNMinus1HOverEMax_);
    photonHistPtrVecHOverEMax.push_back(&HCALIsoNMinus1HOverEMax_);
    photonHistPtrVecHOverEMax.push_back(&HOverENMinus1HOverEMax_);
    photonHistPtrVecHOverEMax.push_back(&trackIsoNMinus1HOverEMax_);
    photonHistPtrVecHOverEMax.push_back(&sigmaIetaIetaNMinus1HOverEMax_);
    /*photonHistPtrVecETMin.push_back(&seedTimeNMinus1HOverEMax_);
      photonHistPtrVecETMin.push_back(&e2OverE9NMinus1HOverEMax_);*/
    photonHistPtrVecTrackIsoMax.push_back(&ETNMinus1TrackIsoMax_);
    photonHistPtrVecTrackIsoMax.push_back(&etaNMinus1TrackIsoMax_);
    photonHistPtrVecTrackIsoMax.push_back(&phiNMinus1TrackIsoMax_);
    photonHistPtrVecTrackIsoMax.push_back(&ECALIsoNMinus1TrackIsoMax_);
    photonHistPtrVecTrackIsoMax.push_back(&HCALIsoNMinus1TrackIsoMax_);
    photonHistPtrVecTrackIsoMax.push_back(&HOverENMinus1TrackIsoMax_);
    photonHistPtrVecTrackIsoMax.push_back(&trackIsoNMinus1TrackIsoMax_);
    photonHistPtrVecTrackIsoMax.push_back(&sigmaIetaIetaNMinus1TrackIsoMax_);
    /*photonHistPtrVecETMin.push_back(&seedTimeNMinus1TrackIsoMax_);
      photonHistPtrVecETMin.push_back(&e2OverE9NMinus1TrackIsoMax_);*/
    photonHistPtrVecSigmaIetaIetaMax.push_back(&ETNMinus1SigmaIetaIetaMax_);
    photonHistPtrVecSigmaIetaIetaMax.push_back(&etaNMinus1SigmaIetaIetaMax_);
    photonHistPtrVecSigmaIetaIetaMax.push_back(&phiNMinus1SigmaIetaIetaMax_);
    photonHistPtrVecSigmaIetaIetaMax.push_back(&ECALIsoNMinus1SigmaIetaIetaMax_);
    photonHistPtrVecSigmaIetaIetaMax.push_back(&HCALIsoNMinus1SigmaIetaIetaMax_);
    photonHistPtrVecSigmaIetaIetaMax.push_back(&HOverENMinus1SigmaIetaIetaMax_);
    photonHistPtrVecSigmaIetaIetaMax.push_back(&trackIsoNMinus1SigmaIetaIetaMax_);
    photonHistPtrVecSigmaIetaIetaMax.push_back(&sigmaIetaIetaNMinus1SigmaIetaIetaMax_);
    /*photonHistPtrVecETMin.push_back(&seedTimeNMinus1sigmaIetaIetaMax_);
      photonHistPtrVecETMin.push_back(&e2OverE9NMinus1sigmaIetaIetaMax_);*/

    /*proceed if the event passed the trigger (no matching of offline reco object to online 
      trigger object)*/
    const edm::TriggerNames& trgNames = iEvent.triggerNames(*pTrgResults);
    for (VSTRING_IT iHLT = realHLT.begin(); iHLT != realHLT.end(); ++iHLT) {
      const unsigned int trgIndex = trgNames.triggerIndex(*iHLT);
      if ((trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex))) {
	
	//fill histograms of event quantities corresponding to the various scan points
	if ((*(pCategory.product()) == GG) || (*(pCategory.product()) == EG) || 
	    (*(pCategory.product()) == EE) || (*(pCategory.product()) == FF)) {
	  fillHistograms(numPassingETMin, evtHistPtrVecETMin, evtQuantityVec, VBOOL(2, true), 
			 "ETMinScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingECALIsoMax, evtHistPtrVecECALIsoMax, evtQuantityVec, 
			 VBOOL(2, true), "ECALIsoMaxScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingHCALIsoMax, evtHistPtrVecHCALIsoMax, evtQuantityVec, 
			 VBOOL(2, true), "HCALIsoMaxScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingHOverEMax, evtHistPtrVecHOverEMax, evtQuantityVec, 
			 VBOOL(2, true), "HOverEMaxScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingTrackIsoMax, evtHistPtrVecTrackIsoMax, evtQuantityVec, 
			 VBOOL(2, true), "trackIsoMaxScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingSigmaIetaIetaMax, evtHistPtrVecSigmaIetaIetaMax, 
			 evtQuantityVec, VBOOL(2, true), "sigmaIetaIetaMaxScan", 
			 HLTPaths->at(iHLT - realHLT.begin()));
	}

	//loop over photons
	for (reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin(); 
	     iPhoton != pPhotons->end(); ++iPhoton) {

	  //plot H/E vs. |eta| of photons in events passing the trigger
	  HOverEVsAbsEta_->Fill(fabs(iPhoton->eta()), iPhoton->hadronicOverEm());

	  //fill the vector of photon quantities
	  VFLOAT photonQuantityVec;
	  photonQuantityVec.push_back((float)(iPhoton->et()));
	  photonQuantityVec.push_back((float)(iPhoton->eta()));
	  photonQuantityVec.push_back((float)(iPhoton->phi()));
	  photonQuantityVec.push_back((float)(iPhoton->ecalRecHitSumEtConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->hcalTowerSumEtConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->hadronicOverEm()));
	  photonQuantityVec.push_back((float)(iPhoton->trkSumPtHollowConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->sigmaIetaIeta()));
	  /*photonQuantityVec.push_back(absSeedTimeMax);
	    photonQuantityVec.push_back(e2OverE9Max);*/

	  /*fill the vector of photon pass flags
	    photon passes if it passes all selection cuts except the one being plotted*/
	  VBOOL passVec;
	  /*edm::Ref<reco::PhotonCollection> ref(pPhotons, iPhoton - pPhotons->begin());
	  passVec.push_back(passNMinus1(passMapVec, ref, &pPassETMin));
	  passVec.push_back(passNMinus1(passMapVec, ref, &pPassAbsEtaMax));
	  passVec.push_back(passNMinus1(passMapVec, ref, NULL));
	  passVec.push_back(passNMinus1(passMapVec, ref, &pPassECALIsoMax));
	  passVec.push_back(passNMinus1(passMapVec, ref, &pPassHCALIsoMax));
	  passVec.push_back(passNMinus1(passMapVec, ref, &pPassHOverEMax));
	  passVec.push_back(passNMinus1(passMapVec, ref, &pPassTrackIsoMax));
	  passVec.push_back(passNMinus1(passMapVec, ref, &pPassSigmaIetaIetaMax));*/
	  /*passVec.push_back(passNMinus1(passMapVec, ref, &pPassAbsSeedTimeMax));
	    passVec.push_back(passNMinus1(passMapVec, ref, &pPassE2OverE9Max));*/

	  //hack to plot distributions with no cuts at all
	  passVec.push_back(true);
	  passVec.push_back(true);
	  passVec.push_back(true);
	  passVec.push_back(true);
	  passVec.push_back(true);
	  passVec.push_back(true);
	  passVec.push_back(true);
	  passVec.push_back(true);

	  //fill histograms of photon quantities corresponding to the various scan points
	  fillHistograms(numPassingETMin, photonHistPtrVecETMin, photonQuantityVec, 
			 passVec, "ETMinScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingECALIsoMax, photonHistPtrVecECALIsoMax, photonQuantityVec, 
			 passVec, "ECALIsoMaxScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingHCALIsoMax, photonHistPtrVecHCALIsoMax, photonQuantityVec, 
			 passVec, "HCALIsoMaxScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingHOverEMax, photonHistPtrVecHOverEMax, photonQuantityVec, 
			 passVec, "HOverEMaxScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingTrackIsoMax, photonHistPtrVecTrackIsoMax, photonQuantityVec, 
			 passVec, "trackIsoMaxScan", HLTPaths->at(iHLT - realHLT.begin()));
	  fillHistograms(numPassingSigmaIetaIetaMax, photonHistPtrVecSigmaIetaIetaMax, 
			 photonQuantityVec, passVec, "sigmaIetaIetaMaxScan", 
			 HLTPaths->at(iHLT - realHLT.begin()));
	}
      }
    }
  }
}

// ---- method called once each run  ---
void SkimAnalyzer::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
  bool changed = false;
  if (!HLTCfg_.init(iRun, iSetup, HLTProcessName_, changed) ) {
    edm::LogError("SkimAnalyzer") << "Error: can't initialize HLTConfigProvider.\n";
    throw cms::Exception("SkimAnalyzer") << "HLTConfigProvider::init() returned non-0.\n";
  }
  if (changed) edm::LogInfo("SkimAnalyzer") << "HLT configuration changed!\n";
}

// ------------ method called once each job just before starting event loop  ------------
void 
SkimAnalyzer::beginJob()
{
  //book histograms, 1 for each plot-category-scan-trigger
  for (VSTRING_IT iGGHLT = HLTPathGG_.begin(); iGGHLT != HLTPathGG_.end(); ++iGGHLT) {
    bookGGMETHistograms(*iGGHLT);
    bookEGMETHistograms(*iGGHLT);
    bookEEMETHistograms(*iGGHLT);
    bookGGDiEMETHistograms(*iGGHLT);
    bookEGDiEMETHistograms(*iGGHLT);
    bookEEDiEMETHistograms(*iGGHLT);
    bookETNMinus1Histograms(*iGGHLT);
    bookEtaNMinus1Histograms(*iGGHLT);
    bookPhiNMinus1Histograms(*iGGHLT);
    bookECALIsoNMinus1Histograms(*iGGHLT);
    bookHCALIsoNMinus1Histograms(*iGGHLT);
    bookHOverENMinus1Histograms(*iGGHLT);
    bookTrackIsoNMinus1Histograms(*iGGHLT);
    bookSigmaIetaIetaNMinus1Histograms(*iGGHLT);
  }
  for (VSTRING_IT iFFHLT = HLTPathFF_.begin(); iFFHLT != HLTPathFF_.end(); ++iFFHLT) {
    bookFFMETHistograms(*iFFHLT);
    bookFFDiEMETHistograms(*iFFHLT);
    bookETNMinus1Histograms(*iFFHLT);
    bookEtaNMinus1Histograms(*iFFHLT);
    bookPhiNMinus1Histograms(*iFFHLT);
    bookECALIsoNMinus1Histograms(*iFFHLT);
    bookHCALIsoNMinus1Histograms(*iFFHLT);
    bookHOverENMinus1Histograms(*iFFHLT);
    bookTrackIsoNMinus1Histograms(*iFFHLT);
    bookSigmaIetaIetaNMinus1Histograms(*iFFHLT);
  }

  //hack
  HOverEVsAbsEta_ = new TH2F("", "", 30, 0.0, 3.0, 20, 0.0, 0.2);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimAnalyzer::endJob() {

  //open output file
  out_ = new TFile(outputFile_.c_str(), "RECREATE");
  if (out_->IsOpen()) {
    out_->cd();

    //write ETMin scan plots
    makeOutputDirectories(ETMinScan_.size(), "ETMinScan");
    fillOutputDirectories(ggMETETMin_, "ETMinScan");
    fillOutputDirectories(ffMETETMin_, "ETMinScan");
    fillOutputDirectories(eeMETETMin_, "ETMinScan");
    fillOutputDirectories(egMETETMin_, "ETMinScan");
    fillOutputDirectories(ggDiEMETETMin_, "ETMinScan");
    fillOutputDirectories(ffDiEMETETMin_, "ETMinScan");
    fillOutputDirectories(eeDiEMETETMin_, "ETMinScan");
    fillOutputDirectories(egDiEMETETMin_, "ETMinScan");
    fillOutputDirectories(ETNMinus1ETMin_, "ETMinScan");
    fillOutputDirectories(etaNMinus1ETMin_, "ETMinScan");
    fillOutputDirectories(phiNMinus1ETMin_, "ETMinScan");
    fillOutputDirectories(ECALIsoNMinus1ETMin_, "ETMinScan");
    fillOutputDirectories(HCALIsoNMinus1ETMin_, "ETMinScan");
    fillOutputDirectories(HOverENMinus1ETMin_, "ETMinScan");
    fillOutputDirectories(trackIsoNMinus1ETMin_, "ETMinScan");
    fillOutputDirectories(sigmaIetaIetaNMinus1ETMin_, "ETMinScan");
    /*fillOutputDirectories(seedTimeNMinus1ETMin_, "ETMinScan");
      fillOutputDirectories(e2OverE9NMinus1ETMin_, "ETMinScan");*/

    //write ECALIsoMax scan plots
    makeOutputDirectories(ECALIsoMaxPTMultiplierScan_.size(), "ECALIsoMaxScan");
    fillOutputDirectories(ggMETECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(ffMETECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(eeMETECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(egMETECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(ggDiEMETECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(ffDiEMETECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(eeDiEMETECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(egDiEMETECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(ETNMinus1ECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(etaNMinus1ECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(phiNMinus1ECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(ECALIsoNMinus1ECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(HCALIsoNMinus1ECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(HOverENMinus1ECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(trackIsoNMinus1ECALIsoMax_, "ECALIsoMaxScan");
    fillOutputDirectories(sigmaIetaIetaNMinus1ECALIsoMax_, "ECALIsoMaxScan");
    /*fillOutputDirectories(seedTimeNMinus1ECALIsoMax_, "ECALIsoMaxScan");
      fillOutputDirectories(e2OverE9NMinus1ECALIsoMax_, "ECALIsoMaxScan");*/

    //write HCALIsoMax scan plots
    makeOutputDirectories(HCALIsoMaxPTMultiplierScan_.size(), "HCALIsoMaxScan");
    fillOutputDirectories(ggMETHCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(ffMETHCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(eeMETHCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(egMETHCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(ggDiEMETHCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(ffDiEMETHCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(eeDiEMETHCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(egDiEMETHCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(ETNMinus1HCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(etaNMinus1HCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(phiNMinus1HCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(ECALIsoNMinus1HCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(HCALIsoNMinus1HCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(HOverENMinus1HCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(trackIsoNMinus1HCALIsoMax_, "HCALIsoMaxScan");
    fillOutputDirectories(sigmaIetaIetaNMinus1HCALIsoMax_, "HCALIsoMaxScan");
    /*fillOutputDirectories(seedTimeNMinus1HCALIsoMax_, "HCALIsoMaxScan");
      fillOutputDirectories(e2OverE9NMinus1HCALIsoMax_, "HCALIsoMaxScan");*/

    //write HOverEMax scan plots
    makeOutputDirectories(HOverEMaxScan_.size(), "HOverEMaxScan");
    fillOutputDirectories(ggMETHOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(ffMETHOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(eeMETHOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(egMETHOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(ggDiEMETHOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(ffDiEMETHOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(eeDiEMETHOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(egDiEMETHOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(ETNMinus1HOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(etaNMinus1HOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(phiNMinus1HOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(ECALIsoNMinus1HOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(HCALIsoNMinus1HOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(HOverENMinus1HOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(trackIsoNMinus1HOverEMax_, "HOverEMaxScan");
    fillOutputDirectories(sigmaIetaIetaNMinus1HOverEMax_, "HOverEMaxScan");
    /*fillOutputDirectories(seedTimeNMinus1HOverEMax_, "HOverEMaxScan");
      fillOutputDirectories(e2OverE9NMinus1HOverEMax_, "HOverEMaxScan");*/

    //write trackIsoMax scan plots
    makeOutputDirectories(trackIsoMaxPTMultiplierScan_.size(), "trackIsoMaxScan");
    fillOutputDirectories(ggMETTrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(ffMETTrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(eeMETTrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(egMETTrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(ggDiEMETTrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(ffDiEMETTrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(eeDiEMETTrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(egDiEMETTrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(ETNMinus1TrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(etaNMinus1TrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(phiNMinus1TrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(ECALIsoNMinus1TrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(HCALIsoNMinus1TrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(HOverENMinus1TrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(trackIsoNMinus1TrackIsoMax_, "trackIsoMaxScan");
    fillOutputDirectories(sigmaIetaIetaNMinus1TrackIsoMax_, "trackIsoMaxScan");
    /*fillOutputDirectories(seedTimeNMinus1TrackIsoMax_, "trackIsoMaxScan");
      fillOutputDirectories(e2OverE9NMinus1TrackIsoMax_, "trackIsoMaxScan");*/

    //write sigmaIetaIeta scan plots
    makeOutputDirectories(sigmaIetaIetaMaxScan_.size(), "sigmaIetaIetaMaxScan");
    fillOutputDirectories(ggMETSigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(ffMETSigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(eeMETSigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(egMETSigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(ggDiEMETSigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(ffDiEMETSigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(eeDiEMETSigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(egDiEMETSigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(ETNMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(etaNMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(phiNMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(ECALIsoNMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(HCALIsoNMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(HOverENMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(trackIsoNMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    fillOutputDirectories(sigmaIetaIetaNMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
    /*fillOutputDirectories(seedTimeNMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");
      fillOutputDirectories(e2OverE9NMinus1SigmaIetaIetaMax_, "sigmaIetaIetaMaxScan");*/

    //hack
    out_->cd("/");
    HOverEVsAbsEta_->Write();

    //write output file
    out_->cd("/"/*extractFileFromFullPath(outputFile_).c_str()*/);
    out_->Write();
    out_->Close();
  }
  else edm::LogError("SkimAnalyzer") << "Error opening file " << outputFile_ << ".\n";
  delete out_;
}

STRING SkimAnalyzer::name(const STRING& category, const STRING& plot, const STRING& scanName, 
			  const unsigned int scanIndex, const STRING& HLT) const
{
  STRINGSTREAM name;
  name << category << "_" << plot << "_" << scanName << scanIndex << "_" << HLT;
  return name.str();
}

void SkimAnalyzer::bookScanPointHistograms(VTH1F& histVec, const VDOUBLE& scanVec, 
					   const STRING& category, const STRING& plot, 
					   const STRING& scanName, const STRING& HLT, 
					   const unsigned int numBins, const float binLowEdge, 
					   const float binHighEdge)
{
  for (VDOUBLE_IT i = scanVec.begin(); i != scanVec.end(); ++i) {
    histVec.push_back(new TH1F(name(category, plot, scanName, i - scanVec.begin(), HLT).c_str(), 
			       "", numBins, binLowEdge, binHighEdge));
  }
}

void SkimAnalyzer::bookGGMETHistograms(const STRING& HLT)
{
  bookScanPointHistograms(ggMETETMin_, ETMinScan_, "gg", "MET", "ETMinScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggMETECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "gg", "MET", 
			  "ECALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggMETHCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "gg", "MET", 
			  "HCALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggMETHOverEMax_, HOverEMaxScan_, "gg", "MET", 
			  "HOverEMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggMETTrackIsoMax_, trackIsoMaxPTMultiplierScan_, "gg", "MET", 
			  "trackIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggMETSigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "gg", "MET", 
			  "sigmaIetaIetaMaxScan", HLT, 30, 0.0, 150.0);
}

void SkimAnalyzer::bookFFMETHistograms(const STRING& HLT)
{
  bookScanPointHistograms(ffMETETMin_, ETMinScan_, "ff", "MET", "ETMinScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffMETECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "ff", "MET", 
			  "ECALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffMETHCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "ff", "MET", 
			  "HCALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffMETHOverEMax_, HOverEMaxScan_, "ff", "MET", 
			  "HOverEMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffMETTrackIsoMax_, trackIsoMaxPTMultiplierScan_, "ff", "MET", 
			  "trackIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffMETSigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "ff", "MET", 
			  "sigmaIetaIetaMaxScan", HLT, 30, 0.0, 150.0);
}

void SkimAnalyzer::bookEGMETHistograms(const STRING& HLT)
{
  bookScanPointHistograms(egMETETMin_, ETMinScan_, "eg", "MET", "ETMinScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egMETECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "eg", "MET", 
			  "ECALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egMETHCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "eg", "MET", 
			  "HCALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egMETHOverEMax_, HOverEMaxScan_, "eg", "MET", 
			  "HOverEMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egMETTrackIsoMax_, trackIsoMaxPTMultiplierScan_, "eg", "MET", 
			  "trackIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egMETSigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "eg", "MET", 
			  "sigmaIetaIetaMaxScan", HLT, 30, 0.0, 150.0);
}

void SkimAnalyzer::bookEEMETHistograms(const STRING& HLT)
{
  bookScanPointHistograms(eeMETETMin_, ETMinScan_, "ee", "MET", "ETMinScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeMETECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "ee", "MET", 
			  "ECALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeMETHCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "ee", "MET", 
			  "HCALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeMETHOverEMax_, HOverEMaxScan_, "ee", "MET", 
			  "HOverEMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeMETTrackIsoMax_, trackIsoMaxPTMultiplierScan_, "ee", "MET", 
			  "trackIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeMETSigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "ee", "MET", 
			  "sigmaIetaIetaMaxScan", HLT, 30, 0.0, 150.0);
}

void SkimAnalyzer::bookGGDiEMETHistograms(const STRING& HLT)
{
  bookScanPointHistograms(ggDiEMETETMin_, ETMinScan_, "gg", "diEMET", "ETMinScan", HLT, 30, 0.0, 
			  150.0);
  bookScanPointHistograms(ggDiEMETECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "gg", "diEMET", 
			  "ECALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggDiEMETHCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "gg", "diEMET", 
			  "HCALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggDiEMETHOverEMax_, HOverEMaxScan_, "gg", "diEMET", 
			  "HOverEMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggDiEMETTrackIsoMax_, trackIsoMaxPTMultiplierScan_, "gg", "diEMET", 
			  "trackIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ggDiEMETSigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "gg", "diEMET", 
			  "sigmaIetaIetaMaxScan", HLT, 30, 0.0, 150.0);
}

void SkimAnalyzer::bookFFDiEMETHistograms(const STRING& HLT)
{
  bookScanPointHistograms(ffDiEMETETMin_, ETMinScan_, "ff", "diEMET", "ETMinScan", HLT, 30, 0.0, 
			  150.0);
  bookScanPointHistograms(ffDiEMETECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "ff", "diEMET", 
			  "ECALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffDiEMETHCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "ff", "diEMET", 
			  "HCALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffDiEMETHOverEMax_, HOverEMaxScan_, "ff", "diEMET", 
			  "HOverEMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffDiEMETTrackIsoMax_, trackIsoMaxPTMultiplierScan_, "ff", "diEMET", 
			  "trackIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(ffDiEMETSigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "ff", "diEMET", 
			  "sigmaIetaIetaMaxScan", HLT, 30, 0.0, 150.0);
}

void SkimAnalyzer::bookEGDiEMETHistograms(const STRING& HLT)
{
  bookScanPointHistograms(egDiEMETETMin_, ETMinScan_, "eg", "diEMET", "ETMinScan", HLT, 30, 0.0, 
			  150.0);
  bookScanPointHistograms(egDiEMETECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "eg", "diEMET", 
			  "ECALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egDiEMETHCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "eg", "diEMET", 
			  "HCALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egDiEMETHOverEMax_, HOverEMaxScan_, "eg", "diEMET", 
			  "HOverEMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egDiEMETTrackIsoMax_, trackIsoMaxPTMultiplierScan_, "eg", "diEMET", 
			  "trackIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(egDiEMETSigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "eg", "diEMET", 
			  "sigmaIetaIetaMaxScan", HLT, 30, 0.0, 150.0);
}

void SkimAnalyzer::bookEEDiEMETHistograms(const STRING& HLT)
{
  bookScanPointHistograms(eeDiEMETETMin_, ETMinScan_, "ee", "diEMET", "ETMinScan", HLT, 30, 0.0, 
			  150.0);
  bookScanPointHistograms(eeDiEMETECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "ee", "diEMET", 
			  "ECALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeDiEMETHCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "ee", "diEMET", 
			  "HCALIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeDiEMETHOverEMax_, HOverEMaxScan_, "ee", "diEMET", 
			  "HOverEMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeDiEMETTrackIsoMax_, trackIsoMaxPTMultiplierScan_, "ee", "diEMET", 
			  "trackIsoMaxScan", HLT, 30, 0.0, 150.0);
  bookScanPointHistograms(eeDiEMETSigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "ee", "diEMET", 
			  "sigmaIetaIetaMaxScan", HLT, 30, 0.0, 150.0);
}

void SkimAnalyzer::bookETNMinus1Histograms(const STRING& HLT)
{
  bookScanPointHistograms(ETNMinus1ETMin_, ETMinScan_, "", "ET", "ETMinScan", HLT, 80, 0.0, 
			  800.0);
  bookScanPointHistograms(ETNMinus1ECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "", "ET", 
			  "ECALIsoMaxScan", HLT, 80, 0.0, 800.0);
  bookScanPointHistograms(ETNMinus1HCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "", "ET", 
			  "HCALIsoMaxScan", HLT, 80, 0.0, 800.0);
  bookScanPointHistograms(ETNMinus1HOverEMax_, HOverEMaxScan_, "", "ET", 
			  "HOverEMaxScan", HLT, 80, 0.0, 800.0);
  bookScanPointHistograms(ETNMinus1TrackIsoMax_, trackIsoMaxPTMultiplierScan_, "", "ET", 
			  "trackIsoMaxScan", HLT, 80, 0.0, 800.0);
  bookScanPointHistograms(ETNMinus1SigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "", "ET", 
			  "sigmaIetaIetaMaxScan", HLT, 80, 0.0, 800.0);
}

void SkimAnalyzer::bookEtaNMinus1Histograms(const STRING& HLT)
{
  bookScanPointHistograms(etaNMinus1ETMin_, ETMinScan_, "", "eta", "ETMinScan", HLT, 60, -3.0, 
			  3.0);
  bookScanPointHistograms(etaNMinus1ECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "", "eta", 
			  "ECALIsoMaxScan", HLT, 60, -3.0, 3.0);
  bookScanPointHistograms(etaNMinus1HCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "", "eta", 
			  "HCALIsoMaxScan", HLT, 60, -3.0, 3.0);
  bookScanPointHistograms(etaNMinus1HOverEMax_, HOverEMaxScan_, "", "eta", 
			  "HOverEMaxScan", HLT, 60, -3.0, 3.0);
  bookScanPointHistograms(etaNMinus1TrackIsoMax_, trackIsoMaxPTMultiplierScan_, "", "eta", 
			  "trackIsoMaxScan", HLT, 60, -3.0, 3.0);
  bookScanPointHistograms(etaNMinus1SigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "", "eta", 
			  "sigmaIetaIetaMaxScan", HLT, 60, -3.0, 3.0);
}

void SkimAnalyzer::bookPhiNMinus1Histograms(const STRING& HLT)
{
  bookScanPointHistograms(phiNMinus1ETMin_, ETMinScan_, "", "phi", "ETMinScan", HLT, 64, -3.2, 
			  3.2);
  bookScanPointHistograms(phiNMinus1ECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "", "phi", 
			  "ECALIsoMaxScan", HLT, 64, -3.2, 3.2);
  bookScanPointHistograms(phiNMinus1HCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "", "phi", 
			  "HCALIsoMaxScan", HLT, 64, -3.2, 3.2);
  bookScanPointHistograms(phiNMinus1HOverEMax_, HOverEMaxScan_, "", "phi", 
			  "HOverEMaxScan", HLT, 64, -3.2, 3.2);
  bookScanPointHistograms(phiNMinus1TrackIsoMax_, trackIsoMaxPTMultiplierScan_, "", "phi", 
			  "trackIsoMaxScan", HLT, 64, -3.2, 3.2);
  bookScanPointHistograms(phiNMinus1SigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "", "phi", 
			  "sigmaIetaIetaMaxScan", HLT, 64, -3.2, 3.2);
}

void SkimAnalyzer::bookECALIsoNMinus1Histograms(const STRING& HLT)
{
  bookScanPointHistograms(ECALIsoNMinus1ETMin_, ETMinScan_, "", "ECALIso", 
			  "ETMinScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(ECALIsoNMinus1ECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "", "ECALIso", 
			  "ECALIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(ECALIsoNMinus1HCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "", "ECALIso", 
			  "HCALIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(ECALIsoNMinus1HOverEMax_, HOverEMaxScan_, "", "ECALIso", 
			  "HOverEMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(ECALIsoNMinus1TrackIsoMax_, trackIsoMaxPTMultiplierScan_, "", "ECALIso", 
			  "trackIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(ECALIsoNMinus1SigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "", "ECALIso", 
			  "sigmaIetaIetaMaxScan", HLT, 20, 0.0, 40.0);
}

void SkimAnalyzer::bookHCALIsoNMinus1Histograms(const STRING& HLT)
{
  bookScanPointHistograms(HCALIsoNMinus1ETMin_, ETMinScan_, "", "HCALIso", 
			  "ETMinScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(HCALIsoNMinus1ECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "", "HCALIso", 
			  "ECALIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(HCALIsoNMinus1HCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "", "HCALIso", 
			  "HCALIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(HCALIsoNMinus1HOverEMax_, HOverEMaxScan_, "", "HCALIso", 
			  "HOverEMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(HCALIsoNMinus1TrackIsoMax_, trackIsoMaxPTMultiplierScan_, "", "HCALIso", 
			  "trackIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(HCALIsoNMinus1SigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "", "HCALIso", 
			  "sigmaIetaIetaMaxScan", HLT, 20, 0.0, 40.0);
}

void SkimAnalyzer::bookHOverENMinus1Histograms(const STRING& HLT)
{
  bookScanPointHistograms(HOverENMinus1ETMin_, ETMinScan_, "", "HOverE", 
			  "ETMinScan", HLT, 20, 0.0, 0.2);
  bookScanPointHistograms(HOverENMinus1ECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "", "HOverE", 
			  "ECALIsoMaxScan", HLT, 20, 0.0, 0.2);
  bookScanPointHistograms(HOverENMinus1HCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "", "HOverE", 
			  "HCALIsoMaxScan", HLT, 20, 0.0, 0.2);
  bookScanPointHistograms(HOverENMinus1HOverEMax_, HOverEMaxScan_, "", "HOverE", 
			  "HOverEMaxScan", HLT, 20, 0.0, 0.2);
  bookScanPointHistograms(HOverENMinus1TrackIsoMax_, trackIsoMaxPTMultiplierScan_, "", "HOverE", 
			  "trackIsoMaxScan", HLT, 20, 0.0, 0.2);
  bookScanPointHistograms(HOverENMinus1SigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "", "HOverE", 
			  "sigmaIetaIetaMaxScan", HLT, 20, 0.0, 0.2);
}

void SkimAnalyzer::bookTrackIsoNMinus1Histograms(const STRING& HLT)
{
  bookScanPointHistograms(trackIsoNMinus1ETMin_, ETMinScan_, "", "trackIso", 
			  "ETMinScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(trackIsoNMinus1ECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "", "trackIso", 
			  "ECALIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(trackIsoNMinus1HCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "", 
			  "trackIso", "HCALIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(trackIsoNMinus1HOverEMax_, HOverEMaxScan_, "", "trackIso", 
			  "HOverEMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(trackIsoNMinus1TrackIsoMax_, trackIsoMaxPTMultiplierScan_, "", 
			  "trackIso", "trackIsoMaxScan", HLT, 20, 0.0, 40.0);
  bookScanPointHistograms(trackIsoNMinus1SigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "", "trackIso", 
			  "sigmaIetaIetaMaxScan", HLT, 20, 0.0, 40.0);
}

void SkimAnalyzer::bookSigmaIetaIetaNMinus1Histograms(const STRING& HLT)
{
  bookScanPointHistograms(sigmaIetaIetaNMinus1ETMin_, ETMinScan_, "", "sigmaIetaIeta", 
			  "ETMinScan", HLT, 20, 0.0, 0.02);
  bookScanPointHistograms(sigmaIetaIetaNMinus1ECALIsoMax_, ECALIsoMaxPTMultiplierScan_, "", 
			  "sigmaIetaIeta", "ECALIsoMaxScan", HLT, 20, 0.0, 0.02);
  bookScanPointHistograms(sigmaIetaIetaNMinus1HCALIsoMax_, HCALIsoMaxPTMultiplierScan_, "", 
			  "sigmaIetaIeta", "HCALIsoMaxScan", HLT, 20, 0.0, 0.02);
  bookScanPointHistograms(sigmaIetaIetaNMinus1HOverEMax_, HOverEMaxScan_, "", "sigmaIetaIeta", 
			  "HOverEMaxScan", HLT, 20, 0.0, 0.02);
  bookScanPointHistograms(sigmaIetaIetaNMinus1TrackIsoMax_, trackIsoMaxPTMultiplierScan_, "", 
			  "sigmaIetaIeta", "trackIsoMaxScan", HLT, 20, 0.0, 0.02);
  bookScanPointHistograms(sigmaIetaIetaNMinus1SigmaIetaIetaMax_, sigmaIetaIetaMaxScan_, "", 
			  "sigmaIetaIeta", "sigmaIetaIetaMaxScan", HLT, 20, 0.0, 0.02);
}

void SkimAnalyzer::incrementCounters(const VDOUBLE& scanVec, VUINT& counterVec, 
					const unsigned int decisionWord) const
{
  for (unsigned int i = 0; i < scanVec.size(); ++i) {
    if (((decisionWord >> i) & 0x1) == 1) ++counterVec[i];
  }
}

int SkimAnalyzer::scanPoint(const STRING& histName, const STRING& scanLabel) const
{
  int scanPointInt = -1;
  size_t scanLabelPos = histName.find(scanLabel);
  size_t underscorePos = histName.find('_', histName.find(scanLabel));
  if ((scanLabelPos != STRING::npos) && (underscorePos != STRING::npos) && 
      (underscorePos > scanLabelPos)) {
    STRING scanPoint = histName.substr(scanLabelPos + scanLabel.length(), 
				       underscorePos - histName.find(scanLabel) - 
				       scanLabel.length());
    STRINGSTREAM scanPointStream;
    scanPointStream << scanPoint;
    scanPointStream >> scanPointInt;
  }
  return scanPointInt;
}

void SkimAnalyzer::fillIndHistogram(VTH1F* histVec, const float quantity, 
				    const bool pass, const STRING& scanLabel, const STRING& HLT, 
				    const unsigned int currentScanPoint)
{
  if (histVec != NULL) {
    for (VTH1F_IT iHist = histVec->begin(); iHist != histVec->end(); ++iHist) {
      STRING histName((*iHist)->GetName());
      int theScanPoint = scanPoint(histName, scanLabel);
      if (theScanPoint >= 0) {

	//if the current histogram corresponds to this trigger and scan point, fill it
	if ((histName.find(HLT) != STRING::npos) && (theScanPoint == (int)currentScanPoint) && 
	    pass) {
	  (*iHist)->Fill(quantity);
	}
      }
      else {
	STRINGSTREAM err;
	err << "Error: unable to parse histogram name " << histName << ".\n";
	edm::LogError("SkimAnalzyer") << err.str();
      }
    }
  }
  else edm::LogInfo("SkimAnalyzerNullPointer") << "In fillIndHistogram, histVec is null.\n";
}

void SkimAnalyzer::fillHistograms(const VUINT& counterVec, 
				  std::vector<VTH1F*>& histPtrVec, 
				  const VFLOAT& quantityVec, 
				  const VBOOL& passVec, const STRING& scanLabel, const STRING& HLT)
{
  /*loop over the vector of scan point counters -- each element holds the number of photons 
    passing the scan point given by the element's index*/
  for (VUINT_IT iCount = counterVec.begin(); iCount != counterVec.end(); ++iCount) {

    //proceed if >=2 photons passed the scan point
    if (*iCount >= 2) {

      //loop over the vector of pointers to histograms of physical quantities
      for (std::vector<VTH1F*>::iterator iHistPtrVec = histPtrVec.begin(); 
	   iHistPtrVec != histPtrVec.end(); ++iHistPtrVec) {

	//fill histograms
	fillIndHistogram(*iHistPtrVec, quantityVec[iHistPtrVec - histPtrVec.begin()], 
			 passVec[iHistPtrVec - histPtrVec.begin()], scanLabel, HLT, 
			 iCount - counterVec.begin());
      }
    }
  }
}

bool SkimAnalyzer::passNMinus1(std::vector<edm::Handle<edm::ValueMap<bool> >*>& passMaps, 
			       edm::Ref<reco::PhotonCollection>& ref, 
			       edm::Handle<edm::ValueMap<bool> >* passMapToIgnore) const
{
  bool passFlag = true;
  for (std::vector<edm::Handle<edm::ValueMap<bool> >*>::const_iterator iPassMap = 
	 passMaps.begin(); iPassMap != passMaps.end(); ++iPassMap) {
    passFlag = passFlag && ((*iPassMap == passMapToIgnore) || (*(*iPassMap)->product())[ref]);
  }
  return passFlag;
}

STRING SkimAnalyzer::extractFileFromFullPath(const STRING& fullPath) const
{
  STRING file = "";
  size_t lastSlashPos = fullPath.rfind('/');
  if (lastSlashPos != STRING::npos) file = fullPath.substr(lastSlashPos + 1);
  return file;
}

void SkimAnalyzer::makeOutputDirectories(const unsigned int numScanPoints, const STRING& scanLabel)
{
  out_->cd();
  out_->mkdir(scanLabel.c_str());
  for (unsigned int i = 0; i < numScanPoints; ++i) {
    STRINGSTREAM dirName;
    dirName << scanLabel << i;
    out_->cd(scanLabel.c_str());
    out_->mkdir(dirName.str().c_str());
    out_->cd(dirName.str().c_str());
  }
}

void SkimAnalyzer::fillOutputDirectories(VTH1F& histVec, const STRING& scanLabel)
{
  for (VTH1F_IT iHist = histVec.begin(); iHist != histVec.end(); ++iHist) {
    STRING histName((*iHist)->GetName());
    int theScanPoint = scanPoint(histName, scanLabel);
    STRINGSTREAM targetDir1;
    targetDir1 << "/" << scanLabel;
    out_->cd(targetDir1.str().c_str());
    STRINGSTREAM targetDir2;
    targetDir2 << scanLabel << theScanPoint;
    out_->cd(targetDir2.str().c_str());
    (*iHist)->Write();
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SkimAnalyzer);
