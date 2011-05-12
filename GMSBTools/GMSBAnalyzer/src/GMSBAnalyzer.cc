// -*- C++ -*-
//
// Package:    GMSBTools
// Class:      GMSBAnalyzer
// 
/*
 Description: 2 photons + MET analysis

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay
//         Created:  Tue May 3 10:56:33 CET 2011
// $Id: GMSBAnalyzer.cc,v 1.1 2011/05/03 10:56:33 yohay Exp $
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
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
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

class GMSBAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GMSBAnalyzer(const edm::ParameterSet&);
      ~GMSBAnalyzer();

   private:
      virtual void beginRun(const edm::Run&, edm::EventSetup const&);
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  //return a properly formatted histogram name
  STRING name(const STRING&, const STRING&, const STRING&) const;

  //book 1D histogram
  void bookHistogram(VTH1F&, const STRING&, const STRING&, const STRING&, const unsigned int, 
		     const float, const float);

  //book 2D histogram
  void bookHistogram(VTH2F&, const STRING&, const STRING&, const STRING&, const unsigned int, 
		     const float, const float, const unsigned int, const float, const float);

  /*book all histograms filled once per event for the samples passing sigmaIetaIeta and track 
    isolation (those with E and G)*/
  void bookEAndGHistograms(const STRING&);

  /*book all histograms filled once per event for the samples failing either sigmaIetaIeta or 
    track isolation (those with F)*/
  void bookFHistograms(const STRING&);

  //book all histograms filled once per photon
  void bookPhotonHistograms(const STRING&);

  //get real trigger names from approximate inputs
  VSTRING getRealTriggerNames(VSTRING*) const;

  //set pointers to the 2 passing photons, ordered by ET
  void getPassingPhotons(edm::Handle<reco::PhotonCollection>&, reco::Photon&, reco::Photon&, 
			 unsigned int&, unsigned int&, edm::Handle<edm::ValueMap<bool> >&);

  //fill individual 1D event or photon histogram corresponding to supplied trigger
  void fillIndHistogram(VTH1F*, const float, const bool, const STRING&);

  //fill individual 2D event or photon histogram corresponding to supplied trigger
  void fillIndHistogram(VTH2F*, const float, const float, const bool, const STRING&);

  //fill 1D histograms of event and photon quantities properly
  void fillHistograms(std::vector<VTH1F*>&, const VFLOAT&, const VBOOL&, const STRING&);

  //fill 2D histograms of event and photon quantities properly
  void fillHistograms(std::vector<VTH2F*>&, const VFLOAT&, const VFLOAT&, const VBOOL&, 
		      const STRING&);

  //return true if the photon should go into the N-1 plot, false otherwise
  bool passNMinus1(std::vector<edm::Handle<edm::ValueMap<bool> >*>&, 
		   edm::Ref<reco::PhotonCollection>&, edm::Handle<edm::ValueMap<bool> >*) const;

  //return true if >=2 photons should go into the N-1 plot, false otherwise
  bool twoPhotonsPassNMinus1(std::vector<edm::Handle<edm::ValueMap<bool> >*>&, 
			     edm::Handle<edm::ValueMap<bool> >*, 
			     const edm::Handle<reco::PhotonCollection>&) const;

  /*return number of PVs passing "good PV criteria" (HEAD of 
    DPGAnalysis/Skims/python/GoodVertex_cfg.py)*/
  unsigned int numValidPVs(edm::Handle<reco::VertexCollection>&) const;

  //if plot should be filled for efficiency calculation, return NPV; else return 0
  unsigned int effDecision(std::vector<edm::Handle<edm::ValueMap<bool> >*>&, 
			   edm::Ref<reco::PhotonCollection>&, edm::Handle<edm::ValueMap<bool> >*, 
			   const unsigned int, 
			   const edm::Handle<reco::PhotonCollection>& pPhotons = 
			   edm::Handle<reco::PhotonCollection>()) const;

  //make ROOT directories in the output file corresponding to the different event samples
  void makeOutputDirectories(const STRING&);

  //calculate efficiencies stored in TGraphAsymmErrors objects
  void calculateEff(VTGRAPHASYMMERRORS&, VTH1F&, VTH1F&);

  //write gg plots
  void writeGGPlots();

  //write eg plots
  void writeEGPlots();

  //write ee plots
  void writeEEPlots();

  //write ff plots
  void writeFFPlots();

  //write photon plots
  void writePhotonPlots();

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

  //fill ROOT directories in the output file with corresponding histograms
  template<typename T>
  void fillOutputDirectories(T& histVec, const STRING& dirName)
  {
    out_->cd();
    out_->cd(dirName.c_str());
    for (typename T::iterator iHist = histVec.begin(); iHist != histVec.end(); ++iHist) { 
      (*iHist)->Write();
    }
  }

      // ----------member data ---------------------------

  //input
  STRING outputFile_;
  edm::InputTag photonTag_;
  edm::InputTag categoryTag_;
  edm::InputTag diEMETTag_;
  edm::InputTag photonSeedTimeTag_;
  edm::InputTag photonE2OverE9Tag_;
  edm::InputTag passAllTag_;
  edm::InputTag trgResultsTag_;
  edm::InputTag tcMETTag_;
  edm::InputTag vertexTag_;
  edm::InputTag passETMinTag_;
  edm::InputTag passAbsEtaMaxTag_;
  edm::InputTag passECALIsoMaxTag_;
  edm::InputTag passHCALIsoMaxTag_;
  edm::InputTag passHOverEMaxTag_;
  edm::InputTag passTrackIsoMaxTag_;
  edm::InputTag passSigmaIetaIetaMaxTag_;
  edm::InputTag passAbsSeedTimeMaxTag_;
  edm::InputTag passE2OverE9MaxTag_;

  //triggers
  STRING HLTProcessName_;
  HLTConfigProvider HLTCfg_;
  VSTRING HLTPathGG_;
  VSTRING HLTPathFF_;
  VSTRING HLTPathAllPhotons_;

  //MET plots (1 for each trigger, filled once per event)
  VTH1F ggMET_;
  VTH1F ffMET_;
  VTH1F eeMET_;
  VTH1F egMET_;

  //di-EM ET plots (1 for each trigger, filled once per event)
  VTH1F ggDiEMET_;
  VTH1F ffDiEMET_;
  VTH1F eeDiEMET_;
  VTH1F egDiEMET_;

  //lead and trail photon plots (1 for each trigger, filled once per event)
  VTH1F ggPhoton1ET_;
  VTH1F ggPhoton1Eta_;
  VTH1F ggPhoton1Phi_;
  VTH1F ggPhoton1ECALIso_;
  VTH1F ggPhoton1HCALIso_;
  VTH1F ggPhoton1HOverE_;
  VTH1F ggPhoton1TrackIso_;
  VTH1F ggPhoton1SigmaIetaIeta_;
  VTH1F ggPhoton1SeedTime_;
  VTH1F ggPhoton1E2OverE9_;
  VTH1F ggPhoton2ET_;
  VTH1F ggPhoton2Eta_;
  VTH1F ggPhoton2Phi_;
  VTH1F ggPhoton2ECALIso_;
  VTH1F ggPhoton2HCALIso_;
  VTH1F ggPhoton2HOverE_;
  VTH1F ggPhoton2TrackIso_;
  VTH1F ggPhoton2SigmaIetaIeta_;
  VTH1F ggPhoton2SeedTime_;
  VTH1F ggPhoton2E2OverE9_;
  VTH1F egPhoton1ET_;
  VTH1F egPhoton1Eta_;
  VTH1F egPhoton1Phi_;
  VTH1F egPhoton1ECALIso_;
  VTH1F egPhoton1HCALIso_;
  VTH1F egPhoton1HOverE_;
  VTH1F egPhoton1TrackIso_;
  VTH1F egPhoton1SigmaIetaIeta_;
  VTH1F egPhoton1SeedTime_;
  VTH1F egPhoton1E2OverE9_;
  VTH1F egPhoton2ET_;
  VTH1F egPhoton2Eta_;
  VTH1F egPhoton2Phi_;
  VTH1F egPhoton2ECALIso_;
  VTH1F egPhoton2HCALIso_;
  VTH1F egPhoton2HOverE_;
  VTH1F egPhoton2TrackIso_;
  VTH1F egPhoton2SigmaIetaIeta_;
  VTH1F egPhoton2SeedTime_;
  VTH1F egPhoton2E2OverE9_;
  VTH1F eePhoton1ET_;
  VTH1F eePhoton1Eta_;
  VTH1F eePhoton1Phi_;
  VTH1F eePhoton1ECALIso_;
  VTH1F eePhoton1HCALIso_;
  VTH1F eePhoton1HOverE_;
  VTH1F eePhoton1TrackIso_;
  VTH1F eePhoton1SigmaIetaIeta_;
  VTH1F eePhoton1SeedTime_;
  VTH1F eePhoton1E2OverE9_;
  VTH1F eePhoton2ET_;
  VTH1F eePhoton2Eta_;
  VTH1F eePhoton2Phi_;
  VTH1F eePhoton2ECALIso_;
  VTH1F eePhoton2HCALIso_;
  VTH1F eePhoton2HOverE_;
  VTH1F eePhoton2TrackIso_;
  VTH1F eePhoton2SigmaIetaIeta_;
  VTH1F eePhoton2SeedTime_;
  VTH1F eePhoton2E2OverE9_;
  VTH1F ffPhoton1ET_;
  VTH1F ffPhoton1Eta_;
  VTH1F ffPhoton1Phi_;
  VTH1F ffPhoton1ECALIso_;
  VTH1F ffPhoton1HCALIso_;
  VTH1F ffPhoton1HOverE_;
  VTH1F ffPhoton1TrackIso_;
  VTH1F ffPhoton1SigmaIetaIeta_;
  VTH1F ffPhoton1SeedTime_;
  VTH1F ffPhoton1E2OverE9_;
  VTH1F ffPhoton2ET_;
  VTH1F ffPhoton2Eta_;
  VTH1F ffPhoton2Phi_;
  VTH1F ffPhoton2ECALIso_;
  VTH1F ffPhoton2HCALIso_;
  VTH1F ffPhoton2HOverE_;
  VTH1F ffPhoton2TrackIso_;
  VTH1F ffPhoton2SigmaIetaIeta_;
  VTH1F ffPhoton2SeedTime_;
  VTH1F ffPhoton2E2OverE9_;

  //N-1 plots (1 for each trigger)
  VTH1F ETNMinus1_;
  VTH1F etaNMinus1_;
  VTH1F phiNMinus1_;
  VTH1F ECALIsoNMinus1_;
  VTH1F HCALIsoNMinus1_;
  VTH1F HOverENMinus1_;
  VTH1F trackIsoNMinus1_;
  VTH1F sigmaIetaIetaNMinus1_;
  VTH1F seedTimeNMinus1_;
  VTH1F e2OverE9NMinus1_;
  VTH1F ECALIsoNMinus1PreselectionAll_;
  VTH1F HCALIsoNMinus1PreselectionAll_;
  VTH1F ECALIsoNMinus1PreselectionDoublePhotonEvts_;
  VTH1F HCALIsoNMinus1PreselectionDoublePhotonEvts_;
  VTH2F seedTimeVsENMinus1_;

  //N-2 plots (1 for each trigger)
  VTH2F ECALIsoVsETNMinus1PreselectionAll_;
  VTH2F HCALIsoVsETNMinus1PreselectionAll_;
  VTH2F ECALIsoVsETNMinus1PreselectionDoublePhotonEvts_;
  VTH2F HCALIsoVsETNMinus1PreselectionDoublePhotonEvts_;

  //efficiency w.r.t acceptance (ET AND |eta| AND timing), all photons (1 for each trigger)
  VTH1F numPassingECALIsoETAbsEtaVsNPVAll_;
  VTGRAPHASYMMERRORS ECALIsoEffWRTAcceptanceVsNPVAll_;
  VTH1F numPassingHCALIsoETAbsEtaVsNPVAll_;
  VTGRAPHASYMMERRORS HCALIsoEffWRTAcceptanceVsNPVAll_;
  VTH1F numPassingHOverEETAbsEtaVsNPVAll_;
  VTGRAPHASYMMERRORS HOverEEffWRTAcceptanceVsNPVAll_;
  VTH1F numPassingPreselectionVsNPVAll_;
  VTGRAPHASYMMERRORS preselectionEffWRTAcceptanceVsNPVAll_;
  VTH1F numPassingETAbsEtaVsNPVAll_;

  /*efficiency w.r.t acceptance (ET AND |eta| AND timing), only photons in 2-photon events (1 for 
    each trigger)*/
  VTH1F numPassingECALIsoETAbsEtaVsNPVDoublePhotonEvts_;
  VTGRAPHASYMMERRORS ECALIsoEffWRTAcceptanceVsNPVDoublePhotonEvts_;
  VTH1F numPassingHCALIsoETAbsEtaVsNPVDoublePhotonEvts_;
  VTGRAPHASYMMERRORS HCALIsoEffWRTAcceptanceVsNPVDoublePhotonEvts_;
  VTH1F numPassingHOverEETAbsEtaVsNPVDoublePhotonEvts_;
  VTGRAPHASYMMERRORS HOverEEffWRTAcceptanceVsNPVDoublePhotonEvts_;
  VTH1F numPassingPreselectionVsNPVDoublePhotonEvts_;
  VTGRAPHASYMMERRORS preselectionEffWRTAcceptanceVsNPVDoublePhotonEvts_;
  VTH1F numPassingETAbsEtaVsNPVDoublePhotonEvts_;

  //efficiency w.r.t acceptance (ET AND |eta| AND timing) AND H/E, all photons (1 for each trigger)
  VTH1F numPassingECALIsoETAbsEtaHOverEVsNPVAll_;
  VTGRAPHASYMMERRORS ECALIsoEffWRTETAbsEtaHOverEVsNPVAll_;
  VTH1F numPassingHCALIsoETAbsEtaHOverEVsNPVAll_;
  VTGRAPHASYMMERRORS HCALIsoEffWRTETAbsEtaHOverEVsNPVAll_;
  VTH1F numPassingETAbsEtaHOverEVsNPVAll_;

  /*efficiency w.r.t acceptance (ET AND |eta| AND timing) AND H/E, only photons in 2-photon events 
    (1 for each trigger)*/
  VTH1F numPassingECALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_;
  VTGRAPHASYMMERRORS ECALIsoEffWRTETAbsEtaHOverEVsNPVDoublePhotonEvts_;
  VTH1F numPassingHCALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_;
  VTGRAPHASYMMERRORS HCALIsoEffWRTETAbsEtaHOverEVsNPVDoublePhotonEvts_;
  VTH1F numPassingETAbsEtaHOverEVsNPVDoublePhotonEvts_;

  //counters
  unsigned int numEvts_;
  unsigned int numGG_;
  unsigned int numEG_;
  unsigned int numFF_;
  unsigned int numEE_;

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
GMSBAnalyzer::GMSBAnalyzer(const edm::ParameterSet& iConfig) :

  //input
  outputFile_(iConfig.getUntrackedParameter<STRING>("outputFile", "skim_analysis.root")),
  photonTag_(iConfig.getUntrackedParameter<edm::InputTag>
	     ("photonTag", edm::InputTag("photons", "", "RECOCleaned"))),
  categoryTag_(iConfig.getUntrackedParameter<edm::InputTag>
	       ("categoryTag", edm::InputTag("EDMCategoryProducer", "eventCategory", 
					     "CATEGORY"))),
  diEMETTag_(iConfig.getUntrackedParameter<edm::InputTag>
	     ("diEMETTag", edm::InputTag("EDMCategoryProducer", "evtDiEMET", "CATEGORY"))),
  photonSeedTimeTag_(iConfig.getUntrackedParameter<edm::InputTag>
		     ("photonSeedTimeTag", 
		      edm::InputTag("EDMCategoryProducer", "photonSeedTime", "OWNPARTICLES"))),
  photonE2OverE9Tag_(iConfig.getUntrackedParameter<edm::InputTag>
		     ("photonE2OverE9Tag", 
		      edm::InputTag("EDMCategoryProducer", "photonE2OverE9", "OWNPARTICLES"))),
  passAllTag_(iConfig.getUntrackedParameter<edm::InputTag>
	      ("passAllTag", edm::InputTag("EDMCategoryProducer", "passingPhotons", "CATEGORY"))),
  trgResultsTag_(iConfig.getUntrackedParameter<edm::InputTag>
		 ("trgResultsTag", edm::InputTag("TriggerResults", "", "HLT"))),
  tcMETTag_(iConfig.getUntrackedParameter<edm::InputTag>
	    ("tcMETTag", edm::InputTag("tcMet", "", "RECO"))),
  vertexTag_(iConfig.getUntrackedParameter<edm::InputTag>
	    ("vertexTag", edm::InputTag("offlinePrimaryVertices", "", "RECO"))),
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
  passAbsSeedTimeMaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
			 ("passAbsSeedTimeMaxTag", edm::InputTag("EDMCategoryProducer", 
								 "passAbsSeedTimeMax", 
								 "OWNPARTICLES"))),
  passE2OverE9MaxTag_(iConfig.getUntrackedParameter<edm::InputTag>
		      ("passE2OverE9MaxTag", edm::InputTag("EDMCategoryProducer", 
		      "passE2OverE9Max", "OWNPARTICLES"))),

  //triggers
  HLTProcessName_(iConfig.getUntrackedParameter<STRING>("HLTProcessName", "HLT")),
  HLTPathGG_(iConfig.getUntrackedParameter<VSTRING>
	     ("HLTPathGG", VSTRING(1, "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1"))),
  HLTPathFF_(iConfig.getUntrackedParameter<VSTRING>
	     ("HLTPathFF", VSTRING(1, "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1"))),
  HLTPathAllPhotons_(iConfig.getUntrackedParameter<VSTRING>
		     ("HLTPathAllPhotons", VSTRING(1, "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1")))
{
}

GMSBAnalyzer::~GMSBAnalyzer()
{}


//
// member functions
//

// ------------ method called on each new Event  ------------
void
GMSBAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //increment event counter
  ++numEvts_;

  //get a bunch of needed products
  edm::Handle<reco::PhotonCollection> pPhotons;
  edm::Handle<int> pCategory;
  edm::Handle<double> pDiEMET;
  edm::Handle<edm::ValueMap<double> > pPhotonSeedTime;
  edm::Handle<edm::ValueMap<double> > pPhotonE2OverE9;
  edm::Handle<edm::ValueMap<bool> > pPassAll;
  edm::Handle<edm::TriggerResults> pTrgResults;
  edm::Handle<reco::METCollection> pTCMET;
  edm::Handle<reco::VertexCollection> pVertices;
  edm::Handle<edm::ValueMap<bool> > pPassETMin;
  edm::Handle<edm::ValueMap<bool> > pPassAbsEtaMax;
  edm::Handle<edm::ValueMap<bool> > pPassECALIsoMax;
  edm::Handle<edm::ValueMap<bool> > pPassHCALIsoMax;
  edm::Handle<edm::ValueMap<bool> > pPassHOverEMax;
  edm::Handle<edm::ValueMap<bool> > pPassTrackIsoMax;
  edm::Handle<edm::ValueMap<bool> > pPassSigmaIetaIetaMax;
  edm::Handle<edm::ValueMap<bool> > pPassAbsSeedTimeMax;
  edm::Handle<edm::ValueMap<bool> > pPassE2OverE9Max;
  if (getCollection_(pPhotons, photonTag_, iEvent) && 
      getCollection_(pCategory, categoryTag_, iEvent) && 
      getCollection_(pDiEMET, diEMETTag_, iEvent) && 
      getCollection_(pPhotonSeedTime, photonSeedTimeTag_, iEvent) && 
      getCollection_(pPhotonE2OverE9, photonE2OverE9Tag_, iEvent) && 
      getCollection_(pPassAll, passAllTag_, iEvent) && 
      getCollection_(pTrgResults, trgResultsTag_, iEvent) && 
      getCollection_(pTCMET, tcMETTag_, iEvent) &&  
      getCollection_(pVertices, vertexTag_, iEvent) &&  
      getCollection_(pPassETMin, passETMinTag_, iEvent) && 
      getCollection_(pPassAbsEtaMax, passAbsEtaMaxTag_, iEvent) && 
      getCollection_(pPassECALIsoMax, passECALIsoMaxTag_, iEvent) && 
      getCollection_(pPassHCALIsoMax, passHCALIsoMaxTag_, iEvent) && 
      getCollection_(pPassHOverEMax, passHOverEMaxTag_, iEvent) && 
      getCollection_(pPassTrackIsoMax, passTrackIsoMaxTag_, iEvent) && 
      getCollection_(pPassSigmaIetaIetaMax, passSigmaIetaIetaMaxTag_, iEvent) && 
      getCollection_(pPassAbsSeedTimeMax, passAbsSeedTimeMaxTag_, iEvent) && 
      getCollection_(pPassE2OverE9Max, passE2OverE9MaxTag_, iEvent)) {

    //pointers to triggers and plots relevant to this event
    VSTRING* HLTPaths = NULL;
    VTH1F* MET = NULL;
    VTH1F* diEMET = NULL;
    VTH1F* photon1ET = NULL;
    VTH1F* photon1Eta = NULL;
    VTH1F* photon1Phi = NULL;
    VTH1F* photon1ECALIso = NULL;
    VTH1F* photon1HCALIso = NULL;
    VTH1F* photon1HOverE = NULL;
    VTH1F* photon1TrackIso = NULL;
    VTH1F* photon1SigmaIetaIeta = NULL;
    VTH1F* photon1SeedTime = NULL;
    VTH1F* photon1E2OverE9 = NULL;
    VTH1F* photon2ET = NULL;
    VTH1F* photon2Eta = NULL;
    VTH1F* photon2Phi = NULL;
    VTH1F* photon2ECALIso = NULL;
    VTH1F* photon2HCALIso = NULL;
    VTH1F* photon2HOverE = NULL;
    VTH1F* photon2TrackIso = NULL;
    VTH1F* photon2SigmaIetaIeta = NULL;
    VTH1F* photon2SeedTime = NULL;
    VTH1F* photon2E2OverE9 = NULL;

    //vector of pointers to the value maps for the straight N-1 plot decisions
    std::vector<edm::Handle<edm::ValueMap<bool> >*> passMapVec;
    passMapVec.push_back(&pPassETMin);
    passMapVec.push_back(&pPassAbsEtaMax);
    passMapVec.push_back(&pPassECALIsoMax);
    passMapVec.push_back(&pPassHCALIsoMax);
    passMapVec.push_back(&pPassHOverEMax);
    passMapVec.push_back(&pPassTrackIsoMax);
    passMapVec.push_back(&pPassSigmaIetaIetaMax);
    passMapVec.push_back(&pPassAbsSeedTimeMax);
    passMapVec.push_back(&pPassE2OverE9Max);

    //vector of pointers to the value maps for the preselection only N-1 plot decisions
    std::vector<edm::Handle<edm::ValueMap<bool> >*> passMapVecPreselection;
    passMapVecPreselection.push_back(&pPassETMin);
    passMapVecPreselection.push_back(&pPassAbsEtaMax);
    passMapVecPreselection.push_back(&pPassECALIsoMax);
    passMapVecPreselection.push_back(&pPassHCALIsoMax);
    passMapVecPreselection.push_back(&pPassHOverEMax);
    passMapVecPreselection.push_back(&pPassAbsSeedTimeMax);

    //vector of pointers to the value maps for the preselection only N-2 plot decisions
    std::vector<edm::Handle<edm::ValueMap<bool> >*> passMapVecPreselectionNMinus2;
    passMapVecPreselectionNMinus2.push_back(&pPassAbsEtaMax);
    passMapVecPreselectionNMinus2.push_back(&pPassECALIsoMax);
    passMapVecPreselectionNMinus2.push_back(&pPassHCALIsoMax);
    passMapVecPreselectionNMinus2.push_back(&pPassHOverEMax);

    /*vector of pointers to the value maps for the ECAL isolation efficiency w.r.t. acceptance (ET 
      and |eta| and timing) plots*/
    std::vector<edm::Handle<edm::ValueMap<bool> >*> passMapVecECALIsoEffWRTAcceptance;
    passMapVecECALIsoEffWRTAcceptance.push_back(&pPassETMin);
    passMapVecECALIsoEffWRTAcceptance.push_back(&pPassAbsEtaMax);
    passMapVecECALIsoEffWRTAcceptance.push_back(&pPassECALIsoMax);
    passMapVecECALIsoEffWRTAcceptance.push_back(&pPassAbsSeedTimeMax);

    /*vector of pointers to the value maps for the HCAL isolation efficiency w.r.t. acceptance (ET 
      and |eta| and timing) plots*/
    std::vector<edm::Handle<edm::ValueMap<bool> >*> passMapVecHCALIsoEffWRTAcceptance;
    passMapVecHCALIsoEffWRTAcceptance.push_back(&pPassETMin);
    passMapVecHCALIsoEffWRTAcceptance.push_back(&pPassAbsEtaMax);
    passMapVecHCALIsoEffWRTAcceptance.push_back(&pPassHCALIsoMax);
    passMapVecHCALIsoEffWRTAcceptance.push_back(&pPassAbsSeedTimeMax);

    /*vector of pointers to the value maps for the H/E efficiency w.r.t. acceptance (ET and |eta| 
      and timing) plots*/
    std::vector<edm::Handle<edm::ValueMap<bool> >*> passMapVecHOverEEffWRTAcceptance;
    passMapVecHOverEEffWRTAcceptance.push_back(&pPassETMin);
    passMapVecHOverEEffWRTAcceptance.push_back(&pPassAbsEtaMax);
    passMapVecHOverEEffWRTAcceptance.push_back(&pPassHOverEMax);
    passMapVecHOverEEffWRTAcceptance.push_back(&pPassAbsSeedTimeMax);

    /*vector of pointers to the value maps for the isolation efficiency w.r.t. preselection 
      (ET, |eta|, and H/E) plots*/
    std::vector<edm::Handle<edm::ValueMap<bool> >*> passMapVecIsoEffWRTPreselection;
    passMapVecIsoEffWRTPreselection.push_back(&pPassETMin);
    passMapVecIsoEffWRTPreselection.push_back(&pPassAbsEtaMax);
    passMapVecIsoEffWRTPreselection.push_back(&pPassHOverEMax);
    passMapVecIsoEffWRTPreselection.push_back(&pPassAbsSeedTimeMax);

    //get the correct plots and triggers for event category
    switch (*(pCategory.product())) {
    case FAIL:
      HLTPaths = &HLTPathAllPhotons_;
      break;
    case GG:
      HLTPaths = &HLTPathGG_;
      MET = &ggMET_;
      diEMET = &ggDiEMET_;
      photon1ET = &ggPhoton1ET_;
      photon1Eta = &ggPhoton1Eta_;
      photon1Phi = &ggPhoton1Phi_;
      photon1ECALIso = &ggPhoton1ECALIso_;
      photon1HCALIso = &ggPhoton1HCALIso_;
      photon1HOverE = &ggPhoton1HOverE_;
      photon1TrackIso = &ggPhoton1TrackIso_;
      photon1SigmaIetaIeta = &ggPhoton1SigmaIetaIeta_;
      photon1SeedTime = &ggPhoton1SeedTime_;
      photon1E2OverE9 = &ggPhoton1E2OverE9_;
      photon2ET = &ggPhoton2ET_;
      photon2Eta = &ggPhoton2Eta_;
      photon2Phi = &ggPhoton2Phi_;
      photon2ECALIso = &ggPhoton2ECALIso_;
      photon2HCALIso = &ggPhoton2HCALIso_;
      photon2HOverE = &ggPhoton2HOverE_;
      photon2TrackIso = &ggPhoton2TrackIso_;
      photon2SigmaIetaIeta = &ggPhoton2SigmaIetaIeta_;
      photon2SeedTime = &ggPhoton2SeedTime_;
      photon2E2OverE9 = &ggPhoton2E2OverE9_;
      break;
    case FF:
      HLTPaths = &HLTPathFF_;
      MET = &ffMET_;
      diEMET = &ffDiEMET_;
      photon1ET = &ffPhoton1ET_;
      photon1Eta = &ffPhoton1Eta_;
      photon1Phi = &ffPhoton1Phi_;
      photon1ECALIso = &ffPhoton1ECALIso_;
      photon1HCALIso = &ffPhoton1HCALIso_;
      photon1HOverE = &ffPhoton1HOverE_;
      photon1TrackIso = &ffPhoton1TrackIso_;
      photon1SigmaIetaIeta = &ffPhoton1SigmaIetaIeta_;
      photon1SeedTime = &ffPhoton1SeedTime_;
      photon1E2OverE9 = &ffPhoton1E2OverE9_;
      photon2ET = &ffPhoton2ET_;
      photon2Eta = &ffPhoton2Eta_;
      photon2Phi = &ffPhoton2Phi_;
      photon2ECALIso = &ffPhoton2ECALIso_;
      photon2HCALIso = &ffPhoton2HCALIso_;
      photon2HOverE = &ffPhoton2HOverE_;
      photon2TrackIso = &ffPhoton2TrackIso_;
      photon2SigmaIetaIeta = &ffPhoton2SigmaIetaIeta_;
      photon2SeedTime = &ffPhoton2SeedTime_;
      photon2E2OverE9 = &ffPhoton2E2OverE9_;
      break;
    case EG:
      HLTPaths = &HLTPathGG_;
      MET = &egMET_;
      diEMET = &egDiEMET_;
      photon1ET = &egPhoton1ET_;
      photon1Eta = &egPhoton1Eta_;
      photon1Phi = &egPhoton1Phi_;
      photon1ECALIso = &egPhoton1ECALIso_;
      photon1HCALIso = &egPhoton1HCALIso_;
      photon1HOverE = &egPhoton1HOverE_;
      photon1TrackIso = &egPhoton1TrackIso_;
      photon1SigmaIetaIeta = &egPhoton1SigmaIetaIeta_;
      photon1SeedTime = &egPhoton1SeedTime_;
      photon1E2OverE9 = &egPhoton1E2OverE9_;
      photon2ET = &egPhoton2ET_;
      photon2Eta = &egPhoton2Eta_;
      photon2Phi = &egPhoton2Phi_;
      photon2ECALIso = &egPhoton2ECALIso_;
      photon2HCALIso = &egPhoton2HCALIso_;
      photon2HOverE = &egPhoton2HOverE_;
      photon2TrackIso = &egPhoton2TrackIso_;
      photon2SigmaIetaIeta = &egPhoton2SigmaIetaIeta_;
      photon2SeedTime = &egPhoton2SeedTime_;
      photon2E2OverE9 = &egPhoton2E2OverE9_;
      break;
    case EE:
      HLTPaths = &HLTPathGG_;
      MET = &eeMET_;
      diEMET = &eeDiEMET_;
      photon1ET = &eePhoton1ET_;
      photon1Eta = &eePhoton1Eta_;
      photon1Phi = &eePhoton1Phi_;
      photon1ECALIso = &eePhoton1ECALIso_;
      photon1HCALIso = &eePhoton1HCALIso_;
      photon1HOverE = &eePhoton1HOverE_;
      photon1TrackIso = &eePhoton1TrackIso_;
      photon1SigmaIetaIeta = &eePhoton1SigmaIetaIeta_;
      photon1SeedTime = &eePhoton1SeedTime_;
      photon1E2OverE9 = &eePhoton1E2OverE9_;
      photon2ET = &eePhoton2ET_;
      photon2Eta = &eePhoton2Eta_;
      photon2Phi = &eePhoton2Phi_;
      photon2ECALIso = &eePhoton2ECALIso_;
      photon2HCALIso = &eePhoton2HCALIso_;
      photon2HOverE = &eePhoton2HOverE_;
      photon2TrackIso = &eePhoton2TrackIso_;
      photon2SigmaIetaIeta = &eePhoton2SigmaIetaIeta_;
      photon2SeedTime = &eePhoton2SeedTime_;
      photon2E2OverE9 = &eePhoton2E2OverE9_;
      break;
    default:
      break;
    }

    //get the real trigger names from the approximate inputs
    VSTRING realHLT;
    VSTRING realHLTAllPhotons;
    try {
      realHLT = getRealTriggerNames(HLTPaths);
      realHLTAllPhotons = getRealTriggerNames(&HLTPathAllPhotons_);
    }
    catch (cms::Exception& ex) {
      STRINGSTREAM err;
      err << "Category: " << *(pCategory.product()) << std::endl;
      throw cms::Exception("GMSBAnalyzer") << ex.what() << err.str();
    }

    //fill the vectors of pointers to vectors of 1D histograms of photon quantities
    std::vector<VTH1F*> photon1DHistPtrVec;
    photon1DHistPtrVec.push_back(&ETNMinus1_);
    photon1DHistPtrVec.push_back(&etaNMinus1_);
    photon1DHistPtrVec.push_back(&phiNMinus1_);
    photon1DHistPtrVec.push_back(&ECALIsoNMinus1_);
    photon1DHistPtrVec.push_back(&HCALIsoNMinus1_);
    photon1DHistPtrVec.push_back(&HOverENMinus1_);
    photon1DHistPtrVec.push_back(&trackIsoNMinus1_);
    photon1DHistPtrVec.push_back(&sigmaIetaIetaNMinus1_);
    photon1DHistPtrVec.push_back(&seedTimeNMinus1_);
    photon1DHistPtrVec.push_back(&e2OverE9NMinus1_);
    photon1DHistPtrVec.push_back(&ECALIsoNMinus1PreselectionAll_);
    photon1DHistPtrVec.push_back(&HCALIsoNMinus1PreselectionAll_);
    photon1DHistPtrVec.push_back(&ECALIsoNMinus1PreselectionDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&HCALIsoNMinus1PreselectionDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&numPassingECALIsoETAbsEtaVsNPVAll_);
    photon1DHistPtrVec.push_back(&numPassingHCALIsoETAbsEtaVsNPVAll_);
    photon1DHistPtrVec.push_back(&numPassingHOverEETAbsEtaVsNPVAll_);
    photon1DHistPtrVec.push_back(&numPassingPreselectionVsNPVAll_);
    photon1DHistPtrVec.push_back(&numPassingETAbsEtaVsNPVAll_);
    photon1DHistPtrVec.push_back(&numPassingECALIsoETAbsEtaVsNPVDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&numPassingHCALIsoETAbsEtaVsNPVDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&numPassingHOverEETAbsEtaVsNPVDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&numPassingPreselectionVsNPVDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&numPassingETAbsEtaVsNPVDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&numPassingECALIsoETAbsEtaHOverEVsNPVAll_);
    photon1DHistPtrVec.push_back(&numPassingHCALIsoETAbsEtaHOverEVsNPVAll_);
    photon1DHistPtrVec.push_back(&numPassingETAbsEtaHOverEVsNPVAll_);
    photon1DHistPtrVec.push_back(&numPassingECALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&numPassingHCALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_);
    photon1DHistPtrVec.push_back(&numPassingETAbsEtaHOverEVsNPVDoublePhotonEvts_);

    //fill the vectors of pointers to vectors of 2D histograms of photon quantities
    std::vector<VTH2F*> photon2DHistPtrVec;
    photon2DHistPtrVec.push_back(&seedTimeVsENMinus1_);
    photon2DHistPtrVec.push_back(&ECALIsoVsETNMinus1PreselectionAll_);
    photon2DHistPtrVec.push_back(&HCALIsoVsETNMinus1PreselectionAll_);
    photon2DHistPtrVec.push_back(&ECALIsoVsETNMinus1PreselectionDoublePhotonEvts_);
    photon2DHistPtrVec.push_back(&HCALIsoVsETNMinus1PreselectionDoublePhotonEvts_);

    //do sample-based processing
    if ((*(pCategory.product()) == GG) || (*(pCategory.product()) == EG) || 
	(*(pCategory.product()) == EE) || (*(pCategory.product()) == FF)) {

      //fill the vectors of pointers to vectors of histograms of event quantities
      std::vector<VTH1F*> evtHistPtrVec;
      evtHistPtrVec.push_back(MET);
      evtHistPtrVec.push_back(diEMET);
      evtHistPtrVec.push_back(photon1ET);
      evtHistPtrVec.push_back(photon1Eta);
      evtHistPtrVec.push_back(photon1Phi);
      evtHistPtrVec.push_back(photon1ECALIso);
      evtHistPtrVec.push_back(photon1HCALIso);
      evtHistPtrVec.push_back(photon1HOverE);
      evtHistPtrVec.push_back(photon1TrackIso);
      evtHistPtrVec.push_back(photon1SigmaIetaIeta);
      evtHistPtrVec.push_back(photon1SeedTime);
      evtHistPtrVec.push_back(photon1E2OverE9);
      evtHistPtrVec.push_back(photon2ET);
      evtHistPtrVec.push_back(photon2Eta);
      evtHistPtrVec.push_back(photon2Phi);
      evtHistPtrVec.push_back(photon2ECALIso);
      evtHistPtrVec.push_back(photon2HCALIso);
      evtHistPtrVec.push_back(photon2HOverE);
      evtHistPtrVec.push_back(photon2TrackIso);
      evtHistPtrVec.push_back(photon2SigmaIetaIeta);
      evtHistPtrVec.push_back(photon2SeedTime);
      evtHistPtrVec.push_back(photon2E2OverE9);

      //fill the vector of event quantities
      VFLOAT evtQuantityVec;
      evtQuantityVec.push_back((float)(pTCMET->begin()->et()));
      evtQuantityVec.push_back((float)(*(pDiEMET.product())));
      reco::Photon photon1;
      reco::Photon photon2;
      unsigned int iPhoton1 = pPhotons->size();
      unsigned int iPhoton2 = pPhotons->size();
      getPassingPhotons(pPhotons, photon1, photon2, iPhoton1, iPhoton2, pPassAll);
      evtQuantityVec.push_back((float)photon1.et());
      evtQuantityVec.push_back((float)photon1.eta());
      evtQuantityVec.push_back((float)photon1.phi());
      evtQuantityVec.push_back(photon1.ecalRecHitSumEtConeDR04());
      evtQuantityVec.push_back(photon1.hcalTowerSumEtConeDR04());
      evtQuantityVec.push_back(photon1.hadronicOverEm());
      evtQuantityVec.push_back(photon1.trkSumPtHollowConeDR04());
      evtQuantityVec.push_back(photon1.sigmaIetaIeta());
      edm::Ref<reco::PhotonCollection> ref1(pPhotons, iPhoton1);
      evtQuantityVec.push_back((*(pPhotonSeedTime.product()))[ref1]);
      evtQuantityVec.push_back((*(pPhotonE2OverE9.product()))[ref1]);
      evtQuantityVec.push_back((float)photon2.et());
      evtQuantityVec.push_back((float)photon2.eta());
      evtQuantityVec.push_back((float)photon2.phi());
      evtQuantityVec.push_back(photon2.ecalRecHitSumEtConeDR04());
      evtQuantityVec.push_back(photon2.hcalTowerSumEtConeDR04());
      evtQuantityVec.push_back(photon2.hadronicOverEm());
      evtQuantityVec.push_back(photon2.trkSumPtHollowConeDR04());
      evtQuantityVec.push_back(photon2.sigmaIetaIeta());
      edm::Ref<reco::PhotonCollection> ref2(pPhotons, iPhoton2);
      evtQuantityVec.push_back((*(pPhotonSeedTime.product()))[ref2]);
      evtQuantityVec.push_back((*(pPhotonE2OverE9.product()))[ref2]);

      /*proceed if the event passed the trigger (no matching of offline reco object to online 
	trigger object)*/
      const edm::TriggerNames& trgNames = iEvent.triggerNames(*pTrgResults);
      for (VSTRING_IT iHLT = realHLT.begin(); iHLT != realHLT.end(); ++iHLT) {
	const unsigned int trgIndex = trgNames.triggerIndex(*iHLT);
	if ((trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex))) {

	  //fill histograms of event quantities
	  fillHistograms(evtHistPtrVec, evtQuantityVec, VBOOL(evtHistPtrVec.size(), true), 
			 HLTPaths->at(iHLT - realHLT.begin()));

	  //increment sample counters
	  if (*(pCategory.product()) == GG) ++numGG_;
	  if (*(pCategory.product()) == EG) ++numEG_;
	  if (*(pCategory.product()) == FF) ++numFF_;
	  if (*(pCategory.product()) == EE) ++numEE_;
	}
      }
    }

    //do photon-based processing
    for (reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin(); 
	 iPhoton != pPhotons->end(); ++iPhoton) {

      /*proceed if the event passed the trigger (no matching of offline reco object to online 
	trigger object)*/
      const edm::TriggerNames& trgNames = iEvent.triggerNames(*pTrgResults);
      for (VSTRING_IT iHLT = realHLTAllPhotons.begin(); iHLT != realHLTAllPhotons.end(); 
	   ++iHLT) {
	const unsigned int trgIndex = trgNames.triggerIndex(*iHLT);
	if ((trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex))) {

	  //fill the vector of photon quantities for 1D histograms
	  VFLOAT photonQuantityVec;
	  photonQuantityVec.push_back((float)(iPhoton->et()));
	  photonQuantityVec.push_back((float)(iPhoton->eta()));
	  photonQuantityVec.push_back((float)(iPhoton->phi()));
	  photonQuantityVec.push_back((float)(iPhoton->ecalRecHitSumEtConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->hcalTowerSumEtConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->hadronicOverEm()));
	  photonQuantityVec.push_back((float)(iPhoton->trkSumPtHollowConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->sigmaIetaIeta()));
	  edm::Ref<reco::PhotonCollection> ref(pPhotons, iPhoton - pPhotons->begin());
	  photonQuantityVec.push_back((*(pPhotonSeedTime.product()))[ref]);
	  photonQuantityVec.push_back((*(pPhotonE2OverE9.product()))[ref]);
	  photonQuantityVec.push_back((float)(iPhoton->ecalRecHitSumEtConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->hcalTowerSumEtConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->ecalRecHitSumEtConeDR04()));
	  photonQuantityVec.push_back((float)(iPhoton->hcalTowerSumEtConeDR04()));
	  unsigned int numPassingPVs = numValidPVs(pVertices);
	  photonQuantityVec.push_back((float)effDecision(passMapVecECALIsoEffWRTAcceptance, ref, 
							 NULL, numPassingPVs));
	  photonQuantityVec.push_back((float)effDecision(passMapVecHCALIsoEffWRTAcceptance, ref, 
							 NULL, numPassingPVs));
	  photonQuantityVec.push_back((float)effDecision(passMapVecHOverEEffWRTAcceptance, ref, 
							 NULL, numPassingPVs));
	  photonQuantityVec.push_back((float)effDecision(passMapVecPreselection, ref, NULL, 
							 numPassingPVs));
	  photonQuantityVec.push_back((float)effDecision(passMapVecHCALIsoEffWRTAcceptance, ref, 
							 &pPassHCALIsoMax, numPassingPVs));
	  photonQuantityVec.push_back((float)effDecision(passMapVecECALIsoEffWRTAcceptance, ref, 
							 NULL, numPassingPVs, pPhotons));
	  photonQuantityVec.push_back((float)effDecision(passMapVecHCALIsoEffWRTAcceptance, ref, 
							 NULL, numPassingPVs, pPhotons));
	  photonQuantityVec.push_back((float)effDecision(passMapVecHOverEEffWRTAcceptance, ref, 
							 NULL, numPassingPVs, pPhotons));
	  photonQuantityVec.push_back((float)effDecision(passMapVecPreselection, ref, NULL, 
							 numPassingPVs, pPhotons));
	  photonQuantityVec.push_back((float)effDecision(passMapVecHCALIsoEffWRTAcceptance, ref, 
							 &pPassHCALIsoMax, numPassingPVs, 
							 pPhotons));
	  photonQuantityVec.push_back((float)effDecision(passMapVecPreselection, ref, 
							 &pPassHCALIsoMax, numPassingPVs));
	  photonQuantityVec.push_back((float)effDecision(passMapVecPreselection, ref, 
							 &pPassECALIsoMax, numPassingPVs));
	  photonQuantityVec.push_back((float)effDecision(passMapVecIsoEffWRTPreselection, ref, 
							 NULL, numPassingPVs));
	  photonQuantityVec.push_back((float)effDecision(passMapVecPreselection, ref, 
							 &pPassHCALIsoMax, numPassingPVs, 
							 pPhotons));
	  photonQuantityVec.push_back((float)effDecision(passMapVecPreselection, ref, 
							 &pPassECALIsoMax, numPassingPVs, 
							 pPhotons));
	  photonQuantityVec.push_back((float)effDecision(passMapVecIsoEffWRTPreselection, ref, 
							 NULL, numPassingPVs, pPhotons));

	  //fill the 2 vectors of photon quantities for 2D histograms
	  VFLOAT photonQuantity1Vec;
	  photonQuantity1Vec.push_back(iPhoton->energy());
	  photonQuantity1Vec.push_back(numPassingPVs);
	  photonQuantity1Vec.push_back(numPassingPVs);
	  photonQuantity1Vec.push_back(numPassingPVs);
	  photonQuantity1Vec.push_back(numPassingPVs);
	  photonQuantity1Vec.push_back(iPhoton->et());
	  photonQuantity1Vec.push_back(iPhoton->et());
	  photonQuantity1Vec.push_back(iPhoton->et());
	  photonQuantity1Vec.push_back(iPhoton->et());
	  VFLOAT photonQuantity2Vec;
	  photonQuantity2Vec.push_back((*(pPhotonSeedTime.product()))[ref]);
	  photonQuantity2Vec.push_back((float)(iPhoton->ecalRecHitSumEtConeDR04()));
	  photonQuantity2Vec.push_back((float)(iPhoton->hcalTowerSumEtConeDR04()));
	  photonQuantity2Vec.push_back((float)(iPhoton->ecalRecHitSumEtConeDR04()));
	  photonQuantity2Vec.push_back((float)(iPhoton->hcalTowerSumEtConeDR04()));
	  photonQuantity2Vec.push_back((float)(iPhoton->ecalRecHitSumEtConeDR04()));
	  photonQuantity2Vec.push_back((float)(iPhoton->hcalTowerSumEtConeDR04()));
	  photonQuantity2Vec.push_back((float)(iPhoton->ecalRecHitSumEtConeDR04()));
	  photonQuantity2Vec.push_back((float)(iPhoton->hcalTowerSumEtConeDR04()));

	  /*fill the vector of photon pass flags for 1D plots
	    photon passes if it passes all selection cuts except the one being plotted*/
	  VBOOL passVec1D;
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassETMin));
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassAbsEtaMax));
	  passVec1D.push_back(true);
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassECALIsoMax));
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassHCALIsoMax));
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassHOverEMax));
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassTrackIsoMax));
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassSigmaIetaIetaMax));
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassAbsSeedTimeMax));
	  passVec1D.push_back(passNMinus1(passMapVec, ref, &pPassE2OverE9Max));
	  passVec1D.push_back(passNMinus1(passMapVecPreselection, ref, &pPassECALIsoMax));
	  passVec1D.push_back(passNMinus1(passMapVecPreselection, ref, &pPassHCALIsoMax));
	  passVec1D.push_back(twoPhotonsPassNMinus1(passMapVecPreselection, &pPassECALIsoMax, 
						    pPhotons));
	  passVec1D.push_back(twoPhotonsPassNMinus1(passMapVecPreselection, &pPassHCALIsoMax, 
						    pPhotons));
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);
	  passVec1D.push_back(true);

	  /*fill the vector of photon pass flags flags for 2D plots
	    photon passes if it passes all selection cuts except the one being plotted*/
	  VBOOL passVec2D;
	  passVec2D.push_back(passNMinus1(passMapVec, ref, &pPassAbsSeedTimeMax));
	  passVec2D.push_back(passNMinus1(passMapVecPreselection, ref, &pPassECALIsoMax));
	  passVec2D.push_back(passNMinus1(passMapVecPreselection, ref, &pPassHCALIsoMax));
	  passVec2D.push_back(twoPhotonsPassNMinus1(passMapVecPreselection, &pPassECALIsoMax, 
						    pPhotons));
	  passVec2D.push_back(twoPhotonsPassNMinus1(passMapVecPreselection, &pPassHCALIsoMax, 
						    pPhotons));
	  passVec2D.push_back(passNMinus1(passMapVecPreselectionNMinus2, ref, &pPassECALIsoMax));
	  passVec2D.push_back(passNMinus1(passMapVecPreselectionNMinus2, ref, &pPassHCALIsoMax));
	  passVec2D.push_back(twoPhotonsPassNMinus1(passMapVecPreselectionNMinus2, 
						    &pPassECALIsoMax, pPhotons));
	  passVec2D.push_back(twoPhotonsPassNMinus1(passMapVecPreselectionNMinus2, 
						    &pPassHCALIsoMax, pPhotons));

	  //fill 1D histograms of photon quantities
	  fillHistograms(photon1DHistPtrVec, photonQuantityVec, passVec1D, 
			 HLTPathAllPhotons_[iHLT - realHLTAllPhotons.begin()]);

	  //fill 2D histograms of photon quantities
	  fillHistograms(photon2DHistPtrVec, photonQuantity1Vec, photonQuantity2Vec, passVec2D, 
			 HLTPathAllPhotons_[iHLT - realHLTAllPhotons.begin()]);
	}
      }
    }
  }
}

// ---- method called once each run  ---
void GMSBAnalyzer::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
  bool changed = false;
  if (!HLTCfg_.init(iRun, iSetup, HLTProcessName_, changed) ) {
    edm::LogError("GMSBAnalyzer") << "Error: can't initialize HLTConfigProvider.\n";
    throw cms::Exception("GMSBAnalyzer") << "HLTConfigProvider::init() returned non-0.\n";
  }
  if (changed) edm::LogInfo("GMSBAnalyzer") << "HLT configuration changed!\n";
}

// ------------ method called once each job just before starting event loop  ------------
void 
GMSBAnalyzer::beginJob()
{
  //book histograms (1 for each trigger)
  for (VSTRING_IT iGGHLT = HLTPathGG_.begin(); iGGHLT != HLTPathGG_.end(); ++iGGHLT) {
    bookEAndGHistograms(*iGGHLT);
  }
  for (VSTRING_IT iFFHLT = HLTPathFF_.begin(); iFFHLT != HLTPathFF_.end(); ++iFFHLT) {
    bookFHistograms(*iFFHLT);
  }
  for (VSTRING_IT iAllPhotonsHLT = HLTPathAllPhotons_.begin(); 
       iAllPhotonsHLT != HLTPathAllPhotons_.end(); ++iAllPhotonsHLT) {
    bookPhotonHistograms(*iAllPhotonsHLT);
  }

  //initialize counters
  numEvts_ = 0;
  numGG_ = 0;
  numEG_ = 0;
  numFF_ = 0;
  numEE_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GMSBAnalyzer::endJob() {

  //open output file
  out_ = new TFile(outputFile_.c_str(), "RECREATE");
  if (out_->IsOpen()) {
    out_->cd();

    //write sample plots
    writeGGPlots();
    writeEGPlots();
    writeEEPlots();
    writeFFPlots();

    //write N-1 histograms
    writePhotonPlots();

    //print counters
    std::cout << "Total number of events: " << numEvts_ << std::endl;
    std::cout << "Number of GG events: " << numGG_ << std::endl;
    std::cout << "Number of EG events: " << numEG_ << std::endl;
    std::cout << "Number of FF events: " << numFF_ << std::endl;
    std::cout << "Number of EE events: " << numEE_ << std::endl;

    //write output file
    out_->cd();
    //out_->Write(); //appears unnecessary if Write() is called on each histogram individually
    out_->Close();
  }
  else edm::LogError("GMSBAnalyzer") << "Error opening file " << outputFile_ << ".\n";
  delete out_;
}

STRING GMSBAnalyzer::name(const STRING& category, const STRING& plot, const STRING& HLT) const
{
  STRINGSTREAM name;
  name << category << "_" << plot << "_" << "_" << HLT;
  return name.str();
}

void GMSBAnalyzer::bookHistogram(VTH1F& histVec, const STRING& category, const STRING& plot, 
				 const STRING& HLT, const unsigned int numBins, 
				 const float binLowEdge, const float binHighEdge)
{
  histVec.push_back(new TH1F(name(category, plot, HLT).c_str(), "", numBins, binLowEdge, 
			     binHighEdge));
}

void GMSBAnalyzer::bookHistogram(VTH2F& histVec, const STRING& category, const STRING& plot, 
				 const STRING& HLT, const unsigned int numBinsX, 
				 const float binLowEdgeX, const float binHighEdgeX, 
				 const unsigned int numBinsY, const float binLowEdgeY, 
				 const float binHighEdgeY)
{
  histVec.push_back(new TH2F(name(category, plot, HLT).c_str(), "", numBinsX, binLowEdgeX, 
			     binHighEdgeX, numBinsY, binLowEdgeY, binHighEdgeY));
}

void GMSBAnalyzer::bookEAndGHistograms(const STRING& HLT)
{
  bookHistogram(ggMET_, "gg", "MET", HLT, 30, 0.0, 150.0);
  bookHistogram(egMET_, "eg", "MET", HLT, 30, 0.0, 150.0);
  bookHistogram(eeMET_, "ee", "MET", HLT, 30, 0.0, 150.0);
  bookHistogram(ggDiEMET_, "gg", "diEMET", HLT, 30, 0.0, 150.0);
  bookHistogram(egDiEMET_, "eg", "diEMET", HLT, 30, 0.0, 150.0);
  bookHistogram(eeDiEMET_, "ee", "diEMET", HLT, 30, 0.0, 150.0);
  bookHistogram(ggPhoton1ET_, "gg", "photon1ET", HLT, 80, 0.0, 800.0);
  bookHistogram(ggPhoton1Eta_, "gg", "photon1Eta", HLT, 60, -3.0, 3.0);
  bookHistogram(ggPhoton1Phi_, "gg", "photon1Phi", HLT, 64, -3.2, 3.2);
  bookHistogram(ggPhoton1ECALIso_, "gg", "photon1ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ggPhoton1HCALIso_, "gg", "photon1HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ggPhoton1HOverE_, "gg", "photon1HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(ggPhoton1TrackIso_, "gg", "photon1TrackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ggPhoton1SigmaIetaIeta_, "gg", "photon1SigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(ggPhoton1SeedTime_, "gg", "photon1SeedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(ggPhoton1E2OverE9_, "gg", "photon1E2OverE9", HLT, 40, 0.0, 1.0);
  bookHistogram(ggPhoton2ET_, "gg", "photon2ET", HLT, 80, 0.0, 800.0);
  bookHistogram(ggPhoton2Eta_, "gg", "photon2Eta", HLT, 60, -3.0, 3.0);
  bookHistogram(ggPhoton2Phi_, "gg", "photon2Phi", HLT, 64, -3.2, 3.2);
  bookHistogram(ggPhoton2ECALIso_, "gg", "photon2ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ggPhoton2HCALIso_, "gg", "photon2HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ggPhoton2HOverE_, "gg", "photon2HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(ggPhoton2TrackIso_, "gg", "photon2TrackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ggPhoton2SigmaIetaIeta_, "gg", "photon2SigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(ggPhoton2SeedTime_, "gg", "photon2SeedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(ggPhoton2E2OverE9_, "gg", "photon2E2OverE9", HLT, 40, 0.0, 1.0);
  bookHistogram(egPhoton1ET_, "eg", "photon1ET", HLT, 80, 0.0, 800.0);
  bookHistogram(egPhoton1Eta_, "eg", "photon1Eta", HLT, 60, -3.0, 3.0);
  bookHistogram(egPhoton1Phi_, "eg", "photon1Phi", HLT, 64, -3.2, 3.2);
  bookHistogram(egPhoton1ECALIso_, "eg", "photon1ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(egPhoton1HCALIso_, "eg", "photon1HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(egPhoton1HOverE_, "eg", "photon1HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(egPhoton1TrackIso_, "eg", "photon1TrackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(egPhoton1SigmaIetaIeta_, "eg", "photon1SigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(egPhoton1SeedTime_, "eg", "photon1SeedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(egPhoton1E2OverE9_, "eg", "photon1E2OverE9", HLT, 40, 0.0, 1.0);
  bookHistogram(egPhoton2ET_, "eg", "photon2ET", HLT, 80, 0.0, 800.0);
  bookHistogram(egPhoton2Eta_, "eg", "photon2Eta", HLT, 60, -3.0, 3.0);
  bookHistogram(egPhoton2Phi_, "eg", "photon2Phi", HLT, 64, -3.2, 3.2);
  bookHistogram(egPhoton2ECALIso_, "eg", "photon2ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(egPhoton2HCALIso_, "eg", "photon2HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(egPhoton2HOverE_, "eg", "photon2HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(egPhoton2TrackIso_, "eg", "photon2TrackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(egPhoton2SigmaIetaIeta_, "eg", "photon2SigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(egPhoton2SeedTime_, "eg", "photon2SeedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(egPhoton2E2OverE9_, "eg", "photon2E2OverE9", HLT, 40, 0.0, 1.0);
  bookHistogram(eePhoton1ET_, "ee", "photon1ET", HLT, 80, 0.0, 800.0);
  bookHistogram(eePhoton1Eta_, "ee", "photon1Eta", HLT, 60, -3.0, 3.0);
  bookHistogram(eePhoton1Phi_, "ee", "photon1Phi", HLT, 64, -3.2, 3.2);
  bookHistogram(eePhoton1ECALIso_, "ee", "photon1ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(eePhoton1HCALIso_, "ee", "photon1HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(eePhoton1HOverE_, "ee", "photon1HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(eePhoton1TrackIso_, "ee", "photon1TrackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(eePhoton1SigmaIetaIeta_, "ee", "photon1SigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(eePhoton1SeedTime_, "ee", "photon1SeedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(eePhoton1E2OverE9_, "ee", "photon1E2OverE9", HLT, 40, 0.0, 1.0);
  bookHistogram(eePhoton2ET_, "ee", "photon2ET", HLT, 80, 0.0, 800.0);
  bookHistogram(eePhoton2Eta_, "ee", "photon2Eta", HLT, 60, -3.0, 3.0);
  bookHistogram(eePhoton2Phi_, "ee", "photon2Phi", HLT, 64, -3.2, 3.2);
  bookHistogram(eePhoton2ECALIso_, "ee", "photon2ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(eePhoton2HCALIso_, "ee", "photon2HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(eePhoton2HOverE_, "ee", "photon2HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(eePhoton2TrackIso_, "ee", "photon2TrackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(eePhoton2SigmaIetaIeta_, "ee", "photon2SigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(eePhoton2SeedTime_, "ee", "photon2SeedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(eePhoton2E2OverE9_, "ee", "photon2E2OverE9", HLT, 40, 0.0, 1.0);
}

void GMSBAnalyzer::bookFHistograms(const STRING& HLT)
{
  bookHistogram(ffMET_, "ff", "MET", HLT, 30, 0.0, 150.0);
  bookHistogram(ffDiEMET_, "ff", "diEMET", HLT, 30, 0.0, 150.0);
  bookHistogram(ffPhoton1ET_, "ff", "photon1ET", HLT, 80, 0.0, 800.0);
  bookHistogram(ffPhoton1Eta_, "ff", "photon1Eta", HLT, 60, -3.0, 3.0);
  bookHistogram(ffPhoton1Phi_, "ff", "photon1Phi", HLT, 64, -3.2, 3.2);
  bookHistogram(ffPhoton1ECALIso_, "ff", "photon1ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ffPhoton1HCALIso_, "ff", "photon1HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ffPhoton1HOverE_, "ff", "photon1HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(ffPhoton1TrackIso_, "ff", "photon1TrackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ffPhoton1SigmaIetaIeta_, "ff", "photon1SigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(ffPhoton1SeedTime_, "ff", "photon1SeedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(ffPhoton1E2OverE9_, "ff", "photon1E2OverE9", HLT, 40, 0.0, 1.0);
  bookHistogram(ffPhoton2ET_, "ff", "photon2ET", HLT, 80, 0.0, 800.0);
  bookHistogram(ffPhoton2Eta_, "ff", "photon2Eta", HLT, 60, -3.0, 3.0);
  bookHistogram(ffPhoton2Phi_, "ff", "photon2Phi", HLT, 64, -3.2, 3.2);
  bookHistogram(ffPhoton2ECALIso_, "ff", "photon2ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ffPhoton2HCALIso_, "ff", "photon2HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ffPhoton2HOverE_, "ff", "photon2HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(ffPhoton2TrackIso_, "ff", "photon2TrackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ffPhoton2SigmaIetaIeta_, "ff", "photon2SigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(ffPhoton2SeedTime_, "ff", "photon2SeedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(ffPhoton2E2OverE9_, "ff", "photon2E2OverE9", HLT, 40, 0.0, 1.0);
}

void GMSBAnalyzer::bookPhotonHistograms(const STRING& HLT)
{
  bookHistogram(ETNMinus1_, "", "ET", HLT, 80, 0.0, 800.0);
  bookHistogram(etaNMinus1_, "", "eta", HLT, 60, -3.0, 3.0);
  bookHistogram(phiNMinus1_, "", "phi", HLT, 64, -3.2, 3.2);
  bookHistogram(ECALIsoNMinus1_, "", "ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(HCALIsoNMinus1_, "", "HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(HOverENMinus1_, "", "HOverE", HLT, 20, 0.0, 0.2);
  bookHistogram(trackIsoNMinus1_, "", "trackIso", HLT, 20, 0.0, 40.0);
  bookHistogram(sigmaIetaIetaNMinus1_, "", "sigmaIetaIeta", HLT, 20, 0.0, 0.02);
  bookHistogram(seedTimeNMinus1_, "", "seedTime", HLT, 35, -3.5, 3.5);
  bookHistogram(e2OverE9NMinus1_, "", "e2OverE9", HLT, 40, 0.0, 1.0);
  bookHistogram(ECALIsoNMinus1PreselectionAll_, "preselectionAll", "ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(HCALIsoNMinus1PreselectionAll_, "preselectionAll", "HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ECALIsoNMinus1PreselectionDoublePhotonEvts_, "preselectionDoublePhotonEvts", 
		"ECALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(HCALIsoNMinus1PreselectionDoublePhotonEvts_, "preselectionDoublePhotonEvts", 
		"HCALIso", HLT, 20, 0.0, 40.0);
  bookHistogram(ECALIsoVsETNMinus1PreselectionAll_, "preselectionAll", "ECALIsoVsET", HLT, 80, 
		0.0, 800.0, 20, 0.0, 40.0);
  bookHistogram(HCALIsoVsETNMinus1PreselectionAll_, "preselectionAll", "HCALIsoVsET", HLT, 80, 
		0.0, 800.0, 20, 0.0, 40.0);
  bookHistogram(ECALIsoVsETNMinus1PreselectionDoublePhotonEvts_, "preselectionDoublePhotonEvts", 
		"ECALIsoVsET", HLT, 80, 0.0, 800.0, 20, 0.0, 40.0);
  bookHistogram(HCALIsoVsETNMinus1PreselectionDoublePhotonEvts_, "preselectionDoublePhotonEvts", 
		"HCALIsoVsET", HLT, 80, 0.0, 800.0, 20, 0.0, 40.0);
  bookHistogram(seedTimeVsENMinus1_, "", "seedTimeVsE", HLT, 80, 0.0, 800.0, 35, -3.5, 3.5);
  bookHistogram(numPassingECALIsoETAbsEtaVsNPVAll_, "all", "numPassingECALIsoETAbsEtaVsNPV", HLT, 
		20, 0.5, 20.5);
  bookHistogram(numPassingHCALIsoETAbsEtaVsNPVAll_, "all", "numPassingHCALIsoETAbsEtaVsNPV", HLT, 
		20, 0.5, 20.5);
  bookHistogram(numPassingHOverEETAbsEtaVsNPVAll_, "all", "numPassingHOverEETAbsEtaVsNPV", HLT, 
		20, 0.5, 20.5);
  bookHistogram(numPassingPreselectionVsNPVAll_, "all", "numPassingPreselectionVsNPV", HLT, 20, 
		0.5, 20.5);
  bookHistogram(numPassingETAbsEtaVsNPVAll_, "all", "numPassingETAbsEtaVsNPV", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingECALIsoETAbsEtaVsNPVDoublePhotonEvts_, "doublePhotonEvts", 
		"numPassingECALIsoETAbsEtaVsNPVDoublePhotonEvts", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingHCALIsoETAbsEtaVsNPVDoublePhotonEvts_, "doublePhotonEvts", 
		"numPassingHCALIsoETAbsEtaVsNPVDoublePhotonEvts", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingHOverEETAbsEtaVsNPVDoublePhotonEvts_, "doublePhotonEvts", 
		"numPassingHOverEETAbsEtaVsNPV", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingPreselectionVsNPVDoublePhotonEvts_, "doublePhotonEvts", 
		"numPassingPreselectionVsNPV", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingETAbsEtaVsNPVDoublePhotonEvts_, "doublePhotonEvts", 
		"numPassingETAbsEtaVsNPVDoublePhotonEvts", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingECALIsoETAbsEtaHOverEVsNPVAll_, "all", 
		"numPassingECALIsoETAbsEtaHOverEVsNPVAll", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingHCALIsoETAbsEtaHOverEVsNPVAll_, "all", 
		"numPassingHCALIsoETAbsEtaHOverEVsNPVAll", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingETAbsEtaHOverEVsNPVAll_, "all", "numPassingETAbsEtaHOverEVsNPVAll", HLT, 
		20, 0.5, 20.5);
  bookHistogram(numPassingECALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_, "doublePhotonEvts", 
		"numPassingECALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingHCALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_, "doublePhotonEvts", 
		"numPassingHCALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts", HLT, 20, 0.5, 20.5);
  bookHistogram(numPassingETAbsEtaHOverEVsNPVDoublePhotonEvts_, "doublePhotonEvts", 
		"numPassingETAbsEtaHOverEVsNPVDoublePhotonEvts", HLT, 20, 0.5, 20.5);
}

VSTRING GMSBAnalyzer::getRealTriggerNames(VSTRING* HLTPaths) const
{
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
    mess << "In getRealTriggerNames(), HLTPaths is null.\n";
    throw cms::Exception("GMSBAnalyzer") << mess.str();
  }
  return realHLT;
}

void GMSBAnalyzer::getPassingPhotons(edm::Handle<reco::PhotonCollection>& pPhotons, 
				     reco::Photon& photon1, reco::Photon& photon2, 
				     unsigned int& iPhoton1, unsigned int& iPhoton2, 
				     edm::Handle<edm::ValueMap<bool> >& pPassAll)
{
  double photon1ET = -1.0;
  double photon2ET = -1.0;
  unsigned int counter = 0;
  for (reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin(); 
       iPhoton != pPhotons->end(); ++iPhoton) {
    edm::Ref<reco::PhotonCollection> ref(pPhotons, iPhoton - pPhotons->begin());
    if ((*(pPassAll.product()))[ref]) {
      ++counter;
      if (iPhoton->et() > photon1ET) {
	photon2ET = photon1ET;
	photon1ET = iPhoton->et();
	photon2 = photon1;
	photon1 = *iPhoton;
	iPhoton2 = iPhoton1;
	iPhoton1 = iPhoton - pPhotons->begin();
      }
      else if (iPhoton->et() > photon2ET) {
	photon2ET = iPhoton->et();
	photon2 = *iPhoton;
	iPhoton2 = iPhoton - pPhotons->begin();
      }
    }
  }
  if (counter != 2) {
    throw cms::Exception("GMSBAnalyzer") << "Error: " << counter << " passing photons found.\n";
  }
}

void GMSBAnalyzer::fillIndHistogram(VTH1F* histVec, const float quantity, const bool pass, 
				    const STRING& HLT)
{
  if (histVec != NULL) {
    VTH1F_IT iHist = histVec->begin();
    bool foundHist = false;
    while ((iHist != histVec->end()) && !foundHist) {
      STRING histName((*iHist)->GetName());
      if ((histName.find(HLT) != STRING::npos) && pass) {
	(*iHist)->Fill(quantity);
	foundHist = true;
      }
      ++iHist;
    }
    if (pass && !foundHist) {
      STRINGSTREAM err;
      err << "Error: unable to find the right histogram for " << HLT << ".\n";
      edm::LogError("GMSBAnalzyer") << err.str();
    }
  }
  else edm::LogInfo("GMSBAnalyzerNullPointer") << "In fillIndHistogram, histVec is null.\n";
}

void GMSBAnalyzer::fillIndHistogram(VTH2F* histVec, const float quantity1, const float quantity2, 
				    const bool pass, const STRING& HLT)
{
  if (histVec != NULL) {
    VTH2F_IT iHist = histVec->begin();
    bool foundHist = false;
    while ((iHist != histVec->end()) && !foundHist) {
      STRING histName((*iHist)->GetName());
      if ((histName.find(HLT) != STRING::npos) && pass) {
	(*iHist)->Fill(quantity1, quantity2);
	foundHist = true;
      }
      ++iHist;
    }
    if (pass && !foundHist) {
      STRINGSTREAM err;
      err << "Error: unable to find the right histogram for " << HLT << ".\n";
      edm::LogError("GMSBAnalzyer") << err.str();
    }
  }
  else edm::LogInfo("GMSBAnalyzerNullPointer") << "In fillIndHistogram, histVec is null.\n";
}

void GMSBAnalyzer::fillHistograms(std::vector<VTH1F*>& histPtrVec, const VFLOAT& quantityVec, 
				  const VBOOL& passVec, const STRING& HLT)
{
  for (std::vector<VTH1F*>::iterator iHistPtrVec = histPtrVec.begin(); 
       iHistPtrVec != histPtrVec.end(); ++iHistPtrVec) {
    fillIndHistogram(*iHistPtrVec, quantityVec[iHistPtrVec - histPtrVec.begin()], 
		     passVec[iHistPtrVec - histPtrVec.begin()], HLT);
  }
}

void GMSBAnalyzer::fillHistograms(std::vector<VTH2F*>& histPtrVec, const VFLOAT& quantity1Vec, 
				  const VFLOAT& quantity2Vec, const VBOOL& passVec, 
				  const STRING& HLT)
{
  for (std::vector<VTH2F*>::iterator iHistPtrVec = histPtrVec.begin(); 
       iHistPtrVec != histPtrVec.end(); ++iHistPtrVec) {
    fillIndHistogram(*iHistPtrVec, quantity1Vec[iHistPtrVec - histPtrVec.begin()], 
		     quantity2Vec[iHistPtrVec - histPtrVec.begin()], 
		     passVec[iHistPtrVec - histPtrVec.begin()], HLT);
  }
}

bool GMSBAnalyzer::passNMinus1(std::vector<edm::Handle<edm::ValueMap<bool> >*>& passMaps, 
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

bool 
GMSBAnalyzer::twoPhotonsPassNMinus1(std::vector<edm::Handle<edm::ValueMap<bool> >*>& passMaps, 
				    edm::Handle<edm::ValueMap<bool> >* passMapToIgnore, 
				    const edm::Handle<reco::PhotonCollection>& pPhotons) const
{
  unsigned int numPassing = 0;
  reco::PhotonCollection::const_iterator iPhoton = pPhotons->begin();
  while ((iPhoton != pPhotons->end()) && (numPassing < 2)) {
    edm::Ref<reco::PhotonCollection> ref(pPhotons, iPhoton - pPhotons->begin());
    if (passNMinus1(passMaps, ref, passMapToIgnore)) ++numPassing;
    ++iPhoton;
  }
  return (numPassing >= 2);
}

unsigned int GMSBAnalyzer::numValidPVs(edm::Handle<reco::VertexCollection>& pVertices) const
{
  unsigned int numPassing = 0;
  for (reco::VertexCollection::const_iterator iVertex = pVertices->begin(); 
       iVertex != pVertices->end(); ++iVertex) {
    if (!(iVertex->isFake()) && (iVertex->ndof() > 4) && (fabs(iVertex->z()) <= 24.0/*cm*/) &&  
	(iVertex->position().Rho() <= 2.0/*cm*/)) ++numPassing;
  }
  return numPassing;
}

unsigned int GMSBAnalyzer::effDecision(std::vector<edm::Handle<edm::ValueMap<bool> >*>& passMaps, 
				       edm::Ref<reco::PhotonCollection>& ref, 
				       edm::Handle<edm::ValueMap<bool> >* passMapToIgnore, 
				       const unsigned int numPassingPVs, 
				       const edm::Handle<reco::PhotonCollection>& pPhotons) const
{
  unsigned int numPassingPVsForEffCalc = 0;
  if (pPhotons.isValid()) {
    if (twoPhotonsPassNMinus1(passMaps, passMapToIgnore, pPhotons)) {
      numPassingPVsForEffCalc = numPassingPVs;
    }
  }
  else if (passNMinus1(passMaps, ref, passMapToIgnore)) numPassingPVsForEffCalc = numPassingPVs;
  return numPassingPVsForEffCalc;
}

void GMSBAnalyzer::makeOutputDirectories(const STRING& dirName)
{
  out_->cd();
  out_->mkdir(dirName.c_str());
}

void GMSBAnalyzer::calculateEff(VTGRAPHASYMMERRORS& effGraphs, VTH1F& numeratorHists, 
				VTH1F& denominatorHists)
{
  //sanity checks
  effGraphs.clear();
  if (numeratorHists.size() != denominatorHists.size()) {
    STRINGSTREAM err;
    err << "Error in calculateEff(VTGRAPHASYMMERRORS& effGraphs, const VTH1F& numeratorHists, ";
    err << "const VTH1F& denominatorHists): numeratorHists.size() = " << numeratorHists.size();
    err << ", denominatorHists.size() = " << denominatorHists.size() << ".\n";
    edm::LogError("GMSBAnalyzer") << err.str();
    return;
  }

  //if bins are unfilled, set bin content to 1 and ignore results from those bins
  for (VTH1F_IT iNumeratorHist = numeratorHists.begin(); iNumeratorHist != numeratorHists.end(); 
       ++iNumeratorHist) {
    const unsigned int iHist = iNumeratorHist - numeratorHists.begin();
    STRINGSTREAM info;
    info << "Passing histogram: "  << (*iNumeratorHist)->GetName() << std::endl;
    info << "Total histogram: "  << denominatorHists[iHist]->GetName() << std::endl;
    edm::LogInfo("GMSBAnalyzerEfficiencyCalculation") << info.str();
    for (Int_t iBin = 0; iBin <= ((*iNumeratorHist)->GetNbinsX() + 2); ++iBin) {
      STRINGSTREAM info;
      info << "Bin " << iBin << ", #passing " << (*iNumeratorHist)->GetBinContent(iBin);
      info << ", #total " << denominatorHists[iHist]->GetBinContent(iBin);
      info << std::endl;
      edm::LogInfo("GMSBAnalyzerEfficiencyCalculation") << info.str();
      if ((iBin == 0) && ((*iNumeratorHist)->GetBinContent(iBin) > 
			  denominatorHists[iHist]->GetBinContent(iBin))) {
	(*iNumeratorHist)->SetBinContent(iBin, 1.0);
	denominatorHists[iNumeratorHist - numeratorHists.begin()]->SetBinContent(iBin, 1.0);
      }
    }
    effGraphs.push_back(new TGraphAsymmErrors(*iNumeratorHist, denominatorHists
					      [iNumeratorHist - numeratorHists.begin()], "cp"));
  }
}

void GMSBAnalyzer::writeGGPlots()
{
  makeOutputDirectories("gg");
  fillOutputDirectories(ggMET_, "gg");
  fillOutputDirectories(ggDiEMET_, "gg");
  fillOutputDirectories(ggPhoton1ET_, "gg");
  fillOutputDirectories(ggPhoton1Eta_, "gg");
  fillOutputDirectories(ggPhoton1Phi_, "gg");
  fillOutputDirectories(ggPhoton1ECALIso_, "gg");
  fillOutputDirectories(ggPhoton1HCALIso_, "gg");
  fillOutputDirectories(ggPhoton1HOverE_, "gg");
  fillOutputDirectories(ggPhoton1TrackIso_, "gg");
  fillOutputDirectories(ggPhoton1SigmaIetaIeta_, "gg");
  fillOutputDirectories(ggPhoton1SeedTime_, "gg");
  fillOutputDirectories(ggPhoton1E2OverE9_, "gg");
  fillOutputDirectories(ggPhoton2ET_, "gg");
  fillOutputDirectories(ggPhoton2Eta_, "gg");
  fillOutputDirectories(ggPhoton2Phi_, "gg");
  fillOutputDirectories(ggPhoton2ECALIso_, "gg");
  fillOutputDirectories(ggPhoton2HCALIso_, "gg");
  fillOutputDirectories(ggPhoton2HOverE_, "gg");
  fillOutputDirectories(ggPhoton2TrackIso_, "gg");
  fillOutputDirectories(ggPhoton2SigmaIetaIeta_, "gg");
  fillOutputDirectories(ggPhoton2SeedTime_, "gg");
  fillOutputDirectories(ggPhoton2E2OverE9_, "gg");
}

void GMSBAnalyzer::writeEGPlots()
{
  makeOutputDirectories("eg");
  fillOutputDirectories(egMET_, "eg");
  fillOutputDirectories(egDiEMET_, "eg");
  fillOutputDirectories(egPhoton1ET_, "eg");
  fillOutputDirectories(egPhoton1Eta_, "eg");
  fillOutputDirectories(egPhoton1Phi_, "eg");
  fillOutputDirectories(egPhoton1ECALIso_, "eg");
  fillOutputDirectories(egPhoton1HCALIso_, "eg");
  fillOutputDirectories(egPhoton1HOverE_, "eg");
  fillOutputDirectories(egPhoton1TrackIso_, "eg");
  fillOutputDirectories(egPhoton1SigmaIetaIeta_, "eg");
  fillOutputDirectories(egPhoton1SeedTime_, "eg");
  fillOutputDirectories(egPhoton1E2OverE9_, "eg");
  fillOutputDirectories(egPhoton2ET_, "eg");
  fillOutputDirectories(egPhoton2Eta_, "eg");
  fillOutputDirectories(egPhoton2Phi_, "eg");
  fillOutputDirectories(egPhoton2ECALIso_, "eg");
  fillOutputDirectories(egPhoton2HCALIso_, "eg");
  fillOutputDirectories(egPhoton2HOverE_, "eg");
  fillOutputDirectories(egPhoton2TrackIso_, "eg");
  fillOutputDirectories(egPhoton2SigmaIetaIeta_, "eg");
  fillOutputDirectories(egPhoton2SeedTime_, "eg");
  fillOutputDirectories(egPhoton2E2OverE9_, "eg");
}

void GMSBAnalyzer::writeEEPlots()
{
  makeOutputDirectories("ee");
  fillOutputDirectories(eeMET_, "ee");
  fillOutputDirectories(eeDiEMET_, "ee");
  fillOutputDirectories(eePhoton1ET_, "ee");
  fillOutputDirectories(eePhoton1Eta_, "ee");
  fillOutputDirectories(eePhoton1Phi_, "ee");
  fillOutputDirectories(eePhoton1ECALIso_, "ee");
  fillOutputDirectories(eePhoton1HCALIso_, "ee");
  fillOutputDirectories(eePhoton1HOverE_, "ee");
  fillOutputDirectories(eePhoton1TrackIso_, "ee");
  fillOutputDirectories(eePhoton1SigmaIetaIeta_, "ee");
  fillOutputDirectories(eePhoton1SeedTime_, "ee");
  fillOutputDirectories(eePhoton1E2OverE9_, "ee");
  fillOutputDirectories(eePhoton2ET_, "ee");
  fillOutputDirectories(eePhoton2Eta_, "ee");
  fillOutputDirectories(eePhoton2Phi_, "ee");
  fillOutputDirectories(eePhoton2ECALIso_, "ee");
  fillOutputDirectories(eePhoton2HCALIso_, "ee");
  fillOutputDirectories(eePhoton2HOverE_, "ee");
  fillOutputDirectories(eePhoton2TrackIso_, "ee");
  fillOutputDirectories(eePhoton2SigmaIetaIeta_, "ee");
  fillOutputDirectories(eePhoton2SeedTime_, "ee");
  fillOutputDirectories(eePhoton2E2OverE9_, "ee");
}

void GMSBAnalyzer::writeFFPlots()
{
  makeOutputDirectories("ff");
  fillOutputDirectories(ffMET_, "ff");
  fillOutputDirectories(ffDiEMET_, "ff");
  fillOutputDirectories(ffPhoton1ET_, "ff");
  fillOutputDirectories(ffPhoton1Eta_, "ff");
  fillOutputDirectories(ffPhoton1Phi_, "ff");
  fillOutputDirectories(ffPhoton1ECALIso_, "ff");
  fillOutputDirectories(ffPhoton1HCALIso_, "ff");
  fillOutputDirectories(ffPhoton1HOverE_, "ff");
  fillOutputDirectories(ffPhoton1TrackIso_, "ff");
  fillOutputDirectories(ffPhoton1SigmaIetaIeta_, "ff");
  fillOutputDirectories(ffPhoton1SeedTime_, "ff");
  fillOutputDirectories(ffPhoton1E2OverE9_, "ff");
  fillOutputDirectories(ffPhoton2ET_, "ff");
  fillOutputDirectories(ffPhoton2Eta_, "ff");
  fillOutputDirectories(ffPhoton2Phi_, "ff");
  fillOutputDirectories(ffPhoton2ECALIso_, "ff");
  fillOutputDirectories(ffPhoton2HCALIso_, "ff");
  fillOutputDirectories(ffPhoton2HOverE_, "ff");
  fillOutputDirectories(ffPhoton2TrackIso_, "ff");
  fillOutputDirectories(ffPhoton2SigmaIetaIeta_, "ff");
  fillOutputDirectories(ffPhoton2SeedTime_, "ff");
  fillOutputDirectories(ffPhoton2E2OverE9_, "ff");
}

void GMSBAnalyzer::writePhotonPlots()
{
  //write non-efficiency plots
  makeOutputDirectories("N-1_plots");  
  fillOutputDirectories(ETNMinus1_, "N-1_plots");
  fillOutputDirectories(etaNMinus1_, "N-1_plots");
  fillOutputDirectories(phiNMinus1_, "N-1_plots");
  fillOutputDirectories(ECALIsoNMinus1_, "N-1_plots");
  fillOutputDirectories(HCALIsoNMinus1_, "N-1_plots");
  fillOutputDirectories(HOverENMinus1_, "N-1_plots");
  fillOutputDirectories(trackIsoNMinus1_, "N-1_plots");
  fillOutputDirectories(sigmaIetaIetaNMinus1_, "N-1_plots");
  fillOutputDirectories(seedTimeNMinus1_, "N-1_plots");
  fillOutputDirectories(e2OverE9NMinus1_, "N-1_plots");
  fillOutputDirectories(ECALIsoNMinus1PreselectionAll_, "N-1_plots");
  fillOutputDirectories(HCALIsoNMinus1PreselectionAll_, "N-1_plots");
  fillOutputDirectories(ECALIsoNMinus1PreselectionDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(HCALIsoNMinus1PreselectionDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(seedTimeVsENMinus1_, "N-1_plots");
  fillOutputDirectories(numPassingECALIsoETAbsEtaVsNPVAll_, "N-1_plots");
  fillOutputDirectories(numPassingHCALIsoETAbsEtaVsNPVAll_, "N-1_plots");
  fillOutputDirectories(numPassingHOverEETAbsEtaVsNPVAll_, "N-1_plots");
  fillOutputDirectories(numPassingPreselectionVsNPVAll_, "N-1_plots");
  fillOutputDirectories(numPassingETAbsEtaVsNPVAll_, "N-1_plots");
  fillOutputDirectories(numPassingECALIsoETAbsEtaVsNPVDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(numPassingHCALIsoETAbsEtaVsNPVDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(numPassingHOverEETAbsEtaVsNPVDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(numPassingPreselectionVsNPVDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(numPassingETAbsEtaVsNPVDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(numPassingECALIsoETAbsEtaHOverEVsNPVAll_, "N-1_plots");
  fillOutputDirectories(numPassingHCALIsoETAbsEtaHOverEVsNPVAll_, "N-1_plots");
  fillOutputDirectories(numPassingETAbsEtaHOverEVsNPVAll_, "N-1_plots");
  fillOutputDirectories(numPassingECALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(numPassingHCALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_, "N-1_plots");
  fillOutputDirectories(numPassingETAbsEtaHOverEVsNPVDoublePhotonEvts_, "N-1_plots");

  //calculate efficiencies
  calculateEff(ECALIsoEffWRTAcceptanceVsNPVAll_, numPassingECALIsoETAbsEtaVsNPVAll_, 
	       numPassingETAbsEtaVsNPVAll_);
  calculateEff(HCALIsoEffWRTAcceptanceVsNPVAll_, numPassingHCALIsoETAbsEtaVsNPVAll_, 
	       numPassingETAbsEtaVsNPVAll_);
  calculateEff(HOverEEffWRTAcceptanceVsNPVAll_, numPassingHOverEETAbsEtaVsNPVAll_, 
	       numPassingETAbsEtaVsNPVAll_);
  calculateEff(preselectionEffWRTAcceptanceVsNPVAll_, numPassingPreselectionVsNPVAll_, 
	       numPassingETAbsEtaVsNPVAll_);
  calculateEff(ECALIsoEffWRTAcceptanceVsNPVDoublePhotonEvts_, 
	       numPassingECALIsoETAbsEtaVsNPVDoublePhotonEvts_, 
	       numPassingETAbsEtaVsNPVDoublePhotonEvts_);
  calculateEff(HCALIsoEffWRTAcceptanceVsNPVDoublePhotonEvts_, 
	       numPassingHCALIsoETAbsEtaVsNPVDoublePhotonEvts_, 
	       numPassingETAbsEtaVsNPVDoublePhotonEvts_);
  calculateEff(HOverEEffWRTAcceptanceVsNPVDoublePhotonEvts_, 
	       numPassingHOverEETAbsEtaVsNPVDoublePhotonEvts_, numPassingETAbsEtaVsNPVAll_);
  calculateEff(preselectionEffWRTAcceptanceVsNPVDoublePhotonEvts_, 
	       numPassingPreselectionVsNPVDoublePhotonEvts_, numPassingETAbsEtaVsNPVAll_);
  calculateEff(ECALIsoEffWRTETAbsEtaHOverEVsNPVAll_, numPassingECALIsoETAbsEtaHOverEVsNPVAll_, 
	       numPassingETAbsEtaHOverEVsNPVAll_);
  calculateEff(HCALIsoEffWRTETAbsEtaHOverEVsNPVAll_, numPassingHCALIsoETAbsEtaHOverEVsNPVAll_, 
	       numPassingETAbsEtaHOverEVsNPVAll_);
  calculateEff(ECALIsoEffWRTETAbsEtaHOverEVsNPVDoublePhotonEvts_, 
	       numPassingECALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_, 
	       numPassingETAbsEtaHOverEVsNPVDoublePhotonEvts_);
  calculateEff(HCALIsoEffWRTETAbsEtaHOverEVsNPVDoublePhotonEvts_, 
	       numPassingHCALIsoETAbsEtaHOverEVsNPVDoublePhotonEvts_, 
	       numPassingETAbsEtaHOverEVsNPVDoublePhotonEvts_);

  //write efficiency plots
  makeOutputDirectories("eff_plots");
  fillOutputDirectories(ECALIsoEffWRTAcceptanceVsNPVAll_, "eff_plots");
  fillOutputDirectories(HCALIsoEffWRTAcceptanceVsNPVAll_, "eff_plots");
  fillOutputDirectories(HOverEEffWRTAcceptanceVsNPVAll_, "eff_plots");
  fillOutputDirectories(preselectionEffWRTAcceptanceVsNPVAll_, "eff_plots");
  fillOutputDirectories(ECALIsoEffWRTAcceptanceVsNPVDoublePhotonEvts_, "eff_plots");
  fillOutputDirectories(HCALIsoEffWRTAcceptanceVsNPVDoublePhotonEvts_, "eff_plots");
  fillOutputDirectories(HOverEEffWRTAcceptanceVsNPVDoublePhotonEvts_, "eff_plots");
  fillOutputDirectories(preselectionEffWRTAcceptanceVsNPVDoublePhotonEvts_, "eff_plots");
  fillOutputDirectories(ECALIsoEffWRTETAbsEtaHOverEVsNPVAll_, "eff_plots");
  fillOutputDirectories(HCALIsoEffWRTETAbsEtaHOverEVsNPVAll_, "eff_plots");
  fillOutputDirectories(ECALIsoEffWRTETAbsEtaHOverEVsNPVDoublePhotonEvts_, "eff_plots");
  fillOutputDirectories(HCALIsoEffWRTETAbsEtaHOverEVsNPVDoublePhotonEvts_, "eff_plots");
}

//define this as a plug-in
DEFINE_FWK_MODULE(GMSBAnalyzer);
