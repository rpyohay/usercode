// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      TauAnalyzer
// 
/**\class TauAnalyzer TauAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/TauAnalyzer.cc

   Description: analyze tau variables that can discriminate signal to background

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: TauAnalyzer.cc,v 1.9 2012/12/06 17:44:24 yohay Exp $
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

//
// class declaration
//

class TauAnalyzer : public edm::EDAnalyzer {
public:
  explicit TauAnalyzer(const edm::ParameterSet&);
  ~TauAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  //delete memory
  void reset(const bool);

  //fill pT histogram with object of arbitrary type
  template<typename U>
  void fillPTHistogramArbitrary(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
  }

  //fill ET histogram
  template<typename U>
  void fillETHistogram(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->et());
  }

  //make visible gen tau pT canvas for 1 decay mode, multiple pT ranks
  void makePTRankCanvas(TCanvas&, TLegend&, const std::string&, std::vector<TH1F*>&);

  //format and draw multiple pT histograms on one canvas
  void drawMultiplePTHistograms(TCanvas&, std::vector<TH1F*>&, const std::vector<unsigned int>&, 
				const std::vector<unsigned int>&, TLegend&, 
				const std::vector<std::string>&, const std::string&);

  //plot histogram of dPhi(muon, MET) for user's choice of muon
  void plotDPhiMuMet(const reco::MuonRef&, const edm::Handle<edm::View<reco::PFMET> >&, TH1F*);

  //plot histogram of MT for user's choice of muon
  void plotMT(const reco::MuonRef&, const edm::Handle<edm::View<reco::PFMET> >&, TH1F*);

  //plot histogram of HT for user's choice of input candidates
  void plotHT(const std::vector<reco::Candidate*>&, TH1F*);

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //tau tag
  edm::InputTag tauTag_;

  //MET tag
  edm::InputTag METTag_;

  //muon tag
  edm::InputTag muonTag_;

  //gen-matched muon tag
  edm::InputTag genMatchedMuonTag_;

  //old jet tag
  edm::InputTag oldJetTag_;

  //new jet tag
  edm::InputTag newJetTag_;

  //jet-muon map tag
  edm::InputTag jetMuonMapTag_;

  //old-new jet map tag
  edm::InputTag oldNewJetMapTag_;

  //gen particle tag
  edm::InputTag genParticleTag_;

  //hadronic tau deltaBeta-corrected isolation energy tag
  edm::InputTag tauHadIsoTag_;

  //dR matching distance
  double dR_;

  //minimum tau pT
  double tauPTMin_;

  //tau decay mode
  reco::PFTau::hadronicDecayMode tauDecayMode_;

  //marker colors for histograms with different pT rank
  std::vector<unsigned int> pTRankColors_;

  //marker styles for histograms with different pT rank
  std::vector<unsigned int> pTRankStyles_;

  //legend entries for histograms with different pT rank
  std::vector<std::string> pTRankEntries_;

  //histogram of MET
  TH1F* MET_;

  //histograms of selected tau associated muon multiplicity
  TH1F* hadTauAssociatedMuMultiplicity_;

  //histogram of mu+had mass
  TH1F* muHadMass_;

  //histogram of mu charge + had charge
  TH1F* muHadCharge_;

  //histogram of W muon transverse mass
  TH1F* WMuMT_;

  //histogram of tau muon transverse mass
  TH1F* tauMuMT_;

  //histogram of dPhi(W muon, MET)
  TH1F* dPhiWMuMET_;

  //histogram of dPhi(tau muon, MET)
  TH1F* dPhiTauMuMET_;

  //histogram of HT (tau muon + hadronic tau + leading distinct corrected jet)
  TH1F* tauMuTauHadJetHT_;

  //histogram of HT (two leading corrected jets)
  TH1F* diJetHT_;

  //histogram of HT (corrected jet associated to hadronic tau + leading distinct corrected jet)
  TH1F* jetTauJetHT_;

  //histogram of HT (tau muon + hadronic tau + leading distinct corrected jet + W muon)
  TH1F* tauMuTauHadJetWMuHT_;

  //histogram of HT (two leading corrected jets + W muon)
  TH1F* diJetWMuHT_;

  /*histogram of HT (corrected jet associated to hadronic tau + leading distinct corrected jet + W 
    muon)*/
  TH1F* jetTauJetWMuHT_;

  //histogram of the parent parton of the jet
  TH1F* jetParentParton_;

  //histogram of dR(W muon, leading soft muon per jet)
  TH1F* dRWMuSoftMu_;

  //histogram of dR(W muon, leading soft gen-matched muon per jet)
  TH1F* dRWMuSoftGenMatchedMu_;

  //histogram of dR(W muon, leading soft muon per jet) when mu+had mass > 2 GeV
  TH1F* dRWMuSoftMuMuHadMassGe2_;

  //histogram of the type of gen tau matches
  TH1F* genTauMatchType_;

  //histogram of tau muon pT
  TH1F* tauMuPT_;

  //histogram of hadronic tau pT
  TH1F* tauHadPT_;

  //histogram of hadronic tau deltaBeta corrected isolation
  TH1F* tauHadIso_;

  /*histogram of the number of removed-muon-matched PF muons found in the isolation candidate 
    collection per tau*/
  TH1F* softMuIsoCandMultiplicity_;

  //histogram of cleaned jet pT vs. cleaned tau pT
  TH2F* cleanedJetPTVsCleanedTauPT_;

  //histogram of mu+had mass vs. dR(tagged soft muon, tau axis)
  TH2F* muHadMassVsDRSoftMuTau_;

  //histogram of mu+had mass vs. tagged soft muon pT for dR(tagged soft muon, tau axis) < 0.5
  TH2F* muHadMassVsSoftMuPTDRSoftMuTauLt05_;
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
TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  genMatchedMuonTag_(iConfig.getParameter<edm::InputTag>("genMatchedMuonTag")),
  oldJetTag_(iConfig.getParameter<edm::InputTag>("oldJetTag")),
  newJetTag_(iConfig.getParameter<edm::InputTag>("newJetTag")),
  jetMuonMapTag_(iConfig.getParameter<edm::InputTag>("jetMuonMapTag")),
  oldNewJetMapTag_(iConfig.getParameter<edm::InputTag>("oldNewJetMapTag")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  tauHadIsoTag_(iConfig.getParameter<edm::InputTag>("tauHadIsoTag")),
  dR_(iConfig.getParameter<double>("dR")),
  tauPTMin_(iConfig.getParameter<double>("tauPTMin")),
  tauDecayMode_(static_cast<reco::PFTau::hadronicDecayMode>
		(iConfig.getParameter<int>("tauDecayMode"))),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries"))
{
  //now do what ever initialization is needed
  reset(false);
}

TauAnalyzer::~TauAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);
  
  //get MET
  edm::Handle<edm::View<reco::PFMET> > pMET;
  iEvent.getByLabel(METTag_, pMET);

  //get muons
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //get gen-matched muons
  edm::Handle<reco::MuonRefVector> pGenMatchedMuons;
  iEvent.getByLabel(genMatchedMuonTag_, pGenMatchedMuons);

  //get old jets
  edm::Handle<reco::PFJetCollection> pOldJets;
  iEvent.getByLabel(oldJetTag_, pOldJets);

  //get new jets
  edm::Handle<reco::PFJetCollection> pNewJets;
  iEvent.getByLabel(newJetTag_, pNewJets);

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByLabel(jetMuonMapTag_, pMuonJetMap);

  //get old-new jet map
  edm::Handle<edm::ValueMap<reco::PFJetRef> > pOldNewJetMap;
  iEvent.getByLabel(oldNewJetMapTag_, pOldNewJetMap);

  //get gen particles
  edm::Handle<reco::GenParticleRefVector> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get hadronic tau deltaBeta-corrected isolation
  edm::Handle<reco::PFTauDiscriminator> pTauHadIso;
  iEvent.getByLabel(tauHadIsoTag_, pTauHadIso);

  //get AK5 PF L1FastL2L3 jet correction service
  const JetCorrector* corrector = JetCorrector::getJetCorrector("ak5PFL1FastL2L3", iSetup);

//   //debug
//   for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
//        ++iTau) {
//     const reco::PFJetRef& tauJetRef = (*iTau)->jetRef();
//     std::cerr << "Tau jet ref key: " << tauJetRef.key() << std::endl;
//     const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];
//     for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
// 	 iMuon != removedMuons.end(); ++iMuon) {
//       std::cerr << "Muon ref key: " << iMuon->key() << std::endl;
//       std::cerr << "Muon pT: " << (*iMuon)->pt() << " GeV\n";
//       std::cerr << "Muon eta: " << (*iMuon)->eta() << std::endl;
//       std::cerr << "Muon phi: " << (*iMuon)->phi() << std::endl;
//     }
//   }

  //plot MET distribution
  fillETHistogram(pMET, MET_);  

  //find the highest pT W muon
  std::vector<reco::MuonRef> WMuonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) { WMuonRefs.push_back(*iMuon); }
  Common::sortByPT(WMuonRefs);

  //plot dPhi(highest pT W muon, MET)
  plotDPhiMuMet(WMuonRefs[WMuonRefs.size() - 1], pMET, dPhiWMuMET_);

  //plot the transverse mass for the highest pT W muon
  plotMT(WMuonRefs[WMuonRefs.size() - 1], pMET, WMuMT_);

  //fill an STL container of pointers to the gen particles
  std::vector<reco::GenParticle*> genParticlePtrs;
  for (reco::GenParticleRefVector::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {
    genParticlePtrs.push_back(const_cast<reco::GenParticle*>((*iGenParticle).get()));
  }

  /*fill collection of all corrected jets
    N.B. "corrected jets" automatically excludes the jet associated to the W muon*/
  std::vector<reco::PFJet> correctedOldJets;
  for (reco::PFJetCollection::const_iterator iJet = pOldJets->begin(); iJet != pOldJets->end(); 
       ++iJet) {
    reco::PFJet correctedJet = *iJet;
    double JEC = corrector->correction(*iJet, iEvent, iSetup);
    correctedJet.scaleEnergy(JEC);
    if (reco::deltaR(*WMuonRefs[WMuonRefs.size() - 1], *iJet) >= dR_) { /*technically you should 
									  check that the W muon 
									  ref key is not 
									  identical to the ref 
									  key of any jet 
									  constituent, but that 
									  would take a much 
									  longer time*/
      correctedOldJets.push_back(correctedJet);
    }
  }

  //loop over selected taus
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) {
    const reco::PFJetRef& tauJetRef = (*iTau)->jetRef();
    const reco::PFJetRef& tauOldJetRef = (*pOldNewJetMap)[tauJetRef];
    const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];

    //find the highest pT associated muon
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	 iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
    Common::sortByPT(removedMuonRefs);

    /*fill collections of
      - corrected jets excluding the jet associated to the hadronic tau
      - the corrected jet associated to the hadronic tau
      N.B. "corrected jets" automatically excludes the jet associated to the W muon*/
    std::vector<reco::PFJet> correctedOldJetsExcludingTau;
    reco::PFJet correctedTauJet;
    for (reco::PFJetCollection::const_iterator iJet = correctedOldJets.begin(); 
	 iJet != correctedOldJets.end(); ++iJet) {
      const unsigned int jetIndex = iJet - pOldJets->begin();
      if (tauOldJetRef.key() != jetIndex) {
	correctedOldJetsExcludingTau.push_back(*iJet);
      }
      else correctedTauJet = *iJet;
    }

    //find the highest pT corrected jet in the original collection distinct from this tau
    std::vector<reco::PFJetRef> oldJetRefsExcludingTau;
    for (reco::PFJetCollection::const_iterator iCorrectedJetExcludingTau = 
	   correctedOldJetsExcludingTau.begin(); 
	 iCorrectedJetExcludingTau != correctedOldJetsExcludingTau.end(); 
	 ++iCorrectedJetExcludingTau) {
      oldJetRefsExcludingTau.
	push_back(reco::PFJetRef(&correctedOldJetsExcludingTau, 
				 iCorrectedJetExcludingTau - 
				 correctedOldJetsExcludingTau.begin()));
    }
    Common::sortByPT(oldJetRefsExcludingTau);

    //find the 2 highest pT corrected jets in the original collection
    std::vector<reco::PFJetRef> oldJetRefs;
    for (reco::PFJetCollection::const_iterator iCorrectedJet = correctedOldJets.begin(); 
	 iCorrectedJet != correctedOldJets.end(); ++iCorrectedJet) {
      oldJetRefs.
	push_back(reco::PFJetRef(&correctedOldJets, iCorrectedJet - correctedOldJets.begin()));
    }
    Common::sortByPT(oldJetRefs);

    /*find number of removed soft muon isolation candidates (corresponding to the highest pT 
      removed muon, so should be either 0 or 1)*/
    const reco::PFCandidateRefVector& PFChargedHadronIsoCands = 
      (*iTau)->isolationPFChargedHadrCands();
    std::vector<reco::PFCandidate*> PFChargedHadronIsoCandPtrs;
    for (reco::PFCandidateRefVector::const_iterator iCand = PFChargedHadronIsoCands.begin(); 
	 iCand != PFChargedHadronIsoCands.end(); ++iCand) {
      PFChargedHadronIsoCandPtrs.push_back(const_cast<reco::PFCandidate*>((*iCand).get()));
    }
    std::vector<double> sortedPFChargedHadronIsoCands = 
      Common::sortByProximity(PFChargedHadronIsoCandPtrs, 
			      *removedMuonRefs[removedMuonRefs.size() - 1]);

    //impose pT and decay mode cut on tau
    if (((*iTau)->pt() > tauPTMin_) && 
	((tauDecayMode_ == reco::PFTau::kNull) || ((*iTau)->decayMode() == tauDecayMode_))//  && 
// 	!((sortedPFChargedHadronIsoCands.size() > 0) && 
// 	  (sortedPFChargedHadronIsoCands[0] < dR_))
	) {

      //plot multiplicity of muons associated to this tau
      hadTauAssociatedMuMultiplicity_->Fill(removedMuons.size());

      //plot the mu + tau invariant mass for the highest pT muon
      const double muHadMass = 
	(removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();
      muHadMass_->Fill(muHadMass);

      //plot the mu + tau charge for the highest pT muon
      muHadCharge_->
	Fill(removedMuonRefs[removedMuonRefs.size() - 1]->charge() + (*iTau)->charge());

      //plot dPhi(highest pT muon, MET)
      plotDPhiMuMet(removedMuonRefs[removedMuonRefs.size() - 1], pMET, dPhiTauMuMET_);

      //plot the transverse mass for the highest pT muon
      plotMT(removedMuonRefs[removedMuonRefs.size() - 1], pMET, tauMuMT_);

      //plot HT (tau muon + hadronic tau + leading distinct corrected jet)
      std::vector<reco::Candidate*> tauMuTauHadJet;
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::Muon*>
				(removedMuonRefs[removedMuonRefs.size() - 1].get())));
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::PFTau*>(((*iTau).get()))));
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::PFJet*>
				(oldJetRefsExcludingTau
				 [oldJetRefsExcludingTau.size() - 1].get())));
      plotHT(tauMuTauHadJet, tauMuTauHadJetHT_);

      //plot HT (tau muon + hadronic tau + leading distinct corrected jet + W muon)
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(tauMuTauHadJet, tauMuTauHadJetWMuHT_);

      //plot HT (two leading corrected jets if both exist)
      std::vector<reco::Candidate*> diJet;
      diJet.push_back(dynamic_cast<reco::Candidate*>
		      (const_cast<reco::PFJet*>(oldJetRefs[oldJetRefs.size() - 1].get())));
      if (oldJetRefs.size() > 1) {
	diJet.push_back(dynamic_cast<reco::Candidate*>
			(const_cast<reco::PFJet*>(oldJetRefs[oldJetRefs.size() - 2].get())));
	plotHT(diJet, diJetHT_);
      }

      //plot HT (two leading corrected jets + W muon)
      diJet.push_back(dynamic_cast<reco::Candidate*>
		      (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(diJet, diJetWMuHT_);

      //plot HT (corrected jet associated to hadronic tau + leading distinct corrected jet)
      std::vector<reco::Candidate*> jetTauJet;
      jetTauJet.push_back(dynamic_cast<reco::Candidate*>
			  (const_cast<reco::PFJet*>(&correctedTauJet)));
      jetTauJet.push_back(dynamic_cast<reco::Candidate*>
			  (const_cast<reco::PFJet*>
			   (oldJetRefsExcludingTau[oldJetRefsExcludingTau.size() - 1].get())));
      plotHT(jetTauJet, jetTauJetHT_);

      /*plot HT (corrected jet associated to hadronic tau + leading distinct corrected jet + W 
	muon)*/
      jetTauJet.push_back(dynamic_cast<reco::Candidate*>
			  (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(jetTauJet, jetTauJetWMuHT_);

      //plot the parent parton of the tau jet
      int nearestPartonIndex = -1;
      const reco::GenParticle* nearestParton = 
	Common::nearestObject(tauOldJetRef, genParticlePtrs, nearestPartonIndex);
      if (nearestParton != NULL) {
	if (reco::deltaR(*nearestParton, *tauOldJetRef) < dR_) {
	  const unsigned int absNearestPartonPDGID = fabs(nearestParton->pdgId());
	  unsigned int bin = absNearestPartonPDGID == GenTauDecayID::G ? 0 : absNearestPartonPDGID;
	  jetParentParton_->Fill(bin);
	}
	else jetParentParton_->Fill(-1);
      }
      else jetParentParton_->Fill(7);

      //plot dR(W muon, leading soft muon per jet) when mu+had mass > 2 GeV
      if (muHadMass > 2.0/*GeV*/) {
	dRWMuSoftMuMuHadMassGe2_->Fill(reco::deltaR(*WMuonRefs[WMuonRefs.size() - 1], 
						    *removedMuonRefs[removedMuonRefs.size() - 1]));
      }

      //plot tau muon pT
      const double tauMuPT = removedMuonRefs[removedMuonRefs.size() - 1]->pt();
      tauMuPT_->Fill(tauMuPT);

      //plot hadronic tau pT
      tauHadPT_->Fill((*iTau)->pt());

      //plot hadronic tau deltaBeta-corrected isolation energy
      tauHadIso_->Fill((*pTauHadIso)[*iTau]);

      //plot soft muon isolation candidate multiplicity per removed muon
      softMuIsoCandMultiplicity_->Fill((sortedPFChargedHadronIsoCands.size() > 0) && 
				       (sortedPFChargedHadronIsoCands[0] < dR_));

      //plot cleaned jet pT vs. cleaned tau pT
      cleanedJetPTVsCleanedTauPT_->Fill((*iTau)->pt(), tauJetRef->pt());

      //plot mu+had mass vs. dR(tagged soft muon, tau axis)
      const double dRSoftMuTau = reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], **iTau);
      muHadMassVsDRSoftMuTau_->Fill(dRSoftMuTau, muHadMass);

      //plot mu+had mass vs. tagged soft muon pT for dR(tagged soft muon, tau axis) < 0.5
      if (dRSoftMuTau < 0.5) muHadMassVsSoftMuPTDRSoftMuTauLt05_->Fill(tauMuPT, muHadMass);

      //     //plot type of gen tau match
      //     if (std::find(genMatchedMuonRefs.begin(), genMatchedMuonRefs.end(), 
      // 		  removedMuonRefs[removedMuonRefs.size() - 1].key()) != genMatchedMuonRefs.end()) {
      //       genTauMatchType_->Fill();
      //     }
    }
  }

  //fill an STL container of gen-matched muon refs
  std::vector<unsigned int> genMatchedMuonRefs;
  for (reco::MuonRefVector::const_iterator iGenMatchedMu = pGenMatchedMuons->begin(); 
       iGenMatchedMu != pGenMatchedMuons->end(); ++iGenMatchedMu) {
    genMatchedMuonRefs.push_back(iGenMatchedMu->key());
  }

  //loop over cleaned jets
  for (reco::PFJetCollection::const_iterator iNewJet = pNewJets->begin(); 
       iNewJet != pNewJets->end(); ++iNewJet) {
    const reco::MuonRefVector& removedMuons = 
      (*pMuonJetMap)[reco::PFJetRef(pNewJets, iNewJet - pNewJets->begin())];

    //find the highest pT associated muon
    std::vector<reco::MuonRef> removedMuonRefs;
    for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	 iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
    Common::sortByPT(removedMuonRefs);
    if (removedMuonRefs.size() > 0) {
      const double dR = reco::deltaR(*WMuonRefs[WMuonRefs.size() - 1], 
				     *removedMuonRefs[removedMuonRefs.size() - 1]);

      //plot dR(W muon, leading soft muon per jet)
      dRWMuSoftMu_->Fill(dR);

      //plot dR(W muon, leading soft gen-matched muon per jet)
      if (std::find(genMatchedMuonRefs.begin(), genMatchedMuonRefs.end(), 
		    removedMuonRefs[removedMuonRefs.size() - 1].key()) != 
	  genMatchedMuonRefs.end()) dRWMuSoftGenMatchedMu_->Fill(dR);
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void TauAnalyzer::beginJob()
{
  //open output files
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book histograms
  MET_ = new TH1F("MET", ";#slash{E}_{T} (GeV);", 20, 0.0, 100.0);
  hadTauAssociatedMuMultiplicity_ = 
    new TH1F("hadTauAssociatedMuMultiplicity", ";N_{#mu}/#tau;", 2, 0.5, 2.5);
  muHadMass_ = new TH1F("muHadMass", ";m_{#mu+had} (GeV);", 20, 0.0, 20.0);
  muHadCharge_ = new TH1F("muHadCharge", ";q_{#mu} + q_{had};", 5, -2.5, 2.5);
  WMuMT_ = new TH1F("WMuMT", ";W muon M_{T} (GeV);", 50, 0.0, 200.0);
  tauMuMT_ = new TH1F("tauMuMT", ";#tau muon M_{T} (GeV);", 50, 0.0, 200.0);
  dPhiWMuMET_ = new TH1F("dPhiWMuMET", ";#Delta#phi(W muon, #slash{E}_{T});", 63, 0.0, 3.15);
  dPhiTauMuMET_ = 
    new TH1F("dPhiTauMuMET", ";#Delta#phi(#tau muon, #slash{E}_{T});", 63, 0.0, 3.15);
  tauMuTauHadJetHT_ = new TH1F("tauMuTauHadJetHT", ";H_{T} (GeV);", 50, 0.0, 200.0);
  diJetHT_ = new TH1F("diJetHT", ";H_{T} (GeV);", 50, 0.0, 200.0);
  jetTauJetHT_ = new TH1F("jetTauJetHT", ";H_{T} (GeV);", 50, 0.0, 200.0);
  tauMuTauHadJetWMuHT_ = new TH1F("tauMuTauHadJetWMuHT", ";H_{T} (GeV);", 50, 0.0, 200.0);
  diJetWMuHT_ = new TH1F("diJetWMuHT", ";H_{T} (GeV);", 50, 0.0, 200.0);
  jetTauJetWMuHT_ = new TH1F("jetTauJetWMuHT", ";H_{T} (GeV);", 50, 0.0, 200.0);
  jetParentParton_ = new TH1F("jetParentParton", ";Jet parent parton;", 7, -0.5, 6.5);
  dRWMuSoftMu_ = new TH1F("dRWMuSoftMu", ";#DeltaR(W muon, soft muon);", 30, 0.0, 3.0);
  dRWMuSoftGenMatchedMu_ = 
    new TH1F("dRWMuSoftGenMatchedMu", ";#DeltaR(W muon, soft muon);", 30, 0.0, 3.0);
  dRWMuSoftMuMuHadMassGe2_ = 
    new TH1F("dRWMuSoftMuMuHadMassGe2", ";#DeltaR(W muon, soft muon);", 30, 0.0, 3.0);
  genTauMatchType_ = new TH1F("genTauMatchType", ";;", 3, -0.5, 2.5);
  tauMuPT_ = new TH1F("tauMuPT", ";p_{T} (GeV);", 20, 0.0, 100.0);
  tauHadPT_ = new TH1F("tauHadPT", ";p_{T} (GeV);", 20, 0.0, 100.0);
  tauHadIso_ = new TH1F("tauHadIso", ";Isolation energy (GeV);", 20, 0.0, 20.0);
  softMuIsoCandMultiplicity_ = new TH1F("softMuIsoCandMultiplicity", 
					";No. removed soft muon isolation candidates;", 
					3, -0.5, 2.5);
  cleanedJetPTVsCleanedTauPT_ = 
    new TH2F("cleanedJetPTVsCleanedTauPT", ";#tau p_{T} (GeV);Jet p_{T} (GeV)", 
	     50, 0.0, 100.0, 50, 0.0, 100.0);
  muHadMassVsDRSoftMuTau_ = 
    new TH2F("muHadMassVsDRSoftMuTau", ";#DeltaR(soft muon, #tau_{had});m_{#mu+had} (GeV)", 
	     30, 0.0, 3.0, 20, 0.0, 20.0);
  muHadMassVsSoftMuPTDRSoftMuTauLt05_ = 
    new TH2F("muHadMassVsSoftMuPTDRSoftMuTauLt05", ";p_{T} (GeV);m_{#mu+had} (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 20.0);

  //set bin labels
  jetParentParton_->GetXaxis()->SetBinLabel(1, "g");
  jetParentParton_->GetXaxis()->SetBinLabel(2, "d");
  jetParentParton_->GetXaxis()->SetBinLabel(3, "u");
  jetParentParton_->GetXaxis()->SetBinLabel(4, "s");
  jetParentParton_->GetXaxis()->SetBinLabel(5, "c");
  jetParentParton_->GetXaxis()->SetBinLabel(6, "b");
  jetParentParton_->GetXaxis()->SetBinLabel(7, "t");
  genTauMatchType_->GetXaxis()->SetBinLabel(1, "#tau_{#mu}#tau_{e}");
  genTauMatchType_->GetXaxis()->SetBinLabel(2, "#tau_{#mu}#tau_{had}");
  genTauMatchType_->GetXaxis()->SetBinLabel(3, "#tau_{#mu}#tau_{#mu}");
}

// ------------ method called once each job just after ending the event loop  ------------
void TauAnalyzer::endJob() 
{
  //make canvases
  out_->cd();
  TCanvas METCanvas("METCanvas", "", 600, 600);
  TCanvas 
    hadTauAssociatedMuMultiplicityCanvas("hadTauAssociatedMuMultiplicityCanvas", "", 600, 600);
  TCanvas muHadMassCanvas("muHadMassCanvas", "", 600, 600);
  TCanvas muHadChargeCanvas("muHadChargeCanvas", "", 600, 600);
  TCanvas WMuMTCanvas("WMuMTCanvas", "", 600, 600);
  TCanvas tauMuMTCanvas("tauMuMTCanvas", "", 600, 600);
  TCanvas dPhiWMuMETCanvas("dPhiWMuMETCanvas", "", 600, 600);
  TCanvas dPhiTauMuMETCanvas("dPhiTauMuMETCanvas", "", 600, 600);
  TCanvas tauMuTauHadJetHTCanvas("tauMuTauHadJetHTCanvas", "", 600, 600);
  TCanvas diJetHTCanvas("diJetHTCanvas", "", 600, 600);
  TCanvas jetTauJetHTCanvas("jetTauJetHTCanvas", "", 600, 600);
  TCanvas tauMuTauHadJetWMuHTCanvas("tauMuTauHadJetWMuHTCanvas", "", 600, 600);
  TCanvas diJetWMuHTCanvas("diJetWMuHTCanvas", "", 600, 600);
  TCanvas jetTauJetWMuHTCanvas("jetTauJetWMuHTCanvas", "", 600, 600);
  TCanvas jetParentPartonCanvas("jetParentPartonCanvas", "", 600, 600);
  TCanvas dRWMuSoftMuCanvas("dRWMuSoftMuCanvas", "", 600, 600);
  TCanvas dRWMuSoftGenMatchedMuCanvas("dRWMuSoftGenMatchedMuCanvas", "", 600, 600);
  TCanvas dRWMuSoftMuMuHadMassGe2Canvas("dRWMuSoftMuMuHadMassGe2Canvas", "", 600, 600);
  TCanvas genTauMatchTypeCanvas("genTauMatchTypeCanvas", "", 600, 600);
  TCanvas tauMuPTCanvas("tauMuPTCanvas", "", 600, 600);
  TCanvas tauHadPTCanvas("tauHadPTCanvas", "", 600, 600);
  TCanvas tauHadIsoCanvas("tauHadIsoCanvas", "", 600, 600);
  TCanvas softMuIsoCandMultiplicityCanvas("softMuIsoCandMultiplicityCanvas", "", 600, 600);
  TCanvas cleanedJetPTVsCleanedTauPTCanvas("cleanedJetPTVsCleanedTauPTCanvas", "", 600, 600);
  TCanvas muHadMassVsDRSoftMuTauCanvas("muHadMassVsDRSoftMuTauCanvas", "", 600, 600);
  TCanvas muHadMassVsSoftMuPTDRSoftMuTauLt05Canvas("muHadMassVsSoftMuPTDRSoftMuTauLt05Canvas", "", 600, 600);

  //format and draw 1D plots
  Common::draw1DHistograms(METCanvas, MET_);
  Common::draw1DHistograms(hadTauAssociatedMuMultiplicityCanvas, hadTauAssociatedMuMultiplicity_);
  Common::draw1DHistograms(muHadMassCanvas, muHadMass_);
  Common::draw1DHistograms(muHadChargeCanvas, muHadCharge_);
  Common::draw1DHistograms(WMuMTCanvas, WMuMT_);
  Common::draw1DHistograms(tauMuMTCanvas, tauMuMT_);
  Common::draw1DHistograms(dPhiWMuMETCanvas, dPhiWMuMET_);
  Common::draw1DHistograms(dPhiTauMuMETCanvas, dPhiTauMuMET_);
  Common::draw1DHistograms(tauMuTauHadJetHTCanvas, tauMuTauHadJetHT_);
  Common::draw1DHistograms(diJetHTCanvas, diJetHT_);
  Common::draw1DHistograms(jetTauJetHTCanvas, jetTauJetHT_);
  Common::draw1DHistograms(tauMuTauHadJetWMuHTCanvas, tauMuTauHadJetWMuHT_);
  Common::draw1DHistograms(diJetWMuHTCanvas, diJetWMuHT_);
  Common::draw1DHistograms(jetTauJetWMuHTCanvas, jetTauJetWMuHT_);
  Common::draw1DHistograms(jetParentPartonCanvas, jetParentParton_);
  Common::draw1DHistograms(dRWMuSoftMuCanvas, dRWMuSoftMu_);
  Common::draw1DHistograms(dRWMuSoftGenMatchedMuCanvas, dRWMuSoftGenMatchedMu_);
  Common::draw1DHistograms(dRWMuSoftMuMuHadMassGe2Canvas, dRWMuSoftMuMuHadMassGe2_);
  Common::draw1DHistograms(tauMuPTCanvas, tauMuPT_);
  Common::draw1DHistograms(tauHadPTCanvas, tauHadPT_);
  Common::draw1DHistograms(tauHadIsoCanvas, tauHadIso_);
  Common::draw1DHistograms(genTauMatchTypeCanvas, genTauMatchType_);
  Common::draw1DHistograms(softMuIsoCandMultiplicityCanvas, softMuIsoCandMultiplicity_);

  //format and draw 2D plots
  Common::draw2DHistograms(cleanedJetPTVsCleanedTauPTCanvas, cleanedJetPTVsCleanedTauPT_);
  Common::draw2DHistograms(muHadMassVsDRSoftMuTauCanvas, muHadMassVsDRSoftMuTau_);
  Common::draw2DHistograms(muHadMassVsSoftMuPTDRSoftMuTauLt05Canvas, muHadMassVsSoftMuPTDRSoftMuTauLt05_);

  //set custom options
  softMuIsoCandMultiplicity_->GetXaxis()->SetTitleSize(0.05);
  softMuIsoCandMultiplicity_->GetXaxis()->SetTitleOffset(1.05);

  //write output files
  out_->cd();
  METCanvas.Write();
  hadTauAssociatedMuMultiplicityCanvas.Write();
  muHadMassCanvas.Write();
  muHadChargeCanvas.Write();
  WMuMTCanvas.Write();
  tauMuMTCanvas.Write();
  dPhiWMuMETCanvas.Write();
  dPhiTauMuMETCanvas.Write();
  tauMuTauHadJetHTCanvas.Write();
  diJetHTCanvas.Write();
  jetTauJetHTCanvas.Write();
  tauMuTauHadJetWMuHTCanvas.Write();
  diJetWMuHTCanvas.Write();
  jetTauJetWMuHTCanvas.Write();
  jetParentPartonCanvas.Write();
  dRWMuSoftMuCanvas.Write();
  dRWMuSoftGenMatchedMuCanvas.Write();
  dRWMuSoftMuMuHadMassGe2Canvas.Write();
  genTauMatchTypeCanvas.Write();
  tauMuPTCanvas.Write();
  tauHadPTCanvas.Write();
  tauHadIsoCanvas.Write();
  softMuIsoCandMultiplicityCanvas.Write();
  cleanedJetPTVsCleanedTauPTCanvas.Write();
  muHadMassVsDRSoftMuTauCanvas.Write();
  muHadMassVsSoftMuPTDRSoftMuTauLt05Canvas.Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void TauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void TauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void TauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
				       edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
				     edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TauAnalyzer::makePTRankCanvas(TCanvas& canvas, TLegend& legend, 
				   const std::string& header, 
				   std::vector<TH1F*>& hists)
{
  drawMultiplePTHistograms(canvas, hists, pTRankColors_, pTRankStyles_, legend, 
			   pTRankEntries_, header);
}

void TauAnalyzer::drawMultiplePTHistograms(TCanvas& canvas, 
					   std::vector<TH1F*>& hists, 
					   const std::vector<unsigned int>& colors, 
					   const std::vector<unsigned int>& styles, 
					   TLegend& legend, 
					   const std::vector<std::string>& entries, 
					   const std::string& header)
{
  Common::setLegendOptions(legend, header.c_str());
  canvas.cd();
  Common::setCanvasOptions(canvas, 1, 0, 0);
  TH1F* pHistWithMaxMaxBin = NULL;
  Double_t maxBinContent = 0.0;
  for (std::vector<TH1F*>::iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    const unsigned int i = iHist - hists.begin();
    Common::setHistogramOptions(*iHist, colors[i], 0.7, styles[i], 1.0, 
				"p_{T} (GeV)", "", 0.04);
    legend.AddEntry(*iHist, entries[i].c_str(), "l");
    Int_t histMaxBin = (*iHist)->GetMaximumBin();
    Double_t histMaxBinContent = (*iHist)->GetBinContent(histMaxBin);
    if (histMaxBinContent >= maxBinContent) {
      maxBinContent = histMaxBinContent;
      pHistWithMaxMaxBin = *iHist;
    }
  }
  pHistWithMaxMaxBin->Draw();
  for (std::vector<TH1F*>::iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    (*iHist)->Draw("SAME");
  }
  legend.Draw();
  canvas.Write();
}

void TauAnalyzer::plotDPhiMuMet(const reco::MuonRef& muonRef, 
				const edm::Handle<edm::View<reco::PFMET> >& pMET, TH1F* hist)
{
  hist->Fill(fabs(reco::deltaPhi(muonRef->phi(), pMET->refAt(0)->phi())));
}

void TauAnalyzer::plotMT(const reco::MuonRef& muonRef, 
			 const edm::Handle<edm::View<reco::PFMET> >& pMET, TH1F* hist)
{
  edm::RefToBase<reco::PFMET> METRefToBase = pMET->refAt(0);
  hist->Fill(sqrt(2*muonRef->pt()*METRefToBase->et()*
		  (1.0 - cos(reco::deltaPhi(muonRef->phi(), METRefToBase->phi())))));
}

void TauAnalyzer::plotHT(const std::vector<reco::Candidate*>& cands, TH1F* hist)
{
  double HT = 0.0;
  for (std::vector<reco::Candidate*>::const_iterator iCand = cands.begin(); iCand != cands.end(); 
       ++iCand) {
    HT+=(*iCand)->pt();
  }
  hist->Fill(HT);
}

void TauAnalyzer::reset(const bool doDelete)
{
  if (doDelete && (out_ != NULL)) delete out_;
  out_ = NULL;
  if (doDelete && (MET_ != NULL)) delete MET_;
  MET_ = NULL;
  if (doDelete && (hadTauAssociatedMuMultiplicity_ != NULL)) {
    delete hadTauAssociatedMuMultiplicity_;
  }
  hadTauAssociatedMuMultiplicity_ = NULL;
  if (doDelete && (muHadMass_ != NULL)) delete muHadMass_;
  muHadMass_ = NULL;
  if (doDelete && (muHadCharge_ != NULL)) delete muHadCharge_;
  muHadCharge_ = NULL;
  if (doDelete && (WMuMT_ != NULL)) delete WMuMT_;
  WMuMT_ = NULL;
  if (doDelete && (tauMuMT_ != NULL)) delete tauMuMT_;
  tauMuMT_ = NULL;
  if (doDelete && (dPhiWMuMET_ != NULL)) delete dPhiWMuMET_;
  dPhiWMuMET_ = NULL;
  if (doDelete && (dPhiTauMuMET_ != NULL)) delete dPhiTauMuMET_;
  dPhiTauMuMET_ = NULL;
  if (doDelete && (tauMuTauHadJetHT_ != NULL)) delete tauMuTauHadJetHT_;
  tauMuTauHadJetHT_ = NULL;
  if (doDelete && (diJetHT_ != NULL)) delete diJetHT_;
  diJetHT_ = NULL;
  if (doDelete && (jetTauJetHT_ != NULL)) delete jetTauJetHT_;
  jetTauJetHT_ = NULL;
  if (doDelete && (tauMuTauHadJetWMuHT_ != NULL)) delete tauMuTauHadJetWMuHT_;
  tauMuTauHadJetWMuHT_ = NULL;
  if (doDelete && (diJetWMuHT_ != NULL)) delete diJetWMuHT_;
  diJetWMuHT_ = NULL;
  if (doDelete && (jetTauJetWMuHT_ != NULL)) delete jetTauJetWMuHT_;
  jetTauJetWMuHT_ = NULL;
  if (doDelete && (jetParentParton_ != NULL)) delete jetParentParton_;
  jetParentParton_ = NULL;
  if (doDelete && (dRWMuSoftMu_ != NULL)) delete dRWMuSoftMu_;
  dRWMuSoftMu_ = NULL;
  if (doDelete && (dRWMuSoftGenMatchedMu_ != NULL)) delete dRWMuSoftGenMatchedMu_;
  dRWMuSoftGenMatchedMu_ = NULL;
  if (doDelete && (dRWMuSoftMuMuHadMassGe2_ != NULL)) delete dRWMuSoftMuMuHadMassGe2_;
  dRWMuSoftMuMuHadMassGe2_ = NULL;
  if (doDelete && (genTauMatchType_ != NULL)) delete genTauMatchType_;
  genTauMatchType_ = NULL;
  if (doDelete && (tauMuPT_ != NULL)) delete tauMuPT_;
  tauMuPT_ = NULL;
  if (doDelete && (tauHadPT_ != NULL)) delete tauHadPT_;
  tauHadPT_ = NULL;
  if (doDelete && (tauHadIso_ != NULL)) delete tauHadIso_;
  tauHadIso_ = NULL;
  if (doDelete && (softMuIsoCandMultiplicity_ != NULL)) delete softMuIsoCandMultiplicity_;
  softMuIsoCandMultiplicity_ = NULL;
  if (doDelete && (cleanedJetPTVsCleanedTauPT_ != NULL)) delete cleanedJetPTVsCleanedTauPT_;
  cleanedJetPTVsCleanedTauPT_ = NULL;
  if (doDelete && (muHadMassVsDRSoftMuTau_ != NULL)) delete muHadMassVsDRSoftMuTau_;
  muHadMassVsDRSoftMuTau_ = NULL;
  if (doDelete && (muHadMassVsSoftMuPTDRSoftMuTauLt05_ != NULL)) delete muHadMassVsSoftMuPTDRSoftMuTauLt05_;
  muHadMassVsSoftMuPTDRSoftMuTauLt05_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
