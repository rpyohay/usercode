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
// $Id: TauAnalyzer.cc,v 1.10 2012/12/12 16:04:07 yohay Exp $
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
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
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
  void fillETHistogram(edm::Handle<edm::View<U> >& pView, TH1F* hist, const double weight)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->et(), weight);
  }

  //make visible gen tau pT canvas for 1 decay mode, multiple pT ranks
  void makePTRankCanvas(TCanvas&, TLegend&, const std::string&, std::vector<TH1F*>&);

  //format and draw multiple pT histograms on one canvas
  void drawMultiplePTHistograms(TCanvas&, std::vector<TH1F*>&, const std::vector<unsigned int>&, 
				const std::vector<unsigned int>&, TLegend&, 
				const std::vector<std::string>&, const std::string&);

  //plot histogram of dPhi(muon, MET) for user's choice of muon
  void plotDPhiMuMet(const reco::MuonRef&, const edm::Handle<edm::View<reco::PFMET> >&, TH1F*, 
		     const double);

  //plot histogram of MT for user's choice of muon
  void plotMT(const reco::MuonRef&, const edm::Handle<edm::View<reco::PFMET> >&, TH1F*, 
	      const double);

  //plot histogram of HT for user's choice of input candidates
  void plotHT(const std::vector<reco::Candidate*>&, TH1F*, const double);

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

  //all muons tag
  edm::InputTag allMuonTag_;

  //muon gen particle tag
  edm::InputTag muonGenParticleTag_;

  //PU info tag
  edm::InputTag PUTag_;

  //vertex tag
  edm::InputTag vtxTag_;

  //dR matching distance
  double dR_;

  //minimum tau pT
  double tauPTMin_;

  //tau decay mode
  reco::PFTau::hadronicDecayMode tauDecayMode_;

  //minimum uncorrected jet pT
  double uncorrJetPTMin_;

  //tau arbitration method when there are multiple selected taus per event
  std::string tauArbitrationMethod_;

  //PU reweighting scenario
  std::string PUScenario_;

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

  //histogram of hadronic tau eta
  TH1F* tauHadEta_;

  //histogram of pT(soft muon)/(mu+had mass)
  TH1F* softMuPTOverMuHadMass_;

  //histogram of pT(mu+had)/(mu+had mass)
  TH1F* muHadPTOverMuHadMass_;

  //histogram of dR(soft muon, nearest gen muon)
  TH1F* dRSoftMuNearestGenMuHist_;

  //histogram of mu+had pT
  TH1F* muHadPT_;

  //histogram of number of good vertices (to check PU reweighting)
  TH1F* nGoodVtx_;

  //histogram of cleaned jet pT vs. cleaned tau pT
  TH2F* cleanedJetPTVsCleanedTauPT_;

  //histogram of uncleaned jet pT vs. cleaned tau pT
  TH2F* uncleanedJetPTVsCleanedTauPT_;

  //histogram of mu+had mass vs. dR(tagged soft muon, tau axis)
  TH2F* muHadMassVsDRSoftMuTau_;

  //histogram of mu+had mass vs. tagged soft muon pT
  TH2F* muHadMassVsSoftMuPT_;

  //histogram of hadronic tau isolation vs. soft muon pT
  TH2F* tauHadIsoVsSoftMuPT_;

  //histogram of properties of nearest muon to soft muon
  TH2F* genMuExistsVsSoftMuNearestMuProperties_;

  //histogram of mu+had multiplicity
  TH1F* muHadMultiplicity_;

  //histogram of mu+had mass vs. hadronic tau eta
  TH2F* muHadMassVsTauHadEta_;

  //histogram of mu+had mass vs. soft muon eta
  TH2F* muHadMassVsSoftMuEta_;

  //histogram of mu+had mass vs. hadronic tau isolation
  TH2F* muHadMassVsTauHadIso_;

  //histogram of mu+had mass vs. hadronic tau pT
  TH2F* muHadMassVsTauHadPT_;

  //histogram of hadronic tau isolation vs. eta
  TH2F* tauHadIsoVsEta_;

  //histogram of hadronic tau eta vs. soft muon eta
  TH2F* tauHadEtaVsSoftMuEta_;

  //histogram of dEta(hadronic tau, soft muon) vs. dPhi(hadronic tau, soft muon)
  TH2F* dEtaTauHadSoftMuVsDPhiTauHadSoftMu_;

  //histogram of pT(hadronic tau)/(mu+had mass) vs. hadronic tau isolation
  TH2F* tauHadPTOverMuHadMassVsTauHadIso_;

  //histogram of pT(soft muon)/(mu+had mass) vs. hadronic tau isolation
  TH2F* softMuPTOverMuHadMassVsTauHadIso_;

  //histogram of (pT(soft muon) + pT(hadronic tau)/(2(mu+had mass)) vs. hadronic tau isolation
  TH2F* avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_;

  //histogram of pT(mu+had)/(mu+had mass) vs. hadronic tau isolation
  TH2F* muHadPTOverMuHadMassVsTauHadIso_;

  //histogram of soft muon pT vs. hadronic tau pT
  TH2F* softMuPTVsTauHadPT_;

  //PU reweighting object
  edm::LumiReWeighting PUReweight_;
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
  allMuonTag_(iConfig.getParameter<edm::InputTag>("allMuonTag")),
  muonGenParticleTag_(iConfig.getParameter<edm::InputTag>("muonGenParticleTag")),
  PUTag_(iConfig.getParameter<edm::InputTag>("PUTag")),
  vtxTag_(iConfig.getParameter<edm::InputTag>("vtxTag")),
  dR_(iConfig.getParameter<double>("dR")),
  tauPTMin_(iConfig.getParameter<double>("tauPTMin")),
  tauDecayMode_(static_cast<reco::PFTau::hadronicDecayMode>
		(iConfig.getParameter<int>("tauDecayMode"))),
  uncorrJetPTMin_(iConfig.getParameter<double>("uncorrJetPTMin")),
  tauArbitrationMethod_(iConfig.getParameter<std::string>("tauArbitrationMethod")),
  PUScenario_(iConfig.getParameter<std::string>("PUScenario")),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries"))
{
  //now do what ever initialization is needed
  reset(false);

  //check that tau arbitration method is valid
  if ((tauArbitrationMethod_ != "pT") && (tauArbitrationMethod_ != "m") && 
      (tauArbitrationMethod_ != "none")) {
    std::stringstream err;
    err << "Error: unsupported tau arbitration method.  Valid choices are:\n";
    err << "---\"pT\": highest pT selected tau chosen per event\n";
    err << "---\"m\": highest mu+had mass selected tau chosen per event\n";
    err << "--\"none\": all selected taus are used\n";
    throw cms::Exception("TauAnalyzer") << err.str();
  }

  //instantiate the PU reweighting object
  float Data2012PUDistArray[60] = {1.226e+04, 3.285e+04, 9.233e+04, 3.395e+05, 6.185e+05, 
				   3.050e+06, 1.772e+07, 5.414e+07, 1.305e+08, 2.590e+08, 
				   4.463e+08, 6.856e+08, 8.816e+08, 9.991e+08, 1.079e+09, 
				   1.138e+09, 1.172e+09, 1.182e+09, 1.177e+09, 1.161e+09, 
				   1.136e+09, 1.105e+09, 1.068e+09, 1.021e+09, 9.556e+08, 
				   8.671e+08, 7.587e+08, 6.389e+08, 5.164e+08, 3.999e+08, 
				   2.963e+08, 2.101e+08, 1.424e+08, 9.205e+07, 5.654e+07, 
				   3.291e+07, 1.815e+07, 9.512e+06, 4.764e+06, 2.300e+06, 
				   1.081e+06, 5.020e+05, 2.337e+05, 1.111e+05, 5.483e+04, 
				   2.840e+04, 1.549e+04, 8.845e+03, 5.236e+03, 3.180e+03, 
				   1.964e+03, 1.225e+03, 7.678e+02, 4.813e+02, 3.006e+02, 
				   1.866e+02, 1.147e+02, 6.969e+01, 4.179e+01, 2.470e+01};
  std::vector<float> Data2012PUDist(Data2012PUDistArray, Data2012PUDistArray + 
				    sizeof(Data2012PUDistArray)/sizeof(float));
  if (PUScenario_ == "S7") {
    float S7PUDistArray[60] = {2.344E-05, 2.344E-05, 2.344E-05, 2.344E-05, 4.687E-04, 4.687E-04, 
			       7.032E-04, 9.414E-04, 1.234E-03, 1.603E-03, 2.464E-03, 3.250E-03, 
			       5.021E-03, 6.644E-03, 8.502E-03, 1.121E-02, 1.518E-02, 2.033E-02, 
			       2.608E-02, 3.171E-02, 3.667E-02, 4.060E-02, 4.338E-02, 4.520E-02, 
			       4.641E-02, 4.735E-02, 4.816E-02, 4.881E-02, 4.917E-02, 4.909E-02, 
			       4.842E-02, 4.707E-02, 4.501E-02, 4.228E-02, 3.896E-02, 3.521E-02, 
			       3.118E-02, 2.702E-02, 2.287E-02, 1.885E-02, 1.508E-02, 1.166E-02, 
			       8.673E-03, 6.190E-03, 4.222E-03, 2.746E-03, 1.698E-03, 9.971E-04, 
			       5.549E-04, 2.924E-04, 1.457E-04, 6.864E-05, 3.054E-05, 1.282E-05, 
			       5.081E-06, 1.898E-06, 6.688E-07, 2.221E-07, 6.947E-08, 2.047E-08};
    std::vector<float> S7PUDist(S7PUDistArray, S7PUDistArray + 
				sizeof(S7PUDistArray)/sizeof(float));
    PUReweight_ = edm::LumiReWeighting(S7PUDist, Data2012PUDist);
  }
  else if (PUScenario_ == "S10") {
    float S10PUDistArray[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 
				1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 
				2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 
				5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 
				4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 
				2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 
				1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 
				2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 
				5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 
				7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06};
    std::vector<float> S10PUDist(S10PUDistArray, S10PUDistArray + 
				 sizeof(S10PUDistArray)/sizeof(float));
    PUReweight_ = edm::LumiReWeighting(S10PUDist, Data2012PUDist);
  }
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

  //get W muons
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

  //get muons
  edm::Handle<reco::MuonCollection> pAllMuons;
  iEvent.getByLabel(allMuonTag_, pAllMuons);

  //get muon gen particles
  edm::Handle<reco::CandidateView> pMuonGenParticles;
  iEvent.getByLabel(muonGenParticleTag_, pMuonGenParticles);

  //get AK5 PF L1FastL2L3 jet correction service
  const JetCorrector* corrector = JetCorrector::getJetCorrector("ak5PFL1FastL2L3", iSetup);

  //get PU info
  edm::Handle<std::vector<PileupSummaryInfo> > pPU;
  iEvent.getByLabel(PUTag_, pPU);

  //get vertices
  edm::Handle<reco::VertexCollection> pVtx;
  iEvent.getByLabel(vtxTag_, pVtx);

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

//get PU weight
  double PUWeight = 1.0;
  float trueNInt = -1;
  if ((PUScenario_ == "S7") || (PUScenario_ == "S10")) {
    std::vector<PileupSummaryInfo>::const_iterator iPU = pPU->begin();
    int BX = 0;
    while ((iPU != pPU->end()) && (BX == 0)) {
      int BX = iPU->getBunchCrossing();
      if (BX == 0) { 
	trueNInt = iPU->getTrueNumInteractions();
	BX = -1;
      }
      ++iPU;
    }
    PUWeight = PUReweight_.weight(trueNInt);
  }

  //plot the number of good vertices
  nGoodVtx_->Fill(/*Common::numGoodVertices(pVtx)*/trueNInt, PUWeight);

  //plot MET distribution
  fillETHistogram(pMET, MET_, PUWeight);  

  //find the highest pT W muon
  std::vector<reco::MuonRef> WMuonRefs;
  for (reco::MuonRefVector::const_iterator iMuon = pMuons->begin(); iMuon != pMuons->end(); 
       ++iMuon) { WMuonRefs.push_back(*iMuon); }
  Common::sortByPT(WMuonRefs);

  //plot dPhi(highest pT W muon, MET)
  plotDPhiMuMet(WMuonRefs[WMuonRefs.size() - 1], pMET, dPhiWMuMET_, PUWeight);

  //plot the transverse mass for the highest pT W muon
  plotMT(WMuonRefs[WMuonRefs.size() - 1], pMET, WMuMT_, PUWeight);

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

//   //debug
//   edm::Handle<reco::MuonCollection> pAllMuons;
//   iEvent.getByLabel("muons", pAllMuons);
//   for (reco::MuonCollection::const_iterator iMu = pAllMuons->begin(); iMu != pAllMuons->end(); 
//        ++iMu) {
//     std::cerr << "Reco muon\n----------\n";
//     std::cerr << "pT = " << iMu->pt() << " GeV\n";
//     std::cerr << "eta = " << iMu->eta() << std::endl;
//     std::cerr << "phi = " << iMu->phi() << std::endl;
//   }

//fill an STL container of muon gen particles
  std::vector<reco::Candidate*> muGenParticlePtrs;
  for (unsigned int iMuGenParticle = 0; iMuGenParticle < pMuonGenParticles->size(); 
       ++iMuGenParticle) {
    muGenParticlePtrs.
      push_back(const_cast<reco::Candidate*>(pMuonGenParticles->refAt(iMuGenParticle).get()));
  }

  //sort selected taus by descending order in pT
  std::vector<reco::PFTauRef> pTSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { pTSortedTaus.push_back(*iTau); }
  std::vector<reco::PFTauRef> taus = pTSortedTaus;
  Common::sortByPT(pTSortedTaus);
  std::reverse(pTSortedTaus.begin(), pTSortedTaus.end());

  //sort selected taus by descending order in mu+had mass
  std::vector<reco::PFTauRef> muHadMassSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { muHadMassSortedTaus.push_back(*iTau); }
  Common::sortByMass(pMuonJetMap, muHadMassSortedTaus);

  //mu+had object counter
  unsigned int nMuHad = 0;

  //loop over selected taus
  std::vector<reco::PFTauRef>::const_iterator iTau = taus.begin();
  std::vector<reco::PFTauRef>::const_iterator endTau = taus.end();
  if (tauArbitrationMethod_ == "pT") {
    iTau = pTSortedTaus.begin();
    endTau = iTau + 1;
  }
  if (tauArbitrationMethod_ == "m") {
    iTau = muHadMassSortedTaus.begin();
    endTau = iTau + 1;
  }
  while (iTau != endTau) {
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
//     const reco::PFCandidateRefVector& PFChargedHadronIsoCands = 
//       (*iTau)->isolationPFCands();
    std::vector<reco::PFCandidate*> PFChargedHadronIsoCandPtrs;
    unsigned int nSoftMuIsoCands = 0;
    for (reco::PFCandidateRefVector::const_iterator iCand = PFChargedHadronIsoCands.begin(); 
	 iCand != PFChargedHadronIsoCands.end(); ++iCand) {
      PFChargedHadronIsoCandPtrs.push_back(const_cast<reco::PFCandidate*>((*iCand).get()));
      const reco::MuonRef isoCandMuRef = (*iCand)->muonRef();
      if (isoCandMuRef.isNonnull() && 
	  (isoCandMuRef.key() == removedMuonRefs[removedMuonRefs.size() - 1].key())) {
	++nSoftMuIsoCands;
// 	std::cerr << "Muon isolation candidate\n------------------------\n";
// 	std::cerr << "pT = " << (*iCand)->pt() << " GeV\n";
// 	std::cerr << "eta = " << (*iCand)->eta() << std::endl;
// 	std::cerr << "phi = " << (*iCand)->phi() << std::endl;
// 	std::cerr << "Product ID = (" << (*iCand).id().processIndex() << ", ";
// 	std::cerr << (*iCand).id().productIndex() << ")\n";
// 	std::cerr << "Muon isolation candidate muon ref\n---------------------------------\n";
// 	std::cerr << "pT = " << isoCandMuRef->pt() << " GeV\n";
// 	std::cerr << "eta = " << isoCandMuRef->eta() << std::endl;
// 	std::cerr << "phi = " << isoCandMuRef->phi() << std::endl;
// 	std::cerr << "Product ID = (" << isoCandMuRef.id().processIndex() << ", ";
// 	std::cerr << isoCandMuRef.id().productIndex() << ")\n";
      }
    }
    std::vector<double> sortedPFChargedHadronIsoCands = /*FIXME use min_element instead of sort 
							  for speed*/
      Common::sortByProximity(PFChargedHadronIsoCandPtrs, 
			      *removedMuonRefs[removedMuonRefs.size() - 1]);

    //debug
    int index = -1;
    const reco::PFCandidate* nearestIsoCand = 
      Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], 
			    PFChargedHadronIsoCandPtrs, index);
    double dR2 = -1.0;
    if (index >= 0) {
      dR2 = reco::deltaR2(*nearestIsoCand, *removedMuonRefs[removedMuonRefs.size() - 1]);
    }
    if ((sortedPFChargedHadronIsoCands.size() > 0) && 
	(sortedPFChargedHadronIsoCands[0] != dR2)) { //sortByProximity returns dR2, not dR!
      std::cerr << "Error: " << sortedPFChargedHadronIsoCands[0] << " != " << dR2 << std::endl;
    }

    //impose pT and decay mode cut on tau
    if (((*iTau)->pt() > tauPTMin_) && 
	((tauDecayMode_ == reco::PFTau::kNull) || ((*iTau)->decayMode() == tauDecayMode_))//  && 
// 	!((sortedPFChargedHadronIsoCands.size() > 0) && 
// 	  (sortedPFChargedHadronIsoCands[0] < dR_))
	) {

      //plot multiplicity of muons associated to this tau
      hadTauAssociatedMuMultiplicity_->Fill(removedMuons.size(), PUWeight);

      //plot the mu + tau invariant mass for the highest pT muon
      const double muHadMass = 
	(removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();
      muHadMass_->Fill(muHadMass, PUWeight);

      //plot the mu + tau charge for the highest pT muon
      muHadCharge_->
	Fill(removedMuonRefs[removedMuonRefs.size() - 1]->charge() + (*iTau)->charge(), PUWeight);

      //plot dPhi(highest pT muon, MET)
      plotDPhiMuMet(removedMuonRefs[removedMuonRefs.size() - 1], pMET, dPhiTauMuMET_, PUWeight);

      //plot the transverse mass for the highest pT muon
      plotMT(removedMuonRefs[removedMuonRefs.size() - 1], pMET, tauMuMT_, PUWeight);

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
      plotHT(tauMuTauHadJet, tauMuTauHadJetHT_, PUWeight);

      //plot HT (tau muon + hadronic tau + leading distinct corrected jet + W muon)
      tauMuTauHadJet.push_back(dynamic_cast<reco::Candidate*>
			       (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(tauMuTauHadJet, tauMuTauHadJetWMuHT_, PUWeight);

      //plot HT (two leading corrected jets if both exist)
      std::vector<reco::Candidate*> diJet;
      diJet.push_back(dynamic_cast<reco::Candidate*>
		      (const_cast<reco::PFJet*>(oldJetRefs[oldJetRefs.size() - 1].get())));
      if (oldJetRefs.size() > 1) {
	diJet.push_back(dynamic_cast<reco::Candidate*>
			(const_cast<reco::PFJet*>(oldJetRefs[oldJetRefs.size() - 2].get())));
	plotHT(diJet, diJetHT_, PUWeight);
      }

      //plot HT (two leading corrected jets + W muon)
      diJet.push_back(dynamic_cast<reco::Candidate*>
		      (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(diJet, diJetWMuHT_, PUWeight);

      //plot HT (corrected jet associated to hadronic tau + leading distinct corrected jet)
      std::vector<reco::Candidate*> jetTauJet;
      jetTauJet.push_back(dynamic_cast<reco::Candidate*>
			  (const_cast<reco::PFJet*>(&correctedTauJet)));
      jetTauJet.push_back(dynamic_cast<reco::Candidate*>
			  (const_cast<reco::PFJet*>
			   (oldJetRefsExcludingTau[oldJetRefsExcludingTau.size() - 1].get())));
      plotHT(jetTauJet, jetTauJetHT_, PUWeight);

      /*plot HT (corrected jet associated to hadronic tau + leading distinct corrected jet + W 
	muon)*/
      jetTauJet.push_back(dynamic_cast<reco::Candidate*>
			  (const_cast<reco::Muon*>(WMuonRefs[WMuonRefs.size() - 1].get())));
      plotHT(jetTauJet, jetTauJetWMuHT_, PUWeight);

      //plot the parent parton of the tau jet
      int nearestPartonIndex = -1;
      const reco::GenParticle* nearestParton = 
	Common::nearestObject(tauOldJetRef, genParticlePtrs, nearestPartonIndex);
      if (nearestParton != NULL) {
	if (reco::deltaR(*nearestParton, *tauOldJetRef) < dR_) {
	  const unsigned int absNearestPartonPDGID = fabs(nearestParton->pdgId());
	  unsigned int bin = absNearestPartonPDGID == GenTauDecayID::G ? 0 : absNearestPartonPDGID;
	  jetParentParton_->Fill(bin, PUWeight);
	}
	else jetParentParton_->Fill(-1, PUWeight);
      }
      else jetParentParton_->Fill(7, PUWeight);

      //plot dR(W muon, leading soft muon per jet) when mu+had mass > 2 GeV
      if (muHadMass > 2.0/*GeV*/) {
	dRWMuSoftMuMuHadMassGe2_->Fill(reco::deltaR(*WMuonRefs[WMuonRefs.size() - 1], 
						    *removedMuonRefs[removedMuonRefs.size() - 1]), 
				       PUWeight);
      }

      //plot tau muon pT
      const double tauMuPT = removedMuonRefs[removedMuonRefs.size() - 1]->pt();
      tauMuPT_->Fill(tauMuPT, PUWeight);

      //plot hadronic tau pT
      tauHadPT_->Fill((*iTau)->pt(), PUWeight);

      //plot hadronic tau deltaBeta-corrected isolation energy
      tauHadIso_->Fill((*pTauHadIso)[*iTau], PUWeight);

      //plot soft muon isolation candidate multiplicity per removed muon
//       softMuIsoCandMultiplicity_->Fill((sortedPFChargedHadronIsoCands.size() > 0) && 
// 				       (sortedPFChargedHadronIsoCands[0] < dR_*dR_));
      softMuIsoCandMultiplicity_->Fill(nSoftMuIsoCands, PUWeight);

      //find nearest muon to soft muon
      int nearestMuIndex = -1;
      std::vector<reco::Muon*> muPtrs;
      for (reco::MuonCollection::const_iterator iMu = pAllMuons->begin(); iMu != pAllMuons->end(); 
	   ++iMu) {
	if ((iMu - pAllMuons->begin()) != removedMuonRefs[removedMuonRefs.size() - 1].key()) {
	  muPtrs.push_back(const_cast<reco::Muon*>(&*iMu));
	}
      }
      const reco::Muon* nearestMu = 
	Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], muPtrs, nearestMuIndex);

      //find nearest gen muon to soft muon
      int nearestGenMuIndex = -1;
      const reco::Candidate* nearestGenMu = 
	Common::nearestObject(removedMuonRefs[removedMuonRefs.size() - 1], muGenParticlePtrs, 
			      nearestGenMuIndex);

      //plot soft muon nearest muon properties
      double dRSoftMuNearestMu = -1.0;
      double dRSoftMuNearestGenMu = -1.0;
      int xBinVal = -1;
      int yBinVal = -1;
      if (nearestMuIndex != -1) {
	dRSoftMuNearestMu = reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], *nearestMu);
      }
      if (nearestGenMuIndex != -1) {
	dRSoftMuNearestGenMu = 
	  reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], *nearestGenMu);
      }
      if ((dRSoftMuNearestMu >= 0.0) && (dRSoftMuNearestMu < dR_)) {
	if (nearestMuIndex == (int)WMuonRefs[WMuonRefs.size() - 1].key()) {
	  xBinVal = 0;
	  if ((dRSoftMuNearestGenMu >= 0.0) && (dRSoftMuNearestGenMu < dR_)) yBinVal = 0;
	  else yBinVal = 1;
	}
	else {
	  xBinVal = 1;
	  if ((dRSoftMuNearestGenMu >= 0.0) && (dRSoftMuNearestGenMu < dR_)) yBinVal = 0;
	  else yBinVal = 1;
	}
      }
      else {
	xBinVal = 2;
	if ((dRSoftMuNearestGenMu >= 0.0) && (dRSoftMuNearestGenMu < dR_)) yBinVal = 0;
	else yBinVal = 1;
      }
      genMuExistsVsSoftMuNearestMuProperties_->Fill(xBinVal, yBinVal, PUWeight);

      //plot hadronic tau eta
      tauHadEta_->Fill((*iTau)->eta(), PUWeight);

      //plot pT(soft muon)/(mu+had mass)
      softMuPTOverMuHadMass_->Fill(removedMuonRefs[removedMuonRefs.size() - 1]->pt()/muHadMass, 
				   PUWeight);

      //plot pT(mu+had)/(mu+had mass) for pT(mu+had) > 30 GeV
      double muHadPT = (removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).pt();
//       if (muHadPT > 30.0/*GeV*/) muHadPTOverMuHadMass_->Fill(muHadPT/muHadMass);
      if (tauOldJetRef->pt() > uncorrJetPTMin_) {
	muHadPTOverMuHadMass_->Fill(tauOldJetRef->pt()/tauOldJetRef->mass(), PUWeight);
      }

      //plot dR(soft muon, nearest gen muon)
      dRSoftMuNearestGenMuHist_->Fill(dRSoftMuNearestGenMu, PUWeight);

      //plot mu+had pT
      muHadPT_->Fill(muHadPT, PUWeight);

      //plot cleaned jet pT vs. cleaned tau pT
      cleanedJetPTVsCleanedTauPT_->Fill((*iTau)->pt(), tauJetRef->pt(), PUWeight);

      //plot uncleaned jet pT vs. cleaned tau pT
      uncleanedJetPTVsCleanedTauPT_->Fill((*iTau)->pt(), tauOldJetRef->pt(), PUWeight);

      //plot mu+had mass vs. dR(tagged soft muon, tau axis)
      const double dRSoftMuTau = 
	reco::deltaR(*removedMuonRefs[removedMuonRefs.size() - 1], **iTau);
      muHadMassVsDRSoftMuTau_->Fill(dRSoftMuTau, muHadMass, PUWeight);

      //plot mu+had mass vs. tagged soft muon pT
      muHadMassVsSoftMuPT_->Fill(tauMuPT, muHadMass, PUWeight);

      //plot hadronic tau isolation vs. soft muon pT
      tauHadIsoVsSoftMuPT_->Fill(tauMuPT, (*pTauHadIso)[*iTau], PUWeight);

      //plot mu+had mass vs. hadronic tau eta
      muHadMassVsTauHadEta_->Fill((*iTau)->eta(), muHadMass, PUWeight);

      //plot mu+had mass vs. soft muon eta
      muHadMassVsSoftMuEta_->Fill(removedMuonRefs[removedMuonRefs.size() - 1]->eta(), muHadMass, 
				  PUWeight);

      //plot mu+had mass vs. hadronic tau isolation
      muHadMassVsTauHadIso_->Fill((*pTauHadIso)[*iTau], muHadMass, PUWeight);

      //plot mu+had mass vs. hadronic tau pT
      muHadMassVsTauHadPT_->Fill((*iTau)->pt(), muHadMass, PUWeight);

      //plot hadronic tau isolation vs. eta
      tauHadIsoVsEta_->Fill((*iTau)->eta(), (*pTauHadIso)[*iTau], PUWeight);

      //plot hadronic tau eta vs. soft muon eta (sanity check)
      tauHadEtaVsSoftMuEta_->Fill(removedMuonRefs[removedMuonRefs.size() - 1]->eta(), 
				  (*iTau)->eta(), PUWeight);

      //plot histogram of dEta(hadronic tau, soft muon) vs. dPhi(hadronic tau, soft muon)
      dEtaTauHadSoftMuVsDPhiTauHadSoftMu_->
	Fill(fabs(reco::deltaPhi(**iTau, *removedMuonRefs[removedMuonRefs.size() - 1])), 
	     fabs((*iTau)->eta() - removedMuonRefs[removedMuonRefs.size() - 1]->eta()), PUWeight);

      //plot histogram of pT(hadronic tau)/(mu+had mass) vs. hadronic tau isolation
      tauHadPTOverMuHadMassVsTauHadIso_->Fill((*pTauHadIso)[*iTau], (*iTau)->pt()/muHadMass, 
					      PUWeight);

      //plot histogram of pT(soft muon)/(mu+had mass) vs. hadronic tau isolation
      softMuPTOverMuHadMassVsTauHadIso_->
	Fill((*pTauHadIso)[*iTau], removedMuonRefs[removedMuonRefs.size() - 1]->pt()/muHadMass, 
	     PUWeight);

      /*plot histogram of (pT(soft muon) + pT(hadronic tau)/(2(mu+had mass)) vs. hadronic tau 
	isolation*/
      avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_->
	Fill((*pTauHadIso)[*iTau], 
	     (removedMuonRefs[removedMuonRefs.size() - 1]->pt() + (*iTau)->pt())/(2*muHadMass), 
	     PUWeight);

      /*plot histogram of pT(mu+had)/(mu+had mass) vs. hadronic tau isolation for pT(mu+had) > 30 
	GeV*/
      /*this plot is currently inconsistent with the 1D plot of pT(mu+had)/(mu+had mass), which 
	uses the parent uncleaned uncorrected jet pT and mass*/
      if (muHadPT > 30.0/*GeV*/) {
	muHadPTOverMuHadMassVsTauHadIso_->Fill((*pTauHadIso)[*iTau], muHadPT/muHadMass, PUWeight);
      }

      //plot soft muon pT vs. hadronic tau pT
      softMuPTVsTauHadPT_->Fill((*iTau)->pt(), removedMuonRefs[removedMuonRefs.size() - 1]->pt(), 
				PUWeight);

      //     //plot type of gen tau match
      //     if (std::find(genMatchedMuonRefs.begin(), genMatchedMuonRefs.end(), 
      // 		  removedMuonRefs[removedMuonRefs.size() - 1].key()) != genMatchedMuonRefs.end()) {
      //       genTauMatchType_->Fill();
      //     }

      //increment mu+had multiplicity counter
      ++nMuHad;
    }

    //increment tau iterator
    ++iTau;
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
      dRWMuSoftMu_->Fill(dR, PUWeight);

      //plot dR(W muon, leading soft gen-matched muon per jet)
      if (std::find(genMatchedMuonRefs.begin(), genMatchedMuonRefs.end(), 
		    removedMuonRefs[removedMuonRefs.size() - 1].key()) != 
	  genMatchedMuonRefs.end()) dRWMuSoftGenMatchedMu_->Fill(dR, PUWeight);
    }
  }

  //plot mu+had multiplicity
  muHadMultiplicity_->Fill(nMuHad, PUWeight);
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
  tauHadEta_ = new TH1F("tauHadEta", ";#eta;", 46, -2.3, 2.3);
  softMuPTOverMuHadMass_ = 
    new TH1F("softMuPTOverMuHadMass", ";#frac{p_{T}^{#mu}}{m_{#mu+had}};", 80, 0.0, 80.0);
  muHadPTOverMuHadMass_ = 
    new TH1F("muHadPTOverMuHadMass", ";#frac{p_{T}^{#mu+had}}{m_{#mu+had}};", 80, 0.0, 80.0);
  dRSoftMuNearestGenMuHist_ = 
    new TH1F("dRSoftMuNearestGenMuHist", ";#DeltaR(soft #mu, gen #mu);", 40, 0.0, 2.0);
  muHadPT_ = new TH1F("muHadPT", ";p_{T}^{#mu+had} (GeV);", 50, 0.0, 100.0);
  muHadMultiplicity_ = new TH1F("muHadMultiplicity", ";N_{#mu+had}/event;", 4, 0.5, 4.5);
  nGoodVtx_ = new TH1F("nGoodVtx", ";No. good vertices;", 60, -0.5, 59.5);
  cleanedJetPTVsCleanedTauPT_ = 
    new TH2F("cleanedJetPTVsCleanedTauPT", ";#tau p_{T} (GeV);Jet p_{T} (GeV)", 
	     50, 0.0, 100.0, 50, 0.0, 100.0);
  uncleanedJetPTVsCleanedTauPT_ = 
    new TH2F("uncleanedJetPTVsCleanedTauPT", ";#tau p_{T} (GeV);Jet p_{T} (GeV)", 
	     50, 0.0, 100.0, 50, 0.0, 100.0);
  muHadMassVsDRSoftMuTau_ = 
    new TH2F("muHadMassVsDRSoftMuTau", ";#DeltaR(soft muon, #tau_{had});m_{#mu+had} (GeV)", 
	     30, 0.0, 3.0, 20, 0.0, 20.0);
  muHadMassVsSoftMuPT_ = 
    new TH2F("muHadMassVsSoftMuPT", ";p_{T} (GeV);m_{#mu+had} (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 20.0);
  tauHadIsoVsSoftMuPT_ = 
    new TH2F("tauHadIsoVsSoftMuPT", ";p_{T} (GeV);Isolation energy (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 20.0);
  genMuExistsVsSoftMuNearestMuProperties_ = 
    new TH2F("genMuExistsVsSoftMuNearestMuProperties", ";Nearest muon property;Gen muon property", 
	     3, -0.5, 2.5, 2, -0.5, 1.5);
  muHadMassVsTauHadEta_ = 
    new TH2F("muHadMassVsTauHadEta", ";#eta;m_{#mu+had} (GeV)", 46, -2.3, 2.3, 20, 0.0, 20.0);
  muHadMassVsSoftMuEta_ = 
    new TH2F("muHadMassVsSoftMuEta", ";#eta;m_{#mu+had} (GeV)", 46, -2.3, 2.3, 20, 0.0, 20.0);
  muHadMassVsTauHadIso_ = 
    new TH2F("muHadMassVsTauHadIso", ";Isolation energy (GeV);m_{#mu+had} (GeV)", 
	     20, 0.0, 20.0, 20, 0.0, 20.0);
  muHadMassVsTauHadPT_ = 
    new TH2F("muHadMassVsTauHadPT", ";p_{T} (GeV);m_{#mu+had} (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 20.0);
  tauHadIsoVsEta_ = 
    new TH2F("tauHadIsoVsEta", ";#eta;Isolation energy (GeV)", 46, -2.3, 2.3, 20, 0.0, 20.0);
  tauHadEtaVsSoftMuEta_ = 
    new TH2F("tauHadEtaVsSoftMuEta", ";#eta_{#mu};#eta_{#tau}", 46, -2.3, 2.3, 46, -2.3, 2.3);
  dEtaTauHadSoftMuVsDPhiTauHadSoftMu_ = 
    new TH2F("dEtaTauHadSoftMuVsDPhiTauHadSoftMu", ";#Delta#phi(#tau, #mu);#Delta#eta(#tau, #mu)", 
	     20, 0.0, 1.0, 20, 0.0, 1.0);
  tauHadPTOverMuHadMassVsTauHadIso_ = 
    new TH2F("tauHadPTOverMuHadMassVsTauHadIso", 
	     ";Isolation energy (GeV);#frac{p_{T}^{#tau}}{m_{#mu+had}}", 
	     20, 0.0, 20.0, 80, 0.0, 80.0);
  softMuPTOverMuHadMassVsTauHadIso_ = 
    new TH2F("softMuPTOverMuHadMassVsTauHadIso", 
	     ";Isolation energy (GeV);#frac{p_{T}^{#mu}}{m_{#mu+had}}", 
	     20, 0.0, 20.0, 80, 0.0, 80.0);
  avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_ = 
    new TH2F("avgTauHadSoftMuPTOverMuHadMassVsTauHadIso", 
	     ";Isolation energy (GeV);#frac{p_{T}^{#tau} + p_{T}^{#mu}}{2m_{#mu+had}}", 
	     20, 0.0, 20.0, 80, 0.0, 80.0);
  muHadPTOverMuHadMassVsTauHadIso_ = 
    new TH2F("muHadPTOverMuHadMassVsTauHadIso", 
	     ";Isolation energy (GeV);#frac{p_{T}^{#mu+had}}{m_{#mu+had}}", 
	     20, 0.0, 20.0, 80, 0.0, 80.0);
  softMuPTVsTauHadPT_ = 
    new TH2F("softMuPTVsTauHadPT", ";p_{T}^{#tau} (GeV);p_{T}^{#mu} (GeV)", 
	     20, 0.0, 100.0, 20, 0.0, 100.0);

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
  std::stringstream dRStream;
  dRStream << dR_;
  genMuExistsVsSoftMuNearestMuProperties_->
    GetXaxis()->SetBinLabel(1, ("#DeltaR < " + dRStream.str() + " and W muon").c_str());
  genMuExistsVsSoftMuNearestMuProperties_->
    GetXaxis()->SetBinLabel(2, ("#DeltaR < " + dRStream.str() + " and not W muon").c_str());
  genMuExistsVsSoftMuNearestMuProperties_->
    GetXaxis()->SetBinLabel(3, ("#DeltaR #geq " + dRStream.str()).c_str());
  genMuExistsVsSoftMuNearestMuProperties_->
    GetYaxis()->SetBinLabel(1, ("#DeltaR < " + dRStream.str()).c_str());
  genMuExistsVsSoftMuNearestMuProperties_->
    GetYaxis()->SetBinLabel(2, ("#DeltaR #geq " + dRStream.str()).c_str());

  //set sumw2
  MET_->Sumw2();
  hadTauAssociatedMuMultiplicity_->Sumw2();
  muHadMass_->Sumw2();
  muHadCharge_->Sumw2();
  WMuMT_->Sumw2();
  tauMuMT_->Sumw2();
  dPhiWMuMET_->Sumw2();
  dPhiTauMuMET_->Sumw2();
  tauMuTauHadJetHT_->Sumw2();
  diJetHT_->Sumw2();
  jetTauJetHT_->Sumw2();
  tauMuTauHadJetWMuHT_->Sumw2();
  diJetWMuHT_->Sumw2();
  jetTauJetWMuHT_->Sumw2();
  jetParentParton_->Sumw2();
  dRWMuSoftMu_->Sumw2();
  dRWMuSoftGenMatchedMu_->Sumw2();
  dRWMuSoftMuMuHadMassGe2_->Sumw2();
  genTauMatchType_->Sumw2();
  tauMuPT_->Sumw2();
  tauHadPT_->Sumw2();
  tauHadIso_->Sumw2();
  softMuIsoCandMultiplicity_->Sumw2();
  tauHadEta_->Sumw2();
  softMuPTOverMuHadMass_->Sumw2();
  muHadPTOverMuHadMass_->Sumw2();
  dRSoftMuNearestGenMuHist_->Sumw2();
  muHadPT_->Sumw2();
  muHadMultiplicity_->Sumw2();
  nGoodVtx_->Sumw2();
  cleanedJetPTVsCleanedTauPT_->Sumw2();
  uncleanedJetPTVsCleanedTauPT_->Sumw2();
  muHadMassVsDRSoftMuTau_->Sumw2();
  muHadMassVsSoftMuPT_->Sumw2();
  tauHadIsoVsSoftMuPT_->Sumw2();
  genMuExistsVsSoftMuNearestMuProperties_->Sumw2();
  muHadMassVsTauHadEta_->Sumw2();
  muHadMassVsSoftMuEta_->Sumw2();
  muHadMassVsTauHadIso_->Sumw2();
  muHadMassVsTauHadPT_->Sumw2();
  tauHadIsoVsEta_->Sumw2();
  tauHadEtaVsSoftMuEta_->Sumw2();
  dEtaTauHadSoftMuVsDPhiTauHadSoftMu_->Sumw2();
  tauHadPTOverMuHadMassVsTauHadIso_->Sumw2();
  softMuPTOverMuHadMassVsTauHadIso_->Sumw2();
  avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_->Sumw2();
  muHadPTOverMuHadMassVsTauHadIso_->Sumw2();
  softMuPTVsTauHadPT_->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void TauAnalyzer::endJob() 
{
  //make canvases
  out_->cd();
  TCanvas METCanvas("METCanvas", "", 600, 600);
  TCanvas hadTauAssociatedMuMultiplicityCanvas("hadTauAssociatedMuMultiplicityCanvas", "", 
					       600, 600);
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
  TCanvas tauHadEtaCanvas("tauHadEtaCanvas", "", 600, 600);
  TCanvas softMuPTOverMuHadMassCanvas("softMuPTOverMuHadMassCanvas", "", 600, 600);
  TCanvas muHadPTOverMuHadMassCanvas("muHadPTOverMuHadMassCanvas", "", 600, 600);
  TCanvas dRSoftMuNearestGenMuHistCanvas("dRSoftMuNearestGenMuHistCanvas", "", 600, 600);
  TCanvas muHadPTCanvas("muHadPTCanvas", "", 600, 600);
  TCanvas muHadMultiplicityCanvas("muHadMultiplicityCanvas", "", 600, 600);
  TCanvas nGoodVtxCanvas("nGoodVtxCanvas", "", 600, 600);
  TCanvas cleanedJetPTVsCleanedTauPTCanvas("cleanedJetPTVsCleanedTauPTCanvas", "", 600, 600);
  TCanvas uncleanedJetPTVsCleanedTauPTCanvas("uncleanedJetPTVsCleanedTauPTCanvas", "", 600, 600);
  TCanvas muHadMassVsDRSoftMuTauCanvas("muHadMassVsDRSoftMuTauCanvas", "", 600, 600);
  TCanvas muHadMassVsSoftMuPTCanvas("muHadMassVsSoftMuPTCanvas", "", 600, 600);
  TCanvas tauHadIsoVsSoftMuPTCanvas("tauHadIsoVsSoftMuPTCanvas", "", 600, 600);
  TCanvas 
    genMuExistsVsSoftMuNearestMuPropertiesCanvas("genMuExistsVsSoftMuNearestMuPropertiesCanvas", 
						 "", 600, 600);
  TCanvas muHadMassVsTauHadEtaCanvas("muHadMassVsTauHadEtaCanvas", "", 600, 600);
  TCanvas muHadMassVsSoftMuEtaCanvas("muHadMassVsSoftMuEtaCanvas", "", 600, 600);
  TCanvas muHadMassVsTauHadIsoCanvas("muHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas muHadMassVsTauHadPTCanvas("muHadMassVsTauHadPTCanvas", "", 600, 600);
  TCanvas tauHadIsoVsEtaCanvas("tauHadIsoVsEtaCanvas", "", 600, 600);
  TCanvas tauHadEtaVsSoftMuEtaCanvas("tauHadEtaVsSoftMuEtaCanvas", "", 600, 600);
  TCanvas 
    dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas("dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas", "", 
					     600, 600);
  TCanvas 
    tauHadPTOverMuHadMassVsTauHadIsoCanvas("tauHadPTOverMuHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas 
    softMuPTOverMuHadMassVsTauHadIsoCanvas("softMuPTOverMuHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas 
    avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas("avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas 
    muHadPTOverMuHadMassVsTauHadIsoCanvas("muHadPTOverMuHadMassVsTauHadIsoCanvas", "", 600, 600);
  TCanvas softMuPTVsTauHadPTCanvas("softMuPTVsTauHadPTCanvas", "", 600, 600);

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
  Common::draw1DHistograms(tauHadEtaCanvas, tauHadEta_);
  Common::draw1DHistograms(softMuPTOverMuHadMassCanvas, softMuPTOverMuHadMass_);
  Common::draw1DHistograms(muHadPTOverMuHadMassCanvas, muHadPTOverMuHadMass_);
  Common::draw1DHistograms(dRSoftMuNearestGenMuHistCanvas, dRSoftMuNearestGenMuHist_);
  Common::draw1DHistograms(muHadPTCanvas, muHadPT_);
  Common::draw1DHistograms(muHadMultiplicityCanvas, muHadMultiplicity_);
  Common::draw1DHistograms(nGoodVtxCanvas, nGoodVtx_);

  //format and draw 2D plots
  Common::draw2DHistograms(cleanedJetPTVsCleanedTauPTCanvas, cleanedJetPTVsCleanedTauPT_);
  Common::draw2DHistograms(uncleanedJetPTVsCleanedTauPTCanvas, uncleanedJetPTVsCleanedTauPT_);
  Common::draw2DHistograms(muHadMassVsDRSoftMuTauCanvas, muHadMassVsDRSoftMuTau_);
  Common::draw2DHistograms(muHadMassVsSoftMuPTCanvas, muHadMassVsSoftMuPT_);
  Common::draw2DHistograms(tauHadIsoVsSoftMuPTCanvas, tauHadIsoVsSoftMuPT_);
  Common::draw2DHistograms(genMuExistsVsSoftMuNearestMuPropertiesCanvas, 
			   genMuExistsVsSoftMuNearestMuProperties_);
  Common::draw2DHistograms(muHadMassVsTauHadEtaCanvas, muHadMassVsTauHadEta_);
  Common::draw2DHistograms(muHadMassVsSoftMuEtaCanvas, muHadMassVsSoftMuEta_);
  Common::draw2DHistograms(muHadMassVsTauHadIsoCanvas, muHadMassVsTauHadIso_);
  Common::draw2DHistograms(muHadMassVsTauHadPTCanvas, muHadMassVsTauHadPT_);
  Common::draw2DHistograms(tauHadIsoVsEtaCanvas, tauHadIsoVsEta_);
  Common::draw2DHistograms(tauHadEtaVsSoftMuEtaCanvas, tauHadEtaVsSoftMuEta_);
  Common::draw2DHistograms(dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas, 
			   dEtaTauHadSoftMuVsDPhiTauHadSoftMu_);
  Common::draw2DHistograms(tauHadPTOverMuHadMassVsTauHadIsoCanvas, 
			   tauHadPTOverMuHadMassVsTauHadIso_);
  Common::draw2DHistograms(softMuPTOverMuHadMassVsTauHadIsoCanvas, 
			   softMuPTOverMuHadMassVsTauHadIso_);
  Common::draw2DHistograms(avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas, 
			   avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_);
  Common::draw2DHistograms(muHadPTOverMuHadMassVsTauHadIsoCanvas, 
			   muHadPTOverMuHadMassVsTauHadIso_);
  Common::draw2DHistograms(softMuPTVsTauHadPTCanvas, softMuPTVsTauHadPT_);

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
  tauHadEtaCanvas.Write();
  softMuPTOverMuHadMassCanvas.Write();
  muHadPTOverMuHadMassCanvas.Write();
  dRSoftMuNearestGenMuHistCanvas.Write();
  muHadPTCanvas.Write();
  muHadMultiplicityCanvas.Write();
  nGoodVtxCanvas.Write();
  cleanedJetPTVsCleanedTauPTCanvas.Write();
  uncleanedJetPTVsCleanedTauPTCanvas.Write();
  muHadMassVsDRSoftMuTauCanvas.Write();
  muHadMassVsSoftMuPTCanvas.Write();
  tauHadIsoVsSoftMuPTCanvas.Write();
  genMuExistsVsSoftMuNearestMuPropertiesCanvas.Write();
  muHadMassVsTauHadEtaCanvas.Write();
  muHadMassVsSoftMuEtaCanvas.Write();
  muHadMassVsTauHadIsoCanvas.Write();
  muHadMassVsTauHadPTCanvas.Write();
  tauHadIsoVsEtaCanvas.Write();
  tauHadEtaVsSoftMuEtaCanvas.Write();
  dEtaTauHadSoftMuVsDPhiTauHadSoftMuCanvas.Write();
  tauHadPTOverMuHadMassVsTauHadIsoCanvas.Write();
  softMuPTOverMuHadMassVsTauHadIsoCanvas.Write();
  avgTauHadSoftMuPTOverMuHadMassVsTauHadIsoCanvas.Write();
  muHadPTOverMuHadMassVsTauHadIsoCanvas.Write();
  softMuPTVsTauHadPTCanvas.Write();
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
				const edm::Handle<edm::View<reco::PFMET> >& pMET, TH1F* hist, 
				const double weight)
{
  hist->Fill(fabs(reco::deltaPhi(muonRef->phi(), pMET->refAt(0)->phi())), weight);
}

void TauAnalyzer::plotMT(const reco::MuonRef& muonRef, 
			 const edm::Handle<edm::View<reco::PFMET> >& pMET, TH1F* hist, 
			 const double weight)
{
  edm::RefToBase<reco::PFMET> METRefToBase = pMET->refAt(0);
  hist->Fill(sqrt(2*muonRef->pt()*METRefToBase->et()*
		  (1.0 - cos(reco::deltaPhi(muonRef->phi(), METRefToBase->phi())))), weight);
}

void TauAnalyzer::plotHT(const std::vector<reco::Candidate*>& cands, TH1F* hist, 
			 const double weight)
{
  double HT = 0.0;
  for (std::vector<reco::Candidate*>::const_iterator iCand = cands.begin(); iCand != cands.end(); 
       ++iCand) {
    HT+=(*iCand)->pt();
  }
  hist->Fill(HT, weight);
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
  if (doDelete && (tauHadEta_ != NULL)) delete tauHadEta_;
  tauHadEta_ = NULL;
  if (doDelete && (softMuPTOverMuHadMass_ != NULL)) delete softMuPTOverMuHadMass_;
  softMuPTOverMuHadMass_ = NULL;
  if (doDelete && (muHadPTOverMuHadMass_ != NULL)) delete muHadPTOverMuHadMass_;
  muHadPTOverMuHadMass_ = NULL;
  if (doDelete && (dRSoftMuNearestGenMuHist_ != NULL)) delete dRSoftMuNearestGenMuHist_;
  dRSoftMuNearestGenMuHist_ = NULL;
  if (doDelete && (muHadPT_ != NULL)) delete muHadPT_;
  muHadPT_ = NULL;
  if (doDelete && (muHadMultiplicity_ != NULL)) delete muHadMultiplicity_;
  muHadMultiplicity_ = NULL;
  if (doDelete && (nGoodVtx_ != NULL)) delete nGoodVtx_;
  nGoodVtx_ = NULL;
  if (doDelete && (cleanedJetPTVsCleanedTauPT_ != NULL)) delete cleanedJetPTVsCleanedTauPT_;
  cleanedJetPTVsCleanedTauPT_ = NULL;
  if (doDelete && (uncleanedJetPTVsCleanedTauPT_ != NULL)) delete uncleanedJetPTVsCleanedTauPT_;
  uncleanedJetPTVsCleanedTauPT_ = NULL;
  if (doDelete && (muHadMassVsDRSoftMuTau_ != NULL)) delete muHadMassVsDRSoftMuTau_;
  muHadMassVsDRSoftMuTau_ = NULL;
  if (doDelete && (muHadMassVsSoftMuPT_ != NULL)) delete muHadMassVsSoftMuPT_;
  muHadMassVsSoftMuPT_ = NULL;
  if (doDelete && (tauHadIsoVsSoftMuPT_ != NULL)) delete tauHadIsoVsSoftMuPT_;
  tauHadIsoVsSoftMuPT_ = NULL;
  if (doDelete && (genMuExistsVsSoftMuNearestMuProperties_ != NULL)) {
    delete genMuExistsVsSoftMuNearestMuProperties_;
  }
  genMuExistsVsSoftMuNearestMuProperties_ = NULL;
  if (doDelete && (muHadMassVsTauHadEta_ != NULL)) delete muHadMassVsTauHadEta_;
  muHadMassVsTauHadEta_ = NULL;
  if (doDelete && (muHadMassVsSoftMuEta_ != NULL)) delete muHadMassVsSoftMuEta_;
  muHadMassVsSoftMuEta_ = NULL;
  if (doDelete && (muHadMassVsTauHadIso_ != NULL)) delete muHadMassVsTauHadIso_;
  muHadMassVsTauHadIso_ = NULL;
  if (doDelete && (muHadMassVsTauHadPT_ != NULL)) delete muHadMassVsTauHadPT_;
  muHadMassVsTauHadPT_ = NULL;
  if (doDelete && (tauHadIsoVsEta_ != NULL)) delete tauHadIsoVsEta_;
  tauHadIsoVsEta_ = NULL;
  if (doDelete && (tauHadEtaVsSoftMuEta_ != NULL)) delete tauHadEtaVsSoftMuEta_;
  tauHadEtaVsSoftMuEta_ = NULL;
  if (doDelete && (dEtaTauHadSoftMuVsDPhiTauHadSoftMu_ != NULL)) {
    delete dEtaTauHadSoftMuVsDPhiTauHadSoftMu_;
  }
  dEtaTauHadSoftMuVsDPhiTauHadSoftMu_ = NULL;
  if (doDelete && (tauHadPTOverMuHadMassVsTauHadIso_ != NULL)) {
    delete tauHadPTOverMuHadMassVsTauHadIso_;
  }
  tauHadPTOverMuHadMassVsTauHadIso_ = NULL;
  if (doDelete && (softMuPTOverMuHadMassVsTauHadIso_ != NULL)) {
    delete softMuPTOverMuHadMassVsTauHadIso_;
  }
  softMuPTOverMuHadMassVsTauHadIso_ = NULL;
  if (doDelete && (avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_ != NULL)) {
    delete avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_;
  }
  avgTauHadSoftMuPTOverMuHadMassVsTauHadIso_ = NULL;
  if (doDelete && (muHadPTOverMuHadMassVsTauHadIso_ != NULL)) {
    delete muHadPTOverMuHadMassVsTauHadIso_;
  }
  muHadPTOverMuHadMassVsTauHadIso_ = NULL;
  if (doDelete && (softMuPTVsTauHadPT_ != NULL)) delete softMuPTVsTauHadPT_;
  softMuPTVsTauHadPT_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
