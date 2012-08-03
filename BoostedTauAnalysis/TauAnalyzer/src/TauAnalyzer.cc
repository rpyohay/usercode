// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      TauAnalyzer
// 
/**\class TauAnalyzer TauAnalyzer.cc BoostedTauAnalysis/TauAnalyzer/src/TauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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

      void reset();
      double dR(const std::vector<double>& dR2, const unsigned int pos = 0) const;
      void fillDRHistogramNearFar(const unsigned int, const bool, 
				  const std::vector<double>&) const;
      void fillDRHistogramLeadTrail(const bool, const unsigned int, 
				    const std::vector<double>&) const;
      void fillDRHistogramHighLowPT(const bool, const std::vector<double>&) const;

      // ----------member data ---------------------------

      //output
      TFile* out_;

      //input
      std::string outFileName_;
      edm::InputTag genParticleTag_;
      edm::InputTag tauTag_;
      std::vector<edm::InputTag> HPSDiscriminatorTags_;
      std::map<std::string, edm::Handle<reco::PFTauDiscriminator> > HPSDiscriminators_;
      unsigned int momPDGID_;

      //histograms
      TH1F* status3TauMultiplicity_;
      TH1F* ZTauTauDecayMultiplicity_;
      TH1F* HPSTauMultiplicity_;
      TH1F* HPSTausPassingEtaCutMultiplicity_;
      TH1F* HPSIsoTausPassingEtaCutMultiplicity_;
      TH1F* dRToHPSTau0_;
      TH1F* dRToHPSTau1_;
      TH1F* dRToHPSIsoTau0_;
      TH1F* dRToHPSIsoTau1_;
      TH2F* HPSTauMultiplicityVsDRToHPSTau0_;
      TH2F* HPSIsoTauMultiplicityVsDRToHPSIsoTau0_;
      TH1F* dRLeadingGenTauToHPSTau0_;
      TH1F* dRTrailingGenTauToHPSTau0_;
      TH1F* dRLeadingGenTauToHPSIsoTau0_;
      TH1F* dRTrailingGenTauToHPSIsoTau0_;
      TH1F* PFTauPT_;
      TH1F* HPSTauPT_;
      TH1F* HPSIsoTauPT_;
      TH1F* HPSIsoTauPTPassingEtaCut_;
      TH1F* genTauPT_;
      TH1F* dRToHPSTau0GenPTGe30_;
      TH1F* dRToHPSIsoTau0GenPTGe30_;

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
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  HPSDiscriminatorTags_(iConfig.getParameter<std::vector<edm::InputTag> >("HPSDiscriminatorTags")),
  momPDGID_(iConfig.getParameter<unsigned int>("momPDGID"))
{
   //now do what ever initialization is needed

}


TauAnalyzer::~TauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   reset();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get GEN particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get HPS PF taus
  edm::Handle<reco::PFTauCollection> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);

  //get HPS discriminators
  for (std::vector<edm::InputTag>::const_iterator 
	 iHPSDiscriminatorTag = HPSDiscriminatorTags_.begin(); 
       iHPSDiscriminatorTag != HPSDiscriminatorTags_.end(); ++iHPSDiscriminatorTag) {
    edm::Handle<reco::PFTauDiscriminator> pHPSDiscriminator;
    iEvent.getByLabel(*iHPSDiscriminatorTag, pHPSDiscriminator);
    HPSDiscriminators_[iHPSDiscriminatorTag->label()] = pHPSDiscriminator;
  }

  //HPS tau multiplicity
  unsigned int nHPSTaus = 0;
  unsigned int nHPSTausPassingEtaCut = 0;
  unsigned int nHPSIsoTausPassingEtaCut = 0;
  std::vector<reco::PFTau*> HPSIsoTaus;
  std::vector<reco::PFTau*> HPSTaus;
  for (reco::PFTauCollection::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) {
    reco::PFTauRef tauRef(pTaus, iTau - pTaus->begin());
    std::map<std::string, edm::Handle<reco::PFTauDiscriminator> >::const_iterator 
      decayModeFinding = HPSDiscriminators_.find("hpsPFTauDiscriminationByDecayModeFinding");
    if ((decayModeFinding != HPSDiscriminators_.end()) && 
	((*(decayModeFinding->second))[tauRef] == 1.0)) {
      ++nHPSTaus;
      std::map<std::string, edm::Handle<reco::PFTauDiscriminator> >::const_iterator 
	looseIsoMVA = HPSDiscriminators_.find("hpsPFTauDiscriminationByLooseIsolationMVA");
      if (fabs(iTau->eta()) < 2.3) {
	++nHPSTausPassingEtaCut;
	if ((looseIsoMVA != HPSDiscriminators_.end()) && 
	    ((*(looseIsoMVA->second))[tauRef] == 1.0)) {
	  ++nHPSIsoTausPassingEtaCut;

	  //fill pT distribution of all isolated HPS taus passing eta cut
	  HPSIsoTauPTPassingEtaCut_->Fill(iTau->pt());
	}
      }

      //save isolated HPS tau
      if ((looseIsoMVA != HPSDiscriminators_.end()) && 
	  ((*(looseIsoMVA->second))[tauRef] == 1.0)) {
	HPSIsoTaus.push_back(const_cast<reco::PFTau*>(&*iTau));

	//fill pT distribution of all isolated HPS taus
	HPSIsoTauPT_->Fill(iTau->pt());
      }

      //save non-isolated HPS tau
      HPSTaus.push_back(const_cast<reco::PFTau*>(&*iTau));

      //fill pT distribution of all HPS taus
      HPSTauPT_->Fill(iTau->pt());
    }

    //fill pT distribution of all PF taus
    PFTauPT_->Fill(iTau->pt());
  }

  //fill multiplicity plots
  HPSTauMultiplicity_->Fill(nHPSTaus);
  HPSTausPassingEtaCutMultiplicity_->Fill(nHPSTausPassingEtaCut);
  HPSIsoTausPassingEtaCutMultiplicity_->Fill(nHPSIsoTausPassingEtaCut);

  //get multiplicity of status 3 taus and Z-->tautau decays
  unsigned int nStatus3TausPerEvt = 0;
  unsigned int nZTauTauDecaysPerEvt = 0;
  reco::Candidate* pMother = NULL;
  std::vector<reco::Candidate*> ZDecayTaus;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {
    if ((fabs(iGenParticle->pdgId()) == 15) && (iGenParticle->status() == 3)) {
      ++nStatus3TausPerEvt;
      size_t numMoms = iGenParticle->numberOfMothers();
      if ((numMoms == 1) && (iGenParticle->mother(0)->pdgId() == (int)momPDGID_)) {

	const reco::Candidate* pZDecayTau = dynamic_cast<const reco::Candidate*>(&*iGenParticle);
	ZDecayTaus.push_back(const_cast<reco::Candidate*>(pZDecayTau));
	if ((pMother == NULL) || 
	    (&*(iGenParticle->mother(0)) != &*pMother)) ++nZTauTauDecaysPerEvt;

	//sort non-isolated HPS taus by proximity to GEN tau
	std::vector<double> dR2MinHPSTau = Common::sortByProximity(HPSTaus, *iGenParticle);
	fillDRHistogramNearFar(0, false, dR2MinHPSTau);
	fillDRHistogramNearFar(1, false, dR2MinHPSTau);

	//sort isolated HPS taus by proximity to GEN tau
	std::vector<double> dR2MinHPSIsoTau = Common::sortByProximity(HPSIsoTaus, *iGenParticle);
	fillDRHistogramNearFar(0, true, dR2MinHPSIsoTau);
	fillDRHistogramNearFar(1, true, dR2MinHPSIsoTau);

	//min dR for gen tau pT > 30 GeV
	if (iGenParticle->pt() > 30.0) {
	  fillDRHistogramHighLowPT(false, dR2MinHPSTau);
	  fillDRHistogramHighLowPT(true, dR2MinHPSIsoTau);
	}

	//fill pT distribution of all GEN taus from desired resonance decay (Z or a)
	genTauPT_->Fill(iGenParticle->pt());
      }
      if (numMoms > 0) pMother = const_cast<reco::Candidate*>(iGenParticle->mother(0));
    }
  }
  status3TauMultiplicity_->Fill(nStatus3TausPerEvt);
  ZTauTauDecayMultiplicity_->Fill(nZTauTauDecaysPerEvt);

  //dR MC matching plots split by leading and trailing GEN tau
  Common::sortByPT(ZDecayTaus);
  if (ZDecayTaus.size() >= 2) {
    std::vector<double> dR2MinLeadingGenTauHPSTau = 
      Common::sortByProximity(HPSTaus, *(ZDecayTaus[0]));
    std::vector<double> dR2MinTrailingGenTauHPSTau = 
      Common::sortByProximity(HPSTaus, *(ZDecayTaus[1]));
    std::vector<double> dR2MinLeadingGenTauHPSIsoTau = 
      Common::sortByProximity(HPSIsoTaus, *(ZDecayTaus[0]));
    std::vector<double> dR2MinTrailingGenTauHPSIsoTau = 
      Common::sortByProximity(HPSIsoTaus, *(ZDecayTaus[1]));
    fillDRHistogramLeadTrail(false, 1, dR2MinLeadingGenTauHPSTau);
    fillDRHistogramLeadTrail(false, 2, dR2MinTrailingGenTauHPSTau);
    fillDRHistogramLeadTrail(true, 1, dR2MinLeadingGenTauHPSIsoTau);
    fillDRHistogramLeadTrail(true, 2, dR2MinTrailingGenTauHPSIsoTau);
  }

  //keep events with >=1 tau above 20 GeV vs. keep only taus above 20 GeV

  //look at tau hadronic decays only, forget about prompt e/mu fakes for the moment
  //apply loose isolation as done by Higgs group
  //how does HPS with/without isolation do for signal MC?

  //check distribution of min dR between gen jet and nearest reco tau
  //check distribution of min dR between gen jet and 2nd nearest reco tau

}


// ------------ method called once each job just before starting event loop  ------------
void 
TauAnalyzer::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book histograms
  status3TauMultiplicity_ = new TH1F("status3TauMultiplicity", "", 7, -0.5, 6.5);
  ZTauTauDecayMultiplicity_ = new TH1F("ZTauTauDecayMultiplicity_", "", 3, -0.5, 2.5);
  HPSTauMultiplicity_ = new TH1F("HPSTauMultiplicity", "", 30, -0.5, 29.5);
  HPSTausPassingEtaCutMultiplicity_ = 
    new TH1F("HPSTausPassingEtaCutMultiplicity", "", 30, -0.5, 29.5);
  HPSIsoTausPassingEtaCutMultiplicity_ = 
    new TH1F("HPSIsoTausPassingEtaCutMultiplicity", "", 30, -0.5, 29.5);
  dRToHPSTau0_ = new TH1F("dRToHPSTau0", "", 500, 0.0, 10.0);
  dRToHPSTau1_ = new TH1F("dRToHPSTau1", "", 500, 0.0, 10.0);
  dRToHPSIsoTau0_ = new TH1F("dRToHPSIsoTau0", "", 500, 0.0, 10.0);
  dRToHPSIsoTau1_ = new TH1F("dRToHPSIsoTau1", "", 500, 0.0, 10.0);
  HPSTauMultiplicityVsDRToHPSTau0_ = 
    new TH2F("HPSTauMultiplicityVsDRToHPSTau0", "", 500, 0.0, 10.0, 30, -0.5, 29.5);
  HPSIsoTauMultiplicityVsDRToHPSIsoTau0_ = 
    new TH2F("HPSIsoTauMultiplicityVsDRToHPSIsoTau0", "", 500, 0.0, 10.0, 30, -0.5, 29.5);
  dRLeadingGenTauToHPSTau0_ = new TH1F("dRLeadingGenTauToHPSTau0", "", 500, 0.0, 10.0);
  dRTrailingGenTauToHPSTau0_ = new TH1F("dRTrailingGenTauToHPSTau0", "", 500, 0.0, 10.0);
  dRLeadingGenTauToHPSIsoTau0_ = new TH1F("dRLeadingGenTauToHPSIsoTau0", "", 500, 0.0, 10.0);
  dRTrailingGenTauToHPSIsoTau0_ = new TH1F("dRTrailingGenTauToHPSIsoTau0", "", 500, 0.0, 10.0);
  PFTauPT_ = new TH1F("PFTauPT", "", 50, 0.0, 100.0);
  HPSTauPT_ = new TH1F("HPSTauPT", "", 50, 0.0, 100.0);
  HPSIsoTauPT_ = new TH1F("HPSIsoTauPT", "", 50, 0.0, 100.0);
  HPSIsoTauPTPassingEtaCut_ = new TH1F("HPSIsoTauPTPassingEtaCut", "", 50, 0.0, 100.0);
  genTauPT_ = new TH1F("genTauPT", "", 50, 0.0, 100.0);
  dRToHPSTau0GenPTGe30_ = new TH1F("dRToHPSTau0GenPTGe30", "", 500, 0.0, 10.0);
  dRToHPSIsoTau0GenPTGe30_ = new TH1F("dRToHPSIsoTau0GenPTGe30", "", 500, 0.0, 10.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauAnalyzer::endJob() 
{
  //write output file
  out_->cd();
  status3TauMultiplicity_->Write();
//   ZTauTauDecayMultiplicity_->Write();
  HPSTauMultiplicity_->Write();
  HPSTausPassingEtaCutMultiplicity_->Write();
  HPSIsoTausPassingEtaCutMultiplicity_->Write();
  dRToHPSTau0_->Write();
  dRToHPSTau1_->Write();
  dRToHPSIsoTau0_->Write();
  dRToHPSIsoTau1_->Write();
  HPSTauMultiplicityVsDRToHPSTau0_->Write();
  HPSIsoTauMultiplicityVsDRToHPSIsoTau0_->Write();
//   dRLeadingGenTauToHPSTau0_->Write();
//   dRTrailingGenTauToHPSTau0_->Write();
//   dRLeadingGenTauToHPSIsoTau0_->Write();
//   dRTrailingGenTauToHPSIsoTau0_->Write();
  PFTauPT_->Write();
  HPSTauPT_->Write();
  HPSIsoTauPT_->Write();
  HPSIsoTauPTPassingEtaCut_->Write();
  genTauPT_->Write();
  dRToHPSTau0GenPTGe30_->Write();
  dRToHPSIsoTau0GenPTGe30_->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
TauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TauAnalyzer::reset()
{
  if (out_ != NULL) delete out_;
  out_ = NULL;
  outFileName_.clear();
  HPSDiscriminatorTags_.clear();
  HPSDiscriminators_.clear();
  momPDGID_ = 0;
  if (status3TauMultiplicity_ != NULL) delete status3TauMultiplicity_;
  status3TauMultiplicity_ = NULL;
  if (ZTauTauDecayMultiplicity_ != NULL) delete ZTauTauDecayMultiplicity_;
  ZTauTauDecayMultiplicity_ = NULL;
  if (HPSTauMultiplicity_ != NULL) delete HPSTauMultiplicity_;
  HPSTauMultiplicity_ = NULL;
  if (HPSTausPassingEtaCutMultiplicity_ != NULL) delete HPSTausPassingEtaCutMultiplicity_;
  HPSTausPassingEtaCutMultiplicity_ = NULL;
  if (HPSIsoTausPassingEtaCutMultiplicity_ != NULL) delete HPSIsoTausPassingEtaCutMultiplicity_;
  HPSIsoTausPassingEtaCutMultiplicity_ = NULL;
  if (dRToHPSTau0_ != NULL) delete dRToHPSTau0_;
  dRToHPSTau0_ = NULL;
  if (dRToHPSTau1_ != NULL) delete dRToHPSTau1_;
  dRToHPSTau1_ = NULL;
  if (dRToHPSIsoTau0_ != NULL) delete dRToHPSIsoTau0_;
  dRToHPSIsoTau0_ = NULL;
  if (dRToHPSIsoTau1_ != NULL) delete dRToHPSIsoTau1_;
  dRToHPSIsoTau1_ = NULL;
  if (HPSTauMultiplicityVsDRToHPSTau0_ != NULL) delete HPSTauMultiplicityVsDRToHPSTau0_;
  HPSTauMultiplicityVsDRToHPSTau0_ = NULL;
  if (HPSIsoTauMultiplicityVsDRToHPSIsoTau0_ != NULL) {
    delete HPSIsoTauMultiplicityVsDRToHPSIsoTau0_;
  }
  HPSIsoTauMultiplicityVsDRToHPSIsoTau0_ = NULL;
  if (dRLeadingGenTauToHPSTau0_ != NULL) delete dRLeadingGenTauToHPSTau0_;
  dRLeadingGenTauToHPSTau0_ = NULL;
  if (dRTrailingGenTauToHPSTau0_ != NULL) delete dRTrailingGenTauToHPSTau0_;
  dRTrailingGenTauToHPSTau0_ = NULL;
  if (dRLeadingGenTauToHPSIsoTau0_ != NULL) delete dRLeadingGenTauToHPSIsoTau0_;
  dRLeadingGenTauToHPSIsoTau0_ = NULL;
  if (dRTrailingGenTauToHPSIsoTau0_ != NULL) delete dRTrailingGenTauToHPSIsoTau0_;
  dRTrailingGenTauToHPSIsoTau0_ = NULL;
  if (PFTauPT_ != NULL) delete PFTauPT_;
  PFTauPT_ = NULL;
  if (HPSTauPT_ != NULL) delete HPSTauPT_;
  HPSTauPT_ = NULL;
  if (HPSIsoTauPT_ != NULL) delete HPSIsoTauPT_;
  HPSIsoTauPT_ = NULL;
  if (HPSIsoTauPTPassingEtaCut_ != NULL) delete HPSIsoTauPTPassingEtaCut_;
  HPSIsoTauPTPassingEtaCut_ = NULL;
  if (genTauPT_ != NULL) delete genTauPT_;
  genTauPT_ = NULL;
  if (dRToHPSTau0GenPTGe30_ != NULL) delete dRToHPSTau0GenPTGe30_;
  dRToHPSTau0GenPTGe30_ = NULL;
  if (dRToHPSIsoTau0GenPTGe30_ != NULL) delete dRToHPSIsoTau0GenPTGe30_;
  dRToHPSIsoTau0GenPTGe30_ = NULL;
}

double TauAnalyzer::dR(const std::vector<double>& dR2, const unsigned int pos) const
{
  double dR = -1.0;
  if (dR2.size() == 0) dR = 100.0;
  else if (pos >= dR2.size()) dR = sqrt(dR2[dR2.size() - 1]);
  else dR = sqrt(dR2[pos]);
  return dR;
}

void TauAnalyzer::fillDRHistogramNearFar(const unsigned int pos, const bool iso, 
					 const std::vector<double>& dR2) const
{
  double deltaR = dR(dR2, pos);
  if (!iso) {
    if (pos == 0) {
      dRToHPSTau0_->Fill(deltaR);
      HPSTauMultiplicityVsDRToHPSTau0_->Fill(deltaR, dR2.size());
    }
    else if (pos == 1) dRToHPSTau1_->Fill(deltaR);
    else {
      std::stringstream err;
      err << "Unsupported position " << pos << " in void TauAnalyzer::fillDRHistogram(const ";
      err << "unsigned int pos, const unsigned int iso, const std::vector<double> dR2) const.\n";
      edm::LogError("TauAnalyzer") << err.str();
    }
  }
  else {
    if (pos == 0) {
      dRToHPSIsoTau0_->Fill(deltaR);
      HPSIsoTauMultiplicityVsDRToHPSIsoTau0_->Fill(deltaR, dR2.size());
    }
    else if (pos == 1) dRToHPSIsoTau1_->Fill(deltaR);
    else {
      std::stringstream err;
      err << "Unsupported position " << pos << " in void TauAnalyzer::fillDRHistogram(const ";
      err << "unsigned int pos, const unsigned int iso, const std::vector<double> dR2) const.\n";
      edm::LogError("TauAnalyzer") << err.str();
    }
  }
}

void TauAnalyzer::fillDRHistogramLeadTrail(const bool iso, const unsigned int lead, 
					   const std::vector<double>& dR2) const
{
  double deltaR = dR(dR2);
  if (!iso) {
    if (lead == 1) dRLeadingGenTauToHPSTau0_->Fill(deltaR);
    else if (lead == 2) dRTrailingGenTauToHPSTau0_->Fill(deltaR);
  }
  else {
    if (lead == 1) dRLeadingGenTauToHPSIsoTau0_->Fill(deltaR);
    else if (lead == 2) dRTrailingGenTauToHPSIsoTau0_->Fill(deltaR);
  }
}

void TauAnalyzer::fillDRHistogramHighLowPT(const bool iso, const std::vector<double>& dR2) const
{
  double deltaR = dR(dR2);
  if (!iso) dRToHPSTau0GenPTGe30_->Fill(deltaR);
  else dRToHPSIsoTau0GenPTGe30_->Fill(deltaR);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
