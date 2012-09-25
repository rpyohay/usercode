// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      EfficiencyAnalyzer
// 
/**\class EfficiencyAnalyzer EfficiencyAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/EfficiencyAnalyzer.cc

   Description: make efficiency plots

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: EfficiencyAnalyzer.cc,v 1.1 2012/09/19 10:57:14 yohay Exp $
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
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

//
// class declaration
//

template<class T>
class EfficiencyAnalyzer : public edm::EDAnalyzer {
public:
  explicit EfficiencyAnalyzer(const edm::ParameterSet&);
  ~EfficiencyAnalyzer();

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

  //fill pT histogram
  void fillPTHistogram(edm::Handle<edm::View<T> >&, TH1F*);

  //fill eta histogram
  void fillEtaHistogram(edm::Handle<edm::View<T> >&, TH1F*);

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //denominator input tag
  edm::InputTag denominatorTag_;

  //numerator input tag
  edm::InputTag numeratorTag_;

  //marker colors for histograms with different pT rank
  std::vector<unsigned int> pTRankColors_;

  //marker colors for histograms with different decay mode
  std::vector<unsigned int> decayModeColors_;

  //marker styles for histograms with different pT rank
  std::vector<unsigned int> pTRankStyles_;

  //marker styles for histograms with different decay mode
  std::vector<unsigned int> decayModeStyles_;

  //legend entries for histograms with different pT rank
  std::vector<std::string> pTRankEntries_;

  //legend entries for histograms with different decay mode
  std::vector<std::string> decayModeEntries_;

  //histogram of denominator pT
  TH1F* denominatorPT_;

  //histogram of numerator pT
  TH1F* numeratorPT_;

  //histogram of denominator eta
  TH1F* denominatorEta_;

  //histogram of numerator eta
  TH1F* numeratorEta_;
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
template<class T>
EfficiencyAnalyzer<T>::EfficiencyAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  denominatorTag_(iConfig.getParameter<edm::InputTag>("denominatorTag")),
  numeratorTag_(iConfig.getParameter<edm::InputTag>("numeratorTag")),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  decayModeColors_(iConfig.getParameter<std::vector<unsigned int> >("decayModeColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  decayModeStyles_(iConfig.getParameter<std::vector<unsigned int> >("decayModeStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries")),
  decayModeEntries_(iConfig.getParameter<std::vector<std::string> >("decayModeEntries"))
{
  //now do what ever initialization is needed
  reset(false);
}

template<class T>
EfficiencyAnalyzer<T>::~EfficiencyAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
template<class T>
void EfficiencyAnalyzer<T>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get denominator collection
  edm::Handle<edm::View<T> > pDenominatorView;
  iEvent.getByLabel(denominatorTag_, pDenominatorView);

  //get numerator collection
  edm::Handle<edm::View<T> > pNumeratorView;
  iEvent.getByLabel(numeratorTag_, pNumeratorView);

  //plot pT distributions
  fillPTHistogram(pDenominatorView, denominatorPT_);
  fillPTHistogram(pNumeratorView, numeratorPT_);

  //plot eta distributions
  fillEtaHistogram(pDenominatorView, denominatorEta_);
  fillEtaHistogram(pNumeratorView, numeratorEta_);
}


// ------------ method called once each job just before starting event loop  ------------
template<class T>
void EfficiencyAnalyzer<T>::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book pT histograms
  const Double_t bins[11] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 100.0};
  denominatorPT_ = new TH1F("denominatorPT", "", 20, 0.0, 100.0/*10, bins*/);
  numeratorPT_ = new TH1F("numeratorPT", "", 20, 0.0, 100.0/*10, bins*/);

  //book eta histograms
  denominatorEta_ = new TH1F("denominatorEta", "", 20, -5.0, 5.0);
  numeratorEta_ = new TH1F("numeratorEta", "", 20, -5.0, 5.0);
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void EfficiencyAnalyzer<T>::endJob() 
{
  //write output file
  out_->cd();
  denominatorPT_->Write();
  numeratorPT_->Write();
  denominatorEta_->Write();
  numeratorEta_->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
template<class T>
void EfficiencyAnalyzer<T>::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
void EfficiencyAnalyzer<T>::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
void EfficiencyAnalyzer<T>::beginLuminosityBlock(edm::LuminosityBlock const&, 
						    edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
void EfficiencyAnalyzer<T>::endLuminosityBlock(edm::LuminosityBlock const&, 
						  edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void EfficiencyAnalyzer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template<class T>
void EfficiencyAnalyzer<T>::fillPTHistogram(edm::Handle<edm::View<T> >& pView, TH1F* hist)
{
  for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
}

template<class T>
void EfficiencyAnalyzer<T>::fillEtaHistogram(edm::Handle<edm::View<T> >& pView, TH1F* hist)
{
  for (unsigned int i = 0; i < pView->size(); ++i) { hist->Fill(pView->refAt(i)->eta()); }
}

template<class T>
void EfficiencyAnalyzer<T>::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
  if ((doDelete) && (denominatorPT_ != NULL)) delete denominatorPT_;
  denominatorPT_ = NULL;
  if ((doDelete) && (numeratorPT_ != NULL)) delete numeratorPT_;
  numeratorPT_ = NULL;
  if ((doDelete) && (denominatorEta_ != NULL)) delete denominatorEta_;
  denominatorEta_ = NULL;
  if ((doDelete) && (numeratorEta_ != NULL)) delete numeratorEta_;
  numeratorEta_ = NULL;
}

//define this as a plug-in
typedef EfficiencyAnalyzer<reco::GenParticle> GenParticleEfficiencyAnalyzer;
typedef EfficiencyAnalyzer<reco::Muon> MuonEfficiencyAnalyzer;
typedef EfficiencyAnalyzer<reco::PFTau> TauEfficiencyAnalyzer;
DEFINE_FWK_MODULE(GenParticleEfficiencyAnalyzer);
DEFINE_FWK_MODULE(MuonEfficiencyAnalyzer);
DEFINE_FWK_MODULE(TauEfficiencyAnalyzer);
