// -*- C++ -*-
//
// Package:    ObjectTriggerAnalyzer
// Class:      ObjectTriggerAnalyzer
// 
/**\class ObjectTriggerAnalyzer ObjectTriggerAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/ObjectTriggerAnalyzer.cc

   Description: analyze objects firing triggers

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: ObjectTriggerAnalyzer.cc,v 1.3 2012/08/27 14:55:29 yohay Exp $
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
#include "TFile.h"
#include "TH1F.h"

//
// class declaration
//

template<class T>
class ObjectTriggerAnalyzer : public edm::EDAnalyzer {
public:
  explicit ObjectTriggerAnalyzer(const edm::ParameterSet&);
  ~ObjectTriggerAnalyzer();

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

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //denominator input tag
  edm::InputTag denominatorTag_;

  //numerator input tag
  edm::InputTag numeratorTag_;

  //histogram of denominator pT
  TH1F* denominatorPT_;

  //histogram of numerator pT
  TH1F* numeratorPT_;

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
ObjectTriggerAnalyzer<T>::ObjectTriggerAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  denominatorTag_(iConfig.getParameter<edm::InputTag>("denominatorTag")),
  numeratorTag_(iConfig.getParameter<edm::InputTag>("numeratorTag"))
{
  //now do what ever initialization is needed
  reset(false);
}

template<class T>
ObjectTriggerAnalyzer<T>::~ObjectTriggerAnalyzer()
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
void ObjectTriggerAnalyzer<T>::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

}


// ------------ method called once each job just before starting event loop  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::beginJob()
{
  //open output file
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book histograms
  denominatorPT_ = new TH1F("denominatorPT", "", 20, 0.0, 100.0);
  numeratorPT_ = new TH1F("numeratorPT", "", 20, 0.0, 100.0);
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::endJob() 
{
  //write output file
  out_->cd();
  denominatorPT_->Write();
  numeratorPT_->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::beginLuminosityBlock(edm::LuminosityBlock const&, 
						    edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
void ObjectTriggerAnalyzer<T>::endLuminosityBlock(edm::LuminosityBlock const&, 
						  edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void ObjectTriggerAnalyzer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

template<class T>
void ObjectTriggerAnalyzer<T>::fillPTHistogram(edm::Handle<edm::View<T> >& pView, TH1F* hist)
{
  for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
}

template<class T>
void ObjectTriggerAnalyzer<T>::reset(const bool doDelete)
{
  if ((doDelete) && (out_ != NULL)) delete out_;
  out_ = NULL;
  if ((doDelete) && (denominatorPT_ != NULL)) delete denominatorPT_;
  denominatorPT_ = NULL;
  if ((doDelete) && (numeratorPT_ != NULL)) delete numeratorPT_;
  numeratorPT_ = NULL;
}

//define this as a plug-in
typedef ObjectTriggerAnalyzer<reco::GenParticle> GenParticleTriggerAnalyzer;
DEFINE_FWK_MODULE(GenParticleTriggerAnalyzer);
