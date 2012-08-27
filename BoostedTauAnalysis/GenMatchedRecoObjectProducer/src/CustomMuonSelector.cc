// -*- C++ -*-
//
// Package:    CustomMuonSelector
// Class:      CustomMuonSelector
// 
/**\class CustomMuonSelector CustomMuonSelector.cc 
   BoostedTauAnalysis/CustomMuonSelector/src/CustomMuonSelector.cc

 Description: create a collection of custom selected muons to put in the event, and stop 
 processing if no muons are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"

//
// class declaration
//

class CustomMuonSelector : public edm::EDFilter {
public:
  explicit CustomMuonSelector(const edm::ParameterSet&);
  ~CustomMuonSelector();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  //input tag for reco muon collection
  edm::InputTag muonTag_;

  //input tag for reco vertex collection
  edm::InputTag vtxTag_;

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
CustomMuonSelector::CustomMuonSelector(const edm::ParameterSet& iConfig) :
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  vtxTag_(iConfig.getParameter<edm::InputTag>("vtxTag"))
{
  produces<reco::MuonCollection>();
}


CustomMuonSelector::~CustomMuonSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool CustomMuonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<reco::MuonCollection> muonColl(new reco::MuonCollection);

  //get muons
  edm::Handle<reco::MuonCollection> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //get vertices
  edm::Handle<reco::VertexCollection> pVertices;
  iEvent.getByLabel(vtxTag_, pVertices);

  //identify the first good vertex (the "primary" (?))
  reco::Vertex* pPV = Common::getPrimaryVertex(pVertices);

  /*fill STL container with muons passing the 2012 tight selection 
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId)*/
  std::vector<reco::Muon*> tightMuons = Common::getTightRecoMuons(pMuons, pPV);

  //fill output collection
  for (std::vector<reco::Muon*>::const_iterator iMuon = tightMuons.begin(); 
       iMuon != tightMuons.end(); ++iMuon) { muonColl->push_back(**iMuon); }
  iEvent.put(muonColl);

  //if no muons passing cuts were found in this event, stop processing
  return (tightMuons.size() > 0);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CustomMuonSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CustomMuonSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
CustomMuonSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
CustomMuonSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
CustomMuonSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
CustomMuonSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CustomMuonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CustomMuonSelector);
