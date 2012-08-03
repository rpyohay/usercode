// -*- C++ -*-
//
// Package:    HadronicTauDecayFinder
// Class:      HadronicTauDecayFinder
// 
/**\class HadronicTauDecayFinder HadronicTauDecayFinder.cc BoostedTauAnalysis/TauSkimmer/plugins/HadronicTauDecayFinder.cc

 Description: return true if at least 1 hadronic tau decay was found at GEN level

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 11:29:40 CEST 2012
// $Id$
//
//


// system include files
#include <memory>
#include <algorithm>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//constants
#define EPDGID 11
#define MUPDGID 13
#define TAUPDGID 15
#define ENEUTRINOPDGID 12
#define MUNEUTRINOPDGID 14
#define TAUNEUTRINOPDGID 16

//
// class declaration
//

class HadronicTauDecayFinder : public edm::EDFilter {
   public:
      explicit HadronicTauDecayFinder(const edm::ParameterSet&);
      ~HadronicTauDecayFinder();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      void reset();

      // ----------member data ---------------------------

      //input
      edm::InputTag genParticleTag_;
      unsigned int momPDGID_;

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
HadronicTauDecayFinder::HadronicTauDecayFinder(const edm::ParameterSet& iConfig) :
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  momPDGID_(iConfig.getParameter<unsigned int>("momPDGID"))
{
   //now do what ever initialization is needed

}


HadronicTauDecayFinder::~HadronicTauDecayFinder()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   reset();

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HadronicTauDecayFinder::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get gen particles
   Handle<reco::GenParticleCollection> pGenParticles;
   iEvent.getByLabel(genParticleTag_, pGenParticles);

   //loop over gen particles
   bool foundHadronicDecay = false;
   reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin();
   while ((iGenParticle != pGenParticles->end()) && !foundHadronicDecay) {

     //look for a status 3 tau from boson decay
     if ((fabs(iGenParticle->pdgId()) == TAUPDGID) && (iGenParticle->status() == 3) && 
	 (iGenParticle->numberOfMothers() == 1) && 
	 (iGenParticle->mother(0)->pdgId() == (int)momPDGID_)) {

       //search for a leptonic decay (stop searching when an e-nu or mu-nu pair is found)
       unsigned int leptonPDGID = 0;
       unsigned int neutrinoPDGID = 0;
       const size_t numDaughters = iGenParticle->numberOfDaughters();
       if (numDaughters == 1) {
	 const reco::Candidate* daughter = iGenParticle->daughter(0);
	 reco::Candidate::const_iterator iDaughter = daughter->begin();
	 while ((iDaughter != daughter->end()) && ((leptonPDGID == 0) || (neutrinoPDGID == 0))) {
	   if (leptonPDGID == 0) {
	     if (((neutrinoPDGID == 0) || (neutrinoPDGID == ENEUTRINOPDGID)) && 
		 (fabs(iDaughter->pdgId()) == EPDGID)) leptonPDGID = EPDGID;
	     if (((neutrinoPDGID == 0) || (neutrinoPDGID == MUNEUTRINOPDGID)) && 
		 (fabs(iDaughter->pdgId()) == MUPDGID)) leptonPDGID = MUPDGID;
	   }
	   if (neutrinoPDGID == 0) {
	     if (((leptonPDGID == 0) || (leptonPDGID == EPDGID)) && 
		 (fabs(iDaughter->pdgId()) == ENEUTRINOPDGID)) neutrinoPDGID = ENEUTRINOPDGID;
	     if (((leptonPDGID == 0) || (leptonPDGID == MUPDGID)) && 
		 (fabs(iDaughter->pdgId()) == MUNEUTRINOPDGID)) neutrinoPDGID = MUNEUTRINOPDGID;
	   }
	   ++iDaughter;
	 }
       }
       else {
	 std::stringstream err;
	 err << "Gen particle " << iGenParticle - pGenParticles->begin() << " has ";
	 err << numDaughters << " daughters.\n";
	 throw cms::Exception("HadronicTauDecayFinder") << err.str();
       }

       //leptonic decay not found ==> found a hadronic decay ==> stop searching, this event passes
       if ((leptonPDGID == 0) || (neutrinoPDGID == 0)) foundHadronicDecay = true;
     }

     //advance to the next gen particle
     ++iGenParticle;
   }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return foundHadronicDecay;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HadronicTauDecayFinder::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HadronicTauDecayFinder::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
HadronicTauDecayFinder::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
HadronicTauDecayFinder::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
HadronicTauDecayFinder::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
HadronicTauDecayFinder::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HadronicTauDecayFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void HadronicTauDecayFinder::reset() { momPDGID_ = 0; }
//define this as a plug-in
DEFINE_FWK_MODULE(HadronicTauDecayFinder);
