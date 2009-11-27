// -*- C++ -*-
//
// Package:    OrbitNumberAnalyzer
// Class:      OrbitNumberAnalyzer
// 
/**\class OrbitNumberAnalyzer OrbitNumberAnalyzer.cc LEDSoakTools/OrbitNumberAnalyzer/src/OrbitNumberAnalyzer.cc

 Description: for a given run, determine the orbit numbers for a specified orbit type

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Rachel Yohay
//         Created:  Wed Nov 25 16:00:20 CET 2009
// $Id$
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
//
// class decleration
//

class OrbitNumberAnalyzer : public edm::EDAnalyzer {
   public:
      explicit OrbitNumberAnalyzer(const edm::ParameterSet&);
      ~OrbitNumberAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  unsigned int orbitType_;
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
OrbitNumberAnalyzer::OrbitNumberAnalyzer(const edm::ParameterSet& iConfig) :
  orbitType_(iConfig.getParameter<unsigned int>("orbitType")

{
   //now do what ever initialization is needed

}


OrbitNumberAnalyzer::~OrbitNumberAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
OrbitNumberAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   

}


// ------------ method called once each job just before starting event loop  ------------
void 
OrbitNumberAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
OrbitNumberAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(OrbitNumberAnalyzer);
