// -*- C++ -*-
//
// Package:    LeakageFilter
// Class:      LeakageFilter
// 
/**\class LeakageFilter LeakageFilter.cc LEDSoakTools/LeakageFilter/src/LeakageFilter.cc

 Description: filter events based on BX

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Rachel Yohay
//         Created:  Tue Nov 24 13:43:53 CET 2009
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
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

//
// class declaration
//

class LeakageFilter : public edm::EDFilter {
   public:
      explicit LeakageFilter(const edm::ParameterSet&);
      ~LeakageFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
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
LeakageFilter::LeakageFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


LeakageFilter::~LeakageFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
LeakageFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //BX
   bool BXDecision;
   if ((iEvent.bunchCrossing() < 10) || 
       ((iEvent.bunchCrossing() < 2010) && (iEvent.bunchCrossing() > 2000)) || 
       ((iEvent.bunchCrossing() < 3564) && (iEvent.bunchCrossing() > 3519))) BXDecision = true;
   else BXDecision = false;

   //L1 bit
   //from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DQMOffline/Ecal/src/EEClusterTaskExtras.cc?revision=1.6&view=markup
   /*ESHandle<L1GtTriggerMenu> menuRcd;
   iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);
   const L1GtTriggerMenu* menu = menuRcd.product();
   Handle<L1GlobalTriggerReadoutRecord> gtRecord;
   iEvent.getByLabel(InputTag("gtDigis"), gtRecord);
   const DecisionWord dWord = gtRecord->decisionWord();
   bool l1SingleEG2 = menu->gtAlgorithmResult("L1_SingleEG2", dWord);
   bool l1SingleEG5 = menu->gtAlgorithmResult("L1_SingleEG5", dWord);
   bool l1SingleEG8 = menu->gtAlgorithmResult("L1_SingleEG8", dWord);
   bool l1SingleEG10 = menu->gtAlgorithmResult("L1_SingleEG10", dWord);
   bool l1SingleEG12 = menu->gtAlgorithmResult("L1_SingleEG12", dWord);
   bool l1SingleEG15 = menu->gtAlgorithmResult("L1_SingleEG15", dWord);
   bool l1SingleEG20 = menu->gtAlgorithmResult("L1_SingleEG20", dWord);
   bool l1SingleEG25 = menu->gtAlgorithmResult("L1_SingleEG25", dWord);
   bool l1DoubleNoIsoEGBTBtight = menu->gtAlgorithmResult("L1_DoubleNoIsoEG_BTB_tight", dWord);
   bool l1DoubleNoIsoEGBTBloose = menu->gtAlgorithmResult("L1_DoubleNoIsoEG_BTB_loose ", dWord);
   bool l1DoubleNoIsoEGTopBottom = menu->gtAlgorithmResult("L1_DoubleNoIsoEGTopBottom", dWord);
   bool l1DoubleNoIsoEGTopBottomCen  = menu->gtAlgorithmResult("L1_DoubleNoIsoEGTopBottomCen", dWord);
   bool l1DoubleNoIsoEGTopBottomCen2  = menu->gtAlgorithmResult("L1_DoubleNoIsoEGTopBottomCen2", dWord);
   bool l1DoubleNoIsoEGTopBottomCenVert  = menu->gtAlgorithmResult("L1_DoubleNoIsoEGTopBottomCenVert", dWord);
   bool L1Decision = l1SingleEG2 || l1SingleEG5 || l1SingleEG8 || l1SingleEG10 || l1SingleEG12 || l1SingleEG15 || l1SingleEG20 || 
     l1SingleEG25 || l1DoubleNoIsoEGBTBtight || l1DoubleNoIsoEGBTBloose || l1DoubleNoIsoEGTopBottom || l1DoubleNoIsoEGTopBottomCen || 
     l1DoubleNoIsoEGTopBottomCen2 || l1DoubleNoIsoEGTopBottomCenVert;

     return (BXDecision && L1Decision);*/

   return BXDecision;
}

// ------------ method called once each job just before starting event loop  ------------
void 
LeakageFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LeakageFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeakageFilter);
