// -*- C++ -*-
//
// Package:    TEFilter
// Class:      TEFilter
// 
/**\class TEFilter TEFilter.cc LEDGapTools/TEFilter/src/TEFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Rachel Yohay
//         Created:  Thu Oct 15 09:54:39 CEST 2009
// $Id$
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"

//for bad user input
#define INVALID_EVENT_TYPE 255

//ECAL partitions
#define EEP 1
#define EBP 2
#define EBM 3
#define EEM 4
#define NO_PARTITION 0

using namespace std;
using namespace edm;
using namespace cms;

//
// class declaration
//

class TEFilter : public EDFilter {
   public:
      explicit TEFilter(const ParameterSet&);
      ~TEFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(Event&, const EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  //input parameters
  InputTag unpacker_;
  string requestedEventTypeString_;
  unsigned int partition1_, partition2_, partition3_, partition4_;

  short requestedEventTypeShort_;
  const bool isInPartition_(const unsigned int);
  const unsigned int getPartition_(const unsigned int);
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
TEFilter::TEFilter(const ParameterSet& iConfig) :
  unpacker_(iConfig.getParameter<InputTag>("unpacker")),
  requestedEventTypeString_(iConfig.getParameter<string>("eventType")),
  partition1_(iConfig.getParameter<unsigned int>("partition1")),
  partition2_(iConfig.getParameter<unsigned int>("partition2")),
  partition3_(iConfig.getParameter<unsigned int>("partition3")),
  partition4_(iConfig.getParameter<unsigned int>("partition4"))
{
   //now do what ever initialization is needed

  //see http://cmslxr.fnal.gov/lxr/source/DataFormats/EcalRawData/interface/EcalDCCHeaderBlock.h?v=CMSSW_3_2_1#020
  if (requestedEventTypeString_ == "laser") requestedEventTypeShort_ = EcalDCCHeaderBlock::LASER_GAP;
  else if (requestedEventTypeString_ == "test pulse") requestedEventTypeShort_ = EcalDCCHeaderBlock::TESTPULSE_GAP;
  else if (requestedEventTypeString_ == "pedestal") requestedEventTypeShort_ = EcalDCCHeaderBlock::PEDESTAL_GAP;
  else if (requestedEventTypeString_ == "LED") requestedEventTypeShort_ = EcalDCCHeaderBlock::LED_GAP;
  else requestedEventTypeShort_ = INVALID_EVENT_TYPE;
}


TEFilter::~TEFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TEFilter::filter(Event& iEvent, const EventSetup& iSetup)
{
  //filter return value
  bool retVal = false;

  //get the event DCC header block
  bool foundDCCHeaders = true;
  Handle<EcalRawDataCollection> pDCCHeaders;
  foundDCCHeaders = iEvent.getByLabel(unpacker_, pDCCHeaders);
  if (pDCCHeaders->size() == 0) {
    foundDCCHeaders = false;
    cerr << "Error getting DCC headers.\n";
  }

  //apply filtering criteria
  if (foundDCCHeaders) {

    //loop over all FED DCC headers
    for (EcalRawDataCollection::const_iterator iHeader = pDCCHeaders->begin(); iHeader != pDCCHeaders->end(); ++iHeader) {

      //get the event type
      int eventType = iHeader->getRunType();

      //get the DCC
      unsigned int DCC = (unsigned int)(iHeader->getDccInTCCCommand());

      //save the event if it passes the filtering criteria
      if ((eventType == requestedEventTypeShort_) && (isInPartition_(DCC))) {
	retVal = true;
      }
    }
  }

  return retVal;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void 
TEFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TEFilter::endJob() {
}

//determine if the given DCC is in the requested partition(s)
const bool TEFilter::isInPartition_(const unsigned int DCC)
{
  unsigned int partition = getPartition_(DCC);
  if (partition == NO_PARTITION) throw DCC;
  if ((partition == partition1_) || (partition == partition2_) || (partition == partition3_) || (partition == partition4_)) return true;
  else return false;
}

//get the partition of the DCC
const unsigned int TEFilter::getPartition_(const unsigned int DCC)
{
  if ((DCC >= 1) && (DCC <= 9)) return EEM;
  else if ((DCC >= 10) && (DCC <= 27)) return EBM;
  else if ((DCC >= 28) && (DCC <= 45)) return EBP;
  else if ((DCC >= 46) && (DCC <= 54)) return EEP;
  else return NO_PARTITION;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TEFilter);
