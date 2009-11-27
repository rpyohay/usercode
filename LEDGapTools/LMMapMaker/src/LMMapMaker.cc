// -*- C++ -*-
//
// Package:    LMMapMaker
// Class:      LMMapMaker
// 
/**\class LMMapMaker LMMapMaker.cc LEDGapTools/LMMapMaker/src/LMMapMaker.cc

 Description: create a text file with the mapping between dee, diffusing sphere, crystal (ix, iy, iz), and crystal hashed index

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Rachel Yohay
//         Created:  Fri Nov  6 15:04:15 CET 2009
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

//
// class decleration
//

using namespace std;
using namespace edm;
using namespace cms;

class LMMapMaker : public edm::EDAnalyzer {
   public:
      explicit LMMapMaker(const edm::ParameterSet&);
      ~LMMapMaker();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  //input parameters
  string mapFileName_;
  string dee1FileName_;
  string dee2FileName_;
  string dee3FileName_;
  string dee4FileName_;

  //file streams
  ofstream mapStream_;
  ifstream dee1Stream_;
  ifstream dee2Stream_;
  ifstream dee3Stream_;
  ifstream dee4Stream_;

  //arrays of information indexed by hashedIndex
  unsigned int dees[EEDetId::kSizeForDenseIndexing];
  unsigned int diffusingSpheres[EEDetId::kSizeForDenseIndexing];
  unsigned int ix[EEDetId::kSizeForDenseIndexing];
  unsigned int iy[EEDetId::kSizeForDenseIndexing];
  int iz[EEDetId::kSizeForDenseIndexing];
};

//
// constants, enums and typedefs
//
const int EEM = -1;
const int EEP = 1;

//
// static data member definitions
//

//
// constructors and destructor
//
LMMapMaker::LMMapMaker(const edm::ParameterSet& iConfig) :
  mapFileName_(iConfig.getParameter<string>("mapFileName")),
  dee1FileName_(iConfig.getParameter<string>("dee1FileName")),
  dee2FileName_(iConfig.getParameter<string>("dee2FileName")),
  dee3FileName_(iConfig.getParameter<string>("dee3FileName")),
  dee4FileName_(iConfig.getParameter<string>("dee4FileName"))
{
   //now do what ever initialization is needed

}


LMMapMaker::~LMMapMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
LMMapMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



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
LMMapMaker::beginJob()
{
  //open the files
  mapStream_.open(mapFileName_.c_str());
  if (!mapStream_.good()) {
    cerr << "Error opening file " << mapFileName_ << ".\n";
    return;
  }
  dee1Stream_.open(dee1FileName_.c_str());
  if (!dee1Stream_.good()) {
    cerr << "Error opening file " << dee1FileName_ << ".\n";
    mapStream_.close();
    return;
  }
  dee2Stream_.open(dee2FileName_.c_str());
  if (!dee2Stream_.good()) {
    cerr << "Error opening file " << dee2FileName_ << ".\n";
    mapStream_.close();
    dee1Stream_.close();
    return;
  }
  dee3Stream_.open(dee3FileName_.c_str());
  if (!dee3Stream_.good()) {
    cerr << "Error opening file " << dee3FileName_ << ".\n";
    mapStream_.close();
    dee1Stream_.close();
    dee2Stream_.close();
    return;
  }
  dee4Stream_.open(dee4FileName_.c_str());
  if (!dee4Stream_.good()) {
    cerr << "Error opening file " << dee4FileName_ << ".\n";
    mapStream_.close();
    dee1Stream_.close();
    dee2Stream_.close();
    dee3Stream_.close();
    return;
  }

  //loop over the dee 1 file
  while (!dee1Stream_.eof()) {

    //get the data
    unsigned int diffusingSphere, x, y;
    dee1Stream_ >> diffusingSphere >> x >> y;
    dee1Stream_.seekg(1, ios::cur);
    const int hashedIndex = (const int)EEDetId(x, y, EEP, EEDetId::XYMODE).hashedIndex();

    //fill the arrays
    dees[hashedIndex] = 1;
    diffusingSpheres[hashedIndex] = diffusingSphere;
    ix[hashedIndex] = x;
    iy[hashedIndex] = y;
    iz[hashedIndex] = EEP;
  }

  //loop over the dee 2 file
  while (!dee2Stream_.eof()) {

    //get the data
    unsigned int diffusingSphere, x, y;
    dee2Stream_ >> diffusingSphere >> x >> y;
    dee2Stream_.seekg(1, ios::cur);
    const int hashedIndex = (const int)EEDetId(x, y, EEP, EEDetId::XYMODE).hashedIndex();

    //fill the arrays
    dees[hashedIndex] = 2;
    diffusingSpheres[hashedIndex] = diffusingSphere;
    ix[hashedIndex] = x;
    iy[hashedIndex] = y;
    iz[hashedIndex] = EEP;
  }

  //loop over the dee 3 file
  while (!dee3Stream_.eof()) {

    //get the data
    unsigned int diffusingSphere, x, y;
    dee3Stream_ >> diffusingSphere >> x >> y;
    dee3Stream_.seekg(1, ios::cur);
    const int hashedIndex = (const int)EEDetId(x, y, EEM, EEDetId::XYMODE).hashedIndex();

    //fill the arrays
    dees[hashedIndex] = 3;
    diffusingSpheres[hashedIndex] = diffusingSphere;
    ix[hashedIndex] = x;
    iy[hashedIndex] = y;
    iz[hashedIndex] = EEM;
  }

  //loop over the dee 4 file
  while (!dee4Stream_.eof()) {

    //get the data
    unsigned int diffusingSphere, x, y;
    dee4Stream_ >> diffusingSphere >> x >> y;
    dee4Stream_.seekg(1, ios::cur);
    const int hashedIndex = (const int)EEDetId(x, y, EEM, EEDetId::XYMODE).hashedIndex();

    //fill the arrays
    dees[hashedIndex] = 4;
    diffusingSpheres[hashedIndex] = diffusingSphere;
    ix[hashedIndex] = x;
    iy[hashedIndex] = y;
    iz[hashedIndex] = EEM;
  }

  //write the values in the array to the map file
  for (unsigned int iCrystal = 0; iCrystal < EEDetId::kSizeForDenseIndexing; ++iCrystal) {
    mapStream_ << dees[iCrystal] << " " << diffusingSpheres[iCrystal] << " ";
    mapStream_ << ix[iCrystal] << " " << iy[iCrystal] << " " << iz[iCrystal] << " ";
    mapStream_ << iCrystal << endl;
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LMMapMaker::endJob() {
  //close all files
  mapStream_.close();
  dee1Stream_.close();
  dee2Stream_.close();
  dee3Stream_.close();
  dee4Stream_.close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(LMMapMaker);
