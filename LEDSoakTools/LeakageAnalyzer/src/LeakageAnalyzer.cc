// -*- C++ -*-
//
// Package:    LeakageAnalyzer
// Class:      LeakageAnalyzer
// 
/**\class LeakageAnalyzer LeakageAnalyzer.cc LEDSoakTools/LeakageAnalyzer/src/LeakageAnalyzer.cc

 Description: look for evidence of leakage of LED pulse tails into the physics triggers

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Rachel YOHAY
//         Created:  Mon Nov 23 11:02:01 CET 2009
// $Id$
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
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

//ROOT include files
#include "TH2F.h"
#include "TFile.h"

//inputs to getSample
#define MIN 0
#define MAX 1

#define PERIOD_CAL_ORBITS 112

using namespace std;
using namespace edm;

//
// class decleration
//

class LeakageAnalyzer : public edm::EDAnalyzer {
   public:
      explicit LeakageAnalyzer(const edm::ParameterSet&);
      ~LeakageAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  //get min/max sample of the digi
  const unsigned int maxSample(EcalDataFrame*);
  const unsigned int minSample(EcalDataFrame*);
  const unsigned int getSample(EcalDataFrame*, const unsigned int);
  const int gainAdjustedADC(EcalMGPASample*);

      // ----------member data ---------------------------
      InputTag EEDigiCollection_;

  //histograms
  TH2F* AVsBX_;
  TH2F* AVsBX2_;
  TH2F* AVsBX3_;
  TH2F* numDigisPerEvtVsBX_;
  TH1F* numDigisPerEvt_;
  TH2F* AVsBXByOrbit_[PERIOD_CAL_ORBITS];

  //output
  string outName_;
  TFile* out_;
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
LeakageAnalyzer::LeakageAnalyzer(const edm::ParameterSet& iConfig) :
  EEDigiCollection_(iConfig.getParameter<InputTag>("EEDigiCollection")),
  outName_(iConfig.getParameter<string>("outName"))
{
   //now do what ever initialization is needed

}


LeakageAnalyzer::~LeakageAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete out_;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
LeakageAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //print the L1 path passed by this event
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
  cout << "Event " << iEvent.id().event();
  if (L1Decision == true) cout << " fired ";
  else cout << " did not fire ";
  cout << "an L1 e/gamma trigger.\n";*/

  //get the EE digis
  bool digisFound = true;
  Handle<EEDigiCollection> pDigis;
  digisFound = iEvent.getByLabel(EEDigiCollection_, pDigis);
  if (pDigis->size() == 0) {
    digisFound = false;
    cerr << "Error getting the EE digis.  Processing of this event will be skipped.\n";
  }
  if (digisFound) {
    
    //sum the EE energy
    int sumEEEnergy = 0;
    int sumEEEnergy2 = 0;
    for (EEDigiCollection::const_iterator iDigi = pDigis->begin(); iDigi != pDigis->end(); ++iDigi) {
      EcalDataFrame* digi = new EcalDataFrame(*iDigi);
      sumEEEnergy+=(maxSample(digi) - minSample(digi));
      sumEEEnergy2+=(gainAdjustedADC(&(digi->sample(5))) - gainAdjustedADC(&(digi->sample(0))));
      AVsBX3_->Fill(iEvent.bunchCrossing(), gainAdjustedADC(&(digi->sample(5))) - gainAdjustedADC(&(digi->sample(0))));
      delete digi;
    }

    //fill the histogram
    AVsBX_->Fill(iEvent.bunchCrossing(), sumEEEnergy);
    AVsBX2_->Fill(iEvent.bunchCrossing(), sumEEEnergy2);
    numDigisPerEvtVsBX_->Fill(iEvent.bunchCrossing(), pDigis->size());
    numDigisPerEvt_->Fill(pDigis->size());
  }

  //get the EB digis (for comparison, because there is absolutely no LED light)
}


// ------------ method called once each job just before starting event loop  ------------
void 
LeakageAnalyzer::beginJob()
{
  //book histogram 1
  AVsBX_ = new TH2F("AVsBX_", "Total EE amplitude (max - min) vs. event BX", 3564, -0.5, 3563.5, 4000, 1500.0, 5500.0);
  AVsBX_->GetXaxis()->SetTitle("Event BX");
  AVsBX_->GetXaxis()->CenterTitle();
  AVsBX_->GetYaxis()->SetTitle("#Sigma_{crystals}(max sample - min sample) (ADC counts)");
  AVsBX_->GetYaxis()->CenterTitle();
  AVsBX_->GetYaxis()->SetTitleOffset(1.25);

  //book histogram 2
  AVsBX2_ = new TH2F("AVsBX2_", "Total EE amplitude (6^{th} sample - 1^{st} sample) vs. event BX", 
		     3564, -0.5, 3563.5, 4000, 1500.0, 5500.0);
  AVsBX2_->GetXaxis()->SetTitle("Event BX");
  AVsBX2_->GetXaxis()->CenterTitle();
  AVsBX2_->GetYaxis()->SetTitle("#Sigma_{crystals}(6^{th} sample - 1^{st} sample) (ADC counts)");
  AVsBX2_->GetYaxis()->CenterTitle();
  AVsBX2_->GetYaxis()->SetTitleOffset(1.2);
  AVsBX2_->GetYaxis()->SetLabelSize(0.03);

  //book histogram 3
  AVsBX3_ = new TH2F("AVsBX3_", "EE amplitude per crystal (6^{th} sample - 1^{st} sample) vs. event BX", 
		     3564, -0.5, 3563.5, 40, -20.0, 20.0);
  AVsBX3_->GetXaxis()->SetTitle("Event BX");
  AVsBX3_->GetXaxis()->CenterTitle();
  AVsBX3_->GetYaxis()->SetTitle("Amplitude/crystal (6^{th} sample - 1^{st} sample) (ADC counts)");
  AVsBX3_->GetYaxis()->CenterTitle();
  AVsBX3_->GetYaxis()->SetTitleOffset(1.1);
  AVsBX3_->GetYaxis()->SetTitleSize(0.03);

  //book histogram 4
  numDigisPerEvtVsBX_ = new TH2F("numDigisPerEvtVsBX_", "Number of EE digis per event vs. event BX", 
		     3564, -0.5, 3563.5, 700, 400.0, 1100.0);
  numDigisPerEvtVsBX_->GetXaxis()->SetTitle("Event BX");
  numDigisPerEvtVsBX_->GetXaxis()->CenterTitle();
  numDigisPerEvtVsBX_->GetYaxis()->SetTitle("Number of EE digis per event");
  numDigisPerEvtVsBX_->GetYaxis()->CenterTitle();
  numDigisPerEvtVsBX_->GetYaxis()->SetTitleOffset(1.2);
  numDigisPerEvtVsBX_->GetYaxis()->SetTitleSize(0.04);

  //book histogram 5
  numDigisPerEvt_ = new TH1F("numDigisPerEvt_", "Number of EE digis per event", 700, 400.0, 1100.0);
  numDigisPerEvt_->GetXaxis()->SetTitle("Number of digis/event");
  numDigisPerEvt_->GetXaxis()->CenterTitle();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LeakageAnalyzer::endJob() {
  //save histogram to output file
  out_ = new TFile(outName_.c_str(), "RECREATE");
  if (out_->IsOpen()) {
    out_->cd();
    AVsBX_->Write();
    AVsBX2_->Write();
    AVsBX3_->Write();
    numDigisPerEvtVsBX_->Write();
    numDigisPerEvt_->Write();
    out_->Close();
  }
  else { cerr << "Error opening file " << outName_ << ".  No output will be saved.\n"; }
}

//get the time of the max sample
const unsigned int LeakageAnalyzer::maxSample(EcalDataFrame* digi) { return getSample(digi, MAX); }

//get the time of the min sample
const unsigned int LeakageAnalyzer::minSample(EcalDataFrame* digi) { return getSample(digi, MIN); }

//get the time of the max or min sample
const unsigned int LeakageAnalyzer::getSample(EcalDataFrame* digi, const unsigned int minOrMax)
{
  unsigned int t = 0;
  int A = gainAdjustedADC(&(digi->sample(0)));
  for (int iBX = 1; iBX < EcalDataFrame::MAXSAMPLES; ++iBX) {
    int ADCVal = gainAdjustedADC(&(digi->sample(iBX)));
    if (((minOrMax == MAX) && (ADCVal > A)) || ((minOrMax == MIN) && (ADCVal < A))) {
      A = ADCVal;
      t = iBX;
    }
  }
  return t;
}

//get the gain-adjusted ADC value of the sample
const int LeakageAnalyzer::gainAdjustedADC(EcalMGPASample* sample)
{
  int A = sample->adc();
  const unsigned int gainID = sample->gainId();
  if (gainID > 1) A = A*(gainID - 1)*6;
  return A;
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeakageAnalyzer);
