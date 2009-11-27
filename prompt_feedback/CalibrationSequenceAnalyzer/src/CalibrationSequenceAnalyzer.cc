// -*- C++ -*-
//
// Package:    CalibrationSequenceAnalyzer
// Class:      CalibrationSequenceAnalyzer
// 
/**\class CalibrationSequenceAnalyzer CalibrationSequenceAnalyzer.cc prompt_feedback/CalibrationSequenceAnalyzer/src/CalibrationSequenceAnalyzer.cc

 Description: analyzes LED events for quality, prompt feedback style

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Rachel Yohay
//         Created:  Wed Sep 30 14:37:25 CEST 2009
// $Id$
//
//

// system include files
#include <memory>
#include <fstream>
//#include <sys/time>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
//#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"

#define NUM_EEM_DCCS 9

using namespace std;
using namespace edm;
using namespace cms;

//
// class decleration
//
class CalibrationSequenceAnalyzer : public edm::EDAnalyzer {
   public:
      explicit CalibrationSequenceAnalyzer(const edm::ParameterSet&);
      ~CalibrationSequenceAnalyzer();

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  void getCrystalsInDCC_(int);
  bool isDigiInActiveFED_(const EEDetId) const;
  void setFitParameters_(TF1*, Double_t, Float_t, Float_t);

      // ----------member data ---------------------------
  //input parameters
      InputTag EEDigiCollection_;
  InputTag unpacker_;
  string outFile_;
  //string outtxtFile_;

  //mapping
  int activeFED_;
  vector<EEDetId> crystalsInActiveFED_;
  EcalElectronicsMapping mapper_;

  //output files
  TFile* out_;
  //ofstream outtxt_;

  //fit parameters
  double alpha_;
  double tRise_;

  //event data tree
  TTree* tree_;
  int evt_, ix_, iy_, iz_, DCC_, hashedIndex_;
  float t_, A_;
  TimeValue_t timestamp_;

  //keep track of which crystals you've drawn pulses for
  int drawnPulses[EEDetId::kSizeForDenseIndexing];
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
CalibrationSequenceAnalyzer::CalibrationSequenceAnalyzer(const edm::ParameterSet& iConfig) :
  EEDigiCollection_(iConfig.getParameter<InputTag>("EEDigiCollection")),
  unpacker_(iConfig.getParameter<InputTag>("unpacker")),
  outFile_(iConfig.getParameter<string>("outFile")),
  //outtxtFile_(iConfig.getParameter<string>("outtxtFile")),
  alpha_(iConfig.getParameter<double>("alpha")),
  tRise_(iConfig.getParameter<double>("tRise"))
{
   //now do what ever initialization is needed
}

CalibrationSequenceAnalyzer::~CalibrationSequenceAnalyzer()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete out_;
}

/*alpha/beta fit (cf. CMS NOTE 2004/025):
A(t) = (((t - t0)/beta)^alpha)*e^(-alpha((t - t0 - beta)/beta)) + pedestal
x[0] = t (time)
par[0] = pedestal
par[1] = amplitude
par[2] = tMax
par[3] = alpha
par[4] = tRise*/
Double_t alphaBeta_(Double_t* x, Double_t* par) {
  if (x[0] < (par[2] - par[4])) return par[0];
  return par[1]*pow((x[0] - par[2] + par[4])/par[4], par[3])*exp(-par[3]*((x[0] - par[2])/par[4])) + par[0];
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
CalibrationSequenceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get the ECAL mapping stuff
  /*ESHandle<EcalElectronicsMapping> ECALMapping;
  iSetup.get<EcalMappingRcd>().get(ECALMapping);
  const EcalElectronicsMapping* pMap = ECALMapping.product();*/

   //get the digis
   bool digisFound = true;
   Handle<EEDigiCollection> pDigis;
   digisFound = iEvent.getByLabel(EEDigiCollection_, pDigis);
   if (pDigis->size() == 0) {
     digisFound = false;
     cerr << "Error getting the EE digis.\n";
   }

   //process the event
   if (digisFound) {

     //get the event DCC header block
     bool foundDCCHeaders = true;
     Handle<EcalRawDataCollection> pDCCHeaders;
     foundDCCHeaders = iEvent.getByLabel(unpacker_, pDCCHeaders);
     if (pDCCHeaders->size() == 0) {
       foundDCCHeaders = false;
       cerr << "Error getting DCC headers.\n";
     }

     //search for the active FED
     if (foundDCCHeaders) {
       EcalRawDataCollection::const_iterator firstHeader = pDCCHeaders->begin();
       activeFED_ = firstHeader->getDccInTCCCommand();
       /*vector<DetId> detIDs = pMap->dccConstituents(activeFED_);
       for (vector<DetId>::const_iterator iID = detIDs.begin(); iID != detIDs.end(); ++iID) {
	 crystalsInActiveFED_.push_back(EEDetId(*iID));
       }
       for (vector<EEDetId>::const_iterator iCrystal = crystalsInActiveFED_.begin(); iCrystal != crystalsInActiveFED_.end(); ++iCrystal) {
	 cout << "crystal " << iCrystal->rawId() << ", ix = " << iCrystal->ix() << ", iy = " << iCrystal->iy() << endl;
	 }*/

       //skip if not EE-
       //if ((activeFED_ >= 1) && (activeFED_ <= 9)) {

	 //initialize event-level tree variables
	 evt_ = iEvent.id().event();
	 timestamp_ = iEvent.time().value();
	 DCC_ = activeFED_;

	 //loop over the digis
	 for (EEDigiCollection::const_iterator iDigi = pDigis->begin(); iDigi != pDigis->end(); ++iDigi) {

	   //create an EcalDataFrame object
	   EcalDataFrame digi(*iDigi);

	   //is this digi in the active FED?
	   const DetId ID = digi.id();
	   EEDetId EEID(ID);
	   //if (isDigiInActiveFED_(EEID)) { //broken

	     //fill arrays of time samples and corresponding ADC counts and their errors
	     //determine the maximum amplitude sample
	     Float_t t[EcalDataFrame::MAXSAMPLES];
	     Float_t A[EcalDataFrame::MAXSAMPLES];
	     Float_t tErr[EcalDataFrame::MAXSAMPLES];
	     Float_t AErr[EcalDataFrame::MAXSAMPLES];
	     Float_t maxA = 0;
	     Float_t tOfMaxA = 0;
	     for (int iBX = 0; iBX < EcalDataFrame::MAXSAMPLES; ++iBX) {
	       t[iBX] = iBX;
	       int rawADCVal = digi.sample(iBX).adc();
	       switch (digi.sample(iBX).gainId()) {
	       case 1:
		 A[iBX] = rawADCVal;
		 break;
	       case 2:
		 A[iBX] = rawADCVal*6;
		 break;
	       case 3:
		 A[iBX] = rawADCVal*12;
		 break;
	       }
	       tErr[iBX] = 0;
	       AErr[iBX] = 2; //estimate of the pedestal noise from many runs
	       if (A[iBX] > maxA) {
		 maxA = A[iBX];
		 tOfMaxA = iBX;
	       }
	     }

	     //make a graph of the pulse shape
	     TGraphErrors* pulse = new TGraphErrors(EcalDataFrame::MAXSAMPLES, t, A, tErr, AErr);

	     //make an alpha/beta fit function
	     TF1* fit = new TF1("fit", alphaBeta_, 0, EcalDataFrame::MAXSAMPLES - 1, 5);
	     Double_t ped = /*(*/A[0]/* + A[1] + A[2])/3.0*/;
	     setFitParameters_(fit, ped, maxA, tOfMaxA);
	     pulse->Fit(fit, "QREB");
	     pulse->SetMarkerStyle(1);

	     //initialize digi-level tree variables
	     ix_ = EEID.ix();
	     iy_ = EEID.iy();
	     iz_ = EEID.zside();
	     hashedIndex_ = EEID.hashedIndex();
	     t_ = fit->GetParameter(2);
	     A_ = fit->GetParameter(1);

	     //draw and save the pulse with the fit overlaid
	     /*if ((t_ < 1.0) && (drawnPulses[hashedIndex_] == 0)) {
	       char canvasName[7];
	       sprintf(canvasName, "%i_%i_%i", ix_, iy_, evt_);
	       TCanvas* canvas = new TCanvas(canvasName, canvasName, 500, 500);
	       canvas->SetWindowSize(1000 - canvas->GetWw(), 1000 - canvas->GetWh());
	       canvas->SetFillStyle(4000);
	       canvas->SetFillColor(0); //white
	       out_->cd();
	       canvas->cd();
	       pulse->Draw("A*");
	       canvas->Update();
	       canvas->Write();
	       drawnPulses[hashedIndex_] = 1;
	       }*/

	     //fill the tree
	     tree_->Fill();

	     //write the data to the output text file
	     /*outtxt_ << ix_ << " " << iy_ << " " << iz_ << " " << hashedIndex_ << " " << DCC_ << " " << evt_ << " " << timestamp_ << endl;
	     for (int i = 0; i < EcalDataFrame::MAXSAMPLES; ++ i) {
	       outtxt_ << A[i];
	       if (i < (EcalDataFrame::MAXSAMPLES - 1)) outtxt_ << " ";
	       else outtxt_ << "\n\n";
	       }*/

	     //delete pointers
	     delete fit;
	     delete pulse;

	     //}//if (isDigiInActiveFED(digi.id())) //broken
	 }//for (EEDigiCollection::const_iterator iDigi = pDigis->begin(); iDigi != pDigis->end(); ++iDigi)
	 //}//if ((activeFED_ >= 1) && (activeFED_ <= 9))
     }//if (foundDCCHeaders)
   }//if (digisFound)
}


// ------------ method called once each job just before starting event loop  ------------
void 
CalibrationSequenceAnalyzer::beginJob()
{
  //initialize private data members
  activeFED_ = 0;
  for (unsigned int i = 0; i < EEDetId::kSizeForDenseIndexing; ++i) { drawnPulses[i] = 0; }

  //open output ROOT file
  out_ = new TFile(outFile_.c_str(), "RECREATE");
  if (!out_->IsOpen()) cerr << "Error opening file " << outFile_ << ".\n";

  //initialize the event data tree
  tree_ = new TTree("tree_", "tree_");
  tree_->Branch("evt_", &evt_, "evt_/I");
  tree_->Branch("timestamp_", &timestamp_, "timestamp_/l");
  tree_->Branch("ix_", &ix_, "ix_/I");
  tree_->Branch("iy_", &iy_, "iy_/I");
  tree_->Branch("iz_", &iz_, "iz_/I");
  tree_->Branch("hashedIndex_", &hashedIndex_, "hashedIndex_/I");
  tree_->Branch("t_", &t_, "t_/F");
  tree_->Branch("A_", &A_, "A_/F");
  tree_->Branch("DCC_", &DCC_, "DCC_/I");

  //open output text file
  /*outtxt_.open(outtxtFile_.c_str());
  if (!outtxt_.good()) cerr << "Error opening file " << outtxtFile_ << ".\n";
  outtxt_ << "ix iy iz hashed_index DCC event timestamp\nt0 t1 t2 t3 t4 t5 t6 t7 t8 t9\n\n";*/
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CalibrationSequenceAnalyzer::endJob() {
  //save and close the output files
  if (out_->IsOpen()) {
    out_->cd();
    tree_->Write();
    out_->Write();
    out_->Close();
  }
  //outtxt_.close();
}

//test whether given crystal is in the active FED
bool CalibrationSequenceAnalyzer::isDigiInActiveFED_(const EEDetId digiID) const
{
  bool inFED = false;
  vector<EEDetId>::const_iterator iID = crystalsInActiveFED_.begin();
  while ((iID != crystalsInActiveFED_.end()) && (!inFED)) {
    if (*iID == digiID) inFED = true;
    else ++iID;
  }
  return inFED;
}

//set alpha/beta fit parameters
void CalibrationSequenceAnalyzer::setFitParameters_(TF1* fit, Double_t ped, Float_t maxA, Float_t tOfMaxA)
{
  fit->SetParLimits(0, 0, 1000);
  fit->SetParLimits(1, 0, 4000);
  fit->SetParLimits(2, 0, EcalDataFrame::MAXSAMPLES - 1);
  fit->SetParLimits(3, 0.01, 10);
  fit->SetParLimits(4, 0, 10);
  fit->SetParameter(0, ped);
  fit->SetParameter(1, maxA - ped);
  fit->SetParameter(2, tOfMaxA);
  fit->SetParameter(3, 1.0);
  fit->SetParameter(4, 3.0);
  fit->FixParameter(0, ped);
  fit->FixParameter(3, alpha_);
  fit->FixParameter(4, tRise_);
  fit->SetParNames("pedestal", "amplitude", "tMax", "alpha", "tRise");
  fit->SetLineColor(6);
  fit->SetLineWidth(6);
  fit->SetRange(tOfMaxA - 1.1, tOfMaxA + 4.1);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CalibrationSequenceAnalyzer);

//money plots:
//average amplitude per crystal
//A/A0 vs. t - t0 per crystal
//timestamp distribution for all events in 1 FED
