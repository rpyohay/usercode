// -*- C++ -*-
//
// Package:    PATVerifier
// Class:      PATVerifier
// 
/**\class PATVerifier PATVerifier.cc PATTools/PATVerifier/src/PATVerifier.cc

 Description: verify PAT recipe for photons

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Rachel YOHAY
//         Created:  Thu Jan 14 14:17:29 CET 2010
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
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

//ROOT include files
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
//
// class decleration
//

using namespace cms;
using namespace edm;
using namespace std;

class PATVerifier : public edm::EDAnalyzer {
   public:
      explicit PATVerifier(const edm::ParameterSet&);
      ~PATVerifier();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
  //output
  string outName_;
  TFile* out_;

  //input
  InputTag photonSrc_;

      // data members for histograms to be filled

      // PhotonID Histograms
      TH1F* h_isoEcalRecHit_;
      TH1F* h_isoHcalRecHit_;
      TH1F* h_trk_pt_solid_;
      TH1F* h_trk_pt_hollow_;
      TH1F* h_ntrk_solid_;
      TH1F* h_ntrk_hollow_;
      TH1F* h_ebgap_;
      TH1F* h_eeGap_;
      TH1F* h_ebeeGap_;
      TH1F* h_r9_;

      // Photon Histograms
      TH1F* h_photonEt_;
      TH1F* h_photonEta_;
      TH1F* h_photonPhi_;
      TH1F* h_hadoverem_;

      // Photon's SuperCluster Histograms
      TH1F* h_photonScEt_;
      TH1F* h_photonScEta_;
      TH1F* h_photonScPhi_;
      TH1F* h_photonScEtaWidth_;

      // Composite or Other Histograms
      TH1F* h_photonInAnyGap_;
      TH1F* h_nPassingPho_;
      TH1F* h_nPho_;

  //cuts
      double minPhotonEt_;       // minimum photon Et
      double minPhotonAbsEta_;   // min and
      double maxPhotonAbsEta_;   // max abs(eta)
      double minPhotonR9_;       // minimum R9 = E(3x3)/E(SuperCluster)
      double maxPhotonHoverE_;   // maximum HCAL / ECAL
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
PATVerifier::PATVerifier(const edm::ParameterSet& iConfig) :
  outName_(iConfig.getParameter<string>("outName")),
  photonSrc_(iConfig.getParameter<InputTag>("photonSrc")),
  minPhotonEt_(iConfig.getParameter<double>("minPhotonEt")),
  minPhotonAbsEta_(iConfig.getParameter<double>("minPhotonAbsEta")),
  maxPhotonAbsEta_(iConfig.getParameter<double>("maxPhotonAbsEta")),
  minPhotonR9_(iConfig.getParameter<double>("minPhotonR9")),
  maxPhotonHoverE_(iConfig.getParameter<double>("maxPhotonHoverE"))

{
   //now do what ever initialization is needed

}


PATVerifier::~PATVerifier()
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
PATVerifier::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get the cleanLayer1Photons
   edm::Handle<edm::View<pat::Photon> > photons;
   bool foundPhotons = iEvent.getByLabel(photonSrc_,photons);
   if ((!foundPhotons) || (photons->size() == 0)) {
     cerr << "No pat::Photon objects found in event " << iEvent.id().event() << ".\n";
     return;
   }

   //loop over the photons
   int photonCounter = 0;
   cerr << "----------\n";
   cerr << "iEvent.id().event() = " << iEvent.id().event() << endl;
   for (edm::View<pat::Photon>::const_iterator iPhoton=photons->begin(); iPhoton!=photons->end(); ++iPhoton){
     float photonEt       = iPhoton->et();
     cerr << "photonEt = " << photonEt << endl;
     float superClusterEt = (iPhoton->superCluster()->energy())/(cosh(iPhoton->superCluster()->position().eta()));

    // Only store photon candidates (SuperClusters) that pass some simple cuts
     cerr << "fabs(iPhoton->eta()) = " << fabs(iPhoton->eta()) << endl;
     cerr << "iPhoton->r9() = " << iPhoton->r9() << endl;
     cerr << "iPhoton->hadronicOverEm() = " << iPhoton->hadronicOverEm() << endl;
    bool passCuts = (              photonEt > minPhotonEt_     ) &&
                    (      fabs(iPhoton->eta()) > minPhotonAbsEta_ ) &&
                    (      fabs(iPhoton->eta()) < maxPhotonAbsEta_ ) &&
                    (          iPhoton->r9() > minPhotonR9_     ) &&
                    ( iPhoton->hadronicOverEm() < maxPhotonHoverE_ ) ;

    if ( passCuts )
    {
      cerr << "Passed\n";
      ///////////////////////////////////////////////////////
      //                fill histograms                    //
      ///////////////////////////////////////////////////////
      // PhotonID Variables
      h_isoEcalRecHit_->Fill(iPhoton->ecalRecHitSumEtConeDR04());
      h_isoHcalRecHit_->Fill(iPhoton->hcalTowerSumEtConeDR04());
      h_trk_pt_solid_ ->Fill(iPhoton->trkSumPtSolidConeDR04());
      h_trk_pt_hollow_->Fill(iPhoton->trkSumPtHollowConeDR04());
      h_ntrk_solid_->   Fill(iPhoton->nTrkSolidConeDR04());
      h_ntrk_hollow_->  Fill(iPhoton->nTrkHollowConeDR04());
      h_ebgap_->        Fill(iPhoton->isEBGap());
      h_eeGap_->        Fill(iPhoton->isEEGap());
      h_ebeeGap_->      Fill(iPhoton->isEBEEGap());
      h_r9_->           Fill(iPhoton->r9());

      // Photon Variables
      h_photonEt_->  Fill(photonEt);
      h_photonEta_-> Fill(iPhoton->eta());
      h_photonPhi_-> Fill(iPhoton->phi());
      h_hadoverem_-> Fill(iPhoton->hadronicOverEm());

      // Photon's SuperCluster Variables
      // eta is with respect to detector (not physics) vertex,
      // thus Et and eta are different from photon.
      h_photonScEt_->      Fill(superClusterEt);
      h_photonScEta_->     Fill(iPhoton->superCluster()->position().eta());
      h_photonScPhi_->     Fill(iPhoton->superCluster()->position().phi());
      h_photonScEtaWidth_->Fill(iPhoton->superCluster()->etaWidth());

      // It passed photon cuts, mark it
      h_nPassingPho_->Fill(1.0);

      // Record whether it was near any module gap.
      // Very convoluted at the moment.
      bool inAnyGap = iPhoton->isEBEEGap() || (iPhoton->isEB()&&iPhoton->isEBGap()) || (iPhoton->isEE()&&iPhoton->isEEGap());
      if (inAnyGap) {
        h_photonInAnyGap_->Fill(1.0);
      } else {
        h_photonInAnyGap_->Fill(0.0);
      }

      photonCounter++;
    }
    else
    {
      cerr << "Failed\n";
      // This didn't pass photon cuts, mark it
      h_nPassingPho_->Fill(0.0);
    }

  } // End Loop over photons
  h_nPho_->Fill(photonCounter);
}


// ------------ method called once each job just before starting event loop  ------------
void 
PATVerifier::beginJob()
{
  // Book Histograms
  // PhotonID Histograms
  h_isoEcalRecHit_ = new TH1F("photonEcalIso",          "Ecal Rec Hit Isolation", 100, 0, 100);
  h_isoHcalRecHit_ = new TH1F("photonHcalIso",          "Hcal Rec Hit Isolation", 100, 0, 100);
  h_trk_pt_solid_  = new TH1F("photonTrackSolidIso",    "Sum of track pT in a cone of #DeltaR" , 100, 0, 100);
  h_trk_pt_hollow_ = new TH1F("photonTrackHollowIso",   "Sum of track pT in a hollow cone" ,     100, 0, 100);
  h_ntrk_solid_    = new TH1F("photonTrackCountSolid",  "Number of tracks in a cone of #DeltaR", 100, 0, 100);
  h_ntrk_hollow_   = new TH1F("photonTrackCountHollow", "Number of tracks in a hollow cone",     100, 0, 100);
  h_ebgap_         = new TH1F("photonInEBgap",          "Ecal Barrel gap flag",  2, -0.5, 1.5);
  h_eeGap_         = new TH1F("photonInEEgap",          "Ecal Endcap gap flag",  2, -0.5, 1.5);
  h_ebeeGap_       = new TH1F("photonInEEgap",          "Ecal Barrel/Endcap gap flag",  2, -0.5, 1.5);
  h_r9_            = new TH1F("photonR9",               "R9 = E(3x3) / E(SuperCluster)", 300, 0, 3);

  // Photon Histograms
  h_photonEt_      = new TH1F("photonEt",     "Photon E_{T}",  500,  0, 1000); //changed from 200 1-GeV bins
  h_photonEta_     = new TH1F("photonEta",    "Photon #eta",   200, -4,   4);
  h_photonPhi_     = new TH1F("photonPhi",    "Photon #phi",   200, -1.*TMath::Pi(), TMath::Pi());
  h_hadoverem_     = new TH1F("photonHoverE", "Hadronic over EM", 200, 0, 1);

  // Photon's SuperCluster Histograms
  h_photonScEt_       = new TH1F("photonScEt",  "Photon SuperCluster E_{T}", 500,  0, 1000); //changed from 200 1-GeV bins
  h_photonScEta_      = new TH1F("photonScEta", "Photon #eta",               200, -4,   4);
  h_photonScPhi_      = new TH1F("photonScPhi", "Photon #phi", 200, -1.*TMath::Pi(), TMath::Pi());
  h_photonScEtaWidth_ = new TH1F("photonScEtaWidth","#eta-width",            100,  0,  .1);

  // Composite or Other Histograms
  h_photonInAnyGap_   = new TH1F("photonInAnyGap",     "Photon in any gap flag",  2, -0.5, 1.5);
  h_nPassingPho_      = new TH1F("photonPassingCount", "Total number photons (0=NotPassing, 1=Passing)", 2, -0.5, 1.5);
  h_nPho_             = new TH1F("photonCount",        "Number of photons passing cuts in event",  10,  0,  10);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATVerifier::endJob() {
  //write histograms to ROOT files
  out_ = new TFile(outName_.c_str(), "RECREATE");
  if (out_->IsOpen()) {
    out_->cd();
      // PhotonID Histograms
  h_isoEcalRecHit_->Write();
  h_isoHcalRecHit_->Write();
  h_trk_pt_solid_-> Write();
  h_trk_pt_hollow_->Write();
  h_ntrk_solid_->   Write();
  h_ntrk_hollow_->  Write();
  h_ebgap_->     Write();
  h_eeGap_->     Write();
  h_ebeeGap_->   Write();
  h_r9_->        Write();

  // Photon Histograms
  h_photonEt_->  Write();
  h_photonEta_-> Write();
  h_photonPhi_-> Write();
  h_hadoverem_-> Write();

  // Photon's SuperCluster Histograms
  h_photonScEt_->      Write();
  h_photonScEta_->     Write();
  h_photonScPhi_->     Write();
  h_photonScEtaWidth_->Write();

  // Composite or Other Histograms
  h_photonInAnyGap_->Write();
  h_nPassingPho_->   Write();
  h_nPho_->          Write();

  // Write the root file (really writes the TTree)
  out_->Write();
  out_->Close();
  }
  else cerr << "Error opening file " << outName_ << ".\n";
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATVerifier);
