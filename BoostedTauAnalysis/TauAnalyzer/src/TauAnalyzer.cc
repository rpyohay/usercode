// -*- C++ -*-
//
// Package:    TauAnalyzer
// Class:      TauAnalyzer
// 
/**\class TauAnalyzer TauAnalyzer.cc 
   BoostedTauAnalysis/TauAnalyzer/src/TauAnalyzer.cc

   Description: analyze tau variables that can discriminate signal to background

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Wed Jul 18 16:40:51 CEST 2012
// $Id: TauAnalyzer.cc,v 1.4 2012/10/04 15:10:13 yohay Exp $
//
//


// system include files
#include <memory>
#include <string>
// #include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

//
// class declaration
//

class TauAnalyzer : public edm::EDAnalyzer {
public:
  explicit TauAnalyzer(const edm::ParameterSet&);
  ~TauAnalyzer();

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

  //fill pT histogram with object of arbitrary type
  template<typename U>
  void fillPTHistogramArbitrary(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->pt());
  }

  //fill ET histogram
  template<typename U>
  void fillETHistogram(edm::Handle<edm::View<U> >& pView, TH1F* hist)
  {
    for (unsigned int i = 0; i < pView->size(); ++i) hist->Fill(pView->refAt(i)->et());
  }

  //make visible gen tau pT canvas for 1 decay mode, multiple pT ranks
  void makePTRankCanvas(TCanvas&, TLegend&, const std::string&, std::vector<TH1F*>&);

  //format and draw multiple pT histograms on one canvas
  void drawMultiplePTHistograms(TCanvas&, std::vector<TH1F*>&, const std::vector<unsigned int>&, 
				const std::vector<unsigned int>&, TLegend&, 
				const std::vector<std::string>&, const std::string&);

  // ----------member data ---------------------------

  //pointer to output file object
  TFile* out_;

  //name of output file
  std::string outFileName_;

  //tau tag
  edm::InputTag tauTag_;

  //MET tag
  edm::InputTag METTag_;

  //muon tag
  edm::InputTag muonTag_;

  //jet-muon map tag
  edm::InputTag jetMuonMapTag_;

  //marker colors for histograms with different pT rank
  std::vector<unsigned int> pTRankColors_;

  //marker styles for histograms with different pT rank
  std::vector<unsigned int> pTRankStyles_;

  //legend entries for histograms with different pT rank
  std::vector<std::string> pTRankEntries_;

  //histogram of MET
  TH1F* MET_;

  //histogram of mu+had mass
  TH1F* muHadMass_;

  //histogram of mu charge + had charge
  TH1F* muHadCharge_;

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
TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  tauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  METTag_(iConfig.getParameter<edm::InputTag>("METTag")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonTag")),
  jetMuonMapTag_(iConfig.getParameter<edm::InputTag>("jetMuonMapTag")),
  pTRankColors_(iConfig.getParameter<std::vector<unsigned int> >("pTRankColors")),
  pTRankStyles_(iConfig.getParameter<std::vector<unsigned int> >("pTRankStyles")),
  pTRankEntries_(iConfig.getParameter<std::vector<std::string> >("pTRankEntries"))
{
  //now do what ever initialization is needed
  reset(false);
}

TauAnalyzer::~TauAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get taus
  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByLabel(tauTag_, pTaus);
  
  //get MET
  edm::Handle<edm::View<reco::PFMET> > pMET;
  iEvent.getByLabel(METTag_, pMET);

  //get muons
  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByLabel(muonTag_, pMuons);

  //get jet-muon map
  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByLabel(jetMuonMapTag_, pMuonJetMap);

//   //debug
//   for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
//        ++iTau) {
//     const reco::PFJetRef& tauJetRef = (*iTau)->jetRef();
//     std::cerr << "Tau jet ref key: " << tauJetRef.key() << std::endl;
//     const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];
//     for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
// 	 iMuon != removedMuons.end(); ++iMuon) {
//       std::cerr << "Muon ref key: " << iMuon->key() << std::endl;
//       std::cerr << "Muon pT: " << (*iMuon)->pt() << " GeV\n";
//       std::cerr << "Muon eta: " << (*iMuon)->eta() << std::endl;
//       std::cerr << "Muon phi: " << (*iMuon)->phi() << std::endl;
//     }
//   }

  //plot MET distribution
  fillETHistogram(pMET, MET_);  
}


// ------------ method called once each job just before starting event loop  ------------
void TauAnalyzer::beginJob()
{
  //open output files
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  //book jet histograms
  MET_ = new TH1F("MET", "", 20, 0.0, 100.0);
}

// ------------ method called once each job just after ending the event loop  ------------
void TauAnalyzer::endJob() 
{
  //make the jet canvases
  out_->cd();

  //write output files
  out_->cd();
  MET_->Write();
  out_->Write();
  out_->Close();
}

// ------------ method called when starting to processes a run  ------------
void TauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void TauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void TauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, 
				       edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void TauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, 
				     edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void TauAnalyzer::makePTRankCanvas(TCanvas& canvas, TLegend& legend, 
				   const std::string& header, 
				   std::vector<TH1F*>& hists)
{
  drawMultiplePTHistograms(canvas, hists, pTRankColors_, pTRankStyles_, legend, 
			   pTRankEntries_, header);
}

void TauAnalyzer::drawMultiplePTHistograms(TCanvas& canvas, 
					   std::vector<TH1F*>& hists, 
					   const std::vector<unsigned int>& colors, 
					   const std::vector<unsigned int>& styles, 
					   TLegend& legend, 
					   const std::vector<std::string>& entries, 
					   const std::string& header)
{
  Common::setLegendOptions(legend, header.c_str());
  canvas.cd();
  Common::setCanvasOptions(canvas, 1, 0, 0);
  TH1F* pHistWithMaxMaxBin = NULL;
  Double_t maxBinContent = 0.0;
  for (std::vector<TH1F*>::iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    const unsigned int i = iHist - hists.begin();
    Common::setHistogramOptions(*iHist, colors[i], 0.7, styles[i], 1.0, 
				"p_{T} (GeV)", "", 0.04);
    legend.AddEntry(*iHist, entries[i].c_str(), "l");
    Int_t histMaxBin = (*iHist)->GetMaximumBin();
    Double_t histMaxBinContent = (*iHist)->GetBinContent(histMaxBin);
    if (histMaxBinContent >= maxBinContent) {
      maxBinContent = histMaxBinContent;
      pHistWithMaxMaxBin = *iHist;
    }
  }
  pHistWithMaxMaxBin->Draw();
  for (std::vector<TH1F*>::iterator iHist = hists.begin(); 
       iHist != hists.end(); ++iHist) {
    (*iHist)->Draw("SAME");
  }
  legend.Draw();
  canvas.Write();
}

void TauAnalyzer::reset(const bool doDelete)
{
  if (doDelete && (out_ != NULL)) delete out_;
  out_ = NULL;
  if (doDelete && (MET_ != NULL)) delete MET_;
  MET_ = NULL;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
