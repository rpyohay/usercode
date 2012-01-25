// -*- C++ -*-
//
// Package:    GMSBTools
// Class:      EventAnalyzer.cc
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: EventAnalyzer.cc,v 1.11 2011/06/07 20:30:39 dwjang Exp $
//

#define EventAnalyzer_cxx

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TObjectTable.h>

#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>

#include "EventAnalyzer.h"
#include "SusyEventPrinter.h"

#include "../jec/JetMETObjects/interface/JetCorrectorParameters.h"
#include "../jec/JetMETObjects/interface/FactorizedJetCorrector.h"


template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}


void EventAnalyzer::InitializePerEvent() {

}


bool EventAnalyzer::isSameObject(TLorentzVector& p1, TLorentzVector& p2) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < 0.5) return true;
  return false;
}


float EventAnalyzer::d0correction(TVector3& beamSpot, susy::Track& track) const {

  float d0 = track.d0() - beamSpot.X()*std::sin(track.phi()) + beamSpot.Y()*std::cos(track.phi());
  return d0;
}


bool EventAnalyzer::PassTrigger(TString path) {
  bool pass = false;
  for(susy::TriggerMap::iterator it = event->hltMap.begin(); it != event->hltMap.end(); it++) {
    if(it->first.Contains(path) && (int(it->second.second)) ) {
      pass = true;
      break;
    }
  }
  return pass;
}


bool EventAnalyzer::PassTriggers() {
  bool pass = false;
  for(std::vector<TString>::iterator it = hltNames.begin(); it != hltNames.end(); it++) {
    if(PassTrigger(*it)) {
      pass = true;
      break;
    }
  }
  return pass;
}

/*this function only works if you believe the prescales from LumiSummary that are saved in the 
  ntuple*/
std::string EventAnalyzer::PassLowestUnprescaledTrigger() const
{
  unsigned int minET = 9999;
  std::string minHLT;
  if (printLevel > 0) std::cout << "---------------------\n";
  for (std::vector<TString>::const_iterator iDesiredHLT = hltNames.begin(); 
       iDesiredHLT != hltNames.end(); ++iDesiredHLT) {
    susy::TriggerMap::const_iterator iHLT = event->hltMap.find(*iDesiredHLT);
    if ((iHLT != event->hltMap.end()) && (iHLT->second.first == 1) && (iHLT->second.second > 0)) {
      std::string HLTName((const char*)iHLT->first);
      if (printLevel > 0) std::cout << "Fired unprescaled trigger " << HLTName << std::endl;
//       size_t secondPhotonPos = HLTName.rfind("Photon");
      size_t firstPhotonPos = HLTName.find("Photon");
//       size_t fourthUnderscorePos = HLTName.rfind("_", HLTName.rfind("_") - 1);
      size_t secondUnderscorePos = HLTName.find("_", HLTName.find("_") + 1);
//       if ((secondPhotonPos != std::string::npos) && (fourthUnderscorePos != std::string::npos)) {
      if ((firstPhotonPos != std::string::npos) && (secondUnderscorePos != std::string::npos)) {
// 	std::string ET(HLTName.substr(secondPhotonPos + 6, 
// 				      fourthUnderscorePos - secondPhotonPos - 6));
	std::string ET(HLTName.substr(firstPhotonPos + 6, 
				      secondUnderscorePos - firstPhotonPos - 6));
	int ETInt = atoi(ET.c_str());
	if (ETInt < (int)minET) {
	  minET = ETInt;
	  minHLT = HLTName;
	}
      }
      else throw HLTName;
    }
  }
  if (printLevel > 0) {
    std::cout << "Lowest ET threshold unprescaled trigger fired by event " << event->eventNumber;
    std::cout << ", run " << event->runNumber << ", lumi section " << event->luminosityBlockNumber;
    std::cout << ": " << ((minHLT == "") ? "none" : minHLT) << std::endl;
  }
  return minHLT;
}

//this function uses foreknowledge of the proper prescales
std::string EventAnalyzer::PassUnprescaledTrigger() const
{
  std::string passedHLT;
  if (printLevel > 0) std::cout << "---------------------\n";
  for (std::vector<TString>::const_iterator iDesiredHLT = hltNames.begin(); 
       iDesiredHLT != hltNames.end(); ++iDesiredHLT) {
    susy::TriggerMap::const_iterator iHLT = event->hltMap.begin();
    bool foundTrigger = false;
    while ((iHLT != event->hltMap.end()) && !foundTrigger) {
      if (iHLT->first.Contains(*iDesiredHLT)) foundTrigger = true;
      else ++iHLT;
    }
    if ((iHLT != event->hltMap.end()) && (iHLT->second.second > 0)) {
      std::string HLTName((const char*)iHLT->first);
      size_t secondPhotonPos = HLTName.rfind("Photon");
      size_t firstPhotonPos = HLTName.find("Photon");
      size_t fourthUnderscorePos = HLTName.rfind("_", HLTName.rfind("_") - 1);
      size_t secondUnderscorePos = HLTName.find("_", HLTName.find("_") + 1);
      if ((firstPhotonPos != std::string::npos) && (secondUnderscorePos != std::string::npos) && 
	  (secondPhotonPos != std::string::npos) && (fourthUnderscorePos != std::string::npos)) {
	std::string trailingET(HLTName.substr(secondPhotonPos + 6, 
					      fourthUnderscorePos - secondPhotonPos - 6));
	std::string leadingET(HLTName.substr(firstPhotonPos + 6, 
					     secondUnderscorePos - firstPhotonPos - 6));
	int trailingETInt = atoi(trailingET.c_str());
	int leadingETInt = atoi(leadingET.c_str());
	/*hack until ntuple prescales are fixed*/
	if (((event->runNumber <= 161176) && (leadingETInt == 32) && (trailingETInt == 26)) || 
	    ((event->runNumber > 161176) && (event->runNumber <= 163869) && 
	     (((leadingETInt == 32) && (trailingETInt == 26)) || 
	      ((leadingETInt == 36) && (trailingETInt == 22)))) || 
	    ((event->runNumber > 163869) && (event->runNumber <= 166967) && 
	     (((leadingETInt == 36) && (trailingETInt == 22)) || 
	      ((leadingETInt == 40) && (trailingETInt == 28)))) || 
	    ((event->runNumber > 166967) && (leadingETInt == 40) && (trailingETInt == 28))) {
	  passedHLT = HLTName;
	}
	/*hack until ntuple prescales are fixed*/
      }
      else throw HLTName;
    }
  }
  if (printLevel > 0) {
    std::cout << "Unprescaled trigger fired by event " << event->eventNumber;
    std::cout << ", run " << event->runNumber << ", lumi section " << event->luminosityBlockNumber;
    std::cout << ": " << ((passedHLT == "") ? "none" : passedHLT) << std::endl;
  }
  return passedHLT;
}

void EventAnalyzer::countTriggers(const std::string& fileLabel)
{
  //counters
  unsigned int nPassingJSON = 0;
  unsigned int nPassingOR = 0;
  unsigned int nPassingLowestUnprescaled = 0;
  unsigned int nPassingEta = 0;
  unsigned int nPassingET = 0;
  unsigned int nPassingEtaOR = 0;
  unsigned int nPassingEtaLowestUnprescaled = 0;
  unsigned int nPassingETOR = 0;
  unsigned int nPassingETLowestUnprescaled = 0;

  //histograms
  //                                   32/26   36/22
  //                         v 32/26 v 36/22 v 40/28 v 40/28 v
  Float_t runRanges[5] = {160404, 161177, 163870, 166968, 167913};
  Float_t triggers[4] = {0, 1, 2, 3};
  TH2F* triggerPassedVsRunRange = new TH2F("triggerPassedVsRunRange", "", 
					   4, runRanges, 3, triggers);
  triggerPassedVsRunRange->GetYaxis()->SetBinLabel(1, "32/26");
  triggerPassedVsRunRange->GetYaxis()->SetBinLabel(2, "36/22");
  triggerPassedVsRunRange->GetYaxis()->SetBinLabel(3, "40/28");
  TH1F* leadingPhotonET = new TH1F("leadingPhotonET", "", 60, 0.0, 300.0);
  TH1F* trailingPhotonET = new TH1F("trailingPhotonET", "", 60, 0.0, 300.0);
  TH1F* leadingPhotonETPassingETORNotLowestUnprescaled = 
    new TH1F("leadingPhotonETPassingETORNotLowestUnprescaled", "", 60, 0.0, 300.0);
  TH1F* trailingPhotonETPassingETORNotLowestUnprescaled = 
    new TH1F("trailingPhotonETPassingETORNotLowestUnprescaled", "", 60, 0.0, 300.0);

  //events that pass by one trigger counting method and not the other
  std::vector<Int_t> passEtaORNotLowestUnprescaledRun;
  std::vector<ULong_t> passEtaORNotLowestUnprescaledEvt;
  std::vector<Int_t> passEtaORNotLowestUnprescaledLumiSection;
  std::vector<Int_t> passEtaLowestUnprescaledNotORRun;
  std::vector<ULong_t> passEtaLowestUnprescaledNotOREvt;
  std::vector<Int_t> passEtaLowestUnprescaledNotORLumiSection;
  std::vector<Int_t> passETORNotLowestUnprescaledRun;
  std::vector<ULong_t> passETORNotLowestUnprescaledEvt;
  std::vector<Int_t> passETORNotLowestUnprescaledLumiSection;
  std::vector<susy::TriggerMap> passETORNotLowestUnprescaledHLTMap;
  std::vector<Int_t> passETLowestUnprescaledNotORRun;
  std::vector<ULong_t> passETLowestUnprescaledNotOREvt;
  std::vector<Int_t> passETLowestUnprescaledNotORLumiSection;

  //remove this memory-intensive step
  //to check duplicated events
  // std::map<int, std::set<int> > allEvents;

  //open files
  ofstream out((fileLabel + ".txt").c_str());
  TFile outROOT((fileLabel + ".root").c_str(), "RECREATE");
  if (!out.is_open() || !outROOT.IsOpen()) {
    std::cout << "Error opening file " << fileLabel << ".txt or " << fileLabel << ".root.\n";
    delete triggerPassedVsRunRange;
    return;
  }
  outROOT.cd();

  //loop over events
  if (fChain == 0) return;
  Long64_t nEvts = (processNEvents != -1) ? processNEvents : fChain->GetEntriesFast();
  for (Long64_t iEvt = 0; iEvt < nEvts; ++iEvt) {
    if (printLevel > 0) std::cout << "Get the tree contents." << std::endl;
    Long64_t localEntry = LoadTree(iEvt);
    if (localEntry < 0) break;

    /*branches you don't need:
      - everything but hltMap and photons
     */
    fChain->GetEntry(iEvt); //fix this to only load the branches you need

    //print debug information
    if (printLevel > 0 || 
	(printInterval > 0 && (/*iEvt >= printInterval && */iEvt % printInterval == 0))) {
      std::cout << (int(iEvt) + 1) << " events processed with run=" << event->runNumber;
      std::cout << ", event=" << event->eventNumber << ", lumi section = ";
      std::cout << event->luminosityBlockNumber << std::endl;
    }

    //initialize global variables
    if(printLevel > 0) {
      std::cout << "Initialize any global variables to be reset per event." << std::endl;
    }
    InitializePerEvent();

    //apply JSON selection
    if (printLevel > 0) std::cout << "Apply good run list." << std::endl;
    if (!isInJson(event->runNumber, event->luminosityBlockNumber)) {
      continue;
    }
    ++nPassingJSON;

//     //uncomment this to print all ntuple variables
//     Print(*event);

    //check for duplicate events
    //get rid of this memory-intensive step
    // if (printLevel > 0) std::cout << "Check duplicated events for data only." << std::endl;
    // bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
    // if (event->isRealData && duplicateEvent) {
    //   continue;
    // }

    //count number of triggers
    bool passOR = useTrigger ? PassTriggers() : true;
    if (passOR) ++nPassingOR;
    std::string HLTPath;
    try { HLTPath = PassUnprescaledTrigger(); }
    catch (std::string badHLTPath) {
      out << "Error: badly formatted trigger name " << badHLTPath << ".  Ending event loop at ";
      out << "event " << iEvt << " (run " << event->runNumber << ", event " << event->eventNumber;
      out << ", lumi section " << event->luminosityBlockNumber << ".\n";
      break;
    }
    if (HLTPath != "") {
      ++nPassingLowestUnprescaled;
      unsigned int bin = 0;
      if ((HLTPath.find("32") != std::string::npos) && 
	  (HLTPath.find("26") != std::string::npos)) bin = 0;
      if ((HLTPath.find("36") != std::string::npos) && 
	  (HLTPath.find("22") != std::string::npos)) bin = 1;
      if ((HLTPath.find("40") != std::string::npos) && 
	  (HLTPath.find("28") != std::string::npos)) bin = 2;
      triggerPassedVsRunRange->Fill(event->runNumber, bin);
    }

    //get the photons
    std::map<TString,susy::PhotonCollection> photonMap = event->photons;
    std::map<TString,susy::PhotonCollection>::const_iterator iPhotonMap = 
      photonMap.find((const char*)tag_);
    unsigned int nPhotonsPassingEta = 0;
    unsigned int nPassingLeadingET = 0;
    unsigned int nPassingTrailingET = 0;
    bool passETAndEta = false;
    if (iPhotonMap != photonMap.end()) {
      susy::PhotonCollection photons = iPhotonMap->second;

      //loop over photons, keeping track of how many pass the ET and eta cuts
      susy::PhotonCollection::const_iterator iPhoton = photons.begin();
      while ((iPhoton != photons.end()) && !passETAndEta) {
	if (fabs(iPhoton->caloPosition.Eta()) < 1.4442) {
	  ++nPhotonsPassingEta;
	  if (iPhoton->momentum.Et() > 45.0) ++nPassingLeadingET;
	  else if (iPhoton->momentum.Et() > 30.0) ++nPassingTrailingET;
	}
	if ((nPassingLeadingET >= 2) || 
	    ((nPassingLeadingET >= 1) && (nPassingTrailingET >= 1))) passETAndEta = true;
	++iPhoton;
      }
    }

    /*count how many photons pass the offline cuts + the trigger, for both ways of counting 
      triggers*/
    if (nPhotonsPassingEta >= 2) {
      ++nPassingEta;
      if (passOR) ++nPassingEtaOR;
      if (HLTPath != "") ++nPassingEtaLowestUnprescaled;
      if (passOR && (HLTPath == "")) {
	passEtaORNotLowestUnprescaledRun.push_back(event->runNumber);
	passEtaORNotLowestUnprescaledEvt.push_back(event->eventNumber);
	passEtaORNotLowestUnprescaledLumiSection.push_back(event->luminosityBlockNumber);
      }
      if (!passOR && (HLTPath != "")) {
	passEtaLowestUnprescaledNotORRun.push_back(event->runNumber);
	passEtaLowestUnprescaledNotOREvt.push_back(event->eventNumber);
	passEtaLowestUnprescaledNotORLumiSection.push_back(event->luminosityBlockNumber);
      }
    }
    if (passETAndEta) {

      //plot leading and trailing photon ET
      float leadingET = 0.0;
      float trailingET = 0.0;
      if (iPhotonMap != photonMap.end()) {
	susy::PhotonCollection photons = iPhotonMap->second;
	for (susy::PhotonCollection::const_iterator iPhoton = photons.begin(); 
	     iPhoton != photons.end(); ++iPhoton) {
	  if (iPhoton->momentum.Et() > leadingET) {
	    trailingET = leadingET;
	    leadingET = iPhoton->momentum.Et();
	  }
	  else if (iPhoton->momentum.Et() > trailingET) trailingET = iPhoton->momentum.Et();
	}
	leadingPhotonET->Fill(leadingET);
	trailingPhotonET->Fill(trailingET);
      }
      ++nPassingET;
      if (passOR) ++nPassingETOR;
      if (HLTPath != "") ++nPassingETLowestUnprescaled;
      if (passOR && (HLTPath == "")) {

	//save run/event/lumi information on these odd cases
	passETORNotLowestUnprescaledRun.push_back(event->runNumber);
	passETORNotLowestUnprescaledEvt.push_back(event->eventNumber);
	passETORNotLowestUnprescaledLumiSection.push_back(event->luminosityBlockNumber);

	/*leading and trailing photon ET -- are there events in the turn-on of the unprescaled 
	  trigger?*/
	leadingPhotonETPassingETORNotLowestUnprescaled->Fill(leadingET);
	trailingPhotonETPassingETORNotLowestUnprescaled->Fill(trailingET);

	//HLT status for these odd cases
	susy::TriggerMap HLTMap;
	for (std::vector<TString>::const_iterator iDesiredHLT = hltNames.begin(); 
	     iDesiredHLT != hltNames.end(); ++iDesiredHLT) {
	  susy::TriggerMap::const_iterator iHLT = event->hltMap.begin();
	  bool foundTrigger = false;
	  while ((iHLT != event->hltMap.end()) && !foundTrigger) {
	    if (iHLT->first.Contains(*iDesiredHLT)) foundTrigger = true;
	    else ++iHLT;
	  }
	  if (iHLT != event->hltMap.end()) {
	    HLTMap.insert(pair<TString, 
			  pair<Int_t, UChar_t> >(iHLT->first, 
						 pair<Int_t, UChar_t>(iHLT->second.first, 
								      iHLT->second.second)));
	  }
	}
	passETORNotLowestUnprescaledHLTMap.push_back(HLTMap);
      }
      if (!passOR && (HLTPath != "")) {
	passETLowestUnprescaledNotORRun.push_back(event->runNumber);
	passETLowestUnprescaledNotOREvt.push_back(event->eventNumber);
	passETLowestUnprescaledNotORLumiSection.push_back(event->luminosityBlockNumber);
      }
    }
  }

  //print counters
  out << "No. passing JSON: " << nPassingJSON << std::endl;
  out << "--> No. passing OR HLT selection: " << nPassingOR << std::endl;
  out << "--> No. passing lowest unprescaled HLT selection: " << nPassingLowestUnprescaled;
  out << std::endl;
  out << "--> No. passing eta: " << nPassingEta << std::endl;
  out << "  --> No. passing OR HLT selection: " << nPassingEtaOR << std::endl;
  out << "  --> No. passing lowest unprescaled HLT selection: " << nPassingEtaLowestUnprescaled;
  out << std::endl;
  out << "  --> No. passing ET: " << nPassingET << std::endl;
  out << "    --> No. passing OR HLT selection: " << nPassingETOR << std::endl;
  out << "    --> No. passing lowest unprescaled HLT selection: " << nPassingETLowestUnprescaled;
  out << std::endl << std::endl;

  //print trigger info
  out << "Events passing eta OR HLT selection but not lowest unprescaled HLT selection:\n";
  out << "Run | Event | Lumi section\n";
  for (std::vector<Int_t>::const_iterator i = passEtaORNotLowestUnprescaledRun.begin(); 
       i < passEtaORNotLowestUnprescaledRun.end(); ++i) {
    const unsigned int index = i - passEtaORNotLowestUnprescaledRun.begin();
    out << *i << " | " << passEtaORNotLowestUnprescaledEvt[index] << " | ";
    out << passEtaORNotLowestUnprescaledLumiSection[index] << std::endl;
  }
  out << "Events passing eta lowest unprescaled HLT selection but not OR HLT selection:\n";
  out << "Run | Event | Lumi section\n";
  for (std::vector<Int_t>::const_iterator i = passEtaLowestUnprescaledNotORRun.begin(); 
       i < passEtaLowestUnprescaledNotORRun.end(); ++i) {
    const unsigned int index = i - passEtaLowestUnprescaledNotORRun.begin();
    out << *i << " | " << passEtaLowestUnprescaledNotOREvt[index] << " | ";
    out << passEtaLowestUnprescaledNotORLumiSection[index] << std::endl;
  }
  out << "Events passing ET OR HLT selection but not lowest unprescaled HLT selection:\n";
  out << "Run | Event | Lumi section\n";
  out << "HLT status\n";
  for (std::vector<Int_t>::const_iterator i = passETORNotLowestUnprescaledRun.begin(); 
       i < passETORNotLowestUnprescaledRun.end(); ++i) {
    const unsigned int index = i - passETORNotLowestUnprescaledRun.begin();
    out << *i << " | " << passETORNotLowestUnprescaledEvt[index] << " | ";
    out << passETORNotLowestUnprescaledLumiSection[index] << std::endl;
    for (susy::TriggerMap::const_iterator iHLT = 
	   passETORNotLowestUnprescaledHLTMap[index].begin(); 
	 iHLT != passETORNotLowestUnprescaledHLTMap[index].end(); ++iHLT) {
      out << iHLT->first << ": prescale " << iHLT->second.first << ", ";
      out << ((iHLT->second.second > 0) ? "passed" : "failed") << std::endl;
    }
  }
  out << "Events passing ET lowest unprescaled HLT selection but not OR HLT selection:\n";
  out << "Run | Event | Lumi section\n";
  for (std::vector<Int_t>::const_iterator i = passETLowestUnprescaledNotORRun.begin(); 
       i < passETLowestUnprescaledNotORRun.end(); ++i) {
    const unsigned int index = i - passETLowestUnprescaledNotORRun.begin();
    out << *i << " | " << passETLowestUnprescaledNotOREvt[index] << " | ";
    out << passETLowestUnprescaledNotORLumiSection[index] << std::endl;
  }
  out.close();

  //write ROOT file
  triggerPassedVsRunRange->Write();
  leadingPhotonET->Write();
  trailingPhotonET->Write();
  leadingPhotonETPassingETORNotLowestUnprescaled->Write();
  trailingPhotonETPassingETORNotLowestUnprescaled->Write();
  outROOT.Write();
  outROOT.Close();
  delete triggerPassedVsRunRange;
}

void EventAnalyzer::Loop(TTree* outTree, Categorizer categorizer, susy::Category* pCategory, 
			 susy::Event* pEvent) {
  if (fChain == 0) return;

  //remove this memory-intensive step
  // to check duplicated events
  // std::map<int, std::set<int> > allEvents;

  // start event looping
  Long64_t nentries = fChain->GetEntriesFast();
  if ((processNEvents < 0) || (processNEvents > nentries)) processNEvents = nentries;
  Long64_t nEvtsProcessed = 0;
  Long64_t nbytes = 0, nb = 0;
  unsigned int nEvtsPassingJSON = 0;
  unsigned int nEvtsPassingHLT = 0;
  unsigned int nEvtsPassingGoodPV = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {
    event->Init();
    pCategory->reset();
    if (pEvent != NULL) pEvent->Init();
    nEvtsProcessed = jentry + 1;
    if(printLevel > 0) std::cout << "Get the tree contents." << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) {
      nEvtsProcessed = jentry;
      break;
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ((printLevel > 0) || ((printInterval > 0) && ((jentry >= printInterval) && 
						     (jentry % printInterval == 0)))) {
      cout << jentry << " events processed with run = " << event->runNumber;
      cout << ", event = " << event->eventNumber << ", lumi section = ";
      cout << event->luminosityBlockNumber << std::endl;
    }

    if(printLevel > 0) std::cout << "Initialize any global variables to be reset per event." << std::endl;
    InitializePerEvent();

    //uncomment if you want to apply JSON criteria
    if(printLevel > 0) std::cout << "Apply good run list." << std::endl;
    bool JSON = isInJson(event->runNumber,event->luminosityBlockNumber);
    if (JSON) ++nEvtsPassingJSON;

    //uncomment if you want to apply HLT criteria
    bool passHLT = (useTrigger ? PassTriggers() : true);
    if (JSON && passHLT) ++nEvtsPassingHLT;

    //decide if the event passed the PV requirement
    bool foundGoodPV = false;
    vector<susy::Vertex>::const_iterator iPV = event->vertices.begin();
    while ((iPV != event->vertices.end()) && !foundGoodPV) {
      if (!iPV->isFake() && (iPV->ndof > 4) && (iPV->position.Z() <= 24.0/*cm*/) && 
	  (iPV->position.Perp() <= 2.0/*cm*/)) foundGoodPV = true;
      else ++iPV;
    }
    if (JSON && passHLT && foundGoodPV) ++nEvtsPassingGoodPV;

    //filter criteria
    bool filter = !enableFilter || (JSON && passHLT && foundGoodPV);

    // uncomment this to use the Json file to flag good data (or bad depending on your outlook)    
    if (!JSON) continue;

    // uncomment this to print all ntuple variables
//     Print(*event);

//     // remove this memory-intensive step
//     if(printLevel > 0) std::cout << "Check duplicated events for data only." << std::endl;
//     bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
//     if(event->isRealData && duplicateEvent) {
//       if (filter) {
// 	if (pEvent != NULL) *pEvent = *event;
// 	outTree->Fill();
//       }
//       continue;
//     }

    //photon variables to pass to the Categorizer object
    VDOUBLE photonET;
    VDOUBLE photonEta;
    VDOUBLE photonECALIso;
    VDOUBLE PUSubtractedPhotonECALIso;
    VDOUBLE photonHCALIso;
    VDOUBLE PUSubtractedPhotonHCALIso;
    VDOUBLE photonHOverE;
    VDOUBLE photonR9;
    VDOUBLE photonTrackIso;
    VDOUBLE photonSigmaIetaIeta;
    VDOUBLE photonSeedTime;
    VDOUBLE photonE2OverE9;
    VINT photonSeedIeta;
    VDOUBLE photonPhi;
    VBOOL photonHasPixelSeed;
    VBOOL photonPasses;
    TLorentzVector diEMP4;
    unsigned int passingPhotonCount = 0;
    double evtInvMass = -1.0;

    // classify photon objects
    std::map<TString,susy::PhotonCollection> photonMap = event->photons;
    std::map<TString,susy::PhotonCollection>::const_iterator iPhotonMap = 
      photonMap.find((const char*)tag_);
    int category = FAIL;
    if (iPhotonMap != photonMap.end()) {
      susy::PhotonCollection photons = iPhotonMap->second;
      for (susy::PhotonCollection::const_iterator iPhoton = photons.begin(); 
	   iPhoton != photons.end(); ++iPhoton) {

	//get isolations
	double ECALIsoOffline = -1.0;
	double HCALIsoOffline = -1.0;
	double trackIso = -1.0;
	double ECALIsoHLT = -1.0;
	double HCALIsoHLT = -1.0;
	switch (isoConeOffline_) {
	case DR03:
	  ECALIsoOffline = iPhoton->ecalRecHitSumEtConeDR03 - (photonECALIsoEffArea_*event->rho);
	  HCALIsoOffline = iPhoton->hcalTowerSumEtConeDR03() - (photonHCALIsoEffArea_*event->rho);
	  trackIso = iPhoton->trkSumPtHollowConeDR03;
	  break;
	case DR04:
	  ECALIsoOffline = iPhoton->ecalRecHitSumEtConeDR04 - (photonECALIsoEffArea_*event->rho);
	  HCALIsoOffline = iPhoton->hcalTowerSumEtConeDR04() - (photonHCALIsoEffArea_*event->rho);
	  trackIso = iPhoton->trkSumPtHollowConeDR04;
	  break;
	default:
	  cerr << "Error: invalid cone size " << isoConeOffline_ << ".  Isolation requirements ";
	  cerr << "will be ignored.\n";
	}
	switch (isoConeHLT_) {
	case DR03:
	  ECALIsoHLT = iPhoton->ecalRecHitSumEtConeDR03;
	  HCALIsoHLT = iPhoton->hcalTowerSumEtConeDR03();
	  break;
	case DR04:
	  ECALIsoHLT = iPhoton->ecalRecHitSumEtConeDR04;
	  HCALIsoHLT = iPhoton->hcalTowerSumEtConeDR04();
	  break;
	default:
	  cerr << "Error: invalid cone size " << isoConeHLT_ << ".  Isolation requirements will ";
	  cerr << "be ignored.\n";
	}

	//fill photon variable vectors
	photonET.push_back(iPhoton->momentum.Et());
	photonEta.push_back(iPhoton->caloPosition.Eta());
	photonECALIso.push_back(ECALIsoHLT);
	PUSubtractedPhotonECALIso.push_back(ECALIsoOffline);
	photonHCALIso.push_back(HCALIsoHLT);
	PUSubtractedPhotonHCALIso.push_back(HCALIsoOffline);
	photonHOverE.push_back((double)iPhoton->hadronicOverEm);
	photonR9.push_back((double)iPhoton->r9);
	photonTrackIso.push_back(trackIso);
	photonSigmaIetaIeta.push_back((double)iPhoton->sigmaIetaIeta);
	photonSeedTime.push_back((double)iPhoton->seedTime);
	photonE2OverE9.push_back((double)((iPhoton->e1x2)/(iPhoton->e3x3)));
	photonSeedIeta.push_back(99); //obsolete
	photonPhi.push_back(iPhoton->caloPosition.Phi());
	photonHasPixelSeed.push_back(iPhoton->nPixelSeeds > 0);
      }

      //categorize the event
      categorizer.setPhotonET(photonET);
      categorizer.setPhotonEta(photonEta);
      categorizer.setPhotonECALIso(photonECALIso);
      categorizer.setPUSubtractedPhotonECALIso(PUSubtractedPhotonECALIso);
      categorizer.setPhotonHCALIso(photonHCALIso);
      categorizer.setPUSubtractedPhotonHCALIso(PUSubtractedPhotonHCALIso);
      categorizer.setPhotonHOverE(photonHOverE);
      categorizer.setPhotonR9(photonR9);
      categorizer.setPhotonTrackIso(photonTrackIso);
      categorizer.setPhotonSigmaIetaIeta(photonSigmaIetaIeta);
      categorizer.setPhotonSeedTime(photonSeedTime);
      categorizer.setPhotonE2OverE9(photonE2OverE9);
      categorizer.setPhotonPhi(photonPhi);
      categorizer.setPhotonHasPixelSeed(photonHasPixelSeed);
      categorizer.checkInput();
      if (!categorizer.initialized()) {
	STRINGSTREAM err;
	err << "Error: Categorizer object not properly initialized at entry " << jentry;
	err << ", run " << event->runNumber << ", event " << event->eventNumber;
	err << ", lumi section " << event->luminosityBlockNumber << ".\n";
	throw err.str();
      }
      try {
	categorizer.decideAll();
	categorizer.findPassingPhotons();
	categorizer.classify();

	//find the 2 photons that decided the event and the di-EM Lorentz vector
	VINT passingPhotons = categorizer.getPassingPhotons();
	for (susy::PhotonCollection::const_iterator iPhoton = photons.begin(); 
	     iPhoton != photons.end(); ++iPhoton) {
	  if (((iPhoton - photons.begin()) == passingPhotons[0]) || 
	      ((iPhoton - photons.begin()) == passingPhotons[1])) {
	    photonPasses.push_back(true);
	    diEMP4+=iPhoton->momentum;
	    ++passingPhotonCount;
	  }
	  else photonPasses.push_back(false);
	}
      }
      catch (STRING& badInput) { throw; }

      if(printLevel > 0) std::cout << "Apply trigger selection in the event." << std::endl;

      //write the category information to the event
      try {
	for (susy::PhotonCollection::const_iterator iPhoton = photons.begin(); 
	     iPhoton != photons.end(); ++iPhoton) {
	  const unsigned int photonIndex = iPhoton - photons.begin();
	  pCategory->setPhotonType(tag_, photonIndex, categorizer.getPhotonType()[photonIndex]);
	  pCategory->setPassETMin1(tag_, photonIndex, 
				   categorizer.getPhotonPassETMin1()[photonIndex]);
	  pCategory->setPassETMin2(tag_, photonIndex, 
				   categorizer.getPhotonPassETMin2()[photonIndex]);
	  pCategory->setPassAbsEtaMax(tag_, photonIndex, 
				      categorizer.getPhotonPassAbsEtaMax()[photonIndex]);
	  pCategory->setPassECALIsoMax(tag_, photonIndex, 
				       categorizer.getPhotonPassECALIsoMax()[photonIndex]);
	  pCategory->
	    setPassPUSubtractedECALIsoMax(tag_, photonIndex, 
					  categorizer.
					  getPhotonPassPUSubtractedECALIsoMax()[photonIndex]);
	  pCategory->setPassHCALIsoMax(tag_, photonIndex, 
				       categorizer.getPhotonPassHCALIsoMax()[photonIndex]);
	  pCategory->
	    setPassPUSubtractedHCALIsoMax(tag_, photonIndex, 
					  categorizer.
					  getPhotonPassPUSubtractedHCALIsoMax()[photonIndex]);
	  pCategory->setPassHOverEMax(tag_, photonIndex, 
				      categorizer.getPhotonPassHOverEMax()[photonIndex]);
	  pCategory->setPassR9Max(tag_, photonIndex, 
				  categorizer.getPhotonPassR9Max()[photonIndex]);
	  pCategory->setPassR9Min(tag_, photonIndex, 
				  categorizer.getPhotonPassR9Min()[photonIndex]);
	  pCategory->setPassTrackIsoMax(tag_, photonIndex, 
					categorizer.getPhotonPassTrackIsoMax()[photonIndex]);
	  pCategory->setPassCombinedIsoMax(tag_, photonIndex, 
					   categorizer.getPhotonPassCombinedIsoMax()[photonIndex]);
	  pCategory->
	    setPassFakeCombinedIsoMax(tag_, photonIndex, 
				      categorizer.getPhotonPassFakeCombinedIsoMax()[photonIndex]);
	  pCategory->
	    setPassSigmaIetaIetaMax(tag_, photonIndex, 
				    categorizer.getPhotonPassSigmaIetaIetaMax()[photonIndex]);
	  pCategory->
	    setPassHLTSigmaIetaIetaMax(tag_, photonIndex, 
				       categorizer.getPhotonPassHLTSigmaIetaIetaMax()[photonIndex]);
	  pCategory->setPassAbsSeedTimeMax(tag_, photonIndex, 
					   categorizer.getPhotonPassAbsSeedTimeMax()[photonIndex]);
	  pCategory->setPassE2OverE9Max(tag_, photonIndex, 
					categorizer.getPhotonPassE2OverE9Max()[photonIndex]);
	  pCategory->setHasPixelSeed(tag_, photonIndex, 
				     categorizer.getPhotonHasPixelSeed()[photonIndex]);
	  pCategory->setPassPreselection(tag_, photonIndex, 
					 categorizer.getPhotonPassPreselection()[photonIndex]);
	  pCategory->setIsDeciding(tag_, photonIndex, photonPasses[photonIndex]);
	  pCategory->setPhotonSeedTime(tag_, photonIndex, photonSeedTime[photonIndex]);
	  pCategory->setPhotonE2OverE9(tag_, photonIndex, photonE2OverE9[photonIndex]);
	  pCategory->setPhotonSeedIeta(tag_, photonIndex, photonSeedIeta[photonIndex]);
	}
	category = categorizer.getCategory();
	pCategory->setEventCategory(tag_, category);
	pCategory->setPassDPhiMin(tag_, categorizer.getEvtPassDPhiMin());
	pCategory->setPassDRMin(tag_, categorizer.getEvtPassDRMin());
	pCategory->setPassAsymmetricETMin(tag_, categorizer.getEvtPassAsymmetricETMin());
	pCategory->setEvtDiEMET(tag_, categorizer.getEvtDiEMET());
	if (passingPhotonCount > 2) {
	  STRINGSTREAM err;
	  err << "Error: " << passingPhotonCount << " passing photons.\n";
	  throw err.str();
	}
	if (passingPhotonCount == 2) evtInvMass = diEMP4.M();
	pCategory->setEvtInvMass(tag_, evtInvMass);
	pCategory->setPassGoodPV(tag_, foundGoodPV);
      }
      catch (STRING& badInput) { throw; }
    }

    //fill the tree
    if (filter) {
      if (pEvent != NULL) *pEvent = *event;
      outTree->Fill();
    }

    //store run, event, and lumi section info for events passing filter criteria (JSON and HLT)
    if (filter) {
      RUNEVTLUMIPAIR runEvtLumiPair(RUNEVTPAIR(event->runNumber, event->eventNumber), 
				    event->luminosityBlockNumber);
      RUNEVTLUMIMASSPAIR runEvtLumiMassPair(runEvtLumiPair, evtInvMass);
      switch (category) {
      case GG:
	ggEvts_.insert(runEvtLumiPair);
	break;
      case FF:
	ffEvts_.insert(runEvtLumiPair);
	break;
      case EG:
	egEvts_.insert(runEvtLumiPair);
	break;
      case EE:
	eeEvts_.insert(runEvtLumiMassPair);
	break;
      default:
	break;
      }
    }

    //print debug info
    if (printLevel == 2) debugPrint(categorizer, event);

  } // for jentry

  //print JSON, HLT, and good vertex info
  cout << "No. events total: " << nEvtsProcessed << endl;
  cout << "No. events passing JSON: " << nEvtsPassingJSON << endl;
  cout << "No. events passing HLT: " << nEvtsPassingHLT << endl;
  cout << "No. events passing good vertex: " << nEvtsPassingGoodPV << endl;
}

void EventAnalyzer::debugPrint(const Categorizer& categorizer, susy::Event* evt) const
{
  //don't print unless event processing is completely finished
  if (!categorizer.done()) {
    std::cerr << "Error: finish event processing.\n";
    return;
  }

  //get all event information
  VDOUBLE photonET = categorizer.getPhotonET();
  VDOUBLE photonEta = categorizer.getPhotonEta();
  VDOUBLE photonECALIso = categorizer.getPhotonECALIso();
  VDOUBLE photonHCALIso = categorizer.getPhotonHCALIso();
  VDOUBLE photonHOverE = categorizer.getPhotonHOverE();
  VDOUBLE photonR9 = categorizer.getPhotonR9();
  VDOUBLE photonTrackIso = categorizer.getPhotonTrackIso();
  VDOUBLE photonSigmaIetaIeta = categorizer.getPhotonSigmaIetaIeta();
  VDOUBLE photonSeedTime = categorizer.getPhotonSeedTime();
  VDOUBLE photonE2OverE9 = categorizer.getPhotonE2OverE9();
  VDOUBLE photonPhi = categorizer.getPhotonPhi();
  VBOOL photonHasPixelSeed = categorizer.getPhotonHasPixelSeed();
  double photon1ETMin = categorizer.getPhoton1ETMin();
  double photon2ETMin = categorizer.getPhoton2ETMin();
  double photonAbsEtaMax = categorizer.getPhotonAbsEtaMax();
  double photonECALIsoMaxPTMultiplier = categorizer.getPhotonECALIsoMaxPTMultiplier();
  double photonECALIsoMaxConstant = categorizer.getPhotonECALIsoMaxConstant();
  double photonHCALIsoMaxPTMultiplier = categorizer.getPhotonHCALIsoMaxPTMultiplier();
  double photonHCALIsoMaxConstant = categorizer.getPhotonHCALIsoMaxConstant();
  double photonHOverEMax = categorizer.getPhotonHOverEMax();
  double photonR9Max = categorizer.getPhotonR9Max();
  double photonR9Min = categorizer.getPhotonR9Min();
  double photonTrackIsoMaxPTMultiplier = categorizer.getPhotonTrackIsoMaxPTMultiplier();
  double photonTrackIsoMaxConstant = categorizer.getPhotonTrackIsoMaxConstant();
  double photonSigmaIetaIetaMax = categorizer.getPhotonSigmaIetaIetaMax();
  double photonHLTSigmaIetaIetaMax = categorizer.getPhotonHLTSigmaIetaIetaMax();
  double photonAbsSeedTimeMax = categorizer.getPhotonAbsSeedTimeMax();
  double photonE2OverE9Max = categorizer.getPhotonE2OverE9Max();
  double photonDPhiMin = categorizer.getPhotonDPhiMin();
  double photonDRMin = categorizer.getPhotonDRMin();
  VBOOL photonPassETMin1 = categorizer.getPhotonPassETMin1();
  VBOOL photonPassETMin2 = categorizer.getPhotonPassETMin2();
  VBOOL photonPassAbsEtaMax = categorizer.getPhotonPassAbsEtaMax();
  VBOOL photonPassECALIsoMax = categorizer.getPhotonPassECALIsoMax();
  VBOOL photonPassHCALIsoMax = categorizer.getPhotonPassHCALIsoMax();
  VBOOL photonPassHOverEMax = categorizer.getPhotonPassHOverEMax();
  VBOOL photonPassR9Max = categorizer.getPhotonPassR9Max();
  VBOOL photonPassR9Min = categorizer.getPhotonPassR9Min();
  VBOOL photonPassTrackIsoMax = categorizer.getPhotonPassTrackIsoMax();
  VBOOL photonPassCombinedIsoMax = categorizer.getPhotonPassCombinedIsoMax();
  VBOOL photonPassFakeCombinedIsoMax = categorizer.getPhotonPassFakeCombinedIsoMax();
  VBOOL photonPassSigmaIetaIetaMax = categorizer.getPhotonPassSigmaIetaIetaMax();
  VBOOL photonPassHLTSigmaIetaIetaMax = categorizer.getPhotonPassHLTSigmaIetaIetaMax();
  VBOOL photonPassAbsSeedTimeMax = categorizer.getPhotonPassAbsSeedTimeMax();
  VBOOL photonPassE2OverE9Max = categorizer.getPhotonPassE2OverE9Max();
  VBOOL photonPassPreselection = categorizer.getPhotonPassPreselection();
  bool evtPassDPhiMin = categorizer.getEvtPassDPhiMin();
  bool evtPassDRMin = categorizer.getEvtPassDRMin();
  bool evtPassAsymmetricETMin = categorizer.getEvtPassAsymmetricETMin();
  VINT passingPhotons = categorizer.getPassingPhotons();
  VINT photonType = categorizer.getPhotonType();
  int category = categorizer.getCategory();

  //loop over photons
  STRINGSTREAM debug;
  debug << "*****************\n";
  debug << "Run " << evt->runNumber << ", event " << evt->eventNumber << ", lumi section ";
  debug << evt->luminosityBlockNumber << std::endl;
  for (VDOUBLE_IT iET = photonET.begin(); iET != photonET.end(); ++iET) {

    //print photon quantity, cut value, and pass flag
    debug << "%%%%%%%%%%%%%%%%%\n";
    const unsigned int i = iET - photonET.begin();
    debug << "Photon index: " << i << std::endl;
    debug << "-----------------\n";
    debug << "Photon ET: " << *iET << " GeV\n";
    debug << "Cut: photon 1 ET > " << photon1ETMin << " GeV\n";
    debug << "Result: ";
    if (photonPassETMin1[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "Cut: photon 2 ET > " << photon2ETMin << " GeV\n";
    debug << "Result: ";
    if (photonPassETMin2[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon eta: " << photonEta[i] << std::endl;
    debug << "Cut: photon |eta| < " << photonAbsEtaMax << std::endl;
    debug << "Result: ";
    if (photonPassAbsEtaMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon ECAL isolation: " << photonECALIso[i] << " GeV\n";
    debug << "Cut: photon ECAL isolation < (" << photonECALIsoMaxPTMultiplier << "ET + ";
    debug << photonECALIsoMaxConstant << " GeV) = ";
    debug << (photonECALIsoMaxPTMultiplier*(*iET) + photonECALIsoMaxConstant) << " GeV\n";
    debug << "Result: ";
    if (photonPassECALIsoMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon HCAL isolation: " << photonHCALIso[i] << " GeV\n";
    debug << "Cut: photon HCAL isolation < (" << photonHCALIsoMaxPTMultiplier << "ET + ";
    debug << photonHCALIsoMaxConstant << " GeV) = ";
    debug << (photonHCALIsoMaxPTMultiplier*(*iET) + photonHCALIsoMaxConstant) << " GeV\n";
    debug << "Result: ";
    if (photonPassHCALIsoMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon H/E: " << photonHOverE[i] << std::endl;
    debug << "Cut: photon H/E < " << photonHOverEMax << std::endl;
    debug << "Result: ";
    if (photonPassHOverEMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon R9: " << photonR9[i] << std::endl;
    debug << "Cut: photon R9 < " << photonR9Max << std::endl;
    debug << "Result: ";
    if (photonPassR9Max[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "Cut: photon R9 > " << photonR9Min << std::endl;
    debug << "Result: ";
    if (photonPassR9Min[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon track isolation: " << photonTrackIso[i] << " GeV\n";
    debug << "Cut: photon track isolation < (" << photonTrackIsoMaxPTMultiplier << "ET + ";
    debug << photonTrackIsoMaxConstant << " GeV) = ";
    debug << (photonTrackIsoMaxPTMultiplier*(*iET) + photonTrackIsoMaxConstant) << " GeV\n";
    debug << "Result: ";
    if (photonPassTrackIsoMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon combined isolation: ";
    debug << (photonECALIso[i] + photonHCALIso[i] + photonTrackIso[i]) << " GeV\n";
    debug << "Cut: photon combined isolation < 6 GeV\n";
    debug << "Result: ";
    if (photonPassCombinedIsoMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Cut: fake combined isolation < 12 GeV\n";
    debug << "Result: ";
    if (photonPassFakeCombinedIsoMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon sigmaIetaIeta: " << photonSigmaIetaIeta[i] << std::endl;
    debug << "Cut: photon sigmaIetaIeta < " << photonSigmaIetaIetaMax << std::endl;
    debug << "Result: ";
    if (photonPassSigmaIetaIetaMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "Cut: photon sigmaIetaIeta < " << photonHLTSigmaIetaIetaMax << std::endl;
    debug << "Result: ";
    if (photonPassHLTSigmaIetaIetaMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon seed time: " << photonSeedTime[i] << " ns\n";
    debug << "Cut: photon |seed time| < " << photonAbsSeedTimeMax << " ns\n";
    debug << "Result: ";
    if (photonPassAbsSeedTimeMax[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Photon E2/E9: " << photonE2OverE9[i] << std::endl;
    debug << "Cut: photon E2/E9 < " << photonE2OverE9Max << std::endl;
    debug << "Result: ";
    if (photonPassE2OverE9Max[i]) debug << "pass\n";
    else debug << "fail\n";
    debug << "-----------------\n";
    debug << "Result: ";
    if (photonPassPreselection[i]) debug << "photon passed preselection\n";
    else debug << "photon failed preselection\n";
    debug << "-----------------\n";
    debug << "Photon has pixel seed: ";
    if (photonHasPixelSeed[i]) debug << "yes\n";
    else debug << "no\n";
    debug << "-----------------\n";
    debug << "Photon type: ";
    switch (photonType[i]) {
    case FAIL:
      debug << "fail\n";
      break;
    case G:
      debug << "photon\n";
      break;
    case E:
      debug << "electron\n";
      break;
    case F:
      debug << "fake\n";
      break;
    default:
      debug << "unknown\n";
      break;
    }
    debug << "%%%%%%%%%%%%%%%%%\n";
  }

  //print event quantities
  debug << "Passing photon indices: " << passingPhotons[0] << ", " << passingPhotons[1];
  debug << std::endl;
  debug << "#################\n";
  debug << "Passing photon phi: ";
  if (passingPhotons[0] >= 0) debug << photonPhi[passingPhotons[0]] << ", ";
  else debug << "N/A" << ", ";
  if (passingPhotons[1] >= 0) debug << photonPhi[passingPhotons[1]] << std::endl;
  else debug << "N/A" << std::endl;
  debug << "Cut: dPhi(photon, photon) > " << photonDPhiMin << std::endl;
  debug << "Result: ";
  if (evtPassDPhiMin) debug << "pass\n";
  else debug << "fail\n";
  debug << "#################\n";
  debug << "Cut: dR(photon, photon) > " << photonDRMin << std::endl;
  debug << "Result: ";
  if (evtPassDRMin) debug << "pass\n";
  else debug << "fail\n";
  debug << "#################\n";
  debug << "Passing photon ET: ";
  if (passingPhotons[0] >= 0) debug << photonET[passingPhotons[0]] << ", ";
  else debug << "N/A" << ", ";
  if (passingPhotons[1] >= 0) debug << photonET[passingPhotons[1]] << std::endl;
  else debug << "N/A" << std::endl;
  debug << "Cut: ET(photon1) > " << photon1ETMin << ", ET(photon2) > " << photon2ETMin << std::endl;
  debug << "Result: ";
  if (evtPassAsymmetricETMin) debug << "pass\n";
  else debug << "fail\n";
  debug << "#################\n";
  debug << "Event category: ";
  switch (category) {
  case FAIL:
    debug << "fail\n";
    break;
  case GG:
    debug << "candidate\n";
    break;
  case EE:
    debug << "ee\n";
    break;
  case FF:
    debug << "ff\n";
    break;
  case EG:
    debug << "eg\n";
    break;
  default:
    debug << "unknown\n";
    break;
  }
  debug << "*****************\n";

  //log
  std::cerr << debug.str();
}
