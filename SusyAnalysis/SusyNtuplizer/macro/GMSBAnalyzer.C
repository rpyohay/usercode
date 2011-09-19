#define GMSBAnalyzer_cxx
#include "GMSBAnalyzer.h"
#include "../../../GMSBTools/Filters/interface/Categorizer.h"
#include "../../../DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void GMSBAnalyzer::Loop(const std::string& outputFile)
{
//   In a ROOT session, you can do:
//      Root > .L GMSBAnalyzer.C
//      Root > GMSBAnalyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //open file
   TFile out(outputFile.c_str(), "RECREATE");
   out.cd();

   //efficiency histograms
   TH1F allBkgVsNPV("allBkgVsNPV", "", 20, 0.5, 20.5);
   TH1F passingBkgVsNPV("passingBkgVsNPV", "", 20, 0.5, 20.5);

   //set user-specified number of entries to process
   Long64_t nentries = fChain->GetEntriesFast();
   if (nEvts_ != -1) nentries = nEvts_;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //loop over photons
      std::map<TString,susy::PhotonCollection>::const_iterator iPhotonMap = 
	susyEvent->photons.find((const char*)tag_);
      if (iPhotonMap != susyEvent->photons.end()) {
	susy::PhotonCollection photons = iPhotonMap->second;
	for (susy::PhotonCollection::const_iterator iPhoton = photons.begin(); 
	     iPhoton != photons.end(); ++iPhoton) {
	  const unsigned int i = iPhoton - photons.begin();

	  std::cout << "Photon " << i << std::endl;

	  //does the photon pass acceptance criteria, and is it matched to a gen quark or gluon?
	  VINT matchPDGIDs;
	  matchPDGIDs.push_back(1);
	  matchPDGIDs.push_back(2);
	  matchPDGIDs.push_back(3);
	  matchPDGIDs.push_back(4);
	  matchPDGIDs.push_back(5);
	  matchPDGIDs.push_back(6);
	  matchPDGIDs.push_back(21);
	  if (passesDenominatorSelection(i, *iPhoton, matchPDGIDs)) {

	    //calculate number of good vertices in the event
	    const unsigned int nPV = numGoodVertices();

	    //denominator photon
	    allBkgVsNPV.Fill(nPV);

	    //does the photon pass the full selection?
	    if (passesNumeratorSelection(i)) {

	      //numerator photon
	      passingBkgVsNPV.Fill(nPV);
	    }
	  }
	}
      }
   }//for (Long64_t jentry=0; jentry<nentries;jentry++)

   //calculate efficiency
   TCanvas bkgEffVsNPVCanvas("bkgEffVsNPVCanvas", "", 600, 600);
   TGraphAsymmErrors bkgEffVsNPV(&passingBkgVsNPV, &allBkgVsNPV, "cp");

   //write histograms
   bkgEffVsNPVCanvas.cd();
   bkgEffVsNPV.Draw("AP");
   bkgEffVsNPV.Write();
   bkgEffVsNPVCanvas.Write();
   out.Write();
   out.Close();
}

bool GMSBAnalyzer::passesDenominatorSelection(const unsigned int photonIndex, 
					      const susy::Photon& photon, 
					      const VINT& matchPDGIDs) const
{
  return (susyCategory->getPassETMin1(tag_, photonIndex) && 
	  susyCategory->getPassAbsEtaMax(tag_, photonIndex) && 
	  matchedToGenParticle(photon, matchPDGIDs, 3));
}

bool GMSBAnalyzer::passesNumeratorSelection(const unsigned int photonIndex) const
{
  return (susyCategory->getPassECALIsoMax(tag_, photonIndex) && 
	  susyCategory->getPassHCALIsoMax(tag_, photonIndex) && 
	  susyCategory->getPassHOverEMax(tag_, photonIndex) && 
	  susyCategory->getPassR9Max(tag_, photonIndex) && 
	  susyCategory->getPassSigmaIetaIetaMax(tag_, photonIndex) && 
	  susyCategory->getPassTrackIsoMax(tag_, photonIndex));
}

bool GMSBAnalyzer::matchedToGenParticle(const susy::Photon& photon, const VINT& matchPDGIDs, 
					const UChar_t status) const
{
  unsigned int numGenQG = 0;

  //flag indicating whether a gen particle match has been found
  bool foundMatch = false;

  //loop over gen particles
//   std::vector<susy::Particle>::const_iterator iGenParticle = susyEvent->genParticles.begin();
//   while ((iGenParticle != susyEvent->genParticles.end()) && !foundMatch) {
  for (std::vector<susy::Particle>::const_iterator iGenParticle = susyEvent->genParticles.begin(); 
       iGenParticle != susyEvent->genParticles.end(); ++iGenParticle) {

    //loop over PDG IDs to match
    bool foundPDGIDMatch = false;
//     VINT_IT iMatchPDGID = matchPDGIDs.begin();
//     while ((iMatchPDGID != matchPDGIDs.end()) && !foundPDGIDMatch) {
    for (VINT_IT iMatchPDGID = matchPDGIDs.begin(); iMatchPDGID != matchPDGIDs.end(); 
	 ++iMatchPDGID) {

      //is the gen particle the right type and status, and is it DR-matched to the photon?
      if ((iGenParticle->pdgId == (Int_t)(*iMatchPDGID)) && (iGenParticle->status == status) && 
	  (reco::deltaR(iGenParticle->momentum.Eta(), iGenParticle->momentum.Phi(), 
			photon.momentum.Eta(), 
			photon.momentum.Phi()) < 0.5)) {

	++numGenQG;

	//break out of all loops once a match is found
	foundMatch = true;
	foundPDGIDMatch = true;
      }

      //increment PDG ID counter
//       ++iMatchPDGID;
    }

    //increment gen particle counter
//     ++iGenParticle;
  }

  std::cout << "numGenQG: " << numGenQG << std::endl;

  //return
  return foundMatch;
}

unsigned int GMSBAnalyzer::numGoodVertices() const
{
  unsigned int nPV = 0;
  for (std::vector<susy::Vertex>::const_iterator iPV = susyEvent->vertices.begin(); 
       iPV != susyEvent->vertices.end(); ++iPV) {
    if (!(iPV->isFake()) && (iPV->ndof > 4) && (iPV->position.z() <= 24.0/*cm*/) && 
	(iPV->position.Perp() <= 2.0/*cm*/)) ++nPV;
  }
  return nPV;
}

void GMSBAnalyzer::countEE(string& outputFile)
{
   if (fChain == 0) return;

   //open files
   TFile out(outputFile.c_str(), "RECREATE");
   out.cd();
   ofstream outTxt((outputFile.replace(outputFile.find("root"), 4, "txt")).c_str());

   //invariant mass histogram
   Float_t bins[8] = {0.0, 71.0, 76.0, 81.0, 101.0, 106.0, 111.0, 200.0};
   TH1F mee("mee", "", 7, bins);

   //set user-specified number of entries to process
   Long64_t nentries = fChain->GetEntriesFast();
   if (nEvts_ != -1) nentries = nEvts_;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //if this event is an ee event, what invariant mass bin does it fall into?
      if (susyCategory->getEventCategory(tag_) == EE) {
	const double invMass = susyCategory->getEvtInvMass(tag_);
	mee.Fill(invMass);
	outTxt << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
	if ((invMass >= 0.0) && (invMass < 71.0)) outTxt << "1\n";
	if ((invMass >= 71.0) && (invMass < 76.0)) outTxt << "2\n";
	if ((invMass >= 76.0) && (invMass < 81.0)) outTxt << "3\n";
	if ((invMass >= 81.0) && (invMass < 101.0)) outTxt << "4\n";
	if ((invMass >= 101.0) && (invMass < 106.0)) outTxt << "5\n";
	if ((invMass >= 106.0) && (invMass < 111.0)) outTxt << "6\n";
	if (invMass >= 111.0) outTxt << "7\n";
      }
   }

   //write histograms
   outTxt.close();
   mee.Write();
   out.Write();
   out.Close();
}

void GMSBAnalyzer::debugPrint(const unsigned int jentry) const
{
  //loop over photons
  cout << "*******************Event " << (jentry + 1) << "*******************\n\n";
  std::map<TString,susy::PhotonCollection>::const_iterator iPhotonMap = 
    susyEvent->photons.find((const char*)tag_);
  if (iPhotonMap != susyEvent->photons.end()) {
    susy::PhotonCollection photons = iPhotonMap->second;
    for (susy::PhotonCollection::const_iterator iPhoton = photons.begin(); 
	 iPhoton != photons.end(); ++iPhoton) {
      const unsigned int i = iPhoton - photons.begin();

      //print photon quantity value and value of the cut on that quantity
      cout << "Photon index: " << i << endl;
      cout << "-------ET   : " << iPhoton->momentum.Et() << endl;
      cout << "--------------";
      cout << (susyCategory->getPassETMin1(tag_, i) ? "PASSED LEADING ET\n" : "FAILED LEADING ET\n");
      cout << (susyCategory->getPassETMin2(tag_, i) ? "PASSED TRAILING ET\n" : "FAILED TRAILING ET\n");
      cout << "-------eta  : " << iPhoton->caloPosition.Eta() << endl;
      cout << "--------------";
      cout << (susyCategory->getPassAbsEtaMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------IECAL: " << iPhoton->ecalRecHitSumEtConeDR04 << endl;
      cout << "--------------";
      cout << (susyCategory->getPassECALIsoMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------IHCAL: " << iPhoton->hcalTowerSumEtConeDR04() << endl;
      cout << "--------------";
      cout << (susyCategory->getPassHCALIsoMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------H/E  : " << iPhoton->hadronicOverEm << endl;
      cout << "--------------";
      cout << (susyCategory->getPassHOverEMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------R9   : " << iPhoton->r9 << endl;
      cout << "--------------";
      cout << (susyCategory->getPassR9Max(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------t    : " << iPhoton->seedTime << endl;
      cout << "--------------";
      cout << (susyCategory->getPassAbsSeedTimeMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------E2/E9: " << (double)((iPhoton->e1x2)/(iPhoton->e3x3)) << endl;
      cout << "--------------";
      cout << (susyCategory->getPassE2OverE9Max(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "----------------------------";
      cout << (susyCategory->getPassPreselection(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------ITRK : " << iPhoton->trkSumPtHollowConeDR04 << endl;
      cout << "--------------";
      cout << (susyCategory->getPassTrackIsoMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------sieie: " << iPhoton->sigmaIetaIeta << endl;
      cout << "--------------";
      cout << (susyCategory->getPassSigmaIetaIetaMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------Seed : " << ((iPhoton->nPixelSeeds) > 0 ? "yes\n" : "no\n");
      cout << "--------------";
      cout << (susyCategory->getHasPixelSeed(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------Phi  : " << iPhoton->caloPosition.Phi() << endl;
      cout << "----------------------------";
      switch (susyCategory->getPhotonType(tag_, i)) {
      case FAIL:
	cout << "FAIL\n";
	break;
      case G:
	cout << "G\n";
	break;
      case E:
	cout << "E\n";
	break;
      case F:
	cout << "F\n";
	break;
      default:
	cout << "ERROR\n";
	break;
      }
      cout << "----------------------------";
      cout << (susyCategory->getIsDeciding(tag_, i) ? "DECIDING\n" : "NOT DECIDING\n");
      cout << endl;
    }
  }

  //print event quantities
  cout << "--------------Category     : ";
  switch(susyCategory->getEventCategory(tag_)) {
  case FAIL:
    cout << "FAIL\n";
    break;
  case GG:
    cout << "GG\n";
    break;
  case EG:
    cout << "EG\n";
    break;
  case EE:
    cout << "EE\n";
    break;
  case FF:
    cout << "FF\n";
    break;
  default:
    cout << "ERROR\n";
    break;
  }
  cout << "--------------Pass DPhiMin : ";
  cout << (susyCategory->getPassDPhiMin(tag_) ? "PASSED\n" : "FAILED\n") << endl;
  cout << "--------------DiEM ET      : " << susyCategory->getEvtDiEMET(tag_) << endl;
  cout << "--------------InvMass      : " << susyCategory->getEvtInvMass(tag_) << endl;
  cout << endl;
}
