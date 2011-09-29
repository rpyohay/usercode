#define GMSBAnalyzer_cxx
#include "GMSBAnalyzer.h"
#include "../../../GMSBTools/Filters/interface/Categorizer.h"
#include "../../../DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "GMSBTools/Filters/interface/Typedefs.h"
#include "SusyAnalysis/SusyNtuplizer/jec/JetMETObjects/interface/JetCorrectorParameters.h"
#include "SusyAnalysis/SusyNtuplizer/jec/JetMETObjects/interface/FactorizedJetCorrector.h"

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

void GMSBAnalyzer::runMETAnalysis(const std::string& outputFile)
{
   if (fChain == 0) return;

   //open file
   TFile out(outputFile.c_str(), "RECREATE");
   out.cd();

   //normalization region
   const int maxNormBin = 4;

   //invariant mass histograms (sanity check)
   TH1F mgg("mgg", "", 150, 0.0, 300.0);
   TH1F meg("meg", "", 150, 0.0, 300.0);
   TH1F mee("mee", "", 150, 0.0, 300.0);
   TH1F mff("mff", "", 150, 0.0, 300.0);
   mgg.GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
   meg.GetXaxis()->SetTitle("M_{e#gamma} (GeV)");
   mee.GetXaxis()->SetTitle("M_{ee} (GeV)");
   mff.GetXaxis()->SetTitle("M_{ff} (GeV)");

   //MET vs. di-EM ET histograms
   Double_t METBins[14] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 70.0, 100.0, 
			   150.0, 500.0};
   TH2F ggMETVsDiEMET("ggMETVsDiEMET", "", 13, METBins, 150, 0.0, 300.0);
   TH2F egMETVsDiEMET("egMETVsDiEMET", "", 13, METBins, 150, 0.0, 300.0);
   TH2F eeMETVsDiEMET("eeMETVsDiEMET", "", 13, METBins, 150, 0.0, 300.0);
   TH2F eeLowSidebandMETVsDiEMET("eeLowSidebandMETVsDiEMET", "", 13, METBins, 150, 0.0, 300.0);
   TH2F eeHighSidebandMETVsDiEMET("eeHighSidebandMETVsDiEMET", "", 13, METBins, 150, 0.0, 300.0);
   TH2F ffMETVsDiEMET("ffMETVsDiEMET", "", 13, METBins, 150, 0.0, 300.0);
   ggMETVsDiEMET.GetXaxis()->SetTitle("ME_{T} (GeV)");
   egMETVsDiEMET.GetXaxis()->SetTitle("ME_{T} (GeV)");
   eeMETVsDiEMET.GetXaxis()->SetTitle("ME_{T} (GeV)");
   eeLowSidebandMETVsDiEMET.GetXaxis()->SetTitle("ME_{T} (GeV)");
   eeHighSidebandMETVsDiEMET.GetXaxis()->SetTitle("ME_{T} (GeV)");
   ffMETVsDiEMET.GetXaxis()->SetTitle("ME_{T} (GeV)");
   ggMETVsDiEMET.GetYaxis()->SetTitle("Di-EM E_{T} (GeV)");
   egMETVsDiEMET.GetYaxis()->SetTitle("Di-EM E_{T} (GeV)");
   eeMETVsDiEMET.GetYaxis()->SetTitle("Di-EM E_{T} (GeV)");
   eeLowSidebandMETVsDiEMET.GetYaxis()->SetTitle("Di-EM E_{T} (GeV)");
   eeHighSidebandMETVsDiEMET.GetYaxis()->SetTitle("Di-EM E_{T} (GeV)");
   ffMETVsDiEMET.GetYaxis()->SetTitle("Di-EM E_{T} (GeV)");
   ggMETVsDiEMET.Sumw2();
   egMETVsDiEMET.Sumw2();
   eeMETVsDiEMET.Sumw2();
   eeLowSidebandMETVsDiEMET.Sumw2();
   eeHighSidebandMETVsDiEMET.Sumw2();
   ffMETVsDiEMET.Sumw2();

   //scaled ee and ff di-EM ET histograms
   TH1F eeDiEMETScaled("eeDiEMETScaled", "", 150, 0.0, 300.0);
   TH1F eeLowSidebandDiEMETScaled("eeLowSidebandDiEMETScaled", "", 150, 0.0, 300.0);
   TH1F eeHighSidebandDiEMETScaled("eeHighSidebandDiEMETScaled", "", 150, 0.0, 300.0);
   TH1F ffDiEMETScaled("ffDiEMETScaled", "", 150, 0.0, 300.0);
   eeDiEMETScaled.GetXaxis()->SetTitle("Di-EM E_{T} (GeV)");
   eeLowSidebandDiEMETScaled.GetXaxis()->SetTitle("Di-EM E_{T} (GeV)");
   eeHighSidebandDiEMETScaled.GetXaxis()->SetTitle("Di-EM E_{T} (GeV)");
   ffDiEMETScaled.GetXaxis()->SetTitle("Di-EM E_{T} (GeV)");

   //weights histograms
   TH1F eeWeights("eeWeights", "", 150, 0.0, 300.0);
   TH1F eeLowSidebandWeights("eeLowSidebandWeights", "", 150, 0.0, 300.0);
   TH1F eeHighSidebandWeights("eeHighSidebandWeights", "", 150, 0.0, 300.0);
   TH1F ffWeights("ffWeights", "", 150, 0.0, 300.0);
   eeWeights.GetXaxis()->SetTitle("Di-EM E_{T} (GeV)");
   eeLowSidebandWeights.GetXaxis()->SetTitle("Di-EM E_{T} (GeV)");
   eeHighSidebandWeights.GetXaxis()->SetTitle("Di-EM E_{T} (GeV)");
   ffWeights.GetXaxis()->SetTitle("Di-EM E_{T} (GeV)");

   //reweighted and normalized ee and ff MET histograms
   TH1F eeFinal("eeFinal", "", 13, METBins);
   TH1F eeLowSidebandFinal("eeLowSidebandFinal", "", 13, METBins);
   TH1F eeHighSidebandFinal("eeHighSidebandFinal", "", 13, METBins);
   TH1F ffFinal("ffFinal", "", 13, METBins);
   eeFinal.GetXaxis()->SetTitle("ME_{T} (GeV)");
   eeLowSidebandFinal.GetXaxis()->SetTitle("ME_{T} (GeV)");
   eeHighSidebandFinal.GetXaxis()->SetTitle("ME_{T} (GeV)");
   ffFinal.GetXaxis()->SetTitle("ME_{T} (GeV)");

   //canvas for the final plot
   TCanvas METCanvas("METCanvas", "", 600, 600);
   METCanvas.SetFillStyle(0);
   METCanvas.SetFillColor(0);
   METCanvas.SetGrid();
   TH1F* frame = METCanvas.DrawFrame(0.0, 0.0, 500.0, 10000.0);
   frame->GetXaxis()->SetTitle("ME_{T} (GeV)");

   //MET and di-EM ET vectors
   VFLOAT eeMETVec;
   VFLOAT eeLowSidebandMETVec;
   VFLOAT eeHighSidebandMETVec;
   VFLOAT ffMETVec;
   VFLOAT eeDiEMETVec;
   VFLOAT eeLowSidebandDiEMETVec;
   VFLOAT eeHighSidebandDiEMETVec;
   VFLOAT ffDiEMETVec;

   //set user-specified number of entries to process
   Long64_t nentries = fChain->GetEntriesFast();
   if (nEvts_ != -1) nentries = nEvts_;

   //loop over events
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (((jentry + 1) % 10000) == 0) cout << "Event " << (jentry + 1) << endl;

      //fill MET vs. di-EM ET and invariant mass histograms
      const double invMass = susyCategory->getEvtInvMass(tag_);
      double MET = -1.0;
      map<TString, susy::MET>::const_iterator iMET = susyEvent->metMap.find("pfMet");
      if (iMET != susyEvent->metMap.end()) MET = iMET->second.met();
      else cerr << "Error: PFMET not found in event " << (jentry + 1) << ".\n";
      const double diEMET = susyCategory->getEvtDiEMET(tag_);
      switch (susyCategory->getEventCategory(tag_)) {
      case GG:
	ggMETVsDiEMET.Fill(MET, diEMET);
	mgg.Fill(invMass);
	break;
      case EG:
	egMETVsDiEMET.Fill(MET, diEMET);
	meg.Fill(invMass);
	break;
      case EE:
	if ((invMass >= 71.0/*GeV*/) && (invMass < 76.0/*GeV*/)) {
	  eeLowSidebandMETVsDiEMET.Fill(MET, diEMET);
	  eeLowSidebandMETVec.push_back(MET);
	  eeLowSidebandDiEMETVec.push_back(diEMET);
	}
	if ((invMass >= 81.0/*GeV*/) && (invMass < 101.0/*GeV*/)) {
	  eeMETVsDiEMET.Fill(MET, diEMET);
	  eeMETVec.push_back(MET);
	  eeDiEMETVec.push_back(diEMET);
	}
	if ((invMass >= 106.0/*GeV*/) && (invMass < 111.0/*126.0*//*GeV*/)) {
	  eeHighSidebandMETVsDiEMET.Fill(MET, diEMET);
	  eeHighSidebandMETVec.push_back(MET);
	  eeHighSidebandDiEMETVec.push_back(diEMET);
	}
	mee.Fill(invMass);
	break;
      case FF:
	ffMETVsDiEMET.Fill(MET, diEMET);
	ffMETVec.push_back(MET);
	ffDiEMETVec.push_back(diEMET);
	mff.Fill(invMass);
	break;
      default:
	break;
      }
   }

   //fill weights histograms
   TH1D* ggMET = ggMETVsDiEMET.ProjectionX();
   TH1D* ffMET = ffMETVsDiEMET.ProjectionX();
   TH1D* ggDiEMET = ggMETVsDiEMET.ProjectionY();
   TH1D* eeDiEMET = eeMETVsDiEMET.ProjectionY();
   TH1D* eeLowSidebandDiEMET = eeLowSidebandMETVsDiEMET.ProjectionY();
   TH1D* eeHighSidebandDiEMET = eeHighSidebandMETVsDiEMET.ProjectionY();
   TH1D* ffDiEMET = ffMETVsDiEMET.ProjectionY();
   const float eeScale = ggMETVsDiEMET.GetEntries()/eeMETVsDiEMET.GetEntries();
   const float eeLowSidebandScale = 
     ggMETVsDiEMET.GetEntries()/eeLowSidebandMETVsDiEMET.GetEntries();
   const float eeHighSidebandScale = 
     ggMETVsDiEMET.GetEntries()/eeHighSidebandMETVsDiEMET.GetEntries();
   const float ffScale = ggMETVsDiEMET.GetEntries()/ffMETVsDiEMET.GetEntries();
   const float ffNorm = ggMET->Integral(1, maxNormBin)/ffMET->Integral(1, maxNormBin);
   eeDiEMET->Scale(eeScale);
   eeLowSidebandDiEMET->Scale(eeLowSidebandScale);
   eeHighSidebandDiEMET->Scale(eeHighSidebandScale);
   ffDiEMET->Scale(ffScale);
   eeDiEMETScaled.Add(eeDiEMET);
   eeLowSidebandDiEMETScaled.Add(eeLowSidebandDiEMET);
   eeHighSidebandDiEMETScaled.Add(eeHighSidebandDiEMET);
   ffDiEMETScaled.Add(ffDiEMET);
   eeWeights.Divide(ggDiEMET, &eeDiEMETScaled);
   eeLowSidebandWeights.Divide(ggDiEMET, &eeLowSidebandDiEMETScaled);
   eeHighSidebandWeights.Divide(ggDiEMET, &eeHighSidebandDiEMETScaled);
   ffWeights.Divide(ggDiEMET, &ffDiEMETScaled);

   //apply weights and normalization to ee and ff MET histograms
   for (VFLOAT_IT i = eeLowSidebandMETVec.begin(); i != eeLowSidebandMETVec.end(); ++i) {
     eeLowSidebandFinal.Fill(*i, eeLowSidebandWeights.GetBinContent((int)(eeLowSidebandDiEMETVec[i - eeLowSidebandMETVec.begin()]/eeLowSidebandWeights.GetBinWidth(1)) + 1));
   }
   for (VFLOAT_IT i = eeMETVec.begin(); i != eeMETVec.end(); ++i) {
     eeFinal.Fill(*i, eeWeights.GetBinContent((int)(eeDiEMETVec[i - eeMETVec.begin()]/
						    eeWeights.GetBinWidth(1)) + 1));
   }
   for (VFLOAT_IT i = eeHighSidebandMETVec.begin(); i != eeHighSidebandMETVec.end(); ++i) {
     eeHighSidebandFinal.Fill(*i, eeHighSidebandWeights.GetBinContent((int)(eeHighSidebandDiEMETVec[i - eeHighSidebandMETVec.begin()]/eeHighSidebandWeights.GetBinWidth(1)) + 1));
   }
   for (VFLOAT_IT i = ffMETVec.begin(); i != ffMETVec.end(); ++i) {
     ffFinal.Fill(*i, ffNorm*ffWeights.GetBinContent((int)(ffDiEMETVec[i - ffMETVec.begin()]/
							   ffWeights.GetBinWidth(1)) + 1));
   }
   eeLowSidebandFinal.Scale(2.0);
   eeHighSidebandFinal.Scale(2.0);
   eeFinal.Add(&eeLowSidebandFinal, -1.0);
   eeFinal.Add(&eeHighSidebandFinal, -1.0);
   const float eeNorm = ggMET->Integral(1, maxNormBin)/eeFinal.Integral(1, maxNormBin);
   eeFinal.Scale(eeNorm);

   //make final canvas
   METCanvas.cd();
   eeFinal.SetLineColor(9);
   ffFinal.SetLineColor(kMagenta);
   eeFinal.SetLineWidth(2);
   ffFinal.SetLineWidth(2);
   eeFinal.SetFillStyle(0);
   ffFinal.SetFillStyle(0);
   ggMET->SetFillStyle(0);
   eeFinal.Draw("HIST");
   ffFinal.Draw("HISTSAME");
   ggMET->Draw("SAME");

   //save
   mgg.Write();
   meg.Write();
   mee.Write();
   mff.Write();
   ggMETVsDiEMET.Write();
   egMETVsDiEMET.Write();
   eeMETVsDiEMET.Write();
   eeLowSidebandMETVsDiEMET.Write();
   eeHighSidebandMETVsDiEMET.Write();
   ffMETVsDiEMET.Write();
   eeDiEMETScaled.Write();
   eeLowSidebandDiEMETScaled.Write();
   eeHighSidebandDiEMETScaled.Write();
   ffDiEMETScaled.Write();
   eeWeights.Write();
   eeLowSidebandWeights.Write();
   eeHighSidebandWeights.Write();
   ffWeights.Write();
   eeFinal.Write();
   eeLowSidebandFinal.Write();
   eeHighSidebandFinal.Write();
   ffFinal.Write();
   METCanvas.Write();
   out.Write();
   out.Close();
}

void GMSBAnalyzer::runEEVsFFAnalysis(const std::string& outputFile)
{
   if (fChain == 0) return;

   //open file
   TFile out(outputFile.c_str(), "RECREATE");
   out.cd();

   //jet histograms filled once per jet per sample per reweighting category
   TH1F ETEENoReweighting("ETEENoReweighting", "", 50, 0.0, 500.0);
   TH1F ETFFNoReweighting("ETFFNoReweighting", "", 50, 0.0, 500.0);
   TH1F ETEEReweighting("ETEEReweighting", "", 50, 0.0, 500.0);
   TH1F ETFFReweighting("ETFFReweighting", "", 50, 0.0, 500.0);
   ETEENoReweighting.Sumw2();
   ETFFNoReweighting.Sumw2();
   ETEEReweighting.Sumw2();
   ETFFReweighting.Sumw2();
   ETEENoReweighting.GetXaxis()->SetTitle("E_{T} (GeV)");
   ETFFNoReweighting.GetXaxis()->SetTitle("E_{T} (GeV)");
   ETEEReweighting.GetXaxis()->SetTitle("E_{T} (GeV)");
   ETFFReweighting.GetXaxis()->SetTitle("E_{T} (GeV)");

   /*jet histograms filled once per jet per sample per reweighting category--1 for each bin of jet 
     ET (30-60, 60-90, >90 GeV)*/
   TH1F eta30GeVEENoReweighting("eta30GeVEENoReweighting", "", 50, -5.0, 5.0);
   TH1F eta60GeVEENoReweighting("eta60GeVEENoReweighting", "", 50, -5.0, 5.0);
   TH1F eta90GeVEENoReweighting("eta90GeVEENoReweighting", "", 50, -5.0, 5.0);
   TH1F eta30GeVFFNoReweighting("eta30GeVFFNoReweighting", "", 50, -5.0, 5.0);
   TH1F eta60GeVFFNoReweighting("eta60GeVFFNoReweighting", "", 50, -5.0, 5.0);
   TH1F eta90GeVFFNoReweighting("eta90GeVFFNoReweighting", "", 50, -5.0, 5.0);
   TH1F eta30GeVEEReweighting("eta30GeVEEReweighting", "", 50, -5.0, 5.0);
   TH1F eta60GeVEEReweighting("eta60GeVEEReweighting", "", 50, -5.0, 5.0);
   TH1F eta90GeVEEReweighting("eta90GeVEEReweighting", "", 50, -5.0, 5.0);
   TH1F eta30GeVFFReweighting("eta30GeVFFReweighting", "", 50, -5.0, 5.0);
   TH1F eta60GeVFFReweighting("eta60GeVFFReweighting", "", 50, -5.0, 5.0);
   TH1F eta90GeVFFReweighting("eta90GeVFFReweighting", "", 50, -5.0, 5.0);
   eta30GeVEENoReweighting.Sumw2();
   eta60GeVEENoReweighting.Sumw2();
   eta90GeVEENoReweighting.Sumw2();
   eta30GeVFFNoReweighting.Sumw2();
   eta60GeVFFNoReweighting.Sumw2();
   eta90GeVFFNoReweighting.Sumw2();
   eta30GeVEEReweighting.Sumw2();
   eta60GeVEEReweighting.Sumw2();
   eta90GeVEEReweighting.Sumw2();
   eta30GeVFFReweighting.Sumw2();
   eta60GeVFFReweighting.Sumw2();
   eta90GeVFFReweighting.Sumw2();
   eta30GeVEENoReweighting.GetXaxis()->SetTitle("#eta");
   eta60GeVEENoReweighting.GetXaxis()->SetTitle("#eta");
   eta90GeVEENoReweighting.GetXaxis()->SetTitle("#eta");
   eta30GeVFFNoReweighting.GetXaxis()->SetTitle("#eta");
   eta60GeVFFNoReweighting.GetXaxis()->SetTitle("#eta");
   eta90GeVFFNoReweighting.GetXaxis()->SetTitle("#eta");
   eta30GeVEEReweighting.GetXaxis()->SetTitle("#eta");
   eta60GeVEEReweighting.GetXaxis()->SetTitle("#eta");
   eta90GeVEEReweighting.GetXaxis()->SetTitle("#eta");
   eta30GeVFFReweighting.GetXaxis()->SetTitle("#eta");
   eta60GeVFFReweighting.GetXaxis()->SetTitle("#eta");
   eta90GeVFFReweighting.GetXaxis()->SetTitle("#eta");

   /*jet histograms filled once per event per sample per reweighting category--1 for each bin of 
     jet ET (30-60, 60-90, >90 GeV)*/
   TH1F n30GeVEENoReweighting("n30GeVEENoReweighting", "", 10, -0.5, 9.5);
   TH1F n60GeVEENoReweighting("n60GeVEENoReweighting", "", 10, -0.5, 9.5);
   TH1F n90GeVEENoReweighting("n90GeVEENoReweighting", "", 10, -0.5, 9.5);
   TH1F n30GeVFFNoReweighting("n30GeVFFNoReweighting", "", 10, -0.5, 9.5);
   TH1F n60GeVFFNoReweighting("n60GeVFFNoReweighting", "", 10, -0.5, 9.5);
   TH1F n90GeVFFNoReweighting("n90GeVFFNoReweighting", "", 10, -0.5, 9.5);
   TH1F n30GeVEEReweighting("n30GeVEEReweighting", "", 10, -0.5, 9.5);
   TH1F n60GeVEEReweighting("n60GeVEEReweighting", "", 10, -0.5, 9.5);
   TH1F n90GeVEEReweighting("n90GeVEEReweighting", "", 10, -0.5, 9.5);
   TH1F n30GeVFFReweighting("n30GeVFFReweighting", "", 10, -0.5, 9.5);
   TH1F n60GeVFFReweighting("n60GeVFFReweighting", "", 10, -0.5, 9.5);
   TH1F n90GeVFFReweighting("n90GeVFFReweighting", "", 10, -0.5, 9.5);
   n30GeVEENoReweighting.Sumw2();
   n60GeVEENoReweighting.Sumw2();
   n90GeVEENoReweighting.Sumw2();
   n30GeVFFNoReweighting.Sumw2();
   n60GeVFFNoReweighting.Sumw2();
   n90GeVFFNoReweighting.Sumw2();
   n30GeVEEReweighting.Sumw2();
   n60GeVEEReweighting.Sumw2();
   n90GeVEEReweighting.Sumw2();
   n30GeVFFReweighting.Sumw2();
   n60GeVFFReweighting.Sumw2();
   n90GeVFFReweighting.Sumw2();
   n30GeVEENoReweighting.GetXaxis()->SetTitle("N_{j}");
   n60GeVEENoReweighting.GetXaxis()->SetTitle("N_{j}");
   n90GeVEENoReweighting.GetXaxis()->SetTitle("N_{j}");
   n30GeVFFNoReweighting.GetXaxis()->SetTitle("N_{j}");
   n60GeVFFNoReweighting.GetXaxis()->SetTitle("N_{j}");
   n90GeVFFNoReweighting.GetXaxis()->SetTitle("N_{j}");
   n30GeVEEReweighting.GetXaxis()->SetTitle("N_{j}");
   n60GeVEEReweighting.GetXaxis()->SetTitle("N_{j}");
   n90GeVEEReweighting.GetXaxis()->SetTitle("N_{j}");
   n30GeVFFReweighting.GetXaxis()->SetTitle("N_{j}");
   n60GeVFFReweighting.GetXaxis()->SetTitle("N_{j}");
   n90GeVFFReweighting.GetXaxis()->SetTitle("N_{j}");
   // vector<vector<vector<TH1F> > > dPhiEM1ToLeading;
   // vector<vector<vector<TH1F> > > dPhiEM2ToLeading;
   // vector<vector<vector<TH1F> > > dPhiDiEMToLeading;
   // vector<vector<vector<TH1F> > > dPhiEM1ToTrailing;
   // vector<vector<vector<TH1F> > > dPhiEM2ToTrailing;
   // vector<vector<vector<TH1F> > > dPhiDiEMToTrailing;
   // vector<vector<vector<TH1F> > > dPhiLeadingToMET;
   // vector<vector<vector<TH1F> > > dPhiTrailingToMET;

   // //MET histograms filled once per event per sample per reweighting category
   // vector<vector<vector<TH1F> > > dPhiEM1ToMET;
   // vector<vector<vector<TH1F> > > dPhiEM2ToMET;
   // vector<vector<vector<TH1F> > > dPhiDiEMToMET;

   // //HT histograms filled once per event per sample per reweighting category
   // vector<vector<vector<TH1F> > > HT;

   //set up on-the-fly jet corrections for PF jets
   vector<JetCorrectorParameters> PFJECs;
   PFJECs.push_back(JetCorrectorParameters("/Users/rachelyohay/RA3/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L1FastJet.txt"));
   PFJECs.push_back(JetCorrectorParameters("/Users/rachelyohay/RA3/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L2Relative.txt"));
   PFJECs.push_back(JetCorrectorParameters("/Users/rachelyohay/RA3/src/SusyAnalysis/SusyNtuplizer/jec/Jec11_V1_AK5PF_L3Absolute.txt"));
   FactorizedJetCorrector PFJetCorrector(PFJECs);
   vector<float> corrDiff;
   vector<float> corrSame;
   unsigned nEvtsWithJetCorrDiff = 0;

   //set user-specified number of entries to process
   Long64_t nentries = fChain->GetEntriesFast();
   if (nEvts_ != -1) nentries = nEvts_;

   //loop over events
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) {
       cout << "No. events processed: " << jentry << endl;
       break;
     }
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if (((jentry + 1) % 10000) == 0) cout << "Event " << (jentry + 1) << endl;

     //proceed if this is an ee or ff event and jets are found
     const int category = susyCategory->getEventCategory(tag_);
     const double invMass = susyCategory->getEvtInvMass(tag_);
     if (((category == EE) && (invMass >= 81.0/*GeV*/) && (invMass < 101.0/*GeV*/)) || 
	 (category == FF)) {
       map<TString, susy::PFJetCollection>::const_iterator iJets = susyEvent->pfJets.find("ak5");
       if (iJets != susyEvent->pfJets.end()) {

	 //count jets with particular ET
	 unsigned int nJets30EE = 0;
	 unsigned int nJets60EE = 0;
	 unsigned int nJets90EE = 0;
	 unsigned int nJets30FF = 0;
	 unsigned int nJets60FF = 0;
	 unsigned int nJets90FF = 0;

	 //loop over PF jets
	 bool noDiffYet = true;
	 for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
	      iJet != iJets->second.end(); ++iJet) {

	   //compute corrected P4
	   PFJetCorrector.setJetEta(iJet->momentum.Eta());
	   PFJetCorrector.setJetPt(iJet->momentum.Pt());
	   PFJetCorrector.setJetA(iJet->jetArea);
	   PFJetCorrector.setRho(susyEvent->rho);
	   double onTheFlyCorr = PFJetCorrector.getCorrection();
	   float storedCorr = -1.0;
	   map<TString, Float_t>::const_iterator iCorr = iJet->jecScaleFactors.find("L1FastL2L3");
	   if (iCorr != iJet->jecScaleFactors.end()) storedCorr = iCorr->second;
	   if (storedCorr != onTheFlyCorr) {
	     // cout << "PF jet stored correction: " << storedCorr << ", on the fly correction: ";
	     // cout << onTheFlyCorr << endl;
	     corrDiff.push_back((storedCorr - onTheFlyCorr)/storedCorr);
	     if (noDiffYet) {
	       ++nEvtsWithJetCorrDiff;
	       noDiffYet = false;
	     }
	   }
	   else corrSame.push_back(0.0);
	   TLorentzVector corrP4 = storedCorr*iJet->momentum;

	   //calculate the jet-EM overlap
	   map<TString, susy::PhotonCollection>::const_iterator iPhotonMap = 
	     susyEvent->photons.find((const char*)tag_);
	   if (iPhotonMap != susyEvent->photons.end()) {
	     float EM1Eta = -9.0;
	     float EM1Phi = -9.0;
	     float EM2Eta = -9.0;
	     float EM2Phi = -9.0;
	     unsigned int count = 0;
	     for (susy::PhotonCollection::const_iterator iPhoton = iPhotonMap->second.begin(); 
		  iPhoton != iPhotonMap->second.end(); ++iPhoton) {
	       if (susyCategory->getIsDeciding(tag_, iPhoton - iPhotonMap->second.begin())) {
		 if ((EM1Eta != -9.0) && (EM1Phi != -9.0)) {
		   EM2Eta = iPhoton->caloPosition.Eta();
		   EM2Phi = iPhoton->caloPosition.Phi();
		   ++count;
		 }
		 else {
		   EM1Eta = iPhoton->caloPosition.Eta();
		   EM1Phi = iPhoton->caloPosition.Phi();
		   ++count;
		 }
	       }
	     }//loop over photons
	     if (count != 2) {
	       cerr << "Error: found " << count << " deciding EM objects in event ";
	       cerr << (jentry + 1) << ".\n";
	     }

	     /*only consider jets...
	       ...that do not overlap with the 2 selected EM objects
	       ...in |eta| < 5.0
	       ...passing PF jet ID*/
	     if ((deltaR(corrP4.Eta(), corrP4.Phi(), EM1Eta, EM1Phi) > 0.8) && 
		 (deltaR(corrP4.Eta(), corrP4.Phi(), EM2Eta, EM2Phi) > 0.8) && 
		 (fabs(corrP4.Eta()) < 5.0)) {
	       bool passed = false;
	       if ((iJet->neutralHadronEnergy/iJet->momentum.Energy() < 0.99) && 
		   (iJet->neutralEmEnergy/iJet->momentum.Energy() < 0.99) && 
		   ((unsigned int)iJet->nConstituents > 1)) {
	       	 if (fabs(iJet->momentum.Eta()) < 2.4) {
	       	   if ((iJet->chargedHadronEnergy > 0.0) && 
	       	       ((int)iJet->chargedMultiplicity > 0) && 
	       	       (iJet->chargedEmEnergy/iJet->momentum.Energy() < 0.99)) passed = true;
	       	 }
	       	 else passed = true;
	       }
	       if (passed) {

	       	 //fill eta and ET histograms
	       	 if ((corrP4.Et() >= 30.0/*GeV*/) && (corrP4.Et() < 60.0/*GeV*/)) {
		   if (category == EE) {
		     eta30GeVEENoReweighting.Fill(corrP4.Eta());
		     ++nJets30EE;
		   }
		   if (category == FF) {
		     eta30GeVFFNoReweighting.Fill(corrP4.Eta());
		     ++nJets30FF;
		   }
	       	 }
	       	 if ((corrP4.Et() >= 60.0/*GeV*/) && (corrP4.Et() < 90.0/*GeV*/)) {
		   if (category == EE) {
		     eta60GeVEENoReweighting.Fill(corrP4.Eta());
		     ++nJets60EE;
		   }
		   if (category == FF) {
		     eta60GeVFFNoReweighting.Fill(corrP4.Eta());
		     ++nJets60FF;
		   }
	       	 }
	       	 if (corrP4.Et() >= 90.0/*GeV*/) {
		   if (category == EE) {
		     eta90GeVEENoReweighting.Fill(corrP4.Eta());
		     ++nJets90EE;
		   }
		   if (category == FF) {
		     eta90GeVFFNoReweighting.Fill(corrP4.Eta());
		     ++nJets90FF;
		   }
	       	 }
		 if (category == EE) ETEENoReweighting.Fill(corrP4.Et());
		 if (category == FF) ETFFNoReweighting.Fill(corrP4.Et());
	       }//if passed jet ID
	     }//if (deltaR...
	   }//if photons are found
	   else cerr << "Error: photons not found in event " << (jentry + 1) << ".\n";
	 }//loop over jets

	 //fill n histogram
	 if (category == EE) {
	   n30GeVEENoReweighting.Fill(nJets30EE);
	   n60GeVEENoReweighting.Fill(nJets60EE);
	   n90GeVEENoReweighting.Fill(nJets90EE);
	 }
	 if (category == FF) {
	   n30GeVFFNoReweighting.Fill(nJets30FF);
	   n60GeVFFNoReweighting.Fill(nJets60FF);
	   n90GeVFFNoReweighting.Fill(nJets90FF);
	 }
       }//if jets found
       else cerr << "Error: PF jets not found in event " << (jentry + 1) << ".\n";
     }//if ee or ff
   }//loop over events

   //print jet corrections information
   sort(corrDiff.begin(), corrDiff.end());
   cout << "No. of jets in which stored correction differs from on-the-fly correction: ";
   cout << corrDiff.size() << endl;
   cout << "No. of events with at least 1 jet correction disagreement: " << nEvtsWithJetCorrDiff;
   cout << endl;
   cout << "Minimum percentage difference from stored correction: " << *(corrDiff.begin()) << endl;
   cout << "Maximum percentage difference from stored correction: " << *(corrDiff.end() - 1);
   cout << endl;
   TH1F corrDiffHist("corrDiffHist", "", 20, -1.0, 1.0);
   for (VFLOAT_IT i = corrDiff.begin(); i != corrDiff.end(); ++i) { corrDiffHist.Fill(*i); }
   cout << corrSame.size() << endl;
   for (VFLOAT_IT i = corrSame.begin(); i != corrSame.end(); ++i) { corrDiffHist.Fill(*i); }

   //write histograms
   out.Write();
   out.Close();
}

void GMSBAnalyzer::runCutFlowAnalysis(/*const std::string& outputFile*/)
{
   if (fChain == 0) return;

//    //open file
//    TFile out(outputFile.c_str(), "RECREATE");
//    out.cd();

   //set user-specified number of entries to process
   Long64_t nentries = fChain->GetEntriesFast();
   if (nEvts_ != -1) nentries = nEvts_;

   //counters
   unsigned int nPassingET = 0;
   unsigned int nPassingPreselection = 0;
   unsigned int nPassingG = 0;
   unsigned int nPassingE = 0;
   unsigned int nPassingF = 0;
   unsigned int nPassingGG = 0;
   unsigned int nPassingGE = 0;
   unsigned int nPassingEE = 0;
   unsigned int nPassingFF = 0;
   unsigned int nPassingGGE = 0;
   unsigned int nPassingGEE = 0;
   unsigned int nPassingGGMu = 0;
   unsigned int nPassingGEMu = 0;
   unsigned int nPassingEEMu = 0;
   unsigned int nPassingFFMu = 0;

   //loop over events
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     susyEvent->Init();
     susyCategory->reset();
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if (((jentry + 1) % /*10000*/1) == 0) {
       cout << "Event " << (jentry + 1) << endl;
       cout << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
       cout << susyEvent->luminosityBlockNumber << endl;
     }

     //loop over photons
     unsigned int nPhotonsPassingET1 = 0;
     unsigned int nPhotonsPassingET2 = 0;
     unsigned int nPhotonsPassingPreselection1 = 0;
     unsigned int nPhotonsPassingPreselection2 = 0;
     unsigned int nPhotonsPassingG1 = 0;
     unsigned int nPhotonsPassingG2 = 0;
     unsigned int nPhotonsPassingE1 = 0;
     unsigned int nPhotonsPassingE2 = 0;
     unsigned int nPhotonsPassingF1 = 0;
     unsigned int nPhotonsPassingF2 = 0;
     const unsigned int nPhotons = 
       susyCategory->getPassETMin1()->find((const char*)tag_)->second.size();
     cout << "nPhotons = " << nPhotons << endl;
     cout << "# photons = " << susyEvent->photons.find((const char*)tag_)->second.size() << endl;
     if ((nPhotons != susyCategory->getPassETMin2()->find((const char*)tag_)->second.size()) || 
	 (nPhotons != 
	  susyCategory->getPassECALIsoMax()->find((const char*)tag_)->second.size()) || 
	 (nPhotons != 
	  susyCategory->getPassHCALIsoMax()->find((const char*)tag_)->second.size()) || 
	 (nPhotons != susyCategory->getPassHOverEMax()->find((const char*)tag_)->second.size()) || 
	 (nPhotons != susyCategory->getPassR9Max()->find((const char*)tag_)->second.size())) {
       cerr << "Error: category information not filled correctly in event " << (jentry + 1);
       cerr << ".\n";
     }
     for (unsigned int iPhoton = 0; iPhoton < nPhotons; ++iPhoton) {

       //determine if photon passed ECAL isolation, HCAL isolation, H/E, and R9
       bool passPreselection = (susyCategory->getPassECALIsoMax(tag_, iPhoton) && 
				susyCategory->getPassHCALIsoMax(tag_, iPhoton) && 
				susyCategory->getPassHOverEMax(tag_, iPhoton) && 
				susyCategory->getPassR9Max(tag_, iPhoton));

       //determine photon type
       int type = susyCategory->getPhotonType(tag_, iPhoton);

       //count how many passing ET
       cout << "Photon " << iPhoton << " ET: ";
       cout << susyEvent->photons.find((const char*)tag_)->second[iPhoton].momentum.Et() << endl;
       if (susyCategory->getPassETMin1(tag_, iPhoton)) {
	 cout << "Photon " << iPhoton << " passed ET1\n";
	 ++nPhotonsPassingET1;

	 //count how many passing ECAL isolation, HCAL isolation, H/E, and R9
	 if (passPreselection) {
	   ++nPhotonsPassingPreselection1;

	   //count how many e, g, and f
	   if (type == G) ++nPhotonsPassingG1;
	   if (type == E) ++nPhotonsPassingE1;
	   if (type == F) ++nPhotonsPassingF1;
	 }
       }
       else if (susyCategory->getPassETMin2(tag_, iPhoton)) {
	 cout << "Photon " << iPhoton << " passed ET2\n";
	 ++nPhotonsPassingET2;

	 //count how many passing ECAL isolation, HCAL isolation, H/E, and R9
	 if (passPreselection) {
	   ++nPhotonsPassingPreselection2;

	   //count how many e, g, and f
	   if (type == G) ++nPhotonsPassingG2;
	   if (type == E) ++nPhotonsPassingE2;
	   if (type == F) ++nPhotonsPassingF2;
	 }
       }
     }

     //loop over muons
     unsigned int nMu = 0;
     for (vector<susy::Muon>::const_iterator iMuon = susyEvent->muons.begin(); 
	  iMuon != susyEvent->muons.end(); ++iMuon) {

       //count muons passing loose selection cuts
       if (iMuon->isGlobalMuon() && iMuon->isTrackerMuon() && (iMuon->nValidMuonHits >= 1) && 
	   (iMuon->nMatches >= 2) && (iMuon->nValidTrackerHits > 10) && 
	   (susyEvent->tracks[(unsigned int)iMuon->combinedTrackIndex].numberOfValidPixelHits >= 
	    1) && 
	   (susyEvent->tracks[(unsigned int)iMuon->combinedTrackIndex].normChi2() < 10) && 
	   (fabs(susyEvent->tracks[(unsigned int)iMuon->combinedTrackIndex].dxy()) < 0.2)) ++nMu;
     }

     //count how many passing ET
     cout << "nPhotonsPassingET1 = " << nPhotonsPassingET1 << endl;
     cout << "nPhotonsPassingET2 = " << nPhotonsPassingET2 << endl;
     if ((nPhotonsPassingET1 >= 2) || ((nPhotonsPassingET1 >= 1) && (nPhotonsPassingET2 >= 1))) {
       ++nPassingET;

       //count how many passing ECAL isolation, HCAL isolation, H/E, and R9
       if ((nPhotonsPassingPreselection1 >= 2) || 
	   ((nPhotonsPassingPreselection1 >= 1) && 
	    (nPhotonsPassingPreselection2 >= 1))) {
	 ++nPassingPreselection;

	 //count how many g
	 if (nPhotonsPassingG1 >= 1) ++nPassingG;

	 //count how many e
	 if (nPhotonsPassingE1 >= 1) ++nPassingE;

	 //count how many f
	 if (nPhotonsPassingF1 >= 1) ++nPassingF;

	 //count how many gg
	 if ((nPhotonsPassingG1 >= 2) || 
	     ((nPhotonsPassingG1 >= 1) && (nPhotonsPassingG2 >= 1))) {
	   ++nPassingGG;

	   //count how many passing ggmu
	   if (nMu >= 1) ++nPassingGGMu;
	 }

	 //count how many ge
	 if (((nPhotonsPassingG1 >= 1) && (nPhotonsPassingE1 >= 1)) || 
	     ((nPhotonsPassingG1 >= 1) && (nPhotonsPassingE2 >= 1)) || 
	     ((nPhotonsPassingG2 >= 1) && (nPhotonsPassingE1 >= 1))) {
	   ++nPassingGE;

	   //count how many passing gemu
	   if (nMu >= 1) ++nPassingGEMu;
	 }

	 //count how many ee
	 if ((nPhotonsPassingE1 >= 1) && (nPhotonsPassingE2 >= 1)) {
	   ++nPassingEE;

	   //count how many passing eemu
	   if (nMu >= 1) ++nPassingEEMu;
	 }

	 //count how many ff
	 if ((nPhotonsPassingF1 >= 1) && (nPhotonsPassingF2 >= 1)) {
	   ++nPassingFF;

	   //count how many passing ffmu
	   if (nMu >= 1) ++nPassingFFMu;
	 }

	 //count how many gge
	 if (((nPhotonsPassingG1 >= 2) && (nPhotonsPassingE1 >= 1)) || 
	     ((nPhotonsPassingG1 >= 2) && (nPhotonsPassingE2 >= 1)) || 
	     ((nPhotonsPassingG2 >= 2) && (nPhotonsPassingE1 >= 1)) || 
	     ((nPhotonsPassingG1 >= 1) && (nPhotonsPassingG2 >= 1) && (nPhotonsPassingE1 >= 1)) || 
	     ((nPhotonsPassingG1 >= 1) && (nPhotonsPassingG2 >= 1) && 
	      (nPhotonsPassingE2 >= 1))) ++nPassingGGE;

	 //count how many gee
	 if (((nPhotonsPassingE1 >= 2) && (nPhotonsPassingG1 >= 1)) || 
	     ((nPhotonsPassingE1 >= 2) && (nPhotonsPassingG2 >= 1)) || 
	     ((nPhotonsPassingE2 >= 2) && (nPhotonsPassingG1 >= 1)) || 
	     ((nPhotonsPassingE1 >= 1) && (nPhotonsPassingE2 >= 1) && (nPhotonsPassingG1 >= 1)) || 
	     ((nPhotonsPassingE1 >= 1) && (nPhotonsPassingE2 >= 1) && 
	      (nPhotonsPassingG2 >= 1))) ++nPassingGEE;
       }
     }
   }

   //print counts
   cout << "No. passing ET: " << nPassingET << endl;
   cout << "No. passing ECAL isolation, HCAL isolation, H/E, and R9: " << nPassingPreselection;
   cout << endl;
   cout << "No. passing g: " << nPassingG << endl;
   cout << "No. passing e: " << nPassingE << endl;
   cout << "No. passing f: " << nPassingF << endl;
   cout << "No. passing gg: " << nPassingGG << endl;
   cout << "No. passing ge: " << nPassingGE << endl;
   cout << "No. passing ee: " << nPassingEE << endl;
   cout << "No. passing ff: " << nPassingFF << endl;
   cout << "No. passing gge: " << nPassingGGE << endl;
   cout << "No. passing gee: " << nPassingGEE << endl;
   cout << "No. passing ggmu: " << nPassingGGMu << endl;
   cout << "No. passing gemu: " << nPassingGEMu << endl;
   cout << "No. passing eemu: " << nPassingEEMu << endl;
   cout << "No. passing ffmu: " << nPassingFFMu << endl;
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
   // ofstream outTxt((outputFile.replace(outputFile.find("root"), 4, "txt")).c_str());

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
	if (((susyEvent->runNumber == 166841) && (susyEvent->eventNumber == 483584716)) || 
	    ((susyEvent->runNumber == 166864) && (susyEvent->eventNumber == 394801663)) || 
	    ((susyEvent->runNumber == 166864) && (susyEvent->eventNumber == 416227374)) || 
	    ((susyEvent->runNumber == 163817) && (susyEvent->eventNumber == 54364321)) || 
	    ((susyEvent->runNumber == 163589) && (susyEvent->eventNumber == 54614720)) || 
	    ((susyEvent->runNumber == 163596) && (susyEvent->eventNumber == 3521103)) || 
	    ((susyEvent->runNumber == 163660) && (susyEvent->eventNumber == 18512069)) || 
	    ((susyEvent->runNumber == 163664) && (susyEvent->eventNumber == 11099963)) || 
	    ((susyEvent->runNumber == 163758) && (susyEvent->eventNumber == 340434605)) || 
	    ((susyEvent->runNumber == 163817) && (susyEvent->eventNumber == 39521312)) || 
	    ((susyEvent->runNumber == 165567) && (susyEvent->eventNumber == 660879884)) || 
	    ((susyEvent->runNumber == 165993) && (susyEvent->eventNumber == 78420443)) || 
	    ((susyEvent->runNumber == 166380) && (susyEvent->eventNumber == 380536201)) || 
	    ((susyEvent->runNumber == 166380) && (susyEvent->eventNumber == 1455639992)) || 
	    ((susyEvent->runNumber == 166408) && (susyEvent->eventNumber == 645967293)) || 
	    ((susyEvent->runNumber == 166554) && (susyEvent->eventNumber == 681773881)) || 
	    ((susyEvent->runNumber == 166784) && (susyEvent->eventNumber == 136897277)) || 
	    ((susyEvent->runNumber == 166787) && (susyEvent->eventNumber == 9069334)) || 
	    ((susyEvent->runNumber == 166565) && (susyEvent->eventNumber == 127439496)) || 
	    ((susyEvent->runNumber == 166701) && (susyEvent->eventNumber == 69918333)) || 
	    ((susyEvent->runNumber == 167039) && (susyEvent->eventNumber == 146433156)) || 
	    ((susyEvent->runNumber == 167041) && (susyEvent->eventNumber == 120702599)) || 
	    ((susyEvent->runNumber == 166946) && (susyEvent->eventNumber == 200822952)) || 
	    ((susyEvent->runNumber == 166950) && (susyEvent->eventNumber == 587885196)) || 
	    ((susyEvent->runNumber == 167284) && (susyEvent->eventNumber == 429291704)) || 
	    ((susyEvent->runNumber == 167282) && (susyEvent->eventNumber == 22617452)) || 
	    ((susyEvent->runNumber == 167898) && (susyEvent->eventNumber == 485562364)) || 
	    ((susyEvent->runNumber == 167675) && (susyEvent->eventNumber == 622183354)) || 
	    ((susyEvent->runNumber == 167675) && (susyEvent->eventNumber == 646543455)) || 
	    ((susyEvent->runNumber == 167746) && (susyEvent->eventNumber == 90136125)) || 
	    ((susyEvent->runNumber == 167786) && (susyEvent->eventNumber == 45391607)) || 
	    ((susyEvent->runNumber == 167786) && (susyEvent->eventNumber == 50156901)) || 
	    ((susyEvent->runNumber == 167786) && (susyEvent->eventNumber == 54739364))) {
	  cout << susyEvent->runNumber << " " << susyEvent->eventNumber << " " << invMass;
	  cout << " Dave yes Rachel no\n";
	}
	if (((susyEvent->runNumber == 166841) && (susyEvent->eventNumber == 279690939)) || 
	    ((susyEvent->runNumber == 166841) && (susyEvent->eventNumber == 501897075)) || 
	    ((susyEvent->runNumber == 166841) && (susyEvent->eventNumber == 779358311)) || 
	    ((susyEvent->runNumber == 166859) && (susyEvent->eventNumber == 131005283)) || 
	    ((susyEvent->runNumber == 166859) && (susyEvent->eventNumber == 150130237)) || 
	    ((susyEvent->runNumber == 166889) && (susyEvent->eventNumber == 249801723)) || 
	    ((susyEvent->runNumber == 163252) && (susyEvent->eventNumber == 70320781)) || 
	    ((susyEvent->runNumber == 163255) && (susyEvent->eventNumber == 548957669)) || 
	    ((susyEvent->runNumber == 163270) && (susyEvent->eventNumber == 523122198)) || 
	    ((susyEvent->runNumber == 163332) && (susyEvent->eventNumber == 293130827)) || 
	    ((susyEvent->runNumber == 163332) && (susyEvent->eventNumber == 359113407)) || 
	    ((susyEvent->runNumber == 163738) && (susyEvent->eventNumber == 172527800)) || 
	    ((susyEvent->runNumber == 163757) && (susyEvent->eventNumber == 4661645)) || 
	    ((susyEvent->runNumber == 163758) && (susyEvent->eventNumber == 405869965)) || 
	    ((susyEvent->runNumber == 163237) && (susyEvent->eventNumber == 85112659)) || 
	    ((susyEvent->runNumber == 165467) && (susyEvent->eventNumber == 262854623)) || 
	    ((susyEvent->runNumber == 165472) && (susyEvent->eventNumber == 244217004)) || 
	    ((susyEvent->runNumber == 165472) && (susyEvent->eventNumber == 845896439)) || 
	    ((susyEvent->runNumber == 165514) && (susyEvent->eventNumber == 586297947)) || 
	    ((susyEvent->runNumber == 165558) && (susyEvent->eventNumber == 9393965)) || 
	    ((susyEvent->runNumber == 165567) && (susyEvent->eventNumber == 69304564)) || 
	    ((susyEvent->runNumber == 165993) && (susyEvent->eventNumber == 1133251955)) || 
	    ((susyEvent->runNumber == 166049) && (susyEvent->eventNumber == 463355641)) || 
	    ((susyEvent->runNumber == 166049) && (susyEvent->eventNumber == 796442372)) || 
	    ((susyEvent->runNumber == 166346) && (susyEvent->eventNumber == 109458583)) || 
	    ((susyEvent->runNumber == 165364) && (susyEvent->eventNumber == 857448405)) || 
	    ((susyEvent->runNumber == 166380) && (susyEvent->eventNumber == 590691805)) || 
	    ((susyEvent->runNumber == 166408) && (susyEvent->eventNumber == 430241492)) || 
	    ((susyEvent->runNumber == 166408) && (susyEvent->eventNumber == 632010396)) || 
	    ((susyEvent->runNumber == 166512) && (susyEvent->eventNumber == 1027686669)) || 
	    ((susyEvent->runNumber == 166514) && (susyEvent->eventNumber == 22711260)) || 
	    ((susyEvent->runNumber == 166701) && (susyEvent->eventNumber == 571056335)) || 
	    ((susyEvent->runNumber == 167041) && (susyEvent->eventNumber == 244226879)) || 
	    ((susyEvent->runNumber == 167043) && (susyEvent->eventNumber == 105723195)) || 
	    ((susyEvent->runNumber == 166950) && (susyEvent->eventNumber == 1183700460)) || 
	    ((susyEvent->runNumber == 166950) && (susyEvent->eventNumber == 1410187165)) || 
	    ((susyEvent->runNumber == 166966) && (susyEvent->eventNumber == 127715395)) || 
	    ((susyEvent->runNumber == 167103) && (susyEvent->eventNumber == 79983187)) || 
	    ((susyEvent->runNumber == 167284) && (susyEvent->eventNumber == 398037294)) || 
	    ((susyEvent->runNumber == 167284) && (susyEvent->eventNumber == 972681846)) || 
	    ((susyEvent->runNumber == 167807) && (susyEvent->eventNumber == 1387950147)) || 
	    ((susyEvent->runNumber == 167830) && (susyEvent->eventNumber == 895644093)) || 
	    ((susyEvent->runNumber == 167830) && (susyEvent->eventNumber == 1232220016)) || 
	    ((susyEvent->runNumber == 167898) && (susyEvent->eventNumber == 551478829)) || 
	    ((susyEvent->runNumber == 167913) && (susyEvent->eventNumber == 329546015)) || 
	    ((susyEvent->runNumber == 167674) && (susyEvent->eventNumber == 182429068)) || 
	    ((susyEvent->runNumber == 167674) && (susyEvent->eventNumber == 188130251)) || 
	    ((susyEvent->runNumber == 167674) && (susyEvent->eventNumber == 212922847)) || 
	    ((susyEvent->runNumber == 167675) && (susyEvent->eventNumber == 881773084))) {
	  cout << susyEvent->runNumber << " " << susyEvent->eventNumber << " " << invMass;
	  cout << " Rachel yes Dave no\n";
	}
	// outTxt << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
	// if ((invMass >= 0.0) && (invMass < 71.0)) outTxt << "1\n";
	// if ((invMass >= 71.0) && (invMass < 76.0)) outTxt << "2\n";
	// if ((invMass >= 76.0) && (invMass < 81.0)) outTxt << "3\n";
	// if ((invMass >= 81.0) && (invMass < 101.0)) outTxt << "4\n";
	// if ((invMass >= 101.0) && (invMass < 106.0)) outTxt << "5\n";
	// if ((invMass >= 106.0) && (invMass < 111.0)) outTxt << "6\n";
	// if (invMass >= 111.0) outTxt << "7\n";
      }
   }

   //write histograms
   // outTxt.close();
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
