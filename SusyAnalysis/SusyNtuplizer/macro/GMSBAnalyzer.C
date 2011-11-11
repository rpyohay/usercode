#define GMSBAnalyzer_cxx
#include "GMSBAnalyzer.h"
#include "../../../GMSBTools/Filters/interface/Categorizer.h"
#include "../../../DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <fstream>
#include <numeric>
#include <TH2.h>
#include "TH3F.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TRandom3.h"
#include "TRegexp.h"
#include "TH3F.h"
#include "RooRealVar.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "PhysicsTools/TagAndProbe/interface/RooCMSShape.h"
#include "PhysicsTools/TagAndProbe/interface/RooErfXGaussian.h"
#include "GMSBTools/Filters/interface/Typedefs.h"
#include "SusyAnalysis/SusyNtuplizer/jec/JetMETObjects/interface/JetCorrectorParameters.h"
#include "SusyAnalysis/SusyNtuplizer/jec/JetMETObjects/interface/FactorizedJetCorrector.h"

using namespace RooFit;

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

// void GMSBAnalyzer::METCorrection(vector<susy::PFJet*>& jets, susy::MET& MET, 
// 				 const Long64_t evt) const
// {
//   //TODO: use L1FastL2L3(Residual) - L1Fast for MC(data)
//   //TODO: implement eta < 4.7 procedure
//   for (vector<susy::PFJet*>::iterator iJet = jets.begin(); iJet != jets.end(); ++iJet) {
//     susy::PFJet* jet = *iJet;
//     map<TString, Float_t>::const_iterator JEC = jet->jecScaleFactors.find("L1FastL2L3");
//     if (JEC != jet->jecScaleFactors.end()) {
//       float scale = JEC->second;
//       TLorentzVector corrP4 = scale*jet->momentum;
//       met.mEt.Set(met.metX() + jet->momentum.Px() - corrP4.Px(),
// 		  met.metY() + jet->momentum.Py() - corrP4.Py());
//       met.sumEt = met.sumEt - jet->momentum.Pt() + corrP4.Pt();
//     }
//     else {
//       cerr << "Error: L1FastL2L3 JEC not found for jet " << (iJet - jets.begin()) << " in event ";
//       cerr << evt << ".\n";
//     }
//   }
// }

// void GMSBAnalyzer::cleanJetCollection(const vector<susy::PFJet*>& allJets, 
// 				      vector<susy::PFJet*>& cleanedJets) const
// {
  
// }

void GMSBAnalyzer::setCanvasOptions(TCanvas& canvas, const string& xAxisTitle, 
				    const string& yAxisTitle, const float xAxisLowerLim, 
				    const float xAxisUpperLim, const float yAxisLowerLim, 
				    const float yAxisUpperLim, const bool setGrid) const
{
  canvas.SetFillStyle(0);
  canvas.SetFillColor(0);
  if (setGrid) canvas.SetGrid();
  TH1F* frame = canvas.DrawFrame(xAxisLowerLim, yAxisLowerLim, xAxisUpperLim, yAxisUpperLim);
  frame->GetXaxis()->SetTitle(xAxisTitle.c_str());
  frame->GetYaxis()->SetTitle(yAxisTitle.c_str());
}

float GMSBAnalyzer::fillWeightsHistograms(const TH1D* ggDiEMET, TH1D* controlDiEMET, 
					  TH1F& controlDiEMETScaled, TH1F& controlWeights) const
{
  float scale = 0.0;
  if (controlDiEMET->Integral() != 0.0) {
    scale = ggDiEMET->Integral()/controlDiEMET->Integral();
  }
  controlDiEMET->Scale(scale);
  controlDiEMETScaled.Add(controlDiEMET);
  controlWeights.Divide(ggDiEMET, &controlDiEMETScaled);
  return scale;
}

string GMSBAnalyzer::histName(const string& part1, const string& part2, 
			      const unsigned int part3) const
{
  stringstream stream;
  stream << part1 << part2 << part3;
  return stream.str();
}

void GMSBAnalyzer::generateBackgroundSubtractedSpectra(TH3F& METVsDiEMETVsInvMass, 
						       TH2F& bkgSubtractedMETVsDiEMET) const
{
  /*mass variable (fit in -10/+20 GeV window around the Z mass to cut off the trigger turn-on at 
    low mass)*/
  RooRealVar m("m", "m", 80.0, 110.0);

  //fit parameters
  /*from fitting ee with background fixed to values obtained by fitting background-only PDF with 
    all parameters floating to gg invariant mass spectrum in 1.14 fb^-1*/
  //fit converged
  //ZBias = 0.18 +/- 0.01
  //ZRes = 2.00 +/- 0.01
  //from fitting ee with all parameters floating in 1.14 fb^-1
  //fit converged
  //ZBias = 0.18 +/- 0.02
  //ZRes = 1.77 +/- 0.02
  RooRealVar ZMass("ZMass", "ZMass", 91.2, "GeV");            //Z mass (fixed to PDG)
  RooRealVar ZWidth("ZWidth", "ZWidth", 2.5, "GeV");          //natural Z width (fixed to PDG)
  RooRealVar ZBias("ZBias", "ZBias", 0.18, -1.0, 1.0);    //CMS energy scale (fixed)
  RooRealVar ZRes("ZRes", "ZRes", 1.77/*, 0.0, 4.0*/, "GeV"); //CMS resolution (fixed)

  // //unnecessary if not using Crystal Ball parametrization
  // RooRealVar alpha("alpha", "alpha", 5.0, 1.0, 7.0);     //Crystal Ball (floating)
  // RooRealVar n("n", "n", 50, 1, 200);                    //Crystal Ball (floating)

  /*from fitting background only, all parameters floating, to gg invariant mass spectrum in 1.14 
    fb^-1*/
  //fit converged
  //a = 82.2 +/- 0.7
  //b = 0.105 +/- 0.005
  //c = 0.033 +/- 0.004
  //peak = 63 +/- 16
  //from fitting ee with all parameters floating in 1.14 fb^-1
  //fit converged
  //aBkg = 86.8 +/- 0.4
  //bBkg = 0.30 +/- 0.02
  //c = 0.14 +/- 0.01
  RooRealVar aBkg("aBkg", "aBkg", /*86.8*/82.2/*, 50.0, 100.0*/, "GeV"); //RooCMSShape (fixed)
  RooRealVar bBkg("bBkg", "bBkg", /*0.3*/0.105/*, 0.1, 0.5*/);            //RooCMSShape (fixed)
  RooRealVar c("c", "c", /*0.14*/0.033/*, 0.01, 0.5*/);                   //RooCMSShape (fixed)
  RooRealVar peak("peak", "peak", 63.0, /*40.0, 90.0, */"GeV");  //RooCMSShape (fixed), but unused

  // //unnecessary if not using RooErfXGaussian parametrization
  // RooRealVar a("a", "a", 0.0, -10.0, 10.0); //RooCMSShape (fixed to tag and probe)
  // RooRealVar b("b", "b", 5.0, 5.0, 6.0);    //RooCMSShape (fixed to tag and probe)

  //signal and background yields (to be fitted)
  RooRealVar signalYield("signalYield", "signalYield", 0.0, 1.0e10);
  RooRealVar backgroundYield("backgroundYield", "backgroundYield", 0.0, 1.0e10);

  //signal PDF: Breit-Wigner (intrinsic) x Gaussian (resolution)
  RooBreitWigner intrinsic("intrinsic", "intrinsic", m, ZMass, ZWidth);
  RooGaussian resolution("resolution", "resolution", m, ZBias, ZRes);
  // RooCBShape resolution("resolution", "resolution", m, ZBias, ZRes, alpha, n);
  // RooErfXGaussian resolution("resolution", "resolution", m, /*a*/ZBias, b, ZBias, ZRes);
  RooFFTConvPdf signal("signal", "signal", m, intrinsic, resolution);

  //background PDF: RooCMSShape with peak position fixed to PDG Z mass
  RooCMSShape background("background", "background", m, aBkg, bBkg, c, ZMass/*peak*/);

  //signal + background PDF
  RooAddPdf total("total", "total", RooArgList(background, signal), 
  		  RooArgList(backgroundYield, signalYield));
  // RooFFTConvPdf total("total", "total", m, signal, background);

  //need to do a first call to fitTo for some reason, because first call always fucks up?
  TH1D* mDistDummy = METVsDiEMETVsInvMass.ProjectionX("dummy", 0, 0, 0, 0, "e");
  RooDataHist mHistDummy("mHistDummy", "mHistDummy", m, mDistDummy);
  total.fitTo(mHistDummy, Extended(true));

  //loop over (MET, di-EM ET) bins
  for (Int_t iDiEMETBin = 1; iDiEMETBin <= METVsDiEMETVsInvMass.GetNbinsY(); ++iDiEMETBin) {
    for (Int_t iMETBin = 1; iMETBin <= METVsDiEMETVsInvMass.GetNbinsZ(); ++iMETBin) {

      //fit invariant mass distribution in this (MET, di-EM ET) bin
      STRINGSTREAM name;
      name << "mDist_diEMETBin" << iDiEMETBin << "_METBin" << iMETBin;
      TH1D* mDist = METVsDiEMETVsInvMass.ProjectionX(name.str().c_str(), iDiEMETBin, iDiEMETBin, 
						     iMETBin, iMETBin, "e");
      RooDataHist mHist("mHist", "mHist", m, mDist);
      total.fitTo(mHist, Extended(true));

      //fill the background subtracted MET vs. di-EM ET plot with the signal yield
      bkgSubtractedMETVsDiEMET.SetBinContent(iDiEMETBin, iMETBin, signalYield.getVal());
      bkgSubtractedMETVsDiEMET.SetBinError(iDiEMETBin, iMETBin, signalYield.getError());

      //plot the signal and background fit over the data
      name.str("");
      name << "mCanvas_diEMETBin" << iDiEMETBin << "_METBin" << iMETBin;
      TCanvas canvas(name.str().c_str(), "", 600, 600);
      setCanvasOptions(canvas, "m (GeV)", "", 70.0, 110.0, 0.0, 10000.0, true);
      canvas.cd();
      name.str("");
      name << "mData_diEMETBin" << iDiEMETBin << "_METBin" << iMETBin;
      RooPlot* mPlot = m.frame(Name(name.str().c_str()));
      total.paramOn(mPlot, Format("NEU", AutoPrecision(1)), Layout(0.6, 0.98, 0.9));
      mHist.plotOn(mPlot, LineColor(kBlack));
      name.str("");
      name << "mBackgroundFit_diEMETBin" << iDiEMETBin << "_METBin" << iMETBin;
      total.plotOn(mPlot, Name(name.str().c_str()), Components("background"), LineColor(kRed));
      name.str("");
      name << "mSignalPlusBackgroundFit_diEMETBin" << iDiEMETBin << "_METBin" << iMETBin;
      total.plotOn(mPlot, Name(name.str().c_str()), LineColor(kBlue));
      mPlot->Draw();
      canvas.Write();
    }
  }
}

float GMSBAnalyzer::normAndErrorSquared(const TH3F& ggMETVsDiEMETVsInvMass, const TH1F& controlMET, 
					const TH1D* egMET, const unsigned int maxNormBin, 
					float& errSquared) const
{
  const float nGGMinusEW = 
    ggMETVsDiEMETVsInvMass.Integral(0, -1, 0, -1, 1, maxNormBin) - egMET->Integral(1, maxNormBin);
  const float nControl = controlMET.Integral(1, maxNormBin);
  float norm = 0.0;
  if (nControl != 0.0) norm = nGGMinusEW/nControl;
  float nGGErrSquared = 0.0;
  float nEGErrSquared = 0.0;
  for (unsigned int iBin = 1; iBin <= maxNormBin; ++iBin) {
    const float ggErr = 
      ggMETVsDiEMETVsInvMass.ProjectionZ("ggMET", 0, -1, 0, -1, "e")->GetBinError(iBin);
    const float egErr = egMET->GetBinError(iBin);
    nGGErrSquared+=(ggErr*ggErr);
    nEGErrSquared+=(egErr*egErr);
  }
  if (nGGMinusEW == 0.0) errSquared = 0.0;
  else {
    errSquared = 
      norm*norm*(((nGGErrSquared + nEGErrSquared)/(nGGMinusEW*nGGMinusEW)) + (1.0/nControl));
  }
  return norm;
}

void GMSBAnalyzer::bookToyDiEMETWeightsHistograms(const vector<TH1F*> controlWeights, 
						  const string& controlSample, 
						  vector<TH1F*>& controlDiEMETToyDistsByBin, 
						  const Int_t nMETBins) const
{
  for (Int_t iMETBin = 1; iMETBin <= nMETBins; ++iMETBin) {
    for (Int_t iBin = 1; iBin <= controlWeights[iMETBin - 1]->GetNbinsX(); ++iBin) {
      STRINGSTREAM name;
      name << controlSample << "DiEMETToyDistBin" << iBin << "_METBin" << iMETBin;
      float typicalVal = controlWeights[iMETBin - 1]->GetBinContent(iBin);
      float sqrtTypicalVal = sqrt(typicalVal);
      TH1F* hist = new TH1F(name.str().c_str(), "", 160, typicalVal - 8.0*sqrtTypicalVal, 
			    typicalVal + 8.0*sqrtTypicalVal);
      hist->Sumw2();
      controlDiEMETToyDistsByBin.push_back(hist);
    }
  }
}

void GMSBAnalyzer::makeToyDiEMETWeightsHistograms(TRandom3& random, TH1F& controlWeightsToy, 
						  const TH1F& controlWeights, 
						  vector<TH1F*>& controlDiEMETToyDistsByBin, 
						  const unsigned int iMETBin) const
{
  Int_t nDiEMETBins = controlWeights.GetNbinsX();
  for (int iBin = 1; iBin <= nDiEMETBins; ++iBin) {
    Double_t controlRandomWeight = 
      random.Gaus(controlWeights.GetBinContent(iBin), controlWeights.GetBinError(iBin));
    /*random.Poisson(controlWeights.GetBinContent(iBin));*/
    controlWeightsToy.SetBinContent(iBin, controlRandomWeight);
    controlDiEMETToyDistsByBin[(iMETBin - 1)*nDiEMETBins + iBin - 1]->Fill(controlRandomWeight);
  }
}

void GMSBAnalyzer::reweightDefault(const VFLOAT& controlMETVec, const VFLOAT& controlDiEMETVec, 
				   const TH1F& controlWeightsToy, TH1F* hist) const
{
  for (VFLOAT_IT i = controlMETVec.begin(); i != controlMETVec.end(); ++i) {
    Int_t iDiEMET = 1;
    bool foundDiEMETBin = false;
    while ((iDiEMET <= controlWeightsToy.GetNbinsX()) && !foundDiEMETBin) {
      const unsigned int iEvt = i - controlMETVec.begin();
      if ((controlDiEMETVec[iEvt] >= controlWeightsToy.GetBinLowEdge(iDiEMET)) && 
	  (controlDiEMETVec[iEvt] < controlWeightsToy.GetBinLowEdge(iDiEMET + 1))) {
	foundDiEMETBin = true;
      }
      else ++iDiEMET;
    }
    if (iDiEMET > controlWeightsToy.GetNbinsX()) {
      cerr << "Error: di-EM ET bin corresponding to event with MET = " << *i;
      cerr << " GeV not found in histogram " << controlWeightsToy.GetName() << ".\n";
    }
    hist->Fill(*i, controlWeightsToy.GetBinContent(iDiEMET));
  }
}

void GMSBAnalyzer::reweightBinned(const TH2F& controlMETVsDiEMET, 
				  const TH1F& controlWeightsToy, TH1F* hist, 
				  const unsigned int iMETBin) const
{
  STRINGSTREAM name;
  name << controlMETVsDiEMET.GetName() << "_MET";
  TH1D* controlMET = controlMETVsDiEMET.ProjectionY(name.str().c_str(), 0, -1, "e");
  const float METBinCenter = controlMET->GetBinCenter(iMETBin);
  const float METBinContent = controlMET->GetBinContent(iMETBin);
  name.str("");
  name << controlMETVsDiEMET.GetName() << "_diEMET_METBin" << iMETBin;
  TH1D* controlDiEMET = controlMETVsDiEMET.ProjectionX(name.str().c_str(), iMETBin, iMETBin, "e");
  const unsigned int nRolls = 1000;
  for (unsigned int iRoll = 1; iRoll <= nRolls; ++iRoll) {
    Double_t diEMETRandom = controlDiEMET->GetRandom();
    const unsigned int nDiEMETBins = controlWeightsToy.GetNbinsX();
    unsigned int iDiEMET = 1;
    bool foundDiEMETBin = false;
    while ((iDiEMET <= nDiEMETBins) && !foundDiEMETBin) {
      if ((diEMETRandom >= controlWeightsToy.GetBinLowEdge(iDiEMET)) && 
	  (diEMETRandom < controlWeightsToy.GetBinLowEdge(iDiEMET + 1))) {
	foundDiEMETBin = true;
      }
      else ++iDiEMET;
    }
    if (iDiEMET > nDiEMETBins) {
      cerr << "Di-EM ET generated larger than maximum di-EM ET in histogram ";
      cerr << controlWeightsToy.GetName() << "; assuming last bin.\n";
      iDiEMET = nDiEMETBins;
    }
    hist->Fill(METBinCenter, controlWeightsToy.GetBinContent(iDiEMET)*(METBinContent/nRolls));
  }
}

void GMSBAnalyzer::generateToys(vector<TH1F*>& controlFinalToy, 
				vector<TH1F*>& controlDiEMETToyDistsByBin, 
				const vector<TH1F*> controlWeights, const unsigned int nToys, 
				const string& controlSample, const Double_t* diEMETBins, 
				const unsigned int nMETBins, const Double_t* METBins, 
				const VFLOAT& controlMETVec, const VFLOAT& controlDiEMETVec) const
{
  bookToyDiEMETWeightsHistograms(controlWeights, controlSample, controlDiEMETToyDistsByBin);
  TRandom3 random;
  for (unsigned int iToy = 1; iToy <= nToys; ++iToy) {
    STRINGSTREAM nameWeights;
    nameWeights << controlSample << "WeightsToy" << iToy;
    TH1F controlWeightsToy(nameWeights.str().c_str(), "", controlWeights[0]->GetNbinsX(), 
			   diEMETBins);
    controlWeightsToy.Sumw2();
    STRINGSTREAM nameFinal;
    nameFinal << controlSample << "FinalToy" << iToy;
    TH1F* hist = new TH1F(nameFinal.str().c_str(), "", nMETBins, METBins);
    hist->Sumw2();
    hist->SetFillStyle(0);
    makeToyDiEMETWeightsHistograms(random, controlWeightsToy, *controlWeights[0], 
				   controlDiEMETToyDistsByBin);
    reweightDefault(controlMETVec, controlDiEMETVec, controlWeightsToy, hist);
    controlFinalToy.push_back(hist);
  }
}

void GMSBAnalyzer::generateToys(vector<TH1F*>& controlFinalToy, 
				vector<TH1F*>& controlDiEMETToyDistsByBin, 
				const vector<TH1F*> controlWeights, const unsigned int nToys, 
				const string& controlSample, const Double_t* diEMETBins, 
				const unsigned int nMETBins, const Double_t* METBins, 
				const TH2F& controlMETVsDiEMET) const
{
  bookToyDiEMETWeightsHistograms(controlWeights, controlSample, controlDiEMETToyDistsByBin, 
				 nMETBins);
  TRandom3 random;
  for (unsigned int iToy = 1; iToy <= nToys; ++iToy) {
    STRINGSTREAM nameFinal;
    nameFinal << controlSample << "FinalToy" << iToy;
    TH1F* hist = new TH1F(nameFinal.str().c_str(), "", nMETBins, METBins);
    hist->Sumw2();
    hist->SetFillStyle(0);
    for (unsigned int iMETBin = 1; iMETBin <= nMETBins; ++iMETBin) {
      STRINGSTREAM nameWeights;
      nameWeights << controlSample << "WeightsToy" << iToy << "_METBin" << iMETBin;
      TH1F controlWeightsToy(nameWeights.str().c_str(), "", 
			     controlWeights[iMETBin - 1]->GetNbinsX(), diEMETBins);
      controlWeightsToy.Sumw2();
      makeToyDiEMETWeightsHistograms(random, controlWeightsToy, *controlWeights[iMETBin - 1], 
				     controlDiEMETToyDistsByBin, iMETBin);
      reweightBinned(controlMETVsDiEMET, controlWeightsToy, hist, iMETBin);
      controlFinalToy.push_back(hist);
    }
  }
}

void GMSBAnalyzer::fillToyDistributions(vector<TH1F*>& controlMETToyDistsByBin, 
					const TH1F& controlFinal, TCanvas& controlToyCanvas, 
					vector<TH1F*>& controlFinalToy, const unsigned int nToys, 
					const string& controlSample, 
					const unsigned int nMETBins) const
{
  for (unsigned int iBin = 1; iBin <= nMETBins; ++iBin) {
    STRINGSTREAM name;
    name << controlSample << "METToyDistBin" << iBin;
    float typicalVal = controlFinal.GetBinContent(iBin);
    float sqrtTypicalVal = sqrt(typicalVal);
    TH1F* hist = new TH1F(name.str().c_str(), "", 150, typicalVal - 15.0*sqrtTypicalVal, 
			  typicalVal + 15.0*sqrtTypicalVal);
    hist->Sumw2();
    controlMETToyDistsByBin.push_back(hist);
  }
  controlToyCanvas.cd();
  for (unsigned int iToy = 0; iToy < nToys; ++iToy) {
    if (iToy == 0) controlFinalToy[iToy]->Draw("HIST");
    else controlFinalToy[iToy]->Draw("HISTSAME");
    for (unsigned int iBin = 1; iBin <= nMETBins; ++iBin) {
      controlMETToyDistsByBin[iBin - 1]->Fill(controlFinalToy[iToy]->GetBinContent(iBin));
    }
  }
}

void GMSBAnalyzer::setMETErrorBars(TH1F& controlFinal, 
				   const vector<TH1F*>& controlMETToyDistsByBin, 
				   const vector<TH1F*>& controlLowSidebandMETToyDistsByBin, 
				   const vector<TH1F*>& controlHighSidebandMETToyDistsByBin, 
				   const float controlNorm, 
				   const float controlNormErrSquared, const bool doDefault) const
{
  for (int iBin = 1; iBin <= controlFinal.GetNbinsX(); ++iBin) {
    const float reweightingErr = controlMETToyDistsByBin[iBin - 1]->GetRMS();
    float reweightingErrSquared = reweightingErr*reweightingErr;
    if ((STRING(controlFinal.GetName()) == "eeFinal") && doDefault) {
      const float lowSidebandReweightingErr = 
	controlLowSidebandMETToyDistsByBin[iBin - 1]->GetRMS();
      const float highSidebandReweightingErr = 
	controlHighSidebandMETToyDistsByBin[iBin - 1]->GetRMS();
      reweightingErrSquared+=4*(lowSidebandReweightingErr*lowSidebandReweightingErr + 
				highSidebandReweightingErr*highSidebandReweightingErr);
    }
    controlFinal.SetBinError(iBin, 
			     //stat
			     sqrt(controlFinal.GetBinError(iBin)*controlFinal.GetBinError(iBin) + 
			     controlNormErrSquared/*norm*/ + 
			     controlNorm*controlNorm*reweightingErrSquared/*reweighting*/));
  }
}

void GMSBAnalyzer::makeFinalCanvas(TH1* hist, const Color_t lineColor, const Width_t lineWidth, 
				   const Style_t fillStyle, const Color_t fillColor, 
				   const Size_t markerSize, const string& drawOption) const
{
  hist->SetLineColor(lineColor);
  hist->SetLineWidth(lineWidth);
  hist->SetFillStyle(fillStyle);
  hist->SetFillColor(fillColor);
  hist->SetMarkerSize(markerSize);
  hist->Draw(drawOption.c_str());
}

void GMSBAnalyzer::deallocateMemory(VTH1F& vec) const
{
  for (VTH1F_IT i = vec.begin(); i != vec.end(); ++i) {
    delete *i;
    *i = NULL;
  }
}

void GMSBAnalyzer::runMETAnalysis(const std::string& outputFile)
{
   if (fChain == 0) return;

   //open file
   TFile out(outputFile.c_str(), "RECREATE");
   out.cd();

   //define constants
   const bool kSumW2 = true;
   const bool kSetGrid = true;
   const unsigned int maxNormBin = 4;     //normalization region
   const float nToys = 1000;              //number of toys for the MET shape error from reweighting
   const unsigned int nMETBins = 13;      //number of MET bins
   const unsigned int nDiEMETBins = 25;   //number of di-EM ET bins
   const unsigned int nInvMassBins = 150; //number of invariant mass bins
   const Double_t METBins[14] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 70.0, 
				 100.0, 150.0, 500.0};            //MET bin boundaries
   const Double_t diEMETBins[26] = {0.0,                          //di-EM ET bin boundaries
				    3.0, 6.0, 9.0, 12.0, 15.0, 
				    18.0, 21.0, 24.0, 27.0, 30.0, 
				    35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 
				    70.0, 80.0, 90.0, 
				    100.0, 150.0, 200.0, 
				    300.0, 400.0, 500.0};
   const Double_t invMassBins[151] = {0.0,                        //invariant mass bin boundaries
				      2.0, 4.0, 6.0, 8.0, 10.0, 
				      12.0, 14.0, 16.0, 18.0, 20.0, 
				      22.0, 24.0, 26.0, 28.0, 30.0, 
				      32.0, 34.0, 36.0, 38.0, 40.0, 
				      42.0, 44.0, 46.0, 48.0, 50.0, 
				      52.0, 54.0, 56.0, 58.0, 60.0, 
				      62.0, 64.0, 66.0, 68.0, 70.0, 
				      72.0, 74.0, 76.0, 78.0, 80.0, 
				      82.0, 84.0, 86.0, 88.0, 90.0, 
				      92.0, 94.0, 96.0, 98.0, 100.0, 
				      102.0, 104.0, 106.0, 108.0, 110.0, 
				      112.0, 114.0, 116.0, 118.0, 120.0, 
				      122.0, 124.0, 126.0, 128.0, 130.0, 
				      132.0, 134.0, 136.0, 138.0, 140.0, 
				      142.0, 144.0, 146.0, 148.0, 150.0, 
				      152.0, 154.0, 156.0, 158.0, 160.0, 
				      162.0, 164.0, 166.0, 168.0, 170.0, 
				      172.0, 174.0, 176.0, 178.0, 180.0, 
				      182.0, 184.0, 186.0, 188.0, 190.0, 
				      192.0, 194.0, 196.0, 198.0, 200.0, 
				      202.0, 204.0, 206.0, 208.0, 210.0, 
				      212.0, 214.0, 216.0, 218.0, 220.0, 
				      222.0, 224.0, 226.0, 228.0, 230.0, 
				      232.0, 234.0, 236.0, 238.0, 240.0, 
				      242.0, 244.0, 246.0, 248.0, 250.0, 
				      252.0, 254.0, 256.0, 258.0, 260.0, 
				      262.0, 264.0, 266.0, 268.0, 270.0, 
				      272.0, 274.0, 276.0, 278.0, 280.0, 
				      282.0, 284.0, 286.0, 288.0, 290.0, 
				      292.0, 294.0, 296.0, 298.0, 300.0};
   const float egMisIDRate = 0.014;      /*take from CMS AN-2010/294 for now; can be computed with 
					   Z*/
   const float egMisIDRateErr = 0.002;   /*take from CMS AN-2010/294 for now; can be computed with 
					   Z*/

   //combined isolation histograms (sanity check)
   TH1F ggCombinedIso("ggCombinedIso", "", 40, 0.0, 20.0);
   TH1F egCombinedIso("egCombinedIso", "", 40, 0.0, 20.0);
   TH1F eeCombinedIso("eeCombinedIso", "", 40, 0.0, 20.0);
   TH1F ffCombinedIso("ffCombinedIso", "", 40, 0.0, 20.0);
   setHistogramOptions(ggCombinedIso, "Combined isolation (GeV)", "", "", kSumW2);
   setHistogramOptions(egCombinedIso, "Combined isolation (GeV)", "", "", kSumW2);
   setHistogramOptions(eeCombinedIso, "Combined isolation (GeV)", "", "", kSumW2);
   setHistogramOptions(ffCombinedIso, "Combined isolation (GeV)", "", "", kSumW2);

   //leading photon ET histograms (sanity check)
   TH1F ggLeadingPhotonET("ggLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F egLeadingPhotonET("egLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F eeLeadingPhotonET("eeLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F ffLeadingPhotonET("ffLeadingPhotonET", "", 150, 0.0, 300.0);
   setHistogramOptions(ggLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(egLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);

   //nPV histograms (sanity check)
   TH1F ggNPV("ggNPV", "", 36, -0.5, 35.5);
   TH1F egNPV("egNPV", "", 36, -0.5, 35.5);
   TH1F eeNPV("eeNPV", "", 36, -0.5, 35.5);
   TH1F ffNPV("ffNPV", "", 36, -0.5, 35.5);
   setHistogramOptions(ggNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(egNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(eeNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(ffNPV, "n_{PV}", "", "", kSumW2);

   //MET vs. di-EM ET vs. invariant mass histograms
   TH3F ggMETVsDiEMETVsInvMass("ggMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F egMETVsDiEMETVsInvMass("egMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F eeMETVsDiEMETVsInvMass("eeMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F eeLowSidebandMETVsDiEMETVsInvMass("eeLowSidebandMETVsDiEMETVsInvMass", "", nInvMassBins, 
					  invMassBins, nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F eeHighSidebandMETVsDiEMETVsInvMass("eeHighSidebandMETVsDiEMETVsInvMass", "", nInvMassBins, 
					   invMassBins, nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F ffMETVsDiEMETVsInvMass("ffMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   setHistogramOptions(ggMETVsDiEMETVsInvMass, "m_{#gamma#gamma} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(egMETVsDiEMETVsInvMass, "m_{e#gamma} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(eeMETVsDiEMETVsInvMass, "m_{ee} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(eeLowSidebandMETVsDiEMETVsInvMass, "m_{ee} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(eeHighSidebandMETVsDiEMETVsInvMass, "m_{ee} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(ffMETVsDiEMETVsInvMass, "m_{ff} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);

   //reweighted and normalized ee and ff MET histograms
   TH1F eeFinal("eeFinal", "", nMETBins, METBins);
   TH1F eeLowSidebandFinal("eeLowSidebandFinal", "", nMETBins, METBins);
   TH1F eeHighSidebandFinal("eeHighSidebandFinal", "", nMETBins, METBins);
   TH1F ffFinal("ffFinal", "", nMETBins, METBins);
   setHistogramOptions(eeFinal, "ME_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeLowSidebandFinal, "ME_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeHighSidebandFinal, "ME_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffFinal, "ME_{T} (GeV)", "", "", kSumW2);

   //ee and ff HT vs. MET histograms
   TH2F eeHTVsMET("eeHTVsMET", "", nMETBins, METBins, 60, 0.0, 600.0);
   TH2F ffHTVsMET("ffHTVsMET", "", nMETBins, METBins, 60, 0.0, 600.0);
   setHistogramOptions(eeHTVsMET, "ME_{T} (GeV)", "H_{T} (GeV)", "", kSumW2);
   setHistogramOptions(ffHTVsMET, "ME_{T} (GeV)", "H_{T} (GeV)", "", kSumW2);

   //ee and ff MET vs. no. jets histograms
   TH2F eeMETVsNJets30("eeMETVsNJets30", "", 6, -0.5, 5.5, nMETBins, METBins);
   TH2F ffMETVsNJets30("ffMETVsNJets30", "", 6, -0.5, 5.5, nMETBins, METBins);
   TH2F eeMETVsNJets60("eeMETVsNJets60", "", 6, -0.5, 5.5, nMETBins, METBins);
   TH2F ffMETVsNJets60("ffMETVsNJets60", "", 6, -0.5, 5.5, nMETBins, METBins);
   setHistogramOptions(eeMETVsNJets30, "n_{j} (30-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
   setHistogramOptions(ffMETVsNJets30, "n_{j} (30-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
   setHistogramOptions(eeMETVsNJets60, "n_{j} (60-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
   setHistogramOptions(ffMETVsNJets60, "n_{j} (60-60 GeV)", "ME_{T} (GeV)", "", kSumW2);

   /*histogram of the percentage difference between Poisson mean of generated di-EM ET weights and 
     input mean (i.e. the value of the measured weight)*/
   TH1F ffDiEMETInputVsOutput("ffDiEMETInputVsOutput", "", nDiEMETBins, diEMETBins);
   TH1F eeDiEMETInputVsOutput("eeDiEMETInputVsOutput", "", nDiEMETBins, diEMETBins);
   TH1F eeLowSidebandDiEMETInputVsOutput("eeLowSidebandDiEMETInputVsOutput", "", nDiEMETBins, 
					 diEMETBins);
   TH1F eeHighSidebandDiEMETInputVsOutput("eeHighSidebandDiEMETInputVsOutput", "", nDiEMETBins, 
					  diEMETBins);
   setHistogramOptions(ffDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);
   setHistogramOptions(eeDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);
   setHistogramOptions(eeLowSidebandDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);
   setHistogramOptions(eeHighSidebandDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);

   //canvases for the MET toy plots
   TCanvas ffToyCanvas("ffToyCanvas", "", 600, 600);
   TCanvas eeToyCanvas("eeToyCanvas", "", 600, 600);
   TCanvas eeLowSidebandToyCanvas("eeLowSidebandToyCanvas", "", 600, 600);
   TCanvas eeHighSidebandToyCanvas("eeHighSidebandToyCanvas", "", 600, 600);
   setCanvasOptions(ffToyCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);
   setCanvasOptions(eeToyCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);
   setCanvasOptions(eeLowSidebandToyCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, 
		    kSetGrid);
   setCanvasOptions(eeHighSidebandToyCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, 
		    kSetGrid);

   //canvas for the final plot
   TCanvas METCanvas("METCanvas", "", 600, 600);
   setCanvasOptions(METCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);

   //MET and di-EM ET vectors
   VFLOAT eeMETVec;
   VFLOAT eeLowSidebandMETVec;
   VFLOAT eeHighSidebandMETVec;
   VFLOAT ffMETVec;
   VFLOAT eeDiEMETVec;
   VFLOAT eeLowSidebandDiEMETVec;
   VFLOAT eeHighSidebandDiEMETVec;
   VFLOAT ffDiEMETVec;

   //PU rho
   TH1F rhoSkim("rhoSkim", "", 100, 0.0, 10.0);
   TH1F rhoPreselected("rhoPreselected", "", 100, 0.0, 10.0);
   setHistogramOptions(rhoSkim, "#rho (GeV)", "", "", kSumW2);
   setHistogramOptions(rhoPreselected, "#rho (GeV)", "", "", kSumW2);

   //set up on-the-fly jet corrections for PF jets
   vector<JetCorrectorParameters> PFJECs;
   PFJECs.push_back(JetCorrectorParameters(L1JECFile_));
   PFJECs.push_back(JetCorrectorParameters(L2JECFile_));
   PFJECs.push_back(JetCorrectorParameters(L3JECFile_));
   FactorizedJetCorrector PFJetCorrector(PFJECs);

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

     //get int. lumi weight for this dataset
     TString fileName(fChain->GetCurrentFile()->GetName());
     TString dataset;
     /*if (fileName.Contains("ntuple_.*-v[0-9]")) */dataset = fileName("ntuple_.*-v[0-9]");
     if (dataset.Length() >= 7) dataset.Remove(0, 7);
     if ((dataset.Length() >= 3) && 
	 (dataset.Contains("v2_.*-v[0-9]"))) dataset.Remove(0, 3); //for TT sample
     string datasetString((const char*)dataset);
     map<string, unsigned int>::const_iterator iDataset = fileMap_.find(datasetString);
     float lumiWeight = 1.0;
     if (iDataset != fileMap_.end()) lumiWeight = weight_[iDataset->second];
     else {
       /*once comfortable with MC dataset names, remove this and just assume the file is data and 
     	 gets weight 1.0*/
       cerr << "Error: dataset " << datasetString << " was not entered into the map.\n";
     }

     //get PU weight for this dataset
     /*in-time reweighting only following 
       https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities as of 27-Oct-11 
       and Tessa's recommendation*/
     int nPV  = -1;
     float PUWeight = 1.0;
     susy::PUSummaryInfoCollection::const_iterator iBX = susyEvent->PU.begin();
     bool foundInTimeBX = false;
     while ((iBX != susyEvent->PU.end()) && !foundInTimeBX) {
       if (iBX->BX == 0) { 
     	 nPV = iBX->numInteractions;
     	 foundInTimeBX = true;
       }
       ++iBX;
     }
     PUWeight = lumiWeights_.ITweight(nPV);

     //fill MET vs. di-EM ET, invariant mass, rho, leading photon ET, and nPV histograms
     double invMass = -1.0;
     if (susyCategory->getEvtInvMass()->size() > 0) invMass = susyCategory->getEvtInvMass(tag_);
     else {
       cerr << "Error: susyCategory->getEvtInvMass()->size() <= 0 in event " << (jentry + 1);
       cerr << ".  Using invariant mass -1.0 GeV.\n";
     }
     double MET = -1.0;
     map<TString, susy::MET>::const_iterator iMET = susyEvent->metMap.find("pfMet");
     if (iMET != susyEvent->metMap.end()) MET = iMET->second.met();
     else cerr << "Error: PFMET not found in event " << (jentry + 1) << ".\n";
     double diEMET = -1.0;
     if (susyCategory->getEvtDiEMET()->size() > 0) diEMET = susyCategory->getEvtDiEMET(tag_);
     else {
       cerr << "Error: susyCategory->getEvtDiEMET()->size() <= 0 in event " << (jentry + 1);
       cerr << ".  Using di-EM ET -1.0 GeV.\n";
     }
     const Float_t rho = susyEvent->rho;
     rhoSkim.Fill(rho, lumiWeight*PUWeight);
     int evtCategory = FAIL;
     if (susyCategory->getEventCategory()->size() > 0) {
       evtCategory = susyCategory->getEventCategory(tag_);
     }
     else {
       cerr << "Error: susyCategory->getEventCategory()->size() <= 0 in event " << (jentry + 1);
       cerr << ".  Using event category FAIL.\n";
     }
     if (evtCategory != FAIL) rhoPreselected.Fill(rho, lumiWeight*PUWeight);
     unsigned int nGoodRecoPV = 0;
     for (vector<susy::Vertex>::const_iterator iPV = susyEvent->vertices.begin(); 
	  iPV != susyEvent->vertices.end(); ++iPV) {
       if (!iPV->isFake() && (iPV->ndof > 4) && (iPV->position.Z() <= 24.0/*cm*/) && 
	   (iPV->position.Perp() <= 2.0/*cm*/)) ++nGoodRecoPV;
     }

     //count jets
     unsigned int nJets = 0;
     unsigned int nJets30 = 0;
     unsigned int nJets60 = 0;
     vector<float> jetET;
     map<TString, susy::PhotonCollection>::const_iterator iPhotonMap = 
       susyEvent->photons.find((const char*)tag_);
     map<TString, susy::PFJetCollection>::const_iterator iJets = susyEvent->pfJets.find("ak5");
     if ((iJets != susyEvent->pfJets.end()) && (iPhotonMap != susyEvent->photons.end())) {

       //loop over PF jets
       for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
	    iJet != iJets->second.end(); ++iJet) {

	 //compute corrected P4
	 PFJetCorrector.setJetEta(iJet->momentum.Eta());
	 PFJetCorrector.setJetPt(iJet->momentum.Pt());
	 PFJetCorrector.setJetA(iJet->jetArea);
	 PFJetCorrector.setRho(susyEvent->rho);
	 float storedCorr = -1.0;
	 map<TString, Float_t>::const_iterator iCorr = iJet->jecScaleFactors.find("L1FastL2L3");
	 if (iCorr != iJet->jecScaleFactors.end()) storedCorr = iCorr->second;
	 TLorentzVector corrP4 = storedCorr*iJet->momentum;

	 //calculate the jet-EM overlap
	 float EM1Eta = -9.0;
	 float EM1Phi = -9.0;
	 float EM2Eta = -9.0;
	 float EM2Phi = -9.0;
	 unsigned int count = 0;
	 for (susy::PhotonCollection::const_iterator iPhoton = iPhotonMap->second.begin(); 
	      iPhoton != iPhotonMap->second.end(); ++iPhoton) {
	   if (susyCategory->getIsDeciding()->size() > 0) {
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
	   }
	   else {
	     cerr << "Error: susyCategory->getIsDeciding()->size() <= 0 in event " << (jentry + 1);
	     cerr << ".  Assuming no deciding photons.\n";
	   }
	 }
	 if ((count != 2) && (evtCategory != FAIL)) {
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
	     ++nJets;
	     jetET.push_back(corrP4.Et());
	     if ((corrP4.Et() >= 30.0/*GeV*/) && (corrP4.Et() < 60.0/*GeV*/)) ++nJets30;
	     else if (corrP4.Et() >= 60.0/*GeV*/) ++nJets60;
	   }
	 }
       }
     }
     else {
       cerr << "Error: " << tag_ << " photon collection or ak5 jet collection not found in ";
       cerr << "event " << (jentry + 1) << ".\n";
     }
     const float HT = accumulate(jetET.begin(), jetET.end(), 0);

     //get combined isolation and photon ET
     vector<float> combinedIso;
     vector<float> ET;
     if (evtCategory != FAIL) {
       if (iPhotonMap != susyEvent->photons.end()) {
	 for (susy::PhotonCollection::const_iterator iPhoton = iPhotonMap->second.begin(); 
	      iPhoton != iPhotonMap->second.end(); ++iPhoton) {
	   if (susyCategory->getIsDeciding()->size() > 0) {
	     if (susyCategory->getIsDeciding(tag_, iPhoton - iPhotonMap->second.begin())) {
	       ET.push_back(iPhoton->momentum.Et());
	       combinedIso.push_back(iPhoton->ecalRecHitSumEtConeDR03 - 0.1474*rho + 
				     iPhoton->hcalTowerSumEtConeDR03() - 0.0467*rho + 
				     iPhoton->trkSumPtHollowConeDR03);
	     }
	   }
	   else {
	     cerr << "Error: susyCategory->getIsDeciding()->size() <= 0 in event " << (jentry + 1);
	     cerr << ".  Assuming no deciding photons.\n";
	   }
	 }
	 if ((combinedIso.size() != 2) || (ET.size() != 2)) {
	   cerr << "Error: combinedIso.size() = " << combinedIso.size() << " and ET.size() = ";
	   cerr << ET.size() << " in event " << (jentry + 1) << ".\n";
	 }
       }
       else {
	 cerr << "Error: " << tag_ << " photon collection not found in event " << (jentry + 1);
	 cerr << ".\n";
       }
     }

     //sort photon ET in ascending order
     sort(ET.begin(), ET.end());

//      if (nJets == 0) {
     bool passTightCombinedIso = true;
       switch (evtCategory) {
       case GG:
	 ggMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	 for (vector<float>::const_iterator iIso = combinedIso.begin(); 
	      iIso != combinedIso.end(); ++iIso) {
	   ggCombinedIso.Fill(*iIso, lumiWeight*PUWeight);
	 }
	 ggLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 ggNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 break;
       case EG:
	 egMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	 for (vector<float>::const_iterator iIso = combinedIso.begin(); 
	      iIso != combinedIso.end(); ++iIso) {
	   egCombinedIso.Fill(*iIso, lumiWeight*PUWeight);
	 }
	 egLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 egNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 break;
       case EE:
	 if ((invMass >= 71.0/*GeV*/) && (invMass < 76.0/*GeV*/)) {
	   eeLowSidebandMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	   eeLowSidebandMETVec.push_back(MET);
	   eeLowSidebandDiEMETVec.push_back(diEMET);
	 }
	 if ((invMass >= 81.0/*GeV*/) && (invMass < 101.0/*GeV*/)) {
	   eeMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	   eeMETVec.push_back(MET);
	   eeDiEMETVec.push_back(diEMET);
	   eeHTVsMET.Fill(MET, HT, lumiWeight*PUWeight);
	   eeMETVsNJets30.Fill(nJets30, MET, lumiWeight*PUWeight);
	   eeMETVsNJets60.Fill(nJets60, MET, lumiWeight*PUWeight);
	 }
	 if ((invMass >= 106.0/*GeV*/) && (invMass < 111.0/*GeV*/)) {
	   eeHighSidebandMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	   eeHighSidebandMETVec.push_back(MET);
	   eeHighSidebandDiEMETVec.push_back(diEMET);
	 }
	 for (vector<float>::const_iterator iIso = combinedIso.begin(); 
	      iIso != combinedIso.end(); ++iIso) {
	   eeCombinedIso.Fill(*iIso, lumiWeight*PUWeight);
	 }
	 eeLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 eeNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 break;
       case FF:
	 for (vector<float>::const_iterator iIso = combinedIso.begin(); 
	      iIso != combinedIso.end(); ++iIso) {
	   ffCombinedIso.Fill(*iIso, lumiWeight*PUWeight);
	   passTightCombinedIso = passTightCombinedIso && (*iIso < 9.0/*GeV*/);
	 }
	 ffMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	 ffMETVec.push_back(MET);
	 ffDiEMETVec.push_back(diEMET);
	 ffLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 ffNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 ffHTVsMET.Fill(MET, HT, lumiWeight*PUWeight);
	 ffMETVsNJets30.Fill(nJets30, MET, lumiWeight*PUWeight);
	 ffMETVsNJets60.Fill(nJets60, MET, lumiWeight*PUWeight);
	 break;
       default:
	 break;
       }
//      }
   }

   //fill weights histograms
   vector<TH1F*> eeDiEMETScaled;
   vector<TH1F*> eeLowSidebandDiEMETScaled;
   vector<TH1F*> eeHighSidebandDiEMETScaled;
   vector<TH1F*> ffDiEMETScaled;
   vector<TH1F*> eeWeights;
   vector<TH1F*> eeLowSidebandWeights;
   vector<TH1F*> eeHighSidebandWeights;
   vector<TH1F*> ffWeights;
   vector<string> controlSamples;
   controlSamples.push_back("ee");
   controlSamples.push_back("eeLowSideband");
   controlSamples.push_back("eeHighSideband");
   controlSamples.push_back("ff");
   vector<TH3F*> TH3Fs;
   TH3Fs.push_back(&eeMETVsDiEMETVsInvMass);
   TH3Fs.push_back(&eeLowSidebandMETVsDiEMETVsInvMass);
   TH3Fs.push_back(&eeHighSidebandMETVsDiEMETVsInvMass);
   TH3Fs.push_back(&ffMETVsDiEMETVsInvMass);
   vector<vector<TH1F*>* > diEMETScaled;
   diEMETScaled.push_back(&eeDiEMETScaled);
   diEMETScaled.push_back(&eeLowSidebandDiEMETScaled);
   diEMETScaled.push_back(&eeHighSidebandDiEMETScaled);
   diEMETScaled.push_back(&ffDiEMETScaled);
   vector<vector<TH1F*>* > weights;
   weights.push_back(&eeWeights);
   weights.push_back(&eeLowSidebandWeights);
   weights.push_back(&eeHighSidebandWeights);
   weights.push_back(&ffWeights);
   for (vector<string>::const_iterator iControlSample = controlSamples.begin(); 
	iControlSample != controlSamples.end(); ++iControlSample) {
     TH1F* histDiEMETScaled = new TH1F(histName(*iControlSample, "DiEMETScaled_METBin", 
						1).c_str(), "", nDiEMETBins, diEMETBins);
     TH1F* histWeights = new TH1F(histName(*iControlSample, "Weights_METBin", 1).c_str(), "", 
				  nDiEMETBins, diEMETBins);
     setHistogramOptions(*histDiEMETScaled, "Di-EM E_{T} (GeV)", "", "", kSumW2);
     setHistogramOptions(*histWeights, "Di-EM E_{T} (GeV)", "", "", kSumW2);
     const unsigned int i = iControlSample - controlSamples.begin();
     fillWeightsHistograms(ggMETVsDiEMETVsInvMass.ProjectionY(histName("ggDiEMET_METBin", "", 
								       1).c_str(), 0, -1, 0, -1, 
							      "e"), 
			   TH3Fs[i]->ProjectionY(histName(*iControlSample, "DiEMET_METBin", 
							  1).c_str(), 0, -1, 0, -1, "e"), 
			   *histDiEMETScaled, *histWeights);
     diEMETScaled[i]->push_back(histDiEMETScaled);
     weights[i]->push_back(histWeights);
   }

   //generate toys for calculating error due to MET shape from reweighting
   vector<TH1F*> ffFinalToy;
   vector<TH1F*> eeFinalToy;
   vector<TH1F*> eeLowSidebandFinalToy;
   vector<TH1F*> eeHighSidebandFinalToy;
   vector<TH1F*> ffDiEMETToyDistsByBin;
   vector<TH1F*> eeDiEMETToyDistsByBin;
   vector<TH1F*> eeLowSidebandDiEMETToyDistsByBin;
   vector<TH1F*> eeHighSidebandDiEMETToyDistsByBin;
   generateToys(ffFinalToy, ffDiEMETToyDistsByBin, ffWeights, nToys, "ff", diEMETBins, nMETBins, 
   		METBins, ffMETVec, ffDiEMETVec);
   generateToys(eeFinalToy, eeDiEMETToyDistsByBin, eeWeights, nToys, "ee", diEMETBins, nMETBins, 
   		METBins, eeMETVec, eeDiEMETVec);
   generateToys(eeLowSidebandFinalToy, eeLowSidebandDiEMETToyDistsByBin, eeLowSidebandWeights, 
   		nToys, "eeLowSideband", diEMETBins, nMETBins, METBins, eeLowSidebandMETVec, 
   		eeLowSidebandDiEMETVec);
   generateToys(eeHighSidebandFinalToy, eeHighSidebandDiEMETToyDistsByBin, eeHighSidebandWeights, 
   		nToys, "eeHighSideband", diEMETBins, nMETBins, METBins, eeHighSidebandMETVec, 
   		eeHighSidebandDiEMETVec);

   //reweight ee and ff MET histograms
   reweightDefault(eeLowSidebandMETVec, eeLowSidebandDiEMETVec, *eeLowSidebandWeights[0], 
   		   &eeLowSidebandFinal);
   reweightDefault(eeMETVec, eeDiEMETVec, *eeWeights[0], &eeFinal);
   reweightDefault(eeHighSidebandMETVec, eeHighSidebandDiEMETVec, *eeHighSidebandWeights[0], 
   		   &eeHighSidebandFinal);
   reweightDefault(ffMETVec, ffDiEMETVec, *ffWeights[0], &ffFinal);

   //make individual histograms of MET toy distributions, 1 per MET bin
   vector<TH1F*> ffMETToyDistsByBin;
   vector<TH1F*> eeMETToyDistsByBin;
   vector<TH1F*> eeLowSidebandMETToyDistsByBin;
   vector<TH1F*> eeHighSidebandMETToyDistsByBin;
   fillToyDistributions(ffMETToyDistsByBin, ffFinal, ffToyCanvas, ffFinalToy, nToys, "ff", 
   			nMETBins);
   fillToyDistributions(eeMETToyDistsByBin, eeFinal, eeToyCanvas, eeFinalToy, nToys, "ee", 
   			nMETBins);
   fillToyDistributions(eeLowSidebandMETToyDistsByBin, eeLowSidebandFinal, 
   			eeLowSidebandToyCanvas, eeLowSidebandFinalToy, nToys, "eeLowSideband", 
   			nMETBins);
   fillToyDistributions(eeHighSidebandMETToyDistsByBin, eeHighSidebandFinal, 
   			eeHighSidebandToyCanvas, eeHighSidebandFinalToy, nToys, "eeHighSideband", 
   			nMETBins);

   //calculate EW contribution
   TH1D* egMET = egMETVsDiEMETVsInvMass.ProjectionZ("egMET", 0, -1, 0, -1, "e");
   const float egScale = egMisIDRate/(1.0 - egMisIDRate);
   egMET->Scale(egScale);
   for (Int_t iBin = 1; iBin <= egMET->GetNbinsX(); ++iBin) {
     const float binContent = egMET->GetBinContent(iBin);
     const float binPoissonError = egMET->GetBinError(iBin);
     if (binContent == 0.0) egMET->SetBinError(iBin, 0.0);
     else {
       egMET->SetBinError(iBin, 
   			  egScale*sqrt(((binPoissonError*binPoissonError)/
   					(binContent*binContent)) + 
   				       ((egMisIDRateErr*egMisIDRateErr)/
   					(egMisIDRate*egMisIDRate))));
     }
   }

   //normalize ee and ff MET histograms
   eeLowSidebandFinal.Scale(2.0);
   eeHighSidebandFinal.Scale(2.0);
   eeFinal.Add(&eeLowSidebandFinal, -1.0);
   eeFinal.Add(&eeHighSidebandFinal, -1.0);
   float eeNormErrSquared = 0.0;
   float ffNormErrSquared = 0.0;
   const float eeNorm = 
     normAndErrorSquared(ggMETVsDiEMETVsInvMass, eeFinal, egMET, maxNormBin, eeNormErrSquared);
   const float ffNorm = 
     normAndErrorSquared(ggMETVsDiEMETVsInvMass, ffFinal, egMET, maxNormBin, ffNormErrSquared);
   eeFinal.Scale(eeNorm);
   ffFinal.Scale(ffNorm);

   //set MET error bars
   setMETErrorBars(ffFinal, ffMETToyDistsByBin, vector<TH1F*>(), vector<TH1F*>(), ffNorm, 
   		   ffNormErrSquared);
   setMETErrorBars(eeFinal, eeMETToyDistsByBin, eeLowSidebandMETToyDistsByBin, 
   		   eeHighSidebandMETToyDistsByBin, eeNorm, eeNormErrSquared);

   //add EW contribution to QCD control samples
   eeFinal.Add(egMET);
   ffFinal.Add(egMET);

   //make final canvas
   METCanvas.cd();
   makeFinalCanvas(dynamic_cast<TH1*>(&eeFinal), 4, 2, 3005, 4, 0, "E2");
   makeFinalCanvas(dynamic_cast<TH1*>(&ffFinal), kMagenta, 2, 3004, kMagenta, 0, "E2SAME");
   makeFinalCanvas(dynamic_cast<TH1*>(egMET), 8, 2, 3003, 8, 1, "HISTSAME");
   makeFinalCanvas(dynamic_cast<TH1*>(ggMETVsDiEMETVsInvMass.ProjectionZ("ggMET", 0, -1, 0, -1, 
									 "e")), 1, 1, 0, 0, 1, 
		   "SAME");

   //save
   out.cd();
   eeToyCanvas.Write();
   eeLowSidebandToyCanvas.Write();
   eeHighSidebandToyCanvas.Write();
   ffToyCanvas.Write();
   METCanvas.Write();
   out.Write();

   //deallocate memory
   deallocateMemory(eeDiEMETScaled);
   deallocateMemory(eeLowSidebandDiEMETScaled);
   deallocateMemory(eeHighSidebandDiEMETScaled);
   deallocateMemory(ffDiEMETScaled);
   deallocateMemory(eeWeights);
   deallocateMemory(eeLowSidebandWeights);
   deallocateMemory(eeHighSidebandWeights);
   deallocateMemory(ffWeights);
   deallocateMemory(ffFinalToy);
   deallocateMemory(eeFinalToy);
   deallocateMemory(eeLowSidebandFinalToy);
   deallocateMemory(eeHighSidebandFinalToy);
   deallocateMemory(ffMETToyDistsByBin);
   deallocateMemory(eeMETToyDistsByBin);
   deallocateMemory(eeLowSidebandMETToyDistsByBin);
   deallocateMemory(eeHighSidebandMETToyDistsByBin);
   deallocateMemory(ffDiEMETToyDistsByBin);
   deallocateMemory(eeDiEMETToyDistsByBin);
   deallocateMemory(eeLowSidebandDiEMETToyDistsByBin);
   deallocateMemory(eeHighSidebandDiEMETToyDistsByBin);
   out.Close();
}

void GMSBAnalyzer::runMETAnalysisWithEEBackgroundFit(const std::string& outputFile)
{
   if (fChain == 0) return;

   //open file
   TFile out(outputFile.c_str(), "RECREATE");
   out.cd();

   //define constants
   const bool kSumW2 = true;
   const bool kSetGrid = true;
   const bool kDefault = false;
   const unsigned int maxNormBin = 4;     //normalization region
   const float nToys = 1000;              //number of toys for the MET shape error from reweighting
   const unsigned int nMETBins = 13;      //number of MET bins
   // const unsigned int nDiEMETBins = 25;   //number of di-EM ET bins
   const unsigned int nDiEMETBins = 7;   //number of di-EM ET bins
   const unsigned int nInvMassBins = 150; //number of invariant mass bins
   const Double_t METBins[14] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 70.0, 
				 100.0, 150.0, 500.0};            //MET bin boundaries
   // const Double_t diEMETBins[26] = {0.0,                          //di-EM ET bin boundaries
   // 				   3.0, 6.0, 9.0, 12.0, 15.0, 
   // 				   18.0, 21.0, 24.0, 27.0, 30.0, 
   // 				   35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 
   // 				   70.0, 80.0, 90.0, 
   // 				   100.0, 150.0, 200.0, 
   // 				   300.0, 400.0, 500.0};
   const Double_t diEMETBins[8] = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0, 500.0};
   const Double_t invMassBins[151] = {0.0,                        //invariant mass bin boundaries
				      2.0, 4.0, 6.0, 8.0, 10.0, 
				      12.0, 14.0, 16.0, 18.0, 20.0, 
				      22.0, 24.0, 26.0, 28.0, 30.0, 
				      32.0, 34.0, 36.0, 38.0, 40.0, 
				      42.0, 44.0, 46.0, 48.0, 50.0, 
				      52.0, 54.0, 56.0, 58.0, 60.0, 
				      62.0, 64.0, 66.0, 68.0, 70.0, 
				      72.0, 74.0, 76.0, 78.0, 80.0, 
				      82.0, 84.0, 86.0, 88.0, 90.0, 
				      92.0, 94.0, 96.0, 98.0, 100.0, 
				      102.0, 104.0, 106.0, 108.0, 110.0, 
				      112.0, 114.0, 116.0, 118.0, 120.0, 
				      122.0, 124.0, 126.0, 128.0, 130.0, 
				      132.0, 134.0, 136.0, 138.0, 140.0, 
				      142.0, 144.0, 146.0, 148.0, 150.0, 
				      152.0, 154.0, 156.0, 158.0, 160.0, 
				      162.0, 164.0, 166.0, 168.0, 170.0, 
				      172.0, 174.0, 176.0, 178.0, 180.0, 
				      182.0, 184.0, 186.0, 188.0, 190.0, 
				      192.0, 194.0, 196.0, 198.0, 200.0, 
				      202.0, 204.0, 206.0, 208.0, 210.0, 
				      212.0, 214.0, 216.0, 218.0, 220.0, 
				      222.0, 224.0, 226.0, 228.0, 230.0, 
				      232.0, 234.0, 236.0, 238.0, 240.0, 
				      242.0, 244.0, 246.0, 248.0, 250.0, 
				      252.0, 254.0, 256.0, 258.0, 260.0, 
				      262.0, 264.0, 266.0, 268.0, 270.0, 
				      272.0, 274.0, 276.0, 278.0, 280.0, 
				      282.0, 284.0, 286.0, 288.0, 290.0, 
				      292.0, 294.0, 296.0, 298.0, 300.0};
   const float egMisIDRate = 0.014;      /*take from CMS AN-2010/294 for now; can be computed with 
					   Z*/
   const float egMisIDRateErr = 0.002;   /*take from CMS AN-2010/294 for now; can be computed with 
					   Z*/

   //combined isolation histograms (sanity check)
   TH1F ggCombinedIso("ggCombinedIso", "", 40, 0.0, 20.0);
   TH1F egCombinedIso("egCombinedIso", "", 40, 0.0, 20.0);
   TH1F eeCombinedIso("eeCombinedIso", "", 40, 0.0, 20.0);
   TH1F ffCombinedIso("ffCombinedIso", "", 40, 0.0, 20.0);
   setHistogramOptions(ggCombinedIso, "Combined isolation (GeV)", "", "", kSumW2);
   setHistogramOptions(egCombinedIso, "Combined isolation (GeV)", "", "", kSumW2);
   setHistogramOptions(eeCombinedIso, "Combined isolation (GeV)", "", "", kSumW2);
   setHistogramOptions(ffCombinedIso, "Combined isolation (GeV)", "", "", kSumW2);

   //leading photon ET histograms (sanity check)
   TH1F ggLeadingPhotonET("ggLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F egLeadingPhotonET("egLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F eeLeadingPhotonET("eeLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F ffLeadingPhotonET("ffLeadingPhotonET", "", 150, 0.0, 300.0);
   setHistogramOptions(ggLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(egLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);

   //nPV histograms (sanity check)
   TH1F ggNPV("ggNPV", "", 36, -0.5, 35.5);
   TH1F egNPV("egNPV", "", 36, -0.5, 35.5);
   TH1F eeNPV("eeNPV", "", 36, -0.5, 35.5);
   TH1F ffNPV("ffNPV", "", 36, -0.5, 35.5);
   setHistogramOptions(ggNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(egNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(eeNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(ffNPV, "n_{PV}", "", "", kSumW2);

   //MET vs. di-EM ET vs. invariant mass histograms
   TH3F ggMETVsDiEMETVsInvMass("ggMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F egMETVsDiEMETVsInvMass("egMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F eeMETVsDiEMETVsInvMass("eeMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F ffMETVsDiEMETVsInvMass("ffMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   setHistogramOptions(ggMETVsDiEMETVsInvMass, "m_{#gamma#gamma} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(egMETVsDiEMETVsInvMass, "m_{e#gamma} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(eeMETVsDiEMETVsInvMass, "m_{ee} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(ffMETVsDiEMETVsInvMass, "m_{ff} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);

   //ee background subtracted MET vs. di-EM ET histogram
   TH2F eeBkgSubtractedMETVsDiEMET("eeBkgSubtractedMETVsDiEMET", "", nDiEMETBins, diEMETBins, 
				   nMETBins, METBins);
   setHistogramOptions(eeBkgSubtractedMETVsDiEMET, "Di-EM E_{T} (GeV)", "ME_{T} (GeV)", "", kSumW2);

   //reweighted and normalized ee and ff MET histograms
   TH1F eeFinal("eeFinal", "", nMETBins, METBins);
   TH1F ffFinal("ffFinal", "", nMETBins, METBins);
   setHistogramOptions(eeFinal, "ME_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffFinal, "ME_{T} (GeV)", "", "", kSumW2);

   //ee and ff HT vs. MET histograms
   TH2F eeHTVsMET("eeHTVsMET", "", nMETBins, METBins, 60, 0.0, 600.0);
   TH2F ffHTVsMET("ffHTVsMET", "", nMETBins, METBins, 60, 0.0, 600.0);
   setHistogramOptions(eeHTVsMET, "ME_{T} (GeV)", "H_{T} (GeV)", "", kSumW2);
   setHistogramOptions(ffHTVsMET, "ME_{T} (GeV)", "H_{T} (GeV)", "", kSumW2);

   //ee and ff MET vs. no. jets histograms
   TH2F eeMETVsNJets30("eeMETVsNJets30", "", 6, -0.5, 5.5, nMETBins, METBins);
   TH2F ffMETVsNJets30("ffMETVsNJets30", "", 6, -0.5, 5.5, nMETBins, METBins);
   TH2F eeMETVsNJets60("eeMETVsNJets60", "", 6, -0.5, 5.5, nMETBins, METBins);
   TH2F ffMETVsNJets60("ffMETVsNJets60", "", 6, -0.5, 5.5, nMETBins, METBins);
   setHistogramOptions(eeMETVsNJets30, "n_{j} (30-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
   setHistogramOptions(ffMETVsNJets30, "n_{j} (30-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
   setHistogramOptions(eeMETVsNJets60, "n_{j} (60-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
   setHistogramOptions(ffMETVsNJets60, "n_{j} (60-60 GeV)", "ME_{T} (GeV)", "", kSumW2);

   /*histogram of the percentage difference between Poisson mean of generated di-EM ET weights and 
     input mean (i.e. the value of the measured weight)*/
   TH1F ffDiEMETInputVsOutput("ffDiEMETInputVsOutput", "", nDiEMETBins, diEMETBins);
   TH1F eeDiEMETInputVsOutput("eeDiEMETInputVsOutput", "", nDiEMETBins, diEMETBins);
   setHistogramOptions(ffDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);
   setHistogramOptions(eeDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);

   //canvases for the MET toy plots
   TCanvas ffToyCanvas("ffToyCanvas", "", 600, 600);
   TCanvas eeToyCanvas("eeToyCanvas", "", 600, 600);
   setCanvasOptions(ffToyCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);
   setCanvasOptions(eeToyCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);

   //canvas for the final plot
   TCanvas METCanvas("METCanvas", "", 600, 600);
   setCanvasOptions(METCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);

   //MET and di-EM ET vectors
   VFLOAT eeMETVec;
   VFLOAT ffMETVec;
   VFLOAT eeDiEMETVec;
   VFLOAT ffDiEMETVec;

   //PU rho
   TH1F rhoSkim("rhoSkim", "", 100, 0.0, 10.0);
   TH1F rhoPreselected("rhoPreselected", "", 100, 0.0, 10.0);
   setHistogramOptions(rhoSkim, "#rho (GeV)", "", "", kSumW2);
   setHistogramOptions(rhoPreselected, "#rho (GeV)", "", "", kSumW2);

   //set up on-the-fly jet corrections for PF jets
   vector<JetCorrectorParameters> PFJECs;
   PFJECs.push_back(JetCorrectorParameters(L1JECFile_));
   PFJECs.push_back(JetCorrectorParameters(L2JECFile_));
   PFJECs.push_back(JetCorrectorParameters(L3JECFile_));
   FactorizedJetCorrector PFJetCorrector(PFJECs);

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

     //get int. lumi weight for this dataset
     TString fileName(fChain->GetCurrentFile()->GetName());
     TString dataset;
     if (fileName.Contains("ntuple_.*-v[0-9]")) dataset = fileName("ntuple_.*-v[0-9]");
     if (dataset.Length() >= 7) dataset.Remove(0, 7);
     string datasetString((const char*)dataset);
     map<string, unsigned int>::const_iterator iDataset = fileMap_.find(datasetString);
     float lumiWeight = 1.0;
     if (iDataset != fileMap_.end()) lumiWeight = weight_[iDataset->second];
     // else {
     //   /*once comfortable with MC dataset names, remove this and just assume the file is data and 
     // 	 gets weight 1.0*/
     //   cerr << "Error: dataset " << datasetString << " was not entered into the map.\n";
     // }

     //get PU weight for this dataset
     /*in-time reweighting only following 
       https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities as of 27-Oct-11 
       and Tessa's recommendation*/
     // int nPV  = -1;
     float PUWeight = 1.0;
     // susy::PUSummaryInfoCollection::const_iterator iBX = susyEvent->PU.begin();
     // bool foundInTimeBX = false;
     // while ((iBX != susyEvent->PU.end()) && !foundInTimeBX) {
     //   if (iBX->BX == 0) { 
     // 	 nPV = iBX->numInteractions;
     // 	 foundInTimeBX = true;
     //   }
     //   ++iBX;
     // }
     // PUWeight = lumiWeights_.ITweight(nPV);
     // PUWeight = 1.0;

     //fill MET vs. di-EM ET vs. invariant mass, rho, leading photon ET, and nPV histograms
     const double invMass = susyCategory->getEvtInvMass(tag_);
     double MET = -1.0;
     map<TString, susy::MET>::const_iterator iMET = susyEvent->metMap.find("pfMet");
     if (iMET != susyEvent->metMap.end()) MET = iMET->second.met();
     else cerr << "Error: PFMET not found in event " << (jentry + 1) << ".\n";
     const double diEMET = susyCategory->getEvtDiEMET(tag_);
     const Float_t rho = susyEvent->rho;
     rhoSkim.Fill(rho, lumiWeight*PUWeight);
     const int evtCategory = susyCategory->getEventCategory(tag_);
     if (evtCategory != FAIL) rhoPreselected.Fill(rho, lumiWeight*PUWeight);
     unsigned int nGoodRecoPV = 0;
     for (vector<susy::Vertex>::const_iterator iPV = susyEvent->vertices.begin(); 
	  iPV != susyEvent->vertices.end(); ++iPV) {
       if (!iPV->isFake() && (iPV->ndof > 4) && (iPV->position.Z() <= 24.0/*cm*/) && 
	   (iPV->position.Perp() <= 2.0/*cm*/)) ++nGoodRecoPV;
     }

     //count jets
     unsigned int nJets = 0;
     unsigned int nJets30 = 0;
     unsigned int nJets60 = 0;
     vector<float> jetET;
     map<TString, susy::PhotonCollection>::const_iterator iPhotonMap = 
       susyEvent->photons.find((const char*)tag_);
     map<TString, susy::PFJetCollection>::const_iterator iJets = susyEvent->pfJets.find("ak5");
     if ((iJets != susyEvent->pfJets.end()) && (iPhotonMap != susyEvent->photons.end())) {

       //loop over PF jets
       for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
	    iJet != iJets->second.end(); ++iJet) {

	 //compute corrected P4
	 PFJetCorrector.setJetEta(iJet->momentum.Eta());
	 PFJetCorrector.setJetPt(iJet->momentum.Pt());
	 PFJetCorrector.setJetA(iJet->jetArea);
	 PFJetCorrector.setRho(susyEvent->rho);
	 float storedCorr = -1.0;
	 map<TString, Float_t>::const_iterator iCorr = iJet->jecScaleFactors.find("L1FastL2L3");
	 if (iCorr != iJet->jecScaleFactors.end()) storedCorr = iCorr->second;
	 TLorentzVector corrP4 = storedCorr*iJet->momentum;

	 //calculate the jet-EM overlap
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
	 }
	 if ((count != 2) && (evtCategory != FAIL)) {
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
	     ++nJets;
	     jetET.push_back(corrP4.Et());
	     if ((corrP4.Et() >= 30.0/*GeV*/) && (corrP4.Et() < 60.0/*GeV*/)) ++nJets30;
	     else if (corrP4.Et() >= 60.0/*GeV*/) ++nJets60;
	   }
	 }
       }
     }
     else {
       cerr << "Error: " << tag_ << " photon collection or ak5 jet collection not found in ";
       cerr << "event " << (jentry + 1) << ".\n";
     }
     const float HT = accumulate(jetET.begin(), jetET.end(), 0);

     //get combined isolation and photon ET
     vector<float> combinedIso;
     vector<float> ET;
     if (evtCategory != FAIL) {
       if (iPhotonMap != susyEvent->photons.end()) {
	 for (susy::PhotonCollection::const_iterator iPhoton = iPhotonMap->second.begin(); 
	      iPhoton != iPhotonMap->second.end(); ++iPhoton) {
	   if (susyCategory->getIsDeciding(tag_, iPhoton - iPhotonMap->second.begin())) {
	     ET.push_back(iPhoton->momentum.Et());
	     combinedIso.push_back(iPhoton->ecalRecHitSumEtConeDR03 - 0.1474*rho + 
				   iPhoton->hcalTowerSumEtConeDR03() - 0.0467*rho + 
				   iPhoton->trkSumPtHollowConeDR03);
	   }
	 }
	 if ((combinedIso.size() != 2) || (ET.size() != 2)) {
	   cerr << "Error: combinedIso.size() = " << combinedIso.size() << " and ET.size() = ";
	   cerr << ET.size() << " in event " << (jentry + 1) << ".\n";
	 }
       }
       else {
	 cerr << "Error: " << tag_ << " photon collection not found in event " << (jentry + 1);
	 cerr << ".\n";
       }
     }

     //sort photon ET in ascending order
     sort(ET.begin(), ET.end());

//      if (nJets == 0) {
     bool passTightCombinedIso = true;
       switch (evtCategory) {
       case GG:
	 ggMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	 for (vector<float>::const_iterator iIso = combinedIso.begin(); 
	      iIso != combinedIso.end(); ++iIso) {
	   ggCombinedIso.Fill(*iIso, lumiWeight*PUWeight);
	 }
	 ggLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 ggNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 break;
       case EG:
	 egMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	 for (vector<float>::const_iterator iIso = combinedIso.begin(); 
	      iIso != combinedIso.end(); ++iIso) {
	   egCombinedIso.Fill(*iIso, lumiWeight*PUWeight);
	 }
	 egLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 egNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 break;
       case EE:
	 eeMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	 eeMETVec.push_back(MET);
	 eeDiEMETVec.push_back(diEMET);
	 eeHTVsMET.Fill(MET, HT, lumiWeight*PUWeight);
	 eeMETVsNJets30.Fill(nJets30, MET, lumiWeight*PUWeight);
	 eeMETVsNJets60.Fill(nJets60, MET, lumiWeight*PUWeight);
	 for (vector<float>::const_iterator iIso = combinedIso.begin(); 
	      iIso != combinedIso.end(); ++iIso) {
	   eeCombinedIso.Fill(*iIso, lumiWeight*PUWeight);
	 }
	 eeLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 eeNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 break;
       case FF:
	 for (vector<float>::const_iterator iIso = combinedIso.begin(); 
	      iIso != combinedIso.end(); ++iIso) {
	   ffCombinedIso.Fill(*iIso, lumiWeight*PUWeight);
	   passTightCombinedIso = passTightCombinedIso && (*iIso < 9.0/*GeV*/);
	 }
	 ffMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	 ffMETVec.push_back(MET);
	 ffDiEMETVec.push_back(diEMET);
	 ffLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 ffNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 ffHTVsMET.Fill(MET, HT, lumiWeight*PUWeight);
	 ffMETVsNJets30.Fill(nJets30, MET, lumiWeight*PUWeight);
	 ffMETVsNJets60.Fill(nJets60, MET, lumiWeight*PUWeight);
	 break;
       default:
	 break;
       }
//      }
   }

   //perform ee sideband subtraction
   generateBackgroundSubtractedSpectra(eeMETVsDiEMETVsInvMass, eeBkgSubtractedMETVsDiEMET);

   //fill weights histograms
   vector<TH1F*> eeDiEMETScaled;
   vector<TH1F*> eeWeights;
   vector<TH1F*> ffDiEMETScaled;
   vector<TH1F*> ffWeights;
   for (unsigned int iMETBin = 1; iMETBin <= nMETBins; ++iMETBin) {
     TH1F* eeHistDiEMETScaled = new TH1F(histName("ee", "DiEMETScaled_METBin", iMETBin).c_str(), 
					 "", nDiEMETBins, diEMETBins);
     TH1F* eeHistWeights = new TH1F(histName("ee", "Weights_METBin", iMETBin).c_str(), "", 
				    nDiEMETBins, diEMETBins);
     setHistogramOptions(*eeHistDiEMETScaled, "Di-EM E_{T} (GeV)", "", "", kSumW2);
     setHistogramOptions(*eeHistWeights, "Di-EM E_{T} (GeV)", "", "", kSumW2);
     fillWeightsHistograms(ggMETVsDiEMETVsInvMass.ProjectionY(histName("ggDiEMET_METBin", "", 
								   iMETBin).c_str(), 
							      0, -1, iMETBin, iMETBin, "e"), 
			   eeBkgSubtractedMETVsDiEMET.ProjectionX(histName("eeDiEMET_METBin", "", 
								       iMETBin).c_str(), 
								  iMETBin, iMETBin, "e"), 
			   *eeHistDiEMETScaled, *eeHistWeights);
     eeDiEMETScaled.push_back(eeHistDiEMETScaled);
     eeWeights.push_back(eeHistWeights);
     if (iMETBin == 1) {
       TH1F* ffHistDiEMETScaled = new TH1F(histName("ff", "DiEMETScaled_METBin", iMETBin).c_str(), 
					   "", nDiEMETBins, diEMETBins);
       TH1F* ffHistWeights = new TH1F(histName("ff", "Weights_METBin", iMETBin).c_str(), "", 
				      nDiEMETBins, diEMETBins);
       setHistogramOptions(*ffHistDiEMETScaled, "Di-EM E_{T} (GeV)", "", "", kSumW2);
       setHistogramOptions(*ffHistWeights, "Di-EM E_{T} (GeV)", "", "", kSumW2);
       fillWeightsHistograms(ggMETVsDiEMETVsInvMass.ProjectionY(histName("ggDiEMET_METBin", "", 
								   iMETBin).c_str(), 
								0, -1, 0, -1, "e"), 
			     ffMETVsDiEMETVsInvMass.ProjectionY(histName("ffDiEMET_METBin", "", 
									 iMETBin).c_str(), 0, -1, 
								0, -1, "e"), 
			     *ffHistDiEMETScaled, *ffHistWeights);
       ffDiEMETScaled.push_back(ffHistDiEMETScaled);
       ffWeights.push_back(ffHistWeights);
     }
   }

   //generate toys for calculating error due to MET shape from reweighting
   vector<TH1F*> ffFinalToy;
   vector<TH1F*> eeFinalToy;
   vector<TH1F*> ffDiEMETToyDistsByBin;
   vector<TH1F*> eeDiEMETToyDistsByBin;
   generateToys(ffFinalToy, ffDiEMETToyDistsByBin, ffWeights, nToys, "ff", diEMETBins, nMETBins, 
		METBins, ffMETVec, ffDiEMETVec);
   generateToys(eeFinalToy, eeDiEMETToyDistsByBin, eeWeights, nToys, "ee", diEMETBins, nMETBins, 
		METBins, eeBkgSubtractedMETVsDiEMET);

   //reweight ee and ff MET histograms
   for (unsigned int iMETBin = 1; iMETBin <= nMETBins; ++iMETBin) {
     reweightBinned(eeBkgSubtractedMETVsDiEMET, *eeWeights[iMETBin - 1], &eeFinal, iMETBin);
   }
   reweightDefault(ffMETVec, ffDiEMETVec, *ffWeights[0], &ffFinal);

   //make individual histograms of MET toy distributions, 1 per MET bin
   vector<TH1F*> ffMETToyDistsByBin;
   vector<TH1F*> eeMETToyDistsByBin;
   fillToyDistributions(ffMETToyDistsByBin, ffFinal, ffToyCanvas, ffFinalToy, nToys, "ff", 
   			nMETBins);
   fillToyDistributions(eeMETToyDistsByBin, eeFinal, eeToyCanvas, eeFinalToy, nToys, "ee", 
   			nMETBins);

   //calculate EW contribution
   TH1D* egMET = egMETVsDiEMETVsInvMass.ProjectionZ("egMET", 0, -1, 0, -1, "e");
   const float egScale = egMisIDRate/(1.0 - egMisIDRate);
   egMET->Scale(egScale);
   for (Int_t iBin = 1; iBin <= egMET->GetNbinsX(); ++iBin) {
     const float binContent = egMET->GetBinContent(iBin);
     const float binPoissonError = egMET->GetBinError(iBin);
     if (binContent == 0.0) egMET->SetBinError(iBin, 0.0);
     else {
       egMET->SetBinError(iBin, 
			  egScale*sqrt(((binPoissonError*binPoissonError)/
					(binContent*binContent)) + 
				       ((egMisIDRateErr*egMisIDRateErr)/
					(egMisIDRate*egMisIDRate))));
     }
   }

   //normalize ee and ff MET histograms
   float eeNormErrSquared = 0.0;
   float ffNormErrSquared = 0.0;
   const float eeNorm = 
     normAndErrorSquared(ggMETVsDiEMETVsInvMass, eeFinal, egMET, maxNormBin, eeNormErrSquared);
   const float ffNorm = 
     normAndErrorSquared(ggMETVsDiEMETVsInvMass, ffFinal, egMET, maxNormBin, ffNormErrSquared);
   eeFinal.Scale(eeNorm);
   ffFinal.Scale(ffNorm);

   //set MET error bars
   setMETErrorBars(ffFinal, ffMETToyDistsByBin, vector<TH1F*>(), vector<TH1F*>(), ffNorm, 
		   ffNormErrSquared);
   setMETErrorBars(eeFinal, eeMETToyDistsByBin, vector<TH1F*>(), vector<TH1F*>(), eeNorm, 
		   eeNormErrSquared, kDefault);

   //add EW contribution to QCD control samples
   eeFinal.Add(egMET);
   ffFinal.Add(egMET);

   //make final canvas
   METCanvas.cd();
   makeFinalCanvas(dynamic_cast<TH1*>(&eeFinal), 4, 2, 3005, 4, 0, "E2");
   makeFinalCanvas(dynamic_cast<TH1*>(&ffFinal), kMagenta, 2, 3004, kMagenta, 0, "E2SAME");
   makeFinalCanvas(dynamic_cast<TH1*>(egMET), 8, 2, 3003, 8, 1, "HISTSAME");
   makeFinalCanvas(dynamic_cast<TH1*>(ggMETVsDiEMETVsInvMass.ProjectionZ("ggMET", 0, -1, 0, -1, 
									 "e")), 1, 1, 0, 0, 1, 
		   "SAME");

   //save
   out.cd();
   eeToyCanvas.Write();
   ffToyCanvas.Write();
   METCanvas.Write();
   out.Write();

   //deallocate memory
   deallocateMemory(eeDiEMETScaled);
   deallocateMemory(eeWeights);
   deallocateMemory(ffDiEMETScaled);
   deallocateMemory(ffWeights);
   deallocateMemory(ffFinalToy);
   deallocateMemory(eeFinalToy);
   deallocateMemory(ffMETToyDistsByBin);
   deallocateMemory(eeMETToyDistsByBin);
   deallocateMemory(ffDiEMETToyDistsByBin);
   deallocateMemory(eeDiEMETToyDistsByBin);
   out.Close();
}

void GMSBAnalyzer::testFitting(const string& inputFile, const string& outputFile) const
{
  //define constants
  const bool kSumW2 = true;
  const bool kSetGrid = true;
  const bool kDefault = false;
  const unsigned int maxNormBin = 4;     //normalization region
  const float nToys = 1000;              //number of toys for the MET shape error from reweighting
  const unsigned int nMETBins = 13;      //number of MET bins
  // const unsigned int nDiEMETBins = 25;   //number of di-EM ET bins
  const unsigned int nDiEMETBins = 7;   //number of di-EM ET bins
  // const unsigned int nInvMassBins = 150; //number of invariant mass bins
  const Double_t METBins[14] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 70.0, 
				100.0, 150.0, 500.0};            //MET bin boundaries
  // const Double_t diEMETBins[26] = {0.0,                          //di-EM ET bin boundaries
  // 				   3.0, 6.0, 9.0, 12.0, 15.0, 
  // 				   18.0, 21.0, 24.0, 27.0, 30.0, 
  // 				   35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 
  // 				   70.0, 80.0, 90.0, 
  // 				   100.0, 150.0, 200.0, 
  // 				   300.0, 400.0, 500.0};
  const Double_t diEMETBins[8] = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0, 500.0};
  // const Double_t invMassBins[151] = {0.0,                        //invariant mass bin boundaries
  // 				     2.0, 4.0, 6.0, 8.0, 10.0, 
  // 				     12.0, 14.0, 16.0, 18.0, 20.0, 
  // 				     22.0, 24.0, 26.0, 28.0, 30.0, 
  // 				     32.0, 34.0, 36.0, 38.0, 40.0, 
  // 				     42.0, 44.0, 46.0, 48.0, 50.0, 
  // 				     52.0, 54.0, 56.0, 58.0, 60.0, 
  // 				     62.0, 64.0, 66.0, 68.0, 70.0, 
  // 				     72.0, 74.0, 76.0, 78.0, 80.0, 
  // 				     82.0, 84.0, 86.0, 88.0, 90.0, 
  // 				     92.0, 94.0, 96.0, 98.0, 100.0, 
  // 				     102.0, 104.0, 106.0, 108.0, 110.0, 
  // 				     112.0, 114.0, 116.0, 118.0, 120.0, 
  // 				     122.0, 124.0, 126.0, 128.0, 130.0, 
  // 				     132.0, 134.0, 136.0, 138.0, 140.0, 
  // 				     142.0, 144.0, 146.0, 148.0, 150.0, 
  // 				     152.0, 154.0, 156.0, 158.0, 160.0, 
  // 				     162.0, 164.0, 166.0, 168.0, 170.0, 
  // 				     172.0, 174.0, 176.0, 178.0, 180.0, 
  // 				     182.0, 184.0, 186.0, 188.0, 190.0, 
  // 				     192.0, 194.0, 196.0, 198.0, 200.0, 
  // 				     202.0, 204.0, 206.0, 208.0, 210.0, 
  // 				     212.0, 214.0, 216.0, 218.0, 220.0, 
  // 				     222.0, 224.0, 226.0, 228.0, 230.0, 
  // 				     232.0, 234.0, 236.0, 238.0, 240.0, 
  // 				     242.0, 244.0, 246.0, 248.0, 250.0, 
  // 				     252.0, 254.0, 256.0, 258.0, 260.0, 
  // 				     262.0, 264.0, 266.0, 268.0, 270.0, 
  // 				     272.0, 274.0, 276.0, 278.0, 280.0, 
  // 				     282.0, 284.0, 286.0, 288.0, 290.0, 
  // 				     292.0, 294.0, 296.0, 298.0, 300.0};
  const float egMisIDRate = 0.014;      /*take from CMS AN-2010/294 for now; can be computed with 
					  Z*/
  const float egMisIDRateErr = 0.002;   /*take from CMS AN-2010/294 for now; can be computed with 
					  Z*/

  TH2F eeBkgSubtractedMETVsDiEMET("eeBkgSubtractedMETVsDiEMET", "", nDiEMETBins, diEMETBins, 
				  nMETBins, METBins);
  setHistogramOptions(eeBkgSubtractedMETVsDiEMET, "", "Di-EM E_{T} (GeV)", "ME_{T} (GeV)", kSumW2);

  //reweighted and normalized ee and ff MET histograms
  TH1F eeFinal("eeFinal", "", nMETBins, METBins);
  TH1F ffFinal("ffFinal", "", nMETBins, METBins);
  setHistogramOptions(eeFinal, "ME_{T} (GeV)", "", "", kSumW2);
  setHistogramOptions(ffFinal, "ME_{T} (GeV)", "", "", kSumW2);

  //ee and ff HT vs. MET histograms
  TH2F eeHTVsMET("eeHTVsMET", "", nMETBins, METBins, 60, 0.0, 600.0);
  TH2F ffHTVsMET("ffHTVsMET", "", nMETBins, METBins, 60, 0.0, 600.0);
  setHistogramOptions(eeHTVsMET, "ME_{T} (GeV)", "H_{T} (GeV)", "", kSumW2);
  setHistogramOptions(ffHTVsMET, "ME_{T} (GeV)", "H_{T} (GeV)", "", kSumW2);

  //ee and ff MET vs. no. jets histograms
  TH2F eeMETVsNJets30("eeMETVsNJets30", "", 6, -0.5, 5.5, nMETBins, METBins);
  TH2F ffMETVsNJets30("ffMETVsNJets30", "", 6, -0.5, 5.5, nMETBins, METBins);
  TH2F eeMETVsNJets60("eeMETVsNJets60", "", 6, -0.5, 5.5, nMETBins, METBins);
  TH2F ffMETVsNJets60("ffMETVsNJets60", "", 6, -0.5, 5.5, nMETBins, METBins);
  setHistogramOptions(eeMETVsNJets30, "n_{j} (30-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
  setHistogramOptions(ffMETVsNJets30, "n_{j} (30-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
  setHistogramOptions(eeMETVsNJets60, "n_{j} (60-60 GeV)", "ME_{T} (GeV)", "", kSumW2);
  setHistogramOptions(ffMETVsNJets60, "n_{j} (60-60 GeV)", "ME_{T} (GeV)", "", kSumW2);

  /*histogram of the percentage difference between Poisson mean of generated di-EM ET weights and 
    input mean (i.e. the value of the measured weight)*/
  TH1F ffDiEMETInputVsOutput("ffDiEMETInputVsOutput", "", nDiEMETBins, diEMETBins);
  TH1F eeDiEMETInputVsOutput("eeDiEMETInputVsOutput", "", nDiEMETBins, diEMETBins);
  setHistogramOptions(ffDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
		      "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);
  setHistogramOptions(eeDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
		      "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);

  //canvases for the MET toy plots
  TCanvas ffToyCanvas("ffToyCanvas", "", 600, 600);
  TCanvas eeToyCanvas("eeToyCanvas", "", 600, 600);
  setCanvasOptions(ffToyCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);
  setCanvasOptions(eeToyCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);

  //canvas for the final plot
  TCanvas METCanvas("METCanvas", "", 600, 600);
  setCanvasOptions(METCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);

  //get TH3Fs
  TFile in(inputFile.c_str());
  TFile out(outputFile.c_str(), "RECREATE");
  TH3F* ggMETVsDiEMETVsInvMass;
  TH3F* egMETVsDiEMETVsInvMass;
  TH3F* eeMETVsDiEMETVsInvMass;
  TH3F* ffMETVsDiEMETVsInvMass;
  in.cd();
  in.GetObject("ggMETVsDiEMETVsInvMass", ggMETVsDiEMETVsInvMass);
  in.GetObject("egMETVsDiEMETVsInvMass", egMETVsDiEMETVsInvMass);
  in.GetObject("eeMETVsDiEMETVsInvMass", eeMETVsDiEMETVsInvMass);
  in.GetObject("ffMETVsDiEMETVsInvMass", ffMETVsDiEMETVsInvMass);
  out.cd();

  //perform ee sideband subtraction
  generateBackgroundSubtractedSpectra(*eeMETVsDiEMETVsInvMass, eeBkgSubtractedMETVsDiEMET);

  //fill weights histograms
  vector<TH1F*> eeDiEMETScaled;
  vector<TH1F*> eeWeights;
  vector<TH1F*> ffDiEMETScaled;
  vector<TH1F*> ffWeights;
  for (unsigned int iMETBin = 1; iMETBin <= nMETBins; ++iMETBin) {
    TH1F* eeHistDiEMETScaled = new TH1F(histName("ee", "DiEMETScaled_METBin", iMETBin).c_str(), 
					"", nDiEMETBins, diEMETBins);
    TH1F* eeHistWeights = new TH1F(histName("ee", "Weights_METBin", iMETBin).c_str(), "", 
				   nDiEMETBins, diEMETBins);
    setHistogramOptions(*eeHistDiEMETScaled, "Di-EM E_{T} (GeV)", "", "", kSumW2);
    setHistogramOptions(*eeHistWeights, "Di-EM E_{T} (GeV)", "", "", kSumW2);
    fillWeightsHistograms(ggMETVsDiEMETVsInvMass->ProjectionY(histName("ggDiEMET_METBin", "", 
								       iMETBin).c_str(), 
							      0, -1, iMETBin, iMETBin, "e"), 
			  eeBkgSubtractedMETVsDiEMET.ProjectionX(histName("eeDiEMET_METBin", "", 
									  iMETBin).c_str(), 
								 iMETBin, iMETBin, "e"), 
			  *eeHistDiEMETScaled, *eeHistWeights);
    eeDiEMETScaled.push_back(eeHistDiEMETScaled);
    eeWeights.push_back(eeHistWeights);
    if (iMETBin == 1) {
      TH1F* ffHistDiEMETScaled = new TH1F(histName("ff", "DiEMETScaled_METBin", iMETBin).c_str(), 
					  "", nDiEMETBins, diEMETBins);
      TH1F* ffHistWeights = new TH1F(histName("ff", "Weights_METBin", iMETBin).c_str(), "", 
				     nDiEMETBins, diEMETBins);
      setHistogramOptions(*ffHistDiEMETScaled, "Di-EM E_{T} (GeV)", "", "", kSumW2);
      setHistogramOptions(*ffHistWeights, "Di-EM E_{T} (GeV)", "", "", kSumW2);
      fillWeightsHistograms(ggMETVsDiEMETVsInvMass->ProjectionY(histName("ggDiEMET_METBin", "", 
									 iMETBin).c_str(), 
								0, -1, 0, -1, "e"), 
			    ffMETVsDiEMETVsInvMass->ProjectionY(histName("ffDiEMET_METBin", "", 
									 iMETBin).c_str(), 
								0, -1, 0, -1, "e"), 
			    *ffHistDiEMETScaled, *ffHistWeights);
      ffDiEMETScaled.push_back(ffHistDiEMETScaled);
      ffWeights.push_back(ffHistWeights);
    }
  }

  //generate toys for calculating error due to MET shape from reweighting
  vector<TH1F*> ffFinalToy;
  vector<TH1F*> eeFinalToy;
  vector<TH1F*> ffDiEMETToyDistsByBin;
  vector<TH1F*> eeDiEMETToyDistsByBin;
  // generateToys(ffFinalToy, ffDiEMETToyDistsByBin, ffWeights, nToys, "ff", diEMETBins, nMETBins, 
  // 	       METBins, ffMETVec, ffDiEMETVec);
  generateToys(eeFinalToy, eeDiEMETToyDistsByBin, eeWeights, nToys, "ee", diEMETBins, nMETBins, 
	       METBins, eeBkgSubtractedMETVsDiEMET);

  //reweight ee and ff MET histograms
  for (unsigned int iMETBin = 1; iMETBin <= nMETBins; ++iMETBin) {
    reweightBinned(eeBkgSubtractedMETVsDiEMET, *eeWeights[iMETBin - 1], &eeFinal, iMETBin);
  }
  // reweightDefault(ffMETVec, ffDiEMETVec, *ffWeights[0], &ffFinal);


  //make individual histograms of MET toy distributions, 1 per MET bin
  vector<TH1F*> ffMETToyDistsByBin;
  vector<TH1F*> eeMETToyDistsByBin;
  // fillToyDistributions(ffMETToyDistsByBin, ffFinal, ffToyCanvas, ffFinalToy, nToys, "ff", 
  // 		       nMETBins);
  fillToyDistributions(eeMETToyDistsByBin, eeFinal, eeToyCanvas, eeFinalToy, nToys, "ee", 
		       nMETBins);

  //calculate EW contribution
  TH1D* egMET = egMETVsDiEMETVsInvMass->ProjectionZ("egMET", 0, -1, 0, -1, "e");
  const float egScale = egMisIDRate/(1.0 - egMisIDRate);
  egMET->Scale(egScale);
  for (Int_t iBin = 1; iBin <= egMET->GetNbinsX(); ++iBin) {
    const float binContent = egMET->GetBinContent(iBin);
    const float binPoissonError = egMET->GetBinError(iBin);
    if (binContent == 0.0) egMET->SetBinError(iBin, 0.0);
    else {
      egMET->SetBinError(iBin, 
			 egScale*sqrt(((binPoissonError*binPoissonError)/
				       (binContent*binContent)) + 
				      ((egMisIDRateErr*egMisIDRateErr)/
				       (egMisIDRate*egMisIDRate))));
    }
  }

  //normalize ee and ff MET histograms
  float eeNormErrSquared = 0.0;
  // float ffNormErrSquared = 0.0;
  const float eeNorm = 
    normAndErrorSquared(*ggMETVsDiEMETVsInvMass, eeFinal, egMET, maxNormBin, eeNormErrSquared);
  // const float ffNorm = 
  //   normAndErrorSquared(*ggMETVsDiEMETVsInvMass, ffFinal, egMET, maxNormBin, ffNormErrSquared);
  eeFinal.Scale(eeNorm);
  // ffFinal.Scale(ffNorm);

  //set MET error bars
  // setMETErrorBars(ffFinal, ffMETToyDistsByBin, vector<TH1F*>(), vector<TH1F*>(), ffNorm, 
  // 		  ffNormErrSquared);
  setMETErrorBars(eeFinal, eeMETToyDistsByBin, vector<TH1F*>(), vector<TH1F*>(), eeNorm, 
		  eeNormErrSquared, kDefault);

  //add EW contribution to QCD control samples
  eeFinal.Add(egMET);
  // ffFinal.Add(egMET);

  //make final canvas
  METCanvas.cd();
  makeFinalCanvas(dynamic_cast<TH1*>(&eeFinal), 4, 2, 3005, 4, 0, "E2");
  // makeFinalCanvas(dynamic_cast<TH1*>(&ffFinal), kMagenta, 2, 3004, kMagenta, 0, "E2SAME");
  makeFinalCanvas(dynamic_cast<TH1*>(egMET), 8, 2, 3003, 8, 1, "HISTSAME");
  makeFinalCanvas(dynamic_cast<TH1*>(ggMETVsDiEMETVsInvMass->ProjectionZ("ggMET", 0, -1, 0, -1, 
									 "e")), 1, 1, 0, 0, 1, 
		  "SAME");

  //save
  out.cd();
  eeToyCanvas.Write();
  ffToyCanvas.Write();
  METCanvas.Write();
  out.Write();

  //deallocate memory
  deallocateMemory(eeDiEMETScaled);
  cout << "Line " << __LINE__ << endl;
  deallocateMemory(eeWeights);
  cout << "Line " << __LINE__ << endl;
  deallocateMemory(ffDiEMETScaled);
  cout << "Line " << __LINE__ << endl;
  deallocateMemory(ffWeights);
  cout << "Line " << __LINE__ << endl;
  deallocateMemory(ffFinalToy);
  cout << "Line " << __LINE__ << endl;
  // deallocateMemory(eeFinalToy);
  cout << "Line " << __LINE__ << endl;
  deallocateMemory(ffMETToyDistsByBin);
  cout << "Line " << __LINE__ << endl;
  deallocateMemory(eeMETToyDistsByBin);
  cout << "Line " << __LINE__ << endl;
  deallocateMemory(ffDiEMETToyDistsByBin);
  cout << "Line " << __LINE__ << endl;
  deallocateMemory(eeDiEMETToyDistsByBin);
  cout << "Line " << __LINE__ << endl;
  out.Close();
  in.Close();
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
   PFJECs.push_back(JetCorrectorParameters(L1JECFile_));
   PFJECs.push_back(JetCorrectorParameters(L2JECFile_));
   PFJECs.push_back(JetCorrectorParameters(L3JECFile_));
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
     double MET = -1.0;
     map<TString, susy::MET>::const_iterator iMET = susyEvent->metMap.find("pfMet");
     if (iMET != susyEvent->metMap.end()) MET = iMET->second.met();
     else cerr << "Error: PFMET not found in event " << (jentry + 1) << ".\n";
     if ((((category == EE) && (invMass >= 81.0/*GeV*/) && (invMass < 101.0/*GeV*/)) || 
	  (category == FF)) && (MET >= 20.0/*GeV*/) && (MET < 50.0/*GeV*/)) {
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
	     cout << "Run " << susyEvent->runNumber << ", event " << susyEvent->eventNumber;
	     cout << ", lumi section " << susyEvent->luminosityBlockNumber << ", jet ";
	     cout << (iJet - iJets->second.begin()) << endl;
	     cout << "PF jet stored correction: " << storedCorr << ", on the fly correction: ";
	     cout << onTheFlyCorr << endl;
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

void GMSBAnalyzer::skim(const string& outputFile, const int evtCategory)
{
  if (fChain == 0) return;

  //open file and clone tree
  TFile out(outputFile.c_str(), "RECREATE");
  out.cd();
  TTree* skimTree = fChain->CloneTree(0);

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

    //skim on event category
    if (susyCategory->getEventCategory(tag_) == evtCategory) skimTree->Fill();

    //hapless attempt to save memory
    susyEvent->Init();
    susyCategory->reset();
  }

  //close files
  out.Write();
  out.Close();
}

void GMSBAnalyzer::stripBranch(const string& outputFile, const string& branch)
{
  if (fChain == 0) return;

  //open file and clone tree
  TFile out(outputFile.c_str(), "RECREATE");
  out.cd();
  TTree* skimTree = fChain->CloneTree(0);

  //set user-specified number of entries to process
  Long64_t nentries = fChain->GetEntriesFast();
  if (nEvts_ != -1) nentries = nEvts_;

  //copy tree to new file
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    TBranch* pBranch = fChain->GetTree()->GetBranch(branch.c_str());
    nb = /*fChain*/pBranch->GetEntry(jentry);   nbytes += nb;
    if (((jentry + 1) % 10000) == 0) cout << "Event " << (jentry + 1) << endl;
    skimTree->Fill();
  }

  //close files
  out.Write();
  out.Close();
}

void GMSBAnalyzer::debugPrint(const unsigned int jentry) const
{
  //loop over photons
  cout << "*******************Event " << (jentry + 1) << "*******************\n\n";
  cout << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
  cout << susyEvent->luminosityBlockNumber << endl;
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
      cout << "--------------";
      cout << (susyCategory->getPassETMin2(tag_, i) ? "PASSED TRAILING ET\n" : "FAILED TRAILING ET\n");
      cout << "-------eta  : " << iPhoton->caloPosition.Eta() << endl;
      cout << "--------------";
      cout << (susyCategory->getPassAbsEtaMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------phi  : " << iPhoton->caloPosition.Phi() << endl;
      cout << "-------IECAL: " << iPhoton->ecalRecHitSumEtConeDR03 << endl;
      cout << "--------------";
      cout << (susyCategory->getPassECALIsoMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------IHCAL: " << iPhoton->hcalTowerSumEtConeDR03() << endl;
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
      cout << "-------ITRK : " << iPhoton->trkSumPtHollowConeDR03 << endl;
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
  cout << "--------------Pass DRMin : ";
  cout << (susyCategory->getPassDRMin(tag_) ? "PASSED\n" : "FAILED\n") << endl;
  cout << "--------------DiEM ET      : " << susyCategory->getEvtDiEMET(tag_) << endl;
  cout << "--------------InvMass      : " << susyCategory->getEvtInvMass(tag_) << endl;
  cout << endl;
}
