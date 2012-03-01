#define GMSBAnalyzer_cxx
#include "GMSBAnalyzer.h"
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
#include "SusyAnalysis/SusyNtuplizer/jec/JetMETObjects/interface/JetCorrectionUncertainty.h"

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
					  TH1F& controlDiEMETScaled, TH1F& controlWeights, 
					  const float scale) const
{
//   float scale = 0.0;
//   if (controlDiEMET->Integral() != 0.0) {
//     scale = ggDiEMET->Integral()/controlDiEMET->Integral();
//   }
//   float scale = 1.0; //scaling undoes Nj-dependent reweighting cf. e-mail with Yueh-Feng 2-Feb-12
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

void GMSBAnalyzer::format(const char* col1Title, const char* col2Title, const char* col3Title, 
			  const char* col4Title, const char* col5Title) const
{
  cout.width(20);
  cout.fill(' ');
  cout << left << col1Title;
  cout.width(20);
  cout << left << col2Title;
  cout.width(20);
  cout << left << col3Title;
  cout.width(20);
  cout << left << col4Title;
  cout.width(20);
  cout << left << col5Title;
  cout << endl;
}

void GMSBAnalyzer::format(const char* bin, const float weight, const float statErr, 
			  const float JESErr) const
{
  cout.width(20);
  cout.fill(' ');
  cout << left << bin;
  cout.width(20);
  cout << left << weight;
  cout.width(20);
  cout << left << statErr;
  cout.width(20);
  cout << left << JESErr;
  cout.width(20);
  cout << left << sqrt(statErr*statErr + JESErr*JESErr);
  cout << endl;
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

float GMSBAnalyzer::normAndErrorSquared(const TH3F& ggMETVsDiEMETVsInvMass, 
					const TH1F& controlMET, const TH1D* egMET, 
					const unsigned int maxNormBin, float& errSquared) const
{
  const float nGGMinusEW = 
    ggMETVsDiEMETVsInvMass.Integral(0, -1, 0, -1, 1, maxNormBin)/* - egMET->Integral(1, maxNormBin)*/;
  const float nControl = controlMET.Integral(1, maxNormBin);
  float norm = 0.0;
  if (nControl != 0.0) norm = nGGMinusEW/nControl;
  float nGGErrSquared = 0.0;
  float nEGErrSquared = 0.0;
  for (unsigned int iBin = 1; iBin <= maxNormBin; ++iBin) {
    const float ggErr = 
      ggMETVsDiEMETVsInvMass.ProjectionZ("ggMET", 0, -1, 0, -1, "e")->GetBinError(iBin);
//     const float egErr = egMET->GetBinError(iBin);
    nGGErrSquared+=(ggErr*ggErr);
    nEGErrSquared+=/*(egErr*egErr)*/0.0;
  }
  if (nGGMinusEW == 0.0) errSquared = 0.0;
  else {
    errSquared = 
      norm*norm*(((nGGErrSquared + nEGErrSquared)/(nGGMinusEW*nGGMinusEW)) + (1.0/nControl));
  }
  return norm;
}

float GMSBAnalyzer::razorNormAndErrorSquared(const TH2F& ggR2VsMR, const TH2F& controlR2VsMR, 
					     const TH2F& egR2VsMR, 
					     const unsigned int maxMRNormBin, 
					     const unsigned int maxR2NormBin, 
					     float& errSquared) const
{
  const float nGGMinusEW = ggR2VsMR.Integral(1, maxMRNormBin, 1, maxR2NormBin) - 
    egR2VsMR.Integral(1, maxMRNormBin, 1, maxR2NormBin);
  const float nControl = controlR2VsMR.Integral(1, maxMRNormBin, 1, maxR2NormBin);
  float norm = 0.0;
  if (nControl != 0.0) norm = nGGMinusEW/nControl;
  float nGGErrSquared = 0.0;
  float nEGErrSquared = 0.0;
  for (unsigned int iMRBin = 1; iMRBin <= maxMRNormBin; ++iMRBin) {
    for (unsigned int iR2Bin = 1; iR2Bin <= maxR2NormBin; ++iR2Bin) {
      const float ggErr = ggR2VsMR.GetBinError(iMRBin, iR2Bin);
      const float egErr = egR2VsMR.GetBinError(iMRBin, iR2Bin);
      nGGErrSquared+=(ggErr*ggErr);
      nEGErrSquared+=(egErr*egErr);
    }
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
      TH1F* hist = new TH1F(name.str().c_str(), "", 320, typicalVal - 8.0*sqrtTypicalVal, 
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
				   const TH1F& controlWeightsToy, TH1F* hist, 
				   const VFLOAT& lumiWeight, const VFLOAT& PUWeight) const
{
  for (VFLOAT_IT i = controlMETVec.begin(); i != controlMETVec.end(); ++i) {
    Int_t iDiEMET = 1;
    bool foundDiEMETBin = false;
    const unsigned int iEvt = i - controlMETVec.begin();
    while ((iDiEMET <= controlWeightsToy.GetNbinsX()) && !foundDiEMETBin) {
      if ((controlDiEMETVec[iEvt] >= controlWeightsToy.GetBinLowEdge(iDiEMET)) && 
	  (controlDiEMETVec[iEvt] < controlWeightsToy.GetBinLowEdge(iDiEMET + 1))) {
	foundDiEMETBin = true;
      }
      else ++iDiEMET;
    }
    if (iDiEMET > controlWeightsToy.GetNbinsX()) {
      cerr << "Error: di-EM ET bin corresponding to event with MET = " << *i;
      cerr << " GeV, di-EM ET = " << controlDiEMETVec[iEvt] << " GeV not found in histogram ";
      cerr << controlWeightsToy.GetName() << ".\n";
    }
    hist->Fill(*i, controlWeightsToy.GetBinContent(iDiEMET)*lumiWeight[iEvt]*PUWeight[iEvt]);
  }
}

void GMSBAnalyzer::reweightDefault(const VFLOAT& controlMETVec, const VFLOAT& controlMRVec, 
				   const VFLOAT& controlR2Vec, const VFLOAT& controlDiEMETVec, 
				   const TH1F& controlWeightsToy, TH1F* hist, TH2F* hist2D, 
				   const VFLOAT& lumiWeight, const VFLOAT& PUWeight) const
{
  for (VFLOAT_IT i = controlMETVec.begin(); i != controlMETVec.end(); ++i) {
    const unsigned int iEvt = i - controlMETVec.begin();
    Int_t iDiEMET = 1;
    bool foundDiEMETBin = false;
    while ((iDiEMET <= controlWeightsToy.GetNbinsX()) && !foundDiEMETBin) {
      if ((controlDiEMETVec[iEvt] >= controlWeightsToy.GetBinLowEdge(iDiEMET)) && 
	  (controlDiEMETVec[iEvt] < controlWeightsToy.GetBinLowEdge(iDiEMET + 1))) {
	foundDiEMETBin = true;
      }
      else ++iDiEMET;
    }
    if (iDiEMET > controlWeightsToy.GetNbinsX()) {
      cerr << "Error: di-EM ET bin corresponding to event with MET = " << *i;
      cerr << " GeV, di-EM ET = " << controlDiEMETVec[iEvt] << " GeV not found in histogram ";
      cerr << controlWeightsToy.GetName() << ".\n";
    }
    hist->Fill(*i, controlWeightsToy.GetBinContent(iDiEMET)*lumiWeight[iEvt]*PUWeight[iEvt]);
    hist2D->Fill(controlMRVec[iEvt], controlR2Vec[iEvt], 
		 controlWeightsToy.GetBinContent(iDiEMET)*lumiWeight[iEvt]*PUWeight[iEvt]);
  }
}

void GMSBAnalyzer::reweightDefault(const VFLOAT& controlMETVec, const VFLOAT& controlMRVec, 
				   const VFLOAT& controlR2Vec, const VFLOAT& controlDiEMETVec, 
				   const VFLOAT& controlNjVec, const unsigned int nNjBins, 
				   const Double_t* NjBins, 
				   const vector<vector<TH1F*> >& controlWeightsToys, TH1F* hist, 
				   TH2F* hist2D, const VFLOAT& lumiWeight, const VFLOAT& PUWeight, 
				   const unsigned int iToy) const
{
  //loop over events
  for (VFLOAT_IT i = controlMETVec.begin(); i != controlMETVec.end(); ++i) {
    const unsigned int iEvt = i - controlMETVec.begin();

    //get the number of jets in the event
    unsigned int iNjBin = nNjBins;
    for (unsigned int iNjBinBoundary = 0; iNjBinBoundary < nNjBins; ++iNjBinBoundary) {
      if ((controlNjVec[iEvt] >= NjBins[iNjBinBoundary]) && 
	  (controlNjVec[iEvt] < NjBins[iNjBinBoundary + 1])) iNjBin = iNjBinBoundary;
    }
    if (iNjBin == nNjBins) {
      cerr << "Error: no Nj bin found for event with MET = " << *i;
      cerr << " GeV, di-EM ET = " << controlDiEMETVec[iEvt] << " GeV, Nj = ";
      cerr << controlNjVec[iEvt] << ".  Assuming Nj = 0.\n";
      iNjBin = 0;
    }

    /*search the 0th toy weight histogram in the event's Nj bin for the di-EM pT bin of this 
      event--the binning is identical WITHIN 1 Nj BIN, so it doesn't matter which one is used

      binning is NOT the same for Nj = 0, Nj = 1, and Nj >= 2*/
    Int_t iDiEMET = 1;
    bool foundDiEMETBin = false;
    while ((iDiEMET <= controlWeightsToys[iNjBin][0]->GetNbinsX()) && !foundDiEMETBin) {
      if ((controlDiEMETVec[iEvt] >= controlWeightsToys[iNjBin][0]->GetBinLowEdge(iDiEMET)) && 
	  (controlDiEMETVec[iEvt] < controlWeightsToys[iNjBin][0]->GetBinLowEdge(iDiEMET + 1))) {
	foundDiEMETBin = true;
      }
      else ++iDiEMET;
    }
    if (iDiEMET > controlWeightsToys[iNjBin][0]->GetNbinsX()) {
      cerr << "Error: di-EM ET bin corresponding to event with MET = " << *i;
      cerr << " GeV, di-EM ET = " << controlDiEMETVec[iEvt] << " GeV not found in histogram ";
      cerr << controlWeightsToys[iNjBin][0]->GetName() << ".  Assuming maximum di-EM ET bin.\n";
      iDiEMET = controlWeightsToys[iNjBin][0]->GetNbinsX();
    }

    //weight by di-EM ET AND Nj
//     if (iNjBin == 1) {
      hist->Fill(*i, controlWeightsToys[iNjBin][iToy]->
		 GetBinContent(iDiEMET)*lumiWeight[iEvt]*PUWeight[iEvt]);
      hist2D->Fill(controlMRVec[iEvt], controlR2Vec[iEvt], 
		   controlWeightsToys[iNjBin][iToy]->
		   GetBinContent(iDiEMET)*lumiWeight[iEvt]*PUWeight[iEvt]);
//     }
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
				const VFLOAT& controlMETVec, const VFLOAT& controlDiEMETVec, 
				const VFLOAT& lumiWeightVec, const VFLOAT& PUWeightVec) const
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
    reweightDefault(controlMETVec, controlDiEMETVec, controlWeightsToy, hist, lumiWeightVec, 
		    PUWeightVec);
    controlFinalToy.push_back(hist);
  }
}

void GMSBAnalyzer::generateToys(vector<TH1F*>& controlFinalToy, vector<TH2F*>& controlFinalToy2D, 
				vector<TH1F*>& controlDiEMETToyDistsByBin, 
				const vector<TH1F*> controlWeights, const unsigned int nToys, 
				const string& controlSample, const Double_t* diEMETBins, 
				const unsigned int nMETBins, const Double_t* METBins, 
				const unsigned int nMRBins, const Double_t* MRBins, 
				const unsigned int nR2Bins, const Double_t* R2Bins, 
				const VFLOAT& controlMETVec, const VFLOAT& controlMRVec, 
				const VFLOAT& controlR2Vec, const VFLOAT& controlDiEMETVec, 
				const VFLOAT& lumiWeightVec, const VFLOAT& PUWeightVec) const
{
  bookToyDiEMETWeightsHistograms(controlWeights, controlSample, controlDiEMETToyDistsByBin);
  TRandom3 random(0); /*guarantees a unique seed cf. 
			http://root.cern.ch/root/html/TRandom3.html#TRandom3:TRandom3*/
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
    nameFinal.str("");
    nameFinal << controlSample << "FinalToyR2VsMR" << iToy;
    TH2F* hist2D = new TH2F(nameFinal.str().c_str(), "", nMRBins, MRBins, nR2Bins, R2Bins);
    hist2D->Sumw2();
    hist2D->SetFillStyle(0);
    makeToyDiEMETWeightsHistograms(random, controlWeightsToy, *controlWeights[0], 
				   controlDiEMETToyDistsByBin);
    reweightDefault(controlMETVec, controlMRVec, controlR2Vec, controlDiEMETVec, 
		    controlWeightsToy, hist, hist2D, lumiWeightVec, PUWeightVec);
    controlFinalToy.push_back(hist);
    controlFinalToy2D.push_back(hist2D);
  }
}

void GMSBAnalyzer::generateToys(vector<TH1F*>& controlDiEMETToyDistsByBin, 
				vector<TH1F*>& controlWeightsToys, 
				const vector<TH1F*> controlWeights, const unsigned int nToys, 
				const string& controlSample, const Double_t* diEMETBins) const
{
  bookToyDiEMETWeightsHistograms(controlWeights, controlSample, controlDiEMETToyDistsByBin);
  TRandom3 random(0); /*guarantees a unique seed cf. 
			http://root.cern.ch/root/html/TRandom3.html#TRandom3:TRandom3*/
  for (unsigned int iToy = 1; iToy <= nToys; ++iToy) {
    STRINGSTREAM nameWeights;
    nameWeights << controlSample << "WeightsToy" << iToy;
    TH1F* controlWeightsToy = 
      new TH1F(nameWeights.str().c_str(), "", controlWeights[0]->GetNbinsX(), diEMETBins);
    controlWeightsToy->Sumw2();
    makeToyDiEMETWeightsHistograms(random, *controlWeightsToy, *controlWeights[0], 
				   controlDiEMETToyDistsByBin);
    controlWeightsToys.push_back(controlWeightsToy);
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

void GMSBAnalyzer::makeToyMETDists(const unsigned int nToys, const string& controlSample, 
				   const unsigned int nMETBins, const Double_t* METBins, 
				   const unsigned int nMRBins, const Double_t* MRBins, 
				   const unsigned int nR2Bins, const Double_t* R2Bins, 
				   const VFLOAT& controlMETVec, const VFLOAT& controlMRVec, 
				   const VFLOAT& controlR2Vec, const VFLOAT& controlDiEMETVec, 
				   const VFLOAT& controlNjVec, const unsigned int nNjBins, 
				   const Double_t* NjBins, 
				   const vector<vector<TH1F*> >& controlWeightsToys, 
				   vector<TH1F*>& controlFinalToy, 
				   vector<TH2F*>& controlFinalToy2D, const VFLOAT& lumiWeightVec, 
				   const VFLOAT& PUWeightVec) const
{
  for (unsigned int iToy = 1; iToy <= nToys; ++iToy) {
    STRINGSTREAM nameFinal;
    nameFinal << controlSample << "FinalToy" << iToy;
    TH1F* hist = new TH1F(nameFinal.str().c_str(), "", nMETBins, METBins);
    hist->Sumw2();
    hist->SetFillStyle(0);
    nameFinal.str("");
    nameFinal << controlSample << "FinalToyR2VsMR" << iToy;
    TH2F* hist2D = new TH2F(nameFinal.str().c_str(), "", nMRBins, MRBins, nR2Bins, R2Bins);
    hist2D->Sumw2();
    hist2D->SetFillStyle(0);
    reweightDefault(controlMETVec, controlMRVec, controlR2Vec, controlDiEMETVec, controlNjVec, 
		    nNjBins, NjBins, controlWeightsToys, hist, hist2D, lumiWeightVec, PUWeightVec, 
		    iToy - 1);
    controlFinalToy.push_back(hist);
    controlFinalToy2D.push_back(hist2D);
  }
}

void GMSBAnalyzer::generateToyDijets(const unsigned int nEvts, TRandom& random, 
				     const vector<float>& j1JES, const vector<float>& j2JES, 
				     const vector<float>& j1JESErr, const vector<float>& j2JESErr, 
				     const vector<TLorentzVector>& j1UncorrP4, 
				     const vector<TLorentzVector>& j2UncorrP4, TH1F& dijetPT) const
{
  for (unsigned int iEvt = 0; iEvt < nEvts; ++iEvt) {
    float toyScale1 = random.Gaus(j1JES[iEvt], j1JESErr[iEvt]);
    float toyScale2 = random.Gaus(j2JES[iEvt], j2JESErr[iEvt]);
    float toyDijetPT = (toyScale1*j1UncorrP4[iEvt] + toyScale2*j2UncorrP4[iEvt]).Pt();
//     cerr << "toyScale1 = " << toyScale1 << endl;
//     cerr << "toyScale2 = " << toyScale2 << endl;
//     cerr << "toyDijetPT = " << toyDijetPT << endl;
    dijetPT.Fill(toyDijetPT);
  }
}

void GMSBAnalyzer::estimateJESError(const unsigned int nToys, 
				    vector<TH1F*>& ffDijetWeightsToyDists, 
				    const vector<TH1F*>& ffWeights, 
				    const unsigned int nDiEMETBins, const Double_t* diEMETBins, 
				    const vector<float>& ffJ1JES, const vector<float>& ffJ2JES, 
				    const vector<float>& ffJ1JESErr, 
				    const vector<float>& ffJ2JESErr, 
				    const vector<TLorentzVector>& ffJ1UncorrP4, 
				    const vector<TLorentzVector>& ffJ2UncorrP4, 
				    const vector<float>& ggJ1JES, const vector<float>& ggJ2JES, 
				    const vector<float>& ggJ1JESErr, 
				    const vector<float>& ggJ2JESErr, 
				    const vector<TLorentzVector>& ggJ1UncorrP4, 
				    const vector<TLorentzVector>& ggJ2UncorrP4) const
{
  //sanity check
  if ((ffJ1UncorrP4.size() != ffJ2UncorrP4.size()) || (ffJ2UncorrP4.size() != ffJ1JES.size()) || 
      (ffJ1JES.size() != ffJ2JES.size()) || (ffJ2JES.size() != ffJ1JESErr.size()) || 
      (ffJ1JESErr.size() != ffJ2JESErr.size())) {
    cerr << "Error: vector size mismatch.\n";
    cerr << "     Size of ffJ1UncorrP4: " << ffJ1UncorrP4.size() << endl;
    cerr << "     Size of ffJ2UncorrP4: " << ffJ2UncorrP4.size() << endl;
    cerr << "     Size of ffJ1JES: " << ffJ1JES.size() << endl;
    cerr << "     Size of ffJ2JES: " << ffJ2JES.size() << endl;
    cerr << "     Size of ffJ1JESErr: " << ffJ1JESErr.size() << endl;
    cerr << "     Size of ffJ2JESErr: " << ffJ2JESErr.size() << endl;
    cerr << "Quitting.\n";
    return;
  }
  if ((ggJ1UncorrP4.size() != ggJ2UncorrP4.size()) || (ggJ2UncorrP4.size() != ggJ1JES.size()) || 
      (ggJ1JES.size() != ggJ2JES.size()) || (ggJ2JES.size() != ggJ1JESErr.size()) || 
      (ggJ1JESErr.size() != ggJ2JESErr.size())) {
    cerr << "Error: vector size mismatch.\n";
    cerr << "     Size of ggJ1UncorrP4: " << ggJ1UncorrP4.size() << endl;
    cerr << "     Size of ggJ2UncorrP4: " << ggJ2UncorrP4.size() << endl;
    cerr << "     Size of ggJ1JES: " << ggJ1JES.size() << endl;
    cerr << "     Size of ggJ2JES: " << ggJ2JES.size() << endl;
    cerr << "     Size of ggJ1JESErr: " << ggJ1JESErr.size() << endl;
    cerr << "     Size of ggJ2JESErr: " << ggJ2JESErr.size() << endl;
    cerr << "Quitting.\n";
    return;
  }
  const unsigned int nFFEvts = ffJ1UncorrP4.size();
  const unsigned int nGGEvts = ggJ1UncorrP4.size();

  //book one histogram per dijet pT bin to hold the toy weights spreads
  bookToyDiEMETWeightsHistograms(ffWeights, "ffDijet", ffDijetWeightsToyDists);

  //loop over toys
  TRandom3 random;
  for (unsigned int iToy = 1; iToy <= nToys; ++iToy) {

    //book dijet pT histograms for this toy
    TH1F ffDijetPT(histName("ffDijetPT", "", iToy).c_str(), "", nDiEMETBins, diEMETBins);
    TH1F ggDijetPT(histName("ggDijetPT", "", iToy).c_str(), "", nDiEMETBins, diEMETBins);

    //generate toy dijets
//     cerr << "ff toy " << iToy << endl;
    generateToyDijets(nFFEvts, random, ffJ1JES, ffJ2JES, ffJ1JESErr, ffJ2JESErr, ffJ1UncorrP4, 
		      ffJ2UncorrP4, ffDijetPT);
//     cerr << "gg toy " << iToy << endl;
    generateToyDijets(nGGEvts, random, ggJ1JES, ggJ2JES, ggJ1JESErr, ggJ2JESErr, ggJ1UncorrP4, 
		      ggJ2UncorrP4, ggDijetPT);

    //normalize integral of ff to integral of gg
    ffDijetPT.Scale(ggDijetPT.GetEntries()/ffDijetPT.GetEntries());

    //fill weights histogram for this toy
    TH1F ffJESToyWeights(histName("ffJESToyWeights", "", iToy).c_str(), "", nDiEMETBins, 
			 diEMETBins);
    ffJESToyWeights.Divide(&ggDijetPT, &ffDijetPT);

    //fill the spread histograms
    for (Int_t iDijetPTBin = 1; iDijetPTBin <= ffJESToyWeights.GetNbinsX(); ++iDijetPTBin) {
//       cerr << "Dijet pT bin " << iDijetPTBin << ", gg bin content ";
//       cerr << ggDijetPT.GetBinContent(iDijetPTBin) << ", ff bin content ";
//       cerr << ffDijetPT.GetBinContent(iDijetPTBin) << ", weight ";
//       cerr << ffJESToyWeights.GetBinContent(iDijetPTBin) << endl;
      ffDijetWeightsToyDists[iDijetPTBin - 1]->Fill(ffJESToyWeights.GetBinContent(iDijetPTBin));
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

void GMSBAnalyzer::fillToyDistributions(vector<TH1F*>& controlMETToyDistsByBin, 
					vector<TH1F*>& controlRazorToyDistsByBin, 
					const TH1F& controlFinal, const TH2F& controlFinalRazor, 
					TCanvas& controlToyCanvas, vector<TH1F*>& controlFinalToy, 
					vector<TH2F*>& controlFinalToyRazor, 
					const unsigned int nToys, const string& controlSample, 
					const unsigned int nMETBins, const unsigned int nMRBins, 
					const unsigned int nR2Bins) const
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
  for (unsigned int iMRBin = 1; iMRBin <= nMRBins; ++iMRBin) {
    for (unsigned int iR2Bin = 1; iR2Bin <= nR2Bins; ++iR2Bin) {
      STRINGSTREAM name;
      name << controlSample << "RazorToyDist_MRBin" << iMRBin << "_R2Bin" << iR2Bin;
      float typicalVal = controlFinalRazor.GetBinContent(iMRBin, iR2Bin);
      float sqrtTypicalVal = sqrt(typicalVal);
      TH1F* hist = new TH1F(name.str().c_str(), "", 150, typicalVal - 15.0*sqrtTypicalVal, 
			    typicalVal + 15.0*sqrtTypicalVal);
      hist->Sumw2();
      controlRazorToyDistsByBin.push_back(hist);
    }
  }
  controlToyCanvas.cd();
  for (unsigned int iToy = 0; iToy < nToys; ++iToy) {
    if (iToy == 0) controlFinalToy[iToy]->Draw("HIST");
    else controlFinalToy[iToy]->Draw("HISTSAME");
    for (unsigned int iBin = 1; iBin <= nMETBins; ++iBin) {
      controlMETToyDistsByBin[iBin - 1]->Fill(controlFinalToy[iToy]->GetBinContent(iBin));
    }
    for (unsigned int iMRBin = 1; iMRBin <= nMRBins; ++iMRBin) {
      for (unsigned int iR2Bin = 1; iR2Bin <= nR2Bins; ++iR2Bin) {
	controlRazorToyDistsByBin[(iMRBin - 1)*nR2Bins + iR2Bin - 1]->
	  Fill(controlFinalToyRazor[iToy]->GetBinContent(iMRBin, iR2Bin));
      }
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
			     //stat--this might be wrong (
			     sqrt(controlFinal.GetBinError(iBin)*controlFinal.GetBinError(iBin) + 
			     controlNormErrSquared/*norm*/ + 
			     controlNorm*controlNorm*reweightingErrSquared/*reweighting*/));
  }
}

void GMSBAnalyzer::setRazorErrorBars(TH2F& controlFinal, 
				     const vector<TH1F*>& controlRazorToyDistsByBin, 
				     const vector<TH1F*>& controlLowSidebandRazorToyDistsByBin, 
				     const vector<TH1F*>& controlHighSidebandRazorToyDistsByBin, 
				     const float controlNorm, 
				     const float controlNormErrSquared, const bool doDefault) const
{
  Int_t nBinsX = controlFinal.GetNbinsX();
  Int_t nBinsY = controlFinal.GetNbinsY();
  for (int iBinX = 1; iBinX <= nBinsX; ++iBinX) {
    for (int iBinY = 1; iBinY <= nBinsY; ++iBinY) {
      const unsigned int binIndex = (iBinX - 1)*nBinsY + iBinY - 1;
      const float reweightingErr = controlRazorToyDistsByBin[binIndex]->GetRMS();
      float reweightingErrSquared = reweightingErr*reweightingErr;
      if ((STRING(controlFinal.GetName()) == "eeR2VsMRFinal") && doDefault) {
	const float lowSidebandReweightingErr = 
	  controlLowSidebandRazorToyDistsByBin[binIndex]->GetRMS();
	const float highSidebandReweightingErr = 
	  controlHighSidebandRazorToyDistsByBin[binIndex]->GetRMS();
	reweightingErrSquared+=4*(lowSidebandReweightingErr*lowSidebandReweightingErr + 
				  highSidebandReweightingErr*highSidebandReweightingErr);
      }
      controlFinal.SetBinError(iBinX, iBinY, 
			       //stat
			       sqrt(controlFinal.GetBinError(iBinX, iBinY)*
				    controlFinal.GetBinError(iBinX, iBinY) + 
				    controlNormErrSquared/*norm*/ + 
				    controlNorm*controlNorm*reweightingErrSquared/*reweighting*/));
    }
  }
}

void GMSBAnalyzer::makeFinalCanvas(TH1* hist, const Color_t lineColor, const Width_t lineWidth, 
				   const Style_t fillStyle, const Color_t fillColor, 
				   const Size_t markerSize, const Color_t markerColor, 
				   const string& drawOption) const
{
  hist->SetLineColor(lineColor);
  hist->SetLineWidth(lineWidth);
  hist->SetFillStyle(fillStyle);
  hist->SetFillColor(fillColor);
  hist->SetMarkerSize(markerSize);
  hist->SetMarkerColor(markerColor);
  hist->Draw(drawOption.c_str());
}

void GMSBAnalyzer::printDiObjectErrorMessage(const unsigned int diObjectSize, 
					     const string& diObjectType, const int evtCategory, 
					     const unsigned int iEvt) const
{
  cerr << "Error: " << diObjectType << "Size = " << diObjectSize << " in run ";
  cerr << susyEvent->runNumber << ", event " << susyEvent->eventNumber << ", lumi section ";
  cerr << susyEvent->luminosityBlockNumber << " (event " << (iEvt + 1) << "), category ";
  cerr << printCategory(evtCategory);
  cerr << ".  Not filling histogram of MET vs. " << diObjectType;
  cerr << "ET vs. invMass for this event.\n";
}

string GMSBAnalyzer::eventFileName(string outputFile, const string& sample) const
{
  return outputFile.replace(outputFile.find(".root"), 5, "_" + sample + ".txt");
}

float GMSBAnalyzer::HT() const
{
  //HT
  float evtHT = 0.0;

  //loop over AK5 PF jets
  map<TString, susy::PFJetCollection>::const_iterator iJets = susyEvent->pfJets.find("ak5");
  if (iJets != susyEvent->pfJets.end()) {
    for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
	 iJet != iJets->second.end(); ++iJet) {

      //increase HT for jets passing criteria
      unsigned int isGoodJet = isJet(*iJet, 5.0);
      if (isGoodJet == 1) {
	evtHT+=iJet->jecScaleFactors.find("L1FastL2L3")->second*iJet->momentum.Et();
      }
    }
  }

  //error
  else cerr << "Error: ak5 PF jets collection not found.\n";

  //return
  return evtHT;
}

float GMSBAnalyzer::MHT() const
{
  //MHT
  TLorentzVector evtMHT;

  //loop over AK5 PF jets
  map<TString, susy::PFJetCollection>::const_iterator iJets = susyEvent->pfJets.find("ak5");
  if (iJets != susyEvent->pfJets.end()) {
    for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
	 iJet != iJets->second.end(); ++iJet) {

      //increase MHT for jets passing criteria
      unsigned int isGoodJet = isJet(*iJet, 5.0);
      if (isGoodJet == 1) {
	evtMHT+=iJet->jecScaleFactors.find("L1FastL2L3")->second*iJet->momentum;
      }
    }
  }

  //error
  else cerr << "Error: ak5 PF jets collection not found.\n";

  //return
  return evtMHT.Et(); //ET or pT?
}

unsigned int GMSBAnalyzer::numJets(const float absEtaMax) const
{
  //jet counter
  unsigned int nJets = 0;

  //loop over AK5 PF jets
  map<TString, susy::PFJetCollection>::const_iterator iJets = susyEvent->pfJets.find("ak5");
  if (iJets != susyEvent->pfJets.end()) {
    for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
	 iJet != iJets->second.end(); ++iJet) {

      //count jets passing criteria
      unsigned int isGoodJet = isJet(*iJet, absEtaMax);
      if (isGoodJet == 1) ++nJets;
    }
  }

  //error
  else cerr << "Error: ak5 PF jets collection not found.\n";

  //return
  return nJets;
}

float GMSBAnalyzer::leadingJetET() const
{
  //leading jet ET
  float ET1 = 0.0;

  //loop over AK5 PF jets
  map<TString, susy::PFJetCollection>::const_iterator iJets = susyEvent->pfJets.find("ak5");
  if (iJets != susyEvent->pfJets.end()) {
    for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
	 iJet != iJets->second.end(); ++iJet) {

      //count jets passing criteria
      unsigned int isGoodJet = isJet(*iJet, 5.0, 0.0);
      if (isGoodJet == 1) {
	const float ET = iJet->jecScaleFactors.find("L1FastL2L3")->second*iJet->momentum.Et();
	if (ET > ET1) ET1 = ET;
      }
    }
  }

  //error
  else cerr << "Error: ak5 PF jets collection not found.\n";

  //return
  return ET1;
}

unsigned int GMSBAnalyzer::isJet(const susy::PFJet& iJet, const float absEtaMax,
				 const float ETMin) const
{
  //pass flag
  unsigned int jet = 0;

  //compute corrected P4
  float storedCorr = -1.0;
  map<TString, Float_t>::const_iterator iCorr = iJet.jecScaleFactors.find("L1FastL2L3");
  if (iCorr != iJet.jecScaleFactors.end()) storedCorr = iCorr->second;
  TLorentzVector corrP4 = storedCorr*iJet.momentum;
 
  //count as a jet if it passes corrected pT, eta, and jet ID cuts
  if ((corrP4.Et() >= ETMin) && (fabs(corrP4.Eta()) <= absEtaMax)) {
    bool passedJetID = false;
    if ((iJet.neutralHadronEnergy/iJet.momentum.Energy() < 0.99) &&
	(iJet.neutralEmEnergy/iJet.momentum.Energy() < 0.99) &&
	((unsigned int)iJet.nConstituents > 1)) {
      if (fabs(iJet.momentum.Eta()) < 2.4) {
	if ((iJet.chargedHadronEnergy > 0.0) &&
	    ((int)iJet.chargedMultiplicity > 0) &&
	    (iJet.chargedEmEnergy/iJet.momentum.Energy() < 0.99)) passedJetID = true;
      }
      else passedJetID = true;
    }
    if (passedJetID) {
 
      //overlap flag false unless a specific overlap is found
      bool overlap = false;

      //loop over electrons
      map<TString, susy::ElectronCollection>::const_iterator iElectronMap =
	susyEvent->electrons.find("gsfElectrons");
      if (iElectronMap != susyEvent->electrons.end()) {
	susy::ElectronCollection::const_iterator iElectron = iElectronMap->second.begin();
	while ((iElectron != iElectronMap->second.end()) && !overlap) {
 
	  //Ulla doesn't use the electron SC position for eta/phi
// 	  //check that SCs were properly stored for this electron
// 	  unsigned int i = (unsigned int)iElectron->superClusterIndex;
// 	  unsigned int size = susyEvent->superClusters.size();
// 	  if (i >= size) {
// 	    cerr << "Error: electron SC index " << i << " >= size of SC collection " << size;
// 	    cerr << ".\n";
// 	    return 2;
// 	  }
 
	  //identify this electron as one which should not be allowed to overlap with a jet
	  float pT = iElectron->momentum.Pt();
	  // 	    float eta = susyEvent->superClusters[i].position.Eta();
	  float eta = iElectron->momentum.Eta();
	  if (iElectron->isPF() &&
	      (pT > 15.0) &&
	      (fabs(eta) < absEtaMax) &&
	      ((iElectron->chargedHadronIso +
		iElectron->photonIso +
		iElectron->neutralHadronIso)/pT < 0.2) &&
	      (deltaR(corrP4.Eta(), corrP4.Phi(),
// 		      eta, susyEvent->superClusters[i].position.Phi()) < 0.5)) {
		      eta, iElectron->momentum.Phi()) < 0.5)) { /*Ulla doesn't use the electron 
								  SC position for eta/phi*/
	    overlap = true;
	  }
	  else ++iElectron;
	}
      }
//       else {
 
// 	//error
// // 	cerr << "Error: gsfElectrons collection not found.\n"; /*GSF electrons missing from a 
// // 								 number of events in MC due to 
// // 								 pT/eta acceptance*/
// 	jet = 3;
//       }
 
      //loop over muons
      vector<susy::Muon>::const_iterator iMuon = susyEvent->muons.begin();
      while ((iMuon != susyEvent->muons.end()) && !overlap) {
 
	//check that tracks were properly stored for this muon
	//in MC sometimes get an unphysical track index
	int i = (int)iMuon->trackIndex;
	int size = (int)susyEvent->tracks.size();
	if (i >= size) {
	  cerr << "Error: muon track index " << i << " >= size of track collection " << size;
	  cerr << ".\n";
	  return 2;
	}
 
	//identify this muon as one which should not be allowed to overlap with a jet
	float pT = iMuon->momentum.Pt();
	float eta = iMuon->momentum.Eta();
	map<TString, UChar_t>::const_iterator iID =
	  iMuon->idPairs.find("muidGlobalMuonPromptTight");
	if ((iID != iMuon->idPairs.end()) && ((unsigned int)iID->second == 1) &&
	    (susyEvent->tracks[iMuon->trackIndex].d0() < 0.02) &&
	    (susyEvent->tracks[iMuon->trackIndex].dz() < 0.5) &&
	    (pT > 15.0) &&
	    (fabs(eta) < absEtaMax) &&
	    ((iMuon->ecalIsoR03 + iMuon->hcalIsoR03 + iMuon->trackIsoR03)/pT < 0.2) &&
	    (deltaR(corrP4.Eta(), corrP4.Phi(), eta, iMuon->momentum.Phi()) < 0.5)) {
	  overlap = true;
	}
	else ++iMuon;
      }
 
      //loop over photons
      map<TString, susy::PhotonCollection>::const_iterator iPhotonMap =
	susyEvent->photons.find((const char*)tag_);
      susy::PhotonCollection::const_iterator iPhoton = iPhotonMap->second.begin();
      if (iPhotonMap != susyEvent->photons.end()) {
	while ((iPhoton != iPhotonMap->second.end()) && !overlap) {
 
	  //identify this photon as one which should not be allowed to overlap with a jet
	  if (susyCategory->getIsDeciding(tag_, iPhoton - iPhotonMap->second.begin()) &&
	      (deltaR(corrP4.Eta(), corrP4.Phi(),
		      iPhoton->caloPosition.Eta(), iPhoton->caloPosition.Phi()) < 0.5)) {
	    overlap = true;
	  }
	  else ++iPhoton;
	}
      }
//       else {
 
// 	//error
// // 	cerr << "Error: " << tag_ << " collection not found.\n";
// 	jet = 3;
//       }
 
      //if the jet didn't overlap with an electron, muon, or primary EM object, it counts
      if (!overlap) jet = 1;
    }
  }
 
  //return
  return jet;
}

void GMSBAnalyzer::runMETAnalysis(const std::string outputFile)
{
   if (fChain == 0) return;

   //output files of event counts
   ofstream ggEvtFile(eventFileName(outputFile, "gg").c_str());
   ofstream egEvtFile(eventFileName(outputFile, "eg").c_str());
   ofstream eeEvtFile(eventFileName(outputFile, "ee").c_str());
   ofstream ffEvtFile(eventFileName(outputFile, "ff").c_str());

   //open file
   TFile out(outputFile.c_str(), "RECREATE");
   out.cd();

   //define constants
   const bool kSumW2 = true;
   const bool kSetGrid = true;
   const unsigned int maxNormBin = 4;     //MET normalization region
   const unsigned int maxMRNormBin = 6;   //MR normalization region
   const unsigned int maxR2NormBin = 10;  //R2 normalization region
   const float nToys = 1000;              //number of toys for the MET shape error from reweighting
   const unsigned int nMETBins = 13;      //number of MET bins
//    const unsigned int nDiEMETBins = 25;   //number of di-EM ET bins
   const unsigned int nDiEMETBins = 15;   //number of di-EM ET bins
//    const unsigned int nDiEMETBins0j = 12;   //number of di-EM ET bins
//    const unsigned int nDiEMETBins1j = 14;   //number of di-EM ET bins
   const unsigned int nDiEMETBins0j = 13;   //number of di-EM ET bins
   const unsigned int nDiEMETBins1j = 20;   //number of di-EM ET bins
   const unsigned int nNjBins = 3;        //number of Nj bins
//    const unsigned int nInvMassBins = 150; //number of invariant mass bins
   const unsigned int nInvMassBins = 500; //number of invariant mass bins
   const unsigned int nMRBins = 20;       //number of MR bins
   const unsigned int nR2Bins = 25;       //number of R2 bins
   const Double_t METBins[14] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 70.0, 
				 100.0, 150.0, 500.0};            //MET bin boundaries
//    const Double_t diEMETBins[26] = {0.0,                          //di-EM ET bin boundaries
// 				    3.0, 6.0, 9.0, 12.0, 15.0, 
// 				    18.0, 21.0, 24.0, 27.0, 30.0, 
// 				    35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 
// 				    70.0, 80.0, 90.0, 
// 				    100.0, 150.0, 200.0, 
// 				    300.0, 400.0, 650.0};
   const Double_t diEMETBins[16] = {0.0,                          //di-EM ET bin boundaries
				    5.0, 10.0, 15.0, 
				    20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 
				    80.0, 100.0, 
				    120.0, 
				    150.0, 
				    200.0, 
				    650.0};
//    const Double_t diEMETBins0j[13] = {0.0,                        //di-EM ET bin boundaries
// 				      5.0, 10.0, 15.0, 
// 				      20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 
// 				      80.0, 100.0, 
// 				      150.0};
//    const Double_t diEMETBins1j[15] = {0.0,                        //di-EM ET bin boundaries
// 				      5.0, 10.0, 15.0, 
// 				      20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 
// 				      80.0, 100.0, 
// 				      120.0, 
// 				      150.0, 
// 				      200.0};
   const Double_t diEMETBins0j[14] = {0.0,                        //di-EM ET bin boundaries
				      5.0, 10.0, 15.0, 
				      20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 
				      80.0, 100.0, 
				      150.0, 650.0};
   const Double_t diEMETBins1j[21] = {0.0,                        //di-EM ET bin boundaries
				      5.0, 10.0, 15.0, 
				      20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 
				      80.0, 100.0, 
				      120.0, 
				      150.0, 200.0, 250.0, 
				      300.0, 400.0, 500.0, 600.0, 700.0};
   const Double_t NjBins[4] = {-0.5, 0.5, 1.5, 99.5};             //Nj bin boundaries
//    const Double_t invMassBins[151] = {0.0,                        //invariant mass bin boundaries
// 				      2.0, 4.0, 6.0, 8.0, 10.0, 
// 				      12.0, 14.0, 16.0, 18.0, 20.0, 
// 				      22.0, 24.0, 26.0, 28.0, 30.0, 
// 				      32.0, 34.0, 36.0, 38.0, 40.0, 
// 				      42.0, 44.0, 46.0, 48.0, 50.0, 
// 				      52.0, 54.0, 56.0, 58.0, 60.0, 
// 				      62.0, 64.0, 66.0, 68.0, 70.0, 
// 				      72.0, 74.0, 76.0, 78.0, 80.0, 
// 				      82.0, 84.0, 86.0, 88.0, 90.0, 
// 				      92.0, 94.0, 96.0, 98.0, 100.0, 
// 				      102.0, 104.0, 106.0, 108.0, 110.0, 
// 				      112.0, 114.0, 116.0, 118.0, 120.0, 
// 				      122.0, 124.0, 126.0, 128.0, 130.0, 
// 				      132.0, 134.0, 136.0, 138.0, 140.0, 
// 				      142.0, 144.0, 146.0, 148.0, 150.0, 
// 				      152.0, 154.0, 156.0, 158.0, 160.0, 
// 				      162.0, 164.0, 166.0, 168.0, 170.0, 
// 				      172.0, 174.0, 176.0, 178.0, 180.0, 
// 				      182.0, 184.0, 186.0, 188.0, 190.0, 
// 				      192.0, 194.0, 196.0, 198.0, 200.0, 
// 				      202.0, 204.0, 206.0, 208.0, 210.0, 
// 				      212.0, 214.0, 216.0, 218.0, 220.0, 
// 				      222.0, 224.0, 226.0, 228.0, 230.0, 
// 				      232.0, 234.0, 236.0, 238.0, 240.0, 
// 				      242.0, 244.0, 246.0, 248.0, 250.0, 
// 				      252.0, 254.0, 256.0, 258.0, 260.0, 
// 				      262.0, 264.0, 266.0, 268.0, 270.0, 
// 				      272.0, 274.0, 276.0, 278.0, 280.0, 
// 				      282.0, 284.0, 286.0, 288.0, 290.0, 
// 				      292.0, 294.0, 296.0, 298.0, 300.0};
   const Double_t invMassBins[501] = {0.0,                        //invariant mass bin boundaries
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
				      292.0, 294.0, 296.0, 298.0, 300.0, 
				      302.0, 304.0, 306.0, 308.0, 310.0, 
				      312.0, 314.0, 316.0, 318.0, 320.0, 
				      322.0, 324.0, 326.0, 328.0, 330.0, 
				      332.0, 334.0, 336.0, 338.0, 340.0, 
				      342.0, 344.0, 346.0, 348.0, 350.0, 
				      352.0, 354.0, 356.0, 358.0, 360.0, 
				      362.0, 364.0, 366.0, 368.0, 370.0, 
				      372.0, 374.0, 376.0, 378.0, 380.0, 
				      382.0, 384.0, 386.0, 388.0, 390.0, 
				      392.0, 394.0, 396.0, 398.0, 400.0, 
				      402.0, 404.0, 406.0, 408.0, 410.0, 
				      412.0, 414.0, 416.0, 418.0, 420.0, 
				      422.0, 424.0, 426.0, 428.0, 430.0, 
				      432.0, 434.0, 436.0, 438.0, 440.0, 
				      442.0, 444.0, 446.0, 448.0, 450.0, 
				      452.0, 454.0, 456.0, 458.0, 460.0, 
				      462.0, 464.0, 466.0, 468.0, 470.0, 
				      472.0, 474.0, 476.0, 478.0, 480.0, 
				      482.0, 484.0, 486.0, 488.0, 490.0, 
				      492.0, 494.0, 496.0, 498.0, 500.0, 
				      502.0, 504.0, 506.0, 508.0, 510.0, 
				      512.0, 514.0, 516.0, 518.0, 520.0, 
				      522.0, 524.0, 526.0, 528.0, 530.0, 
				      532.0, 534.0, 536.0, 538.0, 540.0, 
				      542.0, 544.0, 546.0, 548.0, 550.0, 
				      552.0, 554.0, 556.0, 558.0, 560.0, 
				      562.0, 564.0, 566.0, 568.0, 570.0, 
				      572.0, 574.0, 576.0, 578.0, 580.0, 
				      582.0, 584.0, 586.0, 588.0, 590.0, 
				      592.0, 594.0, 596.0, 598.0, 600.0, 
				      602.0, 604.0, 606.0, 608.0, 610.0, 
				      612.0, 614.0, 616.0, 618.0, 620.0, 
				      622.0, 624.0, 626.0, 628.0, 630.0, 
				      632.0, 634.0, 636.0, 638.0, 640.0, 
				      642.0, 644.0, 646.0, 648.0, 650.0, 
				      652.0, 654.0, 656.0, 658.0, 660.0, 
				      662.0, 664.0, 666.0, 668.0, 670.0, 
				      672.0, 674.0, 676.0, 678.0, 680.0, 
				      682.0, 684.0, 686.0, 688.0, 690.0, 
				      692.0, 694.0, 696.0, 698.0, 700.0, 
				      702.0, 704.0, 706.0, 708.0, 710.0, 
				      712.0, 714.0, 716.0, 718.0, 720.0, 
				      722.0, 724.0, 726.0, 728.0, 730.0, 
				      732.0, 734.0, 736.0, 738.0, 740.0, 
				      742.0, 744.0, 746.0, 748.0, 750.0, 
				      752.0, 754.0, 756.0, 758.0, 760.0, 
				      762.0, 764.0, 766.0, 768.0, 770.0, 
				      772.0, 774.0, 776.0, 778.0, 780.0, 
				      782.0, 784.0, 786.0, 788.0, 790.0, 
				      792.0, 794.0, 796.0, 798.0, 800.0, 
				      802.0, 804.0, 806.0, 808.0, 810.0, 
				      812.0, 814.0, 816.0, 818.0, 820.0, 
				      822.0, 824.0, 826.0, 828.0, 830.0, 
				      832.0, 834.0, 836.0, 838.0, 840.0, 
				      842.0, 844.0, 846.0, 848.0, 850.0, 
				      852.0, 854.0, 856.0, 858.0, 860.0, 
				      862.0, 864.0, 866.0, 868.0, 870.0, 
				      872.0, 874.0, 876.0, 878.0, 880.0, 
				      882.0, 884.0, 886.0, 888.0, 890.0, 
				      892.0, 894.0, 896.0, 898.0, 900.0, 
				      902.0, 904.0, 906.0, 908.0, 910.0, 
				      912.0, 914.0, 916.0, 918.0, 920.0, 
				      922.0, 924.0, 926.0, 928.0, 930.0, 
				      932.0, 934.0, 936.0, 938.0, 940.0, 
				      942.0, 944.0, 946.0, 948.0, 950.0, 
				      952.0, 954.0, 956.0, 958.0, 960.0, 
				      962.0, 964.0, 966.0, 968.0, 970.0, 
				      972.0, 974.0, 976.0, 978.0, 980.0, 
				      982.0, 984.0, 986.0, 988.0, 990.0, 
				      992.0, 994.0, 996.0, 998.0, 1000.0};
   const Double_t MRBins[21] = {0.0, 
				50.0, 100.0, 150.0, 200.0, 250.0, 
				300.0, 350.0, 400.0, 450.0, 500.0, 
				550.0, 600.0, 650.0, 700.0, 750.0, 
				800.0, 850.0, 900.0, 950.0, 1000.0};
   const Double_t R2Bins[26] = {0.0, 
				0.02, 0.04, 0.06, 0.08, 0.1, 
				0.12, 0.14, 0.16, 0.18, 0.2, 
				0.22, 0.24, 0.26, 0.28, 0.3, 
				0.32, 0.34, 0.36, 0.38, 0.4, 
				0.42, 0.44, 0.46, 0.48, 0.5};
   const float egMisIDRate = 0.015;      /*take from CMS AN-2011/515 for now; can be computed with 
					   Z*/
   const float egMisIDRateErr = 0.005;   /*take from CMS AN-2011/515 for now; can be computed with 
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

   //sigmaIetaIeta histograms (sanity check)
   TH1F ggSigmaIetaIeta("ggSigmaIetaIeta", "", 32, 0.0, 0.016);
   TH1F egSigmaIetaIeta("egSigmaIetaIeta", "", 32, 0.0, 0.016);
   TH1F eeSigmaIetaIeta("eeSigmaIetaIeta", "", 32, 0.0, 0.016);
   TH1F ffSigmaIetaIeta("ffSigmaIetaIeta", "", 32, 0.0, 0.016);
   setHistogramOptions(ggSigmaIetaIeta, "#sigma_{i#etai#eta}", "", "", kSumW2);
   setHistogramOptions(egSigmaIetaIeta, "#sigma_{i#etai#eta}", "", "", kSumW2);
   setHistogramOptions(eeSigmaIetaIeta, "#sigma_{i#etai#eta}", "", "", kSumW2);
   setHistogramOptions(ffSigmaIetaIeta, "#sigma_{i#etai#eta}", "", "", kSumW2);

   //R9 histograms (sanity check)
   TH1F ggR9("ggR9", "", 30, 0.7, 1.0);
   TH1F egR9("egR9", "", 30, 0.7, 1.0);
   TH1F eeR9("eeR9", "", 30, 0.7, 1.0);
   TH1F ffR9("ffR9", "", 30, 0.7, 1.0);
   setHistogramOptions(ggR9, "R9", "", "", kSumW2);
   setHistogramOptions(egR9, "R9", "", "", kSumW2);
   setHistogramOptions(eeR9, "R9", "", "", kSumW2);
   setHistogramOptions(ffR9, "R9", "", "", kSumW2);

   //ECAL isolation histograms (sanity check)
   TH1F ggECALIso("ggECALIso", "", 20, 0.0, 20.0);
   TH1F egECALIso("egECALIso", "", 20, 0.0, 20.0);
   TH1F eeECALIso("eeECALIso", "", 20, 0.0, 20.0);
   TH1F ffECALIso("ffECALIso", "", 20, 0.0, 20.0);
   setHistogramOptions(ggECALIso, "I_{ECAL} (GeV)", "", "", kSumW2);
   setHistogramOptions(egECALIso, "I_{ECAL} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeECALIso, "I_{ECAL} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffECALIso, "I_{ECAL} (GeV)", "", "", kSumW2);

   //HCAL isolation histograms (sanity check)
   TH1F ggHCALIso("ggHCALIso", "", 20, 0.0, 20.0);
   TH1F egHCALIso("egHCALIso", "", 20, 0.0, 20.0);
   TH1F eeHCALIso("eeHCALIso", "", 20, 0.0, 20.0);
   TH1F ffHCALIso("ffHCALIso", "", 20, 0.0, 20.0);
   setHistogramOptions(ggHCALIso, "I_{HCAL} (GeV)", "", "", kSumW2);
   setHistogramOptions(egHCALIso, "I_{HCAL} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeHCALIso, "I_{HCAL} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffHCALIso, "I_{HCAL} (GeV)", "", "", kSumW2);

   //track isolation histograms (sanity check)
   TH1F ggTrackIso("ggTrackIso", "", 20, 0.0, 20.0);
   TH1F egTrackIso("egTrackIso", "", 20, 0.0, 20.0);
   TH1F eeTrackIso("eeTrackIso", "", 20, 0.0, 20.0);
   TH1F ffTrackIso("ffTrackIso", "", 20, 0.0, 20.0);
   setHistogramOptions(ggTrackIso, "I_{track} (GeV)", "", "", kSumW2);
   setHistogramOptions(egTrackIso, "I_{track} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeTrackIso, "I_{track} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffTrackIso, "I_{track} (GeV)", "", "", kSumW2);

   //leading photon ET histograms (sanity check)
   TH1F ggLeadingPhotonET("ggLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F egLeadingPhotonET("egLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F eeLeadingPhotonET("eeLeadingPhotonET", "", 150, 0.0, 300.0);
   TH1F ffLeadingPhotonET("ffLeadingPhotonET", "", 150, 0.0, 300.0);
   setHistogramOptions(ggLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(egLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffLeadingPhotonET, "Leading photon E_{T} (GeV)", "", "", kSumW2);

   //trailing photon ET histograms (sanity check)
   TH1F ggTrailingPhotonET("ggTrailingPhotonET", "", 150, 0.0, 300.0);
   TH1F egTrailingPhotonET("egTrailingPhotonET", "", 150, 0.0, 300.0);
   TH1F eeTrailingPhotonET("eeTrailingPhotonET", "", 150, 0.0, 300.0);
   TH1F ffTrailingPhotonET("ffTrailingPhotonET", "", 150, 0.0, 300.0);
   setHistogramOptions(ggTrailingPhotonET, "Trailing photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(egTrailingPhotonET, "Trailing photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeTrailingPhotonET, "Trailing photon E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffTrailingPhotonET, "Trailing photon E_{T} (GeV)", "", "", kSumW2);

   //pThat histograms
   TH1F ggPTHat("ggPTHat", "", 250, 0.0, 500.0);
   TH1F egPTHat("egPTHat", "", 250, 0.0, 500.0);
   TH1F eePTHat("eePTHat", "", 250, 0.0, 500.0);
   TH1F ffPTHat("ffPTHat", "", 250, 0.0, 500.0);
   setHistogramOptions(ggPTHat, "#hat{p}_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(egPTHat, "#hat{p}_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eePTHat, "#hat{p}_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffPTHat, "#hat{p}_{T} (GeV)", "", "", kSumW2);

   //nPV histograms (sanity check)
   TH1F ggNPV("ggNPV", "", 36, -0.5, 35.5);
   TH1F egNPV("egNPV", "", 36, -0.5, 35.5);
   TH1F eeNPV("eeNPV", "", 36, -0.5, 35.5);
   TH1F ffNPV("ffNPV", "", 36, -0.5, 35.5);
   setHistogramOptions(ggNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(egNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(eeNPV, "n_{PV}", "", "", kSumW2);
   setHistogramOptions(ffNPV, "n_{PV}", "", "", kSumW2);

   //rho histograms (sanity check)
   TH1F ggRho("ggRho", "", 200, 0.0, 20.0);
   TH1F egRho("egRho", "", 200, 0.0, 20.0);
   TH1F eeRho("eeRho", "", 200, 0.0, 20.0);
   TH1F ffRho("ffRho", "", 200, 0.0, 20.0);
   setHistogramOptions(ggRho, "#rho (GeV/#eta#cdot#phi)", "", "", kSumW2);
   setHistogramOptions(egRho, "#rho (GeV/#eta#cdot#phi)", "", "", kSumW2);
   setHistogramOptions(eeRho, "#rho (GeV/#eta#cdot#phi)", "", "", kSumW2);
   setHistogramOptions(ffRho, "#rho (GeV/#eta#cdot#phi)", "", "", kSumW2);

   //PF electron multiplicity (sanity check)
   TH1F ggPFElectronMultiplicity("ggPFElectronMultiplicity", "", 5, -0.5, 4.5);
   TH1F egPFElectronMultiplicity("egPFElectronMultiplicity", "", 5, -0.5, 4.5);
   TH1F eePFElectronMultiplicity("eePFElectronMultiplicity", "", 5, -0.5, 4.5);
   TH1F ffPFElectronMultiplicity("ffPFElectronMultiplicity", "", 5, -0.5, 4.5);
   setHistogramOptions(ggPFElectronMultiplicity, "N_{PF e}", "", "", kSumW2);
   setHistogramOptions(egPFElectronMultiplicity, "N_{PF e}", "", "", kSumW2);
   setHistogramOptions(eePFElectronMultiplicity, "N_{PF e}", "", "", kSumW2);
   setHistogramOptions(ffPFElectronMultiplicity, "N_{PF e}", "", "", kSumW2);
 
   //PF muon multiplicity (sanity check)
   TH1F ggMuonMultiplicity("ggMuonMultiplicity", "", 5, -0.5, 4.5);
   TH1F egMuonMultiplicity("egMuonMultiplicity", "", 5, -0.5, 4.5);
   TH1F eeMuonMultiplicity("eeMuonMultiplicity", "", 5, -0.5, 4.5);
   TH1F ffMuonMultiplicity("ffMuonMultiplicity", "", 5, -0.5, 4.5);
   setHistogramOptions(ggMuonMultiplicity, "N_{#mu}", "", "", kSumW2);
   setHistogramOptions(egMuonMultiplicity, "N_{#mu}", "", "", kSumW2);
   setHistogramOptions(eeMuonMultiplicity, "N_{#mu}", "", "", kSumW2);
   setHistogramOptions(ffMuonMultiplicity, "N_{#mu}", "", "", kSumW2);

   //uniform binning ee, gg, and ff di-EM ET histogram
   vector<TH1F*> ggDiEMETUniform;
   vector<TH1F*> ffDiEMETUniform;
   vector<TH1F*> eeDiEMETUniform;
   for (unsigned int iNjBin = 1; iNjBin <= nNjBins; ++iNjBin) {
     ggDiEMETUniform.push_back(new TH1F(histName("ggDiEMETUniform", "NjBin", iNjBin).c_str(), 
					"", 500, 0.0, 500.0));
     ffDiEMETUniform.push_back(new TH1F(histName("ffDiEMETUniform", "NjBin", iNjBin).c_str(), 
					"", 500, 0.0, 500.0));
     eeDiEMETUniform.push_back(new TH1F(histName("eeDiEMETUniform", "NjBin", iNjBin).c_str(), 
					"", 500, 0.0, 500.0));
     setHistogramOptions(*ggDiEMETUniform[ggDiEMETUniform.size() - 1], "Di-EM E_{T} (GeV)", 
			 "", "", kSumW2);
     setHistogramOptions(*ffDiEMETUniform[ffDiEMETUniform.size() - 1], "Di-EM E_{T} (GeV)", 
			 "", "", kSumW2);
     setHistogramOptions(*eeDiEMETUniform[eeDiEMETUniform.size() - 1], "Di-EM E_{T} (GeV)", 
			 "", "", kSumW2);
   }

   //uniform binning ee, gg, and ff dijet ET histogram
   vector<TH1F*> ggDijetETUniform;
   vector<TH1F*> ffDijetETUniform;
   vector<TH1F*> eeDijetETUniform;
   for (unsigned int iNjBin = 1; iNjBin <= nNjBins; ++iNjBin) {
     ggDijetETUniform.push_back(new TH1F(histName("ggDijetETUniform", "NjBin", iNjBin).c_str(), 
					 "", 500, 0.0, 500.0));
     ffDijetETUniform.push_back(new TH1F(histName("ffDijetETUniform", "NjBin", iNjBin).c_str(), 
					 "", 500, 0.0, 500.0));
     eeDijetETUniform.push_back(new TH1F(histName("eeDijetETUniform", "NjBin", iNjBin).c_str(), 
					 "", 500, 0.0, 500.0));
     setHistogramOptions(*ggDijetETUniform[ggDijetETUniform.size() - 1], "Dijet E_{T} (GeV)", 
			 "", "", kSumW2);
     setHistogramOptions(*ffDijetETUniform[ffDijetETUniform.size() - 1], "Dijet E_{T} (GeV)", 
			 "", "", kSumW2);
     setHistogramOptions(*eeDijetETUniform[eeDijetETUniform.size() - 1], "Dijet E_{T} (GeV)", 
			 "", "", kSumW2);
   }

   //sizes of di-objects
   TH1F ggDijetSize("ggDijetSize", "ggDijetSize", 6, -0.5, 5.5);
   TH1F egDijetSize("egDijetSize", "egDijetSize", 6, -0.5, 5.5);
   TH1F eeDijetSize("eeDijetSize", "eeDijetSize", 6, -0.5, 5.5);
   TH1F ffDijetSize("ffDijetSize", "ffDijetSize", 6, -0.5, 5.5);
   TH1F ggDiEMSize("ggDiEMSize", "ggDiEMSize", 6, -0.5, 5.5);
   TH1F egDiEMSize("egDiEMSize", "egDiEMSize", 6, -0.5, 5.5);
   TH1F eeDiEMSize("eeDiEMSize", "eeDiEMSize", 6, -0.5, 5.5);
   TH1F ffDiEMSize("ffDiEMSize", "ffDiEMSize", 6, -0.5, 5.5);
   setHistogramOptions(ggDijetSize, "No. selected EM objects per event", "", "");
   setHistogramOptions(egDijetSize, "No. selected EM objects per event", "", "");
   setHistogramOptions(eeDijetSize, "No. selected EM objects per event", "", "");
   setHistogramOptions(ffDijetSize, "No. selected EM objects per event", "", "");
   setHistogramOptions(ggDiEMSize, "No. selected jets per event", "", "");
   setHistogramOptions(egDiEMSize, "No. selected jets per event", "", "");
   setHistogramOptions(eeDiEMSize, "No. selected jets per event", "", "");
   setHistogramOptions(ffDiEMSize, "No. selected jets per event", "", "");

   //MET of events with fewer than 2 matching jets
   TH1F ggMETUnmatched("ggMETUnmatched", "ggMETUnmatched", nMETBins, METBins);
   TH1F ffMETUnmatched("ffMETUnmatched", "ffMETUnmatched", nMETBins, METBins);
   setHistogramOptions(ggMETUnmatched, "ME_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffMETUnmatched, "ME_{T} (GeV)", "", "", kSumW2);

   //R2 vs. MR histograms
   TH2F ggR2VsMR("ggR2VsMR", "", nMRBins, MRBins, nR2Bins, R2Bins);
   TH2F egR2VsMR("egR2VsMR", "", nMRBins, MRBins, nR2Bins, R2Bins);
   TH2F eeR2VsMR("eeR2VsMR", "", nMRBins, MRBins, nR2Bins, R2Bins);
   TH2F eeLowSidebandR2VsMR("eeLowSidebandR2VsMR", "", nMRBins, MRBins, nR2Bins, R2Bins);
   TH2F eeHighSidebandR2VsMR("eeHighSidebandR2VsMR", "", nMRBins, MRBins, nR2Bins, R2Bins);
   TH2F ffR2VsMR("ffR2VsMR", "", nMRBins, MRBins, nR2Bins, R2Bins);
   setHistogramOptions(ggR2VsMR, "M_{R} (GeV)", "R^{2}", "", kSumW2);
   setHistogramOptions(egR2VsMR, "M_{R} (GeV)", "R^{2}", "", kSumW2);
   setHistogramOptions(eeR2VsMR, "M_{R} (GeV)", "R^{2}", "", kSumW2);
   setHistogramOptions(eeLowSidebandR2VsMR, "M_{R} (GeV)", "R^{2}", "", kSumW2);
   setHistogramOptions(eeHighSidebandR2VsMR, "M_{R} (GeV)", "R^{2}", "", kSumW2);
   setHistogramOptions(ffR2VsMR, "M_{R} (GeV)", "R^{2}", "", kSumW2);

   //MET vs. di-EM ET vs. invariant mass histograms
   //need an additional 3 copies of the gg, ee, and ff histograms for the 3 Nj bins
   vector<TH3F*> ggMETVsDiEMETVsInvMass;
   vector<TH3F*> ggMETVsDijetETVsInvMass;
   vector<TH3F*> eeMETVsDiEMETVsInvMass;
   vector<TH3F*> eeLowSidebandMETVsDiEMETVsInvMass;
   vector<TH3F*> eeHighSidebandMETVsDiEMETVsInvMass;
   vector<TH3F*> ffMETVsDiEMETVsInvMass;
   for (unsigned int iNjBin = 1; iNjBin <= nNjBins; ++iNjBin) {
     const Double_t* pDiEMETBins = diEMETBins1j;
     unsigned int pNDiEMETBins = nDiEMETBins1j;
     if (iNjBin == 1) {
       pDiEMETBins = diEMETBins0j;
       pNDiEMETBins = nDiEMETBins0j;
     }
     ggMETVsDiEMETVsInvMass.
       push_back(new TH3F(histName("ggMETVsDiEMETVsInvMass_", "NjBin", iNjBin).c_str(), "", 
			  nInvMassBins, invMassBins, pNDiEMETBins, pDiEMETBins, nMETBins, 
			  METBins));
     ggMETVsDijetETVsInvMass.
       push_back(new TH3F(histName("ggMETVsDijetETVsInvMass_", "NjBin", iNjBin).c_str(), "", 
			  nInvMassBins, invMassBins, pNDiEMETBins, pDiEMETBins, nMETBins, 
			  METBins));
     eeMETVsDiEMETVsInvMass.
       push_back(new TH3F(histName("eeMETVsDiEMETVsInvMass_", "NjBin", iNjBin).c_str(), "", 
			  nInvMassBins, invMassBins, pNDiEMETBins, pDiEMETBins, nMETBins, 
			  METBins));
     eeLowSidebandMETVsDiEMETVsInvMass.
       push_back(new TH3F(histName("eeLowSidebandMETVsDiEMETVsInvMass_", "NjBin", iNjBin).c_str(), 
			  "", nInvMassBins, invMassBins, pNDiEMETBins, pDiEMETBins, nMETBins, 
			  METBins));
     eeHighSidebandMETVsDiEMETVsInvMass.
       push_back(new TH3F(histName("eeHighSidebandMETVsDiEMETVsInvMass_", "NjBin", iNjBin).
			  c_str(), "", nInvMassBins, invMassBins, pNDiEMETBins, pDiEMETBins, 
			  nMETBins, METBins));
     ffMETVsDiEMETVsInvMass.
       push_back(new TH3F(histName("ffMETVsDiEMETVsInvMass_", "NjBin", iNjBin).c_str(), "", 
			  nInvMassBins, invMassBins, pNDiEMETBins, pDiEMETBins, nMETBins, 
			  METBins));
     setHistogramOptions(*ggMETVsDiEMETVsInvMass[ggMETVsDiEMETVsInvMass.size() - 1], 
			 "m_{#gamma#gamma} (GeV)", "Di-EM E_{T} (GeV)", "ME_{T} (GeV)", kSumW2);
     setHistogramOptions(*ggMETVsDijetETVsInvMass[ggMETVsDijetETVsInvMass.size() - 1], 
			 "m_{#gamma#gamma} (GeV)", "Dijet E_{T} (GeV)", "ME_{T} (GeV)", kSumW2);
     setHistogramOptions(*eeMETVsDiEMETVsInvMass[eeMETVsDiEMETVsInvMass.size() - 1], 
			 "m_{ee} (GeV)", "Di-EM E_{T} (GeV)", "ME_{T} (GeV)", kSumW2);
     setHistogramOptions(*eeLowSidebandMETVsDiEMETVsInvMass[eeLowSidebandMETVsDiEMETVsInvMass.
							   size() - 1], "m_{ee} (GeV)", 
			 "Di-EM E_{T} (GeV)", "ME_{T} (GeV)", kSumW2);
     setHistogramOptions(*eeHighSidebandMETVsDiEMETVsInvMass[eeHighSidebandMETVsDiEMETVsInvMass.
							    size() - 1], "m_{ee} (GeV)", 
			 "Di-EM E_{T} (GeV)", "ME_{T} (GeV)", kSumW2);
     setHistogramOptions(*ffMETVsDiEMETVsInvMass[ffMETVsDiEMETVsInvMass.size() - 1], 
			 "m_{ff} (GeV)", "Di-EM E_{T} (GeV)", "ME_{T} (GeV)", kSumW2);
   }
   TH3F ggMETVsDiEMETVsInvMassTot("ggMETVsDiEMETVsInvMassTot", "", nInvMassBins, invMassBins, 
				  nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F ggMETVsDijetETVsInvMassTot("ggMETVsDijetETVsInvMassTot", "", nInvMassBins, invMassBins, 
				   nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F egMETVsDiEMETVsInvMass("egMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
			       nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F eeMETVsDiEMETVsInvMassTot("eeMETVsDiEMETVsInvMassTot", "", nInvMassBins, invMassBins, 
				  nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F eeLowSidebandMETVsDiEMETVsInvMassTot("eeLowSidebandMETVsDiEMETVsInvMassTot", "", 
					     nInvMassBins, invMassBins, nDiEMETBins, diEMETBins, 
					     nMETBins, METBins);
   TH3F eeHighSidebandMETVsDiEMETVsInvMassTot("eeHighSidebandMETVsDiEMETVsInvMassTot", "", 
					      nInvMassBins, invMassBins, nDiEMETBins, diEMETBins, 
					      nMETBins, METBins);
   TH3F ffMETVsDiEMETVsInvMassTot("ffMETVsDiEMETVsInvMassTot", "", nInvMassBins, invMassBins, 
				  nDiEMETBins, diEMETBins, nMETBins, METBins);
   TH3F eeffMETVsDiEMETVsInvMass("eeffMETVsDiEMETVsInvMass", "", nInvMassBins, invMassBins, 
				 nDiEMETBins, diEMETBins, nMETBins, METBins);
   setHistogramOptions(ggMETVsDiEMETVsInvMassTot, "m_{#gamma#gamma} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(ggMETVsDijetETVsInvMassTot, "m_{#gamma#gamma} (GeV)", "Dijet E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(egMETVsDiEMETVsInvMass, "m_{e#gamma} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(eeMETVsDiEMETVsInvMassTot, "m_{ee} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(eeLowSidebandMETVsDiEMETVsInvMassTot, "m_{ee} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(eeHighSidebandMETVsDiEMETVsInvMassTot, "m_{ee} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(ffMETVsDiEMETVsInvMassTot, "m_{ff} (GeV)", "Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);
   setHistogramOptions(eeffMETVsDiEMETVsInvMass, "m_{ee} (GeV)", "ff Di-EM E_{T} (GeV)", 
		       "ME_{T} (GeV)", kSumW2);

   //reweighted and normalized ee and ff R2 vs. MR histograms
   TH2F eeR2VsMRFinal("eeR2VsMRFinal", "", nMRBins, MRBins, nR2Bins, R2Bins);
   TH2F eeLowSidebandR2VsMRFinal("eeLowSidebandR2VsMRFinal", "", nMRBins, MRBins, nR2Bins, R2Bins);
   TH2F eeHighSidebandR2VsMRFinal("eeHighSidebandR2VsMRFinal", "", nMRBins, MRBins, nR2Bins, 
				  R2Bins);
   TH2F ffR2VsMRFinal("ffR2VsMRFinal", "", nMRBins, MRBins, nR2Bins, R2Bins);
   setHistogramOptions(eeR2VsMRFinal, "M_{R} (GeV)", "R^{2}", "", kSumW2);
   setHistogramOptions(eeLowSidebandR2VsMRFinal, "M_{R} (GeV)", "R^{2}", "", kSumW2);
   setHistogramOptions(eeHighSidebandR2VsMRFinal, "M_{R} (GeV)", "R^{2}", "", kSumW2);
   setHistogramOptions(ffR2VsMRFinal, "M_{R} (GeV)", "R^{2}", "", kSumW2);

   //reweighted and normalized ee and ff MET histograms
   TH1F eeFinal("eeFinal", "", nMETBins, METBins);
   TH1F eeLowSidebandFinal("eeLowSidebandFinal", "", nMETBins, METBins);
   TH1F eeHighSidebandFinal("eeHighSidebandFinal", "", nMETBins, METBins);
   TH1F ffFinal("ffFinal", "", nMETBins, METBins);
   setHistogramOptions(eeFinal, "ME_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeLowSidebandFinal, "ME_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeHighSidebandFinal, "ME_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffFinal, "ME_{T} (GeV)", "", "", kSumW2);

   //clean HT histograms
   TH1F ggCleanHT("ggCleanHT", "", 60, 0.0, 600.0);
   TH1F egCleanHT("egCleanHT", "", 60, 0.0, 600.0);
   TH1F eeCleanHT("eeCleanHT", "", 60, 0.0, 600.0);
   TH1F ffCleanHT("ffCleanHT", "", 60, 0.0, 600.0);
   setHistogramOptions(ggCleanHT, "H_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(egCleanHT, "H_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeCleanHT, "H_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffCleanHT, "H_{T} (GeV)", "", "", kSumW2);

   //clean MHT histograms
   TH1F ggCleanMHT("ggCleanMHT", "", 60, 0.0, 600.0);
   TH1F egCleanMHT("egCleanMHT", "", 60, 0.0, 600.0);
   TH1F eeCleanMHT("eeCleanMHT", "", 60, 0.0, 600.0);
   TH1F ffCleanMHT("ffCleanMHT", "", 60, 0.0, 600.0);
   setHistogramOptions(ggCleanMHT, "MH_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(egCleanMHT, "MH_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeCleanMHT, "MH_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffCleanMHT, "MH_{T} (GeV)", "", "", kSumW2);

   //clean Nj histograms
   TH1F ggCleanNj("ggCleanNj", "", 10, -0.5, 9.5);
   TH1F egCleanNj("egCleanNj", "", 10, -0.5, 9.5);
   TH1F eeCleanNj("eeCleanNj", "", 10, -0.5, 9.5);
   TH1F ffCleanNj("ffCleanNj", "", 10, -0.5, 9.5);
   setHistogramOptions(ggCleanNj, "N_{j}", "", "", kSumW2);
   setHistogramOptions(egCleanNj, "N_{j}", "", "", kSumW2);
   setHistogramOptions(eeCleanNj, "N_{j}", "", "", kSumW2);
   setHistogramOptions(ffCleanNj, "N_{j}", "", "", kSumW2);

   //clean leading jet ET histograms
   TH1F ggCleanLeadingJetET("ggCleanLeadingJetET", "", 50, 0.0, 500.0);
   TH1F egCleanLeadingJetET("egCleanLeadingJetET", "", 50, 0.0, 500.0);
   TH1F eeCleanLeadingJetET("eeCleanLeadingJetET", "", 50, 0.0, 500.0);
   TH1F ffCleanLeadingJetET("ffCleanLeadingJetET", "", 50, 0.0, 500.0);
   setHistogramOptions(ggCleanLeadingJetET, "E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(egCleanLeadingJetET, "E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(eeCleanLeadingJetET, "E_{T} (GeV)", "", "", kSumW2);
   setHistogramOptions(ffCleanLeadingJetET, "E_{T} (GeV)", "", "", kSumW2);

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

//    /*histogram of the percentage difference between Poisson mean of generated di-EM ET weights and 
//      input mean (i.e. the value of the measured weight)*/
//    TH1F ffDiEMETInputVsOutput("ffDiEMETInputVsOutput", "", nDiEMETBins, diEMETBins);
//    TH1F eeDiEMETInputVsOutput("eeDiEMETInputVsOutput", "", nDiEMETBins, diEMETBins);
//    TH1F eeLowSidebandDiEMETInputVsOutput("eeLowSidebandDiEMETInputVsOutput", "", nDiEMETBins, 
// 					 diEMETBins);
//    TH1F eeHighSidebandDiEMETInputVsOutput("eeHighSidebandDiEMETInputVsOutput", "", nDiEMETBins, 
// 					  diEMETBins);
//    setHistogramOptions(ffDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
// 		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);
//    setHistogramOptions(eeDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
// 		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);
//    setHistogramOptions(eeLowSidebandDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
// 		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);
//    setHistogramOptions(eeHighSidebandDiEMETInputVsOutput, "Di-EM E_{T} (GeV)", 
// 		       "#frac{#mu_{in} - #mu_{out}}{#mu_{in}}", "", kSumW2);

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

   //canvases for the final plot
   TCanvas METCanvas("METCanvas", "", 600, 600);
   TCanvas MRCanvas("MRCanvas", "", 600, 600);
   TCanvas R2Canvas("R2Canvas", "", 600, 600);
   setCanvasOptions(METCanvas, "ME_{T} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);
   setCanvasOptions(MRCanvas, "M_{R} (GeV)", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);
   setCanvasOptions(R2Canvas, "R^{2}", "", 0.0, 500.0, 0.0, 10000.0, kSetGrid);

   //MET, di-EM ET, and Nj vectors
   VFLOAT eeMETVec;
   VFLOAT eeLowSidebandMETVec;
   VFLOAT eeHighSidebandMETVec;
   VFLOAT ffMETVec;
   VFLOAT eeDiEMETVec;
   VFLOAT eeLowSidebandDiEMETVec;
   VFLOAT eeHighSidebandDiEMETVec;
   VFLOAT ffDiEMETVec;
   VFLOAT eeNjVec;
   VFLOAT eeLowSidebandNjVec;
   VFLOAT eeHighSidebandNjVec;
   VFLOAT ffNjVec;

   //vectors for the JES systematic
   vector<TLorentzVector> ggJ1UncorrP4;
   vector<TLorentzVector> ggJ2UncorrP4;
   VFLOAT ggJ1JES;
   VFLOAT ggJ2JES;
   VFLOAT ggJ1JESErr;
   VFLOAT ggJ2JESErr;
   vector<TLorentzVector> ffJ1UncorrP4;
   vector<TLorentzVector> ffJ2UncorrP4;
   VFLOAT ffJ1JES;
   VFLOAT ffJ2JES;
   VFLOAT ffJ1JESErr;
   VFLOAT ffJ2JESErr;

   //MR and R2 vectors
   VFLOAT eeMRVec;
   VFLOAT eeLowSidebandMRVec;
   VFLOAT eeHighSidebandMRVec;
   VFLOAT ffMRVec;
   VFLOAT eeR2Vec;
   VFLOAT eeLowSidebandR2Vec;
   VFLOAT eeHighSidebandR2Vec;
   VFLOAT ffR2Vec;

   //lumi and PU weight vectors
   VFLOAT ggLumiWeightVec;
   VFLOAT egLumiWeightVec;
   VFLOAT eeLumiWeightVec;
   VFLOAT eeLowSidebandLumiWeightVec;
   VFLOAT eeHighSidebandLumiWeightVec;
   VFLOAT ffLumiWeightVec;
   VFLOAT ggPUWeightVec;
   VFLOAT egPUWeightVec;
   VFLOAT eePUWeightVec;
   VFLOAT eeLowSidebandPUWeightVec;
   VFLOAT eeHighSidebandPUWeightVec;
   VFLOAT ffPUWeightVec;

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

   //set up on-the-fly JEC uncertainties
   JetCorrectionUncertainty JECErr(JECErrFile_);

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
	 (dataset.Contains("v2_TT"))) dataset.Remove(0, 3); //for TT sample
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
     int nPV = -1;
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

     //apply user trigger requirement and ECAL/HCAL filters
     if (/*passUserHLT() && (susyEvent->PassesHcalNoiseFilter == 1) && 
	   (susyEvent->PassesEcalDeadCellFilter == 1) || */true) {

       //fill MET vs. di-EM ET, invariant mass, rho, leading/trailing photon ET, and nPV histograms
       double invMass = -1.0;
       if (susyCategory->getEvtInvMass()->size() > 0) invMass = susyCategory->getEvtInvMass(tag_);
       else {
	 cerr << "Error: susyCategory->getEvtInvMass()->size() <= 0 in event " << (jentry + 1);
	 cerr << ".  Using invariant mass -1.0 GeV.\n";
       }
       double MET = -1.0;
       TVector2 MET2DVec;
       map<TString, susy::MET>::const_iterator iMET = susyEvent->metMap.find("pfMet");
       if (iMET != susyEvent->metMap.end()) {
	 MET = iMET->second.met();
	 MET2DVec = iMET->second.mEt;
       }
       else cerr << "Error: PFMET not found in event " << (jentry + 1) << ".\n";
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

       //get PF electron multiplicity
       unsigned int nPFElectrons = 0;
       map<TString, susy::ElectronCollection>::const_iterator iElectronMap =
	 susyEvent->electrons.find("gsfElectrons");
       if (iElectronMap != susyEvent->electrons.end()) {
	 for (susy::ElectronCollection::const_iterator iElectron = iElectronMap->second.begin();
	      iElectron != iElectronMap->second.end(); ++iElectron) {
	   if (iElectron->isPF()) ++nPFElectrons;
	 }
       }
 
       //get PF muon multiplicity
       unsigned int nMuons = 0;
       for (vector<susy::Muon>::const_iterator iMuon = susyEvent->muons.begin();
	    iMuon != susyEvent->muons.end(); ++iMuon) {
	 map<TString, UChar_t>::const_iterator iID =
	   iMuon->idPairs.find("muidGlobalMuonPromptTight");
	 if ((iID != iMuon->idPairs.end()) && ((unsigned int)iID->second == 1)) ++nMuons;
       }

       //count jets
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
	       cerr << "Error: susyCategory->getIsDeciding()->size() <= 0 in event ";
	       cerr << (jentry + 1) << ".  Assuming no deciding photons.\n";
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
	   if ((deltaR(corrP4.Eta(), corrP4.Phi(), EM1Eta, EM1Phi) > 0.5) && 
	       (deltaR(corrP4.Eta(), corrP4.Phi(), EM2Eta, EM2Phi) > 0.5) && 
	       (fabs(corrP4.Eta()) <= 2.6)) {
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
	       jetET.push_back(corrP4.Et());
	       if (corrP4.Et() >= 30.0/*GeV*/) ++nJets30;
	       if (corrP4.Et() >= 60.0/*GeV*/) ++nJets60;
	     }
	   }
	 }
       }
       else {
	 cerr << "Error: " << tag_ << " photon collection or ak5 jet collection not found in ";
	 cerr << "event " << (jentry + 1) << ".\n";
       }
       const float oldHT = accumulate(jetET.begin(), jetET.end(), 0);

       /*HT, MHT, Nj, and leading jet ET from function (cleaned of leptons and the 2 primary EM 
	 objects, |eta| <= 5)*/
       const float cleanHT = HT();
       const float cleanMHT = MHT();
       const float cleanNj = numJets(5.0);
       const float cleanLeadingJetET = leadingJetET();

       //get combined isolation, photon ET, sigmaIetaIeta, R9, and ECAL/HCAL/track isolation
       vector<float> combinedIso;
       vector<float> ET;
       vector<float> sigmaIetaIeta;
       vector<float> R9;
       vector<float> ECALIso;
       vector<float> HCALIso;
       vector<float> trackIso;
       vector<TLorentzVector> photonP4;
       vector<TLorentzVector> jP4;
       vector<TLorentzVector> e;
       vector<TLorentzVector> f;
       vector<float> eET;
       vector<float> fET;
       vector<TLorentzVector> uncorrP4;
       vector<float> JES;
       vector<float> JESErr;
       if (evtCategory != FAIL) {
	 if (iPhotonMap != susyEvent->photons.end()) {
	   for (susy::PhotonCollection::const_iterator iPhoton = iPhotonMap->second.begin(); 
		iPhoton != iPhotonMap->second.end(); ++iPhoton) {
	     const unsigned int photonIndex = iPhoton - iPhotonMap->second.begin();
	     if (susyCategory->getIsDeciding()->size() > 0) {
	       if (susyCategory->getIsDeciding(tag_, photonIndex)) {
		 ET.push_back(iPhoton->momentum.Et());
		 combinedIso.push_back(iPhoton->ecalRecHitSumEtConeDR03 - 0.1474*rho + 
				       iPhoton->hcalTowerSumEtConeDR03() - 0.0467*rho + 
				       iPhoton->trkSumPtHollowConeDR03);
		 sigmaIetaIeta.push_back(iPhoton->sigmaIetaIeta);
		 R9.push_back(iPhoton->r9);
		 ECALIso.push_back(iPhoton->ecalRecHitSumEtConeDR03);
		 HCALIso.push_back(iPhoton->hcalTowerSumEtConeDR03());
		 trackIso.push_back(iPhoton->trkSumPtHollowConeDR03);
		 photonP4.push_back(iPhoton->momentum);

		 //get matching jet corrected p4
		 for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
		      iJet != iJets->second.end(); ++iJet) {
		   
		   //JES is L1Fast scale only...
		   float theJES = 1.0;
		   TLorentzVector corrP4 = 
		     correctedJet4Momentum(iJet, "L1FastL2L3", theJES, "L2L3");

		   if (passJetFiducialCuts(corrP4, 10.0/*GeV*/, 2.6)/* && passPFJetID(iJet)*/ && 
		       (deltaR(iJet->momentum.Eta(), iJet->momentum.Phi(), 
			       iPhoton->caloPosition.Eta(), 
			       iPhoton->caloPosition.Phi()) < 0.3)) {
		     uncorrP4.push_back(iJet->momentum);
		     jP4.push_back(corrP4);
		     JES.push_back(theJES);

		     //...but JES error is for L1FastL2L3(Residual) cf. Mikko Voutilanen
		     JECErr.setJetEta(corrP4.Eta());
		     JECErr.setJetPt(corrP4.Pt());
		     JESErr.push_back(JECErr.getUncertainty(true));
		   }
		 }
	       }

	       /*if not a deciding photon in a gg, eg, ee, or ff event, count if it's an e or f 
		 for eeff purposes*/
// 	       else {
// 		 if (susyCategory->getPhotonType(tag_, photonIndex) == E) {
// 		   e.push_back(iPhoton->momentum);
// 		   eET.push_back(iPhoton->momentum.Et());
// 		 }
// 		 if (susyCategory->getPhotonType(tag_, photonIndex) == F) {
// 		   f.push_back(iPhoton->momentum);
// 		   fET.push_back(iPhoton->momentum.Et());
// 		 }
// 	       }
	     }
	     else {
	       cerr << "Error: susyCategory->getIsDeciding()->size() <= 0 in event ";
	       cerr << (jentry + 1) << ".  Assuming no deciding photons.\n";
	     }
	   }
	   if (combinedIso.size() != 2) {
	     cerr << "Error: " << combinedIso.size() << " selected photons in event ";
	     cerr << (jentry + 1) << ".\n";
	   }
	 }
	 else {
	   cerr << "Error: " << tag_ << " photon collection not found in event " << (jentry + 1);
	   cerr << ".\n";
	 }
       }

       //sort photon ET in ascending order
       sort(ET.begin(), ET.end());

       //get di-EM and dijet pT and razor variables
       double dijetPT = -1.0;
       double diEMET = -1.0;
       double MR = -1.0;
       double R2 = -1.0;
       const unsigned int dijetSize = jP4.size();
       const unsigned int diEMSize = photonP4.size();
       if (dijetSize == 2) dijetPT = (jP4[0] + jP4[1]).Pt();
       if (diEMSize == 2) {
	 diEMET = (photonP4[0] + photonP4[1]).Pt();
	 const double MR2 = 
	   (photonP4[0].Energy() + photonP4[1].Energy())*
	   (photonP4[0].Energy() + photonP4[1].Energy()) - 
	   (photonP4[0].Pz() + photonP4[1].Pz())*(photonP4[0].Pz() + photonP4[1].Pz());
	 R2 = ((MET*(photonP4[0].Pt() + photonP4[1].Pt()) - 
		MET2DVec*((photonP4[0] + photonP4[1]).Vect().XYvector()))/2.0)/MR2;
	 MR = sqrt(MR2);
       }

       //determine Nj bin
       unsigned int iNjBin = nNjBins;
       unsigned int nJets = numJets(2.6);
       for (unsigned int iNjBinBoundary = 0; iNjBinBoundary < nNjBins; ++iNjBinBoundary) {
	 if ((nJets >= NjBins[iNjBinBoundary]) && 
	     (nJets < NjBins[iNjBinBoundary + 1])) iNjBin = iNjBinBoundary;
       }
       if (iNjBin == nNjBins) {
	 cerr << "Error: no Nj bin found for nJets = " << nJets << " in event ";
	 cerr << (jentry + 1) << ".  Assuming nJets = 0.\n";
	 iNjBin = 0;
       }

       //get pThat
       float pTHat = -1.0;
       map<TString, Float_t>::const_iterator iPTHat = susyEvent->gridParams.find("ptHat");
       if (iPTHat != susyEvent->gridParams.end()) pTHat = iPTHat->second;

       switch (evtCategory) {
       case GG:
	 if (datasetString.find("DiPhoton") != string::npos) {

	   //fill lumi and PU weight vectors
	   ggLumiWeightVec.push_back(lumiWeight);
	   ggPUWeightVec.push_back(PUWeight);

	   for (vector<float>::const_iterator i = combinedIso.begin(); i != combinedIso.end(); 
		++i) { ggCombinedIso.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = sigmaIetaIeta.begin(); i != sigmaIetaIeta.end(); 
		++i) { ggSigmaIetaIeta.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = R9.begin(); i != R9.end(); 
		++i) { ggR9.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = ECALIso.begin(); i != ECALIso.end(); 
		++i) { ggECALIso.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = HCALIso.begin(); i != HCALIso.end(); 
		++i) { ggHCALIso.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = trackIso.begin(); i != trackIso.end(); 
		++i) { ggTrackIso.Fill(*i, lumiWeight*PUWeight); }
	   ggLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	   ggTrailingPhotonET.Fill(*(ET.begin()), lumiWeight*PUWeight);
	   ggNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	   ggRho.Fill(rho, lumiWeight*PUWeight);
	   ggPFElectronMultiplicity.Fill(nPFElectrons, lumiWeight*PUWeight);
	   ggMuonMultiplicity.Fill(nMuons, lumiWeight*PUWeight);
	   ggDijetSize.Fill(dijetSize);
	   ggDiEMSize.Fill(diEMSize);
	   ggPTHat.Fill(pTHat, lumiWeight*PUWeight);
	   ggCleanHT.Fill(cleanHT, lumiWeight*PUWeight);
	   ggCleanMHT.Fill(cleanMHT, lumiWeight*PUWeight);
	   ggCleanNj.Fill(cleanNj, lumiWeight*PUWeight);
	   ggCleanLeadingJetET.Fill(cleanLeadingJetET, lumiWeight*PUWeight);
	   if (dijetSize == 2) {
	     ggMETVsDijetETVsInvMass[iNjBin]->Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
	     ggMETVsDijetETVsInvMassTot.Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
	     ggDijetETUniform[iNjBin]->Fill(dijetPT, lumiWeight*PUWeight);
	     ggJ1UncorrP4.push_back(uncorrP4[0]);
	     ggJ2UncorrP4.push_back(uncorrP4[1]);
	     ggJ1JES.push_back(JES[0]);
	     ggJ2JES.push_back(JES[1]);
	     ggJ1JESErr.push_back(JESErr[0]);
	     ggJ2JESErr.push_back(JESErr[1]);
	   }
	   else {
	     // 	   printDiObjectErrorMessage(dijetSize, "dijet", evtCategory, jentry);
	     ggMETUnmatched.Fill(MET, lumiWeight*PUWeight);
	     
	     //use di-EM ET when dijet ET is unavailable
	     if (diEMSize == 2) {
	       ggMETVsDijetETVsInvMass[iNjBin]->Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	       ggMETVsDijetETVsInvMassTot.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	     }
	     else printDiObjectErrorMessage(diEMSize, "diEM", evtCategory, jentry);
	   }
	   if (diEMSize == 2) {
	     ggMETVsDiEMETVsInvMass[iNjBin]->Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	     ggMETVsDiEMETVsInvMassTot.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	     ggDiEMETUniform[iNjBin]->Fill(diEMET, lumiWeight*PUWeight);
	     ggR2VsMR.Fill(MR, R2, lumiWeight*PUWeight);
	   }
	   else printDiObjectErrorMessage(diEMSize, "diEM", evtCategory, jentry);
	   
	   //print the event information to a file
	   ggEvtFile << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
	   ggEvtFile << susyEvent->luminosityBlockNumber << endl;
	   
	 }
	 break;
       case EG:

	 //fill lumi and PU weight vectors
	 egLumiWeightVec.push_back(lumiWeight);
	 egPUWeightVec.push_back(PUWeight);

	 for (vector<float>::const_iterator i = combinedIso.begin(); i != combinedIso.end(); 
	      ++i) { egCombinedIso.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = sigmaIetaIeta.begin(); i != sigmaIetaIeta.end(); 
	      ++i) { egSigmaIetaIeta.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = R9.begin(); i != R9.end(); 
	      ++i) { egR9.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = ECALIso.begin(); i != ECALIso.end(); 
	      ++i) { egECALIso.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = HCALIso.begin(); i != HCALIso.end(); 
	      ++i) { egHCALIso.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = trackIso.begin(); i != trackIso.end(); 
	      ++i) { egTrackIso.Fill(*i, lumiWeight*PUWeight); }
	 egLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 egTrailingPhotonET.Fill(*(ET.begin()), lumiWeight*PUWeight);
	 egNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 egRho.Fill(rho, lumiWeight*PUWeight);
	 egPFElectronMultiplicity.Fill(nPFElectrons, lumiWeight*PUWeight);
	 egMuonMultiplicity.Fill(nMuons, lumiWeight*PUWeight);
	 egDijetSize.Fill(dijetSize);
	 egDiEMSize.Fill(diEMSize);
	 egPTHat.Fill(pTHat, lumiWeight*PUWeight);
	 egCleanHT.Fill(cleanHT, lumiWeight*PUWeight);
	 egCleanMHT.Fill(cleanMHT, lumiWeight*PUWeight);
	 egCleanNj.Fill(cleanNj, lumiWeight*PUWeight);
	 egCleanLeadingJetET.Fill(cleanLeadingJetET, lumiWeight*PUWeight);
	 if (diEMSize == 2) {
	   egMETVsDiEMETVsInvMass.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	   egR2VsMR.Fill(MR, R2, lumiWeight*PUWeight);
	 }
	 else printDiObjectErrorMessage(diEMSize, "diEM", evtCategory, jentry);

	 //print the event information to a file
	 egEvtFile << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
	 egEvtFile << susyEvent->luminosityBlockNumber << endl;

	 break;
       case EE:
	 if ((datasetString.find("DYToEE") != string::npos) || 
	     (datasetString.find("TT") != string::npos)) {
	   for (vector<float>::const_iterator i = combinedIso.begin(); i != combinedIso.end(); 
		++i) { eeCombinedIso.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = sigmaIetaIeta.begin(); i != sigmaIetaIeta.end(); 
		++i) { eeSigmaIetaIeta.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = R9.begin(); i != R9.end(); 
		++i) { eeR9.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = ECALIso.begin(); i != ECALIso.end(); 
		++i) { eeECALIso.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = HCALIso.begin(); i != HCALIso.end(); 
		++i) { eeHCALIso.Fill(*i, lumiWeight*PUWeight); }
	   for (vector<float>::const_iterator i = trackIso.begin(); i != trackIso.end(); 
		++i) { eeTrackIso.Fill(*i, lumiWeight*PUWeight); }
	   eeLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	   eeTrailingPhotonET.Fill(*(ET.begin()), lumiWeight*PUWeight);
	   eeDijetSize.Fill(dijetSize);
	   eeDiEMSize.Fill(diEMSize);
	   eePTHat.Fill(pTHat, lumiWeight*PUWeight);
	   // 	 if (f.size() >= 2) {
	   // 	   vector<float>::iterator i1 = max_element(fET.begin(), fET.end());
	   // 	   const unsigned int index1 = i1 - fET.begin();
	   // 	   fET.erase(i1);
	   // 	   const unsigned int index2 = max_element(fET.begin(), fET.end()) - fET.begin();
	   // 	   eeffMETVsDiEMETVsInvMass.Fill(invMass, (f[index1] + f[index2]).Pt(), MET, 
	   // 					 lumiWeight*PUWeight);
	   // 	 }
	   if (diEMSize == 2) {
	     if ((invMass >= 71.0/*GeV*/) && (invMass < 81.0/*GeV*/)) { /*new sidebands from 
									  Yueh-Feng 6-Jan-12*/
	       //fill lumi and PU weight vectors
	       eeLowSidebandLumiWeightVec.push_back(lumiWeight);
	       eeLowSidebandPUWeightVec.push_back(PUWeight);

	       if (dijetSize == 2) {
		 eeLowSidebandMETVsDiEMETVsInvMass[iNjBin]->
		   Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
		 eeLowSidebandMETVsDiEMETVsInvMassTot.
		   Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
		 eeLowSidebandDiEMETVec.push_back(dijetPT);
	       }
	       else {
		 eeLowSidebandMETVsDiEMETVsInvMass[iNjBin]->
		   Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
		 eeLowSidebandMETVsDiEMETVsInvMassTot.
		   Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
		 eeLowSidebandDiEMETVec.push_back(diEMET);
	       }
	       eeLowSidebandMETVec.push_back(MET);
	       eeLowSidebandNjVec.push_back(nJets);
	       eeLowSidebandMRVec.push_back(MR);
	       eeLowSidebandR2Vec.push_back(R2);
	     }
	     if ((invMass >= 81.0/*GeV*/) && (invMass < 101.0/*GeV*/)) {

	       //fill lumi and PU weight vectors
	       eeLumiWeightVec.push_back(lumiWeight);
	       eePUWeightVec.push_back(PUWeight);

	       //fill HT, MHT, Nj, and leading jet ET histograms for lepton/EM cleaned jets
	       eeCleanHT.Fill(cleanHT, lumiWeight*PUWeight);
	       eeCleanMHT.Fill(cleanMHT, lumiWeight*PUWeight);
	       eeCleanNj.Fill(cleanNj, lumiWeight*PUWeight);
	       eeCleanLeadingJetET.Fill(cleanLeadingJetET, lumiWeight*PUWeight);

	       //fill rho and nPV histograms
	       eeRho.Fill(rho, lumiWeight*PUWeight);
	       eeNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);

	       //fill PF electron and muon multiplicity histograms
	       eePFElectronMultiplicity.Fill(nPFElectrons, lumiWeight*PUWeight);
	       eeMuonMultiplicity.Fill(nMuons, lumiWeight*PUWeight);

	       if (dijetSize == 2) {
		 eeMETVsDiEMETVsInvMass[iNjBin]->
		   Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
		 eeMETVsDiEMETVsInvMassTot.
		   Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
		 eeDiEMETVec.push_back(dijetPT);
		 eeDijetETUniform[iNjBin]->Fill(dijetPT, lumiWeight*PUWeight);
	       }
	       else {
		 eeMETVsDiEMETVsInvMass[iNjBin]->
		   Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
		 eeMETVsDiEMETVsInvMassTot.
		   Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
		 eeDiEMETVec.push_back(diEMET);
	       }
	       eeR2VsMR.Fill(MR, R2, lumiWeight*PUWeight);
	       eeMETVec.push_back(MET);
	       eeNjVec.push_back(nJets);
	       eeMRVec.push_back(MR);
	       eeR2Vec.push_back(R2);
	       eeHTVsMET.Fill(MET, oldHT, lumiWeight*PUWeight);
	       eeMETVsNJets30.Fill(nJets30, MET, lumiWeight*PUWeight);
	       eeMETVsNJets60.Fill(nJets60, MET, lumiWeight*PUWeight);
	       eeDiEMETUniform[iNjBin]->Fill(diEMET, lumiWeight*PUWeight);
	     }
	     if ((invMass >= 101.0/*GeV*/) && (invMass < 111.0/*GeV*/)) { /*new sidebands from 
									    Yueh-Feng 6-Jan-12*/
	       //fill lumi and PU weight vectors
	       eeHighSidebandLumiWeightVec.push_back(lumiWeight);
	       eeHighSidebandPUWeightVec.push_back(PUWeight);

	       if (dijetSize == 2) {
		 eeHighSidebandMETVsDiEMETVsInvMass[iNjBin]->
		   Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
		 eeHighSidebandMETVsDiEMETVsInvMassTot.
		   Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
		 eeHighSidebandDiEMETVec.push_back(dijetPT);
	       }
	       else {
		 eeHighSidebandMETVsDiEMETVsInvMass[iNjBin]->
		   Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
		 eeHighSidebandMETVsDiEMETVsInvMassTot.
		   Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
		 eeHighSidebandDiEMETVec.push_back(diEMET);
	       }
	       eeHighSidebandMETVec.push_back(MET);
	       eeHighSidebandNjVec.push_back(nJets);
	       eeHighSidebandMRVec.push_back(MR);
	       eeHighSidebandR2Vec.push_back(R2);
	     }
	   }
	   else printDiObjectErrorMessage(diEMSize, "diEM", evtCategory, jentry);

	   //print the event information to a file
	   eeEvtFile << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
	   eeEvtFile << susyEvent->luminosityBlockNumber << " " << invMass << endl;
	   
	 }
	 break;
       case FF:

	 //fill lumi and PU weight vectors
	 ffLumiWeightVec.push_back(lumiWeight);
	 ffPUWeightVec.push_back(PUWeight);

	 for (vector<float>::const_iterator i = combinedIso.begin(); i != combinedIso.end(); 
	      ++i) { ffCombinedIso.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = sigmaIetaIeta.begin(); i != sigmaIetaIeta.end(); 
	      ++i) { ffSigmaIetaIeta.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = R9.begin(); i != R9.end(); 
	      ++i) { ffR9.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = ECALIso.begin(); i != ECALIso.end(); 
	      ++i) { ffECALIso.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = HCALIso.begin(); i != HCALIso.end(); 
	      ++i) { ffHCALIso.Fill(*i, lumiWeight*PUWeight); }
	 for (vector<float>::const_iterator i = trackIso.begin(); i != trackIso.end(); 
	      ++i) { ffTrackIso.Fill(*i, lumiWeight*PUWeight); }
	 ffLeadingPhotonET.Fill(*(ET.end() - 1), lumiWeight*PUWeight);
	 ffTrailingPhotonET.Fill(*(ET.begin()), lumiWeight*PUWeight);
	 ffNPV.Fill(nGoodRecoPV, lumiWeight*PUWeight);
	 ffRho.Fill(rho, lumiWeight*PUWeight);
	 ffPFElectronMultiplicity.Fill(nPFElectrons, lumiWeight*PUWeight);
	 ffMuonMultiplicity.Fill(nMuons, lumiWeight*PUWeight);
	 ffDijetSize.Fill(dijetSize);
	 ffDiEMSize.Fill(diEMSize);
	 ffPTHat.Fill(pTHat, lumiWeight*PUWeight);
	 ffCleanHT.Fill(cleanHT, lumiWeight*PUWeight);
	 ffCleanMHT.Fill(cleanMHT, lumiWeight*PUWeight);
	 ffCleanNj.Fill(cleanNj, lumiWeight*PUWeight);
	 ffCleanLeadingJetET.Fill(cleanLeadingJetET, lumiWeight*PUWeight);
// 	 if (e.size() >= 2) {
// 	   vector<float>::iterator i1 = max_element(eET.begin(), eET.end());
// 	   const unsigned int index1 = i1 - eET.begin();
// 	   eET.erase(i1);
// 	   const unsigned int index2 = max_element(fET.begin(), fET.end()) - fET.begin();
// 	   eeffMETVsDiEMETVsInvMass.Fill(invMass, (f[index1] + f[index2]).Pt(), MET, 
// 					 lumiWeight*PUWeight);
// 	 }
	 ffHTVsMET.Fill(MET, oldHT, lumiWeight*PUWeight);
	 ffMETVsNJets30.Fill(nJets30, MET, lumiWeight*PUWeight);
	 ffMETVsNJets60.Fill(nJets60, MET, lumiWeight*PUWeight);
	 if (dijetSize == 2) {
	   ffMETVsDiEMETVsInvMass[iNjBin]->Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
	   ffMETVsDiEMETVsInvMassTot.Fill(invMass, dijetPT, MET, lumiWeight*PUWeight);
	   ffDiEMETVec.push_back(dijetPT);
	   ffDijetETUniform[iNjBin]->Fill(dijetPT, lumiWeight*PUWeight);
	   ffJ1UncorrP4.push_back(uncorrP4[0]);
	   ffJ2UncorrP4.push_back(uncorrP4[1]);
	   ffJ1JES.push_back(JES[0]);
	   ffJ2JES.push_back(JES[1]);
	   ffJ1JESErr.push_back(JESErr[0]);
	   ffJ2JESErr.push_back(JESErr[1]);
	 }
	 else {
// 	   printDiObjectErrorMessage(dijetSize, "dijet", evtCategory, jentry);
	   ffMETUnmatched.Fill(MET, lumiWeight*PUWeight);

	   //use di-EM ET when dijet ET is unavailable
	   if (diEMSize == 2) {
	     ffMETVsDiEMETVsInvMass[iNjBin]->Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	     ffMETVsDiEMETVsInvMassTot.Fill(invMass, diEMET, MET, lumiWeight*PUWeight);
	     ffDiEMETVec.push_back(diEMET);
	   }
	   else printDiObjectErrorMessage(diEMSize, "diEM", evtCategory, jentry);
	 }
	 if (diEMSize == 2) {
	   ffMETVec.push_back(MET);
	   ffNjVec.push_back(nJets);
	   ffR2VsMR.Fill(MR, R2, lumiWeight*PUWeight);
	   ffMRVec.push_back(MR);
	   ffR2Vec.push_back(R2);
	   ffDiEMETUniform[iNjBin]->Fill(diEMET, lumiWeight*PUWeight);
	 }
	 else printDiObjectErrorMessage(diEMSize, "diEM", evtCategory, jentry);

	 //print the event information to a file
	 ffEvtFile << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
	 ffEvtFile << susyEvent->luminosityBlockNumber << endl;

	 break;
       default:
	 break;
       }
     }
   }

   //close file with event info
   ggEvtFile.close();
   egEvtFile.close();
   eeEvtFile.close();
   ffEvtFile.close();

   /*
     - the innermost vector of TH1F*s has size 1 in this implementation (corresponds to the use of 
     1 inclusive MET bin)
     - the vector<vector<TH1F*> > is filled with 3 vector<TH1F*> for the 3 Nj bins
   */
   vector<vector<TH1F*> > eeDiEMETScaled(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeLowSidebandDiEMETScaled(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeHighSidebandDiEMETScaled(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > ffDiEMETScaled(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeWeights(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeLowSidebandWeights(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeHighSidebandWeights(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > ffWeights(nNjBins, vector<TH1F*>());

   /*
     - the vector<vector<vector<TH1F*>* > > is filled with 4 vector<vector<TH1F*>* > for the 4 
     control samples ee, eeLowSideband, eeHighSideband, ff
    */
   vector<vector<vector<TH1F*> >* > diEMETScaled;
   diEMETScaled.push_back(&eeDiEMETScaled);
   diEMETScaled.push_back(&eeLowSidebandDiEMETScaled);
   diEMETScaled.push_back(&eeHighSidebandDiEMETScaled);
   diEMETScaled.push_back(&ffDiEMETScaled);
   vector<vector<vector<TH1F*> >* > weights;
   weights.push_back(&eeWeights);
   weights.push_back(&eeLowSidebandWeights);
   weights.push_back(&eeHighSidebandWeights);
   weights.push_back(&ffWeights);

   /*
     - the innermost vector of TH3F*s has size 3 for the 3 Nj bins
     - the vector<vector<TH3F*> > is filled with 4 vector<TH3F*> for the 4 control samples ee, 
     eeLowSideband, eeHighSideband, ff
    */
   vector<vector<TH3F*>* > TH3Fs;
   TH3Fs.push_back(&eeMETVsDiEMETVsInvMass);
   TH3Fs.push_back(&eeLowSidebandMETVsDiEMETVsInvMass);
   TH3Fs.push_back(&eeHighSidebandMETVsDiEMETVsInvMass);
   TH3Fs.push_back(&ffMETVsDiEMETVsInvMass);

   vector<TH3F*> TH3FTots;
   TH3FTots.push_back(&eeMETVsDiEMETVsInvMassTot);
   TH3FTots.push_back(&eeLowSidebandMETVsDiEMETVsInvMassTot);
   TH3FTots.push_back(&eeHighSidebandMETVsDiEMETVsInvMassTot);
   TH3FTots.push_back(&ffMETVsDiEMETVsInvMassTot);

   //names of 4 control samples
   vector<string> controlSamples;
   controlSamples.push_back("ee");
   controlSamples.push_back("eeLowSideband");
   controlSamples.push_back("eeHighSideband");
   controlSamples.push_back("ff");

   //fill weights histograms
   for (vector<string>::const_iterator iControlSample = controlSamples.begin(); 
	iControlSample != controlSamples.end(); ++iControlSample) {
     for (unsigned int iNjBin = 0; iNjBin < nNjBins; ++iNjBin) {
       const Double_t* pDiEMETBins = diEMETBins1j;
       unsigned int pNDiEMETBins = nDiEMETBins1j;
       if (iNjBin == 0) {
	 pDiEMETBins = diEMETBins0j;
	 pNDiEMETBins = nDiEMETBins0j;
       }
       TH1F* histDiEMETScaled = 
	 new TH1F(histName(*iControlSample, "DiEMETScaled_METBin1_NjBin", iNjBin + 1).c_str(), "", 
		  pNDiEMETBins, pDiEMETBins);
       TH1F* histWeights = new TH1F(histName(*iControlSample, "Weights_METBin1_NjBin", 
					     iNjBin + 1).c_str(), "", pNDiEMETBins, pDiEMETBins);
       setHistogramOptions(*histDiEMETScaled, "Dijet p_{T} (GeV)", "", "", kSumW2);
       setHistogramOptions(*histWeights, "Dijet p_{T} (GeV)", "", "", kSumW2);
       const unsigned int i = iControlSample - controlSamples.begin();
       TH3F* pGG3D = ggMETVsDiEMETVsInvMass[iNjBin];
       if (*iControlSample == "ff") pGG3D = ggMETVsDijetETVsInvMass[iNjBin];
       fillWeightsHistograms(pGG3D->
			     ProjectionY(histName("ggDiEMET_METBin1_", "NjBin", 
						  iNjBin + 1).c_str(), 0, -1, 0, -1, "e"), 
			     TH3Fs[i]->at(iNjBin)->
			     ProjectionY(histName(*iControlSample, "DiEMET_METBin1_NjBin", 
						  iNjBin + 1).c_str(), 0, -1, 0, -1, "e"), 
			     *histDiEMETScaled, *histWeights, 
			     ggMETVsDiEMETVsInvMassTot.Integral(0, -1, 0, -1, 0, -1)/
			     TH3FTots[i]->Integral(0, -1, 0, -1, 0, -1));
       diEMETScaled[i]->at(iNjBin).push_back(histDiEMETScaled);
       weights[i]->at(iNjBin).push_back(histWeights);
     }
   }

//    this whole section needs to be updated to account for Nj binning
//    //generate toys for calculating error due to JES for ff reweighting
//    vector<TH1F*> ffDijetWeightsToyDists;
//    estimateJESError(nToys, ffDijetWeightsToyDists, ffWeights, nDiEMETBins, diEMETBins, ffJ1JES, 
// 		    ffJ2JES, ffJ1JESErr, ffJ2JESErr, ffJ1UncorrP4, ffJ2UncorrP4, ggJ1JES, ggJ2JES, 
// 		    ggJ1JESErr, ggJ2JESErr, ggJ1UncorrP4, ggJ2UncorrP4);

//    //add JES error in quadrature to statistical error of weights, bin by bin
//    format("Dijet pT bin", "Weight", "Stat. error", "JES error", "Total error");
//    for (unsigned int iDijetPTBin = 1; iDijetPTBin <= nDiEMETBins; ++iDijetPTBin) {
//      stringstream bin;
//      bin << diEMETBins[iDijetPTBin - 1] << "-" << diEMETBins[iDijetPTBin] << " GeV";
//      float statErrSquared = ffWeights[0]->GetBinError(iDijetPTBin);
//      float JESErrSquared = ffDijetWeightsToyDists[iDijetPTBin - 1]->GetRMS();
//      format(bin.str().c_str(), ffWeights[0]->GetBinContent(iDijetPTBin), statErrSquared, 
// 	    JESErrSquared);
//      statErrSquared*=statErrSquared;
//      JESErrSquared*=JESErrSquared;
// //      ffWeights[0]->SetBinError(iDijetPTBin, sqrt(statErrSquared + JESErrSquared));
//    }

   /*
     - innermost vector<TH1F*> or vector<TH2F*> is filled 1000 times for the 1000 toys
     - vector<vector<TH1F*> > or vector<vector<TH2F*> > is filled 3 times for the 3 Nj bins
    */
   vector<vector<TH1F*> > ffWeightsToys(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeWeightsToys(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeLowSidebandWeightsToys(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeHighSidebandWeightsToys(nNjBins, vector<TH1F*>());

   /*
     - innermost vector<TH1F*> is filled once per di-EM pT bin
     - vector<vector<TH1F*> > is filled 3 times for the 3 Nj bins
    */
   vector<vector<TH1F*> > ffDiEMETToyDistsByBin(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeDiEMETToyDistsByBin(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeLowSidebandDiEMETToyDistsByBin(nNjBins, vector<TH1F*>());
   vector<vector<TH1F*> > eeHighSidebandDiEMETToyDistsByBin(nNjBins, vector<TH1F*>());

   //generate toy di-EM pT weights for calculating error due to MET/razor shape from reweighting
   for (unsigned int iNjBin = 0; iNjBin < nNjBins; ++iNjBin) {
     const Double_t* pDiEMETBins = diEMETBins1j;
     unsigned int pNDiEMETBins = nDiEMETBins1j;
     if (iNjBin == 0) {
       pDiEMETBins = diEMETBins0j;
       pNDiEMETBins = nDiEMETBins0j;
     }
     generateToys(ffDiEMETToyDistsByBin[iNjBin], ffWeightsToys[iNjBin], ffWeights[iNjBin], nToys, 
		  histName("ff", "NjBin", iNjBin + 1), pDiEMETBins);
     generateToys(eeDiEMETToyDistsByBin[iNjBin], eeWeightsToys[iNjBin], eeWeights[iNjBin], nToys, 
		  histName("ee", "NjBin", iNjBin + 1), pDiEMETBins);
     generateToys(eeLowSidebandDiEMETToyDistsByBin[iNjBin], eeLowSidebandWeightsToys[iNjBin], 
		  eeLowSidebandWeights[iNjBin], nToys, 
		  histName("eeLowSideband", "NjBin", iNjBin + 1), pDiEMETBins);
     generateToys(eeHighSidebandDiEMETToyDistsByBin[iNjBin], eeHighSidebandWeightsToys[iNjBin], 
		  eeHighSidebandWeights[iNjBin], nToys, 
		  histName("eeHighSideband", "NjBin", iNjBin + 1), pDiEMETBins);
   }

   //vector<TH1F*> or vector<TH2F*> is filled 1000 times for the 1000 toys
   vector<TH1F*> ffFinalToy;
   vector<TH1F*> eeFinalToy;
   vector<TH1F*> eeLowSidebandFinalToy;
   vector<TH1F*> eeHighSidebandFinalToy;
   vector<TH2F*> ffFinalToyRazor;
   vector<TH2F*> eeFinalToyRazor;
   vector<TH2F*> eeLowSidebandFinalToyRazor;
   vector<TH2F*> eeHighSidebandFinalToyRazor;

   //use the toy di-EM pT weights to generate toy MET distributions
   makeToyMETDists(nToys, "ff", nMETBins, METBins, nMRBins, MRBins, nR2Bins, R2Bins, ffMETVec, 
		   ffMRVec, ffR2Vec, ffDiEMETVec, ffNjVec, nNjBins, NjBins, ffWeightsToys, 
		   ffFinalToy, ffFinalToyRazor, ffLumiWeightVec, ffPUWeightVec);
   makeToyMETDists(nToys, "ee", nMETBins, METBins, nMRBins, MRBins, nR2Bins, R2Bins, eeMETVec, 
		   eeMRVec, eeR2Vec, eeDiEMETVec, eeNjVec, nNjBins, NjBins, eeWeightsToys, 
		   eeFinalToy, eeFinalToyRazor, eeLumiWeightVec, eePUWeightVec);
   makeToyMETDists(nToys, "eeLowSideband", nMETBins, METBins, nMRBins, MRBins, nR2Bins, R2Bins, 
		   eeLowSidebandMETVec, eeLowSidebandMRVec, eeLowSidebandR2Vec, 
		   eeLowSidebandDiEMETVec, eeLowSidebandNjVec, nNjBins, NjBins, 
		   eeLowSidebandWeightsToys, eeLowSidebandFinalToy, eeLowSidebandFinalToyRazor, 
		   eeLowSidebandLumiWeightVec, eeLowSidebandPUWeightVec);
   makeToyMETDists(nToys, "eeHighSideband", nMETBins, METBins, nMRBins, MRBins, nR2Bins, R2Bins, 
		   eeHighSidebandMETVec, eeHighSidebandMRVec, eeHighSidebandR2Vec, 
		   eeHighSidebandDiEMETVec, eeHighSidebandNjVec, nNjBins, NjBins, 
		   eeHighSidebandWeightsToys, eeHighSidebandFinalToy, eeHighSidebandFinalToyRazor, 
		   eeHighSidebandLumiWeightVec, eeHighSidebandPUWeightVec);

   //reweight ee and ff MET/razor histograms
   //changed reweightDefault to no di-EM pT reweighting
   reweightDefault(eeLowSidebandMETVec, eeLowSidebandMRVec, eeLowSidebandR2Vec, 
		   eeLowSidebandDiEMETVec, eeLowSidebandNjVec, nNjBins, NjBins, 
		   eeLowSidebandWeights, &eeLowSidebandFinal, &eeLowSidebandR2VsMRFinal, 
		   eeLowSidebandLumiWeightVec, eeLowSidebandPUWeightVec);
   reweightDefault(eeMETVec, eeMRVec, eeR2Vec, eeDiEMETVec, eeNjVec, nNjBins, NjBins, eeWeights, 
		   &eeFinal, &eeR2VsMRFinal, eeLumiWeightVec, eePUWeightVec);
   reweightDefault(eeHighSidebandMETVec, eeHighSidebandMRVec, eeHighSidebandR2Vec, 
		   eeHighSidebandDiEMETVec, eeHighSidebandNjVec, nNjBins, NjBins, 
		   eeHighSidebandWeights, &eeHighSidebandFinal, &eeHighSidebandR2VsMRFinal, 
		   eeHighSidebandLumiWeightVec, eeHighSidebandPUWeightVec);
   reweightDefault(ffMETVec, ffMRVec, ffR2Vec, ffDiEMETVec, ffNjVec, nNjBins, NjBins, ffWeights, 
		   &ffFinal, &ffR2VsMRFinal, ffLumiWeightVec, ffPUWeightVec);

   //make individual histograms of MET/razor toy distributions, 1 per MET/razor bin
   vector<TH1F*> ffMETToyDistsByBin;
   vector<TH1F*> eeMETToyDistsByBin;
   vector<TH1F*> eeLowSidebandMETToyDistsByBin;
   vector<TH1F*> eeHighSidebandMETToyDistsByBin;
   vector<TH1F*> ffRazorToyDistsByBin;
   vector<TH1F*> eeRazorToyDistsByBin;
   vector<TH1F*> eeLowSidebandRazorToyDistsByBin;
   vector<TH1F*> eeHighSidebandRazorToyDistsByBin;
   fillToyDistributions(ffMETToyDistsByBin, ffRazorToyDistsByBin, ffFinal, ffR2VsMRFinal, 
			ffToyCanvas, ffFinalToy, ffFinalToyRazor, nToys, "ff", nMETBins, nMRBins, 
			nR2Bins);
   fillToyDistributions(eeMETToyDistsByBin, eeRazorToyDistsByBin, eeFinal, eeR2VsMRFinal, 
			eeToyCanvas, eeFinalToy, eeFinalToyRazor, nToys, "ee", nMETBins, nMRBins, 
			nR2Bins);
   fillToyDistributions(eeLowSidebandMETToyDistsByBin, eeLowSidebandRazorToyDistsByBin, 
			eeLowSidebandFinal, eeLowSidebandR2VsMRFinal, eeLowSidebandToyCanvas, 
			eeLowSidebandFinalToy, eeLowSidebandFinalToyRazor, nToys, "eeLowSideband", 
			nMETBins, nMRBins, nR2Bins);
   fillToyDistributions(eeHighSidebandMETToyDistsByBin, eeHighSidebandRazorToyDistsByBin, 
			eeHighSidebandFinal, eeHighSidebandR2VsMRFinal, eeHighSidebandToyCanvas, 
			eeHighSidebandFinalToy, eeHighSidebandFinalToyRazor, nToys, 
			"eeHighSideband", nMETBins, nMRBins, nR2Bins);

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
   egR2VsMR.Scale(egScale);
   for (Int_t iMRBin = 1; iMRBin <= egR2VsMR.GetNbinsX(); ++iMRBin) {
     for (Int_t iR2Bin = 1; iR2Bin <= egR2VsMR.GetNbinsY(); ++iR2Bin) {
       const float binContent = egR2VsMR.GetBinContent(iMRBin, iR2Bin);
       const float binPoissonError = egR2VsMR.GetBinError(iMRBin, iR2Bin);
       if (binContent == 0.0) egR2VsMR.SetBinError(iMRBin, 0.0);
       else {
	 egR2VsMR.SetBinError(iMRBin, 
			      egScale*sqrt(((binPoissonError*binPoissonError)/
					    (binContent*binContent)) + 
					   ((egMisIDRateErr*egMisIDRateErr)/
					    (egMisIDRate*egMisIDRate))));
       }
     }
   }

   //print weights error information
   for (unsigned int iMETBin = 1; iMETBin <= nMETBins; ++iMETBin) {
     cout << "MET bin " << iMETBin << endl;
     cout << "eeLowSideband\n";
     cout << "Mean of toys: " << eeLowSidebandMETToyDistsByBin[iMETBin - 1]->GetMean();
     cout << ", observed no. events: " << eeLowSidebandFinal.GetBinContent(iMETBin);
     cout << ", RMS of toys: " << eeLowSidebandMETToyDistsByBin[iMETBin - 1]->GetRMS() << endl;
     // //make canvases with TLines at the observed point if need be
     // TLine obs(eeLowSidebandFinal.GetBinContent(iMETBin), gPad->getUymin(), 
     // 	       eeLowSidebandFinal.GetBinContent(iMETBin), gPad->GetUymax());
     cout << "eeHighSideband\n";
     cout << "Mean of toys: " << eeHighSidebandMETToyDistsByBin[iMETBin - 1]->GetMean();
     cout << ", observed no. events: " << eeHighSidebandFinal.GetBinContent(iMETBin);
     cout << ", RMS of toys: " << eeHighSidebandMETToyDistsByBin[iMETBin - 1]->GetRMS() << endl;
     cout << "ee\n";
     cout << "Mean of toys: " << eeMETToyDistsByBin[iMETBin - 1]->GetMean();
     cout << ", observed no. events: " << eeFinal.GetBinContent(iMETBin) << ", RMS of toys: ";
     cout << eeMETToyDistsByBin[iMETBin - 1]->GetRMS() << endl;
     cout << "ff\n";
     cout << "Mean of toys: " << ffMETToyDistsByBin[iMETBin - 1]->GetMean();
     cout << ", observed no. events: " << ffFinal.GetBinContent(iMETBin) << ", RMS of toys: ";
     cout << ffMETToyDistsByBin[iMETBin - 1]->GetRMS() << endl;
     cout << "--------------------";
   }

   //normalize ee and ff MET histograms
   eeFinal.Add(&eeLowSidebandFinal, -1.0);
   eeFinal.Add(&eeHighSidebandFinal, -1.0);
   float eeNormErrSquared = 0.0;
   float ffNormErrSquared = 0.0;
   //norm calculated assuming no EW contribution
   const float eeNorm = 
     normAndErrorSquared(ggMETVsDiEMETVsInvMassTot, eeFinal, egMET, maxNormBin, eeNormErrSquared);
//    const float eeNorm = 
//      normAndErrorSquared(*ggMETVsDiEMETVsInvMass[1], eeFinal, egMET, maxNormBin, eeNormErrSquared);
   const float ffNorm = 
     normAndErrorSquared(ggMETVsDiEMETVsInvMassTot, ffFinal, egMET, maxNormBin, ffNormErrSquared);
   eeFinal.Scale(eeNorm);
   ffFinal.Scale(ffNorm);

   //normalize ee and ff razor histograms
   eeLowSidebandR2VsMRFinal.Scale(2.0);
   eeHighSidebandR2VsMRFinal.Scale(2.0);
   eeR2VsMRFinal.Add(&eeLowSidebandR2VsMRFinal, -1.0);
   eeR2VsMRFinal.Add(&eeHighSidebandR2VsMRFinal, -1.0);
   float eeR2VsMRNormErrSquared = 0.0;
   float ffR2VsMRNormErrSquared = 0.0;
   const float eeR2VsMRNorm = 
     razorNormAndErrorSquared(ggR2VsMR, eeR2VsMRFinal, egR2VsMR, maxMRNormBin, maxR2NormBin, 
			      eeR2VsMRNormErrSquared);
   const float ffR2VsMRNorm = 
     razorNormAndErrorSquared(ggR2VsMR, ffR2VsMRFinal, egR2VsMR, maxMRNormBin, maxR2NormBin, 
			      ffR2VsMRNormErrSquared);
   eeR2VsMRFinal.Scale(eeR2VsMRNorm);
   ffR2VsMRFinal.Scale(ffR2VsMRNorm);

   //set MET error bars
   setMETErrorBars(ffFinal, ffMETToyDistsByBin, vector<TH1F*>(), vector<TH1F*>(), ffNorm, 
   		   ffNormErrSquared);
   setMETErrorBars(eeFinal, eeMETToyDistsByBin, eeLowSidebandMETToyDistsByBin, 
   		   eeHighSidebandMETToyDistsByBin, eeNorm, eeNormErrSquared);

   //set razor error bars
   setRazorErrorBars(ffR2VsMRFinal, ffRazorToyDistsByBin, vector<TH1F*>(), vector<TH1F*>(), 
		     ffR2VsMRNorm, ffR2VsMRNormErrSquared);
   setRazorErrorBars(eeR2VsMRFinal, eeRazorToyDistsByBin, eeLowSidebandRazorToyDistsByBin, 
		     eeHighSidebandRazorToyDistsByBin, eeR2VsMRNorm, eeR2VsMRNormErrSquared);

   //add EW contribution to QCD control samples
//    eeFinal.Add(egMET);
//    ffFinal.Add(egMET);
   eeR2VsMRFinal.Add(&egR2VsMR);
   ffR2VsMRFinal.Add(&egR2VsMR);

   //make final MET canvas
   METCanvas.cd();
   makeFinalCanvas(dynamic_cast<TH1*>(&eeFinal), 4, 2, 3005, 4, 0, 1, "E2");
//    makeFinalCanvas(dynamic_cast<TH1*>(&ffFinal), kMagenta, 2, 3004, kMagenta, 0, 1, "E2SAME");
//    makeFinalCanvas(dynamic_cast<TH1*>(egMET), 8, 2, 3003, 8, 1, 1, "HISTSAME");
   makeFinalCanvas(dynamic_cast<TH1*>(ggMETVsDiEMETVsInvMassTot.
				      ProjectionZ("ggMET", 0, -1, 0, -1, "e")), 1, 1, 0, 0, 1, 1, 
		   "SAME");
//    makeFinalCanvas(dynamic_cast<TH1*>(ggMETVsDiEMETVsInvMass[1]->
// 				      ProjectionZ("ggMET", 0, -1, 0, -1, "e")), 1, 1, 0, 0, 1, 1, 
// 		   "SAME");

   //make final razor canvases
   MRCanvas.cd();
   makeFinalCanvas(dynamic_cast<TH1*>(eeR2VsMRFinal.ProjectionX("eeMR", 0, -1, "")), 4, 2, 3005, 
		   4, 0, 1, "E2");
   makeFinalCanvas(dynamic_cast<TH1*>(ffR2VsMRFinal.ProjectionX("ffMR", 0, -1, "")), kMagenta, 2, 
		   3004, kMagenta, 0, 1, "E2SAME");
   makeFinalCanvas(dynamic_cast<TH1*>(egR2VsMR.ProjectionX("egMR", 0, -1, "")), 8, 2, 3003, 8, 1, 
		   1, "HISTSAME");
   makeFinalCanvas(dynamic_cast<TH1*>(ggR2VsMR.ProjectionX("ggMR", 0, -1, "")), 1, 1, 0, 0, 1, 1, 
		   "SAME");
   R2Canvas.cd();
   makeFinalCanvas(dynamic_cast<TH1*>(eeR2VsMRFinal.ProjectionY("eeR2", 0, -1, "")), 4, 2, 3005, 
		   4, 0, 1, "E2");
   makeFinalCanvas(dynamic_cast<TH1*>(ffR2VsMRFinal.ProjectionY("ffR2", 0, -1, "")), kMagenta, 2, 
		   3004, kMagenta, 0, 1, "E2SAME");
   makeFinalCanvas(dynamic_cast<TH1*>(egR2VsMR.ProjectionY("egR2", 0, -1, "")), 8, 2, 3003, 8, 1, 
		   1, "HISTSAME");
   makeFinalCanvas(dynamic_cast<TH1*>(ggR2VsMR.ProjectionY("ggR2", 0, -1, "")), 1, 1, 0, 0, 1, 1, 
		   "SAME");

   //save
   out.cd();
   eeToyCanvas.Write();
   eeLowSidebandToyCanvas.Write();
   eeHighSidebandToyCanvas.Write();
   ffToyCanvas.Write();
   METCanvas.Write();
   MRCanvas.Write();
   R2Canvas.Write();
   out.Write();

   //deallocate memory
   deallocateMemory(ggMETVsDiEMETVsInvMass);
   deallocateMemory(ggMETVsDijetETVsInvMass);
   deallocateMemory(eeMETVsDiEMETVsInvMass);
   deallocateMemory(eeLowSidebandMETVsDiEMETVsInvMass);
   deallocateMemory(eeHighSidebandMETVsDiEMETVsInvMass);
   deallocateMemory(ffMETVsDiEMETVsInvMass);
   deallocateMemory(ggDiEMETUniform);
   deallocateMemory(eeDiEMETUniform);
   deallocateMemory(ffDiEMETUniform);
   deallocateMemory(eeDiEMETScaled);
   deallocateMemory(eeLowSidebandDiEMETScaled);
   deallocateMemory(eeHighSidebandDiEMETScaled);
   deallocateMemory(ffDiEMETScaled);
   deallocateMemory(eeWeights);
   deallocateMemory(eeLowSidebandWeights);
   deallocateMemory(eeHighSidebandWeights);
   deallocateMemory(ffWeights);
   deallocateMemory(ffWeightsToys);
   deallocateMemory(eeWeightsToys);
   deallocateMemory(eeLowSidebandWeightsToys);
   deallocateMemory(eeHighSidebandWeightsToys);
   deallocateMemory(ffFinalToy);
   deallocateMemory(eeFinalToy);
   deallocateMemory(eeLowSidebandFinalToy);
   deallocateMemory(eeHighSidebandFinalToy);
   deallocateMemory(ffFinalToyRazor);
   deallocateMemory(eeFinalToyRazor);
   deallocateMemory(eeLowSidebandFinalToyRazor);
   deallocateMemory(eeHighSidebandFinalToyRazor);
   deallocateMemory(ffMETToyDistsByBin);
   deallocateMemory(eeMETToyDistsByBin);
   deallocateMemory(eeLowSidebandMETToyDistsByBin);
   deallocateMemory(eeHighSidebandMETToyDistsByBin);
   deallocateMemory(ffRazorToyDistsByBin);
   deallocateMemory(eeRazorToyDistsByBin);
   deallocateMemory(eeLowSidebandRazorToyDistsByBin);
   deallocateMemory(eeHighSidebandRazorToyDistsByBin);
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

   //lumi and PU weight vectors
   VFLOAT lumiWeightVec;
   VFLOAT PUWeightVec;

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
     if (evtCategory != FAIL) {
       rhoPreselected.Fill(rho, lumiWeight*PUWeight);
       lumiWeightVec.push_back(lumiWeight);
       PUWeightVec.push_back(PUWeight);
     }
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
     const float oldHT = accumulate(jetET.begin(), jetET.end(), 0);

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
	 eeHTVsMET.Fill(MET, oldHT, lumiWeight*PUWeight);
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
	 ffHTVsMET.Fill(MET, oldHT, lumiWeight*PUWeight);
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
			   *eeHistDiEMETScaled, *eeHistWeights, 1.0);
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
			     *ffHistDiEMETScaled, *ffHistWeights, 1.0);
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
		METBins, ffMETVec, ffDiEMETVec, lumiWeightVec, PUWeightVec);
   generateToys(eeFinalToy, eeDiEMETToyDistsByBin, eeWeights, nToys, "ee", diEMETBins, nMETBins, 
		METBins, eeBkgSubtractedMETVsDiEMET/*, lumiWeightVec, PUWeightVec*/);

   //reweight ee and ff MET histograms
   for (unsigned int iMETBin = 1; iMETBin <= nMETBins; ++iMETBin) {
     reweightBinned(eeBkgSubtractedMETVsDiEMET, *eeWeights[iMETBin - 1], &eeFinal, iMETBin);
   }
   reweightDefault(ffMETVec, ffDiEMETVec, *ffWeights[0], &ffFinal, lumiWeightVec, PUWeightVec);

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
   makeFinalCanvas(dynamic_cast<TH1*>(&eeFinal), 4, 2, 3005, 4, 0, 1, "E2");
   makeFinalCanvas(dynamic_cast<TH1*>(&ffFinal), kMagenta, 2, 3004, kMagenta, 0, 1, "E2SAME");
   makeFinalCanvas(dynamic_cast<TH1*>(egMET), 8, 2, 3003, 8, 1, 1, "HISTSAME");
   makeFinalCanvas(dynamic_cast<TH1*>(ggMETVsDiEMETVsInvMass.ProjectionZ("ggMET", 0, -1, 0, -1, 
									 "e")), 1, 1, 0, 0, 1, 1, 
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
			  *eeHistDiEMETScaled, *eeHistWeights, 1.0);
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
			    *ffHistDiEMETScaled, *ffHistWeights, 1.0);
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
  // reweightDefault(ffMETVec, ffDiEMETVec, *ffWeights[0], &ffFinal, lumiWeightVec, PUWeightVec);


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
  makeFinalCanvas(dynamic_cast<TH1*>(&eeFinal), 4, 2, 3005, 4, 0, 1, "E2");
  // makeFinalCanvas(dynamic_cast<TH1*>(&ffFinal), kMagenta, 2, 3004, kMagenta, 0, 1, "E2SAME");
  makeFinalCanvas(dynamic_cast<TH1*>(egMET), 8, 2, 3003, 8, 1, 1, "HISTSAME");
  makeFinalCanvas(dynamic_cast<TH1*>(ggMETVsDiEMETVsInvMass->ProjectionZ("ggMET", 0, -1, 0, -1, 
									 "e")), 1, 1, 0, 0, 1, 1, 
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

bool GMSBAnalyzer::passJetFiducialCuts(const TLorentzVector& corrP4, const Double_t pTMin, 
				       const Double_t etaMax) const
{
  return ((corrP4.Pt() > pTMin) && (fabs(corrP4.Eta()) < etaMax));
}

bool GMSBAnalyzer::passPFJetID(susy::PFJetCollection::const_iterator& iJet) const
{
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
  return passed;
}

bool GMSBAnalyzer::passCaloJetID(susy::CaloJetCollection::const_iterator& iJet) const
{
  bool passed = false;
  float absEta = fabs(iJet->detectorP4.Eta());
  if (absEta < 2.6) { //HB/HE
    if ((iJet->emEnergyFraction > 0.01) && (iJet->n90Hits > 1) && 
	(iJet->fHPD < 0.98)) passed = true;
  }
  //implement cuts for other eta regions when variables needed are stored in the ntuple
  //cf. http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2010_067_v2.pdf
  return passed;
}

unsigned int GMSBAnalyzer::numPhotonOverlaps(const susy::PhotonCollection& photons, 
					     const TLorentzVector& corrP4, 
					     int& overlappingPhotonIndex, const float dR) const
{
  unsigned int numOverlaps = 0;
  for (susy::PhotonCollection::const_iterator iPhoton = photons.begin(); 
       iPhoton != photons.end(); ++iPhoton) {
    if (deltaR(corrP4.Eta(), corrP4.Phi(), iPhoton->caloPosition.Eta(), 
	       iPhoton->caloPosition.Phi()) < dR) {
      overlappingPhotonIndex = iPhoton - photons.begin();
      ++numOverlaps;
    }
  }
  return numOverlaps;
}

vector<unsigned int> GMSBAnalyzer::getDecidingPhotonIndices(const unsigned int numPhotons) const
{
  vector<unsigned int> decidingPhotonIndices;
  for (unsigned int iPhoton = 0; iPhoton < numPhotons; ++iPhoton) {
    if (susyCategory->getIsDeciding(tag_, iPhoton)) decidingPhotonIndices.push_back(iPhoton);
  }
  return decidingPhotonIndices;
}

void GMSBAnalyzer::runEMFractionAnalysis(const string& outputFile)
{
  if (fChain == 0) return;

  //open file
  TFile out(outputFile.c_str(), "RECREATE");
  out.cd();

  //define constants
  const bool kSumW2 = true;
  const bool kSetGrid = true;

  //EM fraction of PF jets
  //(a) overlapping with e, g, f, and fail objects
  //(b) not overlapping with e, g, f, and fail objects
  TH1F eLeadingEMFractionPF("eLeadingEMFractionPF", "eLeadingEMFractionPF", 50, 0.0, 1.0);
  TH1F gLeadingEMFractionPF("gLeadingEMFractionPF", "gLeadingEMFractionPF", 50, 0.0, 1.0);
  TH1F fLeadingEMFractionPF("fLeadingEMFractionPF", "fLeadingEMFractionPF", 50, 0.0, 1.0);
  TH1F eTrailingEMFractionPF("eTrailingEMFractionPF", "eTrailingEMFractionPF", 50, 0.0, 1.0);
  TH1F gTrailingEMFractionPF("gTrailingEMFractionPF", "gTrailingEMFractionPF", 50, 0.0, 1.0);
  TH1F fTrailingEMFractionPF("fTrailingEMFractionPF", "fTrailingEMFractionPF", 50, 0.0, 1.0);
  TH1F failEMFractionPF("failEMFractionPF", "failEMFractionPF", 50, 0.0, 1.0);
  TH1F otherEMFractionPF("otherEMFractionPF", "otherEMFractionPF", 50, 0.0, 1.0);
  setHistogramOptions(eLeadingEMFractionPF, "EM fraction", "", "", kSumW2);
  setHistogramOptions(gLeadingEMFractionPF, "EM fraction", "", "", kSumW2);
  setHistogramOptions(fLeadingEMFractionPF, "EM fraction", "", "", kSumW2);
  setHistogramOptions(eTrailingEMFractionPF, "EM fraction", "", "", kSumW2);
  setHistogramOptions(gTrailingEMFractionPF, "EM fraction", "", "", kSumW2);
  setHistogramOptions(fTrailingEMFractionPF, "EM fraction", "", "", kSumW2);
  setHistogramOptions(failEMFractionPF, "EM fraction", "", "", kSumW2);
  setHistogramOptions(otherEMFractionPF, "EM fraction", "", "", kSumW2);
  map<string, TH1F*> histograms;
  histograms[string(eLeadingEMFractionPF.GetName())] = &eLeadingEMFractionPF;
  histograms[string(gLeadingEMFractionPF.GetName())] = &gLeadingEMFractionPF;
  histograms[string(fLeadingEMFractionPF.GetName())] = &fLeadingEMFractionPF;
  histograms[string(eTrailingEMFractionPF.GetName())] = &eTrailingEMFractionPF;
  histograms[string(gTrailingEMFractionPF.GetName())] = &gTrailingEMFractionPF;
  histograms[string(fTrailingEMFractionPF.GetName())] = &fTrailingEMFractionPF;
  histograms[string(failEMFractionPF.GetName())] = &failEMFractionPF;
  histograms[string(otherEMFractionPF.GetName())] = &otherEMFractionPF;

  /*(ETj - ETg)/ETj vs. EM fraction of PF jets overlapping with e, g, f, and fail objects*/
  TH2F eLeadingETJetNormVsEMFractionPF("eLeadingETJetNormVsEMFractionPF", 
				       "eLeadingETJetNormVsEMFractionPF", 
				       50, 0.0, 1.0, 100, -1.0, 1.0);
  TH2F gLeadingETJetNormVsEMFractionPF("gLeadingETJetNormVsEMFractionPF", 
				       "gLeadingETJetNormVsEMFractionPF", 
				       50, 0.0, 1.0, 100, -1.0, 1.0);
  TH2F fLeadingETJetNormVsEMFractionPF("fLeadingETJetNormVsEMFractionPF", 
				       "fLeadingETJetNormVsEMFractionPF", 
				       50, 0.0, 1.0, 100, -1.0, 1.0);
  TH2F eTrailingETJetNormVsEMFractionPF("eTrailingETJetNormVsEMFractionPF", 
					"eTrailingETJetNormVsEMFractionPF", 
					50, 0.0, 1.0, 100, -1.0, 1.0);
  TH2F gTrailingETJetNormVsEMFractionPF("gTrailingETJetNormVsEMFractionPF", 
					"gTrailingETJetNormVsEMFractionPF", 
					50, 0.0, 1.0, 100, -1.0, 1.0);
  TH2F fTrailingETJetNormVsEMFractionPF("fTrailingETJetNormVsEMFractionPF", 
					"fTrailingETJetNormVsEMFractionPF", 
					50, 0.0, 1.0, 100, -1.0, 1.0);
  TH2F failETJetNormVsEMFractionPF("failETJetNormVsEMFractionPF", "failETJetNormVsEMFractionPF", 
				   50, 0.0, 1.0, 100, -1.0, 1.0);
  setHistogramOptions(eLeadingETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tj} - p_{Te}}{p_{Tj}}", "", kSumW2);
  setHistogramOptions(gLeadingETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tj} - p_{T#gamma}}{p_{Tj}}", "", kSumW2);
  setHistogramOptions(fLeadingETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tj} - p_{Tf}}{p_{Tj}}", "", kSumW2);
  setHistogramOptions(eTrailingETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tj} - p_{Te}}{p_{Tj}}", "", kSumW2);
  setHistogramOptions(gTrailingETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tj} - p_{T#gamma}}{p_{Tj}}", "", kSumW2);
  setHistogramOptions(fTrailingETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tj} - p_{Tf}}{p_{Tj}}", "", kSumW2);
  setHistogramOptions(failETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tj} - p_{Tfail}}{p_{Tj}}", "", kSumW2);
  map<string, TH2F*> histograms2D;
  histograms2D[string(eLeadingETJetNormVsEMFractionPF.GetName())] = 
    &eLeadingETJetNormVsEMFractionPF;
  histograms2D[string(gLeadingETJetNormVsEMFractionPF.GetName())] = 
    &gLeadingETJetNormVsEMFractionPF;
  histograms2D[string(fLeadingETJetNormVsEMFractionPF.GetName())] = 
    &fLeadingETJetNormVsEMFractionPF;
  histograms2D[string(eTrailingETJetNormVsEMFractionPF.GetName())] = 
    &eTrailingETJetNormVsEMFractionPF;
  histograms2D[string(gTrailingETJetNormVsEMFractionPF.GetName())] = 
    &gTrailingETJetNormVsEMFractionPF;
  histograms2D[string(fTrailingETJetNormVsEMFractionPF.GetName())] = 
    &fTrailingETJetNormVsEMFractionPF;
  histograms2D[string(failETJetNormVsEMFractionPF.GetName())] = &failETJetNormVsEMFractionPF;

  /*(ETjj - ETgg)/ETjj vs. average EM fraction of PF jets overlapping with e, g, and f objects*/
  TH2F eeETJetNormVsEMFractionPF("eeETJetNormVsEMFractionPF", "eeETJetNormVsEMFractionPF", 
				 50, 0.0, 1.0, 100, -1.0, 1.0);
  TH2F ggETJetNormVsEMFractionPF("ggETJetNormVsEMFractionPF", "ggETJetNormVsEMFractionPF", 
				 50, 0.0, 1.0, 100, -1.0, 1.0);
  TH2F ffETJetNormVsEMFractionPF("ffETJetNormVsEMFractionPF", "ffETJetNormVsEMFractionPF", 
				 50, 0.0, 1.0, 100, -1.0, 1.0);
  setHistogramOptions(eeETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tjj} - p_{Tee}}{p_{Tjj}}", "", kSumW2);
  setHistogramOptions(ggETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tjj} - p_{T#gamma#gamma}}{p_{Tjj}}", "", kSumW2);
  setHistogramOptions(ffETJetNormVsEMFractionPF, "EM fraction", 
		      "#frac{p_{Tjj} - p_{Tff}}{p_{Tjj}}", "", kSumW2);

  //canvases to compare EM fractions for the different types of EM objects
  TCanvas leadingEMFractionPF("leadingEMFractionPF", "", 600, 600);
  TCanvas trailingEMFractionPF("trailingEMFractionPF", "", 600, 600);
  TCanvas EMFractionPF("EMFractionPF", "", 600, 600);
  setCanvasOptions(leadingEMFractionPF, "EM fraction", "", 0.0, 1.0, 0.0, 10000.0, kSetGrid);
  setCanvasOptions(trailingEMFractionPF, "EM fraction", "", 0.0, 1.0, 0.0, 10000.0, kSetGrid);
  setCanvasOptions(EMFractionPF, "EM fraction", "", 0.0, 1.0, 0.0, 10000.0, kSetGrid);

  //count events where !=2 matches are found
  unsigned int numLessThan2 = 0;
  unsigned int numGreaterThan2 = 0;

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

    //check that needed maps exist in this event
    bool proceed = true;
    if ((susyCategory->getPhotonType()->size() <= 0) || 
	(susyCategory->getEventCategory()->size() <= 0) || 
	(susyCategory->getEvtInvMass()->size() <= 0)) {
      cerr << "Error: susyCategory->getPhotonType()->size() = 0 or ";
      cerr << "susyCategory->getEventCategory()->size() <= 0 or ";
      cerr << "susyCategory->getEvtInvMass()->size() <= 0 for event " << (jentry + 1);
      cerr << ".  Skipping this event.\n";
      proceed = false;
    }
    if (proceed) {

      //corrected p4s of the jets matched to the 2 EM objects if ee, gg, or ff event
      int category = susyCategory->getEventCategory(tag_);
      TLorentzVector j1;
      TLorentzVector j2;
      TLorentzVector EM1;
      TLorentzVector EM2;
      float EMF1 = 0.0;
      float EMF2 = 0.0;
      unsigned int EMObjectCount = 0;

      //get photons
      map<TString, susy::PhotonCollection>::const_iterator iPhotonMap = 
	susyEvent->photons.find((const char*)tag_);
      if (iPhotonMap != susyEvent->photons.end()) {

	//loop over PF jets
	map<TString, susy::PFJetCollection>::const_iterator iPFJets = 
	  susyEvent->pfJets.find("ak5");
	if (iPFJets != susyEvent->pfJets.end()) {
	  for (susy::PFJetCollection::const_iterator iJet = iPFJets->second.begin(); 
	       iJet != iPFJets->second.end(); ++iJet) {

	    //compute corrected P4
	    float theJES = 1.0;
	    TLorentzVector corrP4 = correctedJet4Momentum(iJet, "L1FastL2L3", theJES, "L2L3");

	    //get the indices of the 2 deciding photons
	    vector<unsigned int> decidingPhotonIndices = 
	      getDecidingPhotonIndices(iPhotonMap->second.size());

	    //calculate the jet-EM overlap
	    int overlappingPhotonIndex03 = -1;
	    unsigned int numOverlaps03 = numPhotonOverlaps(iPhotonMap->second, 
							   iJet->momentum, 
							   overlappingPhotonIndex03, 0.3);
	    
	    //does this jet overlap with a deciding photon?
	    vector<unsigned int>::const_iterator iOverlappingPhoton03 = 
	      find(decidingPhotonIndices.begin(), decidingPhotonIndices.end(), 
		   overlappingPhotonIndex03);

	    //if it does, set photon2Index to the index of the other overlapping photon
	    //otherwise, set photon2Index to -1
	    unsigned int photon2Index03 = -1;
	    if (iOverlappingPhoton03 != decidingPhotonIndices.end()) {
	      switch (iOverlappingPhoton03 - decidingPhotonIndices.begin()) {
	      case 0:
		photon2Index03 = *(iOverlappingPhoton03 + 1);
		break;
	      case 1:
		photon2Index03 = *(iOverlappingPhoton03 - 1);
		break;
	      default:
		cerr << "Error: iOverlappingPhoton03 - decidingPhotonIndices.begin() = ";
		cerr << (iOverlappingPhoton03 - decidingPhotonIndices.begin());
		cerr << ".  Assuming no second deciding photon.\n";
		break;
	      }
	    }
	    
	    //fill single object histograms
	    Float_t EMFraction = 
	      (iJet->neutralEmEnergy + iJet->chargedEmEnergy)/iJet->momentum.Energy();
	    Double_t jetPT = corrP4.Pt();
	    TH1F* pHist1D = EMFractionHistogram(overlappingPhotonIndex03, photon2Index03, 
						iPhotonMap->second, numOverlaps03, histograms, 
						"", "PF");
	    TH2F* pHist2D = EMFractionHistogram(overlappingPhotonIndex03, photon2Index03, 
						iPhotonMap->second, numOverlaps03, 
						histograms2D, "ETJetNormVs", "PF");
	    if (pHist1D != NULL) {
	      pHist1D->Fill(EMFraction, lumiWeight*PUWeight);

	      //save corrected p4 of jets matched to EM objects
	      if ((((string(pHist1D->GetName()).find("eLead") != string::npos) || 
		    (string(pHist1D->GetName()).find("eTrail") != string::npos)) && 
		   (category == EE)) || 
		  (((string(pHist1D->GetName()).find("fLead") != string::npos) || 
		    (string(pHist1D->GetName()).find("fTrail") != string::npos)) && 
		   (category == FF)) || 
		  (((string(pHist1D->GetName()).find("gLead") != string::npos) || 
		    (string(pHist1D->GetName()).find("gTrail") != string::npos)) && 
		   (category == GG))) {
		switch (EMObjectCount) {
		case 0:
		  j1 = corrP4;
		  EM1 = iPhotonMap->second[overlappingPhotonIndex03].momentum;
		  EMF1 = EMFraction;
		  ++EMObjectCount;
		  break;
		case 1:
		  j2 = corrP4;
		  EM2 = iPhotonMap->second[overlappingPhotonIndex03].momentum;
		  EMF2 = EMFraction;
		  ++EMObjectCount;
		  break;
		default:
		  cerr << "Error: " << (EMObjectCount + 1) << " matched objects in event ";
		  cerr << (jentry + 1) << ".\n";
		  ++EMObjectCount;
		  break;
		}
	      }
	    }
	    if (pHist2D != NULL) {
	      pHist2D->Fill(EMFraction, (jetPT - iPhotonMap->second[overlappingPhotonIndex03].
					 momentum.Pt())/jetPT, lumiWeight*PUWeight);
	    }
	  }
	}
	else cerr << "Error: ak5 PF jet collection not found in event " << (jentry + 1) << ".\n";
      }
      else {
	cerr << "Error: " << tag_ << " photon collection not found in ";
	cerr << "event " << (jentry + 1) << ".\n";
      }

      //fill di-EM pT histograms
      float diJetPT = (j1 + j2).Pt();
      float diEMPT = (EM1 + EM2).Pt();
      TH2F* pHistDiEM = NULL;
      switch (category) {
      case EE:
	pHistDiEM = &eeETJetNormVsEMFractionPF;
	break;
      case FF:
	pHistDiEM = &ffETJetNormVsEMFractionPF;
	break;
      case GG:
	pHistDiEM = &ggETJetNormVsEMFractionPF;
	break;
      default:
	break;
      }
      if (pHistDiEM != NULL) {
	if (EMObjectCount == 2) pHistDiEM->Fill((EMF1 + EMF2)/2.0, (diJetPT - diEMPT)/diJetPT, 
						lumiWeight*PUWeight);
	else if (EMObjectCount < 2) ++numLessThan2;
	else ++numGreaterThan2;
      }
    }
  }

  cout << "No. gg/ee/ff events with <2 matched EM objects: " << numLessThan2 << endl;
  cout << "No. gg/ee/ff events with >2 matched EM objects: " << numGreaterThan2 << endl;

  leadingEMFractionPF.cd();
  makeFinalCanvas(&eLeadingEMFractionPF, kRed, 1, 0, 0, 1, kRed, "E1");
  makeFinalCanvas(&gLeadingEMFractionPF, kBlue, 1, 0, 0, 1, kBlue, "E1SAME");
  makeFinalCanvas(&fLeadingEMFractionPF, kMagenta, 1, 0, 0, 1, kMagenta, "E1SAME");
  leadingEMFractionPF.Write();
  trailingEMFractionPF.cd();
  makeFinalCanvas(&eTrailingEMFractionPF, kRed, 1, 0, 0, 1, kRed, "E1");
  makeFinalCanvas(&gTrailingEMFractionPF, kBlue, 1, 0, 0, 1, kBlue, "E1SAME");
  makeFinalCanvas(&fTrailingEMFractionPF, kMagenta, 1, 0, 0, 1, kMagenta, "E1SAME");
  trailingEMFractionPF.Write();
  EMFractionPF.cd();
  makeFinalCanvas(&eLeadingEMFractionPF, kRed, 1, 0, 0, 1, kRed, "E1");
  makeFinalCanvas(&gLeadingEMFractionPF, kBlue, 1, 0, 0, 1, kBlue, "E1SAME");
  makeFinalCanvas(&fLeadingEMFractionPF, kMagenta, 1, 0, 0, 1, kMagenta, "E1SAME");
  makeFinalCanvas(&eTrailingEMFractionPF, kRed - 7, 1, 0, 0, 1, kRed - 7, "E1SAME");
  makeFinalCanvas(&gTrailingEMFractionPF, kBlue - 7, 1, 0, 0, 1, kBlue - 7, "E1SAME");
  makeFinalCanvas(&fTrailingEMFractionPF, kMagenta - 7, 1, 0, 0, 1, kMagenta - 7, "E1SAME");
  EMFractionPF.Write();

  //close
  out.Write();
  out.Close();
}

void GMSBAnalyzer::compareDataToMC(const string& outputFile)
{
  if (fChain == 0) return;
  const bool kSumW2 = true;

  //open file
  TFile out(outputFile.c_str(), "RECREATE");
  out.cd();

  //1D comparison histograms
  TH1F mee("mee", "mee", 150, 60.0, 360.0);
  TH1F nJets("nJets", "nJets", 8, -0.5, 7.5); //in ee events

  //1D histogram options
  setHistogramOptions(mee, "m_{ee} (GeV)", "", "", kSumW2);
  setHistogramOptions(nJets, "n_{j}", "", "", kSumW2);

  //set user-specified number of entries to process
  Long64_t nentries = fChain->GetEntriesFast();
  if (nEvts_ != -1) nentries = nEvts_;

  //loop over events
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    susyEvent->Init();
    susyCategory->reset();
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (((jentry + 1) % 10000) == 0) {
      cout << "Event " << (jentry + 1) << endl;
      cout << susyEvent->runNumber << " " << susyEvent->eventNumber << " ";
      cout << susyEvent->luminosityBlockNumber << endl;
    }

    //get int. lumi weight for this dataset
    TString fileName(fChain->GetCurrentFile()->GetName());
    TString dataset;
    /*if (fileName.Contains("ntuple_.*-v[0-9]")) */dataset = fileName("ntuple_.*-v[0-9]");
    if (dataset.Length() >= 7) dataset.Remove(0, 7);
    if ((dataset.Length() >= 3) && 
	(dataset.Contains("v2_.*-v[0-9]"))) dataset.Remove(0, 3); //for TT sample
    string datasetString((const char*)dataset);
    map<string, unsigned int>::const_iterator iDataset = fileMap_.find(datasetString);
    float weight = 1.0;
    if (iDataset != fileMap_.end()) weight = weight_[iDataset->second];

    //get PU weight for this dataset
    /*in-time reweighting only following 
      https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities as of 27-Oct-11 
      and Tessa's recommendation*/
    int nPV  = -1;
    susy::PUSummaryInfoCollection::const_iterator iBX = susyEvent->PU.begin();
    bool foundInTimeBX = false;
    while ((iBX != susyEvent->PU.end()) && !foundInTimeBX) {
      if (iBX->BX == 0) { 
    	 nPV = iBX->numInteractions;
    	 foundInTimeBX = true;
      }
      ++iBX;
    }
    weight*=lumiWeights_.ITweight(nPV);

    //check for corrupted data
    map<TString, susy::PhotonCollection>::const_iterator iPhotonMap = 
      susyEvent->photons.find((const char*)tag_);
    map<TString, susy::PFJetCollection>::const_iterator iJets = susyEvent->pfJets.find("ak5");
    if ((susyCategory->getEventCategory()->size() > 0) && 
	(susyCategory->getEvtInvMass()->size() > 0) && (iPhotonMap != susyEvent->photons.end()) && 
	(iJets != susyEvent->pfJets.end()) && (susyCategory->getIsDeciding()->size() > 0)) {

      //consider ee events
      if (susyCategory->getEventCategory(tag_) == EE) {

	//plot invariant mass
	mee.Fill(susyCategory->getEvtInvMass(tag_), weight);

	//loop over PF jets
	unsigned int oldNumJets = 0;
	for (susy::PFJetCollection::const_iterator iJet = iJets->second.begin(); 
	     iJet != iJets->second.end(); ++iJet) {

	  //does jet pass ET, eta, and jet ID cuts?
	  float theJES = 1.0;
	  TLorentzVector corrP4 = correctedJet4Momentum(iJet, "L1FastL2L3", theJES);
	  if ((passJetFiducialCuts(corrP4, 30.0/*GeV*/, 2.6)) && passPFJetID(iJet)) {

	    //does jet not overlap with either e?
	    bool nonoverlapping = true;
	    susy::PhotonCollection::const_iterator iPhoton = iPhotonMap->second.begin();
	    while ((iPhoton != iPhotonMap->second.end()) && nonoverlapping) {
	      if (susyCategory->getIsDeciding(tag_, iPhoton - iPhotonMap->second.begin()) && 
		  deltaR(iJet->momentum.Eta(), iJet->momentum.Phi(), iPhoton->momentum.Eta(), 
			 iPhoton->momentum.Phi()) < 0.8) nonoverlapping = false;
	      ++iPhoton;
	    }
	    if (nonoverlapping) ++oldNumJets;
	  }
	}

	//plot nJets
	nJets.Fill(oldNumJets, weight);
      }
    }
    else {
      cerr << "Error: corrupted data in run " << susyEvent->runNumber << ", event ";
      cerr << susyEvent-> eventNumber << ", lumi section " << susyEvent->luminosityBlockNumber;
      cerr << " (event " << (jentry + 1) << ").\n";
      cerr << "     susyCategory->getEventCategory()->size() = ";
      cerr << susyCategory->getEventCategory()->size() << endl;
      cerr << "     susyCategory->getEvtInvMass()->size() = ";
      cerr << susyCategory->getEvtInvMass()->size() << endl;
      cerr << "     iPhotonMap != susyEvent->photons.end() = ";
      cerr << (iPhotonMap != susyEvent->photons.end() ? true : false) << endl;
      cerr << "     iJets != susyEvent->pfJets.end() = ";
      cerr << (iJets != susyEvent->pfJets.end() ? true : false) << endl;
      cerr << "     susyCategory->getIsDeciding()->size() = ";
      cerr << susyCategory->getIsDeciding()->size() << endl;
      cerr << "Skipping event.\n";
    }
  }

  //close file
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

  //MET histograms
  TH1F uncorrMET("uncorrMET", "uncorrMET", 75, 0.0, 150.0);
  TH1F corrMET("corrMET", "corrMET", 75, 0.0, 150.0);
  setHistogramOptions(uncorrMET, "ME_{T} (GeV)", "", "", true);
  setHistogramOptions(corrMET, "ME_{T} (GeV)", "", "", true);

  //MET canvas
  TCanvas METCanvas("METCanvas", "METCanvas", 600, 600);
  setCanvasOptions(METCanvas, "ME_{T} (GeV)", "", 0.0, 150.0, 0.0, 20000.0, true);

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
//     if (susyCategory->getEventCategory(tag_) == EE) {
//       const double invMass = susyCategory->getEvtInvMass(tag_);
//       mee.Fill(invMass);
//     }

    //print debug information
    for (susy::TriggerMap::const_iterator i = susyEvent->hltMap.begin(); 
	 i != susyEvent->hltMap.end(); ++i) {
      cout << i->first << ": " << (int(i->second.second)) << endl;
    }
    for (susy::TriggerMap::const_iterator i = susyEvent->l1Map.begin(); 
	 i != susyEvent->l1Map.end(); ++i) {
      cout << i->first << ": " << (int(i->second.second)) << endl;
    }
    for (susy::PUSummaryInfoCollection::const_iterator i = susyEvent->PU.begin(); 
	 i != susyEvent->PU.end(); ++i) { cout << i->numInteractions << endl; }

    //plot uncorrected and corrected MET
    map<TString, susy::MET>::const_iterator iUncorrMET = susyEvent->metMap.find("pfMet");
//     map<TString, susy::MET>::const_iterator iCorrMET = 
//       susyEvent->metMap.find("pfType1CorrectedMet");
    if ((iUncorrMET != susyEvent->metMap.end())/* && (iCorrMET != susyEvent->metMap.end())*/) {
      cerr << "Uncorrected MET: " << iUncorrMET->second.met() << " GeV\n";
//       cerr << "Corrected MET: " << iCorrMET->second.met() << " GeV\n";
      uncorrMET.Fill(iUncorrMET->second.met());
//       corrMET.Fill(iCorrMET->second.met());
    }
    else cerr << "Error: could not find MET collections.\n";
  }

  //write MET canvas
  METCanvas.cd();
  makeFinalCanvas(&uncorrMET, kBlack, 1, 0, 0, 0.7, kBlack, "E2");
  makeFinalCanvas(&corrMET, kRed, 1, 0, 0, 0.7, kRed, "E2SAME");
  METCanvas.Write();

  //write histograms
  // outTxt.close();
  mee.Write();
  out.Write();
  out.Close();
}

void GMSBAnalyzer::plotTemp(string& outputFile)
{
  if (fChain == 0) return;

  //open files
  TFile out(outputFile.c_str(), "RECREATE");
  out.cd();

  //R9 vs. SC eta histogram
  TH2F R9VsSCEta("R9VsSCEta", "", 100, -1.4442, 1.4442, 100, 0.0, 1.0);
  setHistogramOptions(R9VsSCEta, "SC #eta", "R9", "", true);

  //set user-specified number of entries to process
  Long64_t nentries = fChain->GetEntriesFast();
  if (nEvts_ != -1) nentries = nEvts_;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    //plot R9 vs. SC eta for gs
    if (susyCategory->getEventCategory(tag_) == GG) {
      map<TString, susy::PhotonCollection>::const_iterator iPhotons = susyEvent->photons.find((const char*)tag_);
      if (iPhotons != susyEvent->photons.end()) {
	for (susy::PhotonCollection::const_iterator iPhoton = iPhotons->second.begin(); iPhoton != iPhotons->second.end(); ++iPhoton) {
	  const unsigned int i = iPhoton - iPhotons->second.begin();
	  if (susyCategory->getIsDeciding(tag_, i)) R9VsSCEta.Fill(iPhoton->caloPosition.Eta(), iPhoton->r9);
	}
      }
      else cerr << "Error: no photon collection with tag \"" << tag_ << "\" found in event " << (jentry + 1) << ".\n";
    }
  }

  //write
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
      cout << "-------ITRK : " << iPhoton->trkSumPtHollowConeDR03 << endl;
      cout << "--------------";
      cout << (susyCategory->getPassTrackIsoMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------H/E  : " << iPhoton->hadronicOverEm << endl;
      cout << "--------------";
      cout << (susyCategory->getPassHOverEMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------R9MAX: " << iPhoton->r9 << endl;
      cout << "--------------";
      cout << (susyCategory->getPassR9Max(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------R9MIN: " << iPhoton->r9 << endl;
      cout << "--------------";
      cout << (susyCategory->getPassR9Min(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------t    : " << iPhoton->seedTime << endl;
      cout << "--------------";
      cout << (susyCategory->getPassAbsSeedTimeMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------E2/E9: " << (double)((iPhoton->e1x2)/(iPhoton->e3x3)) << endl;
      cout << "--------------";
      cout << (susyCategory->getPassE2OverE9Max(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "----------------------------";
      cout << (susyCategory->getPassPreselection(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------ICOMB: " << (iPhoton->ecalRecHitSumEtConeDR03 + 
				   iPhoton->trkSumPtHollowConeDR03 + 
				   iPhoton->trkSumPtHollowConeDR03) << endl;
      cout << "--------------";
      cout << (susyCategory->getPassCombinedIsoMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "--------------";
      cout << (susyCategory->getPassFakeCombinedIsoMax(tag_, i) ? "PASSED\n" : "FAILED\n");
      cout << "-------sieie: " << iPhoton->sigmaIetaIeta << endl;
      cout << "--------------";
      cout << (susyCategory->getPassSigmaIetaIetaMax(tag_, i) ? "PASSED TIGHT\n" : "FAILED TIGHT\n");
      cout << (susyCategory->getPassHLTSigmaIetaIetaMax(tag_, i) ? "PASSED LOOSE\n" : "FAILED LOOSE\n");
      cout << "-------Seed : " << ((iPhoton->nPixelSeeds) > 0 ? "yes\n" : "no\n");
      cout << "--------------";
      cout << (susyCategory->getHasPixelSeed(tag_, i) ? "PASSED\n" : "FAILED\n");
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
  cout << "--------------Pass DPhiMin : ";
  cout << (susyCategory->getPassDPhiMin(tag_) ? "PASSED\n" : "FAILED\n") << endl;
  cout << "--------------Pass DRMin : ";
  cout << (susyCategory->getPassDRMin(tag_) ? "PASSED\n" : "FAILED\n") << endl;
  cout << "--------------Pass asymmetric ETMin: ";
  cout << (susyCategory->getPassAsymmetricETMin(tag_) ? "PASSED\n" : "FAILED\n") << endl;
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
  cout << endl;
  cout << "--------------DiEM ET      : " << susyCategory->getEvtDiEMET(tag_) << endl;
  cout << "--------------InvMass      : " << susyCategory->getEvtInvMass(tag_) << endl;
  cout << endl;
}
