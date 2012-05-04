#include <string>
#include <iostream>
#include <utility>
#include <vector>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"

Int_t findBinForParameter(TH1D* parProjHist, const Int_t par)
{
  Int_t bin = 1;
  bool foundBin = false;
  Int_t iBin = 1;
  while ((iBin <= parProjHist->GetNbinsX()) && !foundBin) {
    Double_t lowEdge = parProjHist->GetBinLowEdge(iBin);
    if (((Double_t)par >= lowEdge) && ((Double_t)par < (lowEdge + parProjHist->GetBinWidth(iBin)))) {
      bin = iBin;
      foundBin = true;
    }
    ++iBin;
  }
  if (foundBin == false) return 0;
  else return bin;
}

pair<Float_t, Float_t> acceptance(const string& nGGFile, const int nTot = -1, const bool includeMETCut = true)
{
  //get no. gg events
  TFile nGGFileStream(nGGFile.c_str());
  if (!nGGFileStream.IsOpen()) {
    cerr << "Error opening file " << nGGFile << ".  Quitting.\n";
    return pair<Float_t, Float_t>(-1.0, -1.0);
  }
  TTree* accTree = NULL;
  nGGFileStream.GetObject("accTree", accTree);
  if (accTree == NULL) {
    cerr << "Error getting tree accTree from file " << nGGFile << ".  Quitting.\n";
    nGGFileStream.Close();
    return pair<Float_t, Float_t>(-1.0, -1.0);
  }
  Float_t nGG = -1.0;
  Float_t nGGErr = -1.0;
  Float_t nGGMETGeq50 = -1.0;
  Float_t nGGMETGeq50Err = -1.0;
  Float_t nGG2Neutralino = -1.0;
  Float_t nGG2NeutralinoErr = -1.0;
  UInt_t n2Neutralino = 0;
  accTree->SetBranchAddress("nGG", &nGG);
  accTree->SetBranchAddress("nGGErr", &nGGErr);
  accTree->SetBranchAddress("nGGMETGeq50", &nGGMETGeq50);
  accTree->SetBranchAddress("nGGMETGeq50Err", &nGGMETGeq50Err);
  accTree->SetBranchAddress("nGG2Neutralino", &nGG2Neutralino);
  accTree->SetBranchAddress("nGG2NeutralinoErr", &nGG2NeutralinoErr);
  accTree->SetBranchAddress("n2Neutralino", &n2Neutralino);
  if (accTree->GetEntries() != 1) {
    cerr << "Error: no. entries in accTree is " << accTree->GetEntries() << ", not 1.  Quitting.\n";
    nGGFileStream.Close();
    return pair<Float_t, Float_t>(-1.0, -1.0);
  }
  accTree->GetEntry(0);
  nGGFileStream.Close();

  //return acceptance and error
  if (nTot == -1) {
    return pair<Float_t, Float_t>(nGG2Neutralino/(Float_t)n2Neutralino, (nGG2Neutralino/(Float_t)n2Neutralino)/nGG2NeutralinoErr);
  }
  else if (includeMETCut) {
    cout << "Including MET cut...\n";
    return pair<Float_t, Float_t>(nGGMETGeq50/(Float_t)nTot, (nGGMETGeq50/(Float_t)nTot)/nGGMETGeq50Err);
  }
  else return pair<Float_t, Float_t>(nGG/(Float_t)nTot, (nGG/(Float_t)nTot)/nGGErr);
}

void makeAcceptancePlot(TH2F* accHist, const vector<string>& accFiles, const string& scan, const int nTot = -1, 
			const bool includeMETCut = true)
{
  //loop over files
  for (vector<string>::const_iterator iFile = accFiles.begin(); iFile != accFiles.end(); ++iFile) {

    //get parameters from file name
    string accPhrase("acceptance_");
    size_t acceptancePos = iFile->find(accPhrase);
    if (acceptancePos == string::npos) {
      cerr << "Error: \"" << accPhrase << "\" not found after position 0 in " << *iFile << ".  Quitting.\n";
      return;
    }
    size_t underscore2Pos = iFile->find('_', acceptancePos + accPhrase.length());
    if (underscore2Pos == string::npos) {
      cerr << "Error: \"_\" not found after position " << (acceptancePos + accPhrase.length()) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    const int par1 = atoi(iFile->substr(acceptancePos + accPhrase.length(), underscore2Pos - acceptancePos - accPhrase.length()).c_str());
    size_t underscore3Pos = iFile->find('_', underscore2Pos + 1);
    if (underscore3Pos == string::npos) {
      cerr << "Error: \"_\" not found after position " << (underscore2Pos + 1) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    const int par2 = atoi(iFile->substr(underscore2Pos + 1, underscore3Pos - underscore2Pos - 1).c_str());
    size_t dotPos = iFile->find('.', underscore3Pos + 1);
    if (dotPos == string::npos) {
      cerr << "Error: \".\" not found after position " << (underscore3Pos + 1) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    const int par3 = atoi(iFile->substr(underscore3Pos + 1, dotPos - underscore3Pos - 1).c_str());

    //find bin corresponding to these parameters
    Int_t xBin = -1;
    Int_t yBin = -1;
    if ((scan == "gsq_B") || (scan == "WB")) {
      xBin = findBinForParameter(accHist->ProjectionX(), par1);
      yBin = findBinForParameter(accHist->ProjectionY(), par2);
    }
    else if (scan == "gB") {
      xBin = findBinForParameter(accHist->ProjectionX(), par3);
      yBin = findBinForParameter(accHist->ProjectionY(), par2);
    }
    else {
      cerr << "Error: unrecognized scan \"" << scan << "\".  Quitting.\n";
      return;
    }

    //check that this bin wasn't already filled
    if (accHist->GetBinContent(xBin, yBin) != 0.0) {
      cerr << "Error: bins (" << xBin << ", " << yBin << ") corresponding to parameters (" << par1 << ", " << par2;
      cerr << ") are already filled.  Quitting\n.";
      return;
    }

    //calculate acceptances
    pair<Float_t, Float_t> acc = acceptance(*iFile, nTot);

    //fill acceptance plot
    accHist->SetBinContent(xBin, yBin, acc.first);
    accHist->SetBinError(xBin, yBin, acc.second);
  }
}

void setAxisOptions(TAxis* axis, const Float_t labelSize, const Float_t titleOffset, 
		    const char* title)
{
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetTitleSize(0.06);
  axis->SetTitleOffset(titleOffset);
  axis->SetTitle(title);
}

void setHistogramOptions(TH2F* histogram, const Color_t color, const Size_t size, 
			 const Style_t style, const Float_t yAxisTitleOffset, const Double_t scale, 
			 const char* xAxisTitle, const char* yAxisTitle)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.03, 0.9, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.03, yAxisTitleOffset, yAxisTitle);
  setAxisOptions(histogram->GetZaxis(), 0.03, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

TH2F* fillAcceptanceGrid(const string& par1Latex, const string& par2Latex, const string& par3Latex, const pair<Int_t, Int_t>& nBins, 
			 const pair<Double_t, Double_t>& minBin, const pair<Double_t, Double_t>& maxBin, const vector<string>& accFiles, 
			 const int nTot = -1, const bool includeMETCut = true)
{
  //parse meaning of arguments
  string scan;
  if ((par1Latex == "m_{#tilde{q}}") && (par2Latex == "m_{#tilde{g}}") && (par3Latex == "m_{#tilde{#chi}_{1}^{0}}")) scan = "gsq_B";
  else if ((par1Latex == "m_{#tilde{B}}") && (par2Latex == "m_{#tilde{g}}") && (par3Latex == "m_{#tilde{q}}")) scan = "gB";
  else if ((par1Latex == "m_{#tilde{W}}") && (par2Latex == "m_{#tilde{B}}") && (par3Latex == "")) scan = "WB";
  else {
    cerr << "Error: unrecognized parameter set (\"" << par1Latex << "\", \"" << par2Latex << "\", \"" << par3Latex << "\").  Quitting.\n";
    return NULL;
  }

  //fill acceptance histogram
  TH2F* accHist = 
    new TH2F((scan + "_accHist").c_str(), "", nBins.first, minBin.first, maxBin.first, nBins.second, minBin.second, maxBin.second);
  setHistogramOptions(accHist, kBlack, 0.7, 20, 1.0, 1.0, par1Latex.c_str(), par2Latex.c_str());
  makeAcceptancePlot(accHist, accFiles, scan, nTot);
  return accHist;
}

void setCanvasOptions(TCanvas* canvas, const Int_t grid, const Int_t logY, const Int_t logZ)
{
  canvas->SetFillStyle(0);
  canvas->SetFillColor(0);
  canvas->SetGrid(grid, grid);
  canvas->SetLogy(logY);
  canvas->SetLogz(logZ);
  canvas->SetLeftMargin(0.2);
  canvas->SetTopMargin(0.2);
  canvas->SetRightMargin(0.2);
  canvas->SetBottomMargin(0.2);
}

void writeAcceptanceFile(const string& outputFile, const pair<Int_t, Int_t>& squarkGluinoNBins, const pair<Int_t, Int_t>& binoGluinoNBins, 
			 const pair<Int_t, Int_t>& winoBinoNBins, const pair<Double_t, Double_t>& squarkGluinoMinBin, 
			 const pair<Double_t, Double_t>& binoGluinoMinBin, const pair<Double_t, Double_t>& winoBinoMinBin, 
			 const pair<Double_t, Double_t>& squarkGluinoMaxBin, const pair<Double_t, Double_t>& binoGluinoMaxBin, 
			 const pair<Double_t, Double_t>& winoBinoMaxBin, const vector<string>& squarkGluinoAccFiles, 
			 const vector<string>& binoGluinoAccFiles, const vector<string>& winoBinoAccFiles)
{
  //open output file
  TFile out(outputFile.c_str(), "RECREATE");
  if (!out.IsOpen()) {
    cerr << "Error opening file " << outputFile << ".  Quitting.\n";
    return;
  }

  //make acceptance plots
  TH2F* squarkGluinoAccHist = fillAcceptanceGrid("m_{#tilde{q}}", "m_{#tilde{g}}", "m_{#tilde{#chi}_{1}^{0}}", squarkGluinoNBins, 
						 squarkGluinoMinBin, squarkGluinoMaxBin, squarkGluinoAccFiles, 10000);
  TH2F* binoGluinoAccHist = fillAcceptanceGrid("m_{#tilde{B}}", "m_{#tilde{g}}", "m_{#tilde{q}}", binoGluinoNBins, binoGluinoMinBin, 
					       binoGluinoMaxBin, binoGluinoAccFiles, 10000);
  TH2F* winoBinoAccHist = fillAcceptanceGrid("m_{#tilde{W}}", "m_{#tilde{B}}", "", winoBinoNBins, winoBinoMinBin, winoBinoMaxBin, 
					     winoBinoAccFiles, 10000);

  //make acceptance canvases
  TCanvas squarkGluinoAccCanvas("squarkGluinoAccCanvas", "", 600, 600);
  TCanvas binoGluinoAccCanvas("binoGluinoAccCanvas", "", 600, 600);
  TCanvas winoBinoAccCanvas("winoBinoAccCanvas", "", 600, 600);
  setCanvasOptions(&squarkGluinoAccCanvas, 0, 0, 0);
  setCanvasOptions(&binoGluinoAccCanvas, 0, 0, 0);
  setCanvasOptions(&winoBinoAccCanvas, 0, 0, 0);
  squarkGluinoAccCanvas.cd();
  if (squarkGluinoAccHist != NULL) squarkGluinoAccHist->Draw("COLZ");
  binoGluinoAccCanvas.cd();
  if (binoGluinoAccHist != NULL) binoGluinoAccHist->Draw("COLZ");
  winoBinoAccCanvas.cd();
  if (winoBinoAccHist != NULL) winoBinoAccHist->Draw("COLZ");

  //exit
  out.cd();
  squarkGluinoAccCanvas.Write();
  binoGluinoAccCanvas.Write();
  winoBinoAccCanvas.Write();
  out.Write();
  out.Close();
}
