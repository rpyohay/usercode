#include <string>
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
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

float nGen(const string& nGGFile, const unsigned int nTot)
{
  //get no. gg events
  TFile nGGFileStream(nGGFile.c_str());
  if (!nGGFileStream.IsOpen()) {
    cerr << "Error opening file " << nGGFile << ".  Quitting.\n";
    return 0;
  }
  TTree* accTree = NULL;
  nGGFileStream.GetObject("accTree", accTree);
  if (accTree == NULL) {
    cerr << "Error getting tree accTree from file " << nGGFile << ".  Quitting.\n";
    nGGFileStream.Close();
    return 0;
  }
  Float_t nGG = -1.0;
  Float_t nGGErr = -1.0;
  Float_t nGGMETGeq50 = -1.0;
  Float_t nGGMETGeq50Err = -1.0;
  Float_t nGG2Neutralino = -1.0;
  Float_t nGG2NeutralinoErr = -1.0;
  UInt_t n2Neutralino = 0;
  UInt_t nNtuplized = 0;
  accTree->SetBranchAddress("nGG", &nGG);
  accTree->SetBranchAddress("nGGErr", &nGGErr);
  accTree->SetBranchAddress("nGGMETGeq50", &nGGMETGeq50);
  accTree->SetBranchAddress("nGGMETGeq50Err", &nGGMETGeq50Err);
  accTree->SetBranchAddress("nGG2Neutralino", &nGG2Neutralino);
  accTree->SetBranchAddress("nGG2NeutralinoErr", &nGG2NeutralinoErr);
  accTree->SetBranchAddress("n2Neutralino", &n2Neutralino);
  accTree->SetBranchAddress("nNtuplized", &nNtuplized);
  if (accTree->GetEntries() != 1) {
    cerr << "Error: no. entries in accTree is " << accTree->GetEntries() << ", not 1.  Quitting.\n";
    nGGFileStream.Close();
    return 0;
  }
  accTree->GetEntry(0);
  nGGFileStream.Close();

  //return no. generated events (corrected for processing failures)
  if (nNtuplized > 0) return nTot*n2Neutralino/nNtuplized;
  else return 0;
}

pair<TH1F*, TH1F*> signalMETHists(TFile& nGGFileStream)
{
  TH1F* ggMET = NULL;
  TH1F* ffMET = NULL;
  nGGFileStream.GetObject("ggMET", ggMET);
  nGGFileStream.GetObject("ffMET", ffMET);
  if ((ggMET == NULL) || (ffMET == NULL)) {
    cerr << "Error getting histogram ggMET (address " << ggMET << ") or ffMET (address " << ffMET << ") from file ";
    cerr << nGGFileStream.GetName() << ".  Quitting.\n";
    return pair<TH1F*, TH1F*>(NULL, NULL);
  }
  return pair<TH1F*, TH1F*>((TH1F*)ggMET->Clone(), (TH1F*)ffMET->Clone());
}

void nGGRelStatErr(const string& nGGFile, vector<vector<float> >& relStatErr, const vector<float>& METBins)
{
  //get no. gg events
  TFile nGGFileStream(nGGFile.c_str());
  if (!nGGFileStream.IsOpen()) {
    cerr << "Error opening file " << nGGFile << ".  Quitting.\n";
    return;
  }
  TH1F* ggMET = NULL;
  nGGFileStream.GetObject("ggMET", ggMET);
  if (ggMET == NULL) {
    cerr << "Error getting histogram ggMET from file " << nGGFile << ".  Quitting.\n";
    nGGFileStream.Close();
    return;
  }
  for (vector<float>::const_iterator iMETBin = METBins.begin(); iMETBin != METBins.end(); ++iMETBin) {
    Int_t bin1 = ggMET->FindBin(*iMETBin);
    Int_t bin2 = bin1;
    if (iMETBin == (METBins.end() - 1)) bin2 = -1;
    TArrayD* sumW2Array = ggMET->GetSumw2();
    Float_t nGG = ggMET->Integral(bin1, bin2);
    Float_t nGGErr = sumW2Array->At(bin1);
    if (iMETBin == (METBins.end() - 1)) nGGErr+=sumW2Array->At(bin1 + 1);
    nGGErr = sqrt(nGGErr);

    //add nGG rel. stat. error to vector
    if (nGG == 0) relStatErr[iMETBin - METBins.begin()].push_back(0.0);
    else relStatErr[iMETBin - METBins.begin()].push_back((float)(nGGErr/nGG));

    if ((float)(nGGErr/nGG) > 0.4) {
      cout << (nGGErr/nGG) << " for MET bin " << (iMETBin - METBins.begin()) << " and file " << nGGFileStream.GetName() << endl;
    }
  }
  nGGFileStream.Close();
}

pair<Float_t, Float_t> acceptance(const string& nGGFile, const unsigned int nTot, const bool includeMETCut = true)
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
  UInt_t nNtuplized = 0;
  accTree->SetBranchAddress("nGG", &nGG);
  accTree->SetBranchAddress("nGGErr", &nGGErr);
  accTree->SetBranchAddress("nGGMETGeq50", &nGGMETGeq50);
  accTree->SetBranchAddress("nGGMETGeq50Err", &nGGMETGeq50Err);
  accTree->SetBranchAddress("nGG2Neutralino", &nGG2Neutralino);
  accTree->SetBranchAddress("nGG2NeutralinoErr", &nGG2NeutralinoErr);
  accTree->SetBranchAddress("n2Neutralino", &n2Neutralino);
  accTree->SetBranchAddress("nNtuplized", &nNtuplized);
  if (accTree->GetEntries() != 1) {
    cerr << "Error: no. entries in accTree is " << accTree->GetEntries() << ", not 1.  Quitting.\n";
    nGGFileStream.Close();
    return pair<Float_t, Float_t>(-1.0, -1.0);
  }
  accTree->GetEntry(0);
  nGGFileStream.Close();

  //calculate no. generated events corrected for processing failures
  float nGen = 0;
  if (nNtuplized > 0) nGen = nTot*n2Neutralino/nNtuplized;

  //return acceptance and error
  if (includeMETCut) {
    if ((nGen > 0.0) && (nGGMETGeq50Err > 0.0)) {
      return pair<Float_t, Float_t>(nGGMETGeq50/(Float_t)nGen, (nGGMETGeq50/(Float_t)nGen)/nGGMETGeq50Err);
    }
    else return pair<Float_t, Float_t>(-1.0, -1.0);
  }
  else if ((nGen > 0.0) && (nGGErr > 0.0)) return pair<Float_t, Float_t>(nGG/(Float_t)nGen, (nGG/(Float_t)nGen)/nGGErr);
  else return pair<Float_t, Float_t>(-1.0, -1.0);
}

void getNGen(const vector<string>& accFiles, const string& scan, ofstream& out, const unsigned int nTot)
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
    if ((underscore3Pos == string::npos) && (scan != "WB")) {
      cerr << "Error: \"_\" not found after position " << (underscore2Pos + 1) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    int par2 = -1;
    if (scan != "WB") par2 = atoi(iFile->substr(underscore2Pos + 1, underscore3Pos - underscore2Pos - 1).c_str());
    else underscore3Pos = underscore2Pos;
    size_t dotPos = iFile->find('.', underscore3Pos + 1);
    if (dotPos == string::npos) {
      cerr << "Error: \".\" not found after position " << (underscore3Pos + 1) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    const int par3 = atoi(iFile->substr(underscore3Pos + 1, dotPos - underscore3Pos - 1).c_str());
    if (scan == "WB") par2 = par3;

    //get no. generated events
    float theNGen = nGen(*iFile, nTot);

    //rel. JES errors on no. gg events (derived from 100 toys on model ... )
    float JESErr[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

    //write no. generated events and JES error per MET bin to text file
    if (scan != "WB") {
      out <<  par2 << " " << par1 << " " << par3 << " " << theNGen << " ";
      for (unsigned int i = 0; i < 5; ++i) out << JESErr[i] << " ";
      out << endl;
    }
    else {
      out << par1 << " " << par2 << " " << theNGen << " ";
      for (unsigned int i = 0; i < 5; ++i) out << JESErr[i] << " ";
      out << endl;
    }
  }
}

void getSignalMETHists(const vector<string>& accFiles, const string& scan, const string& jet, TFile& out)
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
    if ((underscore3Pos == string::npos) && (scan != "WB")) {
      cerr << "Error: \"_\" not found after position " << (underscore2Pos + 1) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    int par2 = -1;
    if (scan != "WB") par2 = atoi(iFile->substr(underscore2Pos + 1, underscore3Pos - underscore2Pos - 1).c_str());
    else underscore3Pos = underscore2Pos;
    size_t dotPos = iFile->find('.', underscore3Pos + 1);
    if (dotPos == string::npos) {
      cerr << "Error: \".\" not found after position " << (underscore3Pos + 1) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    const int par3 = atoi(iFile->substr(underscore3Pos + 1, dotPos - underscore3Pos - 1).c_str());
    if (scan == "WB") par2 = par3;

    //open file
    TFile nGGFileStream(iFile->c_str());
    if (!nGGFileStream.IsOpen()) cerr << "Error opening file " << *iFile << ".  Quitting.\n";

    //get the gg and ff MET histograms for this grid point
    pair<TH1F*, TH1F*> hists = signalMETHists(nGGFileStream/**iFile, ggMETClone, ffMETClone*/);

    //change the name of the histograms to comply with the limit setting code
    stringstream ggMETCloneName;
    stringstream ffMETCloneName;
    if ((hists.first != NULL) && (hists.second != NULL)) {
      if (scan != "WB") {
	ggMETCloneName << "gg_met_" << jet << "_mS" << par1 << "_mG" << par2 << "_mN" << par3;
	ffMETCloneName << "ff_met_" << jet << "_mS" << par1 << "_mG" << par2 << "_mN" << par3;
      }
      else {
	ggMETCloneName << "gg_met_" << jet << "_mW" << par1 << "_mB" << par2;
	ffMETCloneName << "ff_met_" << jet << "_mW" << par1 << "_mB" << par2;
      }
      hists.first->SetName(ggMETCloneName.str().c_str());
      hists.second->SetName(ffMETCloneName.str().c_str());
      out.cd();
      hists.first->Write();
      hists.second->Write();
    }

    //close file
    nGGFileStream.Close();
  }
}

void fillNGGRelStatErrVec(vector<vector<float> >& relStatErr, const vector<string>& accFiles, const vector<float>& METBins)
{
  //loop over files
  for (vector<string>::const_iterator iFile = accFiles.begin(); iFile != accFiles.end(); ++iFile) {

    //fill nGG rel. stat. error vector
    nGGRelStatErr(*iFile, relStatErr, METBins);
  }
}

void makeAcceptancePlot(TH2F* accHist, const vector<string>& accFiles, const string& scan, const unsigned int nTot, 
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
    if ((underscore3Pos == string::npos) && (scan != "WB")) {
      cerr << "Error: \"_\" not found after position " << (underscore2Pos + 1) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    int par2 = -1;
    if (scan != "WB") par2 = atoi(iFile->substr(underscore2Pos + 1, underscore3Pos - underscore2Pos - 1).c_str());
    else underscore3Pos = underscore2Pos;
    size_t dotPos = iFile->find('.', underscore3Pos + 1);
    if (dotPos == string::npos) {
      cerr << "Error: \".\" not found after position " << (underscore3Pos + 1) << " in " << *iFile << ".  Quitting.\n";
      return;
    }
    const int par3 = atoi(iFile->substr(underscore3Pos + 1, dotPos - underscore3Pos - 1).c_str());
    if (scan == "WB") par2 = par3;

    //find bin corresponding to these parameters
    Int_t xBin = -1;
    Int_t yBin = -1;
    if ((scan == "gsq_B") || (scan == "gsq_W") || (scan == "WB")) {
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
    pair<Float_t, Float_t> acc = acceptance(*iFile, nTot, includeMETCut);

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
			 const unsigned int nTot, const bool includeMETCut = true)
{
  //parse meaning of arguments
  string scan;
  if ((par1Latex == "m_{#tilde{q}}") && (par2Latex == "M_{3}") && (par3Latex == "M_{1}")) scan = "gsq_B";
  else if ((par1Latex == "m_{#tilde{q}}") && (par2Latex == "M_{3}") && (par3Latex == "M_{2}")) scan = "gsq_W";
  else if ((par1Latex == "M_{1}") && (par2Latex == "M_{3}") && (par3Latex == "m_{#tilde{q}}")) scan = "gB";
  else if ((par1Latex == "M_{2}") && (par2Latex == "M_{1}") && (par3Latex == "")) scan = "WB";
  else {
    cerr << "Error: unrecognized parameter set (\"" << par1Latex << "\", \"" << par2Latex << "\", \"" << par3Latex << "\").  Quitting.\n";
    return NULL;
  }

  //fill acceptance histogram
  string accHistName = scan + "_accHist";
  if ((accFiles.size() > 0) && ((accFiles[0].find("binolikegrid2") != string::npos) || 
				(accFiles[0].find("winolikegrid2") != string::npos) || 
				(accFiles[0].find("binochigrids") != string::npos))) accHistName = scan + "Old_accHist";
  TH2F* accHist = 
    new TH2F(accHistName.c_str(), "", nBins.first, minBin.first, maxBin.first, nBins.second, minBin.second, maxBin.second);
  setHistogramOptions(accHist, kBlack, 0.7, 20, 1.0, 1.0, par1Latex.c_str(), par2Latex.c_str());
  makeAcceptancePlot(accHist, accFiles, scan, nTot, includeMETCut);
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

void makeNGenFiles(const vector<string>& outputFiles, const vector<vector<string> >& accFiles, const vector<string>& scans)
{
  //loop over output files
  for (vector<string>::const_iterator iOut = outputFiles.begin(); iOut != outputFiles.end(); ++iOut) {
    ofstream out(iOut->c_str());
    if (!out.is_open()) {
      cerr << "Error opening file " << *iOut << ".  Quitting.\n";
      return;
    }
    const unsigned int i = iOut - outputFiles.begin();

    //get the number of generated events
    unsigned int nGen = 10000;
    if ((accFiles[i].size() > 0) && (accFiles[i][0].find("winolikegrid2") != string::npos)) nGen = 60000;

    //fill the file of numbers of generated events per grid point
    getNGen(accFiles[i], scans[i], out, nGen);

    //exit
    out.close();
  }
}

void makeFilesForLimitSetting(const vector<string>& outputFiles, const vector<vector<string> >& accFiles, const vector<string>& scans, 
			      const string& jet)
{
  //loop over output files
  for (vector<string>::const_iterator iOut = outputFiles.begin(); iOut != outputFiles.end(); ++iOut) {
    TFile out(iOut->c_str(), "UPDATE");
    if (!out.IsOpen()) {
      cerr << "Error opening file " << *iOut << ".  Quitting.\n";
      return;
    }
    const unsigned int i = iOut - outputFiles.begin();

    //get the histograms needed for limit setting (they should automatically write to the output file)
    getSignalMETHists(accFiles[i], scans[i], jet, out);

    //exit
    out.Write();
    out.Close();
  }
}

void sortNGGRelStatErr(const vector<string>& squarkGluinoAccFiles, const vector<string>& binoGluinoAccFiles, 
		       const vector<string>& winoBinoAccFiles, const vector<float>& METBins)
{
  //fill vectors of nGG rel. stat. error
  vector<vector<float> > squarkGluinoNGGRelStatErr(METBins.size(), vector<float>());
  vector<vector<float> > binoGluinoNGGRelStatErr(METBins.size(), vector<float>());
  vector<vector<float> > winoBinoNGGRelStatErr(METBins.size(), vector<float>());
  fillNGGRelStatErrVec(squarkGluinoNGGRelStatErr, squarkGluinoAccFiles, METBins);
  fillNGGRelStatErrVec(binoGluinoNGGRelStatErr, binoGluinoAccFiles, METBins);
  fillNGGRelStatErrVec(winoBinoNGGRelStatErr, winoBinoAccFiles, METBins);

  //sort vectors from smallest to largest element
  for (vector<float>::const_iterator iMETBin = METBins.begin(); iMETBin != METBins.end(); ++iMETBin) {
    const unsigned int i = iMETBin - METBins.begin();
    sort(squarkGluinoNGGRelStatErr[i].begin(), squarkGluinoNGGRelStatErr[i].end());
    sort(binoGluinoNGGRelStatErr[i].begin(), binoGluinoNGGRelStatErr[i].end());
    sort(winoBinoNGGRelStatErr[i].begin(), winoBinoNGGRelStatErr[i].end());
  }

  //remove elements with value 0
  for (vector<float>::const_iterator iMETBin = METBins.begin(); iMETBin != METBins.end(); ++iMETBin) {
    const unsigned int i = iMETBin - METBins.begin();
    while ((squarkGluinoNGGRelStatErr[i][0] == 0) && (squarkGluinoNGGRelStatErr[i].size() > 0)) {
      squarkGluinoNGGRelStatErr[i].erase(squarkGluinoNGGRelStatErr[i].begin());
    }
    while ((binoGluinoNGGRelStatErr[i][0] == 0) && (binoGluinoNGGRelStatErr[i].size() > 0)) {
      binoGluinoNGGRelStatErr[i].erase(binoGluinoNGGRelStatErr[i].begin());
    }
    while ((winoBinoNGGRelStatErr[i][0] == 0) && (winoBinoNGGRelStatErr[i].size() > 0)) {
      winoBinoNGGRelStatErr[i].erase(winoBinoNGGRelStatErr[i].begin());
    }
  }

  //print smallest and largest elements of each vector
  for (vector<float>::const_iterator iMETBin = METBins.begin(); iMETBin != METBins.end(); ++iMETBin) {
    const unsigned int i = iMETBin - METBins.begin();
    cout << "nGG rel. stat. error ranges from " << *(squarkGluinoNGGRelStatErr[i].begin()) << " to ";
    cout << *(squarkGluinoNGGRelStatErr[i].end() - 1) << " for scan gsq_B and MET bin " << i << endl;
    cout << "nGG rel. stat. error ranges from " << *(binoGluinoNGGRelStatErr[i].begin()) << " to ";
    cout << *(binoGluinoNGGRelStatErr[i].end() - 1) << " for scan gsq_W and MET bin " << i << endl;
    cout << "nGG rel. stat. error ranges from " << *(winoBinoNGGRelStatErr[i].begin()) << " to " << *(winoBinoNGGRelStatErr[i].end() - 1);
    cout << " for scan gB and MET bin " << i << endl;
  }
}

void writeAcceptanceFile(const string& outputFile, const pair<Int_t, Int_t>& squarkGluinoNBins, const pair<Int_t, Int_t>& binoGluinoNBins, 
			 const pair<Int_t, Int_t>& winoBinoNBins, const pair<Int_t, Int_t>& squarkGluinoBinoOldNBins, 
			 const pair<Int_t, Int_t>& squarkGluinoWinoOldNBins, const pair<Int_t, Int_t>& binoGluinoOldNBins, 
			 const pair<Double_t, Double_t>& squarkGluinoMinBin, const pair<Double_t, Double_t>& binoGluinoMinBin, 
			 const pair<Double_t, Double_t>& winoBinoMinBin, const pair<Double_t, Double_t>& squarkGluinoBinoOldMinBin, 
			 const pair<Double_t, Double_t>& squarkGluinoWinoOldMinBin, const pair<Double_t, Double_t>& binoGluinoOldMinBin, 
			 const pair<Double_t, Double_t>& squarkGluinoMaxBin, const pair<Double_t, Double_t>& binoGluinoMaxBin, 
			 const pair<Double_t, Double_t>& winoBinoMaxBin, const pair<Double_t, Double_t>& squarkGluinoBinoOldMaxBin, 
			 const pair<Double_t, Double_t>& squarkGluinoWinoOldMaxBin, const pair<Double_t, Double_t>& binoGluinoOldMaxBin, 
			 const vector<string>& squarkGluinoAccFiles, const vector<string>& binoGluinoAccFiles, 
			 const vector<string>& winoBinoAccFiles, const vector<string>& squarkGluinoBinoOldAccFiles, 
			 const vector<string>& squarkGluinoWinoOldAccFiles, const vector<string>& binoGluinoOldAccFiles)
{
  //open output file
  TFile out(outputFile.c_str(), "RECREATE");
  if (!out.IsOpen()) {
    cerr << "Error opening file " << outputFile << ".  Quitting.\n";
    return;
  }

  //make acceptance plots
  TH2F* squarkGluinoAccHist = fillAcceptanceGrid("m_{#tilde{q}}", "M_{3}", "M_{1}", squarkGluinoNBins, 
						 squarkGluinoMinBin, squarkGluinoMaxBin, squarkGluinoAccFiles, 10000);
  TH2F* binoGluinoAccHist = fillAcceptanceGrid("M_{1}", "M_{3}", "m_{#tilde{q}}", binoGluinoNBins, binoGluinoMinBin, 
					       binoGluinoMaxBin, binoGluinoAccFiles, 10000);
  TH2F* winoBinoAccHist = fillAcceptanceGrid("M_{2}", "M_{1}", "", winoBinoNBins, winoBinoMinBin, winoBinoMaxBin, 
					     winoBinoAccFiles, 10000);
  TH2F* squarkGluinoBinoOldAccHist = fillAcceptanceGrid("m_{#tilde{q}}", "M_{3}", "M_{1}", squarkGluinoBinoOldNBins, 
							squarkGluinoBinoOldMinBin, squarkGluinoBinoOldMaxBin, squarkGluinoBinoOldAccFiles, 
							10000);
  TH2F* squarkGluinoWinoOldAccHist = fillAcceptanceGrid("m_{#tilde{q}}", "M_{3}", "M_{2}", squarkGluinoWinoOldNBins, 
							squarkGluinoWinoOldMinBin, squarkGluinoWinoOldMaxBin, squarkGluinoWinoOldAccFiles, 
							60000);
  TH2F* binoGluinoOldAccHist = fillAcceptanceGrid("M_{1}", "M_{3}", "m_{#tilde{q}}", binoGluinoOldNBins, binoGluinoOldMinBin, 
						  binoGluinoOldMaxBin, binoGluinoOldAccFiles, 10000);

  //make acceptance canvases
  TCanvas squarkGluinoAccCanvas("squarkGluinoAccCanvas", "", 600, 600);
  TCanvas binoGluinoAccCanvas("binoGluinoAccCanvas", "", 600, 600);
  TCanvas winoBinoAccCanvas("winoBinoAccCanvas", "", 600, 600);
  TCanvas squarkGluinoBinoOldAccCanvas("squarkGluinoBinoOldAccCanvas", "", 600, 600);
  TCanvas squarkGluinoWinoOldAccCanvas("squarkGluinoWinoOldAccCanvas", "", 600, 600);
  TCanvas binoGluinoOldAccCanvas("binoGluinoOldAccCanvas", "", 600, 600);
  setCanvasOptions(&squarkGluinoAccCanvas, 0, 0, 0);
  setCanvasOptions(&binoGluinoAccCanvas, 0, 0, 0);
  setCanvasOptions(&winoBinoAccCanvas, 0, 0, 0);
  setCanvasOptions(&squarkGluinoBinoOldAccCanvas, 0, 0, 0);
  setCanvasOptions(&squarkGluinoWinoOldAccCanvas, 0, 0, 0);
  setCanvasOptions(&binoGluinoOldAccCanvas, 0, 0, 0);
  squarkGluinoAccCanvas.cd();
  if (squarkGluinoAccHist != NULL) squarkGluinoAccHist->Draw("COLZ");
  binoGluinoAccCanvas.cd();
  if (binoGluinoAccHist != NULL) binoGluinoAccHist->Draw("COLZ");
  winoBinoAccCanvas.cd();
  if (winoBinoAccHist != NULL) winoBinoAccHist->Draw("COLZ");
  squarkGluinoBinoOldAccCanvas.cd();
  if (squarkGluinoBinoOldAccHist != NULL) squarkGluinoBinoOldAccHist->Draw("COLZ");
  squarkGluinoWinoOldAccCanvas.cd();
  if (squarkGluinoWinoOldAccHist != NULL) squarkGluinoWinoOldAccHist->Draw("COLZ");
  binoGluinoOldAccCanvas.cd();
  if (binoGluinoOldAccHist != NULL) binoGluinoOldAccHist->Draw("COLZ");

  //exit
  out.cd();
  squarkGluinoAccCanvas.Write();
  binoGluinoAccCanvas.Write();
  winoBinoAccCanvas.Write();
  squarkGluinoBinoOldAccCanvas.Write();
  squarkGluinoWinoOldAccCanvas.Write();
  binoGluinoOldAccCanvas.Write();
  out.Write();
  out.Close();
}
