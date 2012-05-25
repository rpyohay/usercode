#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "PhysicsTools/TagAndProbe/interface/RooCMSShape.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace RooFit;

void doIt()
{
  //input file names
  vector<string> inputFileNames;
  inputFileNames.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_diEMPTReweighting.root");
  inputFileNames.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_dijetPTReweighting.root");
  inputFileNames.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_dijetPTReweighting_noPFJetID.root");
  inputFileNames.push_back("/Users/rachelyohay/RA3/data/4684pb-1_MET_dijetPTReweighting_noPFJetID_L1FastOnly.root");
  vector<TFile*> inputFiles;

  //legend entries
  vector<string> legendEntries;
  legendEntries.push_back("Di-EM p_{T} reweighting");
  legendEntries.push_back("Dijet p_{T} reweighting");
  legendEntries.push_back("Dijet p_{T} reweighting, no jet ID");
  legendEntries.push_back("Dijet p_{T} reweighting, no jet ID, L1FastL2L3/L2L3 JEC");

  //sample names
  vector<string> sampleNames;
  sampleNames.push_back("gg");
  sampleNames.push_back("ee");
  sampleNames.push_back("ff");

  //plots
  vector<vector<TH1F*> > PT(3, vector<TH1F*>());
  vector<vector<TH1F*> > MET(3, vector<TH1F*>());

  //loop over input files
  for (vector<string>::const_iterator iIn = inputFileNames.begin(); iIn != inputFileNames.end(); 
       ++iIn) {

    //open input file
    inputFiles.push_back(new TFile(iIn->c_str()));
    if ((*(inputFiles.end() - 1))->IsOpen()) {

      //get plots
      for (vector<string>::const_iterator iSample = sampleNames.begin(); 
	   iSample != sampleNames.end(); ++iSample) {
	TH1F* pPT;
	TH1F* pMET;
	(*(inputFiles.end() - 1))->GetObject((*iSample + "DiEMETUniform").c_str(), pPT);
	string name(*iSample + "Final");
	if (*iSample == "gg") name = *iSample + "MET";
	(*(inputFiles.end() - 1))->GetObject(name.c_str(), pMET);
	PT[iSample - sampleNames.begin()].push_back(pPT);
	MET[iSample - sampleNames.begin()].push_back(pMET);
      }
    }
    else cerr << "Error: could not open file " << *iIn << ".  Skipping this file.\n";
  }

  //make pT canvases
  for (vector<vector<TH1F*> >::iterator iPT = PT.begin(); iPT != PT.end(); ++iPT) {
    string name(sampleNames[iPT - PT.begin()] + "PT");
    TLegend legend(0.4, 0.7, 0.9, 0.9);
    legend.SetFillStyle(0);
    legend.SetLineColor(0);
    TCanvas canvas(name.c_str(), name.c_str(), 600, 600);
    canvas.SetLogy();
    canvas.cd();
    for (vector<TH1F*>::iterator iMethod = iPT->begin(); iMethod != iPT->end(); ++iMethod) {
      const unsigned int index = iMethod - iPT->begin();
      (*iMethod)->SetMarkerColor(index + 1);
      (*iMethod)->SetMarkerSize(0.5);
      (*iMethod)->SetLineColor(index + 1);
      (*iMethod)->SetTitle("");
      (*iMethod)->SetAxisRange(0.0, 200.0);
      legend.AddEntry(*iMethod, legendEntries[index].c_str(), "lp");
      string drawOpt("E1");
      if (index > 0) {
	drawOpt+="SAME";
	(*iMethod)->Scale((*(iPT->begin()))->GetEntries()/(*iMethod)->GetEntries());
      }
      (*iMethod)->Draw(drawOpt.c_str());
    }
    legend.Draw();
    canvas.SaveAs((name + "_log_zoom.pdf").c_str());
    for (vector<TH1F*>::iterator iMethod = iPT->begin(); iMethod != iPT->end(); ++iMethod) {
      *iMethod = NULL;
    }
  }

  //make MET canvases
  for (vector<vector<TH1F*> >::iterator iMET = MET.begin(); iMET != MET.end(); ++iMET) {
    string sampleName = sampleNames[iMET - MET.begin()];
    if (sampleName != "gg") {
      string name(sampleName + "MET");
      TLegend legend(0.4, 0.7, 0.9, 0.9);
      legend.SetFillStyle(0);
      legend.SetLineColor(0);
      TCanvas canvas(name.c_str(), name.c_str(), 600, 600);
      canvas.SetLogy();
      canvas.cd();
      for (vector<TH1F*>::iterator iMethod = iMET->begin(); iMethod != iMET->end(); ++iMethod) {
	const unsigned int index = iMethod - iMET->begin();
	(*iMethod)->SetMarkerColor(index + 1);
	(*iMethod)->SetLineColor(index + 1);
	(*iMethod)->SetTitle("");
	(*iMethod)->GetXaxis()->SetRange(1, 12);
	legend.AddEntry(*iMethod, legendEntries[index].c_str(), "lp");
	string drawOpt("E1");
	if (index > 0) {
	  drawOpt+="SAME";
	  (*iMethod)->Scale((*(iMET->begin()))->Integral()/(*iMethod)->Integral());
	}
	(*iMethod)->Draw(drawOpt.c_str());
      }
      legend.Draw();
      canvas.SaveAs((name + "_log_zoom.pdf").c_str());
    }
    for (vector<TH1F*>::iterator iMethod = iMET->begin(); iMethod != iMET->end(); ++iMethod) {
      *iMethod = NULL;
    }
  }

  //clean up
  for (vector<TFile*>::iterator iIn = inputFiles.begin(); iIn != inputFiles.end(); ++iIn) {
    delete *iIn;
    *iIn = NULL;
  }
}

void plotIt()
{
  string samples[3] = {"gg", "ee", "ff"};
  TFile RPY("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_correctEESidebandSubtraction.root");
  TFile RPYEEDijetPTReweighting("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_eeDijetPTReweighting.root");
  TFile YFL("/Users/rachelyohay/RA3/data/diempt.root");
  TFile out("/Users/rachelyohay/RA3/data/RPYYFLComparison.root", "RECREATE");
  vector<vector<TH1F*> > dijetPTRPY(3, vector<TH1F*>(3, NULL));
  vector<vector<TH1F*> > dijetPTYFL(3, vector<TH1F*>(3, NULL));
  vector<vector<TH1F*> > weightsRPY(3, vector<TH1F*>(3, NULL));
  vector<vector<TCanvas*> > dijetPT(3, vector<TCanvas*>(3, NULL));
  vector<vector<TCanvas*> > dijetPTRatio(3, vector<TCanvas*>(3, NULL));
  vector<vector<TLegend*> > dijetPTLegend(3, vector<TLegend*>(3, NULL));
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      stringstream nameDijetPT;
      nameDijetPT << samples[i] << "DijetPT" << j << "j";
      dijetPT[i][j] = new TCanvas(nameDijetPT.str().c_str(), "", 600, 600);
      dijetPT[i][j]->SetFillStyle(0);
      dijetPT[i][j]->SetFillColor(0);
      dijetPT[i][j]->SetGrid();
      stringstream nameDijetPTRatio;
      nameDijetPTRatio << samples[i] << "DijetPTRatio" << j << "j";
      dijetPTRatio[i][j] = new TCanvas(nameDijetPTRatio.str().c_str(), "", 600, 300);
      dijetPTRatio[i][j]->SetFillStyle(0);
      dijetPTRatio[i][j]->SetFillColor(0);
      dijetPTRatio[i][j]->SetGrid();
      stringstream nameDijetPTLegend;
      nameDijetPTLegend << samples[i] << " " << j << " jets";
      dijetPTLegend[i][j] = new TLegend(0.6, 0.6, 0.9, 0.8);
      dijetPTLegend[i][j]->SetFillColor(0);
      dijetPTLegend[i][j]->SetHeader(nameDijetPTLegend.str().c_str());
      stringstream nameRPY;
      nameRPY << samples[i] << "DiEMETUniformNjBin" << (j + 1);
      RPY.GetObject(nameRPY.str().c_str(), dijetPTRPY[i][j]);
      stringstream nameWeightsRPY;
      nameWeightsRPY << samples[i] << "Weights_METBin1_NjBin" << (j + 1);
      if (i != 0) RPY.GetObject(nameWeightsRPY.str().c_str(), weightsRPY[i][j]);
      stringstream nameYFL;
      nameYFL << "diempt_" << samples[i] << "_" << j << "jet";
      YFL.GetObject(nameYFL.str().c_str(), dijetPTYFL[i][j]);
      dijetPTLegend[i][j]->AddEntry(dijetPTRPY[i][j], "RPY", "lp");
      dijetPTLegend[i][j]->AddEntry(dijetPTYFL[i][j], "YFL", "lp");
      out.cd();
      if (i != 0) weightsRPY[i][j]->Write();
      dijetPT[i][j]->cd();
      dijetPTRPY[i][j]->SetMarkerColor(kRed);
      dijetPTRPY[i][j]->SetMarkerSize(0.5);
      dijetPTRPY[i][j]->SetLineColor(kRed);
      dijetPTRPY[i][j]->GetXaxis()->SetTitle("Dijet p_{T} (GeV)");
      dijetPTYFL[i][j]->SetMarkerColor(kBlue);
      dijetPTYFL[i][j]->SetMarkerSize(0.5);
      dijetPTYFL[i][j]->SetLineColor(kBlue);
      dijetPTYFL[i][j]->GetXaxis()->SetTitle("Dijet p_{T} (GeV)");
      dijetPTYFL[i][j]->Scale(dijetPTRPY[i][j]->Integral()/dijetPTYFL[i][j]->Integral());
      dijetPTRPY[i][j]->Draw();
      dijetPTYFL[i][j]->Draw("SAME");
      dijetPTLegend[i][j]->Draw();
      dijetPT[i][j]->Write();
      dijetPTRatio[i][j]->cd();
      TH1F* dijetPTRatioDist = (TH1F*)dijetPTRPY[i][j]->Clone();
      dijetPTRatioDist->Divide(dijetPTYFL[i][j]);
      dijetPTRatioDist->SetMarkerColor(kBlack);
      dijetPTRatioDist->SetMarkerSize(0.5);
      dijetPTRatioDist->SetLineColor(kBlack);
      dijetPTRatioDist->GetXaxis()->SetTitle("Dijet p_{T} (GeV)");
      dijetPTRatioDist->GetYaxis()->SetTitle("#frac{N_{RPY}}{N_{YFL}}");
      dijetPTRatioDist->Draw();
      dijetPTRatio[i][j]->Write();
    }
  }
  out.Close();
  RPY.Close();
  RPYEEDijetPTReweighting.Close();
  YFL.Close();
}

void compareEEToGG()
{
  TFile in("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_restrictedDiEMPT.root");
  TFile out("/Users/rachelyohay/RA3/data/eeVsGG_data_restrictedDiEMPT.root", "RECREATE");
  TH3F* gg;
  TH1F* eeEGMET;
  TCanvas* canvas;
  in.GetObject("ggMETVsDiEMETVsInvMassTot", gg);
  in.GetObject("eeFinal", eeEGMET);
  in.GetObject("METCanvas", canvas);
  TH1F* egMET = (TH1F*)canvas->GetPrimitive("egMET");
  TH1F* eeMET = (TH1F*)eeEGMET->Clone();
  eeMET->Add(egMET, -1.0);
  eeMET->SetFillStyle(1001);
  TH1D* ggMET = gg->ProjectionZ("ggMET", 0, -1, 0, -1, "e");
  // TH1D* ggMET = gg->ProjectionZ("ggMET", 0, -1, 6, 9, "e");
  // eeMET->Scale(ggMET->Integral(1, 4)/eeMET->Integral(1, 4));
  TH1D* METRatioDist = (TH1D*)ggMET->Clone();
  METRatioDist->Divide(eeMET);
  METRatioDist->SetMarkerColor(kBlack);
  METRatioDist->SetLineColor(kBlack);
  METRatioDist->GetYaxis()->SetTitle("#frac{N_{#gamma#gamma}}{N_{ee}}");
  TCanvas MET("MET", "", 600, 600);
  MET.SetFillStyle(0);
  MET.SetFillColor(0);
  MET.SetGrid();
  TCanvas METRatio("METRatio", "", 600, 300);
  METRatio.SetFillStyle(0);
  METRatio.SetFillColor(0);
  METRatio.SetGrid();
  TLegend METLegend(0.5, 0.6, 0.9, 0.8);
  METLegend.SetFillColor(0);
  METLegend.SetHeader("CMS 4.7 fb^{-1}, no jet requirement");
  // METLegend.SetHeader("CMS simulation, no jet requirement");
  METLegend.AddEntry(ggMET, "#gamma#gamma", "lp");
  METLegend.AddEntry(eeMET, "ee", "f");
  out.cd();
  MET.cd();
  eeMET->Draw("E2");
  ggMET->Draw("SAME");
  METLegend.Draw();
  MET.Write();
  METRatio.cd();
  METRatioDist->Draw();
  METRatio.Write();
  out.Close();
  in.Close();
}

void setCanvasOptions(TVirtualPad* canvas, const Int_t grid, const Int_t logY, const Int_t logZ)
{
  canvas->SetFillStyle(0);
  canvas->SetFillColor(0);
  canvas->SetGrid(grid);
  canvas->SetLogy(logY);
  canvas->SetLogz(logZ);
}

void setCanvasOptions(TCanvas* canvas, const Int_t grid, const Int_t logY, const Int_t logZ)
{
  canvas->SetFillStyle(0);
  canvas->SetFillColor(0);
  canvas->SetGrid(grid, grid);
  canvas->SetLogy(logY);
  canvas->SetLogz(logZ);
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

void setHistogramOptions(TH1* histogram, const Color_t color, const Size_t size, 
			 const Style_t style, const Double_t scale, const char* xAxisTitle, 
			 const char* yAxisTitle)
{
  histogram->SetMarkerColor(color);
  histogram->SetMarkerSize(size);
  histogram->SetMarkerStyle(style);
  histogram->SetLineColor(color);
  histogram->SetLineWidth(1);
  histogram->SetFillStyle(0);
  setAxisOptions(histogram->GetXaxis(), 0.05, 0.9, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.05, 1.05, yAxisTitle);
  histogram->Scale(scale);
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
  setAxisOptions(histogram->GetXaxis(), 0.04, 0.9, xAxisTitle);
  setAxisOptions(histogram->GetYaxis(), 0.05, yAxisTitleOffset, yAxisTitle);
  setAxisOptions(histogram->GetZaxis(), 0.05, histogram->GetZaxis()->GetTitleOffset(), "");
  histogram->Scale(scale);
}

void setGraphOptions(TGraph* graph, const Color_t color, const Size_t size, const Style_t style, 
		     const char* xAxisTitle, const char* yAxisTitle)
{
  graph->SetMarkerColor(color);
  graph->SetMarkerSize(size);
  graph->SetMarkerStyle(style);
  graph->SetLineColor(color);
  graph->SetTitle("");
  setAxisOptions(graph->GetXaxis(), 0.05, 0.9, xAxisTitle);
  setAxisOptions(graph->GetYaxis(), 0.02, 1.05, yAxisTitle);
}

void setLegendOptions(TLegend* legend, const char* header)
{
  legend->SetFillColor(0);
  legend->SetTextFont(42);
  legend->SetHeader(header);
}

float normAndErrorSquared(const TH3F* ggMETVsDiEMETVsInvMass, const TH1F* controlMET, 
			  const TH1F* egMET, const unsigned int maxNormBin, float& errSquared)
{
  const float nGGMinusEW = 
    ggMETVsDiEMETVsInvMass->Integral(0, -1, 0, -1, 1, maxNormBin) - egMET->Integral(1, maxNormBin);
  const float nControl = controlMET->Integral(1, maxNormBin);
  float norm = 0.0;
  if (nControl != 0.0) norm = nGGMinusEW/nControl;
  float nGGErrSquared = 0.0;
  float nEGErrSquared = 0.0;
  TH1D* ggMET = ggMETVsDiEMETVsInvMass->ProjectionZ("ggMET", 0, -1, 0, -1, "e");
  for (unsigned int iBin = 1; iBin <= maxNormBin; ++iBin) {
    const float ggErr = ggMET->GetBinError(iBin);
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

void setMETErrorBars(TH1F* controlFinal, const vector<TH1F*>& controlMETToyDistsByBin, 
		     const vector<TH1F*>& controlLowSidebandMETToyDistsByBin, 
		     const vector<TH1F*>& controlHighSidebandMETToyDistsByBin, 
		     const float controlNorm, const float controlNormErrSquared, 
		     const bool doDefault)
{
  for (int iBin = 1; iBin <= controlFinal->GetNbinsX(); ++iBin) {
    const float reweightingErr = controlMETToyDistsByBin[iBin - 1]->GetRMS();
    float reweightingErrSquared = reweightingErr*reweightingErr;
    if ((string(controlFinal->GetName()) == "eeFinal") && doDefault) {
      const float lowSidebandReweightingErr = 
	controlLowSidebandMETToyDistsByBin[iBin - 1]->GetRMS();
      const float highSidebandReweightingErr = 
	controlHighSidebandMETToyDistsByBin[iBin - 1]->GetRMS();
      reweightingErrSquared+=4*(lowSidebandReweightingErr*lowSidebandReweightingErr + 
				highSidebandReweightingErr*highSidebandReweightingErr);
    }
    controlFinal->
      SetBinError(iBin, 
		  //stat--this might be wrong (
		  sqrt(controlFinal->GetBinError(iBin)*controlFinal->GetBinError(iBin) + 
		       controlNormErrSquared/*norm*/ + 
		       controlNorm*controlNorm*reweightingErrSquared/*reweighting*/));
  }
}

void compareEEToFF(const vector<string>& inFile, const char* outFile)
{
  TFile* in[2] = {NULL, NULL};
  TCanvas* METCanvas[2] = {NULL, NULL};
  TH1F* egMET[2] = {NULL, NULL};
  TH1F* eeOverFF[2] = {NULL, NULL};
  const Style_t style[2] = {22, 21};
  const Color_t color[2] = {kRed, kBlue};
  const Option_t* drawOption[2] = {"", "SAME"};
  // const string entry[2] = {"i-EM p_{T} + N_{j} reweighting", "i-EM p_{T} reweighting only"};
  // const string entry[2] = {"PF p_{T} reweighting", "ECAL-only p_{T} reweighting"};
  // const string entry[2] = {"ull reweighting", "o reweighting"};
  // const string entry[2] = {"1 GeV #leq m_{ee} < 101 GeV", "o m_{ee} cut"};
  // const string entry[2] = {"o m_{ee} cut", "1 GeV #leq m_{ee} < 101 GeV"};
  // const string entry[2] = {"ncorrected ME_{T}", "ype-I corrected ME_{T}"};
  const string entry[2] = {"F electron and muon cleaning", "/g/f cleaning"};
  // const string entry1stLetterLowercase[2] = {"f", "n"};
  // const string entry1stLetterUppercase[2] = {"F", "N"};
  // const string entry1stLetterLowercase[2] = {"8", "n"};
  // const string entry1stLetterUppercase[2] = {"8", "N"};
  // const string entry1stLetterLowercase[2] = {"n", "8"};
  // const string entry1stLetterUppercase[2] = {"N", "8"};
  // const string entry1stLetterLowercase[2] = {"u", "T"};
  // const string entry1stLetterUppercase[2] = {"U", "T"};
  const string entry1stLetterLowercase[2] = {"P", "e"};
  const string entry1stLetterUppercase[2] = {"P", "e"};
  TCanvas* ratioCanvas = new TCanvas("ratio", "", 600, 300);
  setCanvasOptions(ratioCanvas, 1, 0, 0);
  // TLegend* ratioLegend = new TLegend(0.2, 0.6, 0.6, 0.8);
  TLegend* ratioLegend = new TLegend(0.4, 0.65, 0.8, 0.85);
  setLegendOptions(ratioLegend, "CMS 4.7 fb^{-1}, no jet requirement");
  TFile out(outFile, "RECREATE");
  // const string sample[2] = {"ee", "ff"};
  const string sample[2] = {"gg", "ff"};
  TCanvas* canvas[2] = {NULL, NULL};
  TLegend* legend[2] = {NULL, NULL};
  TH1F* controlEGMET[2][2] = {{NULL, NULL}, {NULL, NULL}};
  TH1F* controlMET[2][2] = {{NULL, NULL}, {NULL, NULL}};
  for (unsigned int i = 0; i < 2; ++i) {
    in[i] = new TFile(inFile[i].c_str());
    in[i]->GetObject("METCanvas", METCanvas[i]);
    egMET[i] = (TH1F*)METCanvas[i]->GetPrimitive("egMET");
    for (unsigned int j = 0; j < 2; ++j) {
      if (sample[j] != "gg") in[i]->GetObject((sample[j] + "Final").c_str(), controlEGMET[i][j]);
      else {
	TH3F* ggEGMETVsDiEMETVsInvMass;
	in[i]->GetObject((sample[j] + "METVsDiEMETVsInvMassTot").c_str(), ggEGMETVsDiEMETVsInvMass);
	controlEGMET[i][j] = 
	  (TH1F*)ggEGMETVsDiEMETVsInvMass->ProjectionZ("ggMET", 0, -1, 0, -1, "e");
      }
      controlMET[i][j] = (TH1F*)controlEGMET[i][j]->Clone();
      controlMET[i][j]->Add(egMET[i], -1.0);
      out.cd();
      if (i == 0) {
	canvas[j] = new TCanvas((sample[j] + "Comparison").c_str(), "", 600, 600);
	setCanvasOptions(canvas[j], 1, 1, 0);
	legend[j] = new TLegend(0.5, 0.6, 0.9, 0.8);
	setLegendOptions(legend[j], "CMS 4.7 fb^{-1}, no jet requirement");
      }
      canvas[j]->cd();
      setHistogramOptions(dynamic_cast<TH1*>(controlMET[i][j]), color[i], 0.7, style[i], 1.0, 
			  "ME_{T} (GeV)", "");
      controlMET[i][j]->Draw(drawOption[i]);
      // legend[j]->AddEntry(controlMET[i][j], (sample[j] + " d" + entry[i]).c_str(), "lp");
      // legend[j]->AddEntry(controlMET[i][j], (sample[j] + " " + entry[i]).c_str(), "lp");
      legend[j]->AddEntry(controlMET[i][j], 
      			  (sample[j] + " " + entry1stLetterLowercase[i] + entry[i]).c_str(), "lp");
      if (i == 1) {
	legend[j]->Draw();
	canvas[j]->Write();
	// canvas[j]->SaveAs(("/Users/rachelyohay/Documents/UVa/dissertation/" + sample[j] + "_dijet_pT_and_Nj_vs_dijet_pT_reweighting.pdf").c_str());
	// canvas[j]->SaveAs(("/Users/rachelyohay/Documents/UVa/dissertation/" + sample[j] + "_di-EM_vs_dijet_pT_reweighting.pdf").c_str());
	// canvas[j]->SaveAs(("/Users/rachelyohay/Documents/UVa/dissertation/" + sample[j] + "_dijet_pT_and_Nj_vs_no_reweighting.pdf").c_str());
      }
    }
    out.cd();
    ratioCanvas->cd();
    eeOverFF[i] = (TH1F*)controlMET[i][0]->Clone();
    eeOverFF[i]->Divide(controlMET[i][1]);
    setHistogramOptions(dynamic_cast<TH1*>(eeOverFF[i]), color[i], 0.7, style[i], 1.0, 
			"ME_{T} (GeV)", "#frac{N_{#gamma#gamma}}{N_{ff}}");
    eeOverFF[i]->Draw(drawOption[i]);
    // ratioLegend->AddEntry(eeOverFF[i], ("D" + entry[i]).c_str(), "lp");
    // ratioLegend->AddEntry(eeOverFF[i], entry[i].c_str(), "lp");
    ratioLegend->AddEntry(eeOverFF[i], (entry1stLetterUppercase[i] + entry[i]).c_str(), "lp");
  }
  ratioLegend->Draw();
  ratioCanvas->Write();
  // ratioCanvas->SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/ee_over_ff_dijet_pT_and_Nj_vs_dijet_pT_reweighting.pdf");
  // ratioCanvas->SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/ee_over_ff_di-EM_vs_dijet_pT_reweighting.pdf");
  // ratioCanvas->SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/ee_over_ff_dijet_pT_and_Nj_vs_no_reweighting.pdf");
  eeOverFF[0]->GetXaxis()->SetRange(1, 10);
  eeOverFF[1]->GetXaxis()->SetRange(1, 10);
  // ratioLegend->SetX1NDC(0.5);
  // ratioLegend->SetX2NDC(0.9);
  // ratioLegend->SetY1NDC(0.6);
  // ratioLegend->SetY2NDC(0.8);
  // ratioLegend->Draw();
  ratioCanvas->Update();
  // ratioCanvas->SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/ee_over_ff_dijet_pT_and_Nj_vs_dijet_pT_reweighting_zoom.pdf");
  // ratioCanvas->SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/ee_over_ff_di-EM_vs_dijet_pT_reweighting_zoom.pdf");
  // ratioCanvas->SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/ee_over_ff_dijet_pT_and_Nj_vs_no_reweighting_zoom.pdf");
  // ratioCanvas->SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/ee_over_ff_Type-I_MET_corrections_vs_uncorrected_MET_zoom.pdf");
  ratioCanvas->SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/gg_over_ff_Ulla_cleaning_vs_Yueh-Feng_cleaning_zoom.pdf");
  for (unsigned int i = 0; i < 2; ++i) { in[i]->Close(); }
  out.Close();
}

void compareEEToFFForEERegions(const char* inFile, const char* outFile)
{
  TFile in(inFile);
  TFile out(outFile, "RECREATE");
  TCanvas ratioCanvas("ratio", "", 600, 300);
  setCanvasOptions(&ratioCanvas, 1, 0, 0);
  TLegend ratioLegend(0.16, 0.17, 0.56, 0.37);
  setLegendOptions(&ratioLegend, "CMS 4.7 fb^{-1}, no jet requirement");
  TH1F* eeOverFF[3] = {NULL, NULL, NULL};
  const Style_t style[3] = {20, 22, 21};
  const Color_t color[3] = {kBlack, kRed, kBlue};
  const Option_t* drawOption[3] = {"", "SAME", "SAME"};
  const string entry[3] = {"81 GeV #leq m_{ee} < 101 GeV", "71 GeV #leq m_{ee} < 81 GeV", 
			   "101 GeV #leq m_{ee} < 111 GeV"};
  const string sample[4] = {"ee", "eeLowSideband", "eeHighSideband", "ff"};
  TH1F* controlEGMET[4] = {NULL, NULL, NULL, NULL};
  TH1F* controlMET[4] = {NULL, NULL, NULL, NULL};
  TCanvas* METCanvas = NULL;
  in.GetObject("METCanvas", METCanvas);
  TH1F* egMET = (TH1F*)METCanvas->GetPrimitive("egMET");
  for (unsigned int i = 0; i < 4; ++i) {
    in.GetObject((sample[i] + "Final").c_str(), controlEGMET[i]);
    controlMET[i] = (TH1F*)controlEGMET[i]->Clone();
    if (sample[i].find("Sideband") == string::npos) controlMET[i]->Add(egMET, -1.0);
    else {
      TH3F* ggMETVsDiEMETVsInvMassTot = NULL;
      in.GetObject("ggMETVsDiEMETVsInvMassTot", ggMETVsDiEMETVsInvMassTot);
      float errSquared = 0.0;
      const float norm = 
	normAndErrorSquared(ggMETVsDiEMETVsInvMassTot, controlMET[i], egMET, 4, errSquared);
      controlMET[i]->Scale(norm);
      vector<TH1F*> controlMETToyDistsByBin(13, NULL);
      for (unsigned int j = 1; j <= 13; ++j) {
	stringstream controlMETToyDistName;
	controlMETToyDistName << sample[i] << "METToyDistBin" << j;
	in.GetObject(controlMETToyDistName.str().c_str(), controlMETToyDistsByBin[j - 1]);
      }
      setMETErrorBars(controlMET[i], controlMETToyDistsByBin, vector<TH1F*>(), vector<TH1F*>(), 
		      norm, errSquared, true);
    }
  }
  out.cd();
  ratioCanvas.cd();
  for (unsigned int i = 0; i < 3; ++i) {
    eeOverFF[i] = (TH1F*)controlMET[i]->Clone();
    eeOverFF[i]->Divide(controlMET[3]);
    setHistogramOptions(dynamic_cast<TH1*>(eeOverFF[i]), color[i], 0.7, style[i], 1.0, 
			"ME_{T} (GeV)", "#frac{N_{ee}}{N_{ff}}");
    eeOverFF[i]->GetXaxis()->SetRange(1, 9);
    eeOverFF[i]->Draw(drawOption[i]);
    ratioLegend.AddEntry(eeOverFF[i], entry[i].c_str(), "lp");
  }
  ratioLegend.Draw();
  ratioCanvas.Write();
  ratioCanvas.SaveAs("/Users/rachelyohay/RA3/data/eeOverFF_ZSignalVsSidebands.pdf");
  in.Close();
  out.Close();
}

void replotEMF(const char* inFile, const char* outFile)
{
  TFile in(inFile);
  TFile out(outFile, "RECREATE");
  TCanvas profileCanvas("profileCanvas", "", 600, 600);
  setCanvasOptions(&profileCanvas, 1, 0, 0);
  TLegend legend(0.25, 0.2, 0.65, 0.4);
  setLegendOptions(&legend, "CMS 4.7 fb^{-1}, no jet requirement");
  const string type[3] = {"f", "g", "e"};
  const string label[3] = {"f", "#gamma", "e"};
  const Color_t color[3] = {kBlue, kBlack, kRed};
  const string order[2] = {"Trailing", "Leading"};
  TH2F* histogram[3][2] = {{NULL, NULL}, {NULL, NULL}, {NULL, NULL}};
  TCanvas* canvas[3][2] = {{NULL, NULL}, {NULL, NULL}, {NULL, NULL}};
  const Style_t style[3][2] = {{25, 21}, {24, 20}, {26, 22}};
  const Option_t* drawOption[3][2] = {{"", "SAME"}, {"SAME", "SAME"}, {"SAME", "SAME"}};
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 2; ++j) {
      in.GetObject((type[i] + order[j] + "ETJetNormVsEMFractionPF").c_str(), histogram[i][j]);
      histogram[i][j]->SetTitle("");
      setHistogramOptions(histogram[i][j], color[i], 0.5, style[i][j], 1.6, 1.0, "EMF", 
			  ("(p_{Tj} - p_{T" + label[i] + "})/p_{Tj}").c_str());
      out.cd();
      canvas[i][j] = 
	new TCanvas(("4684pb-1_" + type[i] + order[j] + "_ETBias_log").c_str(), "", 600, 600);
      setCanvasOptions(canvas[i][j], 0, 0, 1);
      canvas[i][j]->cd()->SetLeftMargin(0.2);
      canvas[i][j]->cd()->SetTopMargin(0.2);
      canvas[i][j]->cd()->SetRightMargin(0.2);
      canvas[i][j]->cd()->SetBottomMargin(0.2);
      histogram[i][j]->Draw("COLZ");
      canvas[i][j]->Write();
      canvas[i][j]->SaveAs(("/Users/rachelyohay/Documents/UVa/dissertation/4684pb-1_" + type[i] + 
			    order[j] + "_ETBias_log.pdf").c_str());
      profileCanvas.cd();
      TProfile* profile = histogram[i][j]->ProfileX();
      setHistogramOptions(dynamic_cast<TH1*>(profile), color[i], 0.5, style[i][j], 1.0, "EMF", 
			  "<(p_{Tj} - p_{TEM})/p_{Tj}>");
      legend.AddEntry(profile, (order[j] + " " + label[i]).c_str(), "lp");
      profile->Draw(drawOption[i][j]);
    }
  }
  legend.Draw();
  profileCanvas.Write();
  profileCanvas.SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/4684pb-1_single_ETBias.pdf");
  out.Close();
  in.Close();
}

void replotMET(const char* inFile, const char* outFile)
{
  TFile in(inFile);
  TFile out(outFile, "RECREATE");
  TCanvas METCanvas("METCanvas", "", 600, 600);
  setCanvasOptions(&METCanvas, 1, 1, 0);
  TCanvas egMETCanvas("egMETCanvas", "", 600, 600);
  setCanvasOptions(&egMETCanvas, 1, 1, 0);
  TLegend legend(0.6, 0.7, 0.9, 0.9);
  setLegendOptions(&legend, "CMS 4.7 fb^{-1}, no jet requirement");
  TLegend egLegend(0.6, 0.7, 0.9, 0.9);
  setLegendOptions(&egLegend, "CMS 4.7 fb^{-1}, no jet requirement");
  const string region[4] = {"ee", "eeLowSideband", "eeHighSideband", "ff"};
  const string sample[3] = {"gg", "ee", "ff"};
  const unsigned int jetBin[3] = {1, 2, 3};
  const Color_t color[3] = {kBlack, kRed, kBlue};
  const Style_t style[3] = {20, 22, 21};
  const Option_t* drawOption[3] = {"", "SAME", "SAME"};
  const char* entry[3] = {"81 GeV #leq m_{ee} < 101 GeV", "71 GeV #leq m_{ee} < 81 GeV", 
			  "101 GeV #leq m_{ee} < 111 GeV"};
  const char* dijetPTEntry[3] = {"#gamma#gamma", "ee", "ff"};
  const char* header[3] = {"CMS 4.7 fb^{-1}, 0 jets", "CMS 4.7 fb^{-1}, 1 jet", 
			   "CMS 4.7 fb^{-1}, #geq2 jets"};
  TH1F* METHistogram[3] = {NULL, NULL, NULL};
  TCanvas* dijetPTCanvas[3] = {NULL, NULL, NULL};
  TLegend* dijetPTLegend[3] = {NULL, NULL, NULL};
  TH1F* dijetPTHistogram[3][3] = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
  TH1F* histogram[4][3] = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, 
			   {NULL, NULL, NULL}, {NULL, NULL, NULL}};
  TCanvas* canvas[4][3] = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, 
			   {NULL, NULL, NULL}, {NULL, NULL, NULL}};
  for (unsigned int i = 0; i < 4; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      stringstream histogramName;
      histogramName << region[i] << "Weights_METBin1_NjBin" << jetBin[j];
      in.GetObject(histogramName.str().c_str(), histogram[i][j]);
      histogram[i][j]->SetTitle("");
      setHistogramOptions(histogram[i][j], kBlack, 0.5, 20, 1.0, "Di-EM p_{T} (GeV)", "Weight");
      if (i < 3) {
	stringstream dijetPTHistogramName;
	dijetPTHistogramName << sample[i] << "DiEMETUniformNjBin" << jetBin[j];
	in.GetObject(dijetPTHistogramName.str().c_str(), dijetPTHistogram[i][j]);
	dijetPTHistogram[i][j]->SetTitle("");
	setHistogramOptions(dijetPTHistogram[i][j], color[i], 0.5, style[i], 
			    dijetPTHistogram[0][j]->Integral()/dijetPTHistogram[i][j]->Integral(), 
			    "Di-EM p_{T} (GeV)", "");
	dijetPTHistogram[i][j]->GetXaxis()->SetRange(1, 200);
      }
      stringstream canvasName;
      canvasName << region[i] << "_" << (jetBin[j] - 1) << "_jet_dijet_pT_weights";
      stringstream dijetPTCanvasName;
      dijetPTCanvasName << (jetBin[j] - 1) << "_jet_dijet_pT";
      out.cd();
      canvas[i][j] = new TCanvas(canvasName.str().c_str(), "", 600, 600);
      setCanvasOptions(canvas[i][j], 1, 0, 0);
      histogram[i][j]->Draw();
      canvas[i][j]->Write();
      // canvas[i][j]->SaveAs(("/Users/rachelyohay/Documents/UVa/dissertation/" + canvasName.str() + 
      // 			    ".pdf").c_str());
      if (i == 0) {
	dijetPTCanvas[j] = new TCanvas(dijetPTCanvasName.str().c_str(), "", 600, 600);
	setCanvasOptions(dijetPTCanvas[j], 1, 1, 0);
	dijetPTLegend[j] = new TLegend(0.6, 0.7, 0.9, 0.9);
	setLegendOptions(dijetPTLegend[j], header[j]);
      }
      dijetPTCanvas[j]->cd();
      if (i < 3) {
	dijetPTHistogram[i][j]->Draw(drawOption[i]);
	dijetPTLegend[j]->AddEntry(dijetPTHistogram[i][j], dijetPTEntry[i], "lp");
      }
      if (i == 2) {
	dijetPTLegend[j]->Draw();
	dijetPTCanvas[j]->Write();
	// dijetPTCanvas[j]->SaveAs(("/Users/rachelyohay/Documents/UVa/dissertation/" + 
	// 			  dijetPTCanvasName.str() + ".pdf").c_str());
      }
    }
    if (region[i] == "ee") {
      TH1F* eePlusEGMET = NULL;
      in.GetObject((region[i] + "Final").c_str(), eePlusEGMET);
      TCanvas* intermediateCanvas = NULL;
      in.GetObject("METCanvas", intermediateCanvas);
      TH1F* egMET = (TH1F*)intermediateCanvas->GetPrimitive("egMET");
      setHistogramOptions(egMET, kBlack, 0.5, 20, 1.0, "ME_{T} (GeV)", "");
      egMETCanvas.cd();
      egMET->Draw();
      egLegend.Draw();
      egMETCanvas.Write();
      egMETCanvas.SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/eg_MET.pdf");
      METHistogram[i] = (TH1F*)eePlusEGMET->Clone();
      METHistogram[i]->Add(egMET, -1.0);
    }
    else in.GetObject((region[i] + "Final").c_str(), METHistogram[i]);
    if (i < 3) {
      setHistogramOptions(METHistogram[i], color[i], 0.5, style[i], 
			  METHistogram[0]->Integral()/METHistogram[i]->Integral(), "ME_{T} (GeV)", 
			  "");
      legend.AddEntry(METHistogram[i], entry[i], "lp");
    }
    out.cd();
    METCanvas.cd();
    if (i < 3) METHistogram[i]->Draw(drawOption[i]);
  }
  legend.Draw();
  METCanvas.Write();
  // METCanvas.SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/all_ee_MET_spectra.pdf");
  out.Close();
  in.Close();
}

vector<float> estimateMEEBkgShapeErr()
{
  // TFile nomSBWeights("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_METBinsForLimits.root");
  // TFile lowSBUnityWeight("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_eeLowSidebandOnly_METBinsForLimits.root");
  // TFile highSBUnityWeight("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_eeHighSidebandOnly_METBinsForLimits.root");
  TFile nomSBWeights("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_correctErrors.root");
  TFile lowSBUnityWeight("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_eeLowSidebandOnly.root");
  TFile highSBUnityWeight("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_eeHighSidebandOnly.root");

  TFile out("/Users/rachelyohay/RA3/data/debug.root", "RECREATE");

  TH1F* eeNom;
  TH1F* eeLowSBUnityWeight;
  TH1F* eeHighSBUnityWeight;
  nomSBWeights.GetObject("eeFinal", eeNom);
  lowSBUnityWeight.GetObject("eeFinal", eeLowSBUnityWeight);
  highSBUnityWeight.GetObject("eeFinal", eeHighSBUnityWeight);
  setHistogramOptions(dynamic_cast<TH1*>(eeNom), kBlack, 0.5, 20, 1.0, "ME_{T} (GeV)", "");
  setHistogramOptions(dynamic_cast<TH1*>(eeLowSBUnityWeight), kRed, 0.5, 21, 1.0, "ME_{T} (GeV)", 
		      "");
  setHistogramOptions(dynamic_cast<TH1*>(eeHighSBUnityWeight), kBlue, 0.5, 22, 1.0, "ME_{T} (GeV)", 
		      "");

  out.cd();
  TCanvas MET("MET", "", 400, 600);
  MET.Divide(1, 2);
  MET.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(MET.cd(1), 1, 1, 0);
  MET.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(MET.cd(2), 1, 0, 0);

  MET.cd(1);
  eeNom->Draw();
  eeLowSBUnityWeight->Draw("SAME");
  eeHighSBUnityWeight->Draw("SAME");
  TLegend legend(0.5, 0.6, 0.9, 0.8);
  setLegendOptions(&legend, "CMS 4.7 fb^{-1}, no jet requirement");
  legend.AddEntry(eeNom, "Nominal ee sideband weights", "lp");
  legend.AddEntry(eeLowSBUnityWeight, "ee low sideband only", "lp");
  legend.AddEntry(eeHighSBUnityWeight, "ee high sideband only", "lp");
  legend.Draw();

  TH1F* lowSBOverNom = (TH1F*)eeLowSBUnityWeight->Clone();
  lowSBOverNom->Divide(eeNom);
  setHistogramOptions(dynamic_cast<TH1*>(lowSBOverNom), kRed, 0.5, 21, 1.0, "ME_{T} (GeV)", 
		      "N_{SB}/N_{nom}");
  TH1F* highSBOverNom = (TH1F*)eeHighSBUnityWeight->Clone();
  highSBOverNom->Divide(eeNom);
  setHistogramOptions(dynamic_cast<TH1*>(highSBOverNom), kBlue, 0.5, 22, 1.0, "ME_{T} (GeV)", 
		      "N_{SB}/N_{nom}");
  MET.cd(2);
  lowSBOverNom->Draw();
  highSBOverNom->Draw("SAME");

  vector<float> meeBkgShapeError;
  for (Int_t iBin = 1; iBin <= eeNom->GetNbinsX(); ++iBin) {
    const float nEENom = eeNom->GetBinContent(iBin);
    const float absNomMinusLowSB = fabs(nEENom - eeLowSBUnityWeight->GetBinContent(iBin));
    const float absNomMinusHighSB = fabs(nEENom - eeHighSBUnityWeight->GetBinContent(iBin));
    cout << "Bin " << iBin << ", ";
    if (absNomMinusLowSB > absNomMinusHighSB) {
      cout << absNomMinusLowSB << endl;
      meeBkgShapeError.push_back(absNomMinusLowSB);
    }
    else {
      cout << absNomMinusHighSB << endl;
      meeBkgShapeError.push_back(absNomMinusHighSB);
    }
  }

  MET.Write();
  // MET.SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/mee_bkg_shape_error.pdf");
  out.Close();
  nomSBWeights.Close();
  lowSBUnityWeight.Close();
  highSBUnityWeight.Close();
  return meeBkgShapeError;
}

void makeFinalMETPlots(const char* inFile, const char* outFile)
{
  //mee bkg. shape error (largest of |nom. ee - low SB ee| and |nom. ee - high SB ee| per bin)
  vector<float> meeBkgShapeError = estimateMEEBkgShapeErr();

  //get distributions
  TFile in(inFile);
  if (!in.IsOpen()) {
    cerr << "Error opening file " << inFile << ".\n";
    return;
  }
  TH1F* eeEG = NULL;
  TH1F* ffEG = NULL;
  TCanvas* canvas = NULL;
  in.GetObject("eeFinal", eeEG);
  in.GetObject("ffFinal", ffEG);
  in.GetObject("METCanvas", canvas);
  if ((eeEG == NULL) || (ffEG == NULL) || (canvas == NULL)) {
    cerr << "Error getting histograms eeFinal or ffFinal or canvas METCanvas.\n";
    in.Close();
    return;
  }
  TH1F* eg = (TH1F*)canvas->GetPrimitive("egMET");
  TH1F* gg = (TH1F*)canvas->GetPrimitive("ggMET");
  if ((eg == NULL) || (gg == NULL)) {
    cerr << "Error getting histogram egMET or ggMET.\n";
    in.Close();
    return;
  }

  //make distributions pretty
  setHistogramOptions(dynamic_cast<TH1*>(eeEG), /*kBlue*/857, 0.0, 20, 1.0, "ME_{T} (GeV)", "");
  setHistogramOptions(dynamic_cast<TH1*>(ffEG), /*kMagenta*/619, 0.0, 20, 1.0, "ME_{T} (GeV)", "");
  setHistogramOptions(dynamic_cast<TH1*>(eg), /*kGreen*/834, 0.0, 20, 1.0, "ME_{T} (GeV)", "");
  setHistogramOptions(dynamic_cast<TH1*>(gg), kBlack, 0.7, 20, 1.0, "ME_{T} (GeV)", "");
  eeEG->SetFillStyle(3017);
  ffEG->SetFillStyle(3018);
  eg->SetFillStyle(3017);
  eeEG->SetFillColor(/*kBlue*/857);
  ffEG->SetFillColor(/*kMagenta*/619);
  eg->SetFillColor(/*kGreen*/834);
  eeEG->GetYaxis()->SetRangeUser(4e-1, 40000.0);
  ffEG->GetYaxis()->SetRangeUser(4e-1, 40000.0);
  eg->GetYaxis()->SetRangeUser(4e-1, 40000.0);
  gg->GetYaxis()->SetRangeUser(4e-1, 40000.0);

  //open output file and make canvas
  TFile out(outFile, "RECREATE");
  if (!out.IsOpen()) {
    cerr << "Error opening file " << outFile << ".\n";
    in.Close();
    return;
  }
  TCanvas allMET("allMET", "", 400, 600);
  allMET.Divide(1, 2);
  allMET.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(allMET.cd(1), 1, 1, 0);
  allMET.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(allMET.cd(2), 1, 0, 0);
  TLegend allMETLegend(0.35, 0.5, 0.9, 0.9);
  TLegend allMETRatioLegend(0.35, 0.7, 0.9, 0.9);
  // setLegendOptions(&allMETLegend, "CMS 4.7 fb^{-1}, no jet requirement");
  // setLegendOptions(&allMETRatioLegend, "CMS 4.7 fb^{-1}, no jet requirement");
  setLegendOptions(&allMETLegend, "CMS 4.7 fb^{-1}, #geq1 jet requirement");
  setLegendOptions(&allMETRatioLegend, "CMS 4.7 fb^{-1}, #geq 1 jet requirement");
  TCanvas ffMET("ffMET", "", 400, 600);
  ffMET.Divide(1, 2);
  ffMET.cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
  setCanvasOptions(ffMET.cd(1), 1, 1, 0);
  ffMET.cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
  setCanvasOptions(ffMET.cd(2), 1, 0, 0);
  TLegend ffMETLegend(0.35, 0.5, 0.9, 0.9);
  TLegend ffMETRatioLegend(0.35, 0.7, 0.9, 0.9);
  // setLegendOptions(&ffMETLegend, "CMS 4.7 fb^{-1}, no jet requirement");
  // setLegendOptions(&ffMETRatioLegend, "CMS 4.7 fb^{-1}, no jet requirement");
  setLegendOptions(&ffMETLegend, "CMS 4.7 fb^{-1}, #geq 1 jet requirement");
  setLegendOptions(&ffMETRatioLegend, "CMS 4.7 fb^{-1}, #geq 1 jet requirement");

  //add mee bkg. shape error in quadrature to ee MET error
  const Int_t nEEEGBins = eeEG->GetNbinsX();
  const Int_t nEEBkgShapeErrs = (Int_t)meeBkgShapeError.size();
  if (nEEEGBins != nEEBkgShapeErrs) {
    cerr << "Error: " << nEEEGBins << " bins in ee + eg MET histogram, " << nEEBkgShapeErrs;
    cerr << " ee bkg. shape errors.\n";
    out.Close();
    in.Close();
    return;
  }
  for (Int_t iBin = 1; iBin <= nEEEGBins; ++iBin) {
    Double_t binErr = eeEG->GetBinError(iBin);
    eeEG->SetBinError(iBin, 
    		      sqrt(binErr*binErr + meeBkgShapeError[iBin - 1]*meeBkgShapeError[iBin - 1]));
  }

  //draw MET distributions
  allMET.cd(1);
  eg->Draw("E2");
  ffEG->Draw("E2SAME");
  eeEG->Draw("E2SAME");
  gg->Draw("SAME");

  //add MET legend entries
  allMETLegend.AddEntry(eg, "EW background (from e#gamma)", "f");
  allMETLegend.AddEntry(ffEG, "Total (QCD + EW) background (from ff + e#gamma)", "f");
  allMETLegend.AddEntry(eeEG, "Total (QCD + EW) background (from ee + e#gamma)", "f");
  allMETLegend.AddEntry(gg, "Data (#gamma#gamma)", "lp");
  allMETLegend.Draw();

  //draw MET ratio distributions
  TH1F* ggOverEEEG = (TH1F*)gg->Clone();
  ggOverEEEG->Divide(eeEG);
  setHistogramOptions(dynamic_cast<TH1*>(ggOverEEEG), /*kBlue*/857, 0.5, 20, 1.0, "ME_{T} (GeV)", 
		      "N_{#gamma#gamma}/N_{QCD+EW}");
  ggOverEEEG->GetYaxis()->SetRangeUser(0.0, 2.0);
  for (Int_t iBin = 1; iBin <= nEEEGBins; ++iBin) {
    Double_t ggBinErr = gg->GetBinError(iBin);
    Double_t eeEGBinErr = eeEG->GetBinError(iBin);
    Double_t ggBinContent = gg->GetBinContent(iBin);
    Double_t eeEGBinContent = eeEG->GetBinContent(iBin);
    ggOverEEEG->
      SetBinError(iBin, ggOverEEEG->GetBinContent(iBin)*sqrt(((ggBinErr*ggBinErr)/
							      (ggBinContent*ggBinContent)) + 
							     ((eeEGBinErr*eeEGBinErr)/
							      (eeEGBinContent*eeEGBinContent))));
  }
  TH1F* ggOverFFEG = (TH1F*)gg->Clone();
  ggOverFFEG->Divide(ffEG);
  setHistogramOptions(dynamic_cast<TH1*>(ggOverFFEG), /*kMagenta*/619, 0.5, 20, 1.0, 
		      "ME_{T} (GeV)", "N_{#gamma#gamma}/N_{QCD+EW}");
  ggOverFFEG->GetYaxis()->SetRangeUser(0.0, 2.0);
  allMET.cd(2);
  ggOverEEEG->Draw();
  ggOverFFEG->Draw("SAME");

  //add MET ratio legend entries
  allMETRatioLegend.AddEntry(ggOverEEEG, "N_{#gamma#gamma}/N_{ee + e#gamma}", "lp");
  allMETRatioLegend.AddEntry(ggOverFFEG, "N_{#gamma#gamma}/N_{ff + e#gamma}", "lp");
  allMETRatioLegend.Draw();

  //write and save
  allMET.Write();
  // allMET.SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/MET_final_ee_ff_eg.pdf");
  allMET.SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/MET_final_ee_ff_eg_geq1.pdf");

  //clone ff + eg MET histogram and add |ff - ee| error in quadrature
  TH1F* ffEGTotErr = (TH1F*)ffEG->Clone();
  const Int_t nFFBins = ffEGTotErr->GetNbinsX();
  for (Int_t iBin = 1; iBin <= nFFBins; ++iBin) {
    Double_t binErr = ffEGTotErr->GetBinError(iBin);
    Double_t absFFMinusEE = fabs(ffEGTotErr->GetBinContent(iBin) - eeEG->GetBinContent(iBin));
    ffEGTotErr->SetBinError(iBin, sqrt(binErr*binErr + absFFMinusEE*absFFMinusEE));
  }

  //print ff results for sensitive MET bins
  cout << "ff:\n";
  for (Int_t iBin = 10; iBin <= nFFBins; ++iBin) {
    cout << iBin << " " << ffEGTotErr->GetBinContent(iBin) << " " << ffEGTotErr->GetBinError(iBin);
    cout << endl;
  }

  cout << "gg:\n";
  //print gg results for sensitive MET bins
  for (Int_t iBin = 10; iBin <= nFFBins; ++iBin) {
    cout << iBin << " " << gg->GetBinContent(iBin) << endl;
  }

  //draw ff MET distribution
  ffMET.cd(1);
  eg->Draw("E2");
  ffEGTotErr->Draw("E2SAME");
  gg->Draw("SAME");

  //add ff MET legend entries
  ffMETLegend.AddEntry(eg, "EW background (from e#gamma)", "f");
  ffMETLegend.AddEntry(ffEGTotErr, "Total (QCD + EW) background (from ff + e#gamma)", "f");
  ffMETLegend.AddEntry(gg, "Data (#gamma#gamma)", "lp");
  ffMETLegend.Draw();

  //draw ff MET ratio distribution
  TH1F* ggOverFFEGTotErr = (TH1F*)gg->Clone();
  ggOverFFEGTotErr->Divide(ffEGTotErr);
  setHistogramOptions(dynamic_cast<TH1*>(ggOverFFEGTotErr), /*kMagenta*/619, 0.5, 20, 1.0, 
		      "ME_{T} (GeV)", "N_{#gamma#gamma}/N_{QCD+EW}");
  ggOverFFEGTotErr->GetYaxis()->SetRangeUser(0.0, 2.0);
  for (Int_t iBin = 1; iBin <= nFFBins; ++iBin) {
    Double_t ggBinErr = gg->GetBinError(iBin);
    Double_t ffEGBinErr = ffEGTotErr->GetBinError(iBin);
    Double_t ggBinContent = gg->GetBinContent(iBin);
    Double_t ffEGBinContent = ffEGTotErr->GetBinContent(iBin);
    ggOverFFEG->
      SetBinError(iBin, ggOverFFEGTotErr->
		  GetBinContent(iBin)*sqrt(((ggBinErr*ggBinErr)/
					    (ggBinContent*ggBinContent)) + 
					   ((ffEGBinErr*ffEGBinErr)/
					    (ffEGBinContent*ffEGBinContent))));
  }
  ffMET.cd(2);
  ggOverFFEGTotErr->Draw();

  //add ff MET ratio legend entry
  ffMETRatioLegend.AddEntry(ggOverFFEGTotErr, "N_{#gamma#gamma}/N_{ff + e#gamma}", "lp");
  // ffMETRatioLegend.Draw();

  //write and save
  ffMET.Write();
  // ffMET.SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/MET_final_ff_eg.pdf");
  ffMET.SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/MET_final_ff_eg_geq1.pdf");

  //close
  out.Write();
  out.Close();
  in.Close();
}

void compareDiEMPTBinSizes(const vector<string>& inputFiles, const vector<Color_t>& color, 
			   const vector<Style_t>& markerStyle, const vector<string>& entryLabel, 
			   const vector<Option_t*>& drawOption, const char* outputFile)
{
  //open output file
  TFile output(outputFile, "RECREATE");

  //instantiate ee and ff MET canvases to show the effects of different di-EM pT binning
  TCanvas eeMETCanvas("eeMETCanvas", "", 600, 600);
  setCanvasOptions(&eeMETCanvas, 1, 0, 0);
  TCanvas ffMETCanvas("ffMETCanvas", "", 600, 600);
  setCanvasOptions(&ffMETCanvas, 1, 0, 0);

  //instantiate legends
  TLegend eeMETLegend(0.6, 0.7, 0.9, 0.9);
  setLegendOptions(&eeMETLegend, "CMS 4.7 fb^{-1}, no jet requirement");
  TLegend ffMETLegend(0.6, 0.7, 0.9, 0.9);
  setLegendOptions(&ffMETLegend, "CMS 4.7 fb^{-1}, no jet requirement");

  //allocate all the needed pointers
  vector<TFile*> input(inputFiles.size(), NULL);
  vector<TH1F*> eeEGMET(inputFiles.size(), NULL);
  vector<TH1F*> ffEGMET(inputFiles.size(), NULL);
  vector<TH1F*> egMET(inputFiles.size(), NULL);
  vector<TH1F*> eeMET(inputFiles.size(), NULL);
  vector<TH1F*> ffMET(inputFiles.size(), NULL);
  vector<TCanvas*> METCanvas(inputFiles.size(), NULL);

  //loop over the input files, 1 for each bin definition
  for (vector<string>::const_iterator iInputFile = inputFiles.begin(); 
       iInputFile != inputFiles.end(); ++iInputFile) {
    const unsigned int i = iInputFile - inputFiles.begin();

    //get eeFinal and ffFinal
    input[i] = new TFile(iInputFile->c_str());
    input[i]->GetObject("eeFinal", eeEGMET[i]);
    input[i]->GetObject("ffFinal", ffEGMET[i]);

    //get MET canvas and egMET
    input[i]->GetObject("METCanvas", METCanvas[i]);
    egMET[i] = (TH1F*)METCanvas[i]->GetPrimitive("egMET");

    //subtract egMET from eeFinal and ffFinal
    eeMET[i] = (TH1F*)eeEGMET[i]->Clone();
    eeMET[i]->Add(egMET[i], -1.0);
    ffMET[i] = (TH1F*)ffEGMET[i]->Clone();
    ffMET[i]->Add(egMET[i], -1.0);

    //set histogram options
    setHistogramOptions(dynamic_cast<TH1*>(eeMET[i]), color[i], 0.5, markerStyle[i], 
			1.0, "ME_{T} (GeV)", "");
    setHistogramOptions(dynamic_cast<TH1*>(ffMET[i]), color[i], 0.5, markerStyle[i], 
			1.0, "ME_{T} (GeV)", "");

    //add legend entries
    eeMETLegend.AddEntry(eeMET[i], entryLabel[i].c_str(), "lp");
    ffMETLegend.AddEntry(ffMET[i], entryLabel[i].c_str(), "lp");

    //plot on respective canvases
    output.cd();
    eeMETCanvas.cd();
    eeMET[i]->Draw(drawOption[i]);
    ffMETCanvas.cd();
    ffMET[i]->Draw(drawOption[i]);
  }

  //close input and output files
  eeMETCanvas.cd();
  eeMETLegend.Draw();
  eeMETCanvas.Write();
  ffMETCanvas.cd();
  ffMETLegend.Draw();
  ffMETCanvas.Write();
  output.Close();
  for (vector<TFile*>::iterator iInput = input.begin(); iInput != input.end(); ++iInput) {
    (*iInput)->Close();
  }
}

void readErrorFile(ifstream& inTxt, const unsigned int nBins, vector<float>& err)
{
  string label;
  while (isalpha(inTxt.peek()) || (inTxt.peek() == '+')) {
    inTxt >> label;
    inTxt.seekg(1, ios::cur);
  }
  for (unsigned int iBin = 0; iBin < nBins; ++iBin) {
    inTxt >> err[iBin];
  }
  inTxt.seekg(2, ios::cur);
}

void printTableLine(const char* label, const vector<unsigned int>& binNums, 
		    const vector<float>& err, const vector<float>& bkg)
{
  cout << label;
  for (vector<unsigned int>::const_iterator iBinNum = binNums.begin(); iBinNum != binNums.end(); 
       ++iBinNum) {
    cout.precision(2);
    cout << " & " << ((err[*iBinNum - 1]/bkg[*iBinNum - 1])*100);
  }
  cout << " \\\\\n";
}

void printErrors(const char* inFileROOT, const char* inFileTxt, 
		 const vector<unsigned int>& binNums)
{
  //vectors for the background estimates
  vector<float> eeEGBkg;
  vector<float> ffEGBkg;
  vector<float> eeBkg;
  vector<float> ffBkg;
  vector<float> egBkg;

  //get background estimates
  TFile inROOT(inFileROOT);
  TH1F* eeEGMET = NULL;
  TH1F* ffEGMET = NULL;
  TCanvas* MET = NULL;
  inROOT.GetObject("ffFinal", ffEGMET);
  inROOT.GetObject("eeFinal", eeEGMET);
  inROOT.GetObject("METCanvas", MET);
  const unsigned int nBins = eeEGMET->GetNbinsX();
  TH1D* egMET = (TH1D*)MET->GetPrimitive("egMET");
  for (unsigned int iBin = 1; iBin <= nBins; ++iBin) {
    const float nEG = egMET->GetBinContent(iBin);
    const float nEEEG = eeEGMET->GetBinContent(iBin);
    const float nFFEG = ffEGMET->GetBinContent(iBin);
    eeEGBkg.push_back(nEEEG);
    ffEGBkg.push_back(nFFEG);
    egBkg.push_back(nEG);
    eeBkg.push_back(nEEEG - nEG);
    ffBkg.push_back(nFFEG - nEG);
  }
  inROOT.Close();

  //vectors for the different kinds of errors
  vector<float> ffNormStat(nBins, 0.0);
  vector<float> ffNormFEG(nBins, 0.0);
  vector<float> ffNormReweighting(nBins, 0.0);
  vector<float> ffNormTot(nBins, 0.0);
  vector<float> ffReweighting(nBins, 0.0);
  vector<float> ffReweightingTot(nBins, 0.0);
  vector<float> ffStat(nBins, 0.0);
  vector<float> ffStatTot(nBins, 0.0);
  vector<float> ffStatReweighting(nBins, 0.0);
  vector<float> ffTot(nBins, 0.0);
  vector<float> ffMinusEE(nBins, 0.0);
  vector<float> ffSys(nBins, 0.0);
  vector<float> eeNormStat(nBins, 0.0);
  vector<float> eeNormFEG(nBins, 0.0);
  vector<float> eeNormReweighting(nBins, 0.0);
  vector<float> eeNormTot(nBins, 0.0);
  vector<float> eeReweighting(nBins, 0.0);
  vector<float> eeReweightingTot(nBins, 0.0);
  vector<float> eeStat(nBins, 0.0);
  vector<float> eeStatTot(nBins, 0.0);
  vector<float> eeStatReweighting(nBins, 0.0);
  vector<float> eeTot(nBins, 0.0);
  vector<float> eeMEEBkgShape(nBins, 0.0);
  vector<float> eeSys(nBins, 0.0);
  vector<float> egStat(nBins, 0.0);
  vector<float> egFEG(nBins, 0.0);
  vector<float> egTot(nBins, 0.0);
  vector<float> eeEGTot(nBins, 0.0);
  vector<float> eeEGStatTot(nBins, 0.0);
  vector<float> eeEGSysTot(nBins, 0.0);
  vector<float> eeEGStatSys(nBins, 0.0);
  vector<float> ffEGTot(nBins, 0.0);
  vector<float> ffEGStatTot(nBins, 0.0);
  vector<float> ffEGSysTot(nBins, 0.0);
  vector<float> ffEGStatSys(nBins, 0.0);

  //read the error file
  ifstream inTxt(inFileTxt);
  readErrorFile(inTxt, nBins, ffNormStat);
  readErrorFile(inTxt, nBins, ffNormFEG);
  readErrorFile(inTxt, nBins, ffNormReweighting);
  readErrorFile(inTxt, nBins, ffNormTot);
  readErrorFile(inTxt, nBins, ffReweighting);
  readErrorFile(inTxt, nBins, ffStat);
  readErrorFile(inTxt, nBins, ffTot);
  readErrorFile(inTxt, nBins, eeNormStat);
  readErrorFile(inTxt, nBins, eeNormFEG);
  readErrorFile(inTxt, nBins, eeNormReweighting);
  readErrorFile(inTxt, nBins, eeNormTot);
  readErrorFile(inTxt, nBins, eeReweighting);
  readErrorFile(inTxt, nBins, eeStat);
  readErrorFile(inTxt, nBins, eeTot);
  readErrorFile(inTxt, nBins, egStat);
  readErrorFile(inTxt, nBins, egFEG);
  readErrorFile(inTxt, nBins, egTot);
  readErrorFile(inTxt, nBins, eeEGTot);
  readErrorFile(inTxt, nBins, ffEGTot);
  inTxt.close();

  //import mee bkg. shape error
  //calculate total QCD stat. error, reweighting error, and their quadrature sum
  //calculate ff-ee error and total QCD sys. error
  //calculate total QCD stat. + sys. error
  //calculate total QCD + EW stat. error
  //calculate total QCD + EW sys. error
  vector<float> meeBkgShapeErr = estimateMEEBkgShapeErr();
  for (unsigned int iBin = 0; iBin < nBins; ++iBin) {
    eeMEEBkgShape[iBin] = meeBkgShapeErr[iBin];
    ffStatTot[iBin] = sqrt(ffNormStat[iBin]*ffNormStat[iBin] + ffStat[iBin]*ffStat[iBin]);
    eeStatTot[iBin] = sqrt(eeNormStat[iBin]*eeNormStat[iBin] + eeStat[iBin]*eeStat[iBin]);
    ffReweightingTot[iBin] = sqrt(ffNormReweighting[iBin]*ffNormReweighting[iBin] + 
				  ffReweighting[iBin]*ffReweighting[iBin]);
    eeReweightingTot[iBin] = sqrt(eeNormReweighting[iBin]*eeNormReweighting[iBin] + 
				  eeReweighting[iBin]*eeReweighting[iBin]);
    ffStatReweighting[iBin] = sqrt(ffStatTot[iBin]*ffStatTot[iBin] + 
				   ffReweightingTot[iBin]*ffReweightingTot[iBin]);
    eeStatReweighting[iBin] = sqrt(eeStatTot[iBin]*eeStatTot[iBin] + 
				   eeReweightingTot[iBin]*eeReweightingTot[iBin]);
    ffMinusEE[iBin] = fabs(ffEGBkg[iBin] - eeEGBkg[iBin]);
    ffSys[iBin] = sqrt(ffNormFEG[iBin]*ffNormFEG[iBin] + ffMinusEE[iBin]*ffMinusEE[iBin]);
    eeSys[iBin] = sqrt(eeNormFEG[iBin]*eeNormFEG[iBin] + eeMEEBkgShape[iBin]*eeMEEBkgShape[iBin]);
    ffTot[iBin] = sqrt(ffTot[iBin]*ffTot[iBin] + ffMinusEE[iBin]*ffMinusEE[iBin]);
    eeTot[iBin] = sqrt(eeTot[iBin]*eeTot[iBin] + eeMEEBkgShape[iBin]*eeMEEBkgShape[iBin]);
    eeEGStatSys[iBin] = sqrt(eeEGTot[iBin]*eeEGTot[iBin] + eeMEEBkgShape[iBin]*eeMEEBkgShape[iBin]);
    ffEGStatSys[iBin] = sqrt(ffEGTot[iBin]*ffEGTot[iBin] + ffMinusEE[iBin]*ffMinusEE[iBin]);
    ffEGStatTot[iBin] = sqrt(ffStatTot[iBin]*ffStatTot[iBin] + egStat[iBin]*egStat[iBin]);
    eeEGStatTot[iBin] = sqrt(eeStatTot[iBin]*eeStatTot[iBin] + egStat[iBin]*egStat[iBin]);
    ffEGSysTot[iBin] = sqrt(ffSys[iBin]*ffSys[iBin] + egFEG[iBin]*egFEG[iBin]);
    eeEGSysTot[iBin] = sqrt(eeSys[iBin]*eeSys[iBin] + egFEG[iBin]*egFEG[iBin]);
  }

  //print the ee table
  printTableLine("Total", binNums, eeTot, eeBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Statistics", binNums, eeStatReweighting, eeBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}No. events", binNums, eeStatTot, eeBkg);
  printTableLine("\\hspace{1.5cm}In norm. region", binNums, eeNormStat, eeBkg);
  printTableLine("\\hspace{1.5cm}In this \\MET bin", binNums, eeStat, eeBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}Reweighting", binNums, eeReweightingTot, eeBkg);
  printTableLine("\\hspace{1.5cm}In norm. region", binNums, eeNormReweighting, eeBkg);
  printTableLine("\\hspace{1.5cm}In this \\MET bin", binNums, eeReweighting, eeBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Systematics", binNums, eeSys, eeBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}$f_{e\\rightarrow\\gamma}$ (in norm. region)", binNums, eeNormFEG, 
		 eeBkg);
  printTableLine("\\hspace{1cm}$m_{ee}$ background shape", binNums, eeMEEBkgShape, eeBkg); //update
  cout << "\\hline\n";

  //print the ff table
  printTableLine("Total", binNums, ffTot, ffBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Statistics", binNums, ffStatReweighting, ffBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}No. events", binNums, ffStatTot, ffBkg);
  printTableLine("\\hspace{1.5cm}In norm. region", binNums, ffNormStat, ffBkg);
  printTableLine("\\hspace{1.5cm}In this \\MET bin", binNums, ffStat, ffBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}Reweighting", binNums, ffReweightingTot, ffBkg);
  printTableLine("\\hspace{1.5cm}In norm. region", binNums, ffNormReweighting, ffBkg);
  printTableLine("\\hspace{1.5cm}In this \\MET bin", binNums, ffReweighting, ffBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Systematics", binNums, ffSys, ffBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}$ee$/$\\mathit{ff}$ difference", binNums, ffMinusEE, ffBkg);
  printTableLine("\\hspace{1cm}$f_{e\\rightarrow\\gamma}$ (in norm. region)", binNums, ffNormFEG, 
		 ffBkg);
  cout << "\\hline\n";

  //print the eg table
  printTableLine("Total", binNums, egTot, egBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Statistics", binNums, egStat, egBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Systematics ($f_{e\\rightarrow\\gamma}$)", binNums, egFEG, egBkg);
  cout << "\\hline\n";

  //print the total bkg. table
  printTableLine("Total ($ee$ + $e\\gamma$)", binNums, eeEGStatSys, eeEGBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Statistics", binNums, eeEGStatTot, eeEGBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}QCD", binNums, eeStatTot, eeEGBkg);
  printTableLine("\\hspace{1cm}Electroweak", binNums, egStat, eeEGBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Systematics", binNums, eeEGSysTot, eeEGBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}QCD", binNums, eeSys, eeEGBkg);
  printTableLine("\\hspace{1cm}Electroweak", binNums, egFEG, eeEGBkg);
  cout << "\\hline\n";
  printTableLine("Total ($\\mathit{ff}$ + $e\\gamma$)", binNums, ffEGStatSys, ffEGBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Statistics", binNums, ffEGStatTot, ffEGBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}QCD", binNums, ffStatTot, ffEGBkg);
  printTableLine("\\hspace{1cm}Electroweak", binNums, egStat, ffEGBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{0.5cm}Systematics", binNums, ffEGSysTot, ffEGBkg);
  cout << "\\hline\n";
  printTableLine("\\hspace{1cm}QCD", binNums, ffSys, ffEGBkg);
  printTableLine("\\hspace{1cm}Electroweak", binNums, egFEG, ffEGBkg);
  cout << "\\hline\n";
}

void compareDataToMC()
{
  string samples[3] = {"gg", "ee", "ff"};
  TFile data("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_correctEESidebandSubtraction.root");
  TFile MC("/Users/rachelyohay/RA3/data/4684pb-1_MET_MC_ggDiPhotonJets_noTT.root");
  TFile out("/Users/rachelyohay/RA3/data/dataMCComparisonDiPhotonJetsNoTTMC.root", "RECREATE");
  vector<vector<TH1F*> > diEMPTData(3, vector<TH1F*>(3, NULL));
  vector<vector<TH1F*> > diEMPTMC(3, vector<TH1F*>(3, NULL));
  vector<vector<TH1F*> > weightsData(3, vector<TH1F*>(3, NULL));
  vector<vector<TH1F*> > weightsMC(3, vector<TH1F*>(3, NULL));
  vector<vector<TCanvas*> > diEMPT(3, vector<TCanvas*>(3, NULL));
  vector<vector<TCanvas*> > diEMPTRatio(3, vector<TCanvas*>(3, NULL));
  vector<vector<TCanvas*> > weights(3, vector<TCanvas*>(3, NULL));
  vector<vector<TLegend*> > diEMPTLegend(3, vector<TLegend*>(3, NULL));
  vector<vector<TLegend*> > weightsLegend(3, vector<TLegend*>(3, NULL));
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      stringstream nameDiEMPT;
      nameDiEMPT << samples[i] << "DiEMPT" << j << "j";
      diEMPT[i][j] = new TCanvas(nameDiEMPT.str().c_str(), "", 600, 600);
      diEMPT[i][j]->SetFillStyle(0);
      diEMPT[i][j]->SetFillColor(0);
      diEMPT[i][j]->SetGrid();
      stringstream nameDiEMPTRatio;
      nameDiEMPTRatio << samples[i] << "DiEMPTRatio" << j << "j";
      diEMPTRatio[i][j] = new TCanvas(nameDiEMPTRatio.str().c_str(), "", 600, 300);
      diEMPTRatio[i][j]->SetFillStyle(0);
      diEMPTRatio[i][j]->SetFillColor(0);
      diEMPTRatio[i][j]->SetGrid();
      stringstream nameWeightsCanvas;
      nameWeightsCanvas << samples[i] << "Weights" << j << "j";
      weights[i][j] = new TCanvas(nameWeightsCanvas.str().c_str(), "", 600, 600);
      weights[i][j]->SetFillStyle(0);
      weights[i][j]->SetFillColor(0);
      weights[i][j]->SetGrid();
      stringstream nameDiEMPTLegend;
      nameDiEMPTLegend << samples[i] << ", " << j << " jets";
      diEMPTLegend[i][j] = new TLegend(0.6, 0.6, 0.9, 0.8);
      diEMPTLegend[i][j]->SetFillColor(0);
      diEMPTLegend[i][j]->SetHeader(nameDiEMPTLegend.str().c_str());
      stringstream nameWeightsLegend;
      nameWeightsLegend << samples[i] << ", " << j << " jets";
      weightsLegend[i][j] = new TLegend(0.6, 0.6, 0.9, 0.8);
      weightsLegend[i][j]->SetFillColor(0);
      weightsLegend[i][j]->SetHeader(nameWeightsLegend.str().c_str());
      stringstream name;
      name << samples[i] << "DiEMETUniformNjBin" << (j + 1);
      data.GetObject(name.str().c_str(), diEMPTData[i][j]);
      MC.GetObject(name.str().c_str(), diEMPTMC[i][j]);
      stringstream nameWeights;
      nameWeights << samples[i] << "Weights_METBin1_NjBin" << (j + 1);
      if (i != 0) {
	data.GetObject(nameWeights.str().c_str(), weightsData[i][j]);
	MC.GetObject(nameWeights.str().c_str(), weightsMC[i][j]);
      }
      diEMPTLegend[i][j]->AddEntry(diEMPTData[i][j], "Data", "lp");
      diEMPTLegend[i][j]->AddEntry(diEMPTMC[i][j], "MC", "lp");
      weightsLegend[i][j]->AddEntry(weightsData[i][j], "Data", "lp");
      weightsLegend[i][j]->AddEntry(weightsMC[i][j], "MC", "lp");
      out.cd();
      if (i != 0) {
	weights[i][j]->cd();
	weightsData[i][j]->SetMarkerColor(kRed);
	weightsData[i][j]->SetMarkerSize(0.5);
	weightsData[i][j]->SetLineColor(kRed);
	weightsData[i][j]->GetXaxis()->SetTitle("Di-EM p_{T} (GeV)");
	weightsMC[i][j]->SetMarkerColor(kBlue);
	weightsMC[i][j]->SetMarkerSize(0.5);
	weightsMC[i][j]->SetLineColor(kBlue);
	weightsMC[i][j]->GetXaxis()->SetTitle("Di-EM p_{T} (GeV)");
	weightsData[i][j]->Draw();
	weightsMC[i][j]->Draw("SAME");
	weightsLegend[i][j]->Draw();
	weights[i][j]->Write();
      }
      diEMPT[i][j]->cd();
      diEMPTData[i][j]->SetMarkerColor(kRed);
      diEMPTData[i][j]->SetMarkerSize(0.5);
      diEMPTData[i][j]->SetLineColor(kRed);
      diEMPTData[i][j]->GetXaxis()->SetTitle("Di-EM p_{T} (GeV)");
      diEMPTMC[i][j]->SetMarkerColor(kBlue);
      diEMPTMC[i][j]->SetMarkerSize(0.5);
      diEMPTMC[i][j]->SetLineColor(kBlue);
      diEMPTMC[i][j]->GetXaxis()->SetTitle("Di-EM p_{T} (GeV)");
      diEMPTMC[i][j]->Scale(diEMPTData[i][j]->Integral()/diEMPTMC[i][j]->Integral());
      diEMPTData[i][j]->Draw();
      diEMPTMC[i][j]->Draw("SAME");
      diEMPTLegend[i][j]->Draw();
      diEMPT[i][j]->Write();
      diEMPTRatio[i][j]->cd();
      TH1F* diEMPTRatioDist = (TH1F*)diEMPTData[i][j]->Clone();
      diEMPTRatioDist->Divide(diEMPTMC[i][j]);
      diEMPTRatioDist->SetMarkerColor(kBlack);
      diEMPTRatioDist->SetMarkerSize(0.5);
      diEMPTRatioDist->SetLineColor(kBlack);
      diEMPTRatioDist->GetXaxis()->SetTitle("Di-EM p_{T} (GeV)");
      diEMPTRatioDist->GetYaxis()->SetTitle("#frac{N_{data}}{N_{MC}}");
      diEMPTRatioDist->Draw();
      diEMPTRatio[i][j]->Write();
    }
  }
  out.Close();
  data.Close();
  MC.Close();
}

void compare2Methods()
{
  TFile in1("/Users/rachelyohay/RA3/data/eeVsGGMCDiPhotonJets.root");
  TFile in2("/Users/rachelyohay/RA3/data/eeVsGGMCDiPhotonJetsNoReweighting.root");
  TFile out("/Users/rachelyohay/RA3/data/eeVsGG_reweightingVsNone_MC.root", "RECREATE");
  TCanvas* obs1Canvas = NULL;
  TCanvas* obs2Canvas = NULL;
  TCanvas* obs = new TCanvas("obs", "", 600, 300);
  obs->SetFillStyle(0);
  obs->SetFillColor(0);
  obs->SetGrid();
  TLegend* obsLegend = new TLegend(0.6, 0.6, 0.9, 0.8);
  obsLegend->SetFillColor(0);
  in1.GetObject("METRatio", obs1Canvas);
  in2.GetObject("METRatio", obs2Canvas);
  TH1F* obs1 = (TH1F*)obs1Canvas->GetPrimitive("ggMET");
  TH1F* obs2 = (TH1F*)obs2Canvas->GetPrimitive("ggMET");
  obsLegend->AddEntry(obs1, "With dijet p_{T} reweighting", "lp");
  obsLegend->AddEntry(obs2, "Without dijet p_{T} reweighting", "lp");
  out.cd();
  obs->cd();
  obs1->SetMarkerColor(kRed);
  obs1->SetMarkerSize(0.5);
  obs1->SetMarkerStyle(20);
  obs1->SetLineColor(kRed);
  obs1->GetXaxis()->SetTitle("ME_{T} (GeV)");
  obs2->SetMarkerColor(kBlue);
  obs2->SetMarkerSize(0.5);
  obs2->SetMarkerStyle(20);
  obs2->SetLineColor(kBlue);
  obs2->GetXaxis()->SetTitle("ME_{T} (GeV)");
  obs2->Draw();
  obs1->Draw("SAME");
  obsLegend->Draw();
  obs->Write();
  out.Close();
  in1.Close();
  in2.Close();
}

void compare3Methods()
{
  TFile in1("/Users/rachelyohay/RA3/data/eeVsGG_Pythia_rhoL4.root");
  TFile in2("/Users/rachelyohay/RA3/data/eeVsGG_Pythia_4LeqRhoL8.root");
  TFile in3("/Users/rachelyohay/RA3/data/eeVsGG_Pythia_rhoGeq8.root");
  TFile out("/Users/rachelyohay/RA3/data/eeVsGG_Pythia_rhoL4Vs4LeqRhoL8VsRhoGeq8.root", "RECREATE");
  TCanvas* obs1Canvas = NULL;
  TCanvas* obs2Canvas = NULL;
  TCanvas* obs3Canvas = NULL;
  TCanvas* obs = new TCanvas("obs", "", 600, 300);
  obs->SetFillStyle(0);
  obs->SetFillColor(0);
  obs->SetGrid();
  TLegend* obsLegend = new TLegend(0.6, 0.6, 0.9, 0.8);
  obsLegend->SetFillColor(0);
  in1.GetObject("METRatio", obs1Canvas);
  in2.GetObject("METRatio", obs2Canvas);
  in3.GetObject("METRatio", obs3Canvas);
  TH1F* obs1 = (TH1F*)obs1Canvas->GetPrimitive("ggMET");
  TH1F* obs2 = (TH1F*)obs2Canvas->GetPrimitive("ggMET");
  TH1F* obs3 = (TH1F*)obs3Canvas->GetPrimitive("ggMET");
  obsLegend->AddEntry(obs1, "#rho < 4.0", "lp");
  obsLegend->AddEntry(obs2, "4.0 #leq #rho < 8.0", "lp");
  obsLegend->AddEntry(obs3, "#rho #geq 8.0", "lp");
  // obsLegend->SetHeader("CMS 4.7 fb^{-1}, no jet requirement");
  obsLegend->SetHeader("CMS simulation, no jet requirement");
  out.cd();
  obs->cd();
  obs1->SetMarkerColor(kRed);
  obs1->SetMarkerSize(0.5);
  obs1->SetMarkerStyle(20);
  obs1->SetLineColor(kRed);
  obs1->GetXaxis()->SetTitle("ME_{T} (GeV)");
  obs2->SetMarkerColor(kBlue);
  obs2->SetMarkerStyle(20);
  obs2->SetMarkerSize(0.5);
  obs2->SetLineColor(kBlue);
  obs2->GetXaxis()->SetTitle("ME_{T} (GeV)");
  obs3->SetMarkerSize(0.5);
  obs3->SetMarkerStyle(20);
  obs3->SetLineColor(kBlack);
  obs3->GetXaxis()->SetTitle("ME_{T} (GeV)");
  obs1->Draw();
  obs2->Draw("SAME");
  obs3->Draw("SAME");
  obsLegend->Draw();
  obs->Write();
  out.Close();
  in1.Close();
  in2.Close();
  in3.Close();
}

void compareEEToGG1ToGG2()
{
  // TFile in1("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_trailingPhotonETPlots.root");
  // TFile in2("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_trailingPhotonETPlots.root");
  // TFile out("/Users/rachelyohay/RA3/data/eeVsGGData_rebinned.root", "RECREATE");
  TFile in1("/Users/rachelyohay/RA3/data/4684pb-1_MET_MC_ggPythia25-250GeV.root");
  TFile in2("/Users/rachelyohay/RA3/data/4684pb-1_MET_MC_ggMadGraph.root");
  TFile out("/Users/rachelyohay/RA3/data/eeVsGGPythia25-250GeVVsGGMadGraph_rebinned.root", 
	    "RECREATE");
  const string samples[3] = {"gg", "gg", "ee"};
  TFile* files[3] = {&in1, &in2, &in2};
  vector<TH1F*> leadingPhotonET(3, NULL);
  vector<TH1F*> trailingPhotonET(3, NULL);
  vector<vector<TH1F*> > diEMPT(3, vector<TH1F*>(3, NULL));
  vector<vector<TH1F*> > dijetPT(3, vector<TH1F*>(3, NULL));
  TCanvas* leadingPhotonETCanvas = new TCanvas("leadingPhotonETCanvas", "", 600, 600);
  leadingPhotonETCanvas->SetFillStyle(0);
  leadingPhotonETCanvas->SetFillColor(0);
  leadingPhotonETCanvas->SetGrid();
  leadingPhotonETCanvas->SetLogy();
  TCanvas* trailingPhotonETCanvas = new TCanvas("trailingPhotonETCanvas", "", 600, 600);
  trailingPhotonETCanvas->SetFillStyle(0);
  trailingPhotonETCanvas->SetFillColor(0);
  trailingPhotonETCanvas->SetGrid();
  trailingPhotonETCanvas->SetLogy();
  vector<TCanvas*> diEMPTCanvas(3, NULL);
  vector<TCanvas*> dijetPTCanvas(3, NULL);
  // TLegend* leadingPhotonETLegend = new TLegend(0.5, 0.6, 0.9, 0.9);
  TLegend* leadingPhotonETLegend = new TLegend(0.2, 0.17, 0.6, 0.47);
  leadingPhotonETLegend->SetFillColor(0);
  leadingPhotonETLegend->SetHeader("CMS simulation, no jet requirement");
  // leadingPhotonETLegend->SetHeader("CMS 4.7 fb^{-1}, no jet requirement");
  TLegend* trailingPhotonETLegend = new TLegend(0.5, 0.6, 0.9, 0.9);
  trailingPhotonETLegend->SetFillColor(0);
  trailingPhotonETLegend->SetHeader("CMS simulation, no jet requirement");
  // trailingPhotonETLegend->SetHeader("CMS 4.7 fb^{-1}, no jet requirement");
  vector<TLegend*> diEMPTLegend(3, NULL);
  vector<TLegend*> dijetPTLegend(3, NULL);
  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      string entryName;
      string drawOpt = "SAME";
      switch (i) {
      case 0:
	entryName = "#gamma#gamma (Pythia 25 GeV #leq #hat{p}_{T} < 250 GeV)";
	// entryName = "#gamma#gamma";
	drawOpt = "";
	break;
      case 1:
	entryName = "#gamma#gamma (MadGraph)";
	// entryName = "#gamma#gamma";
	break;
      default:
	entryName = "ee (Pythia)";
	// entryName = "ee";
	break;
      }
      if (i == 0) {
	stringstream nameDiEMPTCanvas;
	nameDiEMPTCanvas << "diEMPT" << j << "j";
	diEMPTCanvas[j] = new TCanvas(nameDiEMPTCanvas.str().c_str(), "", 600, 600);
	diEMPTCanvas[j]->SetFillStyle(0);
	diEMPTCanvas[j]->SetFillColor(0);
	diEMPTCanvas[j]->SetGrid();
	if (j == 0) diEMPTCanvas[j]->SetLogy();
	stringstream nameDijetPTCanvas;
	nameDijetPTCanvas << "dijetPT" << j << "j";
	dijetPTCanvas[j] = new TCanvas(nameDijetPTCanvas.str().c_str(), "", 600, 600);
	dijetPTCanvas[j]->SetFillStyle(0);
	dijetPTCanvas[j]->SetFillColor(0);
	dijetPTCanvas[j]->SetGrid();
	if (j == 0) dijetPTCanvas[j]->SetLogy();
	stringstream nameDiEMPTLegend;
	nameDiEMPTLegend << "CMS simulation, " << j << " jets";
	// nameDiEMPTLegend << "CMS 4.7 fb^{-1}, " << j << " jets";
	// diEMPTLegend[j] = new TLegend(0.6, 0.6, 0.9, 0.8);
	diEMPTLegend[j] = new TLegend(0.3, 0.6, 0.9, 0.8);
	diEMPTLegend[j]->SetFillColor(0);
	diEMPTLegend[j]->SetHeader(nameDiEMPTLegend.str().c_str());
	// dijetPTLegend[j] = new TLegend(0.6, 0.6, 0.9, 0.8);
	dijetPTLegend[j] = new TLegend(0.3, 0.6, 0.9, 0.8);
	dijetPTLegend[j]->SetFillColor(0);
	dijetPTLegend[j]->SetHeader(nameDiEMPTLegend.str().c_str());
      }
      if (j == 0) {
	files[i]->GetObject((samples[i] + "LeadingPhotonET").c_str(), leadingPhotonET[i]);
	files[i]->GetObject((samples[i] + "TrailingPhotonET").c_str(), trailingPhotonET[i]);
      }
      stringstream nameDiEMPT;
      nameDiEMPT << samples[i] << "DiEMETUniformNjBin" << (j + 1);
      files[i]->GetObject(nameDiEMPT.str().c_str(), diEMPT[i][j]);
      stringstream nameDijetPT;
      nameDijetPT << samples[i] << "DijetETUniformNjBin" << (j + 1);
      files[i]->GetObject(nameDijetPT.str().c_str(), dijetPT[i][j]);
      // const float leadingPhotonETScale = 
      // 	leadingPhotonET[0]->Integral()/leadingPhotonET[i]->Integral();
      // const float trailingPhotonETScale = 
      // 	trailingPhotonET[0]->Integral()/trailingPhotonET[i]->Integral();
      // const float diEMPTScale = diEMPT[0][j]->Integral()/diEMPT[i][j]->Integral();
      // const float dijetPTScale = dijetPT[0][j]->Integral()/dijetPT[i][j]->Integral();
      const float leadingPhotonETScale = 1.0/leadingPhotonET[i]->Integral();
      const float trailingPhotonETScale = 1.0/trailingPhotonET[i]->Integral();
      const float diEMPTScale = 1.0/diEMPT[i][j]->Integral();
      const float dijetPTScale = 1.0/dijetPT[i][j]->Integral();
      // if (i != 1) {
      if (j == 0) {
	leadingPhotonETLegend->AddEntry(leadingPhotonET[i], entryName.c_str(), "lp");
	trailingPhotonETLegend->AddEntry(trailingPhotonET[i], entryName.c_str(), "lp");
      }
      diEMPTLegend[j]->AddEntry(diEMPT[i][j], entryName.c_str(), "lp");
      dijetPTLegend[j]->AddEntry(dijetPT[i][j], entryName.c_str(), "lp");
      // }
      out.cd();
      // if (i != 1) {
      if (j == 0) {
	leadingPhotonETCanvas->cd();
	leadingPhotonET[i]->Rebin(3);
	leadingPhotonET[i]->Scale(leadingPhotonETScale);
	leadingPhotonET[i]->SetMarkerColor(i + 1);
	leadingPhotonET[i]->SetMarkerSize(0.5);
	leadingPhotonET[i]->SetLineColor(i + 1);
	leadingPhotonET[i]->GetXaxis()->SetTitle("E_{T} (GeV)");
	leadingPhotonET[i]->Draw(drawOpt.c_str());
	leadingPhotonETLegend->Draw();
	trailingPhotonETCanvas->cd();
	trailingPhotonET[i]->Rebin(3);
	trailingPhotonET[i]->Scale(trailingPhotonETScale);
	trailingPhotonET[i]->SetMarkerColor(i + 1);
	trailingPhotonET[i]->SetMarkerSize(0.5);
	trailingPhotonET[i]->SetLineColor(i + 1);
	trailingPhotonET[i]->GetXaxis()->SetTitle("E_{T} (GeV)");
	trailingPhotonET[i]->Draw(drawOpt.c_str());
	trailingPhotonETLegend->Draw();
      }
      diEMPTCanvas[j]->cd();
      diEMPT[i][j]->Rebin(5);
      diEMPT[i][j]->Scale(diEMPTScale);
      diEMPT[i][j]->SetMarkerColor(i + 1);
      diEMPT[i][j]->SetMarkerSize(0.5);
      diEMPT[i][j]->SetLineColor(i + 1);
      diEMPT[i][j]->GetXaxis()->SetTitle("Di-EM p_{T} (GeV)");
      diEMPT[i][j]->Draw(drawOpt.c_str());
      diEMPTLegend[j]->Draw();
      dijetPTCanvas[j]->cd();
      dijetPT[i][j]->Rebin(5);
      dijetPT[i][j]->Scale(dijetPTScale);
      dijetPT[i][j]->SetMarkerColor(i + 1);
      dijetPT[i][j]->SetMarkerSize(0.5);
      dijetPT[i][j]->SetLineColor(i + 1);
      dijetPT[i][j]->GetXaxis()->SetTitle("Dijet p_{T} (GeV)");
      dijetPT[i][j]->Draw(drawOpt.c_str());
      dijetPTLegend[j]->Draw();
      // }
      if (i == 2) {
	if (j == 0) {
	  leadingPhotonETCanvas->Write();
	  // leadingPhotonETCanvas->SaveAs("/Users/rachelyohay/RA3/data/ee_MC_closure_plots/leadingPhotonET_data_norm1_rebinned.pdf");
	  leadingPhotonETCanvas->SaveAs("/Users/rachelyohay/RA3/data/ee_MC_closure_plots/leadingPhotonET_MCPythia25-250GeV_norm1_rebinned.pdf");
	  trailingPhotonETCanvas->Write();
	  // trailingPhotonETCanvas->SaveAs("/Users/rachelyohay/RA3/data/ee_MC_closure_plots/trailingPhotonET_data_norm1_rebinned.pdf");
	  trailingPhotonETCanvas->SaveAs("/Users/rachelyohay/RA3/data/ee_MC_closure_plots/trailingPhotonET_MCPythia25-250GeV_norm1_rebinned.pdf");
	}
	diEMPTCanvas[j]->Write();
	stringstream diEMPTFileName;
	diEMPTFileName << "/Users/rachelyohay/RA3/data/ee_MC_closure_plots/diEMPT_" << j;
	// diEMPTFileName << "j_data_norm1_rebinned.pdf";
	diEMPTFileName << "j_MCPythia25-250GeV_norm1_rebinned.pdf";
	diEMPTCanvas[j]->SaveAs(diEMPTFileName.str().c_str());
	stringstream dijetPTFileName;
	dijetPTFileName << "/Users/rachelyohay/RA3/data/ee_MC_closure_plots/dijetPT_" << j;
	// dijetPTFileName << "j_data_norm1_rebinned.pdf";
	dijetPTFileName << "j_MCPythia25-250GeV_norm1_rebinned.pdf";
	dijetPTCanvas[j]->SaveAs(dijetPTFileName.str().c_str());
	dijetPTCanvas[j]->Write();
      }
    }
  }
  out.Close();
  in1.Close();
  in2.Close();
}

void compare3Samples()
{
  TFile in("/Users/rachelyohay/RA3/data/4684pb-1_MET_18-Jan-12_skim_jetActivityHists.root");
  TFile out("/Users/rachelyohay/RA3/data/jetActivity_data.root", "RECREATE");
  const string variable[6] = {"CleanHT", "CleanMHT", "CleanNj", "CleanLeadingJetET", "Rho", "NPV"};
  TCanvas* canvas[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
  TLegend* legend[6][2] = {{NULL, NULL}, {NULL, NULL}, {NULL, NULL}, 
			   {NULL, NULL}, {NULL, NULL}, {NULL, NULL}};
  const char* xAxisTitle[6] = {"H_{T} (GeV)", "MH_{T} (GeV)", "N_{j}", "E_{T} (GeV)", 
			       "GeV/(unit #eta#upoint#phi)", "N_{PV}"};
  const string file[6] = {"HT", "MHT", "Nj", "j1ET", "rho", "nPV"};
  const string sample[3] = {"gg", "ee", "ff"};
  const char* legendEntry[3] = {"#gamma#gamma", "ee", "ff"};
  Color_t color[3] = {kBlack, kRed, kBlue};
  Option_t* drawOption[3] = {"", "SAME", "SAME"};
  TH1F* histogram[6][3] = {{NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}, 
			   {NULL, NULL, NULL}, {NULL, NULL, NULL}, {NULL, NULL, NULL}};
  TH1F* ggOverControl[6][2] = {{NULL, NULL}, {NULL, NULL}, {NULL, NULL}, 
			       {NULL, NULL}, {NULL, NULL}, {NULL, NULL}};
  for (unsigned int i = 0; i < 6; ++i) {
    canvas[i] = new TCanvas((variable[i] + "Canvas").c_str(), "", 400, 600);
    canvas[i]->Divide(1, 2);
    canvas[i]->cd(1)->SetPad(0.0, 0.33, 1.0, 1.0);
    setCanvasOptions(canvas[i]->cd(1), 1, 1, 0);
    canvas[i]->cd(2)->SetPad(0.0, 0.0, 1.0, 0.33);
    setCanvasOptions(canvas[i]->cd(2), 1, 0, 0);
    if (i == 4) legend[i][0] = new TLegend(0.4, 0.3, 0.7, 0.5);
    else legend[i][0] = new TLegend(0.6, 0.7, 0.9, 0.9);
    setLegendOptions(legend[i][0], "CMS 4.7 fb^{-1}, no jet requirement");
    legend[i][1] = new TLegend(0.2, 0.67, 0.5, 0.87);
    setLegendOptions(legend[i][1], "CMS 4.7 fb^{-1}, no jet requirement");
    for (unsigned int j = 0; j < 3; ++j) {
      in.GetObject((sample[j] + variable[i]).c_str(), histogram[i][j]);
      out.cd();
      canvas[i]->cd(1);
      setHistogramOptions(dynamic_cast<TH1*>(histogram[i][j]), color[j], 0.5, 1, 
			  histogram[i][0]->Integral()/histogram[i][j]->Integral(), 
			  xAxisTitle[i], "");
      histogram[i][j]->Draw(drawOption[j]);
      legend[i][0]->AddEntry(histogram[i][j], legendEntry[j], "lp");
      if (j > 0) {
	canvas[i]->cd(2);
	Int_t nBins = histogram[i][j]->GetNbinsX();
	ggOverControl[i][j - 1] = 
	  new TH1F(("gg_over_" + sample[j] + "_" + variable[i]).c_str(), "", nBins, 
		   histogram[i][j]->GetBinLowEdge(1), 
		   histogram[i][j]->GetBinLowEdge(nBins) + histogram[i][j]->GetBinWidth(nBins));
	ggOverControl[i][j - 1]->Divide(histogram[i][0], histogram[i][j]);
	setHistogramOptions(dynamic_cast<TH1*>(ggOverControl[i][j - 1]), color[j], 0.5, 1, 1.0, 
			    xAxisTitle[i], 
			    ("#frac{N_{#gamma#gamma}}{N_{" + sample[j] + "}}").c_str());
	ggOverControl[i][j - 1]->Draw(drawOption[j - 1]);
	legend[i][1]->AddEntry(ggOverControl[i][j - 1], legendEntry[j], "lp");
      }
    }
    canvas[i]->cd(1);
    legend[i][0]->Draw();
    canvas[i]->cd(2);
    legend[i][1]->Draw();
    canvas[i]->Write();
    canvas[i]->SaveAs(("/Users/rachelyohay/Documents/UVa/dissertation/hadronic_activity_" + 
		       file[i] + ".pdf").c_str());
  }
  out.Write();
  out.Close();
  in.Close();
}

float fitDiEMInvariantMass(vector<TH3F*>& METVsDiEMETVsInvMass, float& sigFracErr, 
			   map<string, float>& diEMInvMassFitParLowerBounds, 
			   map<string, float>& diEMInvMassFitParUpperBounds, 
			   map<string, float>& diEMInvMassFitPars, 
			   map<string, string>& diEMInvMassFitParUnits, TFile& out)
{
  //vectors to hold the parameter names and units
  vector<string> expectedPars;
  vector<string> expectedParUnits;

  //expected parameters
  expectedPars.push_back("m"); //mass variable
  expectedPars.push_back("muCB"); //Z signal Crystal Ball mean
  expectedPars.push_back("sigma"); //Z signal Crystal Ball width
  expectedPars.push_back("alphaCB"); //Z signal Crystal Ball power law turn-on
  expectedPars.push_back("n"); //Z signal Crystal Ball power
  expectedPars.push_back("alphaCMS"); //mee background RooCMSShape location of erf slope
  expectedPars.push_back("beta"); //mee background RooCMSShape erf slope parameter
  expectedPars.push_back("gamma"); //mee background RooCMSShape exponent
  expectedPars.push_back("muCMS"); //mee background RooCMSShape location of exponential

  //parameter units
  expectedParUnits.push_back("GeV");
  expectedParUnits.push_back("GeV");
  expectedParUnits.push_back("GeV");
  expectedParUnits.push_back("");
  expectedParUnits.push_back("");
  expectedParUnits.push_back("GeV");
  expectedParUnits.push_back("");
  expectedParUnits.push_back("");
  expectedParUnits.push_back("GeV");

  //map of RooRealVar objects constructed from parameters
  map<string, RooRealVar*> realVars;

  //loop over parameters
  for (vector<string>::const_iterator iPar = expectedPars.begin(); iPar != expectedPars.end(); 
       ++iPar) {
    map<string, float>::const_iterator iExpectedParLowerBound = 
      diEMInvMassFitParLowerBounds.find(*iPar);
    map<string, float>::const_iterator iExpectedParUpperBound = 
      diEMInvMassFitParUpperBounds.find(*iPar);
    map<string, float>::const_iterator iExpectedPar = diEMInvMassFitPars.find(*iPar);
    map<string, string>::const_iterator iExpectedParUnit = diEMInvMassFitParUnits.find(*iPar);

    //if the expected parameters are missing, quit
    if (iExpectedParLowerBound == diEMInvMassFitParLowerBounds.end()) {
      cerr << "Error: lower bound for parameter " << *iPar << " not found.\n";
      return -1.0;
    }
    if (iExpectedParUpperBound == diEMInvMassFitParUpperBounds.end()) {
      cerr << "Error: upper bound for parameter " << *iPar << " not found.\n";
      return -1.0;
    }
    if ((iExpectedPar == diEMInvMassFitPars.end()) && (*iPar != "m")) {
      cerr << "Error: parameter " << *iPar << " not found.\n";
      return -1.0;
    }
    if (iExpectedParUnit == diEMInvMassFitParUnits.end()) {
      cerr << "Error: parameter " << *iPar << " not found.\n";
      return -1.0;
    }

    //create the RooRealVars
    const char* realVarName = (*iPar).c_str();
    if (*iPar != "m") {
      realVars[*iPar] = new RooRealVar(realVarName, realVarName, iExpectedPar->second, 
				       iExpectedParLowerBound->second, 
				       iExpectedParUpperBound->second, 
				       iExpectedParUnit->second.c_str());
    }
    else {
      realVars[*iPar] = new RooRealVar(realVarName, realVarName, iExpectedParLowerBound->second, 
				       iExpectedParUpperBound->second, 
				       iExpectedParUnit->second.c_str());
    }
  }

  //signal and background yields (to be fitted)
  RooRealVar signalYield("signalYield", "signalYield", 0.0, 1.0);
  RooRealVar backgroundYield("backgroundYield", "backgroundYield", 0.0, 1.0);

  //signal PDF: Crystal Ball
  RooCBShape signal("signal", "signal", *(realVars["m"]), *(realVars["muCB"]), 
		    *(realVars["sigma"]), *(realVars["alphaCB"]), *(realVars["n"]));

  //background PDF: RooCMSShape
  RooCMSShape background("background", "background", *(realVars["m"]), *(realVars["alphaCMS"]), 
			 *(realVars["beta"]), *(realVars["gamma"]), *(realVars["muCMS"]));

  //signal + background PDF
  RooAddPdf total("total", "total", RooArgList(signal, background), 
  		  RooArgList(signalYield));

  //check that all input histograms have the same number of bins
  const unsigned int nBins = METVsDiEMETVsInvMass[0]->GetNbinsX();
  for (vector<TH3F*>::const_iterator iHist = METVsDiEMETVsInvMass.begin() + 1; 
       iHist != METVsDiEMETVsInvMass.end(); ++iHist) {
    const Int_t nBinsPrev = (*(iHist - 1))->GetNbinsX();
    const Int_t nBinsCurr = (*iHist)->GetNbinsX();
    if (nBinsPrev != nBinsCurr) {
      cerr << "Error: " << nBinsPrev << " bins in histogram ";
      cerr << (iHist - 1 - METVsDiEMETVsInvMass.begin()) << " but " << nBinsCurr;
      cerr << " bins in histogram " << (iHist - METVsDiEMETVsInvMass.begin()) << ".\n";
      return -1.0;
    }
  }

  //make invariant mass projections
  vector<TH1D*> invMassHists;
  for (vector<TH3F*>::const_iterator iHist = METVsDiEMETVsInvMass.begin(); 
       iHist != METVsDiEMETVsInvMass.end(); ++iHist) {
    invMassHists.push_back((*iHist)->ProjectionX((string((*iHist)->GetName()) + "_px").c_str(), 
						 0, -1, 0, -1, "e"));
  }

  //check that all input histograms have the same binning
  vector<float> binBoundaries;
  for (unsigned int iBin = 1; iBin <= nBins + 1; ++iBin) {
    const float binBoundary = invMassHists[0]->GetBinLowEdge(iBin);
    for (vector<TH1D*>::const_iterator iHist = invMassHists.begin() + 1; 
	 iHist != invMassHists.end(); ++iHist) {
      const float binBoundaryPrev = (*(iHist - 1))->GetBinLowEdge(iBin);
      const float binBoundaryCurr = (*iHist)->GetBinLowEdge(iBin);
      if (binBoundaryPrev != binBoundaryCurr) {
	cerr << "Error: for bin " << iBin << ", lower edge is " << binBoundaryPrev;
	cerr << " in histogram " << (iHist - 1 - invMassHists.begin()) << " but ";
	cerr << binBoundaryCurr << " in histogram " << (iHist - invMassHists.begin()) << ".\n";
	return -1.0;
      }
    }
    binBoundaries.push_back(binBoundary);
  }

  //add input histograms if necessary
  string invMassHistName = 
    "m" + string(METVsDiEMETVsInvMass[0]->GetName()).substr(0, 2) + "Hist";
  TH1F invMassHist(invMassHistName.c_str(), "", nBins, &binBoundaries[0]);
  invMassHist.Sumw2();
  for (vector<TH1D*>::const_iterator iHist = invMassHists.begin(); iHist != invMassHists.end(); 
       ++iHist) { invMassHist.Add(*iHist); }

  //fit
  string diEMInvMassName = invMassHistName.substr(0, 3);
  RooDataHist 
    invMassDataHist((diEMInvMassName + "DataHist").c_str(), "", *(realVars["m"]), &invMassHist);
  RooFitResult* fitResult = total.fitTo(invMassDataHist, Save(true));

  //plot result
  out.cd();
  TCanvas canvas((diEMInvMassName + "Canvas").c_str(), "", 600, 600);
  string sample(invMassHistName.substr(1, 2));
  string sampleLatex = "ee";
  if (sample == "eg") sampleLatex = "e#gamma";
  setCanvasOptions(&canvas, 1, 0, 0);
  canvas.cd();
  RooPlot* mPlot = realVars["m"]->frame(Name((diEMInvMassName + "RooPlot").c_str()));
  total.paramOn(mPlot, Format("NEU", AutoPrecision(1)), Layout(0.6, 0.98, 0.9));
  invMassDataHist.plotOn(mPlot, LineColor(kBlack));
  total.plotOn(mPlot, Name((diEMInvMassName + "BkgFit").c_str()), Components("background"), 
	       LineColor(kRed));
  total.plotOn(mPlot, Name((diEMInvMassName + "SigPlusBkgFit").c_str()), LineColor(kBlue));
  mPlot->Draw();
  canvas.Write();

  //set the signal fraction and its error
  const RooArgList& parsFinal = fitResult->floatParsFinal();
  const Int_t iSignalYield = parsFinal.index("signalYield");
  if (iSignalYield == -1) {
    cerr << "Error: parameter \"signalYield\" not found in fitResult->floatParsFinal().\n";
    return -1.0;
  }
  sigFracErr = ((const RooRealVar&)parsFinal[iSignalYield]).getError();
  const float sigFrac = ((const RooRealVar&)parsFinal[iSignalYield]).getVal();

  //deallocate memory
  for (map<string, RooRealVar*>::iterator iVar = realVars.begin(); iVar != realVars.end(); 
       ++iVar) {
    delete iVar->second;
    iVar->second = NULL;
  }
  realVars.clear();

  //return the signal fraction
  return sigFrac;
}

float electronPhotonMisIDRate(vector<TH3F*>& eeMETVsDiEMETVsInvMass, 
			      vector<TH3F*>& egMETVsDiEMETVsInvMass, 
			      map<string, float>& eeInvMassFitParLowerBounds, 
			      map<string, float>& eeInvMassFitParUpperBounds, 
			      map<string, float>& eeInvMassFitPars, 
			      map<string, float>& egInvMassFitParLowerBounds, 
			      map<string, float>& egInvMassFitParUpperBounds, 
			      map<string, float>& egInvMassFitPars, 
			      map<string, string>& diEMInvMassFitParUnits, TFile& out, 
			      float& statErr)
{
  //perform the fits
  float eeZYieldErr = -1.0;
  float egZYieldErr = -1.0;
  float eeZYield = fitDiEMInvariantMass(eeMETVsDiEMETVsInvMass, eeZYieldErr, 
					eeInvMassFitParLowerBounds, 
					eeInvMassFitParUpperBounds, eeInvMassFitPars, 
					diEMInvMassFitParUnits, out);
  float egZYield = fitDiEMInvariantMass(egMETVsDiEMETVsInvMass, egZYieldErr, 
					egInvMassFitParLowerBounds, 
					egInvMassFitParUpperBounds, egInvMassFitPars, 
					diEMInvMassFitParUnits, out);

  //get the histogram normalizations
  unsigned int nEEEvts = 0;
  unsigned int nEGEvts = 0;
  for (vector<TH3F*>::const_iterator iEEHist = eeMETVsDiEMETVsInvMass.begin(); 
       iEEHist != eeMETVsDiEMETVsInvMass.end(); ++iEEHist) {
    nEEEvts+=(*iEEHist)->GetEntries();
  }
  for (vector<TH3F*>::const_iterator iEGHist = egMETVsDiEMETVsInvMass.begin(); 
       iEGHist != egMETVsDiEMETVsInvMass.end(); ++iEGHist) {
    nEGEvts+=(*iEGHist)->Integral(36, 51, 0, -1, 0, -1);
  }

  /*extract the mis-ID rate f_e-->g = N_eg/(2N_ee + N_eg) where N refers to the number of Z signal 
    events*/
  eeZYield*=nEEEvts;
  egZYield*=nEGEvts;
  float egMisIDRate = -1.0;
  const float denominator = 2*eeZYield + egZYield;
  if (denominator != 0.0) egMisIDRate = egZYield/denominator;
  else cerr << "Error: can't calculate electron-->photon mis-ID rate because denominator is 0.\n";

  /*extract the error on f_e-->g using the error on the signal yield fit parameter (due to stats 
    and fit bias, but much smaller than sys. error*/
  eeZYieldErr*=nEEEvts;
  egZYieldErr*=nEGEvts;
  if ((eeZYield > 0.0) && (egZYield > 0.0)) {
    statErr = 
      2*egMisIDRate*(eeZYield/denominator)*sqrt(((egZYieldErr*egZYieldErr)/(egZYield*egZYield)) + 
						((eeZYieldErr*eeZYieldErr)/(eeZYield*eeZYield)));
  }
  else {
    cerr << "Error: can't calculate statistical/fit error on electron-->photon mis-ID rate ";
    cerr << "because ee or eg Z signal yields are <=0.\n";
  }

  //return f_e-->g
  return egMisIDRate;
}

void calculateElectronPhotonMisIDRate(const char* inFile, const char* outFile)
{
  /*ee background parameters (from sequential ee fits with increasing numbers of parameters fixed):
    alphaCMS = 97.0 +/- 0.1
    beta = 0.0922 +/- 0.0003
    gamma = 0.191 +/- 0.001
    muCMS = 58 +/- 2
  */

  /*ee signal parameters (from sequential ee fits with increasing numbers of parameters fixed):
    alphaCB = 1.063 +/- 0.005
    n = 143.16 +/- 0.02
    sigma = [1.0, 5.0] (starting 1.8)
    muCB = [86.2, 96.2] (starting 91.2)
  */

  /*eg background parameters (from eg fit with all parameters floating):
    alphaCMS = 72.02 +/- 0.09
    beta = 0.098 +/- 0.001
    gamma = 0.0375 +/- 0.0004
    muCMS = 56 +/- 9
  */

  /*eg signal parameters (from ee fit):
    alphaCB = 1.063 +/- 0.005
    n = 143.16 +/- 0.02
    sigma = [1.0, 5.0] (starting 1.8)
    muCB = [86.2, 96.2] (starting 91.2)
  */

  //set ee fit parameters
  map<string, float> eeInvMassFitPars;
  eeInvMassFitPars["muCB"] = 91.2; //Z signal Crystal Ball mean
  eeInvMassFitPars["sigma"] = 1.8; //Z signal Crystal Ball width
  eeInvMassFitPars["alphaCB"] = 1.063; //Z signal Crystal Ball power law turn-on
  eeInvMassFitPars["n"] = 143.16; //Z signal Crystal Ball power
  eeInvMassFitPars["alphaCMS"] = 97.0; //mee background RooCMSShape location of erf slope
  eeInvMassFitPars["beta"] = 0.0922; //mee background RooCMSShape erf slope parameter
  eeInvMassFitPars["gamma"] = 0.191; //mee background RooCMSShape exponent
  eeInvMassFitPars["muCMS"] = 58.0; //mee background RooCMSShape location of exponential

  //set ee fit parameter lower bounds
  map<string, float> eeInvMassFitParLowerBounds;
  eeInvMassFitParLowerBounds["m"] = 71.0; //mass variable (no fixed value, only range)
  eeInvMassFitParLowerBounds["muCB"] = 86.2;
  eeInvMassFitParLowerBounds["sigma"] = 1.0;
  eeInvMassFitParLowerBounds["alphaCB"] = 1.063;
  eeInvMassFitParLowerBounds["n"] = 143.16;
  eeInvMassFitParLowerBounds["alphaCMS"] = 97.0;
  eeInvMassFitParLowerBounds["beta"] = 0.0922;
  eeInvMassFitParLowerBounds["gamma"] = 0.191;
  eeInvMassFitParLowerBounds["muCMS"] = 58.0;

  //set ee fit parameter upper bounds
  map<string, float> eeInvMassFitParUpperBounds;
  eeInvMassFitParUpperBounds["m"] = 111.0;
  eeInvMassFitParUpperBounds["muCB"] = 96.2;
  eeInvMassFitParUpperBounds["sigma"] = 5.0;
  eeInvMassFitParUpperBounds["alphaCB"] = 1.063;
  eeInvMassFitParUpperBounds["n"] = 143.16;
  eeInvMassFitParUpperBounds["alphaCMS"] = 97.0;
  eeInvMassFitParUpperBounds["beta"] = 0.0922;
  eeInvMassFitParUpperBounds["gamma"] = 0.191;
  eeInvMassFitParUpperBounds["muCMS"] = 58.0;

  //set eg fit parameters
  map<string, float> egInvMassFitPars;
  egInvMassFitPars["muCB"] = 91.2; //Z signal Crystal Ball mean
  egInvMassFitPars["sigma"] = 1.8; //Z signal Crystal Ball width
  egInvMassFitPars["alphaCB"] = 1.063; //Z signal Crystal Ball power law turn-on
  egInvMassFitPars["n"] = 143.16; //Z signal Crystal Ball power
  egInvMassFitPars["alphaCMS"] = 72.02; //meg background RooCMSShape location of erf slope
  egInvMassFitPars["beta"] = 0.098; //meg background RooCMSShape erf slope parameter
  egInvMassFitPars["gamma"] = 0.0375; //meg background RooCMSShape exponent
  egInvMassFitPars["muCMS"] = 56.0; //meg background RooCMSShape location of exponential

  //set eg fit parameter lower bounds
  map<string, float> egInvMassFitParLowerBounds;
  egInvMassFitParLowerBounds["m"] = 71.0; //mass variable (no fixed value, only range)
  egInvMassFitParLowerBounds["muCB"] = 86.2;
  egInvMassFitParLowerBounds["sigma"] = 1.0;
  egInvMassFitParLowerBounds["alphaCB"] = 1.063;
  egInvMassFitParLowerBounds["n"] = 143.16;
  egInvMassFitParLowerBounds["alphaCMS"] = 72.02;
  egInvMassFitParLowerBounds["beta"] = 0.098;
  egInvMassFitParLowerBounds["gamma"] = 0.0375;
  egInvMassFitParLowerBounds["muCMS"] = 56.0;

  //set eg fit parameter upper bounds
  map<string, float> egInvMassFitParUpperBounds;
  egInvMassFitParUpperBounds["m"] = 111.0;
  egInvMassFitParUpperBounds["muCB"] = 96.2;
  egInvMassFitParUpperBounds["sigma"] = 5.0;
  egInvMassFitParUpperBounds["alphaCB"] = 1.063;
  egInvMassFitParUpperBounds["n"] = 143.16;
  egInvMassFitParUpperBounds["alphaCMS"] = 72.02;
  egInvMassFitParUpperBounds["beta"] = 0.098;
  egInvMassFitParUpperBounds["gamma"] = 0.0375;
  egInvMassFitParUpperBounds["muCMS"] = 56.0;

  //set fit parameter units
  map<string, string> diEMInvMassFitParUnits;
  diEMInvMassFitParUnits["m"] = "GeV";
  diEMInvMassFitParUnits["muCB"] = "GeV";
  diEMInvMassFitParUnits["sigma"] = "GeV";
  diEMInvMassFitParUnits["alphaCB"] = "";
  diEMInvMassFitParUnits["n"] = "";
  diEMInvMassFitParUnits["alphaCMS"] = "GeV";
  diEMInvMassFitParUnits["beta"] = "";
  diEMInvMassFitParUnits["gamma"] = "";
  diEMInvMassFitParUnits["muCMS"] = "GeV";

  //get the histograms
  vector<TH3F*> eeMETVsDiEMETVsInvMass(3, NULL);
  vector<TH3F*> egMETVsDiEMETVsInvMass(1, NULL);
  TFile in(inFile);
  if (!in.IsOpen()) {
    cerr << "Error opening file " << inFile << ".\n";
    return;
  }
  in.GetObject("eeMETVsDiEMETVsInvMassTot", eeMETVsDiEMETVsInvMass[0]);
  in.GetObject("eeLowSidebandMETVsDiEMETVsInvMassTot", eeMETVsDiEMETVsInvMass[1]);
  in.GetObject("eeHighSidebandMETVsDiEMETVsInvMassTot", eeMETVsDiEMETVsInvMass[2]);
  in.GetObject("egMETVsDiEMETVsInvMass", egMETVsDiEMETVsInvMass[0]);
  if (eeMETVsDiEMETVsInvMass[0] == NULL) {
    cerr << "Error getting TH3F eeMETVsDiEMETVsInvMassTot from file " << inFile << ".\n";
    in.Close();
    return;
  }
  if (eeMETVsDiEMETVsInvMass[1] == NULL) {
    cerr << "Error getting TH3F eeLowSidebandMETVsDiEMETVsInvMassTot from file " << inFile << ".\n";
    in.Close();
    return;
  }
  if (eeMETVsDiEMETVsInvMass[2] == NULL) {
    cerr << "Error getting TH3F eeHighSidebandMETVsDiEMETVsInvMassTot from file " << inFile;
    cerr << ".\n";
    in.Close();
    return;
  }
  if (egMETVsDiEMETVsInvMass[0] == NULL) {
    cerr << "Error getting TH3F egMETVsDiEMETVsInvMass from file " << inFile << ".\n";
    in.Close();
    return;
  }

  //open the output file
  TFile out(outFile, "RECREATE");
  if (!out.IsOpen()) {
    cerr << "Error opening file " << outFile << ".\n";
    in.Close();
    return;
  }

  //run
  float egMisIDRateErr = -1.0;
  float egMisIDRate = electronPhotonMisIDRate(eeMETVsDiEMETVsInvMass, egMETVsDiEMETVsInvMass, 
					      eeInvMassFitParLowerBounds, 
					      eeInvMassFitParUpperBounds, eeInvMassFitPars, 
					      egInvMassFitParLowerBounds, 
					      egInvMassFitParUpperBounds, egInvMassFitPars, 
					      diEMInvMassFitParUnits, out, egMisIDRateErr);
  cerr << "e-->g mis-ID rate: " << egMisIDRate << " +/- " << egMisIDRateErr << endl;

  //close
  out.Write();
  out.Close();
  in.Close();
}

void printFileOpenErrMess(const string& file) {
  cerr << "Error opening file " << file << ".  Quitting.\n";
}

void printNullPtrErrMess(const string& ptr, const string& src) {
  cerr << "Error getting " << ptr << " from " << src << ".  Quitting.\n";
}

void plotDependentEGMisIDRate(const string& inFile, const string& outFile, 
			      const vector<string>& vars, const vector<string>& formattedVars)
{
  //open input file
  TFile in(inFile.c_str());
  if (!in.IsOpen()) {
    printFileOpenErrMess(inFile);
    return;
  }

  //open output file
  TFile out(outFile.c_str(), "RECREATE");
  if (!out.IsOpen()) {
    printFileOpenErrMess(outFile);
    in.Close();
    return;
  }

  //loop over variables
  for (vector<string>::const_iterator iVar = vars.begin(); iVar != vars.end(); ++iVar) {

    //get canvas from file
    TCanvas* inEGMisIDRateVsVarCanvas;
    in.GetObject(("egMisIDRateVs" + *iVar).c_str(), inEGMisIDRateVsVarCanvas);
    if (inEGMisIDRateVsVarCanvas == NULL) {
      printNullPtrErrMess("egMisIDRateVs" + *iVar, inFile);
      in.Close();
      return;
    }

    //get graph from canvas
    TGraph* inEGMisIDRateVsVarGraph = (TGraph*)inEGMisIDRateVsVarCanvas->GetPrimitive("Graph");
    if (inEGMisIDRateVsVarGraph == NULL) {
      printNullPtrErrMess("Graph", "egMisIDRateVs" + *iVar);
      in.Close();
      return;
    }

    //make new graph
    int nBins = inEGMisIDRateVsVarGraph->GetN();
    if (nBins <= 0) {
      cerr << "Error: Graph has " << nBins << " points.  Quitting.\n";
      in.Close();
      return;
    }
    vector<double> binCenters;
    vector<double> binHalfWidths;
    vector<double> egMisIDRates;
    vector<double> egMisIDRateErrs;
    if (*iVar == "Eta") {
      for (int iPoint = 1; iPoint <= nBins; ++iPoint) {
	double binLowEdge = 0.0;
	double binHighEdge = 0.0;
	double egMisIDRate = -1.0;
	double dummy = -1.0;
	inEGMisIDRateVsVarGraph->GetPoint(iPoint - 1, binLowEdge, egMisIDRate);
	inEGMisIDRateVsVarGraph->GetPoint(iPoint, binHighEdge, dummy);
	if (iPoint == nBins) binHighEdge = 1.4442;
	binCenters.push_back((binHighEdge + binLowEdge)/2.0);
	binHalfWidths.push_back((binHighEdge - binLowEdge)/2.0);
	egMisIDRates.push_back(egMisIDRate);
	egMisIDRateErrs.push_back(inEGMisIDRateVsVarGraph->GetErrorY(iPoint - 1));
      }
    }
    TGraphErrors outEGMisIDRateVsVarGraph(nBins, &binCenters[0], &egMisIDRates[0], 
					  &binHalfWidths[0], &egMisIDRateErrs[0]);
    if (*iVar == "Eta") {
      setGraphOptions(dynamic_cast<TGraph*>(&outEGMisIDRateVsVarGraph), kBlack, 0.7, 20, 
		      formattedVars[iVar - vars.begin()].c_str(), "f_{e#rightarrow#gamma}");
    }
    else {
      setGraphOptions(dynamic_cast<TGraph*>(inEGMisIDRateVsVarGraph), kBlack, 0.7, 20, 
		      formattedVars[iVar - vars.begin()].c_str(), "f_{e#rightarrow#gamma}");
    }

    //make new canvas
    TCanvas outEGMisIDRateVsVarCanvas(("egMisIDRateVs" + *iVar + "Canvas").c_str(), "", 600, 600);
    setCanvasOptions(&outEGMisIDRateVsVarCanvas, 1, 0, 0);

    //make legend
    TLegend egMisIDRateVsVarLegend(0.6, 0.2, 0.9, 0.3);
    setLegendOptions(&egMisIDRateVsVarLegend, "CMS 4.7 fb^{-1}, no jet requirement");

    //draw new graph on new canvas
    if (*iVar == "Eta") outEGMisIDRateVsVarGraph.Draw("AP");
    else inEGMisIDRateVsVarGraph->Draw("AP");
    egMisIDRateVsVarLegend.Draw();

    //write new canvas
    out.cd();
    outEGMisIDRateVsVarCanvas.Write();
    // if (*iVar == "Eta") {
    //   outEGMisIDRateVsVarCanvas.
    // 	SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/feg_vs_eta.pdf");
    // }
    // else {
    //   outEGMisIDRateVsVarCanvas.
    // 	SaveAs("/Users/rachelyohay/Documents/UVa/dissertation/feg_vs_pT.pdf");
    // }
  }

  //close all files
  out.Close();
  in.Close();
}

void estimateJESUncertaintyOnBkg(const vector<string>& inputFileNames, const vector<float>& METBins, const string& outputFileName)
{
  //find lowest MET bin to define offset
  float minMET = METBins[0];
  for (vector<float>::const_iterator iMETBin = METBins.begin(); iMETBin != METBins.end(); ++iMETBin) {
    if (*iMETBin < minMET) minMET = *iMETBin;
  }

  //loop over JES uncertainty toys
  vector<vector<float> > nEEBkg(METBins.size(), vector<float>());
  vector<vector<float> > nFFBkg(METBins.size(), vector<float>());
  for (vector<string>::const_iterator iInputFileName = inputFileNames.begin(); iInputFileName != inputFileNames.end(); ++iInputFileName) {
    TFile file(iInputFileName->c_str());
    if (!file.IsOpen()) {
      cerr << "Error opening file " << *iInputFileName << ". Quitting.\n";
      return;
    }

    //get bkg. estimates
    TH1F* eeFinal = NULL;
    TH1F* ffFinal = NULL;
    file.GetObject("eeFinal", eeFinal);
    file.GetObject("ffFinal", ffFinal);
    if ((eeFinal == NULL) || (ffFinal == NULL)) {
      cerr << "Error: eeFinal = " << eeFinal << " and ffFinal = " << ffFinal << " in file " << *iInputFileName << ".  Quitting.\n";
      file.Close();
      return;
    }

    //determine offset
    const Int_t eeOffset = eeFinal->FindBin(minMET);
    const Int_t ffOffset = ffFinal->FindBin(minMET);

    //fill vector of toy bkg. estimates for each MET bin
    for (vector<float>::const_iterator iMETBin = METBins.begin(); iMETBin != METBins.end(); ++iMETBin) {
      Int_t eeBin = eeFinal->FindBin(*iMETBin);
      Int_t ffBin = ffFinal->FindBin(*iMETBin);
      nEEBkg[eeBin - eeOffset].push_back(eeFinal->GetBinContent(eeBin));
      nFFBkg[ffBin - ffOffset].push_back(ffFinal->GetBinContent(ffBin));
    }

    //close file
    file.Close();
  }

  //open file to save fits
  TFile out(outputFileName.c_str(), "RECREATE");
  if (!out.IsOpen()) {
    cerr << "Error opening file " << outputFileName << ".  Quitting.\n";
    return;
  }

  //loop over MET bins
  for (vector<float>::const_iterator iMETBin = METBins.begin(); iMETBin != METBins.end(); ++iMETBin) {
    const unsigned int i = iMETBin - METBins.begin();

    //plot distribution of bkg. estimates from toys
    stringstream eeHistName;
    stringstream ffHistName;
    eeHistName << "eeBkgToyDist_METBin" << i;
    ffHistName << "ffBkgToyDist_METBin" << i;
    vector<float> eeBkgToys(nEEBkg[i]);
    vector<float> ffBkgToys(nFFBkg[i]);
    sort(eeBkgToys.begin(), eeBkgToys.end());
    sort(ffBkgToys.begin(), ffBkgToys.end());
    const float eeStep = (*(eeBkgToys.end() - 1) - *(eeBkgToys.begin()))/25.0;
    const float ffStep = (*(ffBkgToys.end() - 1) - *(ffBkgToys.begin()))/25.0;
    TH1F eeBkgToyDist(eeHistName.str().c_str(), "", 26, *(eeBkgToys.begin()), *(eeBkgToys.end() - 1) + eeStep);
    TH1F ffBkgToyDist(ffHistName.str().c_str(), "", 26, *(ffBkgToys.begin()), *(ffBkgToys.end() - 1) + ffStep);
    for (vector<float>::const_iterator iToy = eeBkgToys.begin(); iToy != eeBkgToys.end(); ++iToy) { eeBkgToyDist.Fill(*iToy); }
    for (vector<float>::const_iterator iToy = ffBkgToys.begin(); iToy != ffBkgToys.end(); ++iToy) { ffBkgToyDist.Fill(*iToy); }

    //fit with Gaussian to get relative error
//     //sometimes the Gaussian peak is not in the center of the axis, so constrain the fit to equal amounts around the maximum bin
//     const Int_t eeMaxBin = eeBkgToyDist.GetMaximumBin();
//     const Int_t ffMaxBin = ffBkgToyDist.GetMaximumBin();
//     const float eeRangeLeftToMax = eeBkgToyDist.GetBinLowEdge(eeMaxBin) - eeBkgToyDist.GetBinLowEdge(1);
//     const float ffRangeLeftToMax = ffBkgToyDist.GetBinLowEdge(ffMaxBin) - ffBkgToyDist.GetBinLowEdge(1);
//     const float eeRangeRightToMax = 
//       eeBkgToyDist.GetBinLowEdge(eeBkgToyDist.GetNbinsX()) - (eeBkgToyDist.GetBinLowEdge(eeMaxBin) + eeBkgToyDist.GetBinWidth(eeMaxBin));
//     const float ffRangeRightToMax = 
//       ffBkgToyDist.GetBinLowEdge(ffBkgToyDist.GetNbinsX()) - (ffBkgToyDist.GetBinLowEdge(ffMaxBin) + ffBkgToyDist.GetBinWidth(ffMaxBin));
//     float eeRangeMin = 0.0;
//     float eeRangeMax = 0.0;
//     float ffRangeMin = 0.0;
//     float ffRangeMax = 0.0;
//     if (eeRangeRightToMax > eeRangeLeftToMax) {
//       eeRangeMin = eeBkgToyDist.GetBinLowEdge(1);
//       eeRangeMax = eeBkgToyDist.GetBinLowEdge(eeMaxBin);
//     }
//     else {
//       eeRangeMax = eeBkgToyDist.GetBinLowEdge(eeMaxBin) + eeBkgToyDist.GetBinWidth(eeMaxBin);
//       eeRangeMin = eeBkgToyDist.GetBinLowEdge(eeMaxBin);
//     }
//     if (ffRangeRightToMax > ffRangeLeftToMax) {
//       ffRangeMin = ffBkgToyDist.GetBinLowEdge(1);
//       ffRangeMax = ffBkgToyDist.GetBinLowEdge(ffMaxBin);
//     }
//     else {
//       ffRangeMax = ffBkgToyDist.GetBinLowEdge(ffMaxBin) + ffBkgToyDist.GetBinWidth(ffMaxBin);
//       ffRangeMin = ffBkgToyDist.GetBinLowEdge(ffMaxBin);
//     }
//     TF1* eeFit = new TF1("eeFit", "gaus", eeRangeMin, eeRangeMax);
//     TF1* ffFit = new TF1("ffFit", "gaus", ffRangeMin, ffRangeMax);
    eeBkgToyDist.Fit("gaus", "LQM");
    ffBkgToyDist.Fit("gaus", "LQM");
//     eeBkgToyDist.Fit(eeFit, "LQMR");
//     ffBkgToyDist.Fit(ffFit, "LQMR");
    TF1* eeFit = eeBkgToyDist.GetFunction("gaus");
    TF1* ffFit = ffBkgToyDist.GetFunction("gaus");
    if ((eeFit == NULL) || (ffFit == NULL)) {
      cerr << "Error: eeFit = " << eeFit << " and ffFit = " << ffFit << ".  Quitting.\n";
      return;
    }
    if (i < (METBins.size() - 1)) cout << *iMETBin << " GeV <= MET < " << *(iMETBin + 1);
    else cout << "MET >= " << *iMETBin;
    cout << " GeV\n";
    cout << "\tRel. JES uncertainty on ee bkg.: " << eeFit->GetParameter(2)/eeFit->GetParameter(1) << endl;
    //     cout << "\tChi2/ndof for ee toy fit: " << eeFit->GetChisquare()/eeFit->GetNDF() << endl;
    cout << "\tRel. JES uncertainty on ff bkg.: " << ffFit->GetParameter(2)/ffFit->GetParameter(1) << endl;
    //     cout << "\tChi2/ndof for ff toy fit: " << ffFit->GetChisquare()/ffFit->GetNDF() << endl;
    cout << "\tAbs. JES uncertainty on ee bkg.: " << eeFit->GetParameter(2) << endl;
    cout << "\tAbs. JES uncertainty on ff bkg.: " << ffFit->GetParameter(2) << endl;

    //write
    eeBkgToyDist.Write();
    ffBkgToyDist.Write();
    delete eeFit;
    delete ffFit;
  }

  //close
  out.Write();
  out.Close();
}
