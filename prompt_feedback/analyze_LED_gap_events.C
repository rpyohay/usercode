#include <sys/time.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <ctime>
#include <sstream>
#include "TChain.h"
#include "TRFIOFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"
#include "TLegend.h"
#include "TChainElement.h"

//constants
#define NUM_CRYSTALS 14648
#define LOW_MASK 0xFFFFFFFF
#define IX_MIN 1
#define IX_MAX 100
#define IY_MIN 1
#define IY_MAX 100
#define NUM_EE_DCCS 18
#define NUM_LM_REGIONS 20
#define FIRST_LM_REGION 73
#define NUM_LED_BOXES 16
#define MIN_DEE 1
#define MAX_DEE 4
#define NUM_CHARS_IN_LM_MAP_FILE_NAME 29
#define NUM_DEES 4
#define DAC_ADDR_INTERVAL 2
#define LOWEST_DAC_ADDR 0xA8

using namespace std;

typedef unsigned long long TimeValue_t;

//load the diffusing sphere map in memory
const int loadMapping(unsigned int diffusingSpheres[], unsigned int dees[], const unsigned int size)
{
  //check array integrity
  if (size != NUM_CRYSTALS) return -1;

  //open the mapping file
  ifstream mapStream("diffusing_sphere_index.txt");
  if (!mapStream.good()) {
    cerr << "Error opening file diffusing_sphere_index.txt.\n";
    return -1;
  }

  //loop over the file
  for (unsigned int i = 0; i < NUM_CRYSTALS; ++i) { diffusingSpheres[i] = 0; }
  for (unsigned int i = 0; i < NUM_CRYSTALS; ++i) { dees[i] = 0; }
  while (!mapStream.eof()) {
    unsigned int dee, diffusingSphere, ix, iy, hashedIndex;
    int iz;
    mapStream >> dee >> diffusingSphere >> ix >> iy >> iz >> hashedIndex;
    mapStream.seekg(1, ios::cur);
    diffusingSpheres[hashedIndex] = diffusingSphere;
    dees[hashedIndex] = dee;
  }

  //success
  return 0;
}

//get the DAC address of a given diffusing sphere in a given dee
const int getDACAddr(const unsigned int diffusingSphere, const unsigned int dee)
{
  const int ERROR = -1;
  int retVal;
  switch (dee) {
  case 1:
    if ((diffusingSphere == 10) || ((diffusingSphere >= 12) && (diffusingSphere <= 14)) || (diffusingSphere == 18) || 
	(diffusingSphere == 19)) retVal = 0xA8;
    else if ((diffusingSphere >= 15) && (diffusingSphere <= 17)) retVal = 0xAA;
    else if ((diffusingSphere == 2) || (diffusingSphere == 4) || ((diffusingSphere >= 6) && (diffusingSphere <= 9))) retVal = 0xAC;
    else if ((diffusingSphere == 1) || (diffusingSphere == 3) || (diffusingSphere == 5) || (diffusingSphere == 11)) retVal = 0xAE;
    else retVal = ERROR;
    break;
  case 2:
    if (((diffusingSphere >= 1) && (diffusingSphere <= 5)) || (diffusingSphere == 8)) retVal = 0xA8;
    else if ((diffusingSphere == 10) || (diffusingSphere == 12) || (diffusingSphere == 14)) retVal = 0xAA;
    else if ((diffusingSphere == 13) || (diffusingSphere == 15) || ((diffusingSphere >= 16) && (diffusingSphere <= 19))) retVal = 0xAC;
    else if ((diffusingSphere == 6) || (diffusingSphere == 7) || (diffusingSphere == 9) || (diffusingSphere == 11)) retVal = 0xAE;
    else retVal = ERROR;
    break;
  case 3:
    if ((diffusingSphere >= 1) && (diffusingSphere <= 6)) retVal = 0xA8;
    else if ((diffusingSphere == 10) || (diffusingSphere == 13) || (diffusingSphere == 14)) retVal = 0xAA;
    else if ((diffusingSphere == 12) || ((diffusingSphere >= 15) && (diffusingSphere <= 19))) retVal = 0xAC;
    else if (((diffusingSphere >= 7) && (diffusingSphere <= 9)) || (diffusingSphere == 11)) retVal = 0xAE;
    else retVal = ERROR;
    break;
  case 4:
    if ((diffusingSphere == 13) || ((diffusingSphere >= 15) && (diffusingSphere <= 19))) retVal = 0xA8;
    else if ((diffusingSphere == 10) || (diffusingSphere == 12) || (diffusingSphere == 14)) retVal = 0xAA;
    else if (((diffusingSphere >= 1) && (diffusingSphere <= 5)) || (diffusingSphere == 8)) retVal = 0xAC;
    else if ((diffusingSphere == 6) || (diffusingSphere == 7) || (diffusingSphere == 9) || (diffusingSphere == 11)) retVal = 0xAE;
    else retVal = ERROR;
    break;
  default:
    retVal = ERROR;
    break;
  }
  return retVal;
}

void analyzeLEDGapEvents(vector<string> sourceFiles, const char* outputFile, time_t start)
{
  //constants
  const int EEM = -1;
  const int EEP = 1;
  //const unsigned int colors[NUM_LED_BOXES] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 28, 36, 40, 42, 46, 49}; //bright ROOT colors

  //create a chain of the input trees
  TChain* sourceTrees = new TChain("tree_");
  for (vector<string>::const_iterator iFile = sourceFiles.begin(); iFile != sourceFiles.end(); ++iFile) {
    sourceTrees->Add((*iFile).c_str());
  }
  const Int_t numPoints = sourceTrees->GetEntries();

  //declare/initialize variables and book histograms
  int ix, iy, iz, hashedIndex, DCC;
  unsigned int ixiy[NUM_CRYSTALS];
  unsigned int diffusingSpheres[NUM_CRYSTALS];
  unsigned int dees[NUM_CRYSTALS];
  int izArray[NUM_CRYSTALS];
  int evtCounter[NUM_CRYSTALS];
  int DCCs[NUM_CRYSTALS];
  float A, t;
  float A0[NUM_CRYSTALS];
  double t0[NUM_CRYSTALS];
  double tEnd = 0.0;
  double tBegin;
  TimeValue_t timestamp;
  TCanvas* timingCanvas = new TCanvas("timingCanvas", "Time of max sample (BX)", 600, 600);
  timingCanvas->SetWindowSize(1200 - timingCanvas->GetWw(), 1200 - timingCanvas->GetWh());
  timingCanvas->SetFillStyle(4000);
  timingCanvas->SetFillColor(0); //white
  TH1F* timingDist = new TH1F("timingDist", "Time of max sample (BX)", 1100, -0.5, 10.5);
  timingDist->GetXaxis()->SetTitle("Time of max sample (BX)");
  timingDist->GetXaxis()->CenterTitle();
  TCanvas* avgTimingCanvas = new TCanvas("avgTimingCanvas", "Average timing per crystal (BX)", 800, 400);
  avgTimingCanvas->SetWindowSize(1600 - avgTimingCanvas->GetWw(), 800 - avgTimingCanvas->GetWh());
  avgTimingCanvas->SetFillStyle(4000);
  avgTimingCanvas->SetFillColor(0); //white
  avgTimingCanvas->Divide(2, 1);
  TH2F* avgTimingEEM = new TH2F("avgTimingEEM", "Average timing per crystal (BX) for EE-", 100, 0.5, 100.5, 100, 0.5, 100.5);
  avgTimingEEM->GetXaxis()->SetTitle("ix");
  avgTimingEEM->GetXaxis()->CenterTitle();
  avgTimingEEM->GetYaxis()->SetTitle("iy");
  avgTimingEEM->GetYaxis()->CenterTitle();
  TH2F* avgTimingEEP = new TH2F("avgTimingEEP", "Average timing per crystal (BX) for EE+", 100, 0.5, 100.5, 100, 0.5, 100.5);
  avgTimingEEP->GetXaxis()->SetTitle("ix");
  avgTimingEEP->GetXaxis()->CenterTitle();
  avgTimingEEP->GetYaxis()->SetTitle("iy");
  avgTimingEEP->GetYaxis()->CenterTitle();
  TCanvas* RMSTimingCanvas = new TCanvas("RMSTimingCanvas", "RMS of the average timing per crystal (BX)", 800, 400);
  RMSTimingCanvas->SetWindowSize(1600 - RMSTimingCanvas->GetWw(), 800 - RMSTimingCanvas->GetWh());
  RMSTimingCanvas->SetFillStyle(4000);
  RMSTimingCanvas->SetFillColor(0); //white
  RMSTimingCanvas->Divide(2, 1);
  TH2F* RMSTimingEEM = new TH2F("RMSTimingEEM", "RMS of the average timing per crystal (BX) for EE-", 100, 0.5, 100.5, 100, 0.5, 100.5);
  RMSTimingEEM->GetXaxis()->SetTitle("ix");
  RMSTimingEEM->GetXaxis()->CenterTitle();
  RMSTimingEEM->GetYaxis()->SetTitle("iy");
  RMSTimingEEM->GetYaxis()->CenterTitle();
  RMSTimingEEM->GetZaxis()->SetLabelSize(0.02);
  RMSTimingEEM->GetZaxis()->SetLabelOffset(0.001);
  TH2F* RMSTimingEEP = new TH2F("RMSTimingEEP", "RMS of the average timing per crystal (BX) for EE+", 100, 0.5, 100.5, 100, 0.5, 100.5);
  RMSTimingEEP->GetXaxis()->SetTitle("ix");
  RMSTimingEEP->GetXaxis()->CenterTitle();
  RMSTimingEEP->GetYaxis()->SetTitle("iy");
  RMSTimingEEP->GetYaxis()->CenterTitle();
  RMSTimingEEP->GetZaxis()->SetLabelSize(0.02);
  RMSTimingEEP->GetZaxis()->SetLabelOffset(0.001);
  TCanvas* lowtMaxOccupancyCanvas = new TCanvas("lowtMaxOccupancyCanvas", "Occupancy of pulses with t_{Max} < 0.06", 800, 400);
  lowtMaxOccupancyCanvas->SetWindowSize(1600 - lowtMaxOccupancyCanvas->GetWw(), 800 - lowtMaxOccupancyCanvas->GetWh());
  lowtMaxOccupancyCanvas->SetFillStyle(4000);
  lowtMaxOccupancyCanvas->SetFillColor(0); //white
  lowtMaxOccupancyCanvas->Divide(2, 1);
  TH2F* lowtMaxOccupancyEEM = new TH2F("lowtMaxOccupancyEEM", "Occupancy of pulses with t_{Max} < 1 in EE-", 100, 0.5, 100.5, 100, 0.5, 100.5);
  lowtMaxOccupancyEEM->GetXaxis()->SetTitle("ix");
  lowtMaxOccupancyEEM->GetXaxis()->CenterTitle();
  lowtMaxOccupancyEEM->GetYaxis()->SetTitle("iy");
  lowtMaxOccupancyEEM->GetYaxis()->CenterTitle();
  TH2F* lowtMaxOccupancyEEP = new TH2F("lowtMaxOccupancyEEP", "Occupancy of pulses with t_{Max} < 1 in EE+", 100, 0.5, 100.5, 100, 0.5, 100.5);
  lowtMaxOccupancyEEP->GetXaxis()->SetTitle("ix");
  lowtMaxOccupancyEEP->GetXaxis()->CenterTitle();
  lowtMaxOccupancyEEP->GetYaxis()->SetTitle("iy");
  lowtMaxOccupancyEEP->GetYaxis()->CenterTitle();
  TH1F* timingDistsPerCrystal[NUM_CRYSTALS];
  TCanvas* avgAmpCanvas = new TCanvas("avgAmpCanvas", "Average amplitude per crystal (ADC counts)", 800, 400);
  avgAmpCanvas->SetWindowSize(1600 - avgAmpCanvas->GetWw(), 800 - avgAmpCanvas->GetWh());
  avgAmpCanvas->SetFillStyle(4000);
  avgAmpCanvas->SetFillColor(0); //white
  avgAmpCanvas->Divide(2, 1);
  TH2F* avgAmpEEM = new TH2F("avgAmpEEM", "Average amplitude per crystal (ADC counts) for EE-", 100, 0.5, 100.5, 100, 0.5, 100.5);
  avgAmpEEM->GetXaxis()->SetTitle("ix");
  avgAmpEEM->GetXaxis()->CenterTitle();
  avgAmpEEM->GetYaxis()->SetTitle("iy");
  avgAmpEEM->GetYaxis()->CenterTitle();
  avgAmpEEM->GetZaxis()->SetLabelSize(0.03);
  TH2F* avgAmpEEP = new TH2F("avgAmpEEP", "Average amplitude per crystal (ADC counts) for EE+", 100, 0.5, 100.5, 100, 0.5, 100.5);
  avgAmpEEP->GetXaxis()->SetTitle("ix");
  avgAmpEEP->GetXaxis()->CenterTitle();
  avgAmpEEP->GetYaxis()->SetTitle("iy");
  avgAmpEEP->GetYaxis()->CenterTitle();
  avgAmpEEP->GetZaxis()->SetLabelSize(0.03);
  TCanvas* timestampCanvas = new TCanvas("timestampCanvas", "Timestamp (seconds since epoch)", 600, 600);
  timestampCanvas->SetWindowSize(1200 - timestampCanvas->GetWw(), 1200 - timestampCanvas->GetWh());
  timestampCanvas->SetFillStyle(4000);
  timestampCanvas->SetFillColor(0); //white
  TCanvas* normACanvas = new TCanvas("normACanvas", "A/A_{0} vs. t - t_{0} for each crystal", 600, 600);
  normACanvas->SetWindowSize(1200 - normACanvas->GetWw(), 1200 - normACanvas->GetWh());
  normACanvas->SetFillStyle(4000);
  normACanvas->SetFillColor(0); //white
  TCanvas* AVsTimestampCanvas[NUM_LED_BOXES];
  sourceTrees->SetBranchStatus("*");
  sourceTrees->SetBranchAddress("ix_", &ix);
  sourceTrees->SetBranchAddress("iy_", &iy);
  sourceTrees->SetBranchAddress("iz_", &iz);
  sourceTrees->SetBranchAddress("hashedIndex_", &hashedIndex);
  sourceTrees->SetBranchAddress("DCC_", &DCC);
  sourceTrees->SetBranchAddress("A_", &A);
  sourceTrees->SetBranchAddress("t_", &t);
  sourceTrees->SetBranchAddress("timestamp_", &timestamp);
  for (unsigned int i = 0; i < NUM_CRYSTALS; ++i) {
    ixiy[i] = 0;
    izArray[i] = 0;
    evtCounter[i] = 0;
    A0[i] = 0.0;
    struct timeval currTime;
    gettimeofday(&currTime, 0);
    t0[i] = currTime.tv_sec + currTime.tv_usec*1e-6;
    tBegin = t0[i];
    timingDistsPerCrystal[i] = new TH1F("timingDistsPerCrystal", "timingDistsPerCrystal", 1000, -0.5, 9.5);
  }
  const int success = loadMapping(diffusingSpheres, dees, NUM_CRYSTALS);

  //loop over the tree
  for (Int_t i = 0; i < numPoints; ++i) {

    //get the entry
    sourceTrees->GetEvent(i);

    //plot the amplitude and timing sums (over all events) per crystal
    switch (iz) {
    case EEM:
      avgAmpEEM->Fill(ix, iy, A);
      avgTimingEEM->Fill(ix, iy, t);
      break;
    case EEP:
      avgAmpEEP->Fill(ix, iy, A);
      avgTimingEEP->Fill(ix, iy, t);
      break;
    default:
      cerr << "Error: iz = " << iz << ".\n";
      break;
    }

    //increment the crystal's event counter
    ++(evtCounter[hashedIndex]);

    //store ix, iy, and iz for this hashed index in ixiy[hashedIndex]
    ixiy[hashedIndex] = ((ix << 7) + iy);
    izArray[hashedIndex] = iz;

    //store the DCC for this crystal
    DCCs[hashedIndex] = DCC;

    //find the crystal's t0 and A0
    double timeSinceEpoch = (timestamp >> 32) + (timestamp & LOW_MASK)*1e-6;
    if (timeSinceEpoch < t0[hashedIndex]) {
      t0[hashedIndex] = timeSinceEpoch;
      A0[hashedIndex] = A;
    }

    //find the time of the last event
    if (timeSinceEpoch > tEnd) tEnd = timeSinceEpoch;
    if (timeSinceEpoch < tBegin) tBegin = timeSinceEpoch;

    //plot the timing distribution
    timingDist->Fill(t);
    timingDistsPerCrystal[hashedIndex]->Fill(t);

    //plot the occupancy of pulses with tMax < 1
    if (t < 1.0) {
      switch (iz) {
      case EEM:
	lowtMaxOccupancyEEM->Fill(ix, iy);
	break;
      case EEP:
	lowtMaxOccupancyEEP->Fill(ix, iy);
	break;
      default:
	cerr << "Error: invalid iz (iz = " << iz << ").\n";
	break;
      }
    }
  }
  double ceiling = ceil(tEnd);
  double theFloor = floor(tBegin);

  //plot the amplitude and timing average (over all events) per crystal
  for (unsigned int itx = IX_MIN; itx <= IX_MAX; ++itx) {
    for (unsigned int ity = IY_MIN; ity <= IY_MAX; ++ity) {
      float currBinContentEEM = avgAmpEEM->GetBinContent(itx, ity);
      float currBinContentTimingEEM = avgTimingEEM->GetBinContent(itx, ity);
      if (currBinContentEEM > 0) {
	for (int i = 0; i < NUM_CRYSTALS; ++i) {
	  if ((ixiy[i] == ((itx << 7) + ity)) && (izArray[i] == EEM)) avgAmpEEM->SetBinContent(itx, ity, currBinContentEEM/evtCounter[i]);
	}
      }
      if (currBinContentTimingEEM > 0) {
	for (int i = 0; i < NUM_CRYSTALS; ++i) {
	  if ((ixiy[i] == ((itx << 7) + ity)) && (izArray[i] == EEM)) avgTimingEEM->SetBinContent(itx, ity, currBinContentTimingEEM/evtCounter[i]);
	}
      }
      float currBinContentEEP = avgAmpEEP->GetBinContent(itx, ity);
      float currBinContentTimingEEP = avgTimingEEP->GetBinContent(itx, ity);
      if (currBinContentEEP > 0) {
	for (int i = 0; i < NUM_CRYSTALS; ++i) {
	  if ((ixiy[i] == ((itx << 7) + ity)) && (izArray[i] == EEP)) avgAmpEEP->SetBinContent(itx, ity, currBinContentEEP/evtCounter[i]);
	}
      }
      if (currBinContentTimingEEP > 0) {
	for (int i = 0; i < NUM_CRYSTALS; ++i) {
	  if ((ixiy[i] == ((itx << 7) + ity)) && (izArray[i] == EEP)) avgTimingEEP->SetBinContent(itx, ity, currBinContentTimingEEP/evtCounter[i]);
	}
      }
    }
  }

  //plot the timing RMS over all events per crystal
  for (int i = 0; i < NUM_CRYSTALS; ++i) {
    if (timingDistsPerCrystal[i]->GetEntries() > 0) {
      int z = izArray[i];
      switch (z) {
      case EEM:
	RMSTimingEEM->Fill(ixiy[i] >> 7, ixiy[i] & 0x7F, timingDistsPerCrystal[i]->GetRMS());
	break;
      case EEP:
	RMSTimingEEP->Fill(ixiy[i] >> 7, ixiy[i] & 0x7F, timingDistsPerCrystal[i]->GetRMS());
	break;
      default:
	cerr << "Incorrect z in izArray[" << i << "] -- z = " << z << ".\n";
	break;
      }
    }
    delete timingDistsPerCrystal[i];
  }

  //make pointers to timestamp histograms, one for each FED
  TH1D* timestampHists[NUM_EE_DCCS];
  TH1F* LMRegionHists[NUM_EE_DCCS];
  TCanvas* LMRegionCanvases[NUM_EE_DCCS];

  //initialize time of first event array
  double timeOfFirstEvt[NUM_EE_DCCS];
  for (int i = 0; i < NUM_EE_DCCS; ++i) {
    struct timeval currTime;
    gettimeofday(&currTime, 0);
    timeOfFirstEvt[i] = currTime.tv_sec + currTime.tv_usec*1e-6;
  }

  //construct a TLegend object
  TLegend* legend = new TLegend(0.6, 0.0, 0.9, 1.0);
  legend->SetHeader("DCC");

  //loop over FEDs
  for (unsigned int ipHist = 0; ipHist < NUM_EE_DCCS; ++ipHist) {

    //0-8: DCC 1-9 (EE-); 9-18: DCC 46-54 (EE+)
    int wantedDCC = ipHist + 1;
    if (ipHist >= NUM_EE_DCCS/2) wantedDCC+=36;

    //get the time of the first event in this DCC
    for (unsigned int i = 0; i < NUM_CRYSTALS; ++i) { if ((DCCs[i] == wantedDCC) && (t0[i] < timeOfFirstEvt[ipHist])) timeOfFirstEvt[ipHist] = t0[i]; }
    
    //book histograms
    char histName[100];
    sprintf(histName, "timestamp%i", wantedDCC);
    char histTitle[100];
    sprintf(histTitle, "Timestamp (seconds since epoch)");
    timestampHists[ipHist] = new TH1D(histName, histTitle, (Int_t)(ceiling - theFloor), theFloor, ceiling);
    timestampHists[ipHist]->GetXaxis()->SetTitle("Timestamp (seconds since epoch)");
    timestampHists[ipHist]->GetXaxis()->CenterTitle();
    timestampHists[ipHist]->GetXaxis()->SetLabelSize(0.02);
    timestampHists[ipHist]->GetYaxis()->SetLabelSize(0.03);
    timestampHists[ipHist]->SetLineColor(wantedDCC);
    timestampHists[ipHist]->SetFillColor(wantedDCC);
    timestampHists[ipHist]->SetFillStyle(1001);
    char legendEntry[2];
    sprintf(legendEntry, "%i", wantedDCC);
    legend->AddEntry(timestampHists[ipHist], legendEntry, "f");
    stringstream histNameStream;
    histNameStream << "LM_region_" << wantedDCC;
    stringstream histTitleStream;
    histTitleStream << "Event time, LM region " << wantedDCC;
    LMRegionHists[ipHist] = new TH1F(histNameStream.str().c_str(), histTitleStream.str().c_str(), 10800, -0.5, 10799.5); //3 hours, 1 sec/bin
    LMRegionHists[ipHist]->GetXaxis()->SetTitle("Time since start of log file (sec)");
    LMRegionHists[ipHist]->GetXaxis()->CenterTitle();
    LMRegionHists[ipHist]->GetXaxis()->SetLabelSize(0.02);
    LMRegionCanvases[ipHist] = new TCanvas(histNameStream.str().c_str(), histTitleStream.str().c_str(), 400, 400);
    LMRegionCanvases[ipHist]->SetWindowSize(800 - LMRegionCanvases[ipHist]->GetWw(), 800 - LMRegionCanvases[ipHist]->GetWh());
    LMRegionCanvases[ipHist]->SetFillStyle(4000);
    LMRegionCanvases[ipHist]->SetFillColor(0); //white
  }

  //book histogram for A/A0 vs. t - t0 for all crystals
  TH2D* normAVst = new TH2D("normAVst", "A/A_{0} vs. t - t_{0} for each crystal", Int_t(ceiling - theFloor), 0.0, ceiling - theFloor, 10000, 0.5, 1.5);
  normAVst->GetXaxis()->SetTitle("t - t_{0}");
  normAVst->GetXaxis()->CenterTitle();
  normAVst->GetYaxis()->SetTitle("A/A_{0}");
  normAVst->GetYaxis()->CenterTitle();
  normAVst->GetZaxis()->SetLabelSize(0.02);

  //make 16 histogram pointers for A vs. t - t0, 1 per LED box
  TH2D* AVsTimestamp[NUM_LED_BOXES];
  for (unsigned int ipHist = 0; ipHist < NUM_LED_BOXES; ++ipHist) {

    //mapping between LED box and ipHist:
    //ipHist = 0: dee 1, DAC address 0xA8
    //ipHist = 1: dee 1, DAC address 0xAA
    //ipHist = 2: dee 1, DAC address 0xAC
    //ipHist = 3: dee 1, DAC address 0xAE
    //ipHist = 4: dee 2, DAC address 0xA8
    //ipHist = 5: dee 2, DAC address 0xAA
    //ipHist = 6: dee 2, DAC address 0xAC
    //ipHist = 7: dee 2, DAC address 0xAE
    //ipHist = 8: dee 3, DAC address 0xA8
    //ipHist = 9: dee 3, DAC address 0xAA
    //ipHist = 10: dee 3, DAC address 0xAC
    //ipHist = 11: dee 3, DAC address 0xAE
    //ipHist = 12: dee 4, DAC address 0xA8
    //ipHist = 13: dee 4, DAC address 0xAA
    //ipHist = 14: dee 4, DAC address 0xAC
    //ipHist = 15: dee 4, DAC address 0xAE
    unsigned int dee = (ipHist/NUM_DEES) + 1;
    unsigned int DACAddr = (ipHist % NUM_DEES)*DAC_ADDR_INTERVAL + LOWEST_DAC_ADDR;
    char name[14];
    sprintf(name, "AVsTimestamp%i", ipHist);
    char title[59];
    sprintf(title, "A vs. t - t_{0} for each crystal in dee %i, DAC address 0x%X", dee, DACAddr);

    //create canvas
    AVsTimestampCanvas[ipHist] = new TCanvas(name, title, 600, 600);
    AVsTimestampCanvas[ipHist]->SetWindowSize(1200 - AVsTimestampCanvas[ipHist]->GetWw(), 1200 - AVsTimestampCanvas[ipHist]->GetWh());
    AVsTimestampCanvas[ipHist]->SetFillStyle(4000);
    AVsTimestampCanvas[ipHist]->SetFillColor(0); //white

    //book histogram
    AVsTimestamp[ipHist] = new TH2D(name, title, Int_t(ceiling - theFloor), 0.0, ceiling - theFloor, 1000, 0.0, 1000.0);
    AVsTimestamp[ipHist]->GetXaxis()->SetTitle("t - t_{0} (sec)");
    AVsTimestamp[ipHist]->GetXaxis()->CenterTitle();
    AVsTimestamp[ipHist]->GetXaxis()->SetTitleOffset(1.75);
    AVsTimestamp[ipHist]->GetYaxis()->SetTitle("A (ADC counts)");
    AVsTimestamp[ipHist]->GetYaxis()->CenterTitle();
    AVsTimestamp[ipHist]->GetYaxis()->SetTitleOffset(2.00);
    AVsTimestamp[ipHist]->GetYaxis()->SetLabelSize(0.03);
    AVsTimestamp[ipHist]->GetZaxis()->SetLabelSize(0.02);
  }

  //loop over the tree
  for (Int_t i = 0; i < numPoints; ++i) {

    //get the entry
    sourceTrees->GetEvent(i);

    //fill the timestamp histograms
    double timeSinceEpoch = (timestamp >> 32) + (timestamp & LOW_MASK)*1e-6;
    int iHist = DCC - 1;
    if (iHist >= NUM_EE_DCCS/2) iHist-=36;
    timestampHists[iHist]->Fill(timeSinceEpoch);
    LMRegionHists[iHist]->Fill(timeSinceEpoch - start);

    //plot amplitude (normalized to amplitude at time of first event) vs. time for each crystal
    //exclude crystals in diffusing sphere 13, dee 4 because we know it wasn't on
    /*if ((!((ix >= 21) && (ix <= 25) && (iy >= 36) && (iy <= 40))) && 
	(!((ix >= 11) && (ix <= 20) && (iy >= 31) && (iy <= 40))) && 
	(!((ix >= 11) && (ix <= 15) && (iy >= 41) && (iy <= 45)))) {*/
    /*normA.push_back(A/A0[hashedIndex]);
      tMinust0.push_back(timeSinceEpoch - t0[hashedIndex]);*/
    double displacedTime = timeSinceEpoch - t0[hashedIndex];
      normAVst->Fill(displacedTime, A/A0[hashedIndex]);
      //}

      //plot A vs. t - t0, divided into different histograms for the different LED boxes
      if (success == 0) {
	const unsigned int dee = dees[hashedIndex];
	const unsigned int diffusingSphere = diffusingSpheres[hashedIndex];
	const int DACAddr = getDACAddr(diffusingSphere, dee);
	AVsTimestamp[NUM_DEES*(dee - 1) + (DACAddr - LOWEST_DAC_ADDR)/DAC_ADDR_INTERVAL]->Fill(displacedTime, A);
      }
      else cerr << "Unable to fill histograms because mapping could not be loaded.\n";
  }

  //get timestamp histogram with largest bin content -- this one will be drawn first
  Double_t max = 0.0;
  int iMax = 0;
  for (int i = 0; i < NUM_EE_DCCS; ++i) {
    Double_t histMax = timestampHists[i]->GetMaximum();
    if (histMax > max) {
      max = histMax;
      iMax = i;
    }
  }

  //open output file
  TFile* output = new TFile(outputFile, "RECREATE");
  if (!output->IsOpen()) {
    cerr << "Error opening output file " << outputFile << ".\n";
    delete sourceTrees;
    delete output;
    return;
  }

  //draw and save all histograms
  output->cd();
  timingCanvas->cd();
  //if (timingDist->GetEntries() > 0) timingCanvas->SetLogx();
  timingDist->Draw();
  timingCanvas->Update();
  timingCanvas->Write();
  avgTimingCanvas->cd(1);
  //if (avgTimingEEM->GetEntries() > 0) avgTimingCanvas->cd(1)->SetLogz();
  avgTimingEEM->Draw("COLZ");
  avgTimingCanvas->cd(2);
  //if (avgTimingEEP->GetEntries() > 0) avgTimingCanvas->cd(2)->SetLogz();
  avgTimingEEP->Draw("COLZ");
  avgTimingCanvas->Update();
  avgTimingCanvas->Write();
  RMSTimingCanvas->cd(1);
  if (RMSTimingEEM->GetEntries() > 0) RMSTimingCanvas->cd(1)->SetLogz();
  RMSTimingEEM->Draw("COLZ");
  RMSTimingCanvas->cd(2);
  if (RMSTimingEEP->GetEntries() > 0) RMSTimingCanvas->cd(2)->SetLogz();
  RMSTimingEEP->Draw("COLZ");
  RMSTimingCanvas->Update();
  RMSTimingCanvas->Write();
  lowtMaxOccupancyCanvas->cd(1);
  if (lowtMaxOccupancyEEM->GetEntries() > 0) lowtMaxOccupancyCanvas->cd(1)->SetLogz();
  lowtMaxOccupancyEEM->Draw("COLZ");
  lowtMaxOccupancyCanvas->cd(2);
  if (lowtMaxOccupancyEEP->GetEntries() > 0) lowtMaxOccupancyCanvas->cd(2)->SetLogz();
  lowtMaxOccupancyEEP->Draw("COLZ");
  lowtMaxOccupancyCanvas->Update();
  lowtMaxOccupancyCanvas->Write();
  timestampCanvas->cd();
  timestampHists[iMax]->Draw();
  for (int i = 0; i < NUM_EE_DCCS; ++i) { if (i != iMax) timestampHists[i]->Draw("same"); }
  legend->Draw();
  timestampCanvas->Update();
  timestampCanvas->Write();
  for (int i = 0; i < NUM_EE_DCCS; ++i) {
    LMRegionCanvases[i]->cd();
    LMRegionHists[i]->Draw();
    LMRegionCanvases[i]->Update();
    LMRegionCanvases[i]->Write();
  }
  avgAmpCanvas->cd(1);
  if (avgAmpEEM->GetEntries() > 0) avgAmpCanvas->cd(1)->SetLogz();
  avgAmpEEM->Draw("COLZ");
  avgAmpCanvas->cd(2);
  if (avgAmpEEP->GetEntries() > 0) avgAmpCanvas->cd(2)->SetLogz();
  avgAmpEEP->Draw("COLZ");
  avgAmpCanvas->Update();
  avgAmpCanvas->Write();
  normACanvas->cd();
  if (normAVst->GetEntries() > 0) normACanvas->cd()->SetLogz();
  normAVst->Draw("COLZ");
  normACanvas->Update();
  normACanvas->Write();
  for (unsigned int ipHist = 0; ipHist < NUM_LED_BOXES; ++ipHist) {
    AVsTimestampCanvas[ipHist]->cd();
    AVsTimestamp[ipHist]->Draw("LEGO");
    AVsTimestampCanvas[ipHist]->Update();
    AVsTimestampCanvas[ipHist]->Write();
  }
  output->Write();
  output->Close();
  delete output;
  delete sourceTrees;
}
