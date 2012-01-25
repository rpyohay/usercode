#include <memory>
#include "DongwookCategoryProducer.h"

DongwookCategoryProducer::DongwookCategoryProducer(const ParameterSet& iConfig) :
  photonTag_(iConfig.photonTag),
  photon1ETMin_(iConfig.photon1ETMin),
  photon2ETMin_(iConfig.photon2ETMin),
  photonAbsEtaMax_(iConfig.photonAbsEtaMax),
  photonECALIsoMaxPTMultiplier_(iConfig.photonECALIsoMaxPTMultiplier),
  photonECALIsoMaxConstant_(iConfig.photonECALIsoMaxConstant),
  photonECALIsoEffArea_(iConfig.photonECALIsoEffArea),
  photonHCALIsoMaxPTMultiplier_(iConfig.photonHCALIsoMaxPTMultiplier),
  photonHCALIsoMaxConstant_(iConfig.photonHCALIsoMaxConstant),
  photonHCALIsoEffArea_(iConfig.photonHCALIsoEffArea),
  photonHOverEMax_(iConfig.photonHOverEMax),
  photonR9Max_(iConfig.photonR9Max),
  photonR9Min_(iConfig.photonR9Min),
  photonTrackIsoMaxPTMultiplier_(iConfig.photonTrackIsoMaxPTMultiplier),
  photonTrackIsoMaxConstant_(iConfig.photonTrackIsoMaxConstant),
  photonCombinedIsoMax_(iConfig.photonCombinedIsoMax),
  fakeCombinedIsoMax_(iConfig.fakeCombinedIsoMax),
  isoConeHLT_(iConfig.isoConeHLT),
  isoConeOffline_(iConfig.isoConeOffline),
  photonSigmaIetaIetaMax_(iConfig.photonSigmaIetaIetaMax),
  photonHLTSigmaIetaIetaMax_(iConfig.photonHLTSigmaIetaIetaMax),
  photonAbsSeedTimeMax_(iConfig.photonAbsSeedTimeMax),
  photonE2OverE9Max_(iConfig.photonE2OverE9Max),
  photonDPhiMin_(iConfig.photonDPhiMin),
  photonDRMin_(iConfig.photonDRMin),
  pixelVetoOnFake_(iConfig.pixelVetoOnFake),
  treeName_(iConfig.treeName),
  input_(iConfig.input),
  HLT_(iConfig.HLT),
  nEvts_(iConfig.nEvts),
  JSON_(iConfig.JSON),
  outputFile_(iConfig.outputFile),
  recategorize_(iConfig.recategorize)
{
  //open output file (to hold tree with added branches)
  TFile out(outputFile_.c_str(), "RECREATE");

  //set up the tree reader
  chain_ = new TChain(treeName_.c_str());
  for (VSTRING_IT iIn = input_.begin(); iIn != input_.end(); ++iIn) {
    chain_->Add((*iIn).c_str());
  }
  treeReader_ = new EventAnalyzer(chain_);
  treeReader_->SetPrintInterval(10000);
  treeReader_->SetPrintLevel(0);
  treeReader_->SetUseTrigger(true);
  for (vector<TString>::const_iterator i = HLT_.begin(); i != HLT_.end(); ++i) {
    treeReader_->AddHltName(*i);
  }
  treeReader_->SetFilter(true);
  treeReader_->SetProcessNEvents(nEvts_);
  treeReader_->IncludeAJson(JSON_);
  treeReader_->SetPhotonTag(photonTag_);
  treeReader_->SetPhotonECALIsoEffArea(photonECALIsoEffArea_);
  treeReader_->SetPhotonHCALIsoEffArea(photonHCALIsoEffArea_);
  treeReader_->SetIsoConeHLT(isoConeHLT_);
  treeReader_->SetIsoConeOffline(isoConeOffline_);

  //create the output tree and define new branches
  out.cd();
  if (recategorize_) {
    outTree_ = new TTree("susyTree", "SUSY Event");
    event_ = new susy::Event();
    outTree_->Branch("susyEvent", "susy::Event", &event_);
  }
  else {
    outTree_ = treeReader_->fChain->CloneTree(0);
    event_ = NULL;
  }
  category_ = new susy::Category();
  outTree_->Branch("susyCategory", "susy::Category", &category_);

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
  Categorizer categorizer(photonET, photonEta, photonECALIso, PUSubtractedPhotonECALIso, 
			  photonHCALIso, PUSubtractedPhotonHCALIso, photonHOverE, photonR9, 
			  photonTrackIso, photonSigmaIetaIeta, photonSeedTime, photonE2OverE9, 
			  photonPhi, photonHasPixelSeed, photon1ETMin_, photon2ETMin_, 
			  photonAbsEtaMax_, photonECALIsoMaxPTMultiplier_, 
			  photonECALIsoMaxConstant_, photonHCALIsoMaxPTMultiplier_, 
			  photonHCALIsoMaxConstant_, photonHOverEMax_, photonR9Max_, photonR9Min_, 
			  photonTrackIsoMaxPTMultiplier_, photonTrackIsoMaxConstant_, 
			  photonCombinedIsoMax_, fakeCombinedIsoMax_, photonSigmaIetaIetaMax_, 
			  photonHLTSigmaIetaIetaMax_, photonAbsSeedTimeMax_, photonE2OverE9Max_, 
			  photonDPhiMin_, photonDRMin_, pixelVetoOnFake_);

  //run the category producer
  try { treeReader_->Loop(outTree_, categorizer, category_, event_); }
  catch (STRING& badInput) { throw; }

  //write the output file
  outTree_->Write();
  out.Write();
  out.Close();
  outTree_ = NULL;

  //write events in each sample to file
  STRING outDebugName(outputFile_.replace(outputFile_.find("root"), 4, "txt"));
  ofstream outDebug(outDebugName.c_str());
  if (outDebug.is_open()) {
    outDebug << "Run Event Lumi\n";
    writeEvents(treeReader_->ggEvts_, outDebug, "gg");
    writeEvents(treeReader_->egEvts_, outDebug, "eg");
    writeEvents(treeReader_->eeEvts_, outDebug, "ee");
    writeEvents(treeReader_->ffEvts_, outDebug, "ff");
    outDebug.close();
  }
  else std::cerr << "Error opening file " << outDebugName << ".\n";
}

DongwookCategoryProducer::~DongwookCategoryProducer()
{
  delete chain_;
  chain_ = NULL;
  delete treeReader_;
  treeReader_ = NULL;
  delete category_;
  category_ = NULL;
  if (event_ != NULL) {
    delete event_;
    event_ = NULL;
  }
}

void DongwookCategoryProducer::writeEvents(const RUNEVTLUMIMAP& evts, ofstream& out, 
					   const STRING& label) const
{
  out << label << ": " << evts.size() << " events\n";
  for (RUNEVTLUMIMAP::const_iterator iEvt = evts.begin(); iEvt != evts.end(); ++iEvt) {
    out << (*iEvt).first.first << " " << (*iEvt).first.second << " " << (*iEvt).second;
    out << std::endl;
  }
  out << "------------\n";
}

void DongwookCategoryProducer::writeEvents(const RUNEVTLUMIMASSMAP& evts, ofstream& out, 
					   const STRING& label) const
{
  out << label << ": " << evts.size() << " events\n";
  for (RUNEVTLUMIMASSMAP::const_iterator iEvt = evts.begin(); iEvt != evts.end(); ++iEvt) {
    out << (*iEvt).first.first.first << " " << (*iEvt).first.first.second << " " << (*iEvt).first.second << " " << (*iEvt).second;
    out << std::endl;
  }
  out << "------------\n";
}
