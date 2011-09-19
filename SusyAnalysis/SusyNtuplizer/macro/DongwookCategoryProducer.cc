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
  photonTrackIsoMaxPTMultiplier_(iConfig.photonTrackIsoMaxPTMultiplier),
  photonTrackIsoMaxPTConstant_(iConfig.photonTrackIsoMaxPTConstant),
  photonSigmaIetaIetaMax_(iConfig.photonSigmaIetaIetaMax),
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
  outputFile_(iConfig.outputFile)
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
  treeReader_->SetUseTrigger(false);
  treeReader_->AddHltName(HLT_);
  treeReader_->SetFilter(false);
  treeReader_->SetProcessNEvents(nEvts_);
  treeReader_->IncludeAJson(JSON_);
  treeReader_->SetPhotonTag(photonTag_);
  treeReader_->SetPhotonECALIsoEffArea(photonECALIsoEffArea_);
  treeReader_->SetPhotonHCALIsoEffArea(photonHCALIsoEffArea_);

  //create the output tree and define new branches
  out.cd();
  outTree_ = treeReader_->fChain->CloneTree(0);
  // outTree_ = new TTree("susyTreeCategoryOnly", "");
  category_ = new susy::Category();
  outTree_->Branch("susyCategory", "susy::Category", &category_);

  //photon variables to pass to the Categorizer object
  VDOUBLE photonET;
  VDOUBLE photonEta;
  VDOUBLE photonECALIso;
  VDOUBLE photonHCALIso;
  VDOUBLE photonHOverE;
  VDOUBLE photonR9;
  VDOUBLE photonTrackIso;
  VDOUBLE photonSigmaIetaIeta;
  VDOUBLE photonSeedTime;
  VDOUBLE photonE2OverE9;
  VINT photonSeedIeta;
  VDOUBLE photonPhi;
  VBOOL photonHasPixelSeed;
  Categorizer categorizer(photonET, photonEta, photonECALIso, photonHCALIso, photonHOverE, 
			  photonR9, photonTrackIso, photonSigmaIetaIeta, photonSeedTime, 
			  photonE2OverE9, photonPhi, photonHasPixelSeed, photon1ETMin_, 
			  photon2ETMin_, photonAbsEtaMax_, photonECALIsoMaxPTMultiplier_, 
			  photonECALIsoMaxConstant_, photonHCALIsoMaxPTMultiplier_, 
			  photonHCALIsoMaxConstant_, photonHOverEMax_, photonR9Max_, 
			  photonTrackIsoMaxPTMultiplier_, photonTrackIsoMaxPTConstant_, 
			  photonSigmaIetaIetaMax_, photonAbsSeedTimeMax_, photonE2OverE9Max_, 
			  photonDPhiMin_, photonDRMin_, pixelVetoOnFake_);

  //run the category producer
  try { treeReader_->Loop(outTree_, categorizer, category_); }
  catch (STRING& badInput) { throw; }

  //write the output file
  outTree_->Write();
  out.Write();
  out.Close();

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
  delete treeReader_;
  delete category_;
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
