/*#ifndef SampleMaker_EventSelector_h
  #define SampleMaker_EventSelector_h*/
#ifndef EventSelector_h
#define EventSelector_h

// system include files
/*#include <memory>
  #include <sstream>*/
#include <fstream>
#include <iostream>
#include "format_old.h"

// user include files
/*#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

//reco::Photon classes
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

//ECAL RecHits
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

//for the general track collection
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

//for the HB/HE RecHit collection
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

//HCAL geometry
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

//muon geometry
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"*/

//ROOT
#include "TMath.h"

//sample types
#define GAMMAGAMMA 1
#define EGAMMA 2
#define EE 3
#define FF 4
#define ETRACK 5

//object types
#define HAS_PIXEL_SEED 0
#define LACKS_PIXEL_SEED 1

//ECAL fiducial regions
#define EB 1
#define EEND 2
#define ECAL 3

//muon endcap coverage
/*#define MUON_1_MIN_Z 5.7
#define MUON_1_MAX_Z 6.7
#define MUON_1_MAX_R 2.864
#define MUON_1_MAX_ETA 2.4
#define MUON_2_MIN_Z 6.7
#define MUON_2_MAX_Z 10.88
#define MUON_2_MIN_ETA 0.9
#define MUON_2_MAX_ETA 2.4*/

using namespace std;
/*using namespace reco;
  using namespace edm;*/

class EventSelector {

 public:

  //default constructor
  EventSelector();

  //constructor from a list of cut values
  EventSelector(const unsigned int, const double, const double, const double, const double, const double, const double, const double, const double, 
		const double, const double, const unsigned int, const bool, const double, const bool, const double, const bool, const double, const double, 
		const double, const double, const double, const bool, const double, const unsigned int, const unsigned int, const unsigned int, 
		const unsigned int, const string, const bool);

  //copy constructor
  EventSelector(/*const */EventSelector&);

  //destructor
  ~EventSelector();

  //assignment operator
  EventSelector& operator=(/*const */EventSelector&);

  //getters
  unsigned int getSampleType() const;
  double getECALIsoMaxPTMultiplierEB() const;
  double getECALIsoMaxConstantEB() const;
  double getECALIsoMaxPTMultiplierEE() const;
  double getECALIsoMaxConstantEE() const;
  double getHCALIsoMaxPTMultiplierEB() const;
  double getHCALIsoMaxConstantEB() const;
  double getHCALIsoMaxPTMultiplierEE() const;
  double getHCALIsoMaxConstantEE() const;
  double getHOverEMaxPresel() const;
  double getETMin() const;
  unsigned int getFiducialRegion() const;
  bool getUseHOverE() const;
  double getHOverEMax() const;
  bool getUseSigmaEtaEta() const;
  double getSigmaEtaEtaMax() const;
  bool getUseTrackIso() const;
  double getTrackIsoMaxPTMultiplier() const;
  double getTrackIsoMaxConstant() const;
  double getTrackPTMin() const;
  double getETrackRMin() const;
  double getMinDRPhotons() const;
  bool getUseTimingCut() const;
  double getMaxSeedTime() const;
  unsigned int getNumReqdCands() const;
  unsigned int getRun() const;
  unsigned int getEvt() const;
  unsigned int getLumiSec() const;
  string getDebugFileName() const;
  bool getDebugFlag() const;

  //setters
  void setRun(const unsigned int);
  void setEvt(const unsigned int);
  void setLumiSec(const unsigned int);

  //helper methods
  void printEvtInfo();
  void printDebug(const string);
  void initializeDebugging();
  bool debugFileOpen() const;
  void closeDebugFile();
  string sampleTypeString() const;
  unsigned int ECALFiducialRegion(/*const Photon**/const PhoInfoBranches&, const int) const;

  //calculate dPhi between two objects
  template <typename T>
  T dPhi(T phi1, T phi2) const
  {
    T dPhi = fabs(phi1 - phi2);
    if (dPhi > TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
    return dPhi;
  }

  //calculate dEta between two objects
  template <typename T>
  T dEta(T eta1, T eta2) const { return fabs(eta1 - eta2); }

  //calculate dR between two objects
  template <typename T>
  T dR(T eta1, T eta2, T phi1, T phi2) const { return sqrt(dEta(eta1, eta2)*dEta(eta1, eta2) + dPhi(phi1, phi2)*dPhi(phi1, phi2)); }

  //pass flags
  /*const bool passesDataQualityCuts(const Handle<HBHERecHitCollection>&, const map<unsigned int, Photon*>&, const CaloGeometry*, 
    const Handle<TrackCollection>&);*/
  bool passesDataQualityCuts(const HEHitInfoBranches&, const PhoInfoBranches&, const CosInfoBranches&, const vector<int>&);
  //const bool passesHEBeamHaloTag(const Handle<HBHERecHitCollection>&, const map<unsigned int, Photon*>&, const CaloGeometry*);
  bool passesHEBeamHaloTag(const HEHitInfoBranches&, const PhoInfoBranches&, const vector<int>&);
  //const bool photonIsHEHalo(const Handle<HBHERecHitCollection>&, const Photon*, const CaloGeometry*);
  bool photonIsHEHalo(const HEHitInfoBranches&, const PhoInfoBranches&, const int);
  //const bool passesMuonBeamHaloTag(const Handle<TrackCollection>&, const map<unsigned int, Photon*>&);
  bool passesMuonBeamHaloTag(const CosInfoBranches&, const PhoInfoBranches&, const vector<int>&);
  //const bool photonIsMuonHalo(const Handle<TrackCollection>&, const Photon*);
  bool photonIsMuonHalo(const CosInfoBranches&, const PhoInfoBranches&, const int);
  //const bool passesPreselection(const Photon*, const vector<EcalRecHit*>&, const unsigned int);
  bool passesPreselection(const PhoInfoBranches&, const int);
  //const bool passesCandidateID(const Photon*, unsigned int&);
  bool passesCandidateID(const PhoInfoBranches&, const int, unsigned int&);
  //const bool foundPhotonCandidates(const Handle<PhotonCollection>&, const vector<EcalRecHit*>&, map<unsigned int, Photon*>&, map<unsigned int, Photon*>&);
  bool foundPhotonCandidates(const PhoInfoBranches&, vector<int>&, vector<int>&);
  //const bool foundTrack(const Handle<TrackCollection>&, const map<unsigned int, Photon*>&);
  bool foundTrack(const TrkInfoBranches&, const PhoInfoBranches&, const vector<int>&);
  //const bool electronHasPassingTrack(const Handle<TrackCollection>&, const Photon*, const unsigned int);
  bool electronHasPassingTrack(const TrkInfoBranches&, const PhoInfoBranches&, const int);
  //const bool passDRCut(const map<unsigned int, Photon*>&, const map<unsigned int, Photon*>&);
  bool passDRCut(const vector<int>&, const vector<int>&, const PhoInfoBranches&);
  //const bool photonIsNonoverlapping(const Photon*, const unsigned int, const map<unsigned int, Photon*>&);
  bool photonIsNonoverlapping(const PhoInfoBranches&, const int, const vector<int>&);
  /*const bool foundAllCandidates(const Handle<PhotonCollection>&, const vector<EcalRecHit*>&, const Handle<TrackCollection>&, map<unsigned int, Photon*>&, 
    map<unsigned int, Photon*>&);*/
  bool foundAllCandidates(const PhoInfoBranches&, const TrkInfoBranches&, vector<int>&, vector<int>&);

 private:

  //event selection criteria
  unsigned int sampleType_;
  double ECALIsoMaxPTMultiplierEB_;
  double ECALIsoMaxConstantEB_;
  double ECALIsoMaxPTMultiplierEE_;
  double ECALIsoMaxConstantEE_;
  double HCALIsoMaxPTMultiplierEB_;
  double HCALIsoMaxConstantEB_;
  double HCALIsoMaxPTMultiplierEE_;
  double HCALIsoMaxConstantEE_;
  double HOverEMaxPresel_;
  double ETMin_;
  unsigned int fiducialRegion_;
  bool useHOverE_;
  double HOverEMax_;
  bool useSigmaEtaEta_;
  double sigmaEtaEtaMax_;
  bool useTrackIso_;
  double trackIsoMaxPTMultiplier_;
  double trackIsoMaxConstant_;
  double trackPTMin_;
  double eTrackRMin_;
  double minDRPhotons_;
  bool useTimingCut_;
  double maxSeedTime_;
  unsigned int numReqdCands_;

  //event ID parameters
  unsigned int run_;
  unsigned int evt_;
  unsigned int lumiSec_;

  //debugging output
  string debugFileName_;
  bool debugFlag_;
  ofstream debug_;
};

#endif
