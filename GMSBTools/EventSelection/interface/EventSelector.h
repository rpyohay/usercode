#ifndef SampleMaker_EventSelector_h
#define SampleMaker_EventSelector_h

// system include files
#include <memory>
#include <sstream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

//reco::Photon classes
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

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
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

//ROOT
#include "TMath.h"

//sample types
#define GAMMAGAMMA 1
#define EGAMMA 2
#define EE 3
#define FF 4
#define ETRACK 5

//ECAL fiducial regions
#define EB 1
#define EEND 2
#define ECAL 3

using namespace std;
using namespace reco;
using namespace edm;
//using namespace cms;

class EventSelector {

 public:

  //default constructor
  EventSelector();

  //constructor from a list of cut values
  EventSelector(const unsigned int, const double, const double, const double, const double, const double, const double, const double, const double, 
		const double, const double, const unsigned int, const bool, const double, const bool, const double, const bool, const double, const double, 
		const double, const double, const double, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const string, 
		const bool);

  //copy constructor
  EventSelector(/*const */EventSelector&);

  //destructor
  ~EventSelector();

  //assignment operator
  EventSelector& operator=(/*const */EventSelector&);

  //getters
  const unsigned int getSampleType() const;
  const double getECALIsoMaxPTMultiplierEB() const;
  const double getECALIsoMaxConstantEB() const;
  const double getECALIsoMaxPTMultiplierEE() const;
  const double getECALIsoMaxConstantEE() const;
  const double getHCALIsoMaxPTMultiplierEB() const;
  const double getHCALIsoMaxConstantEB() const;
  const double getHCALIsoMaxPTMultiplierEE() const;
  const double getHCALIsoMaxConstantEE() const;
  const double getHOverEMaxPresel() const;
  const double getETMin() const;
  const unsigned int getFiducialRegion() const;
  const bool getUseHOverE() const;
  const double getHOverEMax() const;
  const bool getUseSigmaEtaEta() const;
  const double getSigmaEtaEtaMax() const;
  const bool getUseTrackIso() const;
  const double getTrackIsoMaxPTMultiplier() const;
  const double getTrackIsoMaxConstant() const;
  const double getTrackPTMin() const;
  const double getETrackRMin() const;
  const double getMinDRPhotons() const;
  const unsigned int getNumReqdCands() const;
  const unsigned int getRun() const;
  const unsigned int getEvt() const;
  const unsigned int getLumiSec() const;
  const string getDebugFileName() const;
  const bool getDebugFlag() const;

  //setters
  void setRun(const unsigned int);
  void setEvt(const unsigned int);
  void setLumiSec(const unsigned int);

  //helper methods
  void printEvtInfo();
  void printDebug(const string);
  void initializeDebugging();
  const bool debugFileOpen() const;
  void closeDebugFile();
  const string sampleTypeString() const;
  const unsigned int ECALFiducialRegion(const Photon*) const;

  //calculate dPhi between two objects
  template <typename T>
  const T dPhi(T phi1, T phi2) const
  {
    T dPhi = fabs(phi1 - phi2);
    if (dPhi > TMath::Pi()) dPhi = TMath::TwoPi() - dPhi;
    return dPhi;
  }

  //calculate dEta between two objects
  template <typename T>
  const T dEta(T eta1, T eta2) const { return fabs(eta1 - eta2); }

  //calculate dR between two objects
  template <typename T>
  const T dR(T eta1, T eta2, T phi1, T phi2) const { return sqrt(dEta(eta1, eta2)*dEta(eta1, eta2) + dPhi(phi1, phi2)*dPhi(phi1, phi2)); }

  //pass flags
  const bool passesDataQualityCuts(const Handle<HBHERecHitCollection>&, const vector<Photon*>&, const CaloGeometry*, const Handle<TrackCollection>&);
  const bool passesHEBeamHaloTag(const Handle<HBHERecHitCollection>&, const vector<Photon*>&, const CaloGeometry*);
  const bool passesMuonBeamHaloTag(const Handle<TrackCollection>&, const vector<Photon*>&);
  const bool passesPreselection(const Photon*, const unsigned int);
  const bool passesCandidateID(const Photon*, bool&, bool&);
  const bool foundPhotonCandidates(const Handle<PhotonCollection>&, vector<Photon*>&);
  const bool foundTrack(const Handle<TrackCollection>&, const vector<Photon*>&);
  const bool passDRCut(const vector<Photon*>&);
  const bool foundAllCandidates(const Handle<PhotonCollection>&, const Handle<TrackCollection>&, vector<Photon*>&);

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
