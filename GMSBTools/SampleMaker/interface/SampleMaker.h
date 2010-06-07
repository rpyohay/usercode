#ifndef SAMPLEMAKER_H
#define SAMPLEMAKER_H

// system include files
#include <memory>
#include <string>
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

#endif
