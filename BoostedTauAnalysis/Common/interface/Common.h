#ifndef BoostedTauAnalysis_Common_interface_Common_h
#define BoostedTauAnalysis_Common_interface_Common_h

#include <vector>
#include <algorithm>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

class minDR {

 public:

  void setCandidate(const reco::Candidate*);

  void deleteCandidate();

  reco::Candidate* getCandidate() const;

  template<typename T>
    bool operator()(const T* cand1, const T* cand2)
    {
      return (deltaR(*cand_, *cand1) < deltaR(*cand_, *cand2));
    }

 private:

  reco::Candidate* cand_;

};

class Common {

 public:

  template<typename T, typename U>
    static std::vector<double> sortByProximity(const std::vector<T*>& objectsToSort, 
					       const U& referenceObject)
    {
      std::vector<double> sortedDR2;
      for (typename std::vector<T*>::const_iterator iObject = objectsToSort.begin(); 
	   iObject != objectsToSort.end(); ++iObject) {
	sortedDR2.push_back(reco::deltaR2(referenceObject, **iObject));
      }
      std::sort(sortedDR2.begin(), sortedDR2.end());
      return sortedDR2;
    }

  static void sortByPT(std::vector<reco::Candidate*>&);

  template<typename T, typename U>
    static const T* nearestObject(const U& obj, const std::vector<T*>& objs, unsigned int& index)
    {
      minDR comp;
      comp.setCandidate(dynamic_cast<const reco::Candidate*>(obj.get()));
      typename std::vector<T*>::const_iterator iMinElement = 
	min_element(objs.begin(), objs.end(), comp);
      const T* nearestObj = *iMinElement;
      index = iMinElement - objs.begin();
      comp.deleteCandidate();
      return nearestObj;
    }

  //identify the first good vertex (the "primary" (?))
  static reco::Vertex* getPrimaryVertex(edm::Handle<reco::VertexCollection>&);

  /*fill STL container with muons passing the 2012 tight selection 
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId)*/
  static std::vector<reco::Muon*> getTightRecoMuons(edm::Handle<reco::MuonCollection>&, 
						    const reco::Vertex*);

 private:

  static bool comparePT(reco::Candidate*, reco::Candidate*);

};

#endif
