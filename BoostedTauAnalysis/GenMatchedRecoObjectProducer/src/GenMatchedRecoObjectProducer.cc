// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      GenMatchedRecoObjectProducer
// 
/**\class GenMatchedRecoObjectProducer GenMatchedRecoObjectProducer.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/GenMatchedRecoObjectProducer.cc

Description: produce a collection of reco objects matched to gen boosted di-tau objects

Implementation:
this module is designed to only match the primary object of the boosted di-tau pair, so it 
has to be run again with the sister designated as the primary to do sister matching
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu Aug 23 11:23:58 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"


//
// class declaration
//

template<class T>
class GenMatchedRecoObjectProducer : public edm::EDFilter {
public:
  explicit GenMatchedRecoObjectProducer(const edm::ParameterSet&);
  ~GenMatchedRecoObjectProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  //input tag for gen particle collection
  edm::InputTag genParticleTag_;

  //input tag reco object collection
  edm::InputTag recoObjTag_;

  //PDG ID of mother of gen matching object
  int momPDGID_;

  //decay type of one member of the tau pair, arbitrarily designated the primary
  GenTauDecayID::DecayType primaryTauDecayType_;

  //decay type of the second member of the tau pair, arbitrarily designated the sister
  GenTauDecayID::DecayType sisterTauDecayType_;

  //minimum pT to be counted as a charged hadron in hadronic tau decay
  double chargedHadronPTMin_;

  //minimum pT to be counted as a neutral hadron in hadronic tau decay
  double neutralHadronPTMin_;

  //minimum pT to be counted as a charged lepton in leptonic tau decay or mother decay
  double chargedLeptonPTMin_;

  //minimum pT of the visible tau decay products or gen lepton from mother decay
  double totalPTMin_;

  //flag indicating whether pT cuts should be applied in determining valid gen objects
  bool applyPTCuts_;

  //flag indicating whether KShorts should be counted as neutral hadrons
  double countKShort_;

  //pointer to the parameter set supplied to this module
  edm::ParameterSet* cfg_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
template<class T>
GenMatchedRecoObjectProducer<T>::GenMatchedRecoObjectProducer(const edm::ParameterSet& iConfig) :
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  recoObjTag_(iConfig.getParameter<edm::InputTag>("recoObjTag")),
  momPDGID_(iConfig.getParameter<int>("momPDGID")),
  primaryTauDecayType_(static_cast<GenTauDecayID::DecayType>
		       (iConfig.getParameter<unsigned int>("primaryTauDecayType"))),
  sisterTauDecayType_(static_cast<GenTauDecayID::DecayType>
		      (iConfig.getParameter<unsigned int>("sisterTauDecayType"))),
  chargedHadronPTMin_(iConfig.getParameter<double>("chargedHadronPTMin")),
  neutralHadronPTMin_(iConfig.getParameter<double>("neutralHadronPTMin")),
  chargedLeptonPTMin_(iConfig.getParameter<double>("chargedLeptonPTMin")),
  totalPTMin_(iConfig.getParameter<double>("totalPTMin")),
  applyPTCuts_(iConfig.getParameter<bool>("applyPTCuts")),
  countKShort_(iConfig.getParameter<bool>("countKShort")),
  cfg_(const_cast<edm::ParameterSet*>(&iConfig))
{
  //register your products
  produces<std::vector<T> >();

  //now do what ever other initialization is needed
  
}


template<class T>
GenMatchedRecoObjectProducer<T>::~GenMatchedRecoObjectProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get GEN particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //get reco object collection
  edm::Handle<std::vector<T> > pRecoObjs;
  iEvent.getByLabel(recoObjTag_, pRecoObjs);

//   //get HPS discriminators
//   for (std::vector<InputTag>::const_iterator 
// 	 iHPSDiscriminatorTag = HPSDiscriminatorTags_.begin(); 
//        iHPSDiscriminatorTag != HPSDiscriminatorTags_.end(); ++iHPSDiscriminatorTag) {
//     edm::Handle<reco::PFTauDiscriminator> pHPSDiscriminator;
//     iEvent.getByLabel(*iHPSDiscriminatorTag, pHPSDiscriminator);
//     HPSDiscriminators_[iHPSDiscriminatorTag->label()] = pHPSDiscriminator;
//   }

//   //fill STL container with HPS taus passing desired discriminators
//   std::vector<reco::PFTau*> HPSTaus;
//   for (reco::PFTauCollection::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
//        ++iTau) {
//     reco::PFTauRef tauRef(pTaus, iTau - pTaus->begin());
//     bool passFlag = true;
//     std::map<std::string, edm::Handle<reco::PFTauDiscriminator> >::const_iterator iHPSDiscriminator = 
//       HPSDiscriminators_.begin();
//     while ((iHPSDiscriminator != HPSDiscriminators_.end()) && passFlag) {
//       if ((*(iHPSDiscriminator->second))[tauRef] != 1.0) passFlag = false;
//       ++iHPSDiscriminator;
//     }
//     if (passFlag) HPSTaus.push_back(const_cast<reco::PFTau>(&*iTau));
//   }

//fill STL container of pointers to reco objects
  std::vector<T*> recoObjPtrs;
  for (typename std::vector<T>::const_iterator iRecoObj = pRecoObjs->begin(); 
       iRecoObj != pRecoObjs->end(); ++iRecoObj) {
    recoObjPtrs.push_back(const_cast<T*>(&*iRecoObj));
  }

  //fill STL container of gen tau decays
  std::vector<GenTauDecayID> tauDecays;
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {
    try {
      GenTauDecayID tauDecay(*cfg_, pGenParticles, iGenParticle - pGenParticles->begin());
      if (tauDecay.tauIsStatus3DecayProduct()) tauDecays.push_back(tauDecay);
    }
    catch (std::string& ex) { throw cms::Exception("GenMatchedRecoObjectProducer") << ex; }
  }

  //fill STL container of gen muons (i.e. from a-->mumu)

  //declare pointers to output collection to produce
  std::auto_ptr<std::vector<T> > genMatchedRecoObjs(new std::vector<T>);

  //loop over gen tau decays
  std::vector<unsigned int> keysToIgnore;
  for (std::vector<GenTauDecayID>::iterator iTau = tauDecays.begin(); iTau != tauDecays.end(); 
       ++iTau) {
    try {

      /*select gen tau decays of the desired type, passing desired pT cut, with the desired sister 
	decay type*/
      if (iTau->tauDecayType(applyPTCuts_, countKShort_).second == primaryTauDecayType_) {
	iTau->findSister();
	const unsigned int iSister = iTau->getSisterIndex();
	if ((std::find(keysToIgnore.begin(), keysToIgnore.end(), iSister) == 
	     keysToIgnore.end()) && /*keysToIgnore keeps track of whether you've looped over the 
				      other half of the boosted di-tau pair before*/
	    (iTau->sisterDecayType(applyPTCuts_, countKShort_).second == sisterTauDecayType_)) {

	  //will the following work if the gen object to match is not a tau?

	  //get the 4-vector of the visible decay products of the tau
	  reco::LeafCandidate::LorentzVector visibleGenP4 = iTau->getVisibleTauP4();

	  //make a dummy LeafCandidate and ref out of the visible 4-vector
	  std::vector<reco::LeafCandidate> 
	    visibleGenTau(1, reco::LeafCandidate(0.0, visibleGenP4));
	  edm::Ref<std::vector<reco::LeafCandidate> > visibleGenTauRef(&visibleGenTau, 0);

	  //find the nearest reco object to the gen tau
	  unsigned int nearestRecoObjKey = 0;
	  const T* nearestRecoObj = 
	    Common::nearestObject(visibleGenTauRef, recoObjPtrs, nearestRecoObjKey);

	  //if nearest reco object is within 0.3 of the gen tau, save this reco object
	  if (reco::deltaR(*nearestRecoObj, *visibleGenTauRef) < 0.3) {
	    genMatchedRecoObjs->push_back(*nearestRecoObj);
	  }
	}

	//add this tau's key to the list of keys to ignore when we get to its sister
	keysToIgnore.push_back(iTau->getTauIndex());
      }
    }
    catch (std::string& ex) { throw cms::Exception("GenMatchedRecoObjectProducer") << ex; }
  }

  //flag indicating whether >0 gen-matched reco objects were found
  const bool foundGenMatchedRecoObject = genMatchedRecoObjs->size() > 0;

  //put output collection into event
  iEvent.put(genMatchedRecoObjs); //this function frees the auto_ptr argument

  //stop processing if no gen-matched objects were found
  return foundGenMatchedRecoObject;
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void GenMatchedRecoObjectProducer<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void GenMatchedRecoObjectProducer<T>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::beginLuminosityBlock(edm::LuminosityBlock&, 
							   edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool GenMatchedRecoObjectProducer<T>::endLuminosityBlock(edm::LuminosityBlock&, 
							 edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void 
GenMatchedRecoObjectProducer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
typedef GenMatchedRecoObjectProducer<reco::Muon> GenMatchedMuonProducer;
DEFINE_FWK_MODULE(GenMatchedMuonProducer);
