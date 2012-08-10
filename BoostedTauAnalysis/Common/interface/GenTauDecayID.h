#ifndef BoostedTauAnalysis_Common_interface_GenTauDecayID_h
#define BoostedTauAnalysis_Common_interface_GenTauDecayID_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class GenTauDecayID {

 public:

  //PDG IDs
  enum PDGIDs {        EPDGID = 11,         MUPDGID = 13,         TAUPDGID = 15, 
	       ENEUTRINOPDGID = 12, MUNEUTRINOPDGID = 14, TAUNEUTRINOPDGID = 16, 
		       ZPDGID = 23,          APDGID = 36};

  //decay types
  enum decayType {HAD = 0, MU, E, UNKNOWN};

  //default constructor
  GenTauDecayID();

  //constructor that takes a parameter set of arguments
  GenTauDecayID(edm::ParameterSet&, const edm::Handle<reco::GenParticleCollection>&, 
		const unsigned int);

  //copy constructor
  GenTauDecayID(const GenTauDecayID&);

  //destructor
  ~GenTauDecayID();

  //assignment operator
  GenTauDecayID& operator=(const GenTauDecayID&);

  //getter for parameter set
  edm::ParameterSet* getParSet() const;

  //getter for gen particle handle
  edm::Handle<reco::GenParticleCollection> getGenParticleHandle() const;

  //getter for tau index
  const unsigned int getTauIndex() const;

  //getter for tau sister index; warning suppressed for copy constructor and assignment operator
  const unsigned int getSisterIndex(const bool warn = true) const;

  //true if parameter set has been unpacked
  bool unpacked() const;

  //true if pGenParticles_ points to a valid handle
  bool validHandle() const;

  //find sister
  void findSister();

  //is tau a status 3 decay product?
  bool tauIsStatus3DecayProduct() const;

  //get tau decay type
  unsigned int tauDecayType() const;

  //get tau sister decay type
  unsigned int sisterDecayType() const;

 private:

  //parameters that might be needed in the decay ID (e.g. eta)
  edm::ParameterSet* pParSet_;

  //gen particle collection handle
  edm::Handle<reco::GenParticleCollection> pGenParticles_;

  //index of the tau
  unsigned int iTau_;

  //index of the tau sister
  unsigned int iSister_;

  //PDG ID of tau mother
  unsigned int momPDGID_;

  //true if parameter set has been unpacked
  bool unpacked_;

  //true if pGenParticles_ points to a valid handle
  bool validHandle_;

  //true if tau sister has been found
  bool foundSister_;

  //unpack parameter set
  void unpackParSet();

  //find sister decay product
  int sister() const;

  //classify decay for given particle index
  unsigned int decayType(const unsigned int) const;

  //warning that sister was never found
  std::string warnSisterNeverFound(const std::string&) const;

  //warning that mother PDG ID could not be found
  std::string warnMomPDGIDNotFound(const std::string&) const;

  //error that sister could not be found
  std::string errorSisterNotFound(const std::string&, const unsigned int) const;

  //error that gen particle handle is invalid
  std::string errorInvalidGenParticleHandle(const std::string&) const;

  //error that parameter set is not unpacked
  std::string errorParameterSetNotUnpacked(const std::string&) const;

  //error that supplied parameter set pointer is null
  std::string errorNullParameterSetPointer(const std::string&) const;

  //error that number of daughters is unexpected
  std::string errorUnexpectedNumDaughters(const std::string&, const unsigned int) const;

};

#endif
