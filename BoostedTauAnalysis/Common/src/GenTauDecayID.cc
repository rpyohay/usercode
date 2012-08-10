#include <string>
#include <sstream>
#include "BoostedTauAnalysis/Common/interface/GenTauDecayID.h"

//default constructor
GenTauDecayID::GenTauDecayID() :
  pParSet_(NULL),
  iTau_(0),
  iSister_(0),
  momPDGID_(0),
  unpacked_(false),
  validHandle_(false),
  foundSister_(false) {}

//constructor that takes a parameter set of arguments
GenTauDecayID::GenTauDecayID(edm::ParameterSet& parSet, 
			     const edm::Handle<reco::GenParticleCollection>& pGenParticles, 
			     const unsigned int iTau) :
  pParSet_(&parSet),
  pGenParticles_(pGenParticles),
  iTau_(iTau),
  iSister_(0),
  momPDGID_(0),
  unpacked_(false),
  validHandle_(pGenParticles_.isValid()),
  foundSister_(false)
{
  try { unpackParSet(); }
  catch (std::string& ex) { throw; }
}

//copy constructor
GenTauDecayID::GenTauDecayID(const GenTauDecayID& other) :
  pParSet_(other.getParSet()),
  pGenParticles_(other.getGenParticleHandle()),
  iTau_(other.getTauIndex()),
  iSister_(other.getSisterIndex(false)),
  momPDGID_(0),
  unpacked_(false),
  validHandle_(pGenParticles_.isValid()),
  foundSister_(false)
{
  try { unpackParSet(); }
  catch (std::string& ex) { throw; }
}

//destructor
GenTauDecayID::~GenTauDecayID()
{
  pParSet_ = NULL;
  pGenParticles_.clear();
  iTau_ = 0;
  iSister_ = 0;
  momPDGID_ = 0;
  unpacked_ = false;
  validHandle_ = false;
  foundSister_ = false;
}

//assignment operator
GenTauDecayID& GenTauDecayID::operator=(const GenTauDecayID& rhs)
{
  if (this != &rhs) {
    pParSet_ = rhs.getParSet();
    pGenParticles_ = rhs.getGenParticleHandle();
    iTau_ = rhs.getTauIndex();
    iSister_ = rhs.getSisterIndex(false);
    momPDGID_ = 0;
    unpacked_ = false;
    validHandle_ = pGenParticles_.isValid();
    foundSister_ = false;
    try { unpackParSet(); }
    catch (std::string& ex) { throw; }
  }
  return *this;
}

//getter for parameter set
edm::ParameterSet* GenTauDecayID::getParSet() const { return pParSet_; }

//getter for gen particle handle
edm::Handle<reco::GenParticleCollection> GenTauDecayID::getGenParticleHandle() const
{
  return pGenParticles_;
}

//getter for tau index
const unsigned int GenTauDecayID::getTauIndex() const { return iTau_; }

//getter for tau sister index
const unsigned int GenTauDecayID::getSisterIndex(const bool warn) const
{
  if (warn && !foundSister_) {
    std::cerr << warnSisterNeverFound("const unsigned int GenTauDecayID::getSisterIndex() const");
  }
  return iSister_;
}

//true if parameter set has been unpacked
bool GenTauDecayID::unpacked() const { return unpacked_; }

//true if pGenParticles_ points to a valid handle
bool GenTauDecayID::validHandle() const { return validHandle_; }

//find sister
void GenTauDecayID::findSister()
{
  if (!foundSister_) {
    int iSister = -1;
    try { iSister = sister(); }
    catch (std::string& ex) { throw; }
    if (iSister != -1) {
      iSister_ = iSister;
      foundSister_ = true;
    }
    else throw errorSisterNotFound("void GenTauDecayID::findSister()", iTau_);
  }
}

//is tau a status 3 decay product?
bool GenTauDecayID::tauIsStatus3DecayProduct() const
{
  std::string fnName("bool GenTauDecayID::tauIsStatus3DecayProduct() const");
  bool ret = false;
  if (unpacked_) {
    if (validHandle_) {
      reco::GenParticleRef tauRef(pGenParticles_, iTau_);
      ret = ((fabs(tauRef->pdgId()) == TAUPDGID) && (tauRef->status() == 3) && 
	     (tauRef->numberOfMothers() == 1) && (tauRef->mother(0)->pdgId() == (int)momPDGID_)); 
    }
    else throw errorInvalidGenParticleHandle(fnName);
  }
  else throw errorParameterSetNotUnpacked(fnName);
  return ret;
}

//get tau decay type
unsigned int GenTauDecayID::tauDecayType() const
{
  unsigned int decayTypeCode = 0;
  try { decayTypeCode = decayType(iTau_); }
  catch (std::string& ex) { throw; }
  return decayTypeCode;
}

//get tau sister decay type
unsigned int GenTauDecayID::sisterDecayType() const
{
  if (!foundSister_) {
    std::cerr << warnSisterNeverFound("unsigned int GenTauDecayID::sisterDecayType() const");
  }
  unsigned int decayTypeCode = 0;
  try { decayTypeCode = decayType(iSister_); }
  catch (std::string& ex) { throw; }
  return decayTypeCode;
}

//unpack parameter set
void GenTauDecayID::unpackParSet()
{
  std::string fnName("void GenTauDecayID::unpackParSet()");
  if (pParSet_ == NULL) {
    throw errorNullParameterSetPointer(fnName);
  }
  if (pParSet_->existsAs<unsigned int>("momPDGID")) {
    momPDGID_ = pParSet_->getParameter<unsigned int>("momPDGID");
  }
  else {
    std::cerr << warnMomPDGIDNotFound(fnName);
    momPDGID_ = 0;
  }
  unpacked_ = true;
}

//find first sister decay product
int GenTauDecayID::sister() const
{
  int iSister = -1;
  if (validHandle_) {
    reco::GenParticleRef tauRef(pGenParticles_, iTau_);
    reco::GenParticleRef momRef = tauRef->motherRef();
    unsigned int iDaughter = 0;
    while ((iDaughter < momRef->numberOfDaughters()) && (iSister == -1)) {
      reco::GenParticleRef childRef = momRef->daughterRef(iDaughter);
      unsigned int childRefKey = childRef.key();
      if ((childRefKey != iTau_) && (childRef->status() == tauRef->status())) {
	iSister = childRefKey;
      }
      ++iDaughter;
    }
  }
  else throw errorInvalidGenParticleHandle("int GenTauDecayID::sister() const");
  return iSister;
}

//classify decay
unsigned int GenTauDecayID::decayType(const unsigned int iParticle) const
{
  std::string fnName("unsigned int GenTauDecayID::decayType() const");
  unsigned int ret = UNKNOWN;
  unsigned int leptonPDGID = 0;
  unsigned int neutrinoPDGID = 0;
  if (validHandle_) {
    reco::GenParticleRef tauRef(pGenParticles_, iParticle);
    const size_t numDaughters = tauRef->numberOfDaughters();
    if (numDaughters == 1) {
      const reco::Candidate* daughter = tauRef->daughter(0);
      reco::Candidate::const_iterator iDaughter = daughter->begin();
      while ((iDaughter != daughter->end()) && 
	     ((leptonPDGID == 0) || (neutrinoPDGID == 0))) {
	if (leptonPDGID == 0) {
	  if (((neutrinoPDGID == 0) || (neutrinoPDGID == ENEUTRINOPDGID)) && 
	      (fabs(iDaughter->pdgId()) == EPDGID)) leptonPDGID = EPDGID;
	  if (((neutrinoPDGID == 0) || (neutrinoPDGID == MUNEUTRINOPDGID)) && 
	      (fabs(iDaughter->pdgId()) == MUPDGID)) leptonPDGID = MUPDGID;
	}
	if (neutrinoPDGID == 0) {
	  if (((leptonPDGID == 0) || (leptonPDGID == EPDGID)) && 
	      (fabs(iDaughter->pdgId()) == ENEUTRINOPDGID)) neutrinoPDGID = ENEUTRINOPDGID;
	  if (((leptonPDGID == 0) || (leptonPDGID == MUPDGID)) && 
	      (fabs(iDaughter->pdgId()) == MUNEUTRINOPDGID)) {
	    neutrinoPDGID = MUNEUTRINOPDGID;
	  }
	}
	++iDaughter;
      }
    }
    else throw errorUnexpectedNumDaughters(fnName, numDaughters);
  }
  else throw errorInvalidGenParticleHandle(fnName);
  if ((leptonPDGID == 0) || (neutrinoPDGID == 0)) ret = HAD;
  if ((leptonPDGID == EPDGID) || (neutrinoPDGID == ENEUTRINOPDGID)) ret = E;
  if ((leptonPDGID == MUPDGID) || (neutrinoPDGID == MUNEUTRINOPDGID)) ret = MU;
  return ret;
}

//warning that sister was never found
std::string GenTauDecayID::warnSisterNeverFound(const std::string& fnName) const
{
  return ("Warning in " + fnName + 
	  ":\nSister index should not be trusted because sister was never found.\n");
}

//warning that mother PDG ID could not be found
std::string GenTauDecayID::warnMomPDGIDNotFound(const std::string& fnName) const
{
  return ("Warning in " + fnName + ":\nint momPDGID not found.\n");
}

//error that sister could not be found
std::string GenTauDecayID::errorSisterNotFound(const std::string& fnName, 
					       const unsigned int iTau) const
{
  std::stringstream err;
  err << "Error in " << fnName << ":\nCould not find sister of tau with index " << iTau_ << ".\n";
  return (err.str());
}

//error that gen particle handle is invalid
std::string GenTauDecayID::errorInvalidGenParticleHandle(const std::string& fnName) const
{
  return ("Error in " + fnName + ":\nInvalid handle.\n");
}

//error that parameter set is not unpacked
std::string GenTauDecayID::errorParameterSetNotUnpacked(const std::string& fnName) const
{
  return ("Error in " + fnName + ":\nParameter set not unpacked.\n");
}

//error that supplied parameter set pointer is null
std::string GenTauDecayID::errorNullParameterSetPointer(const std::string& fnName) const
{
  return ("Error in " + fnName + ":\nedm::ParameterSet* parSet_ is NULL.\n");
}

//error that number of daughters is unexpected
std::string GenTauDecayID::errorUnexpectedNumDaughters(const std::string& fnName, 
						       const unsigned int numDaughters) const
{
  std::stringstream err;
  err << "Error in " << fnName << ":\npTau_ has " << numDaughters << " daughters.\n";
  return err.str();
}
