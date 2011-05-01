#include "GMSBTools/Filters/interface/Categorizer.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include <iostream>

Categorizer::Categorizer() :
  photonET_(VDOUBLE()),
  photonEta_(VDOUBLE()),
  photonECALIso_(VDOUBLE()),
  photonHCALIso_(VDOUBLE()),
  photonHOverE_(VDOUBLE()),
  photonTrackIso_(VDOUBLE()),
  photonSigmaIetaIeta_(VDOUBLE()),
  photonSeedTime_(VDOUBLE()),
  photonE2OverE9_(VDOUBLE()),
  photonPhi_(VDOUBLE()),
  photonHasPixelSeed_(VBOOL()),
  photonETMin_(-1.0),
  photonAbsEtaMax_(-1.0),
  photonECALIsoMaxPTMultiplier_(-1.0),
  photonECALIsoMaxConstant_(-1.0),
  photonHCALIsoMaxPTMultiplier_(-1.0),
  photonHCALIsoMaxConstant_(-1.0),
  photonHOverEMax_(-1.0),
  photonTrackIsoMaxPTMultiplier_(-1.0),
  photonTrackIsoMaxConstant_(-1.0),
  photonSigmaIetaIetaMax_(-1.0),
  photonAbsSeedTimeMax_(-1.0),
  photonE2OverE9Max_(-1.0),
  photonDPhiMin_(-1.0),
  photonPassETMin_(VBOOL()),
  photonPassAbsEtaMax_(VBOOL()),
  photonPassECALIsoMax_(VBOOL()),
  photonPassHCALIsoMax_(VBOOL()),
  photonPassHOverEMax_(VBOOL()),
  photonPassTrackIsoMax_(VBOOL()),
  photonPassSigmaIetaIetaMax_(VBOOL()),
  photonPassAbsSeedTimeMax_(VBOOL()),
  photonPassE2OverE9Max_(VBOOL()),
  photonPassPreselection_(VBOOL()),
  evtPassDPhiMin_(false),
  evtDiEMET_(-1.0),
  passingPhotons_(VINT(2, -1)),
  photonType_(VINT()),
  category_(FAIL),
  initialized_(false),
  decided_(false),
  categorized_(false),
  foundPassingPhotons_(false) {}

Categorizer::Categorizer(const VDOUBLE& photonET, const VDOUBLE& photonEta, 
			 const VDOUBLE& photonECALIso, const VDOUBLE& photonHCALIso,
			 const VDOUBLE& photonHOverE, const VDOUBLE& photonTrackIso,
			 const VDOUBLE& photonSigmaIetaIeta, const VDOUBLE& photonSeedTime,
			 const VDOUBLE& photonE2OverE9, const VDOUBLE& photonPhi, 
			 const VBOOL& photonHasPixelSeed, const double photonETMin, 
			 const double photonAbsEtaMax, const double photonECALIsoMaxPTMultiplier, 
			 const double photonECALIsoMaxConstant, 
			 const double photonHCALIsoMaxPTMultiplier,
			 const double photonHCALIsoMaxConstant, const double photonHOverEMax,
			 const double photonTrackIsoMaxPTMultiplier,
			 const double photonTrackIsoMaxConstant,
			 const double photonSigmaIetaIetaMax, const double photonAbsSeedTimeMax,
			 const double photonE2OverE9Max, const double photonDPhiMin) :
  photonET_(photonET),
  photonEta_(photonEta),
  photonECALIso_(photonECALIso),
  photonHCALIso_(photonHCALIso),
  photonHOverE_(photonHOverE),
  photonTrackIso_(photonTrackIso),
  photonSigmaIetaIeta_(photonSigmaIetaIeta),
  photonSeedTime_(photonSeedTime),
  photonE2OverE9_(photonE2OverE9),
  photonPhi_(photonPhi),
  photonHasPixelSeed_(photonHasPixelSeed),
  photonETMin_(photonETMin),
  photonAbsEtaMax_(photonAbsEtaMax),
  photonECALIsoMaxPTMultiplier_(photonECALIsoMaxPTMultiplier),
  photonECALIsoMaxConstant_(photonECALIsoMaxConstant),
  photonHCALIsoMaxPTMultiplier_(photonHCALIsoMaxPTMultiplier),
  photonHCALIsoMaxConstant_(photonHCALIsoMaxConstant),
  photonHOverEMax_(photonHOverEMax),
  photonTrackIsoMaxPTMultiplier_(photonTrackIsoMaxPTMultiplier),
  photonTrackIsoMaxConstant_(photonTrackIsoMaxConstant),
  photonSigmaIetaIetaMax_(photonSigmaIetaIetaMax),
  photonAbsSeedTimeMax_(photonAbsSeedTimeMax),
  photonE2OverE9Max_(photonE2OverE9Max),
  photonDPhiMin_(photonDPhiMin),
  photonPassETMin_(VBOOL()),
  photonPassAbsEtaMax_(VBOOL()),
  photonPassECALIsoMax_(VBOOL()),
  photonPassHCALIsoMax_(VBOOL()),
  photonPassHOverEMax_(VBOOL()),
  photonPassTrackIsoMax_(VBOOL()),
  photonPassSigmaIetaIetaMax_(VBOOL()),
  photonPassAbsSeedTimeMax_(VBOOL()),
  photonPassE2OverE9Max_(VBOOL()),
  photonPassPreselection_(VBOOL()),
  evtPassDPhiMin_(false),
  evtDiEMET_(-1.0),
  passingPhotons_(VINT(2, -1)),
  photonType_(VINT()),
  category_(FAIL),
  initialized_(false),
  decided_(false),
  categorized_(false),
  foundPassingPhotons_(false) { checkInput(); }

Categorizer::Categorizer(const Categorizer& copy) :
  photonET_(copy.getPhotonET()),
  photonEta_(copy.getPhotonEta()),
  photonECALIso_(copy.getPhotonECALIso()),
  photonHCALIso_(copy.getPhotonHCALIso()),
  photonHOverE_(copy.getPhotonHOverE()),
  photonTrackIso_(copy.getPhotonTrackIso()),
  photonSigmaIetaIeta_(copy.getPhotonSigmaIetaIeta()),
  photonSeedTime_(copy.getPhotonSeedTime()),
  photonE2OverE9_(copy.getPhotonE2OverE9()),
  photonPhi_(copy.getPhotonPhi()),
  photonHasPixelSeed_(copy.getPhotonHasPixelSeed()),
  photonETMin_(copy.getPhotonETMin()),
  photonAbsEtaMax_(copy.getPhotonAbsEtaMax()),
  photonECALIsoMaxPTMultiplier_(copy.getPhotonECALIsoMaxPTMultiplier()),
  photonECALIsoMaxConstant_(copy.getPhotonECALIsoMaxConstant()),
  photonHCALIsoMaxPTMultiplier_(copy.getPhotonHCALIsoMaxPTMultiplier()),
  photonHCALIsoMaxConstant_(copy.getPhotonHCALIsoMaxConstant()),
  photonHOverEMax_(copy.getPhotonHOverEMax()),
  photonTrackIsoMaxPTMultiplier_(copy.getPhotonTrackIsoMaxPTMultiplier()),
  photonTrackIsoMaxConstant_(copy.getPhotonTrackIsoMaxConstant()),
  photonSigmaIetaIetaMax_(copy.getPhotonSigmaIetaIetaMax()),
  photonAbsSeedTimeMax_(copy.getPhotonAbsSeedTimeMax()),
  photonE2OverE9Max_(copy.getPhotonE2OverE9Max()),
  photonDPhiMin_(copy.getPhotonDPhiMin()),
  photonPassETMin_(copy.getPhotonPassETMin()),
  photonPassAbsEtaMax_(copy.getPhotonPassAbsEtaMax()),
  photonPassECALIsoMax_(copy.getPhotonPassECALIsoMax()),
  photonPassHCALIsoMax_(copy.getPhotonPassHCALIsoMax()),
  photonPassHOverEMax_(copy.getPhotonPassHOverEMax()),
  photonPassTrackIsoMax_(copy.getPhotonPassTrackIsoMax()),
  photonPassSigmaIetaIetaMax_(copy.getPhotonPassSigmaIetaIetaMax()),
  photonPassAbsSeedTimeMax_(copy.getPhotonPassAbsSeedTimeMax()),
  photonPassE2OverE9Max_(copy.getPhotonPassE2OverE9Max()),
  photonPassPreselection_(copy.getPhotonPassPreselection()),
  evtPassDPhiMin_(copy.getEvtPassDPhiMin()),
  evtDiEMET_(copy.getEvtDiEMET()),
  passingPhotons_(copy.getPassingPhotons()),
  photonType_(copy.getPhotonType()),
  category_(copy.getCategory()),
  initialized_(copy.initialized()),
  decided_(copy.decided()),
  categorized_(copy.categorized()),
  foundPassingPhotons_(copy.foundPassingPhotons()) { checkInput(); }

Categorizer::~Categorizer() {}

Categorizer& Categorizer::operator=(const Categorizer& copy)
{
  if (this != &copy) {
    photonET_ = copy.getPhotonET();
    photonEta_ = copy.getPhotonEta();
    photonECALIso_ = copy.getPhotonECALIso();
    photonHCALIso_ = copy.getPhotonHCALIso();
    photonHOverE_ = copy.getPhotonHOverE();
    photonTrackIso_ = copy.getPhotonTrackIso();
    photonSigmaIetaIeta_ = copy.getPhotonSigmaIetaIeta();
    photonSeedTime_ = copy.getPhotonSeedTime();
    photonE2OverE9_ = copy.getPhotonE2OverE9();
    photonPhi_ = copy.getPhotonPhi();
    photonHasPixelSeed_ = copy.getPhotonHasPixelSeed();
    photonETMin_ = copy.getPhotonETMin();
    photonAbsEtaMax_ = copy.getPhotonAbsEtaMax();
    photonECALIsoMaxPTMultiplier_ = copy.getPhotonECALIsoMaxPTMultiplier();
    photonECALIsoMaxConstant_ = copy.getPhotonECALIsoMaxConstant();
    photonHCALIsoMaxPTMultiplier_ = copy.getPhotonHCALIsoMaxPTMultiplier();
    photonHCALIsoMaxConstant_ = copy.getPhotonHCALIsoMaxConstant();
    photonHOverEMax_ = copy.getPhotonHOverEMax();
    photonTrackIsoMaxPTMultiplier_ = copy.getPhotonTrackIsoMaxPTMultiplier();
    photonTrackIsoMaxConstant_ = copy.getPhotonTrackIsoMaxConstant();
    photonSigmaIetaIetaMax_ = copy.getPhotonSigmaIetaIetaMax();
    photonAbsSeedTimeMax_ = copy.getPhotonAbsSeedTimeMax();
    photonE2OverE9Max_ = copy.getPhotonE2OverE9Max();
    photonDPhiMin_ = copy.getPhotonDPhiMin();
    photonPassETMin_ = copy.getPhotonPassETMin();
    photonPassAbsEtaMax_ = copy.getPhotonPassAbsEtaMax();
    photonPassECALIsoMax_ = copy.getPhotonPassECALIsoMax();
    photonPassHCALIsoMax_ = copy.getPhotonPassHCALIsoMax();
    photonPassHOverEMax_ = copy.getPhotonPassHOverEMax();
    photonPassTrackIsoMax_ = copy.getPhotonPassTrackIsoMax();
    photonPassSigmaIetaIetaMax_ = copy.getPhotonPassSigmaIetaIetaMax();
    photonPassAbsSeedTimeMax_ = copy.getPhotonPassAbsSeedTimeMax();
    photonPassE2OverE9Max_ = copy.getPhotonPassE2OverE9Max();
    photonPassPreselection_ = copy.getPhotonPassPreselection();
    evtPassDPhiMin_ = copy.getEvtPassDPhiMin();
    evtDiEMET_ = copy.getEvtDiEMET();
    passingPhotons_ = copy.getPassingPhotons();
    photonType_ = copy.getPhotonType();
    category_ = copy.getCategory();
    initialized_ = copy.initialized();
    decided_ = copy.decided();
    categorized_ = copy.categorized();
    foundPassingPhotons_ = foundPassingPhotons();
    checkInput();
  }
  return *this;
}

VDOUBLE Categorizer::getPhotonET() const { return photonET_; }
VDOUBLE Categorizer::getPhotonEta() const { return photonEta_; }
VDOUBLE Categorizer::getPhotonECALIso() const { return photonECALIso_; }
VDOUBLE Categorizer::getPhotonHCALIso() const { return photonHCALIso_; }
VDOUBLE Categorizer::getPhotonHOverE() const { return photonHOverE_; }
VDOUBLE Categorizer::getPhotonTrackIso() const { return photonTrackIso_; }
VDOUBLE Categorizer::getPhotonSigmaIetaIeta() const { return photonSigmaIetaIeta_; }
VDOUBLE Categorizer::getPhotonSeedTime() const { return photonSeedTime_; }
VDOUBLE Categorizer::getPhotonE2OverE9() const { return photonE2OverE9_; }
VDOUBLE Categorizer::getPhotonPhi() const { return photonPhi_; }
VBOOL Categorizer::getPhotonHasPixelSeed() const { return photonHasPixelSeed_; }
double Categorizer::getPhotonETMin() const { return photonETMin_; }
double Categorizer::getPhotonAbsEtaMax() const { return photonAbsEtaMax_; }
double Categorizer::getPhotonECALIsoMaxPTMultiplier() const
{
  return photonECALIsoMaxPTMultiplier_;
}
double Categorizer::getPhotonECALIsoMaxConstant() const { return photonECALIsoMaxConstant_; }
double Categorizer::getPhotonHCALIsoMaxPTMultiplier() const
{
  return photonHCALIsoMaxPTMultiplier_;
}
double Categorizer::getPhotonHCALIsoMaxConstant() const { return photonHCALIsoMaxConstant_; }
double Categorizer::getPhotonHOverEMax() const { return photonHOverEMax_; }
double Categorizer::getPhotonTrackIsoMaxPTMultiplier() const
{
  return photonTrackIsoMaxPTMultiplier_;
}
double Categorizer::getPhotonTrackIsoMaxConstant() const { return photonTrackIsoMaxConstant_; }
double Categorizer::getPhotonSigmaIetaIetaMax() const { return photonSigmaIetaIetaMax_; }
double Categorizer::getPhotonAbsSeedTimeMax() const { return photonAbsSeedTimeMax_; }
double Categorizer::getPhotonE2OverE9Max() const { return photonE2OverE9Max_; }
double Categorizer::getPhotonDPhiMin() const { return photonDPhiMin_; }
VBOOL Categorizer::getPhotonPassETMin() const
{
  if (!decided_) throw STRING("Error in getPhotonPassETMin(): call decideAll() first.\n");
  return photonPassETMin_;
}
VBOOL Categorizer::getPhotonPassAbsEtaMax() const
{
  if (!decided_) throw STRING("Error in getPhotonPassAbsEtaMax(): call decideAll() first.\n");
  return photonPassAbsEtaMax_;
}
VBOOL Categorizer::getPhotonPassECALIsoMax() const
{
  if (!decided_) throw STRING("Error in getPhotonPassECALIsoMax(): call decideAll() first.\n");
  return photonPassECALIsoMax_;
}
VBOOL Categorizer::getPhotonPassHCALIsoMax() const
{
  if (!decided_) throw STRING("Error in getPhotonPassHCALIsoMax(): call decideAll() first.\n");
  return photonPassHCALIsoMax_;
}
VBOOL Categorizer::getPhotonPassHOverEMax() const
{
  if (!decided_) throw STRING("Error in getPhotonPassHOverEMax(): call decideAll() first.\n");
  return photonPassHOverEMax_;
}
VBOOL Categorizer::getPhotonPassTrackIsoMax() const
{
  if (!decided_) throw STRING("Error in getPhotonPassTrackIsoMax(): call decideAll() first.\n");
  return photonPassTrackIsoMax_;
}
VBOOL Categorizer::getPhotonPassSigmaIetaIetaMax() const
{
  if (!decided_) {
    throw STRING("Error in getPhotonPassSigmaIetaIetaMax(): call decideAll() first.\n");
  }
  return photonPassSigmaIetaIetaMax_;
}
VBOOL Categorizer::getPhotonPassAbsSeedTimeMax() const
{
  if (!decided_) {
    throw STRING("Error in getPhotonPassAbsSeedTimeMax(): call decideAll() first.\n");
  }
  return photonPassAbsSeedTimeMax_;
}
VBOOL Categorizer::getPhotonPassE2OverE9Max() const
{
  if (!decided_) throw STRING("Error in getPhotonPassE2OverE9Max(): call decideAll() first.\n");
  return photonPassE2OverE9Max_;
}
VBOOL Categorizer::getPhotonPassPreselection() const
{
  if (!foundPassingPhotons_) {
    throw STRING("Error in getPhotonPassPreselection(): call findPassingPhotons() first.\n");
  }
  return photonPassPreselection_;
}
bool Categorizer::getEvtPassDPhiMin() const
{
  if (!categorized_) throw STRING("Error in getEvtPassDPhiMin(): call classify() first.\n");
  return evtPassDPhiMin_;
}
double Categorizer::getEvtDiEMET() const
{
  if (!categorized_) throw STRING("Error in getEvtDiEMET(): call classify() first.\n");
  return evtDiEMET_;
}
VINT Categorizer::getPassingPhotons() const
{
  if (!categorized_) throw STRING("Error in getPassingPhotons(): call classify() first.\n");
  return passingPhotons_;
}
VINT Categorizer::getPhotonType() const
{
  if (!categorized_) throw STRING("Error in getPhotonType(): call classify() first.\n");
  return photonType_;
}
int Categorizer::getCategory() const
{
  if (!categorized_) throw STRING("Error in getCategory(): call classify() first.\n");
  return category_;
}

void Categorizer::setPhotonET(const VDOUBLE& photonET)
{
  photonET_ = photonET;
  checkInput();
}
void Categorizer::setPhotonEta(const VDOUBLE& photonEta)
{
  photonEta_ = photonEta;
  checkInput();
}
void Categorizer::setPhotonECALIso(const VDOUBLE& photonECALIso)
{
  photonECALIso_ = photonECALIso;
  checkInput();
}
void Categorizer::setPhotonHCALIso(const VDOUBLE& photonHCALIso)
{
  photonHCALIso_ = photonHCALIso;
  checkInput();
}
void Categorizer::setPhotonHOverE(const VDOUBLE& photonHOverE)
{
  photonHOverE_ = photonHOverE;
  checkInput();
}
void Categorizer::setPhotonTrackIso(const VDOUBLE& photonTrackIso)
{
  photonTrackIso_ = photonTrackIso;
  checkInput();
}
void Categorizer::setPhotonSigmaIetaIeta(const VDOUBLE& photonSigmaIetaIeta)
{
  photonSigmaIetaIeta_ = photonSigmaIetaIeta;
  checkInput();
}
void Categorizer::setPhotonSeedTime(const VDOUBLE& photonSeedTime)
{
  photonSeedTime_ = photonSeedTime;
  checkInput();
}
void Categorizer::setPhotonE2OverE9(const VDOUBLE& photonE2OverE9)
{
  photonE2OverE9_ = photonE2OverE9;
  checkInput();
}
void Categorizer::setPhotonPhi(const VDOUBLE& photonPhi)
{
  photonPhi_ = photonPhi;
  checkInput();
}
void Categorizer::setPhotonHasPixelSeed(const VBOOL& photonHasPixelSeed)
{
  photonHasPixelSeed_ = photonHasPixelSeed;
  checkInput();
}
void Categorizer::setPhotonETMin(const double photonETMin) { photonETMin_ = photonETMin; }
void Categorizer::setPhotonAbsEtaMax(const double photonAbsEtaMax)
{
  photonAbsEtaMax_ = photonAbsEtaMax;
}
void Categorizer::setPhotonECALIsoMaxPTMultiplier(const double photonECALIsoMaxPTMultiplier)
{
  photonECALIsoMaxPTMultiplier_ = photonECALIsoMaxPTMultiplier;
}
void Categorizer::setPhotonECALIsoMaxConstant(const double photonECALIsoMaxConstant)
{
  photonECALIsoMaxConstant_ = photonECALIsoMaxConstant;
}
void Categorizer::setPhotonHCALIsoMaxPTMultiplier(const double photonHCALIsoMaxPTMultiplier)
{
  photonHCALIsoMaxPTMultiplier_ = photonHCALIsoMaxPTMultiplier;
}
void Categorizer::setPhotonHCALIsoMaxConstant(const double photonHCALIsoMaxConstant)
{
  photonHCALIsoMaxConstant_ = photonHCALIsoMaxConstant;
}
void Categorizer::setPhotonHOverEMax(const double photonHOverEMax)
{
  photonHOverEMax_ = photonHOverEMax;
}
void Categorizer::setPhotonTrackIsoMaxPTMultiplier(const double photonTrackIsoMaxPTMultiplier)
{
  photonTrackIsoMaxPTMultiplier_ = photonTrackIsoMaxPTMultiplier;
}
void Categorizer::setPhotonTrackIsoMaxConstant(const double photonTrackIsoMaxConstant)
{
  photonTrackIsoMaxConstant_ = photonTrackIsoMaxConstant;
}
void Categorizer::setPhotonSigmaIetaIetaMax(const double photonSigmaIetaIetaMax)
{
  photonSigmaIetaIetaMax_ = photonSigmaIetaIetaMax;
}
void Categorizer::setPhotonAbsSeedTimeMax(const double photonAbsSeedTimeMax)
{
  photonAbsSeedTimeMax_ = photonAbsSeedTimeMax;
}
void Categorizer::setPhotonE2OverE9Max(const double photonE2OverE9Max)
{
  photonE2OverE9Max_ = photonE2OverE9Max;
}
void Categorizer::setPhotonDPhiMin(const double photonDPhiMin) { photonDPhiMin_ = photonDPhiMin; }

void Categorizer::decideAll()
{
  try {
    decideETMin();
    decideAbsEtaMax();
    decideECALIsoMax();
    decideHCALIsoMax();
    decideHOverEMax();
    decideTrackIsoMax();
    decideSigmaIetaIetaMax();
    decideAbsSeedTimeMax();
    decideE2OverE9Max();
  }
  catch (STRING& vecSizeMismatch) { throw; }
  decided_ = true;
}

void Categorizer::findPassingPhotons()
{
  if (!decided_) throw STRING("Error in findPassingPhotons(): call decideAll() first.\n");
  photonPassPreselection_.clear();
  for (VBOOL_IT iPassETMin = photonPassETMin_.begin(); iPassETMin != photonPassETMin_.end(); 
       ++iPassETMin) {
    const unsigned int i = iPassETMin - photonPassETMin_.begin();
    photonPassPreselection_.push_back(*iPassETMin && photonPassAbsEtaMax_[i] && 
				      photonPassECALIsoMax_[i] && photonPassHCALIsoMax_[i] && 
				      photonPassHOverEMax_[i] && photonPassAbsSeedTimeMax_[i] && 
				      photonPassE2OverE9Max_[i]);
  }
  foundPassingPhotons_ = true;
}

void Categorizer::classify()
{
  if (!foundPassingPhotons_) {
    throw STRING("Error in classify(): call findPassingPhotons() first.\n");
  }
  double photon1ET = -1.0;
  double photon2ET = -1.0;
  photonType_.clear();
  for (VBOOL_IT iPassETMin = photonPassETMin_.begin(); iPassETMin != photonPassETMin_.end(); 
       ++iPassETMin) {
    const unsigned int i = iPassETMin - photonPassETMin_.begin();
    if (photonPassPreselection_[i]) {
      if (photonET_[i] > photon1ET) {
	photon2ET = photon1ET;
	photon1ET = photonET_[i];
	passingPhotons_[1] = passingPhotons_[0];
	passingPhotons_[0] = i;
      }
      else if (photonET_[i] > photon2ET) {
	passingPhotons_[1] = i;
	photon2ET = photonET_[i];
      }
      if (photonPassSigmaIetaIetaMax_[i] && 
	  photonPassTrackIsoMax_[i]) {
	if (photonHasPixelSeed_[i]) photonType_.push_back(E);
	else photonType_.push_back(G);
      }
      else photonType_.push_back(F);
    }
    else photonType_.push_back(FAIL);
  }
  if ((passingPhotons_[0] >= 0) && (passingPhotons_[1] >= 0)) {
    if ((photonType_[passingPhotons_[0]] == G) && 
	(photonType_[passingPhotons_[1]] == G)) category_ = GG;
    if (((photonType_[passingPhotons_[0]] == G) && (photonType_[passingPhotons_[1]] == E)) || 
	((photonType_[passingPhotons_[0]] == E) && 
	 (photonType_[passingPhotons_[1]] == G))) category_ = EG;
    if ((photonType_[passingPhotons_[0]] == E) && 
	(photonType_[passingPhotons_[1]] == E)) category_ = EE;
    if ((photonType_[passingPhotons_[0]] == F) && 
	(photonType_[passingPhotons_[1]] == F)) category_ = FF;
    if (deltaPhi(photonPhi_[passingPhotons_[0]], photonPhi_[passingPhotons_[1]]) < 
	photonDPhiMin_) {
      category_ = FAIL;
    }
    else evtPassDPhiMin_ = true;
    evtDiEMET_ = sqrt(photonET_[passingPhotons_[0]]*photonET_[passingPhotons_[0]] + 
		      photonET_[passingPhotons_[1]]*photonET_[passingPhotons_[1]] + 
		      2*photonET_[passingPhotons_[0]]*photonET_[passingPhotons_[1]]*
		      cos(deltaPhi(photonPhi_[passingPhotons_[0]], 
				   photonPhi_[passingPhotons_[1]])));
  }
  else category_ = FAIL;
  categorized_ = true;
}

bool Categorizer::done() const
{
  return initialized_ && decided_ && foundPassingPhotons_ && categorized_;
}

bool Categorizer::initialized() const { return initialized_; }

bool Categorizer::decided() const { return decided_; }

bool Categorizer::categorized() const { return categorized_; }

bool Categorizer::foundPassingPhotons() const { return foundPassingPhotons_; }

void Categorizer::checkInput()
{
  initialized_ = 
    ((photonET_.size() == photonEta_.size()) && (photonEta_.size() == photonECALIso_.size()) && 
     (photonECALIso_.size() == photonHCALIso_.size()) && 
     (photonHCALIso_.size() == photonHOverE_.size()) && 
     (photonHOverE_.size() == photonTrackIso_.size()) && 
     (photonTrackIso_.size() == photonSigmaIetaIeta_.size()) && 
     (photonSigmaIetaIeta_.size() == photonSeedTime_.size()) && 
     (photonSeedTime_.size() == photonE2OverE9_.size()) && 
     (photonE2OverE9_.size() == photonPhi_.size()) && 
     (photonPhi_.size() == photonHasPixelSeed_.size()));
}

void Categorizer::decide(VBOOL& pass, const VDOUBLE& quantity, const double cut0, 
			 const double cut1, const bool min, const bool abs)
{
  if (!initialized_) {
    STRINGSTREAM err;
    err << "Error in decide(VBOOL& pass, const VDOUBLE& quantity, const double cut0, ";
    err << "const double cut1, const bool min, const bool abs): input vector size mismatch.\n";
    throw STRING(err.str());
  }
  pass.clear();
  for (VDOUBLE_IT i = quantity.begin(); i != quantity.end(); ++i) {
    double val = *i;
    if (abs) val = fabs(*i);
    bool decision = (cut0 == -1.0) || (cut1 == -1.0);
    if (min) decision = decision || (val > (cut1*photonET_[i - quantity.begin()] + cut0));
    else decision = decision || (val < (cut1*photonET_[i - quantity.begin()] + cut0));
    pass.push_back(decision);
  }
}

void Categorizer::decide1D(VBOOL& pass, const VDOUBLE& quantity, const double cut, 
			   const bool min, const bool abs)
{
  try { decide(pass, quantity, cut, 0.0, min, abs); }
  catch (STRING& vecSizeMismatch) { throw; }
}

void Categorizer::decide1DMin(VBOOL& pass, const VDOUBLE& quantity, const double cut)
{
  try { decide1D(pass, quantity, cut, true, false); }
  catch (STRING& vecSizeMismatch) { throw; }
}

void Categorizer::decide1DMax(VBOOL& pass, const VDOUBLE& quantity, const double cut)
{
  try { decide1D(pass, quantity, cut, false, false); }
  catch (STRING& vecSizeMismatch) { throw; }
}

void Categorizer::decide1DMinAbs(VBOOL& pass, const VDOUBLE& quantity, const double cut)
{
  try { decide1D(pass, quantity, cut, true, true); }
  catch (STRING& vecSizeMismatch) { throw; }
}

void Categorizer::decide1DMaxAbs(VBOOL& pass, const VDOUBLE& quantity, const double cut)
{
  try { decide1D(pass, quantity, cut, false, true); }
  catch (STRING& vecSizeMismatch) { throw; }
}

void Categorizer::decide2DMax(VBOOL& pass, const VDOUBLE& quantity, const double cut0, 
			      const double cut1)
{
  try { decide(pass, quantity, cut0, cut1, false, false); }
  catch (STRING& vecSizeMismatch) { throw; }
}

void Categorizer::decideETMin()
{
  try { decide1DMin(photonPassETMin_, photonET_, photonETMin_); }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideAbsEtaMax()
{
  try { decide1DMaxAbs(photonPassAbsEtaMax_, photonEta_, photonAbsEtaMax_); }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideECALIsoMax()
{
  try {
    decide2DMax(photonPassECALIsoMax_, photonECALIso_, photonECALIsoMaxConstant_, 
		photonECALIsoMaxPTMultiplier_);
  }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideHCALIsoMax()
{
  try {
    decide2DMax(photonPassHCALIsoMax_, photonHCALIso_, photonHCALIsoMaxConstant_, 
		photonHCALIsoMaxPTMultiplier_);
  }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideHOverEMax()
{
  try { decide1DMax(photonPassHOverEMax_, photonHOverE_, photonHOverEMax_); }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideTrackIsoMax()
{
  try {
    decide2DMax(photonPassTrackIsoMax_, photonTrackIso_, photonTrackIsoMaxConstant_, 
		photonTrackIsoMaxPTMultiplier_);
  }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideSigmaIetaIetaMax()
{
  try { decide1DMax(photonPassSigmaIetaIetaMax_, photonSigmaIetaIeta_, photonSigmaIetaIetaMax_); }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideAbsSeedTimeMax()
{
  try { decide1DMaxAbs(photonPassAbsSeedTimeMax_, photonSeedTime_, photonAbsSeedTimeMax_); }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideE2OverE9Max()
{
  try { decide1DMax(photonPassE2OverE9Max_, photonE2OverE9_, photonE2OverE9Max_); }
  catch (STRING& vecSizeMismatch) { throw; }
}
