#include "GMSBTools/Filters/interface/Categorizer.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>
#include <utility>

Categorizer::Categorizer() :
  photonET_(VDOUBLE()),
  photonEta_(VDOUBLE()),
  photonECALIso_(VDOUBLE()),
  photonHCALIso_(VDOUBLE()),
  photonHOverE_(VDOUBLE()),
  photonR9_(VDOUBLE()),
  photonTrackIso_(VDOUBLE()),
  photonSigmaIetaIeta_(VDOUBLE()),
  photonSeedTime_(VDOUBLE()),
  photonE2OverE9_(VDOUBLE()),
  photonPhi_(VDOUBLE()),
  photonHasPixelSeed_(VBOOL()),
  photon1ETMin_(-1.0),
  photon2ETMin_(-1.0),
  photonAbsEtaMax_(-1.0),
  photonECALIsoMaxPTMultiplier_(-1.0),
  photonECALIsoMaxConstant_(-1.0),
  photonHCALIsoMaxPTMultiplier_(-1.0),
  photonHCALIsoMaxConstant_(-1.0),
  photonHOverEMax_(-1.0),
  photonR9Max_(-1.0),
  photonTrackIsoMaxPTMultiplier_(-1.0),
  photonTrackIsoMaxConstant_(-1.0),
  photonSigmaIetaIetaMax_(-1.0),
  photonAbsSeedTimeMax_(-1.0),
  photonE2OverE9Max_(-1.0),
  photonDPhiMin_(-1.0),
  photonDRMin_(-1.0),
  pixelVetoOnFake_(true),
  photonPassETMin1_(VBOOL()),
  photonPassETMin2_(VBOOL()),
  photonPassAbsEtaMax_(VBOOL()),
  photonPassECALIsoMax_(VBOOL()),
  photonPassHCALIsoMax_(VBOOL()),
  photonPassHOverEMax_(VBOOL()),
  photonPassR9Max_(VBOOL()),
  photonPassTrackIsoMax_(VBOOL()),
  photonPassSigmaIetaIetaMax_(VBOOL()),
  photonPassAbsSeedTimeMax_(VBOOL()),
  photonPassE2OverE9Max_(VBOOL()),
  photonPassPreselection_(VBOOL()),
  photonPairsPassingDPhiMin_(PHOTON_PAIR_MAP()),
  photonPairsPassingDRMin_(PHOTON_PAIR_MAP()),
  photonPairsPassingAsymmetricETMin_(PHOTON_ET_PAIR_MAP()),
  evtPassDPhiMin_(false),
  evtPassDRMin_(false),
  evtPassAsymmetricETMin_(false),
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
			 const VDOUBLE& photonHOverE, const VDOUBLE& photonR9, 
			 const VDOUBLE& photonTrackIso, const VDOUBLE& photonSigmaIetaIeta, 
			 const VDOUBLE& photonSeedTime, const VDOUBLE& photonE2OverE9, 
			 const VDOUBLE& photonPhi, const VBOOL& photonHasPixelSeed, 
			 const double photon1ETMin, const double photon2ETMin, 
			 const double photonAbsEtaMax, const double photonECALIsoMaxPTMultiplier, 
			 const double photonECALIsoMaxConstant, 
			 const double photonHCALIsoMaxPTMultiplier,
			 const double photonHCALIsoMaxConstant, const double photonHOverEMax, 
			 const double photonR9Max, const double photonTrackIsoMaxPTMultiplier,
			 const double photonTrackIsoMaxConstant,
			 const double photonSigmaIetaIetaMax, const double photonAbsSeedTimeMax,
			 const double photonE2OverE9Max, const double photonDPhiMin, 
			 const double photonDRMin, const bool pixelVetoOnFake) :
  photonET_(photonET),
  photonEta_(photonEta),
  photonECALIso_(photonECALIso),
  photonHCALIso_(photonHCALIso),
  photonHOverE_(photonHOverE),
  photonR9_(photonR9),
  photonTrackIso_(photonTrackIso),
  photonSigmaIetaIeta_(photonSigmaIetaIeta),
  photonSeedTime_(photonSeedTime),
  photonE2OverE9_(photonE2OverE9),
  photonPhi_(photonPhi),
  photonHasPixelSeed_(photonHasPixelSeed),
  photon1ETMin_(photon1ETMin),
  photon2ETMin_(photon2ETMin),
  photonAbsEtaMax_(photonAbsEtaMax),
  photonECALIsoMaxPTMultiplier_(photonECALIsoMaxPTMultiplier),
  photonECALIsoMaxConstant_(photonECALIsoMaxConstant),
  photonHCALIsoMaxPTMultiplier_(photonHCALIsoMaxPTMultiplier),
  photonHCALIsoMaxConstant_(photonHCALIsoMaxConstant),
  photonHOverEMax_(photonHOverEMax),
  photonR9Max_(photonR9Max),
  photonTrackIsoMaxPTMultiplier_(photonTrackIsoMaxPTMultiplier),
  photonTrackIsoMaxConstant_(photonTrackIsoMaxConstant),
  photonSigmaIetaIetaMax_(photonSigmaIetaIetaMax),
  photonAbsSeedTimeMax_(photonAbsSeedTimeMax),
  photonE2OverE9Max_(photonE2OverE9Max),
  photonDPhiMin_(photonDPhiMin),
  photonDRMin_(photonDRMin),
  pixelVetoOnFake_(pixelVetoOnFake),
  photonPassETMin1_(VBOOL()),
  photonPassETMin2_(VBOOL()),
  photonPassAbsEtaMax_(VBOOL()),
  photonPassECALIsoMax_(VBOOL()),
  photonPassHCALIsoMax_(VBOOL()),
  photonPassHOverEMax_(VBOOL()),
  photonPassR9Max_(VBOOL()),
  photonPassTrackIsoMax_(VBOOL()),
  photonPassSigmaIetaIetaMax_(VBOOL()),
  photonPassAbsSeedTimeMax_(VBOOL()),
  photonPassE2OverE9Max_(VBOOL()),
  photonPassPreselection_(VBOOL()),
  photonPairsPassingDPhiMin_(PHOTON_PAIR_MAP()),
  photonPairsPassingDRMin_(PHOTON_PAIR_MAP()),
  photonPairsPassingAsymmetricETMin_(PHOTON_ET_PAIR_MAP()),
  evtPassDPhiMin_(false),
  evtPassDRMin_(false),
  evtPassAsymmetricETMin_(false),
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
  photonR9_(copy.getPhotonR9()),
  photonTrackIso_(copy.getPhotonTrackIso()),
  photonSigmaIetaIeta_(copy.getPhotonSigmaIetaIeta()),
  photonSeedTime_(copy.getPhotonSeedTime()),
  photonE2OverE9_(copy.getPhotonE2OverE9()),
  photonPhi_(copy.getPhotonPhi()),
  photonHasPixelSeed_(copy.getPhotonHasPixelSeed()),
  photon1ETMin_(copy.getPhoton1ETMin()),
  photon2ETMin_(copy.getPhoton2ETMin()),
  photonAbsEtaMax_(copy.getPhotonAbsEtaMax()),
  photonECALIsoMaxPTMultiplier_(copy.getPhotonECALIsoMaxPTMultiplier()),
  photonECALIsoMaxConstant_(copy.getPhotonECALIsoMaxConstant()),
  photonHCALIsoMaxPTMultiplier_(copy.getPhotonHCALIsoMaxPTMultiplier()),
  photonHCALIsoMaxConstant_(copy.getPhotonHCALIsoMaxConstant()),
  photonHOverEMax_(copy.getPhotonHOverEMax()),
  photonR9Max_(copy.getPhotonR9Max()),
  photonTrackIsoMaxPTMultiplier_(copy.getPhotonTrackIsoMaxPTMultiplier()),
  photonTrackIsoMaxConstant_(copy.getPhotonTrackIsoMaxConstant()),
  photonSigmaIetaIetaMax_(copy.getPhotonSigmaIetaIetaMax()),
  photonAbsSeedTimeMax_(copy.getPhotonAbsSeedTimeMax()),
  photonE2OverE9Max_(copy.getPhotonE2OverE9Max()),
  photonDPhiMin_(copy.getPhotonDPhiMin()),
  photonDRMin_(copy.getPhotonDRMin()),
  pixelVetoOnFake_(copy.getPixelVetoOnFake()),
  photonPassETMin1_(copy.getPhotonPassETMin1(false)),
  photonPassETMin2_(copy.getPhotonPassETMin2(false)),
  photonPassAbsEtaMax_(copy.getPhotonPassAbsEtaMax(false)),
  photonPassECALIsoMax_(copy.getPhotonPassECALIsoMax(false)),
  photonPassHCALIsoMax_(copy.getPhotonPassHCALIsoMax(false)),
  photonPassHOverEMax_(copy.getPhotonPassHOverEMax(false)),
  photonPassR9Max_(copy.getPhotonPassR9Max(false)),
  photonPassTrackIsoMax_(copy.getPhotonPassTrackIsoMax(false)),
  photonPassSigmaIetaIetaMax_(copy.getPhotonPassSigmaIetaIetaMax(false)),
  photonPassAbsSeedTimeMax_(copy.getPhotonPassAbsSeedTimeMax(false)),
  photonPassE2OverE9Max_(copy.getPhotonPassE2OverE9Max(false)),
  photonPassPreselection_(copy.getPhotonPassPreselection(false)),
  photonPairsPassingDPhiMin_(copy.getPhotonPairsPassingDPhiMin(false)),
  photonPairsPassingDRMin_(copy.getPhotonPairsPassingDRMin(false)),
  photonPairsPassingAsymmetricETMin_(copy.getPhotonPairsPassingAsymmetricETMin(false)),
  evtPassDPhiMin_(copy.getEvtPassDPhiMin(false)),
  evtPassDRMin_(copy.getEvtPassDRMin(false)),
  evtPassAsymmetricETMin_(copy.getEvtPassAsymmetricETMin(false)),
  evtDiEMET_(copy.getEvtDiEMET(false)),
  passingPhotons_(copy.getPassingPhotons(false)),
  photonType_(copy.getPhotonType(false)),
  category_(copy.getCategory(false)),
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
    photonR9_ = copy.getPhotonR9();
    photonTrackIso_ = copy.getPhotonTrackIso();
    photonSigmaIetaIeta_ = copy.getPhotonSigmaIetaIeta();
    photonSeedTime_ = copy.getPhotonSeedTime();
    photonE2OverE9_ = copy.getPhotonE2OverE9();
    photonPhi_ = copy.getPhotonPhi();
    photonHasPixelSeed_ = copy.getPhotonHasPixelSeed();
    photon1ETMin_ = copy.getPhoton1ETMin();
    photon2ETMin_ = copy.getPhoton2ETMin();
    photonAbsEtaMax_ = copy.getPhotonAbsEtaMax();
    photonECALIsoMaxPTMultiplier_ = copy.getPhotonECALIsoMaxPTMultiplier();
    photonECALIsoMaxConstant_ = copy.getPhotonECALIsoMaxConstant();
    photonHCALIsoMaxPTMultiplier_ = copy.getPhotonHCALIsoMaxPTMultiplier();
    photonHCALIsoMaxConstant_ = copy.getPhotonHCALIsoMaxConstant();
    photonHOverEMax_ = copy.getPhotonHOverEMax();
    photonR9Max_ = copy.getPhotonR9Max();
    photonTrackIsoMaxPTMultiplier_ = copy.getPhotonTrackIsoMaxPTMultiplier();
    photonTrackIsoMaxConstant_ = copy.getPhotonTrackIsoMaxConstant();
    photonSigmaIetaIetaMax_ = copy.getPhotonSigmaIetaIetaMax();
    photonAbsSeedTimeMax_ = copy.getPhotonAbsSeedTimeMax();
    photonE2OverE9Max_ = copy.getPhotonE2OverE9Max();
    photonDPhiMin_ = copy.getPhotonDPhiMin();
    photonDRMin_ = copy.getPhotonDRMin();
    pixelVetoOnFake_ = copy.getPixelVetoOnFake();
    photonPassETMin1_ = copy.getPhotonPassETMin1(false);
    photonPassETMin2_ = copy.getPhotonPassETMin2(false);
    photonPassAbsEtaMax_ = copy.getPhotonPassAbsEtaMax(false);
    photonPassECALIsoMax_ = copy.getPhotonPassECALIsoMax(false);
    photonPassHCALIsoMax_ = copy.getPhotonPassHCALIsoMax(false);
    photonPassHOverEMax_ = copy.getPhotonPassHOverEMax(false);
    photonPassR9Max_ = copy.getPhotonPassR9Max(false);
    photonPassTrackIsoMax_ = copy.getPhotonPassTrackIsoMax(false);
    photonPassSigmaIetaIetaMax_ = copy.getPhotonPassSigmaIetaIetaMax(false);
    photonPassAbsSeedTimeMax_ = copy.getPhotonPassAbsSeedTimeMax(false);
    photonPassE2OverE9Max_ = copy.getPhotonPassE2OverE9Max(false);
    photonPassPreselection_ = copy.getPhotonPassPreselection(false);
    photonPairsPassingDPhiMin_ = copy.getPhotonPairsPassingDPhiMin(false);
    photonPairsPassingDRMin_ = copy.getPhotonPairsPassingDRMin(false);
    photonPairsPassingAsymmetricETMin_ = copy.getPhotonPairsPassingAsymmetricETMin(false);
    evtPassDPhiMin_ = copy.getEvtPassDPhiMin(false);
    evtPassDRMin_ = copy.getEvtPassDRMin(false);
    evtPassAsymmetricETMin_ = copy.getEvtPassAsymmetricETMin(false);
    evtDiEMET_ = copy.getEvtDiEMET(false);
    passingPhotons_ = copy.getPassingPhotons(false);
    photonType_ = copy.getPhotonType(false);
    category_ = copy.getCategory(false);
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
VDOUBLE Categorizer::getPhotonR9() const { return photonR9_; }
VDOUBLE Categorizer::getPhotonTrackIso() const { return photonTrackIso_; }
VDOUBLE Categorizer::getPhotonSigmaIetaIeta() const { return photonSigmaIetaIeta_; }
VDOUBLE Categorizer::getPhotonSeedTime() const { return photonSeedTime_; }
VDOUBLE Categorizer::getPhotonE2OverE9() const { return photonE2OverE9_; }
VDOUBLE Categorizer::getPhotonPhi() const { return photonPhi_; }
VBOOL Categorizer::getPhotonHasPixelSeed() const { return photonHasPixelSeed_; }
double Categorizer::getPhoton1ETMin() const { return photon1ETMin_; }
double Categorizer::getPhoton2ETMin() const { return photon2ETMin_; }
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
double Categorizer::getPhotonR9Max() const { return photonR9Max_; }
double Categorizer::getPhotonTrackIsoMaxPTMultiplier() const
{
  return photonTrackIsoMaxPTMultiplier_;
}
double Categorizer::getPhotonTrackIsoMaxConstant() const { return photonTrackIsoMaxConstant_; }
double Categorizer::getPhotonSigmaIetaIetaMax() const { return photonSigmaIetaIetaMax_; }
double Categorizer::getPhotonAbsSeedTimeMax() const { return photonAbsSeedTimeMax_; }
double Categorizer::getPhotonE2OverE9Max() const { return photonE2OverE9Max_; }
double Categorizer::getPhotonDPhiMin() const { return photonDPhiMin_; }
double Categorizer::getPhotonDRMin() const { return photonDRMin_; }
bool Categorizer::getPixelVetoOnFake() const { return pixelVetoOnFake_; }
VBOOL Categorizer::getPhotonPassETMin1(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassETMin1(): call decideAll() first.\n");
  }
  return photonPassETMin1_;
}
VBOOL Categorizer::getPhotonPassETMin2(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassETMin2(): call decideAll() first.\n");
  }
  return photonPassETMin2_;
}
VBOOL Categorizer::getPhotonPassAbsEtaMax(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassAbsEtaMax(): call decideAll() first.\n");
  }
  return photonPassAbsEtaMax_;
}
VBOOL Categorizer::getPhotonPassECALIsoMax(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassECALIsoMax(): call decideAll() first.\n");
  }
  return photonPassECALIsoMax_;
}
VBOOL Categorizer::getPhotonPassHCALIsoMax(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassHCALIsoMax(): call decideAll() first.\n");
  }
  return photonPassHCALIsoMax_;
}
VBOOL Categorizer::getPhotonPassHOverEMax(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassHOverEMax(): call decideAll() first.\n");
  }
  return photonPassHOverEMax_;
}
VBOOL Categorizer::getPhotonPassR9Max(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassR9Max(): call decideAll() first.\n");
  }
  return photonPassR9Max_;
}
VBOOL Categorizer::getPhotonPassTrackIsoMax(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassTrackIsoMax(): call decideAll() first.\n");
  }
  return photonPassTrackIsoMax_;
}
VBOOL Categorizer::getPhotonPassSigmaIetaIetaMax(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassSigmaIetaIetaMax(): call decideAll() first.\n");
  }
  return photonPassSigmaIetaIetaMax_;
}
VBOOL Categorizer::getPhotonPassAbsSeedTimeMax(const bool throwException) const
{
  if (throwException && !decided_) {
    throw STRING("Error in getPhotonPassAbsSeedTimeMax(): call decideAll() first.\n");
  }
  return photonPassAbsSeedTimeMax_;
}
VBOOL Categorizer::getPhotonPassE2OverE9Max(const bool throwException) const
{
  if (throwException && !decided_) throw STRING("Error in getPhotonPassE2OverE9Max(): call decideAll() first.\n");
  return photonPassE2OverE9Max_;
}
VBOOL Categorizer::getPhotonPassPreselection(const bool throwException) const
{
  if (throwException && !foundPassingPhotons_) {
    throw STRING("Error in getPhotonPassPreselection(): call findPassingPhotons() first.\n");
  }
  return photonPassPreselection_;
}
PHOTON_PAIR_MAP Categorizer::getPhotonPairsPassingDPhiMin(const bool throwException) const
{
  if (throwException && !foundPassingPhotons_) {
    throw STRING("Error in getPhotonPairsPassingDPhiMin(): call findPassingPhotons() first.\n");
  }
  return photonPairsPassingDPhiMin_;
}
PHOTON_PAIR_MAP Categorizer::getPhotonPairsPassingDRMin(const bool throwException) const
{
  if (throwException && !foundPassingPhotons_) {
    throw STRING("Error in getPhotonPairsPassingDRMin(): call findPassingPhotons() first.\n");
  }
  return photonPairsPassingDRMin_;
}
PHOTON_ET_PAIR_MAP
Categorizer::getPhotonPairsPassingAsymmetricETMin(const bool throwException) const
{
  if (throwException && !foundPassingPhotons_) {
    STRINGSTREAM err;
    err << "Error in getPhotonPairsPassingAsymmetricETMin(): call findPassingPhotons() first.\n";
    throw err.str();
  }
  return photonPairsPassingAsymmetricETMin_;
}
bool Categorizer::getEvtPassDPhiMin(const bool throwException) const
{
  if (throwException && !categorized_) {
    throw STRING("Error in getEvtPassDPhiMin(): call classify() first.\n");
  }
  return evtPassDPhiMin_;
}
bool Categorizer::getEvtPassDRMin(const bool throwException) const
{
  if (throwException && !categorized_) {
    throw STRING("Error in getEvtPassDRMin(): call classify() first.\n");
  }
  return evtPassDRMin_;
}
bool Categorizer::getEvtPassAsymmetricETMin(const bool throwException) const
{
  if (throwException && !categorized_) {
    throw STRING("Error in getEvtPassAsymmetricETMin(): call classify() first.\n");
  }
  return evtPassAsymmetricETMin_;
}
bool Categorizer::getEvtPassAsymmetricETMin(const unsigned int iPhoton1, 
					    const unsigned int iPhoton2) const
{
  bool pass = false;
  if ((photonPassETMin1_[iPhoton1] && photonPassETMin2_[iPhoton2]) || 
      (photonPassETMin1_[iPhoton2] && photonPassETMin2_[iPhoton1])) pass = true;
  return pass;
}
double Categorizer::getEvtDiEMET(const bool throwException) const
{
  if (throwException && !categorized_) {
    throw STRING("Error in getEvtDiEMET(): call classify() first.\n");
  }
  return evtDiEMET_;
}
VINT Categorizer::getPassingPhotons(const bool throwException) const
{
  if (throwException && !categorized_) {
    throw STRING("Error in getPassingPhotons(): call classify() first.\n");
  }
  return passingPhotons_;
}
VINT Categorizer::getPhotonType(const bool throwException) const
{
  if (throwException && !categorized_) {
    throw STRING("Error in getPhotonType(): call classify() first.\n");
  }
  return photonType_;
}
int Categorizer::getCategory(const bool throwException) const
{
  if (throwException && !categorized_) {
    throw STRING("Error in getCategory(): call classify() first.\n");
  }
  return category_;
}
int Categorizer::getCategory(const unsigned int iPhoton1, const unsigned int iPhoton2) const
{
  int category = FAIL;
  if ((photonType_[iPhoton1] == G) && (photonType_[iPhoton2] == G)) category = GG;
  if (((photonType_[iPhoton1] == G) && (photonType_[iPhoton2] == E)) || 
      ((photonType_[iPhoton1] == E) && (photonType_[iPhoton2] == G))) category = EG;
  if ((photonType_[iPhoton1] == E) && (photonType_[iPhoton2] == E)) category = EE;
  if ((photonType_[iPhoton1] == F) && (photonType_[iPhoton2] == F)) category = FF;
  return category;
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
void Categorizer::setPhotonR9(const VDOUBLE& photonR9)
{
  photonR9_ = photonR9;
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
void Categorizer::setPhoton1ETMin(const double photon1ETMin) { photon1ETMin_ = photon1ETMin; }
void Categorizer::setPhoton2ETMin(const double photon2ETMin) { photon2ETMin_ = photon2ETMin; }
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
void Categorizer::setPhotonR9Max(const double photonR9Max)
{
  photonR9Max_ = photonR9Max;
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
void Categorizer::setPhotonDRMin(const double photonDRMin) { photonDRMin_ = photonDRMin; }
void Categorizer::setPixelVetoOnFake(const bool pixelVetoOnFake)
{
  pixelVetoOnFake_ = pixelVetoOnFake;
}

void Categorizer::decideAll()
{
  try {
    decideETMin1();
    decideETMin2();
    decideAbsEtaMax();
    decideECALIsoMax();
    decideHCALIsoMax();
    decideHOverEMax();
    decideR9Max();
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
  for (VBOOL_IT iPassETMin2 = photonPassETMin2_.begin(); iPassETMin2 != photonPassETMin2_.end(); 
       ++iPassETMin2) {
    const unsigned int i = iPassETMin2 - photonPassETMin2_.begin();
    photonPassPreselection_.push_back(*iPassETMin2 && photonPassAbsEtaMax_[i] && 
				      photonPassECALIsoMax_[i] && photonPassHCALIsoMax_[i] && 
				      photonPassHOverEMax_[i] && photonPassR9Max_[i] && 
				      photonPassAbsSeedTimeMax_[i] && photonPassE2OverE9Max_[i]);
  }
  foundPassingPhotons_ = true;
}

void Categorizer::classify()
{
  //sanity check
  if (!foundPassingPhotons_) {
    throw STRING("Error in classify(): call findPassingPhotons() first.\n");
  }

  //sort the preselected photons by type
  photonType_.clear();
  VUINT categorizedPhotons;
  for (VBOOL_IT iPassETMin1 = photonPassETMin1_.begin(); iPassETMin1 != photonPassETMin1_.end(); 
       ++iPassETMin1) {
    const unsigned int i = iPassETMin1 - photonPassETMin1_.begin();
    if (photonPassPreselection_[i]) {
      if (photonPassSigmaIetaIetaMax_[i] && 
	  photonPassTrackIsoMax_[i]) {
	if (photonHasPixelSeed_[i]) photonType_.push_back(E);
	else photonType_.push_back(G);
	categorizedPhotons.push_back(i);
      }
      else if (!photonHasPixelSeed_[i] || 
	       (photonHasPixelSeed_[i] && !pixelVetoOnFake_)) {
	photonType_.push_back(F);
	categorizedPhotons.push_back(i);
      }
      else photonType_.push_back(FAIL);
    }
    else photonType_.push_back(FAIL);
  }

  /*for each viable (e.g. ee, gg, ff, or eg) pair of categorized photons, determine if it passes 
    dR, dPhi, or asymmetric ET cuts*/
  photonPairsPassingDPhiMin_.clear();
  photonPairsPassingDRMin_.clear();
  photonPairsPassingAsymmetricETMin_.clear();
  // std::cout << "Size of categorized photons: " << categorizedPhotons.size() << std::endl;
  for (VUINT_IT iCategorizedPhoton1 = categorizedPhotons.begin(); 
       iCategorizedPhoton1 != categorizedPhotons.end(); ++iCategorizedPhoton1) {
    for (VUINT_IT iCategorizedPhoton2 = iCategorizedPhoton1 + 1; 
	 iCategorizedPhoton2 != categorizedPhotons.end(); ++iCategorizedPhoton2) {
      // std::cout << "photon 1: " << *iCategorizedPhoton1 << ", photon 2: " << *iCategorizedPhoton2;
      // std::cout << std::endl;
      if (getCategory(*iCategorizedPhoton1, *iCategorizedPhoton2) != FAIL) {
	// std::cout << "category not FAIL\n";
	const double dPhi = 
	  fabs(deltaPhi(photonPhi_[*iCategorizedPhoton1], photonPhi_[*iCategorizedPhoton2]));
	// std::cout << "dPhi = " << dPhi << std::endl;
	// std::cout << "dPhiMin = " << photonDPhiMin_ << std::endl;
	if (dPhi > photonDPhiMin_) {
	  photonPairsPassingDPhiMin_[std::make_pair(*iCategorizedPhoton1, *iCategorizedPhoton2)] = 
	    dPhi;
	}
	const double dR = 
	  deltaR(photonEta_[*iCategorizedPhoton1], photonPhi_[*iCategorizedPhoton1], 
		 photonEta_[*iCategorizedPhoton2], photonPhi_[*iCategorizedPhoton2]);
	if (dR > photonDRMin_) {
	  photonPairsPassingDRMin_[std::make_pair(*iCategorizedPhoton1, *iCategorizedPhoton2)] = dR;
	}
	if (getEvtPassAsymmetricETMin(*iCategorizedPhoton1, *iCategorizedPhoton2)) {
	  photonPairsPassingAsymmetricETMin_[std::make_pair(*iCategorizedPhoton1, 
							    *iCategorizedPhoton2)] = 
	    std::make_pair(photonET_[*iCategorizedPhoton1], photonET_[*iCategorizedPhoton2]);
	}
      }
    }
  }

  //find the pair with the highest ET photons that passes dPhi, dR, and asymmetric ET
  double photon1ETMax = -1.0;
  double photon2ETMax = -1.0;
  PHOTON_PAIR_MAP_IT iMaxETPair = photonPairsPassingDRMin_.end();
  for (PHOTON_PAIR_MAP_IT iPairPassingDPhiMin = photonPairsPassingDPhiMin_.begin(); 
       iPairPassingDPhiMin != photonPairsPassingDPhiMin_.end(); ++iPairPassingDPhiMin) {
    // std::cout << "Pair " << iPairPassingDPhiMin->first.first << ", ";
    // std::cout << iPairPassingDPhiMin->first.second << std::endl;
    PHOTON_PAIR_MAP_IT iPairPassingDR = photonPairsPassingDRMin_.find(iPairPassingDPhiMin->first);
    // std::cout << "Found a dR match?  ";
    // std::cout << ((iPairPassingDR != photonPairsPassingDRMin_.end()) ? "Yes" : "No") << std::endl;
    PHOTON_ET_PAIR_MAP_IT iPairPassingAsymmetricET = 
      photonPairsPassingAsymmetricETMin_.find(iPairPassingDPhiMin->first);
    // std::cout << "Found an ET match?  ";
    // std::cout << ((iPairPassingAsymmetricET != photonPairsPassingAsymmetricETMin_.end()) ? "Yes" : "No") << std::endl;
    if ((iPairPassingDR != photonPairsPassingDRMin_.end()) && 
	(iPairPassingAsymmetricET != photonPairsPassingAsymmetricETMin_.end())) {
      // std::cout << "This pair: " << iPairPassingAsymmetricET->first.first << " ";
      // std::cout << iPairPassingAsymmetricET->first.second << " ";
      // std::cout << iPairPassingAsymmetricET->second.first << " ";
      // std::cout << iPairPassingAsymmetricET->second.second << std::endl;
      double* pPhotonMax = NULL;
      double* pPhotonMin = NULL;
      if (std::max(photon1ETMax, photon2ETMax) == photon1ETMax) {
	pPhotonMax = &photon1ETMax;
	pPhotonMin = &photon2ETMax;
      }
      else {
	pPhotonMax = &photon2ETMax;
	pPhotonMin = &photon1ETMax;
      }
      double thisPhotonMax = -1.0;
      double thisPhotonMin = -1.0;
      if (std::max(photonET_[iPairPassingDR->first.first], 
		   photonET_[iPairPassingDR->first.second]) == 
	  photonET_[iPairPassingDR->first.first]) {
	thisPhotonMax = photonET_[iPairPassingDR->first.first];
	thisPhotonMin = photonET_[iPairPassingDR->first.second];
      }
      else {
	thisPhotonMax = photonET_[iPairPassingDR->first.second];
	thisPhotonMin = photonET_[iPairPassingDR->first.first];
      }
      if (thisPhotonMax > *pPhotonMax) {
      	*pPhotonMax = thisPhotonMax;
      	*pPhotonMin = thisPhotonMin;
	iMaxETPair = iPairPassingDR;
      }
      else if ((thisPhotonMax == *pPhotonMax) && (thisPhotonMin > *pPhotonMin)) {
      	*pPhotonMin = thisPhotonMin;
	iMaxETPair = iPairPassingDR;
      }
      // std::cout << "Max ET pair: " << iMaxETPair->first.first << " " << iMaxETPair->first.second;
      // std::cout << " " << iMaxETPair->second << std::endl;
    }
  }

  //now that the pair is chosen, set the event category and some event variables
  evtPassDPhiMin_ = false;
  evtPassDRMin_ = false;
  evtPassAsymmetricETMin_ = false;
  passingPhotons_[0] = -1;
  passingPhotons_[1] = -1;
  evtDiEMET_ = -1.0;
  category_ = FAIL;
  if (iMaxETPair != photonPairsPassingDRMin_.end()) {
    evtPassDPhiMin_ = true;
    evtPassDRMin_ = true;
    evtPassAsymmetricETMin_ = true;
    passingPhotons_[0] = iMaxETPair->first.first;
    passingPhotons_[1] = iMaxETPair->first.second;
    evtDiEMET_ = sqrt(photonET_[passingPhotons_[0]]*photonET_[passingPhotons_[0]] + 
		      photonET_[passingPhotons_[1]]*photonET_[passingPhotons_[1]] + 
		      2*photonET_[passingPhotons_[0]]*photonET_[passingPhotons_[1]]*
		      cos(deltaPhi(photonPhi_[passingPhotons_[0]], 
				   photonPhi_[passingPhotons_[1]])));
    category_ = getCategory(iMaxETPair->first.first, iMaxETPair->first.second);
  }

  //finished
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
     (photonHOverE_.size() == photonR9_.size()) && 
     (photonR9_.size() == photonTrackIso_.size()) && 
     (photonTrackIso_.size() == photonSigmaIetaIeta_.size()) && 
     (photonSigmaIetaIeta_.size() == photonSeedTime_.size()) && 
     (photonSeedTime_.size() == photonE2OverE9_.size()) && 
     (photonE2OverE9_.size() == photonPhi_.size()) && 
     (photonPhi_.size() == photonHasPixelSeed_.size()) && (photon1ETMin_ >= photon2ETMin_));
  decided_ = false;
  foundPassingPhotons_ = false;
  categorized_ = false;
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

void Categorizer::decideETMin1()
{
  try { decide1DMin(photonPassETMin1_, photonET_, photon1ETMin_); }
  catch (STRING& vecSizeMismatch) { throw; }
}
void Categorizer::decideETMin2()
{
  try { decide1DMin(photonPassETMin2_, photonET_, photon2ETMin_); }
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
void Categorizer::decideR9Max()
{
  try { decide1DMax(photonPassR9Max_, photonR9_, photonR9Max_); }
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
