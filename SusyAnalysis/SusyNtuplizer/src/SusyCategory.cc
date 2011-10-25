#include "SusyCategory.h"
#include "../../../GMSBTools/Filters/interface/Categorizer.h"

susy::Category::Category() :
  photonType_(susy::TAG_VALUEMAP_INT()),
  passETMin1_(susy::TAG_VALUEMAP_BOOL()),
  passETMin2_(susy::TAG_VALUEMAP_BOOL()),
  passAbsEtaMax_(susy::TAG_VALUEMAP_BOOL()),
  passECALIsoMax_(susy::TAG_VALUEMAP_BOOL()),
  passPUSubtractedECALIsoMax_(susy::TAG_VALUEMAP_BOOL()),
  passHCALIsoMax_(susy::TAG_VALUEMAP_BOOL()),
  passPUSubtractedHCALIsoMax_(susy::TAG_VALUEMAP_BOOL()),
  passHOverEMax_(susy::TAG_VALUEMAP_BOOL()),
  passR9Max_(susy::TAG_VALUEMAP_BOOL()),
  passTrackIsoMax_(susy::TAG_VALUEMAP_BOOL()),
  passCombinedIsoMax_(susy::TAG_VALUEMAP_BOOL()),
  passFakeCombinedIsoMax_(susy::TAG_VALUEMAP_BOOL()),
  passSigmaIetaIetaMax_(susy::TAG_VALUEMAP_BOOL()),
  passAbsSeedTimeMax_(susy::TAG_VALUEMAP_BOOL()),
  passE2OverE9Max_(susy::TAG_VALUEMAP_BOOL()),
  hasPixelSeed_(susy::TAG_VALUEMAP_BOOL()),
  passPreselection_(susy::TAG_VALUEMAP_BOOL()),
  isDeciding_(susy::TAG_VALUEMAP_BOOL()),
  eventCategory_(susy::TAG_INT()),
  passDPhiMin_(susy::TAG_BOOL()),
  passDRMin_(susy::TAG_BOOL()),
  passAsymmetricETMin_(susy::TAG_BOOL()),
  evtDiEMET_(susy::TAG_DOUBLE()),
  evtInvMass_(susy::TAG_DOUBLE()),
  photonSeedTime_(susy::TAG_VALUEMAP_DOUBLE()),
  photonE2OverE9_(susy::TAG_VALUEMAP_DOUBLE()),
  photonSeedIeta_(susy::TAG_VALUEMAP_INT()) {}

susy::Category::Category(const susy::Category& other) :
  photonType_(*(other.getPhotonType())),
  passETMin1_(*(other.getPassETMin1())),
  passETMin2_(*(other.getPassETMin2())),
  passAbsEtaMax_(*(other.getPassAbsEtaMax())),
  passECALIsoMax_(*(other.getPassECALIsoMax())),
  passPUSubtractedECALIsoMax_(*(other.getPassPUSubtractedECALIsoMax())),
  passHCALIsoMax_(*(other.getPassHCALIsoMax())),
  passPUSubtractedHCALIsoMax_(*(other.getPassPUSubtractedHCALIsoMax())),
  passHOverEMax_(*(other.getPassHOverEMax())),
  passR9Max_(*(other.getPassR9Max())),
  passTrackIsoMax_(*(other.getPassTrackIsoMax())),
  passCombinedIsoMax_(*(other.getPassCombinedIsoMax())),
  passFakeCombinedIsoMax_(*(other.getPassFakeCombinedIsoMax())),
  passSigmaIetaIetaMax_(*(other.getPassSigmaIetaIetaMax())),
  passAbsSeedTimeMax_(*(other.getPassAbsSeedTimeMax())),
  passE2OverE9Max_(*(other.getPassE2OverE9Max())),
  hasPixelSeed_(*(other.getHasPixelSeed())),
  passPreselection_(*(other.getPassPreselection())),
  isDeciding_(*(other.getIsDeciding())),
  eventCategory_(*(other.getEventCategory())),
  passDPhiMin_(*(other.getPassDPhiMin())),
  passDRMin_(*(other.getPassDRMin())),
  passAsymmetricETMin_(*(other.getPassAsymmetricETMin())),
  evtDiEMET_(*(other.getEvtDiEMET())),
  evtInvMass_(*(other.getEvtInvMass())),
  photonSeedTime_(*(other.getPhotonSeedTime())),
  photonE2OverE9_(*(other.getPhotonE2OverE9())),
  photonSeedIeta_(*(other.getPhotonSeedIeta())) {}

susy::Category::~Category() { reset(); }

susy::Category& susy::Category::operator=(const susy::Category& other)
{
  if (this != &other) {
    photonType_ = *(other.getPhotonType());
    passETMin1_ = *(other.getPassETMin1());
    passETMin2_ = *(other.getPassETMin2());
    passAbsEtaMax_ = *(other.getPassAbsEtaMax());
    passECALIsoMax_ = *(other.getPassECALIsoMax());
    passPUSubtractedECALIsoMax_ = *(other.getPassPUSubtractedECALIsoMax());
    passHCALIsoMax_ = *(other.getPassHCALIsoMax());
    passPUSubtractedHCALIsoMax_ = *(other.getPassPUSubtractedHCALIsoMax());
    passHOverEMax_ = *(other.getPassHOverEMax());
    passR9Max_ = *(other.getPassR9Max());
    passTrackIsoMax_ = *(other.getPassTrackIsoMax());
    passCombinedIsoMax_ = *(other.getPassCombinedIsoMax());
    passFakeCombinedIsoMax_ = *(other.getPassFakeCombinedIsoMax());
    passSigmaIetaIetaMax_ = *(other.getPassSigmaIetaIetaMax());
    passAbsSeedTimeMax_ = *(other.getPassAbsSeedTimeMax());
    passE2OverE9Max_ = *(other.getPassE2OverE9Max());
    hasPixelSeed_ = *(other.getHasPixelSeed());
    passPreselection_ = *(other.getPassPreselection());
    isDeciding_ = *(other.getIsDeciding());
    eventCategory_ = *(other.getEventCategory());
    passDPhiMin_ = *(other.getPassDPhiMin());
    passDRMin_ = *(other.getPassDRMin());
    passAsymmetricETMin_ = *(other.getPassAsymmetricETMin());
    evtDiEMET_ = *(other.getEvtDiEMET());
    evtInvMass_ = *(other.getEvtInvMass());
    photonSeedTime_ = *(other.getPhotonSeedTime());
    photonE2OverE9_ = *(other.getPhotonE2OverE9());
    photonSeedIeta_ = *(other.getPhotonSeedIeta());
  }
  return *this;
}

int susy::Category::getPhotonType(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_INT, susy::VALUEMAP_INT>(photonType_, tag, 
								  photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassETMin1(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passETMin1_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassETMin2(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passETMin2_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassAbsEtaMax(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passAbsEtaMax_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassECALIsoMax(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passECALIsoMax_, tag, 
							photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassPUSubtractedECALIsoMax(const TString& tag, 
						   const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passPUSubtractedECALIsoMax_, 
								    tag, photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassHCALIsoMax(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passHCALIsoMax_, tag, 
							photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassPUSubtractedHCALIsoMax(const TString& tag, 
						   const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passPUSubtractedHCALIsoMax_, 
								    tag, photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassHOverEMax(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passHOverEMax_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassR9Max(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passR9Max_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassTrackIsoMax(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passTrackIsoMax_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassCombinedIsoMax(const TString& tag, 
					   const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passCombinedIsoMax_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassFakeCombinedIsoMax(const TString& tag, 
					       const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passFakeCombinedIsoMax_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassSigmaIetaIetaMax(const TString& tag, 
					     const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passSigmaIetaIetaMax_, tag, 
							photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassAbsSeedTimeMax(const TString& tag, 
					   const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passAbsSeedTimeMax_, tag, 
							photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassE2OverE9Max(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passE2OverE9Max_, tag, 
							photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getHasPixelSeed(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(hasPixelSeed_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassPreselection(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(passPreselection_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getIsDeciding(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL>(isDeciding_, tag, 
								    photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

int susy::Category::getEventCategory(const TString& tag) const
{
  try { return (*primitiveType<susy::TAG_INT>(eventCategory_, tag)).second; }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassDPhiMin(const TString& tag) const
{
  try { return (*primitiveType<susy::TAG_BOOL>(passDPhiMin_, tag)).second; }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassDRMin(const TString& tag) const
{
  try { return (*primitiveType<susy::TAG_BOOL>(passDRMin_, tag)).second; }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassAsymmetricETMin(const TString& tag) const
{
  try { return (*primitiveType<susy::TAG_BOOL>(passAsymmetricETMin_, tag)).second; }
  catch (STRING& ex) { throw; }
}

double susy::Category::getEvtDiEMET(const TString& tag) const
{
  try { return (*primitiveType<susy::TAG_DOUBLE>(evtDiEMET_, tag)).second; }
  catch (STRING& ex) { throw; }
}

double susy::Category::getEvtInvMass(const TString& tag) const
{
  try { return (*primitiveType<susy::TAG_DOUBLE>(evtInvMass_, tag)).second; }
  catch (STRING& ex) { throw; }
}

double susy::Category::getPhotonSeedTime(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_DOUBLE, susy::VALUEMAP_DOUBLE>(photonSeedTime_, tag, 
									photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

double susy::Category::getPhotonE2OverE9(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_DOUBLE, susy::VALUEMAP_DOUBLE>(photonE2OverE9_, tag, 
									photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

int susy::Category::getPhotonSeedIeta(const TString& tag, const unsigned int photonIndex) const
{
  try {
    return (*valueMap<susy::TAG_VALUEMAP_INT, susy::VALUEMAP_INT>(photonSeedIeta_, tag, 
								  photonIndex)).second;
  }
  catch (STRING& ex) { throw; }
}

bool susy::Category::getPassGoodPV(const TString& tag) const
{
  try { return (*primitiveType<susy::TAG_BOOL>(passGoodPV_, tag)).second; }
  catch (STRING& ex) { throw; }
}

void susy::Category::setPhotonType(const TString& tag, const unsigned int photonIndex, 
				   const int type)
{
  setValueMap<susy::TAG_VALUEMAP_INT, susy::VALUEMAP_INT, int>(photonType_, tag, photonIndex, 
							       type);
}

void susy::Category::setPassETMin1(const TString& tag, const unsigned int photonIndex, 
				   const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passETMin1_, tag, photonIndex, 
								  pass);
}

void susy::Category::setPassETMin2(const TString& tag, const unsigned int photonIndex, 
				   const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passETMin2_, tag, photonIndex, 
								  pass);
}

void susy::Category::setPassAbsEtaMax(const TString& tag, const unsigned int photonIndex, 
				      const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passAbsEtaMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassECALIsoMax(const TString& tag, const unsigned int photonIndex, 
				       const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passECALIsoMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassPUSubtractedECALIsoMax(const TString& tag, 
						   const unsigned int photonIndex, const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passPUSubtractedECALIsoMax_, 
								  tag, photonIndex, pass);
}

void susy::Category::setPassHCALIsoMax(const TString& tag, const unsigned int photonIndex, 
				       const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passHCALIsoMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassPUSubtractedHCALIsoMax(const TString& tag, 
						   const unsigned int photonIndex, const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passPUSubtractedHCALIsoMax_, 
								  tag, photonIndex, pass);
}

void susy::Category::setPassHOverEMax(const TString& tag, const unsigned int photonIndex, 
				      const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passHOverEMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassR9Max(const TString& tag, const unsigned int photonIndex, 
				  const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passR9Max_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassTrackIsoMax(const TString& tag, const unsigned int photonIndex, 
					const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passTrackIsoMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassCombinedIsoMax(const TString& tag, const unsigned int photonIndex, 
					   const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passCombinedIsoMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassFakeCombinedIsoMax(const TString& tag, const unsigned int photonIndex, 
					       const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passFakeCombinedIsoMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassSigmaIetaIetaMax(const TString& tag, const unsigned int photonIndex, 
					     const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passSigmaIetaIetaMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassAbsSeedTimeMax(const TString& tag, const unsigned int photonIndex, 
					   const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passAbsSeedTimeMax_, tag, 
								  photonIndex, pass);
}

void susy::Category::setPassE2OverE9Max(const TString& tag, const unsigned int photonIndex, 
					const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passE2OverE9Max_, tag, 
								  photonIndex, pass);
}

void susy::Category::setHasPixelSeed(const TString& tag, const unsigned int photonIndex, 
				     const bool hasPixelSeed)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(hasPixelSeed_, tag, photonIndex, 
								  hasPixelSeed);
}

void susy::Category::setPassPreselection(const TString& tag, const unsigned int photonIndex, 
				       const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(passPreselection_, tag, 
								  photonIndex, pass);
}

void susy::Category::setIsDeciding(const TString& tag, const unsigned int photonIndex, 
				   const bool pass)
{
  setValueMap<susy::TAG_VALUEMAP_BOOL, susy::VALUEMAP_BOOL, bool>(isDeciding_, tag, 
								  photonIndex, pass);
}

void susy::Category::setEventCategory(const TString& tag, const int category)
{
  setPrimitiveType<susy::TAG_INT, int>(eventCategory_, tag, category);
}

void susy::Category::setPassDPhiMin(const TString& tag, const bool pass)
{
  setPrimitiveType<susy::TAG_BOOL, bool>(passDPhiMin_, tag, pass);
}

void susy::Category::setPassDRMin(const TString& tag, const bool pass)
{
  setPrimitiveType<susy::TAG_BOOL, bool>(passDRMin_, tag, pass);
}

void susy::Category::setPassAsymmetricETMin(const TString& tag, const bool pass)
{
  setPrimitiveType<susy::TAG_BOOL, bool>(passAsymmetricETMin_, tag, pass);
}

void susy::Category::setEvtDiEMET(const TString& tag, const double diEMET)
{
  setPrimitiveType<susy::TAG_DOUBLE, double>(evtDiEMET_, tag, diEMET);
}

void susy::Category::setEvtInvMass(const TString& tag, const double invMass)
{
  setPrimitiveType<susy::TAG_DOUBLE, double>(evtInvMass_, tag, invMass);
}

void susy::Category::setPhotonSeedTime(const TString& tag, const unsigned int photonIndex, 
				       const double seedTime)
{
  setValueMap<susy::TAG_VALUEMAP_DOUBLE, susy::VALUEMAP_DOUBLE, double>(photonSeedTime_, tag, 
									photonIndex, seedTime);
}

void susy::Category::setPhotonE2OverE9(const TString& tag, const unsigned int photonIndex, 
				       const double e2OverE9)
{
  setValueMap<susy::TAG_VALUEMAP_DOUBLE, susy::VALUEMAP_DOUBLE, double>(photonE2OverE9_, tag, 
									photonIndex, e2OverE9);
}

void susy::Category::setPhotonSeedIeta(const TString& tag, const unsigned int photonIndex, 
				       const int seedIeta)
{
  setValueMap<susy::TAG_VALUEMAP_INT, susy::VALUEMAP_INT, int>(photonSeedIeta_, tag, photonIndex, 
							       seedIeta);
}

void susy::Category::setPassGoodPV(const TString& tag, const bool pass)
{
  setPrimitiveType<susy::TAG_BOOL, bool>(passGoodPV_, tag, pass);
}

const susy::TAG_VALUEMAP_INT* susy::Category::getPhotonType() const { return &photonType_; }

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassETMin1() const { return &passETMin1_; }

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassETMin2() const { return &passETMin2_; }

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassAbsEtaMax() const { return &passAbsEtaMax_; }

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassECALIsoMax() const
{
  return &passECALIsoMax_; 
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassPUSubtractedECALIsoMax() const
{
  return &passPUSubtractedECALIsoMax_; 
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassHCALIsoMax() const
{
  return &passHCALIsoMax_;
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassPUSubtractedHCALIsoMax() const
{
  return &passPUSubtractedHCALIsoMax_; 
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassHOverEMax() const { return &passHOverEMax_; }

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassR9Max() const { return &passR9Max_; }

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassTrackIsoMax() const
{
  return &passTrackIsoMax_;
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassCombinedIsoMax() const
{
  return &passCombinedIsoMax_;
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassFakeCombinedIsoMax() const
{
  return &passFakeCombinedIsoMax_;
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassSigmaIetaIetaMax() const
{
  return &passSigmaIetaIetaMax_;
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassAbsSeedTimeMax() const
{
  return &passAbsSeedTimeMax_;
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassE2OverE9Max() const
{
  return &passE2OverE9Max_;
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getHasPixelSeed() const { return &hasPixelSeed_; }

const susy::TAG_VALUEMAP_BOOL* susy::Category::getPassPreselection() const
{
  return &passPreselection_;
}

const susy::TAG_VALUEMAP_BOOL* susy::Category::getIsDeciding() const { return &isDeciding_; }

const susy::TAG_INT* susy::Category::getEventCategory() const { return &eventCategory_; }

const susy::TAG_BOOL* susy::Category::getPassDPhiMin() const { return &passDPhiMin_; }

const susy::TAG_BOOL* susy::Category::getPassDRMin() const { return &passDRMin_; }

const susy::TAG_BOOL* susy::Category::getPassAsymmetricETMin() const { return &passAsymmetricETMin_; }

const susy::TAG_DOUBLE* susy::Category::getEvtDiEMET() const { return &evtDiEMET_; }

const susy::TAG_DOUBLE* susy::Category::getEvtInvMass() const { return &evtInvMass_; }

const susy::TAG_VALUEMAP_DOUBLE* susy::Category::getPhotonSeedTime() const
{
  return &photonSeedTime_;
}

const susy::TAG_VALUEMAP_DOUBLE* susy::Category::getPhotonE2OverE9() const
{
  return &photonE2OverE9_;
}

const susy::TAG_VALUEMAP_INT* susy::Category::getPhotonSeedIeta() const
{
  return &photonSeedIeta_;
}

const susy::TAG_BOOL* susy::Category::getPassGoodPV() const { return &passGoodPV_; }

void susy::Category::setPhotonType(const susy::TAG_VALUEMAP_INT& photonType)
{
  photonType_ = photonType;
}

void susy::Category::setPassETMin1(const susy::TAG_VALUEMAP_BOOL& passETMin1)
{
  passETMin1_ = passETMin1;
}

void susy::Category::setPassETMin2(const susy::TAG_VALUEMAP_BOOL& passETMin2)
{
  passETMin2_ = passETMin2;
}

void susy::Category::setPassAbsEtaMax(const susy::TAG_VALUEMAP_BOOL& passAbsEtaMax)
{
  passAbsEtaMax_ = passAbsEtaMax;
}

void susy::Category::setPassECALIsoMax(const susy::TAG_VALUEMAP_BOOL& passECALIsoMax)
{
  passECALIsoMax_ = passECALIsoMax;
}

void susy::Category::setPassPUSubtractedECALIsoMax(const susy::TAG_VALUEMAP_BOOL& 
						   passPUSubtractedECALIsoMax)
{
  passPUSubtractedECALIsoMax_ = passPUSubtractedECALIsoMax;
}

void susy::Category::setPassHCALIsoMax(const susy::TAG_VALUEMAP_BOOL& passHCALIsoMax)
{
  passHCALIsoMax_ = passHCALIsoMax;
}

void susy::Category::setPassPUSubtractedHCALIsoMax(const susy::TAG_VALUEMAP_BOOL& 
						   passPUSubtractedHCALIsoMax)
{
  passPUSubtractedHCALIsoMax_ = passPUSubtractedHCALIsoMax;
}

void susy::Category::setPassHOverEMax(const susy::TAG_VALUEMAP_BOOL& passHOverEMax)
{
  passHOverEMax_ = passHOverEMax;
}

void susy::Category::setPassR9Max(const susy::TAG_VALUEMAP_BOOL& passR9Max)
{
  passR9Max_ = passR9Max;
}

void susy::Category::setPassTrackIsoMax(const susy::TAG_VALUEMAP_BOOL& passTrackIsoMax)
{
  passTrackIsoMax_ = passTrackIsoMax;
}

void susy::Category::setPassCombinedIsoMax(const susy::TAG_VALUEMAP_BOOL& passCombinedIsoMax)
{
  passCombinedIsoMax_ = passCombinedIsoMax;
}

void 
susy::Category::setPassFakeCombinedIsoMax(const susy::TAG_VALUEMAP_BOOL& passFakeCombinedIsoMax)
{
  passFakeCombinedIsoMax_ = passFakeCombinedIsoMax;
}

void susy::Category::setPassSigmaIetaIetaMax(const susy::TAG_VALUEMAP_BOOL& passSigmaIetaIetaMax)
{
  passSigmaIetaIetaMax_ = passSigmaIetaIetaMax;
}

void susy::Category::setPassAbsSeedTimeMax(const susy::TAG_VALUEMAP_BOOL& passAbsSeedTimeMax)
{
  passAbsSeedTimeMax_ = passAbsSeedTimeMax;
}

void susy::Category::setPassE2OverE9Max(const susy::TAG_VALUEMAP_BOOL& passE2OverE9Max)
{
  passE2OverE9Max_ = passE2OverE9Max;
}

void susy::Category::setHasPixelSeed(const susy::TAG_VALUEMAP_BOOL& hasPixelSeed)
{
  hasPixelSeed_ = hasPixelSeed;
}

void susy::Category::setPassPreselection(const susy::TAG_VALUEMAP_BOOL& passingPhotons)
{
  passPreselection_ = passingPhotons;
}

void susy::Category::setIsDeciding(const susy::TAG_VALUEMAP_BOOL& isDeciding)
{
  isDeciding_ = isDeciding;
}

void susy::Category::setEventCategory(const susy::TAG_INT& eventCategory)
{
  eventCategory_ = eventCategory;
}

void susy::Category::setPassDPhiMin(const susy::TAG_BOOL& passDPhiMin)
{
  passDPhiMin_ = passDPhiMin;
}

void susy::Category::setPassDRMin(const susy::TAG_BOOL& passDRMin)
{
  passDRMin_ = passDRMin;
}

void susy::Category::setPassAsymmetricETMin(const susy::TAG_BOOL& passAsymmetricETMin)
{
  passAsymmetricETMin_ = passAsymmetricETMin;
}

void susy::Category::setEvtDiEMET(const susy::TAG_DOUBLE& evtDiEMET)
{
  evtDiEMET_ = evtDiEMET;
}

void susy::Category::setEvtInvMass(const susy::TAG_DOUBLE& evtInvMass)
{
  evtInvMass_ = evtInvMass;
}

void susy::Category::setPhotonSeedTime(const susy::TAG_VALUEMAP_DOUBLE& photonSeedTime)
{
  photonSeedTime_ = photonSeedTime;
}

void susy::Category::setPhotonE2OverE9(const susy::TAG_VALUEMAP_DOUBLE& photonE2OverE9)
{
  photonE2OverE9_ = photonE2OverE9;
}

void susy::Category::setPhotonSeedIeta(const susy::TAG_VALUEMAP_INT& photonSeedIeta)
{
  photonSeedIeta_ = photonSeedIeta;
}

void susy::Category::setPassGoodPV(const susy::TAG_BOOL& passGoodPV)
{
  passGoodPV_ = passGoodPV;
}

void susy::Category::reset()
{
  photonType_.clear();
  passETMin1_.clear();
  passETMin2_.clear();
  passAbsEtaMax_.clear();
  passECALIsoMax_.clear();
  passPUSubtractedECALIsoMax_.clear();
  passHCALIsoMax_.clear();
  passPUSubtractedHCALIsoMax_.clear();
  passHOverEMax_.clear();
  passR9Max_.clear();
  passTrackIsoMax_.clear();
  passCombinedIsoMax_.clear();
  passFakeCombinedIsoMax_.clear();
  passSigmaIetaIetaMax_.clear();
  passAbsSeedTimeMax_.clear();
  passE2OverE9Max_.clear();
  hasPixelSeed_.clear();
  passPreselection_.clear();
  isDeciding_.clear();
  eventCategory_.clear();
  passDPhiMin_.clear();
  passDRMin_.clear();
  passAsymmetricETMin_.clear();
  evtDiEMET_.clear();
  evtInvMass_.clear();
  photonSeedTime_.clear();
  photonE2OverE9_.clear();
  photonSeedIeta_.clear();
  passGoodPV_.clear();
}
