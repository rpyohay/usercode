#ifndef SusyCategory_h
#define SusyCategory_h

#include <map>
#include "TString.h"
#include "../../../GMSBTools/Filters/interface/Typedefs.h"

namespace susy {

  //typedefs
  typedef std::map<TString, std::map<unsigned int, int> > TAG_VALUEMAP_INT;
  typedef std::map<unsigned int, int> VALUEMAP_INT;
  typedef std::map<TString, std::map<unsigned int, bool> > TAG_VALUEMAP_BOOL;
  typedef std::map<unsigned int, bool> VALUEMAP_BOOL;
  typedef std::map<TString, int> TAG_INT;
  typedef std::map<TString, bool> TAG_BOOL;
  typedef std::map<TString, double> TAG_DOUBLE;
  typedef std::map<TString, std::map<unsigned int, double> > TAG_VALUEMAP_DOUBLE;
  typedef std::map<unsigned int, double> VALUEMAP_DOUBLE;

  class Category {
  public:

    //constructors and destructor
    Category();
    Category(const Category&);
    ~Category();

    //assignment operator
    Category& operator=(const Category&);

    //getters
    int getPhotonType(const TString&, const unsigned int) const;
    bool getPassETMin1(const TString&, const unsigned int) const;
    bool getPassETMin2(const TString&, const unsigned int) const;
    bool getPassAbsEtaMax(const TString&, const unsigned int) const;
    bool getPassECALIsoMax(const TString&, const unsigned int) const;
    bool getPassPUSubtractedECALIsoMax(const TString&, const unsigned int) const;
    bool getPassHCALIsoMax(const TString&, const unsigned int) const;
    bool getPassPUSubtractedHCALIsoMax(const TString&, const unsigned int) const;
    bool getPassHOverEMax(const TString&, const unsigned int) const;
    bool getPassR9Max(const TString&, const unsigned int) const;
    bool getPassR9Min(const TString&, const unsigned int) const;
    bool getPassTrackIsoMax(const TString&, const unsigned int) const;
    bool getPassCombinedIsoMax(const TString&, const unsigned int) const;
    bool getPassFakeCombinedIsoMax(const TString&, const unsigned int) const;
    bool getPassSigmaIetaIetaMax(const TString&, const unsigned int) const;
    bool getPassHLTSigmaIetaIetaMax(const TString&, const unsigned int) const;
    bool getPassAbsSeedTimeMax(const TString&, const unsigned int) const;
    bool getPassE2OverE9Max(const TString&, const unsigned int) const;
    bool getHasPixelSeed(const TString&, const unsigned int) const;
    bool getPassPreselection(const TString&, const unsigned int) const;
    bool getIsDeciding(const TString&, const unsigned int) const;
    int getEventCategory(const TString&) const;
    bool getPassDPhiMin(const TString&) const;
    bool getPassDRMin(const TString&) const;
    bool getPassAsymmetricETMin(const TString&) const;
    double getEvtDiEMET(const TString&) const;
    double getEvtInvMass(const TString&) const;
    double getPhotonSeedTime(const TString&, const unsigned int) const;
    double getPhotonE2OverE9(const TString&, const unsigned int) const;
    int getPhotonSeedIeta(const TString&, const unsigned int) const;
    bool getPassGoodPV(const TString&) const;

    //setters
    void setPhotonType(const TString&, const unsigned int, const int);
    void setPassETMin1(const TString&, const unsigned int, const bool);
    void setPassETMin2(const TString&, const unsigned int, const bool);
    void setPassAbsEtaMax(const TString&, const unsigned int, const bool);
    void setPassECALIsoMax(const TString&, const unsigned int, const bool);
    void setPassPUSubtractedECALIsoMax(const TString&, const unsigned int, const bool);
    void setPassHCALIsoMax(const TString&, const unsigned int, const bool);
    void setPassPUSubtractedHCALIsoMax(const TString&, const unsigned int, const bool);
    void setPassHOverEMax(const TString&, const unsigned int, const bool);
    void setPassR9Max(const TString&, const unsigned int, const bool);
    void setPassR9Min(const TString&, const unsigned int, const bool);
    void setPassTrackIsoMax(const TString&, const unsigned int, const bool);
    void setPassCombinedIsoMax(const TString&, const unsigned int, const bool);
    void setPassFakeCombinedIsoMax(const TString&, const unsigned int, const bool);
    void setPassSigmaIetaIetaMax(const TString&, const unsigned int, const bool);
    void setPassHLTSigmaIetaIetaMax(const TString&, const unsigned int, const bool);
    void setPassAbsSeedTimeMax(const TString&, const unsigned int, const bool);
    void setPassE2OverE9Max(const TString&, const unsigned int, const bool);
    void setHasPixelSeed(const TString&, const unsigned int, const bool);
    void setPassPreselection(const TString&, const unsigned int, const bool);
    void setIsDeciding(const TString&, const unsigned int, const bool);
    void setEventCategory(const TString&, const int);
    void setPassDPhiMin(const TString&, const bool);
    void setPassDRMin(const TString&, const bool);
    void setPassAsymmetricETMin(const TString&, const bool);
    void setEvtDiEMET(const TString&, const double);
    void setEvtInvMass(const TString&, const double);
    void setPhotonSeedTime(const TString&, const unsigned int, const double);
    void setPhotonE2OverE9(const TString&, const unsigned int, const double);
    void setPhotonSeedIeta(const TString&, const unsigned int, const int);
    void setPassGoodPV(const TString&, const bool);

    //get the entire map
    const TAG_VALUEMAP_INT* getPhotonType() const;
    const TAG_VALUEMAP_BOOL* getPassETMin1() const;
    const TAG_VALUEMAP_BOOL* getPassETMin2() const;
    const TAG_VALUEMAP_BOOL* getPassAbsEtaMax() const;
    const TAG_VALUEMAP_BOOL* getPassECALIsoMax() const;
    const TAG_VALUEMAP_BOOL* getPassPUSubtractedECALIsoMax() const;
    const TAG_VALUEMAP_BOOL* getPassHCALIsoMax() const;
    const TAG_VALUEMAP_BOOL* getPassPUSubtractedHCALIsoMax() const;
    const TAG_VALUEMAP_BOOL* getPassHOverEMax() const;
    const TAG_VALUEMAP_BOOL* getPassR9Max() const;
    const TAG_VALUEMAP_BOOL* getPassR9Min() const;
    const TAG_VALUEMAP_BOOL* getPassTrackIsoMax() const;
    const TAG_VALUEMAP_BOOL* getPassCombinedIsoMax() const;
    const TAG_VALUEMAP_BOOL* getPassFakeCombinedIsoMax() const;
    const TAG_VALUEMAP_BOOL* getPassSigmaIetaIetaMax() const;
    const TAG_VALUEMAP_BOOL* getPassHLTSigmaIetaIetaMax() const;
    const TAG_VALUEMAP_BOOL* getPassAbsSeedTimeMax() const;
    const TAG_VALUEMAP_BOOL* getPassE2OverE9Max() const;
    const TAG_VALUEMAP_BOOL* getHasPixelSeed() const;
    const TAG_VALUEMAP_BOOL* getPassPreselection() const;
    const TAG_VALUEMAP_BOOL* getIsDeciding() const;
    const TAG_INT* getEventCategory() const;
    const TAG_BOOL* getPassDPhiMin() const;
    const TAG_BOOL* getPassDRMin() const;
    const TAG_BOOL* getPassAsymmetricETMin() const;
    const TAG_DOUBLE* getEvtDiEMET() const;
    const TAG_DOUBLE* getEvtInvMass() const;
    const TAG_VALUEMAP_DOUBLE* getPhotonSeedTime() const;
    const TAG_VALUEMAP_DOUBLE* getPhotonE2OverE9() const;
    const TAG_VALUEMAP_INT* getPhotonSeedIeta() const;
    const TAG_BOOL* getPassGoodPV() const;

    //set the entire map
    void setPhotonType(const TAG_VALUEMAP_INT&);
    void setPassETMin1(const TAG_VALUEMAP_BOOL&);
    void setPassETMin2(const TAG_VALUEMAP_BOOL&);
    void setPassAbsEtaMax(const TAG_VALUEMAP_BOOL&);
    void setPassECALIsoMax(const TAG_VALUEMAP_BOOL&);
    void setPassPUSubtractedECALIsoMax(const TAG_VALUEMAP_BOOL&);
    void setPassHCALIsoMax(const TAG_VALUEMAP_BOOL&);
    void setPassPUSubtractedHCALIsoMax(const TAG_VALUEMAP_BOOL&);
    void setPassHOverEMax(const TAG_VALUEMAP_BOOL&);
    void setPassR9Max(const TAG_VALUEMAP_BOOL&);
    void setPassR9Min(const TAG_VALUEMAP_BOOL&);
    void setPassTrackIsoMax(const TAG_VALUEMAP_BOOL&);
    void setPassCombinedIsoMax(const TAG_VALUEMAP_BOOL&);
    void setPassFakeCombinedIsoMax(const TAG_VALUEMAP_BOOL&);
    void setPassSigmaIetaIetaMax(const TAG_VALUEMAP_BOOL&);
    void setPassHLTSigmaIetaIetaMax(const TAG_VALUEMAP_BOOL&);
    void setPassAbsSeedTimeMax(const TAG_VALUEMAP_BOOL&);
    void setPassE2OverE9Max(const TAG_VALUEMAP_BOOL&);
    void setHasPixelSeed(const TAG_VALUEMAP_BOOL&);
    void setPassPreselection(const TAG_VALUEMAP_BOOL&);
    void setIsDeciding(const TAG_VALUEMAP_BOOL&);
    void setEventCategory(const TAG_INT&);
    void setPassDPhiMin(const TAG_BOOL&);
    void setPassDRMin(const TAG_BOOL&);
    void setPassAsymmetricETMin(const TAG_BOOL&);
    void setEvtDiEMET(const TAG_DOUBLE&);
    void setEvtInvMass(const TAG_DOUBLE&);
    void setPhotonSeedTime(const TAG_VALUEMAP_DOUBLE&);
    void setPhotonE2OverE9(const TAG_VALUEMAP_DOUBLE&);
    void setPhotonSeedIeta(const TAG_VALUEMAP_INT&);
    void setPassGoodPV(const TAG_BOOL&);

    //clear all maps
    void reset();

  private:

    /*each map is between a string containing the photon collection name and the corresponding 
      value map*/
    TAG_VALUEMAP_INT photonType_;
    TAG_VALUEMAP_BOOL passETMin1_;
    TAG_VALUEMAP_BOOL passETMin2_;
    TAG_VALUEMAP_BOOL passAbsEtaMax_;
    TAG_VALUEMAP_BOOL passECALIsoMax_;
    TAG_VALUEMAP_BOOL passPUSubtractedECALIsoMax_;
    TAG_VALUEMAP_BOOL passHCALIsoMax_;
    TAG_VALUEMAP_BOOL passPUSubtractedHCALIsoMax_;
    TAG_VALUEMAP_BOOL passHOverEMax_;
    TAG_VALUEMAP_BOOL passR9Max_;
    TAG_VALUEMAP_BOOL passR9Min_;
    TAG_VALUEMAP_BOOL passTrackIsoMax_;
    TAG_VALUEMAP_BOOL passCombinedIsoMax_;
    TAG_VALUEMAP_BOOL passFakeCombinedIsoMax_;
    TAG_VALUEMAP_BOOL passSigmaIetaIetaMax_;
    TAG_VALUEMAP_BOOL passHLTSigmaIetaIetaMax_;
    TAG_VALUEMAP_BOOL passAbsSeedTimeMax_;
    TAG_VALUEMAP_BOOL passE2OverE9Max_;
    TAG_VALUEMAP_BOOL hasPixelSeed_;
    TAG_VALUEMAP_BOOL passPreselection_;
    TAG_VALUEMAP_BOOL isDeciding_;
    TAG_INT eventCategory_;
    TAG_BOOL passDPhiMin_;
    TAG_BOOL passDRMin_;
    TAG_BOOL passAsymmetricETMin_;
    TAG_DOUBLE evtDiEMET_;
    TAG_DOUBLE evtInvMass_;
    TAG_VALUEMAP_DOUBLE photonSeedTime_;
    TAG_VALUEMAP_DOUBLE photonE2OverE9_;
    TAG_VALUEMAP_INT photonSeedIeta_;
    TAG_BOOL passGoodPV_;

    //generic getter for a value map
    template<typename T, typename U>
      typename U::const_iterator valueMap(const T& map, const TString& tag, 
					  const unsigned int photonIndex) const
      {
	typename T::const_iterator iMap = map.find(tag);
	if (iMap == map.end()) {
	  STRINGSTREAM err;
	  err << "Error: no tag \"" << tag << "\" found.\n";
	  throw err.str();
	}
	typename U::const_iterator iPhoton = (*iMap).second.find(photonIndex);
	if (iPhoton == (*iMap).second.end()) {
	  STRINGSTREAM err;
	  err << "Error: no photon " << photonIndex << " with tag \"" << tag << "\" found.\n";
	  throw err.str();
	}
	return iPhoton;
      }

    //generic getter for a primitive type
    template<typename T>
      typename T::const_iterator primitiveType(const T& map, const TString& tag) const
      {
	typename T::const_iterator iMap = map.find(tag);
	if (iMap == map.end()) {
	  STRINGSTREAM err;
	  err << "Error: no tag \"" << tag << "\" found.\n";
	  throw err.str();
	}
	return iMap;
      }

    //generic setter for a value map
    template<typename T, typename U, typename V>
      void setValueMap(T& map, const TString& tag, const unsigned int photonIndex, const V val)
      {
	typename T::iterator i = map.find(tag);
	if (i != map.end()) {
	  typename U::iterator j = (*i).second.find(photonIndex);
	  if (j != (*i).second.end()) (*j).second = val;
	  else (*i).second[photonIndex] = val;
	}
	else {
	  map[tag] = U();
	  map[tag][photonIndex] = val;
	}
      }

    //generic setter for a primitive type
    template<typename T, typename V>
      void setPrimitiveType(T& map, const TString& tag, const V val)
      {
	typename T::iterator i = map.find(tag);
	if (i != map.end()) (*i).second = val;
	else map[tag] = val;
      }
  };
}

#endif
