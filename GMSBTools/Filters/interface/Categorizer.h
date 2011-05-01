#ifndef GMSBTools_Filters_Categorizer_h
#define GMSBTools_Filters_Categorizer_h

#include "GMSBTools/Filters/interface/Typedefs.h"

//event categories
#define FAIL -1
#define GG 0
#define EG 1
#define EE 2
#define FF 3

//photon types
#define G 0
#define E 1
#define F 2

class Categorizer {

 public:

  //default constructor
  Categorizer();

  //constructor from a list of arguments
  Categorizer(const VDOUBLE&, const VDOUBLE&, const VDOUBLE&, const VDOUBLE&, const VDOUBLE&, 
	      const VDOUBLE&, const VDOUBLE&, const VDOUBLE&, const VDOUBLE&, const VDOUBLE&, 
	      const VBOOL&, const double, const double, const double, const double, 
	      const double, const double, const double, const double, const double, const double, 
	      const double, const double, const double);

  //copy constructor
  Categorizer(const Categorizer&);

  //destructor
  ~Categorizer();

  //assignment operator
  Categorizer& operator=(const Categorizer&);

  //getters
  VDOUBLE getPhotonET() const;
  VDOUBLE getPhotonEta() const;
  VDOUBLE getPhotonECALIso() const;
  VDOUBLE getPhotonHCALIso() const;
  VDOUBLE getPhotonHOverE() const;
  VDOUBLE getPhotonTrackIso() const;
  VDOUBLE getPhotonSigmaIetaIeta() const;
  VDOUBLE getPhotonSeedTime() const;
  VDOUBLE getPhotonE2OverE9() const;
  VDOUBLE getPhotonPhi() const;
  VBOOL getPhotonHasPixelSeed() const;
  double getPhotonETMin() const;
  double getPhotonAbsEtaMax() const;
  double getPhotonECALIsoMaxPTMultiplier() const;
  double getPhotonECALIsoMaxConstant() const;
  double getPhotonHCALIsoMaxPTMultiplier() const;
  double getPhotonHCALIsoMaxConstant() const;
  double getPhotonHOverEMax() const;
  double getPhotonTrackIsoMaxPTMultiplier() const;
  double getPhotonTrackIsoMaxConstant() const;
  double getPhotonSigmaIetaIetaMax() const;
  double getPhotonAbsSeedTimeMax() const;
  double getPhotonE2OverE9Max() const;
  double getPhotonDPhiMin() const;
  VBOOL getPhotonPassETMin() const;
  VBOOL getPhotonPassAbsEtaMax() const;
  VBOOL getPhotonPassECALIsoMax() const;
  VBOOL getPhotonPassHCALIsoMax() const;
  VBOOL getPhotonPassHOverEMax() const;
  VBOOL getPhotonPassTrackIsoMax() const;
  VBOOL getPhotonPassSigmaIetaIetaMax() const;
  VBOOL getPhotonPassAbsSeedTimeMax() const;
  VBOOL getPhotonPassE2OverE9Max() const;
  VBOOL getPhotonPassPreselection() const;
  bool getEvtPassDPhiMin() const;
  double getEvtDiEMET() const;
  VINT getPassingPhotons() const;
  VINT getPhotonType() const;
  int getCategory() const;

  //setters
  void setPhotonET(const VDOUBLE&);
  void setPhotonEta(const VDOUBLE&);
  void setPhotonECALIso(const VDOUBLE&);
  void setPhotonHCALIso(const VDOUBLE&);
  void setPhotonHOverE(const VDOUBLE&);
  void setPhotonTrackIso(const VDOUBLE&);
  void setPhotonSigmaIetaIeta(const VDOUBLE&);
  void setPhotonSeedTime(const VDOUBLE&);
  void setPhotonE2OverE9(const VDOUBLE&);
  void setPhotonPhi(const VDOUBLE&);
  void setPhotonHasPixelSeed(const VBOOL&);
  void setPhotonETMin(const double);
  void setPhotonAbsEtaMax(const double);
  void setPhotonECALIsoMaxPTMultiplier(const double);
  void setPhotonECALIsoMaxConstant(const double);
  void setPhotonHCALIsoMaxPTMultiplier(const double);
  void setPhotonHCALIsoMaxConstant(const double);
  void setPhotonHOverEMax(const double);
  void setPhotonTrackIsoMaxPTMultiplier(const double);
  void setPhotonTrackIsoMaxConstant(const double);
  void setPhotonSigmaIetaIetaMax(const double);
  void setPhotonAbsSeedTimeMax(const double);
  void setPhotonE2OverE9Max(const double);
  void setPhotonDPhiMin(const double);

  /*fill the vectors of booleans with true or false depending on whether the photon passes or 
    fails the cut*/
  void decideAll();

  //set photonPassPreselection() (true if it passes, false otherwise)
  void findPassingPhotons();

  /*find 2 highest ET photons passing preselection and classify event as belonging to candidate, 
    control, or no sample*/
  void classify();

  //true if all event processing is finished
  bool done() const;

 private:

  //photon quantities
  VDOUBLE photonET_;
  VDOUBLE photonEta_;
  VDOUBLE photonECALIso_;
  VDOUBLE photonHCALIso_;
  VDOUBLE photonHOverE_;
  VDOUBLE photonTrackIso_;
  VDOUBLE photonSigmaIetaIeta_;
  VDOUBLE photonSeedTime_;
  VDOUBLE photonE2OverE9_;
  VDOUBLE photonPhi_;
  VBOOL photonHasPixelSeed_;

  //cuts
  double photonETMin_;
  double photonAbsEtaMax_;
  double photonECALIsoMaxPTMultiplier_;
  double photonECALIsoMaxConstant_;
  double photonHCALIsoMaxPTMultiplier_;
  double photonHCALIsoMaxConstant_;
  double photonHOverEMax_;
  double photonTrackIsoMaxPTMultiplier_;
  double photonTrackIsoMaxConstant_;
  double photonSigmaIetaIetaMax_;
  double photonAbsSeedTimeMax_;
  double photonE2OverE9Max_;
  double photonDPhiMin_;

  //photon pass flags
  VBOOL photonPassETMin_;
  VBOOL photonPassAbsEtaMax_;
  VBOOL photonPassECALIsoMax_;
  VBOOL photonPassHCALIsoMax_;
  VBOOL photonPassHOverEMax_;
  VBOOL photonPassTrackIsoMax_;
  VBOOL photonPassSigmaIetaIetaMax_;
  VBOOL photonPassAbsSeedTimeMax_;
  VBOOL photonPassE2OverE9Max_;
  VBOOL photonPassPreselection_;

  //event pass flags
  bool evtPassDPhiMin_;

  //di-EM ET
  double evtDiEMET_;

  /*photons on which the categorization decision was based (2 highest ET photons passing 
    preselection)*/
  VINT passingPhotons_;

  //photon type (G, F, E, or FAIL)
  VINT photonType_;

  //category
  int category_;

  //true if all inputs are good, false otherwise
  bool initialized_;

  //true if all photons have been classified as passing or failing all cuts, false otherwise
  bool decided_;

  //true if event has been categorized
  bool categorized_;

  //true if photons passing preselection have been found
  bool foundPassingPhotons_;

  //true if intialized_ is true
  bool initialized() const;

  //true if decided_ is true
  bool decided() const;

  //true if categorized_ is true
  bool categorized() const;

  //true if foundPassingPhotons_ is true
  bool foundPassingPhotons() const;

  //check that input is good
  void checkInput();

  /*fill the vector of pass flags for a 2D cut as specified in the args
    @param pass vector of pass flags, 1 element per photon, for a specific cut
    @param quantity vector of values of the quantity relating to the specific cut, 1 element per 
    photon
    @param cut0 y-intercept of the cut
    @param cut1 slope of the cut (x variable is photon ET)
    @param min true if the cut should be decided val > cut, false if it should be decided val < cut
    @param abs true if the val should be fabs(val), false otherwise
   */
  void decide(VBOOL&, const VDOUBLE&, const double, const double, const bool, const bool);

  /*fill the vector of pass flags for a 1D cut as specified in the args
    @param pass vector of pass flags, 1 element per photon, for a specific cut
    @param quantity vector of values of the quantity relating to the specific cut, 1 element per 
    photon
    @param cut value of the cut
    @param min true if the cut should be decided val > cut, false if it should be decided val < cut
    @param abs true if the val should be fabs(val), false otherwise
   */
  void decide1D(VBOOL&, const VDOUBLE&, const double, const bool, const bool);

  /*fill the vector of pass flags for a 1D cut, val > cut passes
    @param pass vector of pass flags, 1 element per photon, for a specific cut
    @param quantity vector of values of the quantity relating to the specific cut, 1 element per 
    photon
    @param cut value of the cut
   */
  void decide1DMin(VBOOL&, const VDOUBLE&, const double);

  /*fill the vector of pass flags for a 1D cut, val < cut passes
    @param pass vector of pass flags, 1 element per photon, for a specific cut
    @param quantity vector of values of the quantity relating to the specific cut, 1 element per 
    photon
    @param cut value of the cut
   */
  void decide1DMax(VBOOL&, const VDOUBLE&, const double);

  /*fill the vector of pass flags for a 1D cut, fabs(val) > cut passes
    @param pass vector of pass flags, 1 element per photon, for a specific cut
    @param quantity vector of values of the quantity relating to the specific cut, 1 element per 
    photon
    @param cut value of the cut
   */
  void decide1DMinAbs(VBOOL&, const VDOUBLE&, const double);

  /*fill the vector of pass flags for a 1D cut, fabs(val) < cut passes
    @param pass vector of pass flags, 1 element per photon, for a specific cut
    @param quantity vector of values of the quantity relating to the specific cut, 1 element per 
    photon
    @param cut value of the cut
   */
  void decide1DMaxAbs(VBOOL&, const VDOUBLE&, const double);

  /*fill the vector of pass flags for a 2D cut, val < cut1*ET + cut0 passes
    @param pass vector of pass flags, 1 element per photon, for a specific cut
    @param quantity vector of values of the quantity relating to the specific cut, 1 element per 
    photon
    @param cut0 y-intercept of the cut
    @param cut1 slope of the cut (x variable is photon ET)
   */
  void decide2DMax(VBOOL&, const VDOUBLE&, const double, const double);

  /*set the appropriate element in the pass vectors (true = photon passed cut, false = photon 
    failed cut)*/
  void decideETMin();
  void decideAbsEtaMax();
  void decideECALIsoMax();
  void decideHCALIsoMax();
  void decideHOverEMax();
  void decideTrackIsoMax();
  void decideSigmaIetaIetaMax();
  void decideAbsSeedTimeMax();
  void decideE2OverE9Max();

};

typedef std::vector<Categorizer> CategorizerCollection;

#endif
