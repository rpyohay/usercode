#include "Skims/SkimTypes/interface/CutScanKey.h"

CutScanKey::CutScanKey() :
  cutName_(""),
  cutVal1_(std::vector<double>()),
  cutVal2_(std::vector<double>()) {}

CutScanKey::CutScanKey(const std::string& cutName, const std::vector<double>& cutVal1, 
		       const std::vector<double>& cutVal2) :
  cutName_(cutName),
  cutVal1_(cutVal1),
  cutVal2_(cutVal2) {}

//LCG dictionary doesn't like copy constructors
/*CutScanKey::CutScanKey(const CutScanKey& copy) :
  cutName_(copy.getCutName()),
  cutVal1_(copy.getCutVal1()),
  cutVal2_(copy.getCutVal2()) {}*/

CutScanKey::~CutScanKey() {}

CutScanKey& CutScanKey::operator=(const CutScanKey& copy)
{
  if (this != &copy) {
    cutName_ = copy.getCutName();
    cutVal1_ = copy.getCutVal1();
    cutVal2_ = copy.getCutVal2();
  }
  return *this;
}

void CutScanKey::setCutName(const std::string& cutName) { cutName_ = cutName; }

void CutScanKey::setCutVal1(const std::vector<double>& cutVal1) { cutVal1_ = cutVal1; }

void CutScanKey::setCutVal2(const std::vector<double>& cutVal2) { cutVal2_ = cutVal2; }

std::string CutScanKey::getCutName() const { return cutName_; }

std::vector<double> CutScanKey::getCutVal1() const { return cutVal1_; }

std::vector<double> CutScanKey::getCutVal2() const { return cutVal2_; }
