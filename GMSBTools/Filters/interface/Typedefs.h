#ifndef GMSBTools_Filters_Typedefs_h
#define GMSBTools_Filters_Typedefs_h

#include <vector>
#include <sstream>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "FWCore/Framework/interface/Event.h"

typedef std::string STRING;
typedef std::vector<STRING> VSTRING;
typedef VSTRING::const_iterator VSTRING_IT;
typedef std::vector<unsigned int> VUINT;
typedef VUINT::const_iterator VUINT_IT;
typedef std::vector<int> VINT;
typedef VINT::const_iterator VINT_IT;
typedef std::vector<double> VDOUBLE;
typedef VDOUBLE::const_iterator VDOUBLE_IT;
typedef std::stringstream STRINGSTREAM;
typedef std::vector<bool> VBOOL;
typedef VBOOL::const_iterator VBOOL_IT;
typedef std::vector<TH1F*> VTH1F;
typedef VTH1F::iterator VTH1F_IT;
typedef VTH1F::const_iterator VTH1F_CONST_IT;
typedef std::vector<float> VFLOAT;
typedef std::vector<TH2F*> VTH2F;
typedef VTH2F::iterator VTH2F_IT;
typedef std::vector<TGraphAsymmErrors*> VTGRAPHASYMMERRORS;
typedef VTGRAPHASYMMERRORS::iterator VTGRAPHASYMMERRORS_IT;
typedef std::pair<edm::RunNumber_t, edm::EventNumber_t> RUNEVTPAIR;
typedef std::pair<RUNEVTPAIR, edm::LuminosityBlockNumber_t> RUNEVTLUMIPAIR;
typedef std::map<RUNEVTPAIR, edm::LuminosityBlockNumber_t> RUNEVTLUMIMAP;
typedef std::vector<TFile> VTFILE;
typedef VTFILE::iterator VTFILE_IT;

#endif
