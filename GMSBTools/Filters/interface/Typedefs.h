#ifndef GMSBTools_Filters_Typedefs_h
#define GMSBTools_Filters_Typedefs_h

#include <vector>
#include <sstream>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"

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
typedef std::vector<TH1F*>::iterator VTH1F_IT;
typedef std::vector<TH1F*>::const_iterator VTH1F_CONST_IT;
typedef std::vector<float> VFLOAT;
typedef std::vector<float>::const_iterator VFLOAT_IT;
typedef std::vector<TH2F*> VTH2F;
typedef std::vector<TH2F*>::iterator VTH2F_IT;
typedef std::vector<TGraphAsymmErrors*> VTGRAPHASYMMERRORS;
typedef std::vector<TGraphAsymmErrors*>::iterator VTGRAPHASYMMERRORS_IT;
typedef std::pair<Int_t, ULong_t> RUNEVTPAIR;
typedef std::pair<RUNEVTPAIR, Int_t> RUNEVTLUMIPAIR;
typedef std::map<RUNEVTPAIR, Int_t> RUNEVTLUMIMAP;
typedef std::vector<TFile> VTFILE;
typedef std::vector<TFile>::iterator VTFILE_IT;
typedef std::pair<unsigned int, unsigned int> PHOTON_PAIR;
typedef std::map<std::pair<unsigned int, unsigned int>, double> PHOTON_PAIR_MAP;
typedef std::map<std::pair<unsigned int, unsigned int>, double>::const_iterator PHOTON_PAIR_MAP_IT;
typedef std::map<std::pair<unsigned int, unsigned int>, std::pair<float, float> > 
  PHOTON_ET_PAIR_MAP;
typedef std::map<std::pair<unsigned int, unsigned int>, std::pair<float, float> >::const_iterator 
  PHOTON_ET_PAIR_MAP_IT;

#endif
