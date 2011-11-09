#pragma GCC diagnostic ignored "-Wwrite-strings" //needed to get rid of pesky "deprecated conversion from string constant to char *" compilation error

#ifdef __CINT__

//never even gets here...
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
//#pragma GCC diagnostic ignored "-Wformat"
// #pragma GCC diagnostic warning "-Wwrite-strings"

#include "PhysicsTools/TagAndProbe/interface/RooCMSShape.h"
#include "PhysicsTools/TagAndProbe/interface/RooErfXGaussian.h"
#include "TVirtualFFT.h"

#pragma link C++ class RooCMSShape;
#pragma link C++ class RooErfXGaussian;


#pragma link C++ global gROOT;
#pragma link C++ global gEnv;


#endif
