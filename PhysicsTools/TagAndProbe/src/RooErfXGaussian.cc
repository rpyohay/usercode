/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 *
 *
 * Authors:
 *   Rachel Yohay, UVa - rpy3y@virginia.edu
 *
 * Description:
 *   Defines a probability density function which has Gaussian form around the Z peak  
 *   but turns over (i.e., error function) at low mass due to threshold 
 *   effect. We use this to model the background shape in Z->ll invariant 
 *   mass.
 * History:
 *   
 *
 *****************************************************************************/

#include "PhysicsTools/TagAndProbe/interface/RooErfXGaussian.h"

ClassImp(RooErfXGaussian) 

 RooErfXGaussian::RooErfXGaussian(const char *name, const char *title, 
				  RooAbsReal& _x,
				  RooAbsReal& _alpha,
				  RooAbsReal& _beta,
				  RooAbsReal& _mu,
				  RooAbsReal& _sigma) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   alpha("alpha","alpha",this,_alpha),
   beta("beta","beta",this,_beta),
   mu("mu","mu",this,_mu),
   sigma("sigma","sigma",this,_sigma)
 { } 


 RooErfXGaussian::RooErfXGaussian(const RooErfXGaussian& other, const char* name):
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   alpha("alpha",this,other.alpha),
   beta("beta",this,other.beta),
   mu("mu",this,other.mu),
   sigma("sigma",this,other.sigma)
 { } 



 Double_t RooErfXGaussian::evaluate() const 
 { 
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 

  Double_t rval=0;
  //Double_t erf = TMath::Erfc((alpha - x) * beta);
  Double_t arg = (x - alpha) * beta;
  Double_t erf = RooMath::erfc(arg);
  Double_t t = (x-mu)/sigma;
  if (arg > 0) rval= erf*exp(-0.5*t*t);
  else rval = 0;
  return rval;
 } 
