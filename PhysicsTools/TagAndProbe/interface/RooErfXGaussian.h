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

#ifndef ROO_ERF_X_GAUSSIAN
#define ROO_ERF_X_GAUSSIAN

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "TMath.h"
#include "RooMath.h"

class RooErfXGaussian : public RooAbsPdf {
public:
  RooErfXGaussian() {};
  RooErfXGaussian(const char *name, const char *title,
		  RooAbsReal& _x,
		  RooAbsReal& _alpha,
		  RooAbsReal& _beta,
		  RooAbsReal& _mu,
		  RooAbsReal& _sigma);

  RooErfXGaussian(const RooErfXGaussian& other, const char* name);
  inline virtual TObject* clone(const char* newname) const
  {
    return new RooErfXGaussian(*this,newname);
  }
  inline ~RooErfXGaussian() {}
  Double_t evaluate() const ;
  

  ClassDef(RooErfXGaussian,1);

protected:

  RooRealProxy x ;
  RooRealProxy alpha ;
  RooRealProxy beta ;
  RooRealProxy mu ;
  RooRealProxy sigma ;
  
};
 
#endif
