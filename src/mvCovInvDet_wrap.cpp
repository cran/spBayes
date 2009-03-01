#include <iostream>
#include <string>
using namespace std;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"
#include "covmodel.h"
#include "covInvDet.h"

extern "C" {
  
  SEXP mvCovInvDet_wrap(SEXP coordsD_r, SEXP C_r, SEXP n_r, SEXP m_r, SEXP Psi_r, SEXP V_r, SEXP theta_r, 
			SEXP tmp_mm_r, SEXP tmp_mm1_r, SEXP tmp_mm2_r, SEXP covModel_r){
    
  int nPramPtr = 1;
  string covModel = CHAR(STRING_ELT(covModel_r,0));
  
  void (covmodel::*cov1ParamPtr)(double, double &, double &) = NULL; 
  void (covmodel::*cov2ParamPtr)(double, double, double &, double&) = NULL;
  
  if(covModel == "exponential"){
    cov1ParamPtr = &covmodel::exponential;
  }else if(covModel == "spherical"){
    cov1ParamPtr = &covmodel::spherical;
  }else if(covModel == "gaussian"){
    cov1ParamPtr = &covmodel::gaussian;
  }else if(covModel == "matern"){
    cov2ParamPtr = &covmodel::matern;
    nPramPtr = 2;
  }else{
    error("c++ error: cov.model is not correctly specified");
  }
  
  //my covmodel object for calling cov function
  covmodel *covModelObj = new covmodel;
  
  SEXP logDet_r;
  PROTECT(logDet_r = allocVector(REALSXP, 1));
  
  REAL(logDet_r)[0] = mvCovInvDet(REAL(coordsD_r), REAL(C_r), INTEGER(n_r)[0], INTEGER(m_r)[0], REAL(Psi_r), REAL(V_r), REAL(theta_r), 
				  REAL(tmp_mm_r), REAL(tmp_mm1_r), REAL(tmp_mm2_r),
				  covModel, nPramPtr, covModelObj, cov1ParamPtr, cov2ParamPtr);
  
  UNPROTECT(1);
  
  return(logDet_r);
  }
  
}
