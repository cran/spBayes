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
  
  SEXP mvMPPCovInvDet_wrap(SEXP knotsD_r, SEXP coordsKnotsD_r, SEXP C_r, SEXP C_str_r, SEXP ct_r, 
			   SEXP E_r, SEXP n_r, SEXP m_r, SEXP q_r, SEXP Psi_r, SEXP V_r, SEXP theta_r, 
			   SEXP tmp_mm_r, SEXP tmp_mm1_r, SEXP tmp_mm2_r, SEXP tmp_nmqm_r, 
			   SEXP tmp_nmqm1_r, SEXP tmp_qmqm_r, SEXP tmp_qmqm1_r, SEXP lwork_r, SEXP work_r,
			   SEXP covModel_r, SEXP brute_r){
  
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
  
  REAL(logDet_r)[0] = mvMPPCovInvDet(REAL(knotsD_r), REAL(coordsKnotsD_r), REAL(C_r), REAL(C_str_r), REAL(ct_r), 
					      REAL(E_r), INTEGER(n_r)[0], INTEGER(m_r)[0], INTEGER(q_r)[0], REAL(Psi_r), REAL(V_r), REAL(theta_r), 
					      REAL(tmp_mm_r), REAL(tmp_mm1_r), REAL(tmp_mm2_r), REAL(tmp_nmqm_r), 
					      REAL(tmp_nmqm1_r), REAL(tmp_qmqm_r), REAL(tmp_qmqm1_r), INTEGER(lwork_r)[0], REAL(work_r),
					      covModel, nPramPtr, covModelObj, cov1ParamPtr, cov2ParamPtr, INTEGER(brute_r)[0]);
  
  UNPROTECT(1);
  
  return(logDet_r);
  }
  
}
