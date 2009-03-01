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

extern "C" {
  
  SEXP pPCov(SEXP knotsD_r, SEXP coordsKnotsD_r, SEXP n_r, SEXP m_r, 
	      SEXP tauSq_r, SEXP sigmaSq_r, SEXP theta_r, SEXP covModel_r){

    //BLAS and LAPACK vars
    char const *lower = "L";
    char const *upper = "U";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *rside = "R";
    char const *lside = "L";
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    const int incOne = 1;

    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int nm = n*m;
    int mm = m*m;
    int i, info;

    double *ct = (double *) R_alloc(nm, sizeof(double));
    double *C_str = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    
    SEXP C_r;
    PROTECT(C_r = allocMatrix(REALSXP, n, n)); 
    
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
    
    //make the correlation matrix
    for(i = 0; i < mm; i++){
      if(nPramPtr == 1)
	(covModelObj->*cov1ParamPtr)(REAL(theta_r)[0], C_str[i], REAL(knotsD_r)[i]);
      else //i.e., 2 parameter matern
	(covModelObj->*cov2ParamPtr)(REAL(theta_r)[0], REAL(theta_r)[1], C_str[i], REAL(knotsD_r)[i]);
    }
    
    for(i = 0; i < nm; i++){
      if(nPramPtr == 1)
	(covModelObj->*cov1ParamPtr)(REAL(theta_r)[0], ct[i], REAL(coordsKnotsD_r)[i]);
      else //i.e., 2 parameter matern
	(covModelObj->*cov2ParamPtr)(REAL(theta_r)[0], REAL(theta_r)[1], ct[i], REAL(coordsKnotsD_r)[i]);
    }
    
    //scale by sigma^2
    F77_NAME(dscal)(&mm, &REAL(sigmaSq_r)[0], C_str, &incOne);	
    F77_NAME(dscal)(&nm, &REAL(sigmaSq_r)[0], ct, &incOne);
    
    F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("pPCov: Cholesky failed (1)\n");}
    F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("pPCov: Cholesky failed (2)\n");}
    
    F77_NAME(dsymm)(rside, upper, &n, &m, &one, C_str, &m, ct, &n, &zero, tmp_nm, &n);
    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, REAL(C_r), &n);
    
    for(i = 0; i < n; i++) REAL(C_r)[i*n+i] += REAL(tauSq_r)[0];
    
    UNPROTECT(1);

    return(C_r);
  }
  
}
