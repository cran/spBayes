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

double mPPCovInvDet(double *knotsD, double *coordsKnotsD, double *C, double *C_str, double *ct, 
		    double *E, double *Einv, int n, int m, double tauSq, double sigmaSq, double* theta, 
		    double *tmp_mm, double *tmp_nm, double *tmp_nm1, double *tmp_nn,
		    string covModel, int nPramPtr, covmodel *covModelObj, 
		    void (covmodel::*cov1ParamPtr)(double, double &, double &),
		    void (covmodel::*cov2ParamPtr)(double, double, double &, double&), int brute){
  
  //E = (double *) R_alloc(n, sizeof(double));
  //Einv = (double *) R_alloc(n, sizeof(double));

  int nm = n*m;
  int mm = m*m;
  int nn = n*n;
  int i, j, info;
  double tauSqInv, negTauSqInv, logDet = 0.0, logDetC_str = 0.0, logDetE = 0.0;

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
  
  //make the correlation matrix
  for(i = 0; i < mm; i++){
    if(nPramPtr == 1)
      (covModelObj->*cov1ParamPtr)(theta[0], C_str[i], knotsD[i]);
    else //i.e., 2 parameter matern
      (covModelObj->*cov2ParamPtr)(theta[0], theta[1], C_str[i], knotsD[i]);
  }
  
  for(i = 0; i < nm; i++){
    if(nPramPtr == 1)
      (covModelObj->*cov1ParamPtr)(theta[0], ct[i], coordsKnotsD[i]);
    else //i.e., 2 parameter matern
      (covModelObj->*cov2ParamPtr)(theta[0], theta[1], ct[i], coordsKnotsD[i]);
  }
  
  //scale by sigma^2
  F77_NAME(dscal)(&mm, &sigmaSq, C_str, &incOne);	
  F77_NAME(dscal)(&nm, &sigmaSq, ct, &incOne);
  
  if(!brute){

  /******************************
     Sherman-Woodbury-Morrison
  *******************************/
  //make ct C
  F77_NAME(dcopy)(&mm, C_str, &incOne, tmp_mm, &incOne);
  F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in mPPCovInvDet (1)\n");}
  
  logDetC_str = 0;
  for(i = 0; i < m; i++) 
    logDetC_str += 2.0*log(tmp_mm[i*m+i]);
  
  F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in mPPCovInvDet (2)\n");}
  
  F77_NAME(dsymm)(rside, upper, &n, &m, &one, tmp_mm, &m, ct, &n, &zero, tmp_nm, &n);
  F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, tmp_nn, &n);
  
  logDetE = 0;
  for(i = 0; i < n; i++){ 
    E[i] = tauSq+sigmaSq-tmp_nn[i*n+i];
    Einv[i] = 1.0/E[i];
    logDetE += 2.0*log(sqrt(E[i]));
  }
  
  //make Einv Ct
  //F77_NAME(dsymm)(lside, upper, &n, &m, &one, Einv, &n, ct, &n, &zero, tmp_nm, &n);
  diagmm(n, m, Einv, ct, tmp_nm);
  
  //make C* + t(Ct) E.inv Ct
  F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &one, ct, &n, tmp_nm, &n, &zero, tmp_mm, &m);
  F77_NAME(daxpy)(&mm, &one, C_str, &incOne, tmp_mm, &incOne);
  
  //get log(|tmp_mm|) then tmp_mm^{-1}
  logDet = 0.0;
  F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in mPPCovInvDet (3)\n");}
  for(i = 0; i < m; i++) logDet += 2.0*log(tmp_mm[i*m+i]);
  F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in mPPCovInvDet (4)\n");}
  
  logDet = logDetE+logDet-logDetC_str;
  
  //C = Einv - Einv ct (C* + t(ct) Einv ct)^{-1} t(ct) Einv
  F77_NAME(dsymm)(rside, upper, &n, &m, &one, tmp_mm, &m, tmp_nm, &n, &zero, tmp_nm1, &n);
  F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm1, &n, tmp_nm, &n, &zero, C, &n);
  
  F77_NAME(dscal)(&nn, &negOne, C, &incOne);
  for(i = 0; i < n; i++) 
    C[i*n+i] = Einv[i]+C[i*n+i];
  }else{
    F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (5)\n");}
    F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (6)\n");}
    
    F77_NAME(dsymm)(rside, upper, &n, &m, &one, C_str, &m, ct, &n, &zero, tmp_nm, &n);
    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, C, &n);
    
    for(i = 0; i < n; i++) C[i*n+i] = tauSq+sigmaSq;
    
    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (7)\n");}
    
    logDet = 0;
    for(i = 0; i < n; i++) logDet += 2*log(C[i*n+i]);
    
    F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (8)\n");}
  }

  return logDet;

}
