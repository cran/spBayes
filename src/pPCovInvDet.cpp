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

double pPCovInvDet(double *knotsD, double *coordsKnotsD, double *C, double *C_str, double *ct, 
		   int n, int m, double tauSq, double sigmaSq, double* theta, 
		   double *tmp_mm, double *tmp_nm,
		   string covModel, int nPramPtr, covmodel *covModelObj, 
		   void (covmodel::*cov1ParamPtr)(double, double &, double &),
		   void (covmodel::*cov2ParamPtr)(double, double, double &, double&), int brute){
  
  int nm = n*m;
  int mm = m*m;
  int i, j, info;
  double tauSqInv, negTauSqInv, logDet = 0.0;

  if(tauSq == 0) error("c++ error: tauSq == 0 in pPCovInvDet\n");

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
      
    tauSqInv = 1.0/tauSq;
    negTauSqInv = -1.0*tauSqInv;
    
    //t(ct) 1/tau^2  ct = tmp_mm
    F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &tauSqInv, ct, &n, ct, &n, &zero, tmp_mm, &m);
    
    //[C* + t(ct) 1/tau^2  ct]^{-1} = tmp_mm, and get the det on the way
    F77_NAME(daxpy)(&mm, &one, C_str, &incOne, tmp_mm, &incOne);
    
    F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (1)\n");}
    
    //get log det cov
    logDet = 0;
    for(i = 0; i < m; i++) logDet += 2.0*log(tmp_mm[i*m+i]);
    
    F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (2)\n");}
    
    //-1/tau^2 ct  tmp_mm = tmp_nm
    F77_NAME(dsymm)(rside, upper, &n, &m, &negTauSqInv, tmp_mm, &m, ct, &n, &zero, tmp_nm, &n);
    
    //diag(1/tau^2) + tmp_nm ct 1/tau^2 = C
    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &tauSqInv, tmp_nm, &n, ct, &n, &zero, C, &n);
    for(i = 0; i < n; i++) C[i*n+i] = tauSqInv+C[i*n+i];
    
    //finish getting the log det cov 
    logDet += n*log(tauSq);
    F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (3)\n");}
    for(i = 0; i < m; i++) logDet -= 2.0*log(C_str[i*m+i]);
    
  }else{

    F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (4)\n");}
    F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (5)\n");}
    
    F77_NAME(dsymm)(rside, upper, &n, &m, &one, C_str, &m, ct, &n, &zero, tmp_nm, &n);
    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, C, &n);
    
    for(i = 0; i < n; i++) C[i*n+i] += tauSq;
    
    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (6)\n");}
    
    logDet = 0;
    for(i = 0; i < n; i++) logDet += 2*log(C[i*n+i]);
    
    F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("pPCovInvDet: Cholesky failed (7)\n");}
  }

  return logDet;

}
