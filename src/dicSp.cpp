// Andrew O. Finley
// Dept. of Forest Resources
// University of Minnesota
// afinley@stat.umn.edu 
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2004 Andrew O. Finley

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <math.h>
using namespace std;

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "covmodel.h"

extern "C" {
  
  SEXP dicSp(SEXP covmod, SEXP noTausq, SEXP X, SEXP Y, SEXP coords, 
		  SEXP theta, SEXP sigmasq, SEXP tausq, SEXP phi, SEXP nu, SEXP verb)
    
  {
    int i,j, verbose, nProtect= 0;

    /*****************************************************************
       get the args and setup
    ******************************************************************/
    /*********************
       verbose
    **********************/
    if(!isInteger(verb))
      error("c++ error: verbose must be numeric");
    
    verbose = INTEGER(verb)[0];

    /*********************
       covariance model
    **********************/
    if(!isString(covmod))
      error("c++ error: cov.model needs to be a string");
    
    string covModel = CHAR(STRING_ELT(covmod, 0));
    bool onePramPtr = true;

    void (covmodel::*cov1ParamPtr)(double, double *r, double *d, int &length) = NULL; 
    void (covmodel::*cov2ParamPtr)(double, double, double *r, double *d, int &length) = NULL;
    
    if(covModel == "exponential"){
      cov1ParamPtr = &covmodel::exponential;
    }else if(covModel == "spherical"){
      cov1ParamPtr = &covmodel::spherical;
    }else if(covModel == "gaussian"){
      cov1ParamPtr = &covmodel::gaussian;
    }else if(covModel == "matern"){
      cov2ParamPtr = &covmodel::matern;
      onePramPtr = false;
    }else{
      error("c++ error: cov.model is not correctly specified");
    }

    //my covmodel object for calling cov function
    covmodel *covModelObj  = new covmodel;
    
    /*********************
       data and samples
    **********************/
    SEXP xdims;
    int xnrow, xncol, nsamples, incOne = 1;
    bool useTausq = false;
    
    //check the X matrix
    if(!isMatrix(X))
      error("c++ error: 'X' is not a matrix \n");
    
    PROTECT(xdims = getAttrib(X, R_DimSymbol)); nProtect++;      
    xnrow = INTEGER(xdims)[0]; 
    xncol = INTEGER(xdims)[1];
    
    //check the Y matrix
    if(!isMatrix(Y))
      error("c++ error: 'Y' is not a matrix \n");
    
    //check theta
    if(!isMatrix(theta))
      error("c++ error: 'theta' is not a matrix \n");
    
    //just make sure I transposed theta for easy access to samples
    if(INTEGER(getAttrib(theta, R_DimSymbol))[0] != xncol)
      error("c++ error: dims of 'theta' are wrong \n");
    
    //number of samples
    nsamples = INTEGER(getAttrib(theta, R_DimSymbol))[1];

    //check var components
    if(!isMatrix(sigmasq))
      error("c++ error: 'sigmasq' is not a matrix \n");
    
    if(INTEGER(noTausq)[0] == 1)
      useTausq = false;
    else if(INTEGER(noTausq)[0] == 0)
      useTausq = true;
    else
      error("c++ error: something wrong with the no.tausq directive \n");
    
    if(useTausq){
      if(!isMatrix(tausq))
	error("c++ error: 'tausq' is not a matrix \n");
    }
    
    if(!isMatrix(phi))
      error("c++ error: 'phi' is not a matrix \n"); 
    
    if(covModel == "matern"){
      if(!isMatrix(nu))
	error("c++ error: 'nu' is not a matrix \n");
    }


    /*********************
     get coords and mk D
    **********************/
    if(!isMatrix(coords))
      error("c++ error: 'coords' is not a matrix \n");
    
    if(INTEGER(getAttrib(coords, R_DimSymbol))[0] != xnrow)
      error("c++ error: nrows of 'coords' is wrong \n");  
    if(INTEGER(getAttrib(coords, R_DimSymbol))[1] != 2)
      error("c++ error: ncols of 'coords' is wrong \n");
    
    int D_length = xnrow*xnrow;
    double* D = (double *) R_alloc(D_length, sizeof(double));
    
    double dmax = 0;
    for(i = 0; i < xnrow; i++){
      for(j = 0; j < xnrow; j++){
	D[i*xnrow+j] = sqrt(pow(REAL(coords)[i]-REAL(coords)[j],2)+pow(REAL(coords)[xnrow+i]-REAL(coords)[xnrow+j],2));
	if(D[i*xnrow+j] > dmax)
	  dmax = D[i*xnrow+j];
      }
    }
     


    /*********************
         get DIC
    **********************/
    if(verbose){
      Rprintf("\n------------------------------------------------------------------------\n");
      Rprintf("\tCalculating Deviance Information Criterion (DIC)");
      Rprintf("\n------------------------------------------------------------------------\n");
      Rprintf("Percent complete: ");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    double *DD; DD = (double *) R_alloc(nsamples, sizeof(double));
    
    double *R; R = (double *) R_alloc(D_length, sizeof(double));
    double *RInvCurrent; RInvCurrent = (double *) R_alloc(D_length, sizeof(double));
    double *tmpXRow; tmpXRow = (double *) R_alloc(xnrow, sizeof(double));
    double *tmpXRow1; tmpXRow1 = (double *) R_alloc(xnrow, sizeof(double));
    int info = 0;
    const char lower = 'L';
    const char upper = 'U';
    const char ntran = 'N';
    const char ytran = 'Y';
    const char rside = 'R';
    const char lside = 'L';
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    double logDetCov = 0;
    int s = 0, status=0, rtnStatus=0;
    double sigmasqTmp, tausqTmp;
    double DDBar, DDBarOmega, pD, DIC;
    DDBar = DDBarOmega = pD = DIC = 0.0;
    double junk;
    int job = 01;

    for(s = 0; s < nsamples; s++){
      
      //make the correlation matrix
      if(onePramPtr)
	(covModelObj->*cov1ParamPtr)(REAL(phi)[s], R, D, D_length);
      else //i.e., 2 parameter matern
	(covModelObj->*cov2ParamPtr)(REAL(phi)[s], REAL(nu)[s], R, D, D_length);

      //scale correlation matrix with sigmasq
      sigmasqTmp = REAL(sigmasq)[s];
      F77_NAME(dscal)(&D_length, &sigmasqTmp, R, &incOne);
      
      //use tausq
      if(useTausq)
	for(i = 0; i < xnrow; i++) R[i*xnrow+i] = R[i*xnrow+i]+REAL(tausq)[s];

      //chol decom
      F77_NAME(dpofa)(R, &xnrow, &xnrow, &info);
      if(info != 0){error("c++ error: Cholesky failed (1), see DIC.sp documentation\n");}
      
      //get logDet
      logDetCov = 0;
      for(i = 0; i < xnrow; i++)
	logDetCov += log(R[i*xnrow+i]);
      logDetCov = 2*logDetCov;
      
      //finish the invert
      F77_NAME(dpodi)(R, &xnrow, &xnrow, &junk, &job);
      if(info != 0){error("c++ error: Cholesky inverse failed (1), see DIC.sp documentation\n");}
   
      //tmp = y - x*theta
      F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &one, &REAL(X)[0], &xnrow, &REAL(theta)[s*xncol], &incOne, &zero, tmpXRow, &incOne);
      F77_NAME(daxpy)(&xnrow, &negOne, &REAL(Y)[0], &incOne, tmpXRow, &incOne);

      //(-1/2) * tmp` * R^{-1} * tmp
      F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
      DD[s] = logDetCov + F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne); //DD[s] = -2.0*(-0.5*logDetCov - 0.5*dotResult);
  
      //report
      if(verbose){
	if(status == 100){
	  Rprintf("%.0f...", (100.0*s)/nsamples);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
	if(rtnStatus == 800){
	  Rprintf("\n\t");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  rtnStatus = 0;
	}
	rtnStatus++;
	status++;
      }
      
      R_CheckUserInterrupt();
    }
    if(verbose){
      Rprintf("\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    for(i = 0; i < nsamples; i++){
      DDBar += DD[i];
    }
    DDBar = DDBar/nsamples;

    //now for DDBarOmega
    //get the means of the params
    double *thetaMean; thetaMean = (double *) R_alloc(xncol, sizeof(double));
    double sigmasqMean, tausqMean, phiMean, nuMean;
    sigmasqMean = tausqMean = phiMean = nuMean = 0.0;

    //zero
    for(i = 0; i < xncol; i++){thetaMean[i] = 0.0;}

    for(i = 0; i < xncol; i++){
      for(s = 0; s < nsamples; s++){
	thetaMean[i] += REAL(theta)[s*xncol+i];
      }
      thetaMean[i] = thetaMean[i]/nsamples;
    }

    for(s = 0; s < nsamples; s++){
      sigmasqMean += REAL(sigmasq)[s];
      if(useTausq)
	tausqMean += REAL(tausq)[s];
      if(onePramPtr){
	phiMean += REAL(phi)[s];
      }else{ //i.e., 2 parameter matern
	phiMean += REAL(phi)[s];
	nuMean += REAL(nu)[s];
      }
    }
    sigmasqMean = sigmasqMean/nsamples;
    if(useTausq)
      tausqMean = tausqMean/nsamples;
    if(onePramPtr){
      phiMean = phiMean/nsamples;
    }else{ //i.e., 2 parameter matern
      phiMean = phiMean/nsamples;
      nuMean = nuMean/nsamples;
    }
    
    //make the correlation matrix
    if(onePramPtr)
      (covModelObj->*cov1ParamPtr)(phiMean, R, D, D_length);
    else //i.e., 2 parameter matern
      (covModelObj->*cov2ParamPtr)(phiMean, nuMean, R, D, D_length);
    
    //scale correlation matrix with sigmasq
    F77_NAME(dscal)(&D_length, &sigmasqMean, R, &incOne);
    
    //use tausq
    if(useTausq)
      for(i = 0; i < xnrow; i++) R[i*xnrow+i] = R[i*xnrow+i]+tausqMean;
    
    //chol decom
    F77_NAME(dpofa)(R, &xnrow, &xnrow, &info);
    if(info != 0){error("c++ error: Cholesky failed (2), see DIC.sp documentation\n");}
    
    //get logDet
    logDetCov = 0;
    for(i = 0; i < xnrow; i++)
      logDetCov += log(R[i*xnrow+i]);
    logDetCov = 2*logDetCov;
    
    //finish the invert
    F77_NAME(dpodi)(R, &xnrow, &xnrow, &junk, &job);
    if(info != 0){error("c++ error: Cholesky inverse failed (2), see DIC.sp documentation\n");}
   
    //tmp = y - x*theta
    F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &one, &REAL(X)[0], &xnrow, &thetaMean[0], &incOne, &zero, tmpXRow, &incOne);
    F77_NAME(daxpy)(&xnrow, &negOne, &REAL(Y)[0], &incOne, tmpXRow, &incOne);
    
    //(-1/2) * tmp` * R^{-1} * tmp
    F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
    DDBarOmega = logDetCov + F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne);;

    /*********************
     make return object
    **********************/
    SEXP dicResults, dicResultsRowNames, dicResultsColNames, dicResultDimName;
    
    PROTECT(dicResults = allocMatrix(REALSXP, 4, 1)); nProtect++; //for Dbar, DbarOmega, pD, and DIC
    REAL(dicResults)[0] = DDBar;
    REAL(dicResults)[1] = DDBarOmega;
    REAL(dicResults)[2] = DDBar - DDBarOmega;
    REAL(dicResults)[3] = DDBar + DDBar - DDBarOmega;   
    
    PROTECT(dicResultsRowNames = allocVector(STRSXP, 4)); nProtect++;
    SET_VECTOR_ELT(dicResultsRowNames, 0, mkChar("bar.D"));
    SET_VECTOR_ELT(dicResultsRowNames, 1, mkChar("D.bar.Theta"));
    SET_VECTOR_ELT(dicResultsRowNames, 2, mkChar("pD"));
    SET_VECTOR_ELT(dicResultsRowNames, 3, mkChar("DIC"));
    
    PROTECT(dicResultsColNames = allocVector(STRSXP, 1)); nProtect++;
    SET_VECTOR_ELT(dicResultsColNames, 0, mkChar("value"));
    
    PROTECT(dicResultDimName = allocVector(VECSXP, 2)); nProtect++;
    SET_VECTOR_ELT(dicResultDimName, 0, dicResultsRowNames);
    SET_VECTOR_ELT(dicResultDimName, 1, dicResultsColNames);
    setAttrib(dicResults, R_DimNamesSymbol, dicResultDimName);

    /*********************
     clean-up and return
    **********************/

    delete covModelObj;

    UNPROTECT(nProtect);
    
    return(dicResults);
  }
  
}
