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

  SEXP predictSp(SEXP covmod, SEXP noTausq, SEXP X, SEXP Y, SEXP coords, 
		 SEXP theta, SEXP sigmasq, SEXP tausq, SEXP phi, SEXP nu, 
		 SEXP D, SEXP gammaD, SEXP predD, SEXP predX, SEXP predJnt, SEXP verb)
    
  {
    int i,j, verbose, predJoint, nProtect= 0;

    /*****************************************************************
       get the args and setup
    ******************************************************************/
    /*********************
       verbose
    **********************/
    if(!isInteger(verb))
      error("c++ error: verbose must be integer");
    
    verbose = INTEGER(verb)[0];

    /*********************
       joint prediction
    **********************/
    if(!isInteger(predJnt))
      error("c++ error: jointPred must be integer");
    
    predJoint = INTEGER(predJnt)[0];


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
       get predictions
    **********************/
    if(verbose){
      Rprintf("\n------------------------------------------------------------------------\n");
      Rprintf("\t\t\tPredicting points");
      Rprintf("\n------------------------------------------------------------------------\n");
      Rprintf("Percent complete: ");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }
    
    
    //idot check 
    if(INTEGER(getAttrib(gammaD, R_DimSymbol))[0] != INTEGER(getAttrib(predD, R_DimSymbol))[0]
       ||INTEGER(getAttrib(gammaD, R_DimSymbol))[0] != INTEGER(getAttrib(predD, R_DimSymbol))[1])
      error("c++ error: something wrong with the predD or gammaD mtrx\n");

    int npredPts = INTEGER(getAttrib(predD, R_DimSymbol))[0];
    int covPred_length = npredPts*npredPts;
    double *predR =  (double *) R_alloc(covPred_length, sizeof(double));

    double *muPred = (double *) R_alloc(npredPts, sizeof(double));
    double *covPred = (double *) R_alloc(covPred_length, sizeof(double));
    
    int gamma_length = npredPts*xnrow;
    double *gamma = (double *) R_alloc(gamma_length, sizeof(double));
    double *tmpGamma = (double *) R_alloc(gamma_length, sizeof(double));
    double *tmpPredRow = (double *) R_alloc(npredPts, sizeof(double));
    double *tmpPredRow1 = (double *) R_alloc(npredPts, sizeof(double));

    //for the independent prediction
    //this is just a row of the gamma matrix as we actually have t(gamma)
    //which is a bit awk for the ind case but nice and efficient for the 
    //joint case, which will should be used more
    double *indGammaD  =  (double *) R_alloc(xnrow, sizeof(double)); 
    double *indGammaR  =  (double *) R_alloc(xnrow, sizeof(double));
    double *tmpIndGamma  =  (double *) R_alloc(xnrow, sizeof(double));
    double indSD = 0.0;
    double indMu = 0.0;

    
    SEXP pred;
    PROTECT(pred = allocMatrix(REALSXP, npredPts, nsamples)); nProtect++;

    int D_nrow = INTEGER(getAttrib(D, R_DimSymbol))[0];
    int D_length = D_nrow*D_nrow;
    double *R = (double *) R_alloc(D_length, sizeof(double));
    double *tmpXRow = (double *) R_alloc(xnrow, sizeof(double));
    double *tmpXRow1 = (double *) R_alloc(xnrow, sizeof(double));
    int info = 0;
    int incZero = 0;
    const char lower = 'L';
    const char upper = 'U';
    const char ntran = 'N';
    const char ytran = 'T';
    const char rside = 'R';
    const char lside = 'L';
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    double logDetCov = 0;
    int s = 0, status=0, rtnStatus=0;
    double sigmasqTmp, tausqTmp;
    double junk;
    int job = 01;

    GetRNGstate();
    
    //for joint prediction
    if(predJoint){
      
      for(s = 0; s < nsamples; s++){
	
	//make the correlation matrix
	if(onePramPtr)
	  (covModelObj->*cov1ParamPtr)(REAL(phi)[s], R, &REAL(D)[0], D_length);
	else //i.e., 2 parameter matern
	  (covModelObj->*cov2ParamPtr)(REAL(phi)[s], REAL(nu)[s], R, &REAL(D)[0], D_length);
	
	//scale correlation matrix with sigmasq
	F77_NAME(dscal)(&D_length, &REAL(sigmasq)[s], R, &incOne);
	
	//use tausq
	if(useTausq)
	  for(i = 0; i < xnrow; i++) R[i*xnrow+i] = R[i*xnrow+i]+REAL(tausq)[s];
	
	//chol decom
	F77_NAME(dpofa)(R, &xnrow, &xnrow, &info);
	if(info != 0){error("c++ error: Cholesky failed (1), see predict.sp documentation\n");}
	
	//finish the invert
	F77_NAME(dpodi)(R, &xnrow, &xnrow, &junk, &job);
	if(info != 0){error("c++ error: Cholesky inverse failed (2), see predict.sp documentation\n");}
	
	//tmpPredRow = X_o*theta
	F77_NAME(dgemv)(&ntran, &npredPts, &xncol, &one, &REAL(predX)[0], &npredPts, &REAL(theta)[s*xncol], &incOne, &zero, tmpPredRow, &incOne);
	
	//make gamma
	if(onePramPtr)
	  (covModelObj->*cov1ParamPtr)(REAL(phi)[s], gamma, &REAL(gammaD)[0], gamma_length);
	else //i.e., 2 parameter matern
	  (covModelObj->*cov2ParamPtr)(REAL(phi)[s], REAL(nu)[s], gamma, &REAL(gammaD)[0], gamma_length);
	
	//scale gamma by sigma^2
	F77_NAME(dscal)(&gamma_length, &REAL(sigmasq)[s], gamma, &incOne);
	
	//tmpGamma = gamma * R{-1}
	F77_NAME(dsymm)(&rside, &upper, &npredPts, &xnrow, &one, R, &xnrow, gamma, &npredPts, &zero, tmpGamma, &npredPts);      
	
	//tmpXrow = y - x*theta
	F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, &REAL(X)[0], &xnrow, &REAL(theta)[s*xncol], &incOne, &zero, tmpXRow, &incOne);
	F77_NAME(daxpy)(&xnrow, &one, &REAL(Y)[0], &incOne, tmpXRow, &incOne);
	
	//tmpGamma * tmpXrow
	F77_NAME(dgemv)(&ntran, &npredPts, &xnrow, &one, tmpGamma, &npredPts, tmpXRow, &incOne, &zero, muPred, &incOne);
	
	//muPred
	F77_NAME(daxpy)(&npredPts, &one, tmpPredRow, &incOne, muPred, &incOne);    
	
	//covPred
	F77_NAME(dgemm)(&ntran, &ytran, &npredPts, &npredPts, &xnrow, &one, tmpGamma, &npredPts, gamma, &npredPts, &zero, covPred, &npredPts);
	
	//make the correlation matrix for the prediction point submatrix
	if(onePramPtr)
	  (covModelObj->*cov1ParamPtr)(REAL(phi)[s], predR, &REAL(predD)[0], covPred_length);
	else //i.e., 2 parameter matern
	  (covModelObj->*cov2ParamPtr)(REAL(phi)[s], REAL(nu)[s], predR, &REAL(predD)[0], covPred_length);
	
	//scale predR by sigma^2
	F77_NAME(dscal)(&covPred_length, &REAL(sigmasq)[s], predR, &incOne);
	
	//use tausq
	if(useTausq)
	  for(i = 0; i < npredPts; i++) predR[i*npredPts+i] = predR[i*npredPts+i]+REAL(tausq)[s];
	
	//predR-covPred
	for(i = 0; i < covPred_length; i++) covPred[i] = predR[i]-covPred[i];
	
	//need the cholesky for mvn
	F77_NAME(dpofa)(covPred, &npredPts, &npredPts, &info);
	if(info != 0){error("c++ error: Cholesky failed (3), see predict.sp documentation\n");}      
	
	//make the draw and put the resulting sample of preds directly in the output matrix
	mvrnorm(&REAL(pred)[s*npredPts], muPred, covPred, npredPts, true); 
	
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

    }else{//independent predictions

      for(s = 0; s < nsamples; s++){//for each sample
	
	//make the correlation matrix
	if(onePramPtr)
	  (covModelObj->*cov1ParamPtr)(REAL(phi)[s], R, &REAL(D)[0], D_length);
	else //i.e., 2 parameter matern
	  (covModelObj->*cov2ParamPtr)(REAL(phi)[s], REAL(nu)[s], R, &REAL(D)[0], D_length);
	
	//scale correlation matrix with sigmasq
	F77_NAME(dscal)(&D_length, &REAL(sigmasq)[s], R, &incOne);
	
	//use tausq
	if(useTausq)
	  for(i = 0; i < xnrow; i++) R[i*xnrow+i] = R[i*xnrow+i]+REAL(tausq)[s];
	
	//chol decom
	F77_NAME(dpofa)(R, &xnrow, &xnrow, &info);
	if(info != 0){error("c++ error: Cholesky failed (1), see predict.sp documentation\n");}
	
	//finish the invert
	F77_NAME(dpodi)(R, &xnrow, &xnrow, &junk, &job);
	if(info != 0){error("c++ error: Cholesky inverse failed (2), see predict.sp documentation\n");}
	

	//tmpPredRow = X_o*theta now tmpPredRow[j] is the expected mean of the jth point 
	F77_NAME(dgemv)(&ntran, &npredPts, &xncol, &one, &REAL(predX)[0], &npredPts, &REAL(theta)[s*xncol], &incOne, &zero, tmpPredRow, &incOne);
	
	//for each pred point
	for(j = 0; j < npredPts; j++){
	  
	  //this is a bit inefficient but works for now, the problem is that I use gammaD = t(gammaD)
	  //so that the joint prediction is more efficient (i.e., no copy required).
	  //So here indGammaD is just the single vector of distance between the given 
	  //point and each observation
	  F77_NAME(dcopy)(&xnrow, &REAL(gammaD)[j], &npredPts, indGammaD, &incOne);
	  
	  //make the gamma vector
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(REAL(phi)[s], indGammaR, indGammaD, xnrow);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(REAL(phi)[s], REAL(nu)[s], indGammaR, indGammaD, xnrow);
	  
	  //scale gamma vector with sigmasq
	  F77_NAME(dscal)(&xnrow, &REAL(sigmasq)[s], indGammaR, &incOne);
	  
	  //tmpGamma = gamma * R{-1}
	  F77_NAME(dsymm)(&rside, &upper, &incOne, &xnrow, &one, R, &xnrow, indGammaR, &incOne, &zero, tmpIndGamma, &incOne);      
	  
	  //get var
	  indSD = sqrt(REAL(sigmasq)[s] + REAL(tausq)[s] - F77_NAME(ddot)(&xnrow, tmpIndGamma, &incOne, indGammaR, &incOne));
	  
	  //get mu
	  //tmpXrow = y - x*theta
	  F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, &REAL(X)[0], &xnrow, &REAL(theta)[s*xncol], &incOne, &zero, tmpXRow, &incOne);
	  F77_NAME(daxpy)(&xnrow, &one, &REAL(Y)[0], &incOne, tmpXRow, &incOne);
	  
	  indMu = tmpPredRow[j] + F77_NAME(ddot)(&xnrow, tmpIndGamma, &incOne, tmpXRow, &incOne);
	  
	  //make the draw
	  REAL(pred)[s*npredPts+j] = rnorm(indMu, indSD);
 	  
	}
	
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
      
    }

    PutRNGstate();

    if(verbose){
      Rprintf("\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }
    
    delete covModelObj;
    
    UNPROTECT(nProtect);
    
    return(pred);
  }
  
}
