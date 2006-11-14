#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "vprior.h"
#include "priors.h"
#include "covmodel.h"
#include "util.h"



extern "C" {
  
  SEXP spPredict(SEXP args){
    int i,j,k,l,s,nProtect= 0;
    
    /*****************************************
                    get stuff
    *****************************************/

    double *X = REAL(getListElement(args, "X"));
    double *Y = REAL(getListElement(args, "Y"));
    double *obsD = REAL(getListElement(args, "obs.D"));
    int nObs = INTEGER(getListElement(args, "n.obs"))[0];
    double *predX = REAL(getListElement(args, "pred.X"));
    double *predD = REAL(getListElement(args, "pred.D"));
    int nPred = INTEGER(getListElement(args, "n.pred"))[0];
    double *predObsD = REAL(getListElement(args, "predObs.D"));
    double *ASamples  = REAL(getListElement(args, "K"));
    int ACase = INTEGER(getListElement(args, "K.case"))[0]; 
    double *PsiSamples;
    bool noPsi = INTEGER(getListElement(args, "no.Psi"))[0];
    int PsiCase;
    if(!noPsi){
      PsiCase = INTEGER(getListElement(args, "Psi.case"))[0];
      PsiSamples = REAL(getListElement(args, "Psi"));
    }
    double *phiSamples = REAL(getListElement(args, "phi"));
    int phiCase = INTEGER(getListElement(args, "phi.case"))[0];  
    double *nuSamples;
    string covModel = CHAR(STRING_ELT(getListElement(args, "cov.model"), 0));
    if(covModel == "matern")
      nuSamples = REAL(getListElement(args, "nu"));
    double *betaSamples = REAL(getListElement(args, "beta"));
    int nSamples = INTEGER(getListElement(args, "n.samples"))[0];
    int m = INTEGER(getListElement(args, "m"))[0];
    int xnrow = INTEGER(getListElement(args, "xrows"))[0]; 
    int xncol = INTEGER(getListElement(args, "xcols"))[0]; 
    bool verbose = INTEGER(getListElement(args, "verbose"))[0];

    //just check
    if(xnrow != m*nObs)
      error("c++ error: xnrow != m*nObs");

    double *w;
    SEXP w_r;
    bool haveW = INTEGER(getListElement(args, "sp.effect"))[0];
    if(haveW){
      w = REAL(getListElement(args, "w"));
    }else{
      PROTECT(w_r = allocMatrix(REALSXP, xnrow, nSamples)); nProtect++;
      zeros(REAL(w_r), xnrow*nSamples);
      w = REAL(w_r);
    }

    
    /*****************************************
        set-up cov. model function pointer
    *****************************************/
    bool onePramPtr = true;
    
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
      onePramPtr = false;
    }else{
      error("c++ error: cov.model is not correctly specified");
    }
   
    //my covmodel object for calling cov function
    covmodel *covModelObj = new covmodel;


    /*****************************************
        calculate w if needed
    *****************************************/
    bool linpack = false;
    int info;
    const char lower = 'L';
    const char upper = 'U';
    const char ntran = 'N';
    const char ytran = 'T';
    const char rside = 'R';
    const char lside = 'L';
    const double one = 1.0;
    const int incOne = 1;
    const double negOne = -1.0;
    const double zero = 0.0;
    double junk;
    int job = 01;
    double sigmasqTmp;
    double tausqTmp;
    double logDetR = 0;
    double logDetRCurrent = 0;
    double r = 0;
    double rUnif = 0;
    int dnrow = nObs;
    int dlength = dnrow*dnrow;
    int rnrow = dnrow*m;
    int rlength = rnrow*rnrow;
    double *R = (double *) R_alloc(rlength, sizeof(double));
    double *RCurrent = (double *) R_alloc(rlength, sizeof(double));
    double *IKPsi = (double *) R_alloc(rlength, sizeof(double)); zeros(IKPsi, rlength);
    double *wMu = (double *) R_alloc(rnrow, sizeof(double));
    int Alength = m*m;
    int nLowerTriA = (m*m-m)/2+m;
    double *mmblk = (double *) R_alloc(Alength, sizeof(double));
    double *tmp = (double *) R_alloc(Alength, sizeof(double));
    double *A = (double *) R_alloc(Alength, sizeof(double));
    double *Psi = (double *) R_alloc(Alength, sizeof(double));
    double *tmpXRow = (double *) R_alloc(xnrow, sizeof(double));
    double *tmpXRow1 = (double *) R_alloc(xnrow, sizeof(double));
    double *tmpXCol = (double *) R_alloc(xncol, sizeof(double));
    double *tmpXCol1 = (double *) R_alloc(xncol, sizeof(double));
    double *tmpXRowCol = (double *) R_alloc(xnrow*xncol, sizeof(double));
    double *tmpXColCol = (double *) R_alloc(xncol*xncol, sizeof(double));
    int status = 0;

    GetRNGstate();
    if(!haveW){
      if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\t\tRecovering spatial effects\n\t\tfor observed points\n");
	Rprintf("-------------------------------------------------\n");
        #ifdef Win32
        R_FlushConsole();
        #endif
      }      
      
      for(s = 0; s < nSamples; s++){

	if(m == 1){
	  if(!noPsi){
	    
	    //make the correlation matrix
	    for(i = 0; i < dlength; i++){
	      if(onePramPtr)
		(covModelObj->*cov1ParamPtr)(phiSamples[s], R[i], obsD[i]);
	      else //i.e., 2 parameter matern
		(covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], R[i], obsD[i]);
	    }
	    
	    //chol decom
	    if(!linpack){
	      F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed\n");}
	    }else{
	      F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed\n");}
	    }
	    
	    //finish the invert
	    if(!linpack){
	      F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky inverse failed\n");}
	    }else{
	      F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	      if(info != 0){error("c++ error: Cholesky inverse failed\n");}
	    }
	    
	    //scale correlation matrix with sigmasq
	    sigmasqTmp = 1.0/ASamples[s];
	    F77_NAME(dscal)(&dlength, &sigmasqTmp, R, &incOne);
	    
	    for(i = 0; i < rnrow; i++) 
	      R[i*rnrow+i] = R[i*rnrow+i]+1.0/PsiSamples[s];
	    
	    //chol decom
	    if(!linpack){
	      F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed\n");}
	    }else{
	      F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed\n");}
	    }
	    
	    //finish the invert
	    if(!linpack){
	      F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky inverse failed\n");}
	    }else{
	      F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	      if(info != 0){error("c++ error: Cholesky inverse failed\n");}
	    }
	    
	    //make w mu
	    F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, &betaSamples[s*xncol], 
			    &incOne, &zero, tmpXRow, &incOne);
	    
	    F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
	    
	    tausqTmp = 1.0/PsiSamples[s];
	    F77_NAME(dscal)(&xnrow, &tausqTmp, tmpXRow, &incOne);
	    
	    F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, wMu, &incOne);
	    
	    
	    //chol decom for the mvnorm
	    if(!linpack){
	      F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed\n");}
	    }else{
	      F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed\n");}
	    }
	    
	    mvrnorm(&REAL(w_r)[s*xnrow], wMu, R, xnrow, true);
	  }else{
	    F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, &betaSamples[s*xncol], 
			    &incOne, &zero, &REAL(w_r)[s*xnrow], &incOne);
	    
	    F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, &REAL(w_r)[s*xnrow], &incOne);
	  }

	}else{ //m > 1
	  if(!noPsi){
	    
	    //get A sample
	    zeros(A, Alength);
	    if(ACase == 1){
	      for(i = 0; i < m; i++){A[i*m+i] = sqrt(ASamples[s]);}
	    }else if(ACase == 2){
	      for(i = 0; i < m; i++){A[i*m+i] = sqrt(ASamples[s*m+i]);}
	    }else{
	      setLowerChol(A, &ASamples[s*nLowerTriA], m);
	    }
	    
	    //make the correlation matrix
	    zeros(mmblk, Alength);
	    for(i = 0; i < dnrow; i++){
	      for(j = 0; j < dnrow; j++){
		
		for(k = 0; k < m; k++){
		  
		  if(phiCase == 1){//separable
		    if(onePramPtr)
		      (covModelObj->*cov1ParamPtr)(phiSamples[s], mmblk[k*m+k], obsD[j*dnrow+i]);
		    else //i.e., 2 parameter matern
		      (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], mmblk[k*m+k], obsD[j*dnrow+i]);
		  }else{//require NuCase == PhiCase == 2 separable
		    if(onePramPtr)
		      (covModelObj->*cov1ParamPtr)(phiSamples[s*m+k], mmblk[k*m+k], obsD[j*dnrow+i]);
		    else //i.e., 2 parameter matern
		      (covModelObj->*cov2ParamPtr)(phiSamples[s*m+k], nuSamples[s*m+k], mmblk[k*m+k], obsD[j*dnrow+i]);
		  }
		}
		
		F77_NAME(dgemm)(&ntran, &ntran, &m, &m, &m, &one, A, &m, mmblk, &m, &zero, tmp, &m);
		F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, tmp, &m, A, &m, &zero, mmblk, &m);
		
		for(k = 0; k < m; k++){
		  for(l = 0; l < m; l++){
		    R[((j*m+l)*rnrow)+(i*m+k)] = mmblk[l*m+k];
		    mmblk[l*m+k] = 0.0; //zero out
		  }
		}
	      }
	    }
	    
	    //chol decom
	    if(!linpack){
	      F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	    }else{
	      F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	    }
	    
	    //finish the invert
	    if(!linpack){
	      F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	    }else{
	      F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	      if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	    }
	    
	    
	    //get Psi sample
	    zeros(Psi, Alength);
	    if(PsiCase == 1){
	      for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiSamples[s]);}
	    }else if(PsiCase == 2){
	      for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiSamples[s*m+i]);}
	    }else{
	      setLowerChol(Psi, &PsiSamples[s*nLowerTriA], m);
	    }
	    
	    //already chol in lower tri
	    F77_NAME(dpotri)(&lower, &m, Psi, &m, &info);
	    
	    //transpose fill upper tri
	    for(i = 0; i < m; i++)
	      for(j = i; j < m; j++)
		Psi[j*m+i] = Psi[i*m+j];
	    
	    //block diag (aka I_n \otimes Psi^{-1}
	    for(i = 0; i < dnrow; i++){
	      for(k = 0; k < m; k++){
		for(l = 0; l < m; l++){
		  IKPsi[((i*m+l)*rnrow)+(i*m+k)] = Psi[l*m+k];
		}
	      }
	    }
	    
	    F77_NAME(daxpy)(&rlength, &one, IKPsi, &incOne, R, &incOne);
	    
	    
	    //chol decom
	    if(!linpack){
	      F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	    }else{
	      F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	    }
	    
	    //finish the invert
	    if(!linpack){
	      F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	    }else{
	      F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	      if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	    }
	    
	    //make w mu
	    F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow,  &betaSamples[s*xncol], 
			    &incOne, &zero, tmpXRow, &incOne);
	    F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
	    
	    F77_NAME(dsymv)(&upper, &xnrow, &one, IKPsi, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
	    F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow1, &incOne, &zero, wMu, &incOne);
	    
	    //chol decom for the mvnorm
	    if(!linpack){
	      F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	    }else{
	      F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	      if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	    }
	    
	    mvrnorm(&REAL(w_r)[s*xnrow], wMu, R, xnrow, true);
	  }else{
	    F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow,  &betaSamples[s*xncol], 
			    &incOne, &zero, &REAL(w_r)[s*xnrow], &incOne);
	    F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, &REAL(w_r)[s*xnrow], &incOne);
	  }
	    
	}
	
	if(verbose){
	  if(status == 100){
	    Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
            #ifdef Win32
	    R_FlushConsole();
            #endif
	    status = 0;
	  }
	  status++;
	}
	R_CheckUserInterrupt();
      }//end sample loop
      
      if(verbose){
	Rprintf("Sampled: %i of %i, %3.2f%%\n", nSamples, nSamples, 100.0);
        #ifdef Win32
	R_FlushConsole();
        #endif
	status = 0;
	R_CheckUserInterrupt();
      }
    }//fi get w
      
    /*****************************************
       predict w and y
    *****************************************/

    int nSObs = m*nObs;
    int SObsLength = nSObs*nSObs;
    double *obsS = (double *) R_alloc(SObsLength, sizeof(double));

    int nSPred = m*nPred;
    int SPredLength = nSPred*nSPred;
    double *predS = (double *) R_alloc(SPredLength, sizeof(double));

    int SPredObsLength = nSPred*nSObs;
    double *predObsS  = (double *) R_alloc(SPredObsLength, sizeof(double));

    double *tmpSPredSObs = (double *) R_alloc(nSPred*nSObs, sizeof(double));   
    double *wPredMu = (double *) R_alloc(nSPred, sizeof(double));
    double *wPredS = (double *) R_alloc(nSPred*nSPred, sizeof(double));   
    
    SEXP wPred;
    PROTECT(wPred = allocMatrix(REALSXP, nSPred, nSamples)); nProtect++;
    zeros(REAL(wPred), nSPred*nSamples);

    SEXP yPred;
    PROTECT(yPred = allocMatrix(REALSXP, nSPred, nSamples)); nProtect++;
    zeros(REAL(yPred), nSPred*nSamples);
    status=0;
    
    
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tPredicting points\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }      

    for(s = 0; s < nSamples; s++){
      if(m == 1){ 

	//
	//obsS inv
	//
	//make the correlation matrix
	for(i = 0; i < SObsLength; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiSamples[s], obsS[i], obsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], obsS[i], obsD[i]);
	}
	
	//scale correlation matrix with sigmasq
	F77_NAME(dscal)(&SObsLength, &ASamples[s], obsS, &incOne);

	//chol decom
	if(!linpack){
	  F77_NAME(dpotrf)(&upper, &nSObs, obsS, &nSObs, &info);
	  if(info != 0){error("c++ error: Cholesky failed 1\n");}
	}else{
	  F77_NAME(dpofa)(obsS, &nSObs, &nSObs, &info);
	  if(info != 0){error("c++ error: Cholesky failed 1\n");}
	}
	
	//finish the invert
	if(!linpack){
	  F77_NAME(dpotri)(&upper, &nSObs, obsS, &nSObs, &info);
	  if(info != 0){error("c++ error: Cholesky inverse failed\n");}
	}else{
	  F77_NAME(dpodi)(obsS, &nSObs, &nSObs, &junk, &job);
	  if(info != 0){error("c++ error: Cholesky inverse failed\n");}
	}

	//
	//predS
	//
	for(i = 0; i < SPredLength ; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiSamples[s], predS[i], predD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], predS[i], predD[i]);
	}
	
	//scale correlation matrix with sigmasq
	F77_NAME(dscal)(&SPredLength, &ASamples[s], predS, &incOne);	


	//
	//predObsS
	//
	//make the correlation matrix
	for(i = 0; i < SPredObsLength; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiSamples[s], predObsS[i], predObsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], predObsS[i], predObsD[i]);
	}
	
	//scale correlation matrix with sigmasq
	F77_NAME(dscal)(&SPredObsLength, &ASamples[s], predObsS, &incOne);

	//get predObsS*obsS^{-1}
	F77_NAME(dsymm)(&rside, &upper, &nSPred, &nSObs, &one, obsS, &nSObs, predObsS, 
			&nSPred, &zero, tmpSPredSObs, &nSPred);      
	  
	//get mu tmpSPredSObs * w
	F77_NAME(dgemv)(&ntran, &nSPred, &nSObs, &one, tmpSPredSObs, &nSPred, 
			&w[s*xnrow], &incOne, &zero, wPredMu, &incOne);

	//get Sigma presS - predObsS*obsS^{-1}*t(predObsS)
	F77_NAME(dgemm)(&ntran, &ytran, &nSPred, &nSPred, &nSObs, &negOne, tmpSPredSObs, &nSPred, predObsS, &nSPred, 
			&zero, wPredS, &nSPred);
	F77_NAME(daxpy)(&SPredLength, &one, predS, &incOne, wPredS, &incOne);

	//chol decom for the mvnorm
	if(!linpack){
	  F77_NAME(dpotrf)(&upper, &nSPred, wPredS, &nSPred, &info);
	  if(info != 0){error("c++ error: Cholesky failed 2\n");}
	}else{
	  F77_NAME(dpofa)(wPredS, &nSPred, &nSPred, &info);
	  if(info != 0){error("c++ error: Cholesky failed 2\n");}
	}

	//make the draw	
	mvrnorm(&REAL(wPred)[s*nSPred], wPredMu, wPredS, nSPred, true);

	//get prediction XB+w, just use wPredMu again for yPredMu
	F77_NAME(dcopy)(&nSPred, &REAL(wPred)[s*nSPred], &incOne, &REAL(yPred)[s*nSPred], &incOne);

	F77_NAME(dgemv)(&ntran, &nSPred, &xncol, &one, predX, 
			&nSPred, &betaSamples[s*xncol], &incOne, &zero, wPredMu, &incOne);

	F77_NAME(daxpy)(&nSPred, &one, wPredMu, &incOne, &REAL(yPred)[s*nSPred], &incOne);
	

      }else{ //m > 1
	
	//
	//obsS inv
	//
	//get A sample
	zeros(A, Alength);
	if(ACase == 1){
	  for(i = 0; i < m; i++){A[i*m+i] = sqrt(ASamples[s]);}
	}else if(ACase == 2){
	  for(i = 0; i < m; i++){A[i*m+i] = sqrt(ASamples[s*m+i]);}
	}else{
	  setLowerChol(A, &ASamples[s*nLowerTriA], m);
	}
	
	//make the correlation matrix
	zeros(mmblk, Alength);
	for(i = 0; i < nObs; i++){
	  for(j = 0; j < nObs; j++){
	    
	    for(k = 0; k < m; k++){
	      
	      if(phiCase == 1){//separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiSamples[s], mmblk[k*m+k], obsD[j*nObs+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], mmblk[k*m+k], obsD[j*nObs+i]);
	      }else{//require NuCase == PhiCase == 2 separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiSamples[s*m+k], mmblk[k*m+k], obsD[j*nObs+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiSamples[s*m+k], nuSamples[s*m+k], mmblk[k*m+k], obsD[j*nObs+i]);
	      }
	    }
	    
	    F77_NAME(dgemm)(&ntran, &ntran, &m, &m, &m, &one, A, &m, mmblk, &m, &zero, tmp, &m);
	    F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, tmp, &m, A, &m, &zero, mmblk, &m);
	    
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
		obsS[((j*m+l)*nSObs)+(i*m+k)] = mmblk[l*m+k];
		mmblk[l*m+k] = 0.0; //zero out
	      }
	    }
	  }
	}
	

	//chol decom
	if(!linpack){
	  F77_NAME(dpotrf)(&upper, &nSObs, obsS, &nSObs, &info);
	  if(info != 0){error("c++ error: Cholesky failed 1\n");}
	}else{
	  F77_NAME(dpofa)(obsS, &nSObs, &nSObs, &info);
	  if(info != 0){error("c++ error: Cholesky failed 1\n");}
	}
	
	//finish the invert
	if(!linpack){
	  F77_NAME(dpotri)(&upper, &nSObs, obsS, &nSObs, &info);
	  if(info != 0){error("c++ error: Cholesky inverse failed\n");}
	}else{
	  F77_NAME(dpodi)(obsS, &nSObs, &nSObs, &junk, &job);
	  if(info != 0){error("c++ error: Cholesky inverse failed\n");}
	}


	//
	//predS
	//
	//get A sample
	//make the correlation matrix
	zeros(mmblk, Alength);
	for(i = 0; i < nPred; i++){
	  for(j = 0; j < nPred; j++){
	    
	    for(k = 0; k < m; k++){
	      
	      if(phiCase == 1){//separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiSamples[s], mmblk[k*m+k], predD[j*nPred+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], mmblk[k*m+k], predD[j*nPred+i]);
	      }else{//require NuCase == PhiCase == 2 separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiSamples[s*m+k], mmblk[k*m+k], predD[j*nPred+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiSamples[s*m+k], nuSamples[s*m+k], mmblk[k*m+k], predD[j*nPred+i]);
	      }
	    }
	    
	    F77_NAME(dgemm)(&ntran, &ntran, &m, &m, &m, &one, A, &m, mmblk, &m, &zero, tmp, &m);
	    F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, tmp, &m, A, &m, &zero, mmblk, &m);
	    
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
		predS[((j*m+l)*nSPred)+(i*m+k)] = mmblk[l*m+k];
		mmblk[l*m+k] = 0.0; //zero out
	      }
	    }
	  }
	}

	//
	//predObsS
	//
	//make the correlation matrix
	zeros(mmblk, Alength);
	for(i = 0; i < nPred; i++){
	  for(j = 0; j < nObs; j++){
	    
	    for(k = 0; k < m; k++){
	      
	      if(phiCase == 1){//separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiSamples[s], mmblk[k*m+k], predObsD[j*nPred+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], mmblk[k*m+k], predObsD[j*nPred+i]);
	      }else{//require NuCase == PhiCase == 2 separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiSamples[s*m+k], mmblk[k*m+k], predObsD[j*nPred+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiSamples[s*m+k], nuSamples[s*m+k], mmblk[k*m+k], predObsD[j*nPred+i]);
	      }
	    }
	    
	    F77_NAME(dgemm)(&ntran, &ntran, &m, &m, &m, &one, A, &m, mmblk, &m, &zero, tmp, &m);
	    F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, tmp, &m, A, &m, &zero, mmblk, &m);
	    
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
		predObsS[((j*m+l)*nSPred)+(i*m+k)] = mmblk[l*m+k];
		mmblk[l*m+k] = 0.0; //zero out
	      }
	    }
	  }
	}

	//get predObsS*obsS^{-1}
	F77_NAME(dsymm)(&rside, &upper, &nSPred, &nSObs, &one, obsS, &nSObs, predObsS, 
			&nSPred, &zero, tmpSPredSObs, &nSPred);      
	  
	//get mu tmpSPredSObs * w
	F77_NAME(dgemv)(&ntran, &nSPred, &nSObs, &one, tmpSPredSObs, &nSPred, 
			&w[s*xnrow], &incOne, &zero, wPredMu, &incOne);

	//get Sigma presS - predObsS*obsS^{-1}*t(predObsS)
	F77_NAME(dgemm)(&ntran, &ytran, &nSPred, &nSPred, &nSObs, &negOne, tmpSPredSObs, &nSPred, predObsS, &nSPred, 
			&zero, wPredS, &nSPred);
	F77_NAME(daxpy)(&SPredLength, &one, predS, &incOne, wPredS, &incOne);

	//chol decom for the mvnorm
	if(!linpack){
	  F77_NAME(dpotrf)(&upper, &nSPred, wPredS, &nSPred, &info);
	  if(info != 0){error("c++ error: Cholesky failed 2\n");}
	}else{
	  F77_NAME(dpofa)(wPredS, &nSPred, &nSPred, &info);
	  if(info != 0){error("c++ error: Cholesky failed 2\n");}
	}

	//make the draw	
	mvrnorm(&REAL(wPred)[s*nSPred], wPredMu, wPredS, nSPred, true);

	//get prediction XB+w, just use wPredMu again for yPredMu
	F77_NAME(dcopy)(&nSPred, &REAL(wPred)[s*nSPred], &incOne, &REAL(yPred)[s*nSPred], &incOne);

	F77_NAME(dgemv)(&ntran, &nSPred, &xncol, &one, predX, 
			&nSPred, &betaSamples[s*xncol], &incOne, &zero, wPredMu, &incOne);

	F77_NAME(daxpy)(&nSPred, &one, wPredMu, &incOne, &REAL(yPred)[s*nSPred], &incOne);

      }

      if(verbose){
	if(status == 100){
	  Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
	status++;
      }
      R_CheckUserInterrupt();

    }//end sample loop

    if(verbose){
      Rprintf("Sampled: %i of %i, %3.2f%%\n", nSamples, nSamples, 100.0);
      #ifdef Win32
      R_FlushConsole();
      #endif
      status = 0;
    }

    PutRNGstate();    
    /*********************
         return
    **********************/
    
    //make the result list object
    int nResultListObjs = 2;
    if(!haveW)
      nResultListObjs = 3;
         
    SEXP result, resultNames;
    
    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //set result list elements
    //pred w predicted
    SET_VECTOR_ELT(result, 0, wPred);
    SET_VECTOR_ELT(resultNames, 0, mkChar("pred.sp.effects")); 

    //pred y predicted
    SET_VECTOR_ELT(result, 1, yPred);
    SET_VECTOR_ELT(resultNames, 1, mkChar("pred.y")); 

    if(!haveW){    
      //w observed
      SET_VECTOR_ELT(result, 2, w_r);
      SET_VECTOR_ELT(resultNames, 2, mkChar("sp.effects"));    
    }

    namesgets(result, resultNames);

   //unprotect
   UNPROTECT(nProtect);
   
   return(result);

  }
}

