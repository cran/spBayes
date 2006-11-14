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
  
  SEXP dic(SEXP args){
    int i,j,k,l,s,nProtect= 0;
    
    /*****************************************
                    get stuff
    *****************************************/
    bool DICMarg = INTEGER(getListElement(args, "DIC.marg"))[0];
    bool DICUnmarg = INTEGER(getListElement(args, "DIC.unmarg"))[0];
    
    double *X = REAL(getListElement(args, "X"));
    double *Y = REAL(getListElement(args, "Y"));
    double *D = REAL(getListElement(args, "D"));
    int n = INTEGER(getListElement(args, "n"))[0];
    double *ASamples  = REAL(getListElement(args, "K"));
    double *AMeans  = REAL(getListElement(args, "K.means"));
    int ACase = INTEGER(getListElement(args, "K.case"))[0]; 
    double *PsiSamples;
    bool noPsi = INTEGER(getListElement(args, "no.Psi"))[0];
    double *PsiMeans;
    int PsiCase;
    if(!noPsi){
      PsiCase = INTEGER(getListElement(args, "Psi.case"))[0];
      PsiSamples = REAL(getListElement(args, "Psi"));
      PsiMeans = REAL(getListElement(args, "Psi.means"));
    }
    double *phiSamples = REAL(getListElement(args, "phi"));
    double *phiMeans = REAL(getListElement(args, "phi.means"));
    int phiCase = INTEGER(getListElement(args, "phi.case"))[0];  
    double *nuSamples;
    double *nuMeans;
    string covModel = CHAR(STRING_ELT(getListElement(args, "cov.model"), 0));
    if(covModel == "matern"){
      nuSamples = REAL(getListElement(args, "nu"));
      nuMeans = REAL(getListElement(args, "nu.means"));
    }
    double *betaSamples = REAL(getListElement(args, "beta"));
    double *betaMeans = REAL(getListElement(args, "beta.means"));
    int nSamples = INTEGER(getListElement(args, "n.samples"))[0];
    int m = INTEGER(getListElement(args, "m"))[0];
    int xnrow = INTEGER(getListElement(args, "xrows"))[0]; 
    int xncol = INTEGER(getListElement(args, "xcols"))[0]; 
    bool verbose = INTEGER(getListElement(args, "verbose"))[0];
    
    //just check
    if(xnrow != m*n)
      error("c++ error: xnrow != m*n");
    
    double *w;
    SEXP w_r;
    bool haveW = INTEGER(getListElement(args, "sp.effect"))[0];
    
    if(DICUnmarg){
      if(haveW){
	w = REAL(getListElement(args, "w"));
      }else{
	PROTECT(w_r = allocMatrix(REALSXP, xnrow, nSamples)); nProtect++;
	zeros(REAL(w_r), xnrow*nSamples);
	w = REAL(w_r);
      }
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
    double logDetIKPsi = 0;
    double logDetR = 0;
    double logDetRCurrent = 0;
    double r = 0;
    double rUnif = 0;
    int dnrow = n;
    int dlength = dnrow*dnrow;
    int rnrow = dnrow*m;
    int rlength = rnrow*rnrow;
    double *R = (double *) R_alloc(rlength, sizeof(double));
    double *IKPsi = (double *) R_alloc(rlength, sizeof(double)); zeros(IKPsi, rlength);//for w recovery
    double *wMu = (double *) R_alloc(rnrow, sizeof(double));//for w recovery
    double *wMeans;//for DIC
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
    double *In = (double *) R_alloc(dlength, sizeof(double)); identity(In, dnrow);

    
    if(!haveW  && DICUnmarg){
      GetRNGstate();
      if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\t\tRecovering spatial effects\n");
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
		(covModelObj->*cov1ParamPtr)(phiSamples[s], R[i], D[i]);
	      else //i.e., 2 parameter matern
		(covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], R[i], D[i]);
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
		      (covModelObj->*cov1ParamPtr)(phiSamples[s], mmblk[k*m+k], D[j*dnrow+i]);
		    else //i.e., 2 parameter matern
		      (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], mmblk[k*m+k], D[j*dnrow+i]);
		  }else{//require NuCase == PhiCase == 2 separable
		    if(onePramPtr)
		      (covModelObj->*cov1ParamPtr)(phiSamples[s*m+k], mmblk[k*m+k], D[j*dnrow+i]);
		    else //i.e., 2 parameter matern
		      (covModelObj->*cov2ParamPtr)(phiSamples[s*m+k], nuSamples[s*m+k], mmblk[k*m+k], D[j*dnrow+i]);
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

      PutRNGstate(); 
      
      if(verbose){
	Rprintf("Sampled: %i of %i, %3.2f%%\n", nSamples, nSamples, 100.0);
        #ifdef Win32
	R_FlushConsole();
        #endif
	status = 0;
	R_CheckUserInterrupt();
      }
      

    }//fi get w

    
    if(DICUnmarg){
      //get w row means
      wMeans = (double *) R_alloc(rnrow, sizeof(double));
      zeros(wMeans,rnrow);
      
      for(s = 0; s < nSamples; s++)
	for(j = 0; j < rnrow; j++)
	  wMeans[j] += w[s*xnrow+j];
      
      for(j = 0; j < rnrow; j++)
	wMeans[j] = wMeans[j]/nSamples;
    }
    
    /*********************
            DICs
    **********************/
    SEXP DICMargResult;
    if(DICMarg){
      PROTECT(DICMargResult = allocMatrix(REALSXP, 4, 1)); nProtect++; //for Dbar, DbarOmega, pD, and DIC
    }
    
    SEXP DICUnmargResult;
    if(DICUnmarg){  
      PROTECT(DICUnmargResult = allocMatrix(REALSXP, 4, 1)); nProtect++; //for Dbar, DbarOmega, pD, and DIC
    }
    
    double *DMarg;
    double *DUnmarg;
    if(DICMarg){DMarg = (double *) R_alloc(nSamples, sizeof(double));}
    if(DICUnmarg){DUnmarg = (double *) R_alloc(nSamples, sizeof(double));}
  
    double DBarMarg, DBarMargOmega, pDMarg, DBarUnmarg, DBarUnmargOmega, pDUnmarg;
    DBarMarg=DBarMargOmega=pDMarg=DBarUnmarg=DBarUnmargOmega=pDUnmarg = 0.0;

    /*********************
           DICMarg DBar
    **********************/
  
    if(DICMarg){
      if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\t\tCalculating marginalized DIC\n");
	Rprintf("-------------------------------------------------\n");
        #ifdef Win32
	R_FlushConsole();
        #endif
      }
    
      for(s = 0; s < nSamples; s++){
	
	//form covariance matrix
	if(m == 1){
	  
	  //make the correlation matrix
	  for(i = 0; i < dlength; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phiSamples[s], R[i], D[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], R[i], D[i]);
	  }
	  
	  //scale correlation matrix with sigmasq
	  sigmasqTmp = ASamples[s];
	  F77_NAME(dscal)(&dlength, &sigmasqTmp, R, &incOne);
	  
	  //use Psi
	  if(!noPsi){
	    for(i = 0; i < rnrow; i++) 
	      R[i*rnrow+i] = R[i*rnrow+i]+PsiSamples[s];
	  }	 
	  
	}else{ //m > 1
	  
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
		    (covModelObj->*cov1ParamPtr)(phiSamples[s], mmblk[k*m+k], D[j*dnrow+i]);
		  else //i.e., 2 parameter matern
		    (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], mmblk[k*m+k], D[j*dnrow+i]);
		}else{//require NuCase == PhiCase == 2 separable
		  if(onePramPtr)
		    (covModelObj->*cov1ParamPtr)(phiSamples[s*m+k], mmblk[k*m+k], D[j*dnrow+i]);
		  else //i.e., 2 parameter matern
		    (covModelObj->*cov2ParamPtr)(phiSamples[s*m+k], nuSamples[s*m+k], mmblk[k*m+k], D[j*dnrow+i]);
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
	  
	  //use Psi
	  if(!noPsi){
	    
	    //get Psi sample
	    zeros(Psi, Alength);
	    if(PsiCase == 1){
	      for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiSamples[s]);}
	    }else if(PsiCase == 2){
	      for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiSamples[s*m+i]);}
	    }else{
	      setLowerChol(Psi, &PsiSamples[s*nLowerTriA], m);
	    }
	    
	    F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, Psi, &m, Psi, &m, &zero, tmp, &m);	
	    
	    for(i = 0; i < dnrow; i++){
	      for(k = 0; k < m; k++){
		for(l = 0; l < m; l++){
		  R[((i*m+l)*rnrow)+(i*m+k)] = R[((i*m+l)*rnrow)+(i*m+k)]+tmp[l*m+k];
		}
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
	
	//get logDet
	logDetR = 0;
	for(i = 0; i < rnrow; i++)
	  logDetR += log(R[i*rnrow+i]);
	logDetR = 2*logDetR;
	
	//finish the invert
	if(!linpack){
	  F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	  if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	}else{
	  F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	  if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	}
	
	for(i = 0; i < rlength; i++){
	  if(isnan(R[i]))
	    error("c++ error: Lapack dpotri routine produced an undiscovered NAN in the R^{-1} computation. This is a rogue bug we have yet to track down.  Rerun with linpack=TRUE in the run.control list.  If this error persists when linpack=TRUE please report to package maintainer.\n");
	}
	
	F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, &betaSamples[s*xncol], 
			&incOne, &zero, tmpXRow, &incOne);
	F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
	
	//(-1/2) * tmp` * R^{-1} * tmp
	F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
	DMarg[s] = logDetR+F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne); //i.e., DMarg = -2.0*(-0.5*logDetCov - 0.5*dotResult);

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
      
      for(i = 0; i < nSamples; i++)
	DBarMarg += DMarg[i];
      
      DBarMarg = DBarMarg/nSamples;
      
      //
      //Now DBarMargOmega
      //
      //form covariance matrix
      if(m == 1){
	
	//make the correlation matrix
	for(i = 0; i < dlength; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiMeans[0], R[i], D[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiMeans[0], nuMeans[0], R[i], D[i]);
	}
	
	//scale correlation matrix with sigmasq
	sigmasqTmp = AMeans[0];
	F77_NAME(dscal)(&dlength, &sigmasqTmp, R, &incOne);
	
	//use Psi
	if(!noPsi){
	  for(i = 0; i < rnrow; i++) 
	    R[i*rnrow+i] = R[i*rnrow+i]+PsiMeans[0];
	}	 
	
      }else{ //m > 1
	
	//get A sample
	zeros(A, Alength);
	if(ACase == 1){
	  for(i = 0; i < m; i++){A[i*m+i] = sqrt(AMeans[0]);}
	}else if(ACase == 2){
	  for(i = 0; i < m; i++){A[i*m+i] = sqrt(AMeans[i]);}
	}else{
	  setLowerChol(A, &AMeans[0], m);
	}
	
	//make the correlation matrix
	zeros(mmblk, Alength);
	for(i = 0; i < dnrow; i++){
	  for(j = 0; j < dnrow; j++){
	    
	    for(k = 0; k < m; k++){
	      
	      if(phiCase == 1){//separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiMeans[0], mmblk[k*m+k], D[j*dnrow+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiMeans[0], nuMeans[0], mmblk[k*m+k], D[j*dnrow+i]);
	      }else{//require NuCase == PhiCase == 2 separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiMeans[k], mmblk[k*m+k], D[j*dnrow+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiMeans[k], nuMeans[k], mmblk[k*m+k], D[j*dnrow+i]);
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
	
	//use Psi
	if(!noPsi){
	  
	  //get Psi sample
	  zeros(Psi, Alength);
	  if(PsiCase == 1){
	    for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiMeans[0]);}
	  }else if(PsiCase == 2){
	    for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiMeans[i]);}
	  }else{
	    setLowerChol(Psi, &PsiMeans[0], m);
	  }
	  
	  F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, Psi, &m, Psi, &m, &zero, tmp, &m);	
	  
	  for(i = 0; i < dnrow; i++){
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
		R[((i*m+l)*rnrow)+(i*m+k)] = R[((i*m+l)*rnrow)+(i*m+k)]+tmp[l*m+k];
	      }
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
      
      //get logDet
      logDetR = 0;
      for(i = 0; i < rnrow; i++)
	logDetR += log(R[i*rnrow+i]);
      logDetR = 2*logDetR;
      
      //finish the invert
      if(!linpack){
	F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
      }else{
	F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
      }
      
      for(i = 0; i < rlength; i++){
	if(isnan(R[i]))
	  error("c++ error: Lapack dpotri routine produced an undiscovered NAN in the R^{-1} computation. This is a rogue bug we have yet to track down.  Rerun with linpack=TRUE in the run.control list.  If this error persists when linpack=TRUE please report to package maintainer.\n");
      }
      
      F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, &betaMeans[0], 
		      &incOne, &zero, tmpXRow, &incOne);
      F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
      
      //(-1/2) * tmp` * R^{-1} * tmp
      F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
      DBarMargOmega = logDetR+F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne); //i.e., DMarg = -2.0*(-0.5*logDetCov - 0.5*dotResult);
      
      if(verbose){
	Rprintf("Sampled: %i of %i, %3.2f%%\n", nSamples, nSamples, 100.0);
        #ifdef Win32
	R_FlushConsole();
        #endif	
	R_CheckUserInterrupt();
      }
      
    }
  


    /*********************
        DICUnmarg DBar
    **********************/
    
    if(DICUnmarg){
      if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\t\tCalculating unmarginalized DIC\n");
	Rprintf("-------------------------------------------------\n");
        #ifdef Win32
	R_FlushConsole();
        #endif
      }
      
      //if no Psi then Psi = I and logDetIKPsi = 0.
      identity(IKPsi, rnrow);
      logDetIKPsi = 0.0;
      
      for(s = 0; s < nSamples; s++){
	
	//form covariance matrix
	if(m == 1){
	  
	  //make the correlation matrix
	  for(i = 0; i < dlength; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phiSamples[s], R[i], D[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], R[i], D[i]);
	  }
	  
 	  //scale correlation matrix with sigmasq
	  sigmasqTmp = ASamples[s];
	  F77_NAME(dscal)(&dlength, &sigmasqTmp, R, &incOne);
	  
	  //use Psi
	  if(!noPsi){
	    logDetIKPsi = dnrow*log(PsiSamples[s]);
	    
	    for(i = 0; i < dnrow; i++)
	      IKPsi[i*dnrow+i] = 1.0/PsiSamples[s];
	  }//else IKPsi = I and logDetIKPsi = 0.0 set above
	  
	  
	}else{ //m > 1
	  
	  
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
		    (covModelObj->*cov1ParamPtr)(phiSamples[s], mmblk[k*m+k], D[j*dnrow+i]);
		  else //i.e., 2 parameter matern
		    (covModelObj->*cov2ParamPtr)(phiSamples[s], nuSamples[s], mmblk[k*m+k], D[j*dnrow+i]);
		}else{//require NuCase == PhiCase == 2 separable
		  if(onePramPtr)
		    (covModelObj->*cov1ParamPtr)(phiSamples[s*m+k], mmblk[k*m+k], D[j*dnrow+i]);
		  else //i.e., 2 parameter matern
		    (covModelObj->*cov2ParamPtr)(phiSamples[s*m+k], nuSamples[s*m+k], mmblk[k*m+k], D[j*dnrow+i]);
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
	  
	  
	  //use Psi
	  if(!noPsi){
	    //get Psi sample
	    zeros(Psi, Alength);
	    if(PsiCase == 1){
	      for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiSamples[s]);}
	    }else if(PsiCase == 2){
	      for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiSamples[s*m+i]);}
	    }else{
	      setLowerChol(Psi, &PsiSamples[s*nLowerTriA], m);
	    }
	    
	    //get det
	    logDetIKPsi = 0.0;
	    for(i = 0; i < m; i++)
	      logDetIKPsi += log(Psi[i*m+i]);
	    logDetIKPsi = dnrow*2.0*logDetIKPsi;
	    
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
	    
	  }//else IKPsi = I and logDetIKPsi = 0.0 set above
	}
	
	//chol decom
	if(!linpack){
	  F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	  if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	}else{
	  F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	  if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	}
	
	//get logDet
	logDetR = 0;
	for(i = 0; i < rnrow; i++)
	  logDetR += log(R[i*rnrow+i]);
	logDetR = 2*logDetR;
	
	
	//finish the invert
	if(!linpack){
	  F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	  if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	}else{
	  F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	  if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	}
	
	
	F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, &betaSamples[s*xncol], 
			&incOne, &zero, tmpXRow, &incOne);
	F77_NAME(daxpy)(&xnrow, &negOne, &w[s*xnrow], &incOne, tmpXRow, &incOne);
	F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
	
	//(-1/2) * tmp` * Psi^{-1} * tmp
	F77_NAME(dsymv)(&upper, &xnrow, &one, IKPsi, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
	DUnmarg[s] = logDetIKPsi+F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne);


	//get w likelihood
	//(-1/2) * w` * R`^{-1} * w
	F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, &w[s*xnrow], &incOne, &zero, tmpXRow1, &incOne);
	DUnmarg[s] += logDetR+F77_NAME(ddot)(&xnrow, &w[s*xnrow], &incOne, tmpXRow1, &incOne);

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
      
      for(i = 0; i < nSamples; i++){
	DBarUnmarg += DUnmarg[i];
      }
      
      DBarUnmarg = DBarUnmarg/nSamples;
      
    
      //
      //Now DBarUnmargOmega
      //
      //form covariance matrix
      if(m == 1){
	
	//make the correlation matrix
	for(i = 0; i < dlength; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiMeans[0], R[i], D[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiMeans[0], nuMeans[0], R[i], D[i]);
	}
	
	//scale correlation matrix with sigmasq
	sigmasqTmp = AMeans[0];
	F77_NAME(dscal)(&dlength, &sigmasqTmp, R, &incOne);
	
	//use Psi
	if(!noPsi){
	  logDetIKPsi = dnrow*log(PsiMeans[0]);
	  
	  for(i = 0; i < dnrow; i++)
	    IKPsi[i*dnrow+i] = 1.0/PsiMeans[0];
	}//else IKPsi = I and logDetIKPsi = 0.0 set above
	
	
      }else{ //m > 1
	
	//get A sample
	zeros(A, Alength);
	if(ACase == 1){
	  for(i = 0; i < m; i++){A[i*m+i] = sqrt(AMeans[0]);}
	}else if(ACase == 2){
	  for(i = 0; i < m; i++){A[i*m+i] = sqrt(AMeans[i]);}
	}else{
	  setLowerChol(A, &AMeans[0], m);
	}
	
	//make the correlation matrix
	zeros(mmblk, Alength);
	for(i = 0; i < dnrow; i++){
	  for(j = 0; j < dnrow; j++){
	    
	    for(k = 0; k < m; k++){
	      
	      if(phiCase == 1){//separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiMeans[0], mmblk[k*m+k], D[j*dnrow+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiMeans[0], nuMeans[0], mmblk[k*m+k], D[j*dnrow+i]);
	      }else{//require NuCase == PhiCase == 2 separable
		if(onePramPtr)
		  (covModelObj->*cov1ParamPtr)(phiMeans[k], mmblk[k*m+k], D[j*dnrow+i]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(phiMeans[k], nuMeans[k], mmblk[k*m+k], D[j*dnrow+i]);
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
	
	
	//use Psi
	if(!noPsi){
	  //get Psi sample
	  zeros(Psi, Alength);
	  if(PsiCase == 1){
	    for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiMeans[0]);}
	  }else if(PsiCase == 2){
	    for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiMeans[i]);}
	  }else{
	    setLowerChol(Psi, &PsiMeans[0], m);
	  }
	  
	  //get det
	  logDetIKPsi = 0.0;
	  for(i = 0; i < m; i++)
	    logDetIKPsi += log(Psi[i*m+i]);
	  logDetIKPsi = dnrow*2.0*logDetIKPsi;
	  
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
	  
	}//else IKPsi = I and logDetIKPsi = 0.0 set above
	
      }
      
      //chol decom
      if(!linpack){
	F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
      }else{
	F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
      }
      
      //get logDet
      logDetR = 0;
      for(i = 0; i < rnrow; i++)
	logDetR += log(R[i*rnrow+i]);
      logDetR = 2*logDetR;
      
      //finish the invert
      if(!linpack){
	F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
      }else{
	F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
      }
      
      for(i = 0; i < rlength; i++){
	if(isnan(R[i]))
	  error("c++ error: Lapack dpotri routine produced an undiscovered NAN in the R^{-1} computation. This is a rogue bug we have yet to track down.  Rerun with linpack=TRUE in the run.control list.  If this error persists when linpack=TRUE please report to package maintainer.\n");
      }
      
      
      F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, &betaMeans[0], 
		      &incOne, &zero, tmpXRow, &incOne);
      F77_NAME(daxpy)(&xnrow, &negOne, wMeans, &incOne, tmpXRow, &incOne);
      F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
      
      //(-1/2) * tmp` * Psi^{-1} * tmp
      F77_NAME(dsymv)(&upper, &xnrow, &one, IKPsi, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
      DBarUnmargOmega = logDetIKPsi+F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne);
      
      //get w likelihood
      //(-1/2) * w` * R`^{-1} * w
      F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, wMeans, &incOne, &zero, tmpXRow1, &incOne);
      DBarUnmargOmega += logDetR+F77_NAME(ddot)(&xnrow, wMeans, &incOne, tmpXRow1, &incOne);
      
      if(verbose){
	Rprintf("Sampled: %i of %i, %3.2f%%\n", nSamples, nSamples, 100.0);
        #ifdef Win32
	R_FlushConsole();
        #endif
	R_CheckUserInterrupt();
      }
      
    }
    
    
    /*********************
         return
    **********************/
    
    //make the result list object
    int nResultListObjs = 0;
    if(DICMarg) nResultListObjs++;
    if(DICUnmarg) nResultListObjs++;
    if(!haveW  && DICUnmarg) nResultListObjs++;
    
    SEXP dicResults, dicResultsRowNames, dicResultsColNames, dicResultDimName;
    PROTECT(dicResultsRowNames = allocVector(VECSXP, 4)); nProtect++;
    PROTECT(dicResultsColNames = allocVector(VECSXP, 1)); nProtect++;
    SET_VECTOR_ELT(dicResultsRowNames, 0, mkChar("bar.D"));
    SET_VECTOR_ELT(dicResultsRowNames, 1, mkChar("D.bar.Omega"));
    SET_VECTOR_ELT(dicResultsRowNames, 2, mkChar("pD"));
    SET_VECTOR_ELT(dicResultsRowNames, 3, mkChar("DIC"));
    SET_VECTOR_ELT(dicResultsColNames, 0, mkChar("value"));
    PROTECT(dicResultDimName = allocVector(VECSXP, 2)); nProtect++;
    SET_VECTOR_ELT(dicResultDimName, 0, dicResultsRowNames);
    SET_VECTOR_ELT(dicResultDimName, 1, dicResultsColNames);
    
    
    SEXP result, resultNames;
    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //set result list elements
    int resultListIndx = 0;
    if(DICMarg){
      SET_VECTOR_ELT(result, resultListIndx, DICMargResult);
      REAL(DICMargResult)[0] = DBarMarg;
      REAL(DICMargResult)[1] = DBarMargOmega;
      REAL(DICMargResult)[2] = DBarMarg - DBarMargOmega;
      REAL(DICMargResult)[3] = DBarMarg + DBarMarg - DBarMargOmega; 
      setAttrib(DICMargResult, R_DimNamesSymbol, dicResultDimName);
      
      SET_VECTOR_ELT(resultNames, resultListIndx, mkChar("DIC.marg")); 
      resultListIndx++;
    }
    
    if(DICUnmarg){
      SET_VECTOR_ELT(result, resultListIndx, DICUnmargResult);
      REAL(DICUnmargResult)[0] = DBarUnmarg;
      REAL(DICUnmargResult)[1] = DBarUnmargOmega;
      REAL(DICUnmargResult)[2] = DBarUnmarg - DBarUnmargOmega;
      REAL(DICUnmargResult)[3] = DBarUnmarg + DBarUnmarg - DBarUnmargOmega; 
      setAttrib(DICUnmargResult, R_DimNamesSymbol, dicResultDimName);

      SET_VECTOR_ELT(resultNames, resultListIndx, mkChar("DIC.unmarg")); 
      resultListIndx++;
    }
    
    if(!haveW  && DICUnmarg){//calculated it here
      //w observed
      SET_VECTOR_ELT(result, resultListIndx, w_r);
      SET_VECTOR_ELT(resultNames, resultListIndx, mkChar("sp.effects"));    
    }
    
    namesgets(result, resultNames);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    
  }
  
}
