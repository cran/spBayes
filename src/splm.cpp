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

  SEXP spLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
	    SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	    SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r,
	    SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP nuTuning_r, 
	    SEXP nugget_r, SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r){


    /*****************************************
                Common variables
    *****************************************/
    int i,j,k,l,info,nProtect= 0;
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

    /*****************************************
                     Set-up
    *****************************************/

    double *Y = REAL(Y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];

    double *coordsD = REAL(coordsD_r);

    //covariance model
    string covModel = CHAR(STRING_ELT(covModel_r,0));

    //priors and starting
    double *sigmaSqIG = REAL(sigmaSqIG_r);
    double *phiUnif = REAL(phiUnif_r);

    double phiStarting = REAL(phiStarting_r)[0];
    double sigmaSqStarting = REAL(sigmaSqStarting_r)[0];
    double *betaStarting = REAL(betaStarting_r);

    double sigmaSqIGa = sigmaSqIG[0]; double sigmaSqIGb = sigmaSqIG[1];
    double phiUnifa = phiUnif[0]; double phiUnifb = phiUnif[1];

    //if nugget
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    double *tauSqIG = NULL;
    double tauSqStarting = 0;
    double tauSqIGa = 0, tauSqIGb = 0;
    if(nugget){
      tauSqIG = REAL(tauSqIG_r);
      tauSqStarting = REAL(tauSqStarting_r)[0];
      tauSqIGa = tauSqIG[0]; tauSqIGb = tauSqIG[1]; 
    }

    //if matern
    double *nuUnif = NULL;
    double nuStarting = 0;
    double nuUnifa = 0, nuUnifb = 0;
    if(covModel == "matern"){
      nuUnif = REAL(nuUnif_r);
      nuStarting = REAL(nuStarting_r)[0];
      nuUnifa = nuUnif[0]; nuUnifb = nuUnif[1]; 
    }

    //tuning
    double phiTuning = sqrt(REAL(phiTuning_r)[0]);
    double sigmaSqTuning = sqrt(REAL(sigmaSqTuning_r)[0]);
    double tauSqTuning = 0;
    double nuTuning = 0;
    if(nugget)
      tauSqTuning = sqrt(REAL(tauSqTuning_r)[0]);

    if(covModel == "matern")
      nuTuning = sqrt(REAL(nuTuning_r)[0]);

    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];


    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      

      Rprintf("Number of MCMC samples %i.\n\n", nSamples);

      Rprintf("Priors and hyperpriors:\n");
      Rprintf("\tbeta flat.\n");   
      Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);
      if(nugget){
	Rprintf("\ttau.sq IG hyperpriors shape=%.5f and scale=%.5f\n", tauSqIGa, tauSqIGb); 
      }
      Rprintf("\tphi Unif hyperpriors a=%.5f and b=%.5f\n", phiUnifa, phiUnifb);
      if(covModel == "matern"){
	Rprintf("\tnu Unif hyperpriors a=%.5f and b=%.5f\n", nuUnifa, nuUnifb);	  
      }
    } 

    /*****************************************
        Set-up cov. model function pointer
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
         Set-up MCMC sample matrices etc.
    *****************************************/
    int nn = n*n;

    //spatial parameters
    int nSpParams, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;

    if(!nugget && covModel != "matern"){
      nSpParams = 2;//sigma^2, phi
      sigmaSqIndx = 0; phiIndx = 1;
    }else if(nugget && covModel != "matern"){
      nSpParams = 3;//sigma^2, tau^2, phi
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
    }else if(!nugget && covModel == "matern"){
      nSpParams = 3;//sigma^2, phi, nu
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
    }else{
      nSpParams = 4;//sigma^2, tau^2, phi, nu
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;//sigma^2, tau^2, phi, nu
    }
    
    double *spParams = (double *) R_alloc(nSpParams, sizeof(double));
    
    //set starting
    spParams[sigmaSqIndx] = log(sigmaSqStarting);

    if(nugget) spParams[tauSqIndx] = log(tauSqStarting);

    spParams[phiIndx] = logit(phiStarting, phiUnifa, phiUnifb);

    if(covModel == "matern") 
      spParams[nuIndx] = logit(nuStarting, nuUnifa, nuUnifb);

    //Beta parameter and set starting
    double *beta = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, betaStarting, &incOne, beta, &incOne);

    //samples and random effects
    int nParams = p+nSpParams;

    SEXP w_r, samples_r, accept_r;

    PROTECT(w_r = allocMatrix(REALSXP, n, nSamples)); nProtect++; 
    double *w = REAL(w_r); zeros(w, n*nSamples);

    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    double *samples = REAL(samples_r);

    PROTECT(accept_r = allocMatrix(REALSXP, 1, 1)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, rtnStatus=0, accept=0, batchAccept = 0;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0, det = 0;
    bool first = true, accepted = true;

    double *C = (double *) R_alloc(nn, sizeof(double)); 
    double *wMu = (double *) R_alloc(n, sizeof(double));
    
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_n1 = (double *) R_alloc(n, sizeof(double));

    double *candSpParams = (double *) R_alloc(nSpParams, sizeof(double));
    double sigmaSq, tauSq, phi, nu, sigmaSqTmp, tauSqTmp;
    double logMHRatio, HastAdj;
    
    double *S_beta = (double *) R_alloc(p*p, sizeof(double));
    double *Mu_beta = (double *) R_alloc(p, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double)); 
    double *tmp_np = (double *) R_alloc(n*p, sizeof(double)); 

    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    GetRNGstate();
    for(s = 0; s < nSamples; s++){
      
      //
      //Current
      //
      sigmaSq = exp(spParams[sigmaSqIndx]);
      
      if(nugget){
	tauSq = exp(spParams[tauSqIndx]);
      }      

      phi = logitInv(spParams[phiIndx], phiUnifa, phiUnifb);

      if(covModel == "matern"){
	nu = logitInv(spParams[nuIndx], nuUnifa, nuUnifb);
      }

      
      //make the correlation matrix
      for(i = 0; i < nn; i++){
	if(onePramPtr)
	  (covModelObj->*cov1ParamPtr)(phi, C[i], coordsD[i]);
	else //i.e., 2 parameter matern
	  (covModelObj->*cov2ParamPtr)(phi, nu, C[i], coordsD[i]);
      }
      
      F77_NAME(dscal)(&nn, &sigmaSq, C, &incOne);
      
      if(nugget){
	for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+tauSq;
      }
      
      //invert C and log det cov
      F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
      
      det = 0;
      for(i = 0; i < n; i++) det += 2*log(C[i*n+i]);

      F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
      
      //
      //Update Beta
      //
      //finish the Gibbs
      //C^-1 X = tmp_np
      F77_NAME(dsymm)(lside, upper, &n, &p, &one, C, &n, X, &n, &zero, tmp_np, &n);
      
      //(t(X) tmp_np)^{-1} = S_beta
      F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, tmp_np, &n, &zero, S_beta, &p);
      
      F77_NAME(dpotrf)(upper, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
      F77_NAME(dpotri)(upper, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky inverse failed\n" << endl;}
      
      //S_beta tmp_np Y = Mu_beta
      F77_NAME(dsymv)(upper, &n, &one, C, &n, Y, &incOne, &zero, tmp_n, &incOne);
      F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, tmp_n, &incOne, &zero, tmp_p, &incOne);
      F77_NAME(dsymv)(upper, &p, &one, S_beta, &p, tmp_p, &incOne, &zero, Mu_beta, &incOne);
      
      //Gibbs draw
      //take upper for the chol for the mv draw
      F77_NAME(dpotrf)(upper, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
      
      mvrnorm(beta, Mu_beta, S_beta, p, true);

      //
      //Likelihood with Jacobian   
      //
      logPostCurrent = -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq);
      
      if(nugget){
	logPostCurrent += -1.0*(1.0+tauSqIGa)*log(tauSq)-tauSqIGb/tauSq+log(tauSq);
      }
      
      logPostCurrent += log(phi - phiUnifa) + log(phiUnifb - phi); 

      if(covModel == "matern"){
	logPostCurrent += log(nu - nuUnifa) + log(nuUnifb - nu);   
      }

      //Y-XB
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
      
      //(-1/2) * tmp_n` *  C^-1 * tmp_n
      F77_NAME(dsymv)(upper, &n, &one,  C, &n, tmp_n, &incOne, &zero, tmp_n1, &incOne);
      logPostCurrent += -0.5*det-0.5*F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n1, &incOne);
    
      
      //
      //Candidate
      //

      //propose   
      candSpParams[sigmaSqIndx] = rnorm(spParams[sigmaSqIndx], sigmaSqTuning);
      sigmaSq = exp(candSpParams[sigmaSqIndx]);

      if(nugget){
	candSpParams[tauSqIndx] = rnorm(spParams[tauSqIndx], tauSqTuning);
	tauSq = exp(candSpParams[tauSqIndx]);
      }      

      candSpParams[phiIndx] = rnorm(spParams[phiIndx], phiTuning);
      phi = logitInv(candSpParams[phiIndx], phiUnifa, phiUnifb);

      if(covModel == "matern"){
	candSpParams[nuIndx] = rnorm(spParams[nuIndx], nuTuning);
	nu = logitInv(candSpParams[nuIndx], nuUnifa, nuUnifb);
      }

      
      //make the correlation matrix
      for(i = 0; i < nn; i++){
	if(onePramPtr)
	  (covModelObj->*cov1ParamPtr)(phi, C[i], coordsD[i]);
	else //i.e., 2 parameter matern
	  (covModelObj->*cov2ParamPtr)(phi, nu, C[i], coordsD[i]);
      }
      
      F77_NAME(dscal)(&nn, &sigmaSq, C, &incOne);
      
      if(nugget){
	for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+tauSq;
      }
      
      //invert C and log det cov
      F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
      
      detCand = 0;
      for(i = 0; i < n; i++) detCand += 2*log(C[i*n+i]);

      F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
      

      //Likelihood with Jacobian   
      logPostCand = -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq);
      
      if(nugget){
	logPostCand += -1.0*(1.0+tauSqIGa)*log(tauSq)-tauSqIGb/tauSq+log(tauSq);
      }
      
      logPostCand += log(phi - phiUnifa) + log(phiUnifb - phi); 

      if(covModel == "matern"){
	logPostCand += log(nu - nuUnifa) + log(nuUnifb - nu);   
      }

      //Y-XB
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
      
      //(-1/2) * tmp_n` *  C^-1 * tmp_n
      F77_NAME(dsymv)(upper, &n, &one,  C, &n, tmp_n, &incOne, &zero, tmp_n1, &incOne);
      logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n1, &incOne);
      
      //
      //MH accept/reject	
      //      
  
      //MH ratio with adjustment
      logMHRatio = logPostCand - logPostCurrent;
      
      if(runif(0.0,1.0) <= exp(logMHRatio)){
	F77_NAME(dcopy)(&nSpParams, candSpParams, &incOne, spParams, &incOne);
	accept++;
	batchAccept++;
      }
              
      /******************************
         Recover w and w* if needed
      *******************************/
      
      //Get previously accepted sample
      sigmaSq = exp(spParams[sigmaSqIndx]);
      
      if(nugget)
	tauSq = exp(spParams[tauSqIndx]);
      
      phi = logitInv(spParams[phiIndx], phiUnifa, phiUnifb);
      
      if(covModel == "matern")
	nu = logitInv(spParams[nuIndx], nuUnifa, nuUnifb);
      
      if(nugget){
	//make the correlation matrix
	for(i = 0; i < nn; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phi, C[i], coordsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phi, nu, C[i], coordsD[i]);
	}
	
	//invert C
	F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	
	//scale correlation matrix with 1/sigmasq and add 1/nugget to diag
	sigmaSqTmp = 1.0/sigmaSq;
	F77_NAME(dscal)(&nn, &sigmaSqTmp, C, &incOne);
	
	for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+1.0/tauSq;
	
	//invert C
	F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	
	//make w mu
	F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
	F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	
	tauSqTmp = 1.0/tauSq;
	F77_NAME(dscal)(&n, &tauSqTmp, tmp_n, &incOne);
	
	F77_NAME(dsymv)(upper, &n, &one, C, &n, tmp_n, &incOne, &zero, wMu, &incOne);
	
	//chol for the mvnorm and draw
	F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	
	mvrnorm(&w[s*n], wMu, C, n, true);
	
      }else{//no nugget so w is just resids
	F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, &w[s*n], &incOne);
	F77_NAME(daxpy)(&n, &one, Y, &incOne, &w[s*n], &incOne);
      }
      
      
      /******************************
          Save samples and report
      *******************************/
      F77_NAME(dcopy)(&p, beta, &incOne, &samples[s*nParams], &incOne);
      F77_NAME(dcopy)(&nSpParams, spParams, &incOne, &samples[s*nParams+p], &incOne);
      
      //report
      if(verbose){
	if(status == nReport){
	  Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
	  Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
	  Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept/s);
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	  batchAccept = 0;
	}
      }
      status++;
      
      
      R_CheckUserInterrupt();
    }//end sample loop
    PutRNGstate();
    
    //final status report
    if(verbose){
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
    }
    Rprintf("-------------------------------------------------\n");
    #ifdef Win32
    R_FlushConsole();
    #endif
    
    //untransform variance variables
    for(s = 0; s < nSamples; s++){
      samples[s*nParams+p+sigmaSqIndx] = exp(samples[s*nParams+p+sigmaSqIndx]);
     
      if(nugget)
	samples[s*nParams+p+tauSqIndx] = exp(samples[s*nParams+p+tauSqIndx]);

      samples[s*nParams+p+phiIndx] = logitInv(samples[s*nParams+p+phiIndx], phiUnifa, phiUnifb);

      if(covModel == "matern")
	samples[s*nParams+p+nuIndx] = logitInv(samples[s*nParams+p+nuIndx], nuUnifa, nuUnifb);
    }
   
    //calculate acceptance rate
    REAL(accept_r)[0] = 100.0*accept/s;

    //make return object
    SEXP result, resultNames;
    
    int nResultListObjs = 3;

    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;

   //samples
    SET_VECTOR_ELT(result, 0, samples_r);
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));
    
    SET_VECTOR_ELT(result, 2, w_r);
    SET_VECTOR_ELT(resultNames, 2, mkChar("sp.effects"));
  
    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    
  }
}
