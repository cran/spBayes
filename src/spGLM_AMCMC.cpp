#include <algorithm>
#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP spGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r, SEXP family_r, SEXP weights_r,
		   SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
		   SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP wStarting_r,
		   SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP nuTuning_r, SEXP betaTuning_r, SEXP wTuning_r,
		   SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int i,j,k,l,b,info,nProtect= 0;
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
    int pp = p*p;
    int n = INTEGER(n_r)[0];

    double *coordsD = REAL(coordsD_r);

    std::string family = CHAR(STRING_ELT(family_r,0));

    int *weights = INTEGER(weights_r);

    //covariance model
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));

    //priors and starting
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));

    double *betaMu = NULL;
    double *betaSd = NULL;
    
    if(betaPrior == "normal"){
      betaMu = REAL(VECTOR_ELT(betaNorm_r, 0)); 
      betaSd = REAL(VECTOR_ELT(betaNorm_r, 1));
    }
    
    double *sigmaSqIG = REAL(sigmaSqIG_r);
    double *phiUnif = REAL(phiUnif_r);

    double phiStarting = REAL(phiStarting_r)[0];
    double sigmaSqStarting = REAL(sigmaSqStarting_r)[0];
    double *betaStarting = REAL(betaStarting_r);
    double *wStarting = REAL(wStarting_r);

    double sigmaSqIGa = sigmaSqIG[0]; double sigmaSqIGb = sigmaSqIG[1];
    double phiUnifa = phiUnif[0]; double phiUnifb = phiUnif[1];

    //if matern
    double *nuUnif = NULL;
    double nuStarting = 0;
    double nuUnifa = 0, nuUnifb = 0;
    if(covModel == "matern"){
      nuUnif = REAL(nuUnif_r);
      nuStarting = REAL(nuStarting_r)[0];
      nuUnifa = nuUnif[0]; nuUnifb = nuUnif[1]; 
    }

    int nBatch = INTEGER(nBatch_r)[0];
    int batchLength = INTEGER(batchLength_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    int nSamples = nBatch*batchLength;

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

      if(betaPrior == "flat"){
	Rprintf("\tbeta flat.\n");
      }else{
	Rprintf("\tbeta normal:\n");
	Rprintf("\t\tmu:"); printVec(betaMu, p);
	Rprintf("\t\tsd:"); printVec(betaSd, p);Rprintf("\n");
      }
      Rprintf("\n");
  
      Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);
      Rprintf("\n");
      
      Rprintf("\tphi Unif hyperpriors a=%.5f and b=%.5f\n", phiUnifa, phiUnifb);
      Rprintf("\n");

      if(covModel == "matern"){
	Rprintf("\tnu Unif hyperpriors a=%.5f and b=%.5f\n", nuUnifa, nuUnifb);	  
	Rprintf("\n");
      }

      Rprintf("Adaptive Metropolis with target acceptance rate: %.1f\n", 100*acceptRate);
 
    } 

    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/
    int nn = n*n;

    //spatial parameters
    int nParams, betaIndx, sigmaSqIndx, phiIndx, nuIndx;

    if(covModel != "matern"){
      nParams = p+2;//sigma^2, phi
      betaIndx = 0; sigmaSqIndx = betaIndx+p; phiIndx = sigmaSqIndx+1;
    }else{
      nParams = p+3;//sigma^2, phi, nu
      betaIndx = 0; sigmaSqIndx = betaIndx+p; phiIndx = sigmaSqIndx+1; nuIndx = phiIndx+1;
    }

    //set starting
    double *spParams = (double *) R_alloc(nParams, sizeof(double));
    
    F77_NAME(dcopy)(&p, betaStarting, &incOne, &spParams[betaIndx], &incOne);

    spParams[sigmaSqIndx] = log(sigmaSqStarting);

    spParams[phiIndx] = logit(phiStarting, phiUnifa, phiUnifb);

    if(covModel == "matern") 
      spParams[nuIndx] = logit(nuStarting, nuUnifa, nuUnifb);

    double *w = (double *) R_alloc(n, sizeof(double));
    F77_NAME(dcopy)(&n, wStarting, &incOne, w, &incOne);

    //set tuning
    double *spTuning = (double *) R_alloc(nParams, sizeof(double));

    F77_NAME(dcopy)(&p, REAL(betaTuning_r), &incOne, &spTuning[betaIndx], &incOne);

    spTuning[sigmaSqIndx] = REAL(sigmaSqTuning_r)[0];

    spTuning[phiIndx] = REAL(phiTuning_r)[0];

    if(covModel == "matern") 
      spTuning[nuIndx] = REAL(nuTuning_r)[0];

    for(i = 0; i < nParams; i++)
      spTuning[i] = log(spTuning[i]);
    
    double *wTuning = (double *) R_alloc(n, sizeof(double));
    F77_NAME(dcopy)(&n, REAL(wTuning_r), &incOne, wTuning, &incOne);

    for(i = 0; i < n; i++)
      wTuning[i] = log(wTuning[i]);
    
    //return stuff  
    SEXP w_r, samples_r, accept_r, accept_w_r, tuning_r, tuning_w_r;
    
    PROTECT(w_r = allocMatrix(REALSXP, n, nSamples)); nProtect++;  
    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    PROTECT(accept_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++;
    PROTECT(accept_w_r = allocMatrix(REALSXP, n, nBatch)); nProtect++;
    PROTECT(tuning_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++;
    PROTECT(tuning_w_r = allocMatrix(REALSXP, n, nBatch)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0, spParamsjCurrent, wjCurrent;

    double *accept = (double *) R_alloc(nParams, sizeof(double)); zeros(accept, nParams);
    double *accept_w = (double *) R_alloc(n, sizeof(double)); zeros(accept_w, n);

    double *C = (double *) R_alloc(nn, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_n1 = (double *) R_alloc(n, sizeof(double));

    double sigmaSq, phi, nu;
    double *beta = (double *) R_alloc(p, sizeof(double));
    double *theta = (double *) R_alloc(3, sizeof(double)); //phi, nu, and perhaps more in the future

    double logMHRatio;

    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    logPostCurrent = R_NegInf;

    GetRNGstate();
    for(b = 0, s = 0; b < nBatch; b++){
      
      for(i = 0; i < batchLength; i++, s++){
	
	for(j = 0; j < nParams; j++){
	  
	  //propose
	  spParamsjCurrent = spParams[j];
	  spParams[j] = rnorm(spParamsjCurrent, exp(spTuning[j]));
 
	  //extract and transform
	  F77_NAME(dcopy)(&p, &spParams[betaIndx], &incOne, beta, &incOne);
	  sigmaSq = theta[0] = exp(spParams[sigmaSqIndx]);
	  phi = theta[1] = logitInv(spParams[phiIndx], phiUnifa, phiUnifb);
	  
	  if(covModel == "matern"){
	    nu = theta[2] = logitInv(spParams[nuIndx], nuUnifa, nuUnifb);
	  }
	  
	  //construct covariance matrix
	  spCovLT(coordsD, n, theta, covModel, C);

	  //invert C and log det cov
	  detCand = 0;
	  F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in spGLM\n");}
	  for(k = 0; k < n; k++) detCand += 2*log(C[k*n+k]);
	  F77_NAME(dpotri)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spGLM\n");}
	  
	  //Likelihood with Jacobian  
	  logPostCand = 0.0;
	  
	  if(betaPrior == "normal"){
	    for(k = 0; k < p; k++){
	      logPostCand += dnorm(beta[k], betaMu[k], betaSd[k], 1);
	    }
	  }
	  
	  logPostCand += -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq); 
	  
	  logPostCand += log(phi - phiUnifa) + log(phiUnifb - phi); 
	  
	  if(covModel == "matern"){
	    logPostCand += log(nu - nuUnifa) + log(nuUnifb - nu);   
	  }
	  
	  F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
	  
	  if(family == "binomial"){
	    logPostCand += binomial_logpost(n, Y, tmp_n, w, weights);
	  }else if(family == "poisson"){
	    logPostCand += poisson_logpost(n, Y, tmp_n, w, weights);
	  }else{
	    error("c++ error: family misspecification in spGLM\n");
	  }
	  
	  //(-1/2) * tmp_n` *  C^-1 * tmp_n
	  F77_NAME(dsymv)(lower, &n, &one,  C, &n, w, &incOne, &zero, tmp_n1, &incOne);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&n, w, &incOne, tmp_n1, &incOne);

	  //
	  //MH accept/reject	
	  //      
	  
	  //MH ratio with adjustment
	  logMHRatio = logPostCand - logPostCurrent;
	  
	  if(runif(0.0,1.0) <= exp(logMHRatio)){
	    logPostCurrent = logPostCand;
	    accept[j]++;
	  }else{
	    spParams[j] = spParamsjCurrent;
	  }

	}//end spParams

	/************************/
	//update the sp. matrices
	/************************/

	//extract and transform
	F77_NAME(dcopy)(&p, &spParams[betaIndx], &incOne, beta, &incOne);
	sigmaSq = theta[0] = exp(spParams[sigmaSqIndx]);
	phi = theta[1] = logitInv(spParams[phiIndx], phiUnifa, phiUnifb);
	
	if(covModel == "matern"){
	  nu = theta[2] = logitInv(spParams[nuIndx], nuUnifa, nuUnifb);
	}
	
	//construct covariance matrix
	spCovLT(coordsD, n, theta, covModel, C);

	//invert C and log det cov
	detCand = 0;
	F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in spGLM\n");}
	for(k = 0; k < n; k++) detCand += 2*log(C[k*n+k]);
	F77_NAME(dpotri)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spGLM\n");}
	
	F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
	  
	for(j = 0; j < n; j++){
	  
	  //propose
	  wjCurrent = w[j];
	  w[j] = rnorm(wjCurrent, exp(wTuning[j]));

	  //Likelihood with Jacobian  
	  logPostCand = 0.0;
	  
	  if(betaPrior == "normal"){
	    for(k = 0; k < p; k++){
	      logPostCand += dnorm(beta[k], betaMu[k], betaSd[k], 1);
	    }
	  }
	  
	  logPostCand += -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq); 
	  
	  logPostCand += log(phi - phiUnifa) + log(phiUnifb - phi); 
	  
	  if(covModel == "matern"){
	    logPostCand += log(nu - nuUnifa) + log(nuUnifb - nu);   
	  }
	  	  
	  if(family == "binomial"){
	    logPostCand += binomial_logpost(n, Y, tmp_n, w, weights);
	  }else if(family == "poisson"){
	    logPostCand += poisson_logpost(n, Y, tmp_n, w, weights);
	  }else{
	    error("c++ error: family misspecification in spGLM\n");
	  }
	  
	  //(-1/2) * tmp_n` *  C^-1 * tmp_n
	  F77_NAME(dsymv)(lower, &n, &one,  C, &n, w, &incOne, &zero, tmp_n1, &incOne);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&n, w, &incOne, tmp_n1, &incOne);
	  
	  //
	  //MH accept/reject	
	  //      
	  
	  //MH ratio with adjustment
	  logMHRatio = logPostCand - logPostCurrent;
	  
	  if(runif(0.0,1.0) <= exp(logMHRatio)){
	    logPostCurrent = logPostCand;
	    accept_w[j]++;
	  }else{
	    w[j] = wjCurrent;
	  }
	  
	}//end w
	
	
	/******************************
               Save samples
	*******************************/
	F77_NAME(dcopy)(&nParams, spParams, &incOne, &REAL(samples_r)[s*nParams], &incOne);
	F77_NAME(dcopy)(&n, w, &incOne, &REAL(w_r)[s*n], &incOne);
	      
	R_CheckUserInterrupt();

      }//end batch
      
      //adjust tuning
      for(j = 0; j < nParams; j++){
	REAL(accept_r)[b*nParams+j] = accept[j]/batchLength;
	REAL(tuning_r)[b*nParams+j] = spTuning[j];
	
	if(accept[j]/batchLength > acceptRate){
	  spTuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}else{
	  spTuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}
	accept[j] = 0.0;
      }
      
      for(j = 0; j < n; j++){
	REAL(accept_w_r)[b*n+j] = accept_w[j]/batchLength;
	REAL(tuning_w_r)[b*n+j] = wTuning[j];
	
	if(accept_w[j]/batchLength > acceptRate){
	  wTuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}else{
	  wTuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}
	accept_w[j] = 0.0;
      }
    
      //report
      if(verbose){
	if(status == nReport){
	  Rprintf("Batch: %i of %i, %3.2f%%\n", b, nBatch, 100.0*b/nBatch);
	  Rprintf("\tparameter\tacceptance\ttuning\n");	  
	  for(j = 0; j < p; j++){
	    Rprintf("\tbeta[%i]\t\t%3.1f\t\t%1.5f\n", j, 100.0*REAL(accept_r)[b*nParams+betaIndx+j], exp(spTuning[betaIndx+j]));
	  }
	  Rprintf("\tsigma.sq\t%3.1f\t\t%1.5f\n", 100.0*REAL(accept_r)[b*nParams+sigmaSqIndx], exp(spTuning[sigmaSqIndx]));
	  Rprintf("\tphi\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(accept_r)[b*nParams+phiIndx], exp(spTuning[phiIndx]));
	  if(covModel == "matern")
	    Rprintf("\tnu\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(accept_r)[b*nParams+nuIndx], exp(spTuning[nuIndx]));
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      
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
      REAL(samples_r)[s*nParams+sigmaSqIndx] = exp(REAL(samples_r)[s*nParams+sigmaSqIndx]);

      REAL(samples_r)[s*nParams+phiIndx] = logitInv(REAL(samples_r)[s*nParams+phiIndx], phiUnifa, phiUnifb);

      if(covModel == "matern")
	REAL(samples_r)[s*nParams+nuIndx] = logitInv(REAL(samples_r)[s*nParams+nuIndx], nuUnifa, nuUnifb);
    }
   
    //make return object
    SEXP result, resultNames;
    
    int nResultListObjs = 6;
    
    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;

    //samples
    SET_VECTOR_ELT(result, 0, samples_r);
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.beta.theta.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));

    SET_VECTOR_ELT(result, 2, accept_w_r);
    SET_VECTOR_ELT(resultNames, 2, mkChar("acceptance.w"));
    
    SET_VECTOR_ELT(result, 3, w_r);
    SET_VECTOR_ELT(resultNames, 3, mkChar("p.w.samples"));

    SET_VECTOR_ELT(result, 4, tuning_r);
    SET_VECTOR_ELT(resultNames, 4, mkChar("tuning"));

    SET_VECTOR_ELT(result, 5, tuning_w_r);
    SET_VECTOR_ELT(resultNames, 5, mkChar("tuning.w"));
  
    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);

    
  }
}
