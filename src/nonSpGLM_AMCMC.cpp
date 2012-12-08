#include <algorithm>
#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP nonSpGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r, SEXP weights_r,
		      SEXP betaPrior_r, SEXP betaNorm_r, SEXP betaStarting_r, SEXP betaTuning_r, 
		      SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r, SEXP amcmc_r){
    
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

    std::string family = CHAR(STRING_ELT(family_r,0));
    int *weights = INTEGER(weights_r);

    //priors and starting
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));

    double *betaMu = NULL;
    double *betaSd = NULL;
    
    if(betaPrior == "normal"){
      betaMu = REAL(VECTOR_ELT(betaNorm_r, 0)); 
      betaSd = REAL(VECTOR_ELT(betaNorm_r, 1));
    }
    
    double *betaStarting = REAL(betaStarting_r);
  
    int nBatch = INTEGER(nBatch_r)[0];
    int batchLength = INTEGER(batchLength_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    int nSamples = nBatch*batchLength;

    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    bool amcmc = static_cast<bool>(INTEGER(amcmc_r)[0]);

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
  
      Rprintf("Priors and hyperpriors:\n");

      if(betaPrior == "flat"){
	Rprintf("\tbeta flat.\n");
      }else{
	Rprintf("\tbeta normal:\n");
	Rprintf("\t\tmu:"); printVec(betaMu, p);
	Rprintf("\t\tsd:"); printVec(betaSd, p);Rprintf("\n");
      }
      Rprintf("\n");
  
      if(amcmc){
	Rprintf("Using adaptive MCMC.\n\n");
	Rprintf("\tNumber of batches %i.\n", nBatch);
	Rprintf("\tBatch length %i.\n", batchLength);
	Rprintf("\tTarget acceptance rate %.5f.\n", acceptRate);
	Rprintf("\n");
      }else{
	Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      }
 
    } 

    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/  
    //set starting
    double *params = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, betaStarting, &incOne, params, &incOne);

    //set tuning
    double *tuning = NULL;
    
    if(amcmc){
      tuning =(double *) R_alloc(p, sizeof(double));
      F77_NAME(dcopy)(&p, REAL(betaTuning_r), &incOne, tuning, &incOne);
      
      for(i = 0; i < p; i++){
	tuning[i] = log(tuning[i]);
      }
    }else{
      tuning =(double *) R_alloc(p*p, sizeof(double));
      F77_NAME(dcopy)(&pp, REAL(betaTuning_r), &incOne, tuning, &incOne);
    }

    //return stuff  
    SEXP samples_r, accept_r, tuning_r;
    
    PROTECT(samples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++; 
    PROTECT(accept_r = allocMatrix(REALSXP, p, nBatch)); nProtect++;
    PROTECT(tuning_r = allocMatrix(REALSXP, p, nBatch)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, batchAccept=0;
    double logPostCurrent = R_NegInf, logPostCand = 0,  paramsjCurrent = 0;

    double *paramsCurrent = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, params, &incOne, paramsCurrent, &incOne);

    double *accept = (double *) R_alloc(p, sizeof(double)); zeros(accept, p);
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *beta = (double *) R_alloc(p, sizeof(double));
   
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
	
	for(j = 0; j < p; j++){
	  
	  //propose
	  if(amcmc){
	    paramsjCurrent = params[j];
	    params[j] = rnorm(paramsjCurrent, exp(tuning[j]));
	  }else{
	    mvrnorm(params, paramsCurrent, tuning, p, false);
	  }
	  
	  F77_NAME(dcopy)(&p, params, &incOne, beta, &incOne);
	  
	  //Likelihood with Jacobian  
	  logPostCand = 0.0;
	  
	  if(betaPrior == "normal"){
	    for(k = 0; k < p; k++){
	      logPostCand += dnorm(beta[k], betaMu[k], betaSd[k], 1);
	    }
	  }
	  	  
	  F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
	  
	  if(family == "binomial"){
	    logPostCand += binomial_logpost(n, Y, tmp_n, weights);
	  }else if(family == "poisson"){
	    logPostCand += poisson_logpost(n, Y, tmp_n, weights);
	  }else{
	    error("c++ error: family misspecification in spGLM\n");
	  }
	  
	  //
	  //MH accept/reject	
	  //      
	  
	  //MH ratio with adjustment
	  logMHRatio = logPostCand - logPostCurrent;
	  
	  if(runif(0.0,1.0) <= exp(logMHRatio)){
	    logPostCurrent = logPostCand;
	    
	    if(amcmc){
	      accept[j]++;
	    }else{
	      accept[0]++;
	      batchAccept++;
	    }

	  }else{

	    if(amcmc){
	      params[j] = paramsjCurrent;
	    }else{
	      F77_NAME(dcopy)(&p, paramsCurrent, &incOne, params, &incOne);
	    }
	  }
	  
	  if(!amcmc){
	    break;
	  }

	}//end params
	
	/******************************
               Save samples
	*******************************/
	F77_NAME(dcopy)(&p, params, &incOne, &REAL(samples_r)[s*p], &incOne);
	
	R_CheckUserInterrupt();
	
      }//end batch
      
      //adjust tuning      
      if(amcmc){
	for(j = 0; j < p; j++){
	  REAL(accept_r)[b*p+j] = accept[j]/batchLength;
	  REAL(tuning_r)[b*p+j] = tuning[j];
	  
	  if(accept[j]/batchLength > acceptRate){
	    tuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	  }else{
	    tuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	  }
	  accept[j] = 0.0;
	}
      }
      
      //report
      if(verbose){
	if(status == nReport){
	  if(amcmc){
	    Rprintf("Batch: %i of %i, %3.2f%%\n", b, nBatch, 100.0*b/nBatch);
	    Rprintf("\tparameter\tacceptance\ttuning\n");	  
	    for(j = 0; j < p; j++){
	      Rprintf("\tbeta[%i]\t\t%3.1f%\t\t%1.5f\n", j, 100.0*REAL(accept_r)[b*p+j], exp(tuning[j]));
	    }
	  }else{
	    Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
      	    Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
      	    Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept[0]/s);
      	  }
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
      	  batchAccept = 0;
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

    //make return object
    SEXP result, resultNames;
    
    int nResultListObjs = 3;
    
    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;

    //samples
    SET_VECTOR_ELT(result, 0, samples_r);
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.beta.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));

    SET_VECTOR_ELT(result, 2, tuning_r);
    SET_VECTOR_ELT(resultNames, 2, mkChar("tuning"));
  
    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
  }
}
