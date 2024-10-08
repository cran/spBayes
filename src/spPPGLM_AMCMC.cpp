#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#define USE_FC_LEN_T
#include <algorithm>
#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"
#ifndef FCONE
# define FCONE
#endif

extern "C" {

  SEXP spPPGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP family_r, SEXP weights_r,
		     SEXP m_r, SEXP knotsD_r, SEXP knotsCoordsD_r, 
		     SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
		     SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP w_strStarting_r,
		     SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP nuTuning_r, SEXP betaTuning_r, SEXP w_strTuning_r,
		     SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int i, j, k, b, info, nProtect= 0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *lside = "L";
    const double one = 1.0;
    const double zero = 0.0;
    const int incOne = 1;

    /*****************************************
                     Set-up
    *****************************************/

    double *Y = REAL(Y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];

    std::string family = CHAR(STRING_ELT(family_r,0));

    int *weights = INTEGER(weights_r);

    //covariance model
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));

    int m = INTEGER(m_r)[0];
    double *knotsD = REAL(knotsD_r);
    double *knotsCoordsD = REAL(knotsCoordsD_r);

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
    double *w_strStarting = REAL(w_strStarting_r);

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
    double acceptRate  = REAL(acceptRate_r)[0];
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
      
      Rprintf("Using non-modified predictive process with %i knots.\n\n", m);
    
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
    int nm = n*m, mm = m*m;

    //spatial parameters
    int nParams, betaIndx, sigmaSqIndx, phiIndx, nuIndx = 0;

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

    double *w_str = (double *) R_alloc(m, sizeof(double));
    F77_NAME(dcopy)(&m, w_strStarting, &incOne, w_str, &incOne);

    //set tuning
    double *spTuning = (double *) R_alloc(nParams, sizeof(double));

    F77_NAME(dcopy)(&p, REAL(betaTuning_r), &incOne, &spTuning[betaIndx], &incOne);

    spTuning[sigmaSqIndx] = REAL(sigmaSqTuning_r)[0];

    spTuning[phiIndx] = REAL(phiTuning_r)[0];

    if(covModel == "matern") 
      spTuning[nuIndx] = REAL(nuTuning_r)[0];

    for(i = 0; i < nParams; i++)
      spTuning[i] = log(spTuning[i]);
    
    double *w_strTuning = (double *) R_alloc(m, sizeof(double));
    F77_NAME(dcopy)(&m, REAL(w_strTuning_r), &incOne, w_strTuning, &incOne);

    for(i = 0; i < m; i++)
      w_strTuning[i] = log(w_strTuning[i]);
    
    //return stuff  
    SEXP w_r, w_str_r, samples_r, accept_r, accept_w_str_r, tuning_r, tuning_w_str_r;
    
    PROTECT(w_r = Rf_allocMatrix(REALSXP, n, nSamples)); nProtect++; 
    PROTECT(w_str_r = Rf_allocMatrix(REALSXP, m, nSamples)); nProtect++; 
    PROTECT(samples_r = Rf_allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    PROTECT(accept_r = Rf_allocMatrix(REALSXP, nParams, nBatch)); nProtect++;
    PROTECT(accept_w_str_r = Rf_allocMatrix(REALSXP, m, nBatch)); nProtect++;
    PROTECT(tuning_r = Rf_allocMatrix(REALSXP, nParams, nBatch)); nProtect++;
    PROTECT(tuning_w_str_r = Rf_allocMatrix(REALSXP, m, nBatch)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0, spParamsjCurrent, w_strjCurrent;

    double *accept = (double *) R_alloc(nParams, sizeof(double)); zeros(accept, nParams);
    double *accept_w_str = (double *) R_alloc(m, sizeof(double)); zeros(accept_w_str, m);

    double *K = (double *) R_alloc(mm, sizeof(double));
    double *P = (double *) R_alloc(nm, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
     
    double sigmaSq, phi, nu = 0;
    double *beta = (double *) R_alloc(p, sizeof(double));
    double *w = (double *) R_alloc(n, sizeof(double)); zeros(w, n);
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
	  
	  //construct covariance matrices 
	  spCovLT(knotsD, m, theta, covModel, K);
	  spCov(knotsCoordsD, nm, theta, covModel, P);
  
	  //invert C and log det cov
	  detCand = 0;
	  F77_NAME(dpotrf)(lower, &m, K, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: Cholesky failed in spGLM\n");}
	  for(k = 0; k < m; k++) detCand += 2*log(K[k*m+k]);
	  F77_NAME(dpotri)(lower, &m, K, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: Cholesky inverse failed in spGLM\n");}
	  
	  //make \tild{w}
	  F77_NAME(dsymv)(lower, &m, &one, K, &m, w_str, &incOne, &zero, tmp_m, &incOne FCONE);     
	  F77_NAME(dgemv)(ytran, &m, &n, &one, P, &m, tmp_m, &incOne, &zero, w, &incOne FCONE);
	  
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
	  
	  F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &incOne, &zero, tmp_n, &incOne FCONE);
	  
	  if(family == "binomial"){
	    logPostCand += binomial_logpost(n, Y, tmp_n, w, weights);
	  }else if(family == "poisson"){
	    logPostCand += poisson_logpost(n, Y, tmp_n, w, weights);
	  }else{
	    Rf_error("c++ Rf_error: family misspecification in spGLM\n");
	  }
	  
	  //(-1/2) * tmp_n` *  C^-1 * tmp_n
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&m, w_str, &incOne, tmp_m, &incOne);

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
	
	//construct covariance matrices 
	spCovLT(knotsD, m, theta, covModel, K);
	spCov(knotsCoordsD, nm, theta, covModel, P);
	
	//invert C and log det cov
	detCand = 0;
	F77_NAME(dpotrf)(lower, &m, K, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: Cholesky failed in spGLM\n");}
	for(k = 0; k < m; k++) detCand += 2*log(K[k*m+k]);
	F77_NAME(dpotri)(lower, &m, K, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: Cholesky inverse failed in spGLM\n");}
		
	F77_NAME(dsymm)(lside, lower, &m, &n, &one, K, &m, P, &m, &zero, tmp_nm, &m FCONE FCONE);

	F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &incOne, &zero, tmp_n, &incOne FCONE);

	for(j = 0; j < m; j++){
	  
	  //propose
	  w_strjCurrent = w_str[j];
	  w_str[j] = rnorm(w_strjCurrent, exp(w_strTuning[j]));

	  //make \tild{w}
	  F77_NAME(dgemv)(ytran, &m, &n, &one, tmp_nm, &m, w_str, &incOne, &zero, w, &incOne FCONE);
	  
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
	    Rf_error("c++ Rf_error: family misspecification in spGLM\n");
	  }
	  
	  //(-1/2) * tmp_n` *  C^-1 * tmp_n
	  F77_NAME(dsymv)(lower, &m, &one, K, &m, w_str, &incOne, &zero, tmp_m, &incOne FCONE);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&m, w_str, &incOne, tmp_m, &incOne);
	  
	  //
	  //MH accept/reject	
	  //      
	  
	  //MH ratio with adjustment
	  logMHRatio = logPostCand - logPostCurrent;
	  
	  if(runif(0.0,1.0) <= exp(logMHRatio)){
	    logPostCurrent = logPostCand;
	    accept_w_str[j]++;
	  }else{
	    w_str[j] = w_strjCurrent;
	  }
	  
	}//end w_str
	
	
	/******************************
               Save samples
	*******************************/
	F77_NAME(dcopy)(&nParams, spParams, &incOne, &REAL(samples_r)[s*nParams], &incOne);
	F77_NAME(dcopy)(&n, w, &incOne, &REAL(w_r)[s*n], &incOne);
	F77_NAME(dcopy)(&m, w_str, &incOne, &REAL(w_str_r)[s*m], &incOne);
	
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
      
      for(j = 0; j < m; j++){
	REAL(accept_w_str_r)[b*m+j] = accept_w_str[j]/batchLength;
	REAL(tuning_w_str_r)[b*m+j] = w_strTuning[j];
	
	if(accept_w_str[j]/batchLength > acceptRate){
	  w_strTuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}else{
	  w_strTuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}
	accept_w_str[j] = 0.0;
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
    
    int nResultListObjs = 7;
    
    PROTECT(result = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result, 0, samples_r);
    SET_VECTOR_ELT(resultNames, 0, Rf_mkChar("p.beta.theta.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, Rf_mkChar("acceptance"));

    SET_VECTOR_ELT(result, 2, accept_w_str_r);
    SET_VECTOR_ELT(resultNames, 2, Rf_mkChar("acceptance.w.knots"));
    
    SET_VECTOR_ELT(result, 3, w_r);
    SET_VECTOR_ELT(resultNames, 3, Rf_mkChar("p.w.samples"));

    SET_VECTOR_ELT(result, 4, w_str_r);
    SET_VECTOR_ELT(resultNames, 4, Rf_mkChar("p.w.knots.samples"));

    SET_VECTOR_ELT(result, 5, tuning_r);
    SET_VECTOR_ELT(resultNames, 5, Rf_mkChar("tuning"));

    SET_VECTOR_ELT(result, 6, tuning_w_str_r);
    SET_VECTOR_ELT(resultNames, 6, Rf_mkChar("tuning.w.knots"));
  
    Rf_namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);

  }
}
