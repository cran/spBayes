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

  SEXP spPPLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP knotsD_r, SEXP knotsCoordsD_r, 
	      SEXP modPP_r, SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	      SEXP betaStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r,
	      SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP nuTuning_r, 
	      SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r){

    /*****************************************
                Common variables
    *****************************************/
    int i, j, k, l, b, s, info, nProtect=0;
    char const *lower = "L";
    char const *nUnit = "N";
    char const *ntran = "N";
    char const *ytran = "T";
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
    int np = n*p;
    int m = INTEGER(m_r)[0];
    int nm = n*m;
    int mm = m*m;
    int mp = m*p;

    double *knotsD = REAL(knotsD_r);
    double *knotsCoordsD = REAL(knotsCoordsD_r);

    bool modPP = static_cast<bool>(INTEGER(modPP_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));

    //priors
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));
    double *betaMu = NULL;
    double *betaC = NULL;
    
    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(p, sizeof(double));
      F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);

      betaC = (double *) R_alloc(pp, sizeof(double)); 
      F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaC, &incOne);
    }

    double sigmaSqIGa = REAL(sigmaSqIG_r)[0]; double sigmaSqIGb = REAL(sigmaSqIG_r)[1];
    double phiUnifa = REAL(phiUnif_r)[0]; double phiUnifb = REAL(phiUnif_r)[1];

    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    double tauSqIGa = 0, tauSqIGb = 0;
    if(nugget){
      tauSqIGa = REAL(tauSqIG_r)[0]; tauSqIGb = REAL(tauSqIG_r)[1]; 
    }

    //matern
    double nuUnifa = 0, nuUnifb = 0;
    if(covModel == "matern"){
      nuUnifa = REAL(nuUnif_r)[0]; nuUnifb = REAL(nuUnif_r)[1]; 
    }

    bool amcmc = static_cast<bool>(INTEGER(amcmc_r)[0]);
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
      
      if(modPP){
	Rprintf("Using modified predictive process with %i knots.\n\n", m);
      }else{
	Rprintf("Using non-modified predictive process with %i knots.\n\n", m);
      }

      if(amcmc){
	Rprintf("Using adaptive MCMC.\n\n");
	Rprintf("\tNumber of batches %i.\n", nBatch);
	Rprintf("\tBatch Rf_length %i.\n", batchLength);
	Rprintf("\tTarget acceptance rate %.5f.\n", acceptRate);
	Rprintf("\n");
      }else{
	Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      }

      if(!nugget){
	Rprintf("tau.sq not included in the model (i.e., no nugget model).\n\n");
      }

      Rprintf("Priors and hyperpriors:\n");
      
      if(betaPrior == "flat"){
	Rprintf("\tbeta flat.\n");
      }else{
	Rprintf("\tbeta normal:\n");
	Rprintf("\tmu:"); printVec(betaMu, p);
	Rprintf("\tcov:\n"); printMtrx(betaC, p, p);
	Rprintf("\n");
      }
 
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
         Set-up MCMC sample matrices etc.
    *****************************************/ 
    //parameters
    int nParams, sigmaSqIndx, tauSqIndx = 0, phiIndx, nuIndx = 0;

    if(!nugget && covModel != "matern"){
      nParams = 2;//sigma^2, phi
      sigmaSqIndx = 0; phiIndx = 1;
    }else if(nugget && covModel != "matern"){
      nParams = 3;//sigma^2, tau^2, phi
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
    }else if(!nugget && covModel == "matern"){
      nParams = 3;//sigma^2, phi, nu
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
    }else{
      nParams = 4;//sigma^2, tau^2, phi, nu
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;//sigma^2, tau^2, phi, nu
    }
    
    double *params = (double *) R_alloc(nParams, sizeof(double));

    //starting
    double *beta = (double *) R_alloc(p, sizeof(double)); 
    F77_NAME(dcopy)(&p, REAL(betaStarting_r), &incOne, beta, &incOne);

    params[sigmaSqIndx] = log(REAL(sigmaSqStarting_r)[0]);

    if(nugget){
      params[tauSqIndx] = log(REAL(tauSqStarting_r)[0]);
    }

    params[phiIndx] = logit(REAL(phiStarting_r)[0], phiUnifa, phiUnifb);

    if(covModel == "matern"){
      params[nuIndx] = logit(REAL(nuStarting_r)[0], nuUnifa, nuUnifb);
    }

    //tuning and fixed
    double *tuning = (double *) R_alloc(nParams, sizeof(double));
    int *fixed = (int *) R_alloc(nParams, sizeof(int)); zeros(fixed, nParams);
    
    tuning[sigmaSqIndx] = REAL(sigmaSqTuning_r)[0];
    if(tuning[sigmaSqIndx] == 0){
      fixed[sigmaSqIndx] = 1;
    }
          
    if(nugget){
      tuning[tauSqIndx] = REAL(tauSqTuning_r)[0];
      if(tuning[tauSqIndx] == 0){
	fixed[tauSqIndx] = 1;
      }
    }
    
    tuning[phiIndx] = REAL(phiTuning_r)[0];
    if(tuning[phiIndx] == 0){
      fixed[phiIndx] = 1;
    }
    
    if(covModel == "matern"){
      tuning[nuIndx] = REAL(nuTuning_r)[0];
      if(tuning[nuIndx] == 0){
	fixed[nuIndx] = 1;
      }
    }

    for(i = 0; i < nParams; i++){
      tuning[i] = log(sqrt(tuning[i]));
    }

    //return stuff  
    SEXP samples_r, accept_r, tuning_r, betaSamples_r;
    PROTECT(samples_r = Rf_allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    PROTECT(accept_r = Rf_allocMatrix(REALSXP, nParams, nBatch)); nProtect++; 
    PROTECT(tuning_r = Rf_allocMatrix(REALSXP, nParams, nBatch)); nProtect++;
    PROTECT(betaSamples_r = Rf_allocMatrix(REALSXP, p, nSamples)); nProtect++; 
    
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status=1, batchAccept=0;
    double logMHRatio =0, logPostCurrent = R_NegInf, logPostCand = 0, paramsjCurrent = 0;
    double det = 0, detCurrent = 0;
    double priors = 0, priorsCurrent = 0;
    bool updateBeta = false;
   
    double *paramsCurrent = (double *) R_alloc(nParams, sizeof(double));
    double *accept = (double *) R_alloc(nParams, sizeof(double)); zeros(accept, nParams);

    double *P = (double *) R_alloc(nm, sizeof(double)); 
    double *K = (double *) R_alloc(mm, sizeof(double)); 
    double *D = (double *) R_alloc(n, sizeof(double)); 
    double *H = (double *) R_alloc(nm, sizeof(double)); 

    double *u = (double *) R_alloc(n, sizeof(double)); 
    F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, u, &incOne FCONE);
    F77_NAME(daxpy)(&n, &one, Y, &incOne, u, &incOne);

    double *DCurrent = (double *) R_alloc(n, sizeof(double)); 
    double *HCurrent = (double *) R_alloc(nm, sizeof(double)); 
   
    double sigmaSq, phi, tauSq = 0, nu = 0, Q;
    double *theta = (double *) R_alloc(3, sizeof(double)); //phi, nu, and perhaps more in the future

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); 
    double *tmp_np = (double *) R_alloc(np, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mp = (double *) R_alloc(mp, sizeof(double));
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_pp2 = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    
    double *betaCInv = NULL;
    double *betaCInvMu = NULL;
    
    if(betaPrior == "normal"){
      betaCInv = (double *) R_alloc(pp, sizeof(double));
      betaCInvMu = (double *) R_alloc(p, sizeof(double));
      
      F77_NAME(dcopy)(&pp, betaC, &incOne, betaCInv, &incOne);
      F77_NAME(dpotrf)(lower, &p, betaCInv, &p, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &p, betaCInv, &p, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotri failed\n");}
      
      F77_NAME(dsymv)(lower, &p, &one, betaCInv, &p, betaMu, &incOne, &zero, betaCInvMu, &incOne FCONE);      
    }

    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }
    
    GetRNGstate();
    
    for(b = 0, s = 0; b < nBatch; b++){    
      for(i = 0; i < batchLength; i++, s++){
	for(j = 0; j < nParams; j++){
	  
	  //propose
	  if(amcmc){
	    if(fixed[j] == 1){
	      paramsjCurrent = params[j];
	    }else{
	      paramsjCurrent = params[j];
	      params[j] = rnorm(paramsjCurrent, exp(tuning[j]));
	    }
	  }else{
	    F77_NAME(dcopy)(&nParams, params, &incOne, paramsCurrent, &incOne);
	    
	    for(j = 0; j < nParams; j++){
	      if(fixed[j] == 1){
		params[j] = params[j];
	      }else{
		params[j] = rnorm(params[j], exp(tuning[j]));
	      }
	    }
	  }
	  
	  //extract and transform
	  sigmaSq = theta[0] = exp(params[sigmaSqIndx]);
	  phi = theta[1] = logitInv(params[phiIndx], phiUnifa, phiUnifb);
	  
	  if(nugget){
	    tauSq = exp(params[tauSqIndx]);
	  }

	  if(covModel == "matern"){
	    nu = theta[2] = logitInv(params[nuIndx], nuUnifa, nuUnifb);
	  }
	  
	  //construct covariance matrices 
	  spCovLT(knotsD, m, theta, covModel, K);
	  spCov(knotsCoordsD, nm, theta, covModel, P);
	  
	  //get D
	  if(modPP){
	    
	    for(k = 0; k < m; k++){
	      for(l = k; l < m; l++){
		tmp_mm[k*m+l] = K[k*m+l];
	      }
	    }
	    
	    F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed\n");}//L_K
	    F77_NAME(dcopy)(&nm, P, &incOne, H, &incOne);
	    
	    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &m, &n, &one, tmp_mm, &m, H, &m FCONE FCONE FCONE FCONE);
	    
	    for(k = 0; k < n; k++){
	      D[k] = sigmaSq - F77_NAME(ddot)(&m, &H[k*m], &incOne, &H[k*m], &incOne);
	      if(nugget){
		D[k] += tauSq;
	      }
	    }
	    
	  }else{
	    for(k = 0; k < n; k++){
	      D[k] = tauSq;
	    }
	  }	 
	  
	  for(k = 0; k < n; k++){
	    for(l = 0; l < m; l++){
	      H[k*m+l] = P[k*m+l]/sqrt(D[k]);//W'
	    }
	  }
	  
	  F77_NAME(dgemm)(ntran, ytran, &m, &m, &n, &one, H, &m, H, &m, &zero, tmp_mm, &m FCONE FCONE);//W'W
	  
	  for(k = 0; k < m; k++){
	    for(l = k; l < m; l++){
	      K[k*m+l] += tmp_mm[k*m+l];
	    }
	  }
	  
	  F77_NAME(dpotrf)(lower, &m, K, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed\n");}//L
	  
	  F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &m, &n, &one, K, &m, H, &m FCONE FCONE FCONE FCONE);//LH = W'
	  
	  for(k = 0; k < n; k++){
	    tmp_n[k] = u[k]/sqrt(D[k]);//v
	  }
	  
	  F77_NAME(dgemv)(ntran, &m, &n, &one, H, &m, tmp_n, &incOne, &zero, tmp_m, &incOne FCONE); //w = Hv
	  
	  Q = F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n, &incOne) - F77_NAME(ddot)(&m, tmp_m, &incOne, tmp_m, &incOne);
	  
	  F77_NAME(dgemm)(ntran, ytran, &m, &m, &n, &negOne, H, &m, H, &m, &zero, tmp_mm, &m FCONE FCONE);//-HH'
	  
	  for(k = 0; k < m; k++){
	    tmp_mm[k*m+k] += 1.0; //J
	  }
	  
	  F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed\n");}//L_J
	  
	  det = 0;
	  for(k = 0; k < n; k++){
	    det += log(D[k]);
	  }
	  
	  for(k = 0; k < m; k++){
	    det -= 2*log(tmp_mm[k*m+k]);
	  }
	  
	  //Priors, jacobian adjustments, and likelihood
	  priors = -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq);
	  
	  if(nugget){
	    priors += -1.0*(1.0+tauSqIGa)*log(tauSq)-tauSqIGb/tauSq+log(tauSq);
	  }

	  priors += log(phi - phiUnifa) + log(phiUnifb - phi); 
	  
	  if(covModel == "matern"){
	    priors += log(nu - nuUnifa) + log(nuUnifb - nu);   
	  }
	  
	  logPostCand = priors-0.5*det-0.5*Q;
	  
	  //
	  //MH accept/reject	
	  //      
	  logMHRatio = logPostCand - logPostCurrent;

	  if(runif(0.0,1.0) <= exp(logMHRatio)){
	    logPostCurrent = logPostCand;

	    F77_NAME(dcopy)(&n, D, &incOne, DCurrent, &incOne);
	    F77_NAME(dcopy)(&nm, H, &incOne, HCurrent, &incOne);
	    detCurrent = det;
	    priorsCurrent = priors;

	    //set to true so beta's mu and var are updated
	    updateBeta = true;

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
	      F77_NAME(dcopy)(&nParams, paramsCurrent, &incOne, params, &incOne);
	    }
	  }
	  
	  if(!amcmc){
	    break;
	  }
	}//end params
	
	//only update beta's mu and var if theta has changed
	if(updateBeta){
	  
	  for(k = 0; k < n; k++){
	    tmp_n[k] = Y[k]/sqrt(DCurrent[k]);//\tilda{y}
	    for(l = 0; l < p; l++){
	      tmp_np[l*n+k] = X[l*n+k]/sqrt(DCurrent[k]);//V
	    }
	  }
	  
	  F77_NAME(dgemm)(ntran, ntran, &m, &p, &n, &one, HCurrent, &m, tmp_np, &n, &zero, tmp_mp, &m FCONE FCONE);//\tilda{V}
	  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, tmp_np, &n, tmp_np, &n, &zero, tmp_pp, &p FCONE FCONE);//V'V
	  F77_NAME(dgemm)(ytran, ntran, &p, &p, &m, &one, tmp_mp, &m, tmp_mp, &m, &zero, tmp_pp2, &p FCONE FCONE);//\tilda{V}'\tilda{V}
	  
	  for(k = 0; k < p; k++){
	    for(l = k; l < p; l++){
	      tmp_pp[k*p+l] -= tmp_pp2[k*p+l];//Z
	      
	      if(betaPrior == "normal"){
		tmp_pp[k*p+l] += betaCInv[k*p+l];
	      }
	      
	    }
	  }
 	  
	  //B
	  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: Cholesky failed\n");}
	  F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: Cholesky inverse failed\n");}
   
	  //b
	  F77_NAME(dgemv)(ytran, &n, &p, &one, tmp_np, &n, tmp_n, &incOne, &zero, tmp_p, &incOne FCONE); //V'\tilda{y}
	  F77_NAME(dgemv)(ntran, &m, &n, &one, HCurrent, &m, tmp_n, &incOne, &zero, tmp_m, &incOne FCONE); //H\tilda{y}
	  F77_NAME(dgemv)(ytran, &m, &p, &one, tmp_mp, &m, tmp_m, &incOne, &zero, tmp_pp2, &incOne FCONE); //\tilda{V}'(H\tilda{y})
	  
	  for(k = 0; k < p; k++){
	    tmp_p[k] -= tmp_pp2[k]; 
	    
	    if(betaPrior == "normal"){
	      tmp_p[k] += betaCInvMu[k];
	    }
	  }
	  
	  F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &incOne, &zero, tmp_p2, &incOne FCONE); //Bb
	  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed\n");}

	}//end updateBeta

	//set to false so beta's mu and var are only updated when theta has changed
       	updateBeta = false;
	
	//draw beta
	mvrnorm(beta, tmp_p2, tmp_pp, p, false);//note, tmp_p2 and tmp_pp will carry over when theta is not updated
	
	//update logPostCurrent (beta changes on every iteration so logPostCurrent must be updated)
	F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, u, &incOne FCONE);
	F77_NAME(daxpy)(&n, &one, Y, &incOne, u, &incOne);//u
	  
	for(k = 0; k < n; k++){
	  tmp_n[k] = u[k]/sqrt(DCurrent[k]);//v
	}
	  
	F77_NAME(dgemv)(ntran, &m, &n, &one, HCurrent, &m, tmp_n, &incOne, &zero, tmp_m, &incOne FCONE); //w = Hv
	  
	Q = F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n, &incOne) - F77_NAME(ddot)(&m, tmp_m, &incOne, tmp_m, &incOne);

	logPostCurrent = priorsCurrent-0.5*detCurrent-0.5*Q;

	/******************************
               Save samples
	*******************************/
	F77_NAME(dcopy)(&p, beta, &incOne, &REAL(betaSamples_r)[s*p], &incOne);
	F77_NAME(dcopy)(&nParams, params, &incOne, &REAL(samples_r)[s*nParams], &incOne);

	R_CheckUserInterrupt();
      }//end batch
      
      //adjust tuning
      if(amcmc){
	for(j = 0; j < nParams; j++){
	  REAL(accept_r)[b*nParams+j] = accept[j]/batchLength;
	  REAL(tuning_r)[b*nParams+j] = tuning[j];
	  
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
	    Rprintf("Batch: %i of %i, %3.2f%%\n", b+1, nBatch, 100.0*(b+1)/nBatch);
	    Rprintf("\tparameter\tacceptance\ttuning\n");	  
	    Rprintf("\tsigma.sq\t%3.1f\t\t%1.5f\n", 100.0*REAL(accept_r)[b*nParams+sigmaSqIndx], exp(tuning[sigmaSqIndx]));
	    if(nugget){
	      Rprintf("\ttau.sq\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(accept_r)[b*nParams+tauSqIndx], exp(tuning[tauSqIndx]));
	    }
	    Rprintf("\tphi\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(accept_r)[b*nParams+phiIndx], exp(tuning[phiIndx]));
	    if(covModel == "matern"){
	      Rprintf("\tnu\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(accept_r)[b*nParams+nuIndx], exp(tuning[nuIndx]));
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

    //untransform variance variables
    for(s = 0; s < nSamples; s++){
      REAL(samples_r)[s*nParams+sigmaSqIndx] = exp(REAL(samples_r)[s*nParams+sigmaSqIndx]);
      if(nugget){
	REAL(samples_r)[s*nParams+tauSqIndx] = exp(REAL(samples_r)[s*nParams+tauSqIndx]);
      }
      REAL(samples_r)[s*nParams+phiIndx] = logitInv(REAL(samples_r)[s*nParams+phiIndx], phiUnifa, phiUnifb);
      
      if(covModel == "matern")
	REAL(samples_r)[s*nParams+nuIndx] = logitInv(REAL(samples_r)[s*nParams+nuIndx], nuUnifa, nuUnifb);
    }
    
    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 4;

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    //samples
    SET_VECTOR_ELT(result_r, 0, samples_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("p.theta.samples")); 

    SET_VECTOR_ELT(result_r, 1, accept_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("acceptance"));

    SET_VECTOR_ELT(result_r, 2, tuning_r);
    SET_VECTOR_ELT(resultName_r, 2, Rf_mkChar("tuning"));
    
    SET_VECTOR_ELT(result_r, 3, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 3, Rf_mkChar("p.beta.samples"));
  
    Rf_namesgets(result_r, resultName_r);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
