#include <algorithm>
#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP spLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
	    SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	    SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r,
	    SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP nuTuning_r, 
	    SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r){

    /*****************************************
                Common variables
    *****************************************/
    int i, j, k, l, b, s, info, nProtect=0;
    char const *lower = "L";
    char const *upper = "U";
    char const *nUnit = "N";
    char const *yUnit = "U";
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
    int nn = n*n;
    int np = n*p;
    
    double *coordsD = REAL(coordsD_r);

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
      
      if(amcmc){
	Rprintf("Using adaptive MCMC.\n\n");
	Rprintf("\tNumber of batches %i.\n", nBatch);
	Rprintf("\tBatch length %i.\n", batchLength);
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
    int nParams, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;

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
    SEXP samples_r, accept_r, tuning_r;
    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 

    if(amcmc){
      PROTECT(accept_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++; 
      PROTECT(tuning_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++; 
    }else{
      PROTECT(accept_r = allocMatrix(REALSXP, 1, floor(static_cast<double>(nSamples/nReport)))); nProtect++; 
    }

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status=1, batchAccept=0, reportCnt=0;
    double logMHRatio =0, logPostCurrent = R_NegInf, logPostCand = 0, det = 0, paramsjCurrent = 0;
   
    double *paramsCurrent = (double *) R_alloc(nParams, sizeof(double));
    double *accept = (double *) R_alloc(nParams, sizeof(double)); zeros(accept, nParams);

    double *C = (double *) R_alloc(nn, sizeof(double)); zeros(C, nn);
    double sigmaSq, phi, tauSq, nu, Q;
    double *theta = (double *) R_alloc(3, sizeof(double)); //phi, nu, and perhaps more in the future

    int p1 = p+1;
    double *vU = (double *) R_alloc(n*p1, sizeof(double));

    double *z = (double *) R_alloc(n, sizeof(double)); 
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_nn = NULL;
    double *Cbeta = NULL;

    if(betaPrior == "normal"){
      tmp_nn = (double *) R_alloc(nn, sizeof(double));
      Cbeta = (double *) R_alloc(nn, sizeof(double));

      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, betaMu, &incOne, &zero, z, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, z, &incOne);

      F77_NAME(dsymm)(rside, lower, &n, &p, &one, betaC, &p, X, &n, &zero, vU, &n);
      F77_NAME(dgemm)(ntran, ytran, &n, &n, &p, &one, vU, &n, X, &n, &zero, tmp_nn, &n);
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
	  
	  //construct covariance matrix
	  spCovLT(coordsD, n, theta, covModel, C);
	  
	  if(nugget){
	    for(k = 0; k < n; k++){
	      C[k*n+k] += tauSq;
	    }
	  }
	  
	  if(betaPrior == "normal"){

	    for(k = 0; k < n; k++){
	      for(l = k; l < n; l++){
	    	Cbeta[k*n+l] = C[k*n+l]+tmp_nn[k*n+l];
	      }
	    }

	    det = 0;
	    F77_NAME(dpotrf)(lower, &n, Cbeta, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	    for(k = 0; k < n; k++) det += 2*log(Cbeta[k*n+k]);
	    
	    F77_NAME(dcopy)(&n, z, &incOne, tmp_n, &incOne);
	    F77_NAME(dtrsv)(lower, ntran, nUnit, &n, Cbeta, &n, tmp_n, &incOne);//L[u] = (y-X'beta)

	    Q = pow(F77_NAME(dnrm2)(&n, tmp_n, &incOne),2);
	  }else{//beta flat
	    det = 0;
	    F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	    for(k = 0; k < n; k++) det += 2*log(C[k*n+k]);
	    
	    F77_NAME(dcopy)(&n, Y, &incOne, vU, &incOne);
	    F77_NAME(dcopy)(&np, X, &incOne, &vU[n], &incOne);
	    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p1, &one, C, &n, vU, &n);//L[v:U] = [y:X]

	    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, &vU[n], &n, &vU[n], &n, &zero, tmp_pp, &p); //U'U
	    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	    for(k = 0; k < p; k++) det += 2*log(tmp_pp[k*p+k]);
	    
	    F77_NAME(dgemv)(ytran, &n, &p, &one, &vU[n], &n, vU, &incOne, &zero, tmp_p, &incOne); //U'v
	    F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p, &incOne);

	    Q = pow(F77_NAME(dnrm2)(&n, vU, &incOne),2) - pow(F77_NAME(dnrm2)(&p, tmp_p, &incOne),2) ;
	  }
	  
	  //
	  //priors, jacobian adjustments, and likelihood
	  //
 	  logPostCand = -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq);
	  
	  if(nugget){
	    logPostCand += -1.0*(1.0+tauSqIGa)*log(tauSq)-tauSqIGb/tauSq+log(tauSq);
	  }
	  
	  logPostCand += log(phi - phiUnifa) + log(phiUnifb - phi); 
	  
	  if(covModel == "matern"){
	    logPostCand += log(nu - nuUnifa) + log(nuUnifb - nu);   
	  }
	  
	  logPostCand += -0.5*det-0.5*Q;
	  
	  //
	  //MH accept/reject	
	  //      
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
	      F77_NAME(dcopy)(&nParams, paramsCurrent, &incOne, params, &incOne);
	    }
	  }
	  
	  if(!amcmc){
	    break;
	  }
	}//end params
	
	/******************************
               Save samples
	*******************************/
	F77_NAME(dcopy)(&nParams, params, &incOne, &REAL(samples_r)[s*nParams], &incOne);
	
	R_CheckUserInterrupt();
      }//end batch
      
      //adjust tuning
      if(amcmc){
      	for(j = 0; j < nParams; j++){
      	  REAL(accept_r)[b*nParams+j] = accept[j]/batchLength;
      	  REAL(tuning_r)[b*nParams+j] = exp(tuning[j]);
	  
      	  if(accept[j]/batchLength > acceptRate){
      	    tuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
      	  }else{
      	    tuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
      	  }
      	  accept[j] = 0.0;
      	}
      }
      
      //report
      if(status == nReport){

	if(verbose){
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
      	}

	if(!amcmc){
	  REAL(accept_r)[reportCnt] = 100.0*batchAccept/nReport;
	  reportCnt++;
	}

	batchAccept = 0;
	status = 0;

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
    int nResultListObjs = 2;

    if(amcmc){
      nResultListObjs++;
    }
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    //samples
    SET_VECTOR_ELT(result_r, 0, samples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.theta.samples")); 

    SET_VECTOR_ELT(result_r, 1, accept_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("acceptance"));

    if(amcmc){
      SET_VECTOR_ELT(result_r, 2, tuning_r);
      SET_VECTOR_ELT(resultName_r, 2, mkChar("tuning"));
    }

    namesgets(result_r, resultName_r);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
