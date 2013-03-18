#include <algorithm>
#include <string>
// #ifdef _OPENMP
// #include <omp.h>
// #endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"


extern "C" {

  SEXP spMvGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,SEXP family_r, SEXP weights_r,
		     SEXP betaPrior_r, SEXP betaNorm_r, SEXP KIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
		     SEXP phiStarting_r, SEXP AStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP wStarting_r,
		     SEXP phiTuning_r, SEXP ATuning_r, SEXP nuTuning_r , SEXP betaTuning_r, SEXP wTuning_r,
		     SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, b, s, ii, jj, info, nProtect= 0;
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
    int m = INTEGER(m_r)[0];
    int mm = m*m;
    int nn = n*n;
    int nm = n*m;
    int nLTr = m*(m-1)/2+m;

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
    
    double *phiUnif = REAL(phiUnif_r);
    double KIW_df = REAL(VECTOR_ELT(KIW_r, 0))[0]; double *KIW_S = REAL(VECTOR_ELT(KIW_r, 1));

    double *phiStarting = REAL(phiStarting_r);
    double *AStarting = REAL(AStarting_r);
    double *betaStarting = REAL(betaStarting_r);
    double *wStarting = REAL(wStarting_r);

    //if matern
    double *nuUnif = NULL;
    double *nuStarting = NULL;

    if(covModel == "matern"){
      nuUnif = REAL(nuUnif_r);
      nuStarting = REAL(nuStarting_r);
    }

    int nBatch = INTEGER(nBatch_r)[0];
    int batchLength = INTEGER(batchLength_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    int nSamples = nBatch*batchLength;

    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    double *A = (double *) R_alloc(mm, sizeof(double));

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
  
      Rprintf("\tK IW hyperpriors df=%.5f, S=\n", KIW_df);
      printMtrx(KIW_S, m, m);
      Rprintf("\n"); 

      Rprintf("\tphi Unif hyperpriors\n");
      Rprintf("\t");   
      for(i = 0; i < m; i++){
	Rprintf("(%.5f, %.5f) ", phiUnif[i*2], phiUnif[i*2+1]);
      }
      Rprintf("\n\n");   

      if(covModel == "matern"){
	Rprintf("\tnu Unif hyperpriors\n");
	Rprintf("\t");
	for(i = 0; i < m; i++){
	  Rprintf("(%.5f, %.5f) ", nuUnif[i*2], nuUnif[i*2+1]);
	}
	Rprintf("\n\n");   
      }

      Rprintf("Adaptive Metropolis with target acceptance rate: %.1f\n", 100*acceptRate);
      Rprintf("\n");   

      Rprintf("Metropolis starting values:\n");
  
      Rprintf("\tbeta starting:\n");
      Rprintf("\t"); printVec(betaStarting, p);
      Rprintf("\n"); 

      covExpand(AStarting, A, m);
      Rprintf("\tA starting\n");
      printMtrx(A, m, m);
      Rprintf("\n"); 

      Rprintf("\tphi starting\n");
      Rprintf("\t"); printVec(phiStarting, m);
      Rprintf("\n");   

      if(covModel == "matern"){
	Rprintf("\tnu starting\n");
	Rprintf("\t"); printVec(nuStarting, m);
      }

    }
 
    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/
    //spatial parameters
    int nParams, betaIndx, AIndx, phiIndx, nuIndx;

    if(covModel != "matern"){
      nParams = p+nLTr+m;//A, phi
      betaIndx = 0; AIndx = betaIndx+p; phiIndx = AIndx+nLTr;
    }else {
      nParams = p+nLTr+2*m;//A, phi, nu
      betaIndx = 0; AIndx = betaIndx+p; phiIndx = AIndx+nLTr, nuIndx = phiIndx+m;
    }
    
    double *spParams = (double *) R_alloc(nParams, sizeof(double));
    
    //set starting
    F77_NAME(dcopy)(&p, betaStarting, &incOne, &spParams[betaIndx], &incOne);

    covTrans(AStarting, &spParams[AIndx], m);

    for(i = 0; i < m; i++){
      spParams[phiIndx+i] = logit(phiStarting[i], phiUnif[i*2], phiUnif[i*2+1]);
      
      if(covModel == "matern"){
	spParams[nuIndx+i] = logit(nuStarting[i], nuUnif[i*2], nuUnif[i*2+1]);
      }
    }

    double *w = (double *) R_alloc(nm, sizeof(double));
    F77_NAME(dcopy)(&nm, wStarting, &incOne, w, &incOne);

    //set tuning
    double *spTuning = (double *) R_alloc(nParams, sizeof(double));
    
    F77_NAME(dcopy)(&p, REAL(betaTuning_r), &incOne, &spTuning[betaIndx], &incOne);

    for(i = 0; i < nLTr; i++)
      spTuning[AIndx+i] = REAL(ATuning_r)[i];
    
    for(i = 0; i < m; i++){
      spTuning[phiIndx+i] = REAL(phiTuning_r)[i];
      
      if(covModel == "matern"){
	spTuning[nuIndx+i] = REAL(nuTuning_r)[i];
      }
    }
    
    for(i = 0; i < nParams; i++)
      spTuning[i] = log(spTuning[i]);

    double *wTuning = (double *) R_alloc(nm, sizeof(double));
    F77_NAME(dcopy)(&nm, REAL(wTuning_r), &incOne, wTuning, &incOne);

    for(i = 0; i < nm; i++)
      wTuning[i] = log(wTuning[i]);

    //samples and random effects
    SEXP w_r, samples_r, accept_r, accept_w_r, tuning_r, tuning_w_r;

    PROTECT(w_r = allocMatrix(REALSXP, nm, nSamples)); nProtect++;  
    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    PROTECT(accept_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++;
    PROTECT(accept_w_r = allocMatrix(REALSXP, nm, nBatch)); nProtect++;
    PROTECT(tuning_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++;
    PROTECT(tuning_w_r = allocMatrix(REALSXP, nm, nBatch)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status=0;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0, logDetK, SKtrace;
    double spParamsjCurrent, wjCurrent;

    double *accept = (double *) R_alloc(nParams, sizeof(double)); zeros(accept, nParams);
    double *accept_w = (double *) R_alloc(nm, sizeof(double)); zeros(accept_w, nm);

    double *C = (double *) R_alloc(nm*nm, sizeof(double)); 
        
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm1 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm2 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *tmp_nm1 = (double *) R_alloc(nm, sizeof(double));

    double *beta = (double *) R_alloc(p, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));

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
	  
	  covTransInvExpand(&spParams[AIndx], A, m);
	  
	  for(k = 0; k < m; k++){
	    phi[k] = logitInv(spParams[phiIndx+k], phiUnif[k*2], phiUnif[k*2+1]);
	    
	    if(covModel == "matern"){
	      nu[k] = logitInv(spParams[nuIndx+k], nuUnif[k*2], nuUnif[k*2+1]);
	    }
	  }

	  //construct covariance matrix
          // #pragma omp parallel 
	  // {
          // #pragma omp for private(ii, k, l, h)
	  for(jj = 0; jj < n; jj++){
	    for(ii = jj; ii < n; ii++){	
	      for(k = 0; k < m; k++){
		for(l = 0; l < m; l++){
		  C[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
		  for(h = 0; h < m; h++){
		    C[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(coordsD[jj*n+ii], phi[h], nu[h], covModel);
		  }
		}
	      }
	    }
	  }
	  // } //parallel for

	  //invert C and log det cov
	  detCand = 0;
	  F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: Cholesky failed in spGLM\n");}
	  for(k = 0; k < nm; k++) detCand += 2*log(C[k*nm+k]);
	  F77_NAME(dpotri)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spGLM\n");}

	  //Likelihood with Jacobian   
	  logPostCand = 0.0;
	  
	  if(betaPrior == "normal"){
	    for(k = 0; k < p; k++){
	      logPostCand += dnorm(beta[k], betaMu[k], betaSd[k], 1);
	    }
	  }      
	  
	  //Jacobian and IW priors for K = A'A
	  logDetK = 0.0;
	  SKtrace = 0.0;
	  
	  for(k = 0; k < m; k++) logDetK += 2*log(A[k*m+k]);
	  
	  //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	  for(k = 0; k < m; k++) logPostCand += (m-k)*log(A[k*m+k])+log(A[k*m+k]);
	  
	  //get S*K^-1, already have the chol of K (i.e., A)
	  F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	  F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
	  for(k = 0; k < m; k++){SKtrace += tmp_mm[k*m+k];}
	  logPostCand += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;
	  
	  for(k = 0; k < m; k++){
	    logPostCand += log(phi[k] - phiUnif[k*2]) + log(phiUnif[k*2+1] - phi[k]); 
	    
	    if(covModel == "matern"){
	      logPostCand += log(nu[k] - nuUnif[k*2]) + log(nuUnif[k*2+1] - nu[k]);  
	    }
	  }
	  
	  F77_NAME(dgemv)(ntran, &nm, &p, &one, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
	  
	  if(family == "binomial"){
	    logPostCand += binomial_logpost(nm, Y, tmp_nm, w, weights);
	  }else if(family == "poisson"){
	    logPostCand += poisson_logpost(nm, Y, tmp_nm, w, weights);
	  }else{
	    error("c++ error: family misspecification in spGLM\n");
	  }
	  
	  //(-1/2) * tmp_n` *  C^-1 * tmp_n
	  F77_NAME(dsymv)(lower, &nm, &one, C, &nm, w, &incOne, &zero, tmp_nm1, &incOne);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&nm, w, &incOne, tmp_nm1, &incOne);
	  
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
	
	covTransInvExpand(&spParams[AIndx], A, m);
	
	for(k = 0; k < m; k++){
	  phi[k] = logitInv(spParams[phiIndx+k], phiUnif[k*2], phiUnif[k*2+1]);
	  
	  if(covModel == "matern"){
	    nu[k] = logitInv(spParams[nuIndx+k], nuUnif[k*2], nuUnif[k*2+1]);
	  }
	}
	
	//construct covariance matrix
        // #pragma omp parallel 
	// {
        // #pragma omp for private(ii, k, l, h)
	  for(jj = 0; jj < n; jj++){
	    for(ii = jj; ii < n; ii++){	
	      for(k = 0; k < m; k++){
		for(l = 0; l < m; l++){
		  C[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
		  for(h = 0; h < m; h++){
		    C[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(coordsD[jj*n+ii], phi[h], nu[h], covModel);
		  }
		}
	      }
	    }
	  }
	// } //parallel for
	
	//invert C and log det cov
	detCand = 0;
	F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	for(k = 0; k < nm; k++) detCand += 2*log(C[k*nm+k]);
	F77_NAME(dpotri)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotri failed\n");}

	F77_NAME(dgemv)(ntran, &nm, &p, &one, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);

	for(j = 0; j < nm; j++){
	  
	  //propose
	  wjCurrent = w[j];
	  w[j] = rnorm(wjCurrent, exp(wTuning[j]));
	  
	  //
	  //Likelihood with Jacobian   
	  //
	  logPostCand = 0.0;
	  
	  if(betaPrior == "normal"){
	    for(k = 0; k < p; k++){
	      logPostCand += dnorm(beta[k], betaMu[k], betaSd[k], 1);
	    }
	  }      
	  
	  //Jacobian and IW priors for K = A'A
	  covTransInvExpand(&spParams[AIndx], A, m);	  

	  logDetK = 0.0;
	  SKtrace = 0.0;
	  
	  for(k = 0; k < m; k++) logDetK += 2*log(A[k*m+k]);
	  
	  //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	  for(k = 0; k < m; k++) logPostCand += (m-k)*log(A[k*m+k])+log(A[k*m+k]);
	  
	  //get S*K^-1, already have the chol of K (i.e., A)
	  F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	  F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
	  for(k = 0; k < m; k++){SKtrace += tmp_mm[k*m+k];}
	  logPostCand += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;
	  
	  for(k = 0; k < m; k++){
	    logPostCand += log(phi[k] - phiUnif[k*2]) + log(phiUnif[k*2+1] - phi[k]); 
	    
	    if(covModel == "matern"){
	      logPostCand += log(nu[k] - nuUnif[k*2]) + log(nuUnif[k*2+1] - nu[k]);  
	    }
	  }
	  
	  if(family == "binomial"){
	    logPostCand += binomial_logpost(nm, Y, tmp_nm, w, weights);
	  }else if(family == "poisson"){
	    logPostCand += poisson_logpost(nm, Y, tmp_nm, w, weights);
	  }else{
	    error("c++ error: family misspecification in spGLM\n");
	  }
	  
	  //(-1/2) * tmp_n` *  C^-1 * tmp_n
	  F77_NAME(dsymv)(lower, &nm, &one, C, &nm, w, &incOne, &zero, tmp_nm1, &incOne);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&nm, w, &incOne, tmp_nm1, &incOne);
	  
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
          Save samples and report
	*******************************/
	F77_NAME(dcopy)(&nParams, spParams, &incOne, &REAL(samples_r)[s*nParams], &incOne);
	F77_NAME(dcopy)(&nm, w, &incOne, &REAL(w_r)[s*nm], &incOne);
	
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
      
      for(j = 0; j < nm; j++){
	REAL(accept_w_r)[b*nm+j] = accept_w[j]/batchLength;
	REAL(tuning_w_r)[b*nm+j] = wTuning[j];
	
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
	  
	  covTransInvExpand(&spParams[AIndx], A, m);
	  
	  for(k = 0, i = 0; k < m; k++){
	    for(j = k; j < m; j++, i++){
	      Rprintf("\tA[%i,%i]\t\t%3.1f\t\t%1.5f\n",j, k, 100.0*REAL(accept_r)[b*nParams+AIndx+j], exp(spTuning[AIndx+i]));
	    }
	  }
	  
	  for(j = 0; j < m; j++){
	    Rprintf("\tphi[%i]\t\t%3.1f\t\t%1.5f\n", j, 100.0*REAL(accept_r)[b*nParams+phiIndx+j], exp(spTuning[phiIndx+j]));	
	    if(covModel == "matern")
	      Rprintf("\tnu[%i]\t\t%3.1f\t\t%1.5f\n", j, 100.0*REAL(accept_r)[b*nParams+nuIndx+j], exp(spTuning[nuIndx+j]));
	  } 
	  
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
      
      covTransInv(&REAL(samples_r)[s*nParams+AIndx], &REAL(samples_r)[s*nParams+AIndx], m);
      
      for(i = 0; i < m; i++){
    	REAL(samples_r)[s*nParams+phiIndx+i] = logitInv(REAL(samples_r)[s*nParams+phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
    	if(covModel == "matern"){
    	  REAL(samples_r)[s*nParams+nuIndx+i] = logitInv(REAL(samples_r)[s*nParams+nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
    	}
      }
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
