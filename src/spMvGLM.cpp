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

   SEXP spMvGLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,SEXP family_r, SEXP weights_r,
		SEXP betaPrior_r, SEXP betaNorm_r, SEXP KIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
		SEXP phiStarting_r, SEXP AStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP wStarting_r,
		SEXP phiTuning_r, SEXP ATuning_r, SEXP nuTuning_r , SEXP betaTuning_r, SEXP wTuning_r,
		SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r){


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

    //tuning
    double *betaTuning = (double *) R_alloc(p*p, sizeof(double)); 
    F77_NAME(dcopy)(&pp, REAL(betaTuning_r), &incOne, betaTuning, &incOne);
    double *phiTuning = REAL(phiTuning_r);
    double *ATuning = REAL(ATuning_r);
    double *wTuning = REAL(wTuning_r);
    double *nuTuning = NULL;

    if(covModel == "matern")
      nuTuning = REAL(nuTuning_r);

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

      Rprintf("Metropolis tuning values:\n");

      Rprintf("\tbeta tuning:\n");
      printMtrx(betaTuning, p, p);
      Rprintf("\n"); 
  
      Rprintf("\tA tuning:\n");
      Rprintf("\t"); printVec(ATuning, nLTr);
      Rprintf("\n"); 

      Rprintf("\tphi tuning\n");
      Rprintf("\t"); printVec(phiTuning, m);
      Rprintf("\n");   

      if(covModel == "matern"){
	Rprintf("\tnu tuning\n");
	Rprintf("\t"); printVec(nuTuning, m);
	Rprintf("\n");
      }

      Rprintf("Metropolis starting values:\n");
  
      Rprintf("\tbeta starting:\n");
      Rprintf("\t"); printVec(betaStarting, p);
      Rprintf("\n"); 

      Rprintf("\tA starting:\n");
      Rprintf("\t"); printVec(AStarting, nLTr);
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
    int nn = n*n;
    int mm = m*m;
    int nm = n*m;

    //spatial parameters
    int nParams, betaIndx, AIndx, LIndx, phiIndx, nuIndx;

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
    }

    if(covModel == "matern"){
      for(i = 0; i < m; i++){
	spParams[nuIndx+i] = logit(nuStarting[i], nuUnif[i*2], nuUnif[i*2+1]);
      }
    }

    double *wCurrent = (double *) R_alloc(nm, sizeof(double));
    F77_NAME(dcopy)(&nm, wStarting, &incOne, wCurrent, &incOne);

   
    //samples and random effects
    SEXP w_r, samples_r, accept_r;

    PROTECT(w_r = allocMatrix(REALSXP, nm, nSamples)); nProtect++; 
    double *w = REAL(w_r); zeros(w, nm*nSamples);

    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    double *samples = REAL(samples_r);

    PROTECT(accept_r = allocMatrix(REALSXP, 1, 1)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status=0, rtnStatus=0, accept=0, batchAccept = 0;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0, logDetK, SKtrace;

    double *C = (double *) R_alloc(nm*nm, sizeof(double)); 
        
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm1 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm2 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *tmp_nm1 = (double *) R_alloc(nm, sizeof(double));

    double *candSpParams = (double *) R_alloc(nParams, sizeof(double));
    double *beta = (double *) R_alloc(p, sizeof(double));
    double *A = (double *) R_alloc(mm, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));

    double *wCand = (double *) R_alloc(nm, sizeof(double));
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
    for(s = 0; s < nSamples; s++){
 
      //propose   
      mvrnorm(&candSpParams[betaIndx], &spParams[betaIndx], betaTuning, p, false);
      F77_NAME(dcopy)(&p, &candSpParams[betaIndx], &incOne, beta, &incOne);

      for(i = 0; i < nLTr; i++){
	candSpParams[AIndx+i] = rnorm(spParams[AIndx+i], sqrt(ATuning[i]));
      }

      covTransInvExpand(&candSpParams[AIndx], A, m);

      for(i = 0; i < m; i++){
	candSpParams[phiIndx+i] = rnorm(spParams[phiIndx+i], sqrt(phiTuning[i]));
	phi[i] = logitInv(candSpParams[phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
	if(covModel == "matern"){
	  candSpParams[nuIndx+i] = rnorm(spParams[nuIndx+i], sqrt(nuTuning[i]));
	  nu[i] = logitInv(candSpParams[nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
	}
      }

      for(i = 0; i < nm; i++){
	wCand[i] = rnorm(wCurrent[i], sqrt(wTuning[i]));
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
      for(i = 0; i < nm; i++) detCand += 2*log(C[i*nm+i]);
      F77_NAME(dpotri)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotri failed\n");}
      
      //Likelihood with Jacobian   
      logPostCand = 0.0;

      if(betaPrior == "normal"){
	for(i = 0; i < p; i++){
	  logPostCand += dnorm(beta[i], betaMu[i], betaSd[i], 1);
	}
      }      

      //Jacobian and IW priors for K = A'A
      logDetK = 0.0;
      SKtrace = 0.0;

      for(i = 0; i < m; i++) logDetK += 2*log(A[i*m+i]);
      
      //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
      for(i = 0; i < m; i++) logPostCand += (m-i)*log(A[i*m+i])+log(A[i*m+i]);
      
      //get S*K^-1, already have the chol of K (i.e., A)
      F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
      F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
      for(i = 0; i < m; i++){SKtrace += tmp_mm[i*m+i];}
      logPostCand += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;

      for(i = 0; i < m; i++){
	logPostCand += log(phi[i] - phiUnif[i*2]) + log(phiUnif[i*2+1] - phi[i]); 
      
	if(covModel == "matern"){
	  logPostCand += log(nu[i] - nuUnif[i*2]) + log(nuUnif[i*2+1] - nu[i]);  
	}
      }
   
      F77_NAME(dgemv)(ntran, &nm, &p, &one, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
     
      if(family == "binomial"){
	logPostCand += binomial_logpost(nm, Y, tmp_nm, wCand, weights);
      }else if(family == "poisson"){
	logPostCand += poisson_logpost(nm, Y, tmp_nm, wCand, weights);
      }else{
	error("c++ error: family misspecification in spGLM\n");
      }

      //(-1/2) * tmp_n` *  C^-1 * tmp_n
      F77_NAME(dsymv)(lower, &nm, &one,  C, &nm, wCand, &incOne, &zero, tmp_nm, &incOne);
      logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&nm, wCand, &incOne, tmp_nm, &incOne);

      //
      //MH accept/reject	
      //      
  
      //MH ratio with adjustment
      logMHRatio = logPostCand - logPostCurrent;
      
      if(runif(0.0,1.0) <= exp(logMHRatio)){
	F77_NAME(dcopy)(&nParams, candSpParams, &incOne, spParams, &incOne);
	F77_NAME(dcopy)(&nm, wCand, &incOne, wCurrent, &incOne);
	logPostCurrent = logPostCand;
	accept++;
	batchAccept++;
      }
      
      /******************************
          Save samples and report
      *******************************/
      F77_NAME(dcopy)(&nParams, spParams, &incOne, &samples[s*nParams], &incOne);
      F77_NAME(dcopy)(&nm, wCurrent, &incOne, &w[s*nm], &incOne);
      
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
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    //untransform variance variables
    for(s = 0; s < nSamples; s++){
 
      covTransInv(&samples[s*nParams+AIndx], &samples[s*nParams+AIndx], m);
   	
      for(i = 0; i < m; i++){
	samples[s*nParams+phiIndx+i] = logitInv(samples[s*nParams+phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
	if(covModel == "matern"){
	  samples[s*nParams+nuIndx+i] = logitInv(samples[s*nParams+nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
	}
      }
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
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.beta.theta.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));
    
    SET_VECTOR_ELT(result, 2, w_r);
    SET_VECTOR_ELT(resultNames, 2, mkChar("p.w.samples"));
  
    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    
  }
}
