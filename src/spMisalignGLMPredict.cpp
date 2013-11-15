#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  
  SEXP spMisalignGLMPredict(SEXP family_r, SEXP Y_r, SEXP X_r, SEXP m_r, SEXP n_r, SEXP p_r, 
			    SEXP Z_r, SEXP nPred_r, SEXP pPred_r, 
			    SEXP obsD_r, SEXP predObsD_r,
			    SEXP samples_r, SEXP w_r, SEXP nSamples_r, 
			    SEXP covModel_r, 
			    SEXP verbose_r, SEXP nReport_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, s, ii, jj, kk, info, nProtect=0;
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
    std::string family = CHAR(STRING_ELT(family_r,0));
    
    double *Y = REAL(Y_r);
    double *X = REAL(X_r);
    int *p = INTEGER(p_r);//number of X columns
    int *n = INTEGER(n_r);//number of observations
    int m = INTEGER(m_r)[0];//number of outcomes
    int nLTr = m*(m-1)/2+m;

    int N = 0;
    int P = 0;
    for(i = 0; i < m; i++){
      N += n[i];
      P += p[i];
    }
  
    int mm = m*m;
    int NN = N*N;
    int NP = N*P;
    int PP = P*P;
    
    int *pPred = INTEGER(pPred_r);
    int *nPred = INTEGER(nPred_r);
    double *Z = REAL(Z_r);
    
    int NPred = 0;
    int PPred = 0;
    for(i = 0; i < m; i++){
      NPred += nPred[i];
      PPred += pPred[i];
    }

    int NPredN = NPred*N;
 
    double *obsD = REAL(obsD_r);
    double *predObsD = REAL(predObsD_r);
    
    double *samples = REAL(samples_r);
    double *w = REAL(w_r);
    int nSamples = INTEGER(nSamples_r)[0];
    
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    
    int nParams, betaIndx, AIndx, phiIndx, nuIndx;
  
    if(covModel != "matern"){
      nParams = P+nLTr+m;//A, phi
      betaIndx = 0; AIndx = betaIndx+P; phiIndx = AIndx+nLTr;
    }else {
      nParams = P+nLTr+2*m;//A, phi, nu
      betaIndx = 0; AIndx = betaIndx+P; phiIndx = AIndx+nLTr, nuIndx = phiIndx+m;
    }
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      if(family == "binomial"){
	Rprintf("Binomial regression.\n\n");
      }else if(family == "poisson"){
	Rprintf("Poisson regression.\n\n");
      }

      Rprintf("Model fit with %i outcome variables.\n\n", m);
      Rprintf("Number of observations within each outcome:"); printVec(n, m);
      Rprintf("\nNumber of covariates for each outcome (including intercept if specified):"); printVec(p, m);
      Rprintf("\nTotal number of observations: %i\n\n", N);
      Rprintf("Total number of covariates (including intercept if specified): %i\n\n", P);

      Rprintf("Number of prediction for each outcome:"); printVec(nPred, m);   
      Rprintf("\nTotal number of predictions: %i\n\n", NPred);

      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
    
    }
    
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    SEXP predSamples_r, predWSamples_r;
    PROTECT(predSamples_r = allocMatrix(REALSXP, NPred, nSamples)); nProtect++;  
    PROTECT(predWSamples_r = allocMatrix(REALSXP, NPred, nSamples)); nProtect++;  
 
    int status=1;
    
    double *C = (double *) R_alloc(NN, sizeof(double));
    double *c = (double *) R_alloc(NPredN, sizeof(double));
    double *z = (double *) R_alloc(NPred, sizeof(double)); new double[NPred];
    double *u = (double *) R_alloc(N, sizeof(double)); new double[N];
    double *v = (double *) R_alloc(N, sizeof(double)); new double[N];
    double muPred, varPred;
    double *beta = (double *) R_alloc(P, sizeof(double));
    double *A = (double *) R_alloc(mm, sizeof(double)); 
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
    int sl, sk;

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
      F77_NAME(dcopy)(&P, &samples[s*nParams+betaIndx], &incOne, beta, &incOne);
      covExpand(&samples[s*nParams+AIndx], A, m);//note this is K, so we need chol
      F77_NAME(dpotrf)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");} 
      clearUT(A, m); //make sure upper tri is clear
   
      for(k = 0; k < m; k++){
	phi[k] = samples[s*nParams+(phiIndx+k)];
	
	if(covModel == "matern"){
	  nu[k] = samples[s*nParams+(nuIndx+k)];
	}

      }

      sl = sk = 0;
      for(k = 0; k < m; k++){
	sl = 0;
	for(l = 0; l < m; l++){
	  
	  for(i = 0; i < n[k]; i++){
	    for(j = 0; j < n[l]; j++){
	      
	      C[(sl+j)*N+(sk+i)] = 0.0;
	      for(ii = 0; ii < m; ii++){
		C[(sl+j)*N+(sk+i)] += A[k+m*ii]*A[l+m*ii]*spCor(obsD[(sl+j)*N+(sk+i)], phi[ii], nu[ii], covModel);
	      }
	    }
	  }
	  sl += n[l];
	}
	sk += n[k];
      }
      
      //c
      sl = sk = 0;
      for(k = 0; k < m; k++){
	sl = 0;
	for(l = 0; l < m; l++){
	  
	  for(i = 0; i < nPred[k]; i++){
	    for(j = 0; j < n[l]; j++){
	      
	      c[(sl+j)*NPred+(sk+i)] = 0.0;
	      for(ii = 0; ii < m; ii++){
		c[(sl+j)*NPred+(sk+i)] += A[k+m*ii]*A[l+m*ii]*spCor(predObsD[(sl+j)*NPred+(sk+i)], phi[ii], nu[ii], covModel);
	      }
	    }
	  }
	  sl += n[l];
	}
	sk += nPred[k];
      }
            
      F77_NAME(dpotrf)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotri failed\n");}
          
      sl = 0;
      for(l = 0; l < m; l++){
      	for(j = 0; j < nPred[l]; j++){
      	  z[sl+j] = F77_NAME(ddot)(&m, &A[l], &m, &A[l], &m);
      	}
      	sl += nPred[l];
      } 

      //prediction location loop
      for(k = 0; k < NPred; k++){
	
	F77_NAME(dsymv)(lower, &N, &one, C, &N, &c[k], &NPred, &zero, u, &incOne);
	
	muPred = F77_NAME(ddot)(&N, u, &incOne, &w[N*s], &incOne);
	varPred = z[k] - F77_NAME(ddot)(&N, u, &incOne, &c[k], &NPred);
	
	for(i = 0; i < N; i++){
	  if(predObsD[i*NPred+k] == 0){
	    varPred = 0; //fitted value (assuming rnorm returns mu when sd is 0)
	  }
	}

	REAL(predWSamples_r)[s*NPred+k] = rnorm(muPred, sqrt(varPred));
	if(family == "binomial"){
	  REAL(predSamples_r)[s*NPred+k] = 1.0/(1.0+exp(-1.0*(F77_NAME(ddot)(&P, &Z[k], &NPred, beta, &incOne)+REAL(predWSamples_r)[s*NPred+k])));
	}else if(family == "poisson"){
	  REAL(predSamples_r)[s*NPred+k] = exp(F77_NAME(ddot)(&P, &Z[k], &NPred, beta, &incOne)+REAL(predWSamples_r)[s*NPred+k]);
	}else{
	  error("c++ error: family misspecification\n");
	}

	R_CheckUserInterrupt();
      }//end prediction location loop
  
      //report
      if(verbose){
	if(status == nReport){
	  Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
    }
    
    PutRNGstate();
    
    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 2;
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, predSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.y.predictive.samples")); 

    SET_VECTOR_ELT(result_r, 1, predWSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.w.predictive.samples")); 
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
