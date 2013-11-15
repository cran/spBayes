#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  
  SEXP spMisalignPredict(SEXP Y_r, SEXP X_r, SEXP m_r, SEXP n_r, SEXP p_r, 
			 SEXP Z_r, SEXP nPred_r, SEXP pPred_r, 
			 SEXP obsD_r, SEXP predObsD_r,
			 SEXP samples_r, SEXP beta_r, SEXP nSamples_r, 
			 SEXP betaPrior_r, SEXP betaNorm_r, 
			 SEXP nugget_r, SEXP covModel_r, 
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
    double *beta = REAL(beta_r);
    int nSamples = INTEGER(nSamples_r)[0];
  
    //priors
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));
    double *betaMu = NULL;
    double *betaC = NULL;
    
    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(P, sizeof(double));
      F77_NAME(dcopy)(&P, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);
      
      betaC = (double *) R_alloc(PP, sizeof(double)); 
      F77_NAME(dcopy)(&PP, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaC, &incOne);
    }
  
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    int nParams, AIndx, PsiIndx, phiIndx, nuIndx;
    
    if(!nugget && covModel != "matern"){
      nParams = nLTr+m;//A, phi
      AIndx = 0; phiIndx = nLTr;
    }else if(nugget && covModel != "matern"){
      nParams = nLTr+m+m;//A, diag(Psi), phi
      AIndx = 0; PsiIndx = nLTr; phiIndx = PsiIndx+m;
    }else if(!nugget && covModel == "matern"){
      nParams = nLTr+2*m;//A, phi, nu
      AIndx = 0; phiIndx = nLTr, nuIndx = phiIndx+m;
    }else{
      nParams = nLTr+3*m;//A, diag(Psi), phi, nu
      AIndx = 0; PsiIndx = nLTr, phiIndx = PsiIndx+m, nuIndx = phiIndx+m;
     }
    
  if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i outcome variables.\n\n", m);
      Rprintf("Number of observations within each outcome:"); printVec(n, m);
      Rprintf("\nNumber of covariates for each outcome (including intercept if specified):"); printVec(p, m);
      Rprintf("\nTotal number of observations: %i\n\n", N);
      Rprintf("Total number of covariates (including intercept if specified): %i\n\n", P);

      Rprintf("Number of prediction for each outcome:"); printVec(nPred, m);   
      Rprintf("\nTotal number of predictions: %i\n\n", NPred);

      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
    
      if(!nugget){
	Rprintf("Psi not included in the model (i.e., no nugget model).\n\n");
      }

    }
 
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    SEXP predSamples_r;
    PROTECT(predSamples_r = allocMatrix(REALSXP, NPred, nSamples)); nProtect++;   
    
    int status=1;
    
    double *C = (double *) R_alloc(NN, sizeof(double));
    double *c = (double *) R_alloc(NPredN, sizeof(double));
    double *z = (double *) R_alloc(NPred, sizeof(double)); new double[NPred];
    double *u = (double *) R_alloc(N, sizeof(double)); new double[N];
    double *v = (double *) R_alloc(N, sizeof(double)); new double[N];
    double muPred, varPred;
    double *A = (double *) R_alloc(mm, sizeof(double)); 
    double *Psi = (double *) R_alloc(m, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
    int *fitted = (int *) R_alloc(NPred, sizeof(int)); zeros(fitted, NPred);
    int sl, sk;

    if(betaPrior == "normal"){
      //fill this in
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
    
    for(s = 0; s < nSamples; s++){
      
      covExpand(&samples[s*nParams+AIndx], A, m);//note this is K, so we need chol
      F77_NAME(dpotrf)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");} 
      clearUT(A, m); //make sure upper tri is clear
   
      for(k = 0; k < m; k++){
	phi[k] = samples[s*nParams+(phiIndx+k)];
	
	if(covModel == "matern"){
	  nu[k] = samples[s*nParams+(nuIndx+k)];
	}

      }

      if(nugget){
	for(k = 0; k < m; k++){
	  Psi[k] = samples[s*nParams+(PsiIndx+k)];
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
      
      sl = 0;
      for(l = 0; l < m; l++){
      	for(j = 0; j < n[l]; j++){
      	  C[(sl+j)*N+(sl+j)] += Psi[l];
      	}
      	sl += n[l];
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

	      if(predObsD[(sl+j)*NPred+(sk+i)] == 0 && l == k){
	      	c[(sl+j)*NPred+(sk+i)] += Psi[l];
	      	fitted[sk+i] = 1;
	      }

	    }
	  }
	  sl += n[l];
	}
	sk += nPred[k];
      }
            
      F77_NAME(dpotrf)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotri failed\n");}
      
      F77_NAME(dgemv)(ntran, &N, &P, &negOne, X, &N, &beta[s*P], &incOne, &zero, v, &incOne);
      F77_NAME(daxpy)(&N, &one, Y, &incOne, v, &incOne);
      
      sl = 0;
      for(l = 0; l < m; l++){
      	for(j = 0; j < nPred[l]; j++){
      	  z[sl+j] = F77_NAME(ddot)(&m, &A[l], &m, &A[l], &m) + Psi[l];
      	}
      	sl += nPred[l];
      } 

      //prediction location loop
      for(k = 0; k < NPred; k++){
	
	F77_NAME(dsymv)(lower, &N, &one, C, &N, &c[k], &NPred, &zero, u, &incOne);
	
	muPred = F77_NAME(ddot)(&P, &Z[k], &NPred, &beta[s*P], &incOne) + F77_NAME(ddot)(&N, u, &incOne, v, &incOne);
	varPred = z[k] - F77_NAME(ddot)(&N, u, &incOne, &c[k], &NPred);
	
	if(fitted[k] == 0){
	  REAL(predSamples_r)[s*NPred+k] = rnorm(muPred, sqrt(varPred));
	}else{//fitted
      	  REAL(predSamples_r)[s*NPred+k] = muPred;
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
    int nResultListObjs = 1;
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, predSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.y.predictive.samples")); 
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
