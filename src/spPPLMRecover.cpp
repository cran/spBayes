#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"


extern "C" {

  SEXP spPPLMRecover(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP p_r, SEXP m_r,
		     SEXP samples_r, SEXP nSamples_r, 
		     SEXP beta_r, SEXP sigmaSqIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		     SEXP nugget_r, SEXP knotsD_r, SEXP knotsObsD_r, SEXP covModel_r, SEXP modPP_r, SEXP verbose_r, SEXP nReport_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int i, j, k, l, s, info, nProtect=0;
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
    int m = INTEGER(m_r)[0];
    int nm = n*m;
    int mm = m*m;
    int mp = m*p;
    
    int nSamples = INTEGER(nSamples_r)[0];

    double *beta = REAL(beta_r);
    int sigmaSqIndx = INTEGER(sigmaSqIndx_r)[0]; 
    int tauSqIndx  = INTEGER(tauSqIndx_r)[0]; 
    int phiIndx = INTEGER(phiIndx_r)[0]; 
    int nuIndx  = INTEGER(nuIndx_r)[0]; 

    double *knotsD = REAL(knotsD_r);
    double *knotsObsD = REAL(knotsObsD_r);

    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    bool modPP = static_cast<bool>(INTEGER(modPP_r)[0]);
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status=1;
    SEXP wStrSamples_r; 
    SEXP wSamples_r;
    PROTECT(wStrSamples_r = allocMatrix(REALSXP, m, nSamples)); nProtect++; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, n, nSamples)); nProtect++; 
  
    double *theta = (double *) R_alloc(3, sizeof(double)); //phi, nu, and perhaps more in the future
    double *P = (double *) R_alloc(nm, sizeof(double)); 
    double *K = (double *) R_alloc(mm, sizeof(double)); 

    double *sigmaSq = &REAL(samples_r)[sigmaSqIndx*nSamples];
    double *tauSq = NULL;
    if(nugget){
      tauSq = &REAL(samples_r)[tauSqIndx*nSamples];
    }
    double *phi = &REAL(samples_r)[phiIndx*nSamples];
    double *nu = NULL;
    if(covModel == "matern"){
      nu = &REAL(samples_r)[nuIndx*nSamples];
    }

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); 
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double)); 
    double *tmp_nm2 = (double *) R_alloc(nm, sizeof(double)); 
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double *tmp_m2 = (double *) R_alloc(m, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));

    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tRecovering w\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }
    
    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){    

      theta[0] = sigmaSq[s];
      theta[1] = phi[s];
      
      if(covModel == "matern"){
	theta[2] = nu[s];
      }
      
      spCovLT(knotsD, m, theta, covModel, K);
      spCov(knotsObsD, nm, theta, covModel, P);
      
      F77_NAME(dpotrf)(lower, &m, K, &m, &info); if(info != 0){error("c++ error: Cholesky failed 1\n");}
      F77_NAME(dpotri)(lower, &m, K, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed\n");}
      
      //P K^{-1}
      F77_NAME(dsymm)(lside, lower, &m, &n, &one, K, &m, P, &m, &zero, tmp_nm, &m);
      
      if(!modPP){
	for(i = 0; i < n; i++){
	  tmp_n[i] = 1.0/tauSq[s];
	}
      }else{
	for(i = 0; i < n; i++){
	  if(nugget){
	    tmp_n[i] = 1.0/(tauSq[s]+sigmaSq[s]-F77_NAME(ddot)(&m, &P[m*i], &incOne, &tmp_nm[m*i], &incOne));
	  }else{
	    tmp_n[i] = 1.0/(sigmaSq[s]-F77_NAME(ddot)(&m, &P[m*i], &incOne, &tmp_nm[m*i], &incOne));
	  }
	}

      }

      for(j = 0; j < n; j++){
	for(i = 0; i < m; i++){
	  tmp_nm2[j*m+i] = tmp_nm[j*m+i]*tmp_n[j];
	}
      }
 
      //(C^{*-1} c) (1/E ct C^{*-1})
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &n, &one, tmp_nm2, &m, tmp_nm, &m, &zero, tmp_mm, &m);

      for(i = 0; i < mm; i++){
	K[i] += tmp_mm[i];
      }
     
      //invert C_str
      F77_NAME(dpotrf)(lower, &m, K, &m, &info); if(info != 0){error("c++ error: Cholesky failed 2\n");}
      F77_NAME(dpotri)(lower, &m, K, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed\n");}
      
      //make w* mu
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s], &nSamples, &zero, tmp_n, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
      
      //(1/E ct C^{*-1})'(Y-XB)
      F77_NAME(dgemv)(ntran, &m, &n, &one, tmp_nm2, &m, tmp_n, &incOne, &zero, tmp_m, &incOne);     
      F77_NAME(dsymv)(lower, &m, &one, K, &m, tmp_m, &incOne, &zero, tmp_m2, &incOne);

      //chol for the mvnorm and draw
      F77_NAME(dpotrf)(lower, &m, K, &m, &info); if(info != 0){error("c++ error: Cholesky failed 3\n");}
      mvrnorm(&REAL(wStrSamples_r)[s*m], tmp_m2, K, m, false);
      
      //make \tild{w}
      F77_NAME(dgemv)(ytran, &m, &n, &one, tmp_nm, &m, &REAL(wStrSamples_r)[s*m], &incOne, &zero, &REAL(wSamples_r)[s*n], &incOne);
   
      R_CheckUserInterrupt();
      
      if(verbose){
	if(status == nReport){
	  Rprintf("Sampled: %i of %i, %3.2f%%\n", s+1, nSamples, 100.0*s/nSamples);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      R_CheckUserInterrupt();
    }//end sample loop
    
    PutRNGstate();
    
    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 2;
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, wSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.w.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, wStrSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.wStr.samples")); 
    
    namesgets(result_r, resultName_r);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
