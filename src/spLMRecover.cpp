#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP spLMRecover(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
		   SEXP samples_r, SEXP nSamples_r, 
		   SEXP sigmaSqIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		   SEXP betaPrior_r, SEXP betaNorm_r, 	   
		   SEXP nugget_r, SEXP covModel_r, 
		   SEXP beta_r, SEXP w_r,
		   SEXP verbose_r, SEXP nReport_r){

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
    
    double *coordsD = REAL(coordsD_r);

    double *samples = REAL(samples_r);
    int nSamples = INTEGER(nSamples_r)[0];

    int sigmaSqIndx = INTEGER(sigmaSqIndx_r)[0]; 
    int phiIndx = INTEGER(phiIndx_r)[0]; 
    int nuIndx  = INTEGER(nuIndx_r)[0]; 
    int tauSqIndx  = INTEGER(tauSqIndx_r)[0]; 

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
     
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    bool getBeta = static_cast<bool>(INTEGER(beta_r)[0]);
    bool getW = static_cast<bool>(INTEGER(w_r)[0]);
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    SEXP betaSamples_r, wSamples_r;
    if(getBeta){
      PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++;
    }
    
    if(getW){ 
      PROTECT(wSamples_r = allocMatrix(REALSXP, n, nSamples)); nProtect++; 
    }

    int status=1;
    double *C = (double *) R_alloc(nn, sizeof(double)); zeros(C, nn);
    double *theta = (double *) R_alloc(3, sizeof(double)); //phi, nu, and perhaps more in the future
    double *B = (double *) R_alloc(pp, sizeof(double));
    double *bb = (double *) R_alloc(p, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double tauSq;
  
    double *W = (double *) R_alloc(nn, sizeof(double)); zeros(W, nn);

    int p1 = p+1;
    double *vU = (double *) R_alloc(n*p1, sizeof(double));
    
    double *betaCInv = NULL;
    double *betaCInvMu = NULL;
    double *Sbeta = NULL;
  
    if(betaPrior == "normal"){
      betaCInv = (double *) R_alloc(pp, sizeof(double));
      betaCInvMu = (double *) R_alloc(p, sizeof(double));
      
      F77_NAME(dcopy)(&pp, betaC, &incOne, betaCInv, &incOne);
      F77_NAME(dpotrf)(lower, &p, betaCInv, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &p, betaCInv, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
      
      F77_NAME(dsymv)(lower, &p, &one, betaCInv, &p, betaMu, &incOne, &zero, betaCInvMu, &incOne);  
    }
    
    if(verbose){
      if(getW){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\t\tRecovering beta and w\n");
	Rprintf("-------------------------------------------------\n");
      }else{
	Rprintf("-------------------------------------------------\n");
	Rprintf("\t\tRecovering beta\n");
	Rprintf("-------------------------------------------------\n");
      }
      #ifdef Win32
      R_FlushConsole();
      #endif
    }
    
    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){
      
      theta[0] = samples[sigmaSqIndx*nSamples+s];
      theta[1] = samples[phiIndx*nSamples+s];
    
      if(covModel == "matern"){
	theta[2] = samples[nuIndx*nSamples+s];
      }
    
      tauSq = samples[tauSqIndx*nSamples+s];
  
      //construct covariance matrix
      spCovLT(coordsD, n, theta, covModel, C);
        
      if(nugget){
	for(k = 0; k < n; k++){
	  C[k*n+k] += tauSq;
	}
      }

      F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}

      F77_NAME(dcopy)(&n, Y, &incOne, vU, &incOne);
      F77_NAME(dcopy)(&np, X, &incOne, &vU[n], &incOne);
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p1, &one, C, &n, vU, &n);//L[v:U] = [y:X]
      
      //B
      F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, &vU[n], &n, &vU[n], &n, &zero, B, &p); //U'U
      
      if(betaPrior == "normal"){
	for(k = 0; k < p; k++){
	  for(l = k; l < p; l++){
	    B[k*p+l] += betaCInv[k*p+l];
	  }
	}
      }
      
      F77_NAME(dpotrf)(lower, &p, B, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &p, B, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	      
      //bb
      F77_NAME(dgemv)(ytran, &n, &p, &one, &vU[n], &n, vU, &incOne, &zero, tmp_p, &incOne); //U'v
      
      if(betaPrior == "normal"){
	for(k = 0; k < p; k++){
	  tmp_p[k] += betaCInvMu[k];
	}
      }
      
      F77_NAME(dsymv)(lower, &p, &one, B, &p, tmp_p, &incOne, &zero, bb, &incOne); 
      F77_NAME(dpotrf)(lower, &p, B, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}

      mvrnorm(&REAL(betaSamples_r)[s*p], bb, B, p, false);

      //get w
      if(getW){
	
	if(nugget){

	  spCovLT(coordsD, n, theta, covModel, C);

	  //v3
	  zeros(W, nn);

	  for(k = 0; k < n; k++){
	    C[k*n+k] += tauSq;
	    W[k*n+k] = tauSq;
	  }
	  
	  //L
	  F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}

	  //W
	  F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &n, &one, C, &n, W, &n);

	  //L_B
	  F77_NAME(dgemm)(ytran, ntran, &n, &n, &n, &negOne, W, &n, W, &n, &zero, C, &n);

	  for(k = 0; k < n; k++){
	    C[k*n+k] += tauSq;
	  }

	  F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}


	  F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &REAL(betaSamples_r)[s*p], &incOne, &zero, vU, &incOne);
	  F77_NAME(daxpy)(&n, &one, Y, &incOne, vU, &incOne);

	  for(k = 0; k < n; k++){
	    vU[k] *= 1.0/tauSq;
	  }

	  F77_NAME(dtrmv)(lower, ytran, nUnit, &n, C, &n, vU, &incOne);
	  F77_NAME(dtrmv)(lower, ntran, nUnit, &n, C, &n, vU, &incOne);

	  for(k = 0; k < n; k++){
	    vU[n+k] = rnorm(0, 1);
	  }
	  
	  F77_NAME(dtrmv)(lower, ntran, nUnit, &n, C, &n, &vU[n], &incOne);

	  for(k = 0; k < n; k++){
	    REAL(wSamples_r)[s*n+k] = vU[k] + vU[n+k];
	  }

	  // //v2
	  // for(k = 0; k < n; k++){
	  //   C[k*n+k] += tauSq;
	  // }

	  // F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	  // F77_NAME(dpotri)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotri failed\n");}

	  // for(i = 0; i < nn; i++){
	  //   C[i] = -1.0*tauSq*C[i]*tauSq;
	  // }

	  // for(i = 0; i < n; i++){
	  //   C[i*n+i] += tauSq;
	  // }

	  // //v
	  // F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &REAL(betaSamples_r)[s*p], &incOne, &zero, vU, &incOne);
	  // F77_NAME(daxpy)(&n, &one, Y, &incOne, vU, &incOne);

	  // for(i = 0; i < n; i++){
	  //   vU[i] *= 1.0/tauSq;
	  // }

	  // F77_NAME(dsymv)(lower, &n, &one, C, &n, vU, &incOne, &zero, &vU[n], &incOne);
	  
	  // F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	  // mvrnorm(&REAL(wSamples_r)[s*n], &vU[n], C, n, false);

	  //v1
	  // for(k = 0; k < n; k++){
	  //   C[k*n+k] += 1.0/tauSq;
	  // }

	  // F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	  // F77_NAME(dpotri)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	  
	  // F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &REAL(betaSamples_r)[s*p], &incOne, &zero, vU, &incOne);
	  // F77_NAME(daxpy)(&n, &one, Y, &incOne, vU, &incOne);
	
	  // for(k = 0; k < n; k++){
	  //   vU[k] *= 1.0/tauSq;
	  // }
	
	  // F77_NAME(dsymv)(lower, &n, &one, C, &n, vU, &incOne, &zero, &vU[n], &incOne);
	
	  // F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	  // mvrnorm(&REAL(wSamples_r)[s*n], &vU[n], C, n, false);
	  
	}else{
	  F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &REAL(betaSamples_r)[s*p], &incOne, &zero, &REAL(wSamples_r)[s*n], &incOne);
	  F77_NAME(daxpy)(&n, &one, Y, &incOne, &REAL(wSamples_r)[s*n], &incOne);
	}

      }

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
      R_CheckUserInterrupt();
    }//end sample loop

    PutRNGstate();
    
    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 1;

    if(getW){
      nResultListObjs++;
    }
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.samples")); 
    
    if(getW){
      SET_VECTOR_ELT(result_r, 1, wSamples_r);
      SET_VECTOR_ELT(resultName_r, 1, mkChar("p.w.samples"));
    }

    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
