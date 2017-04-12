#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  SEXP spLMPredict(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP p_r, SEXP Z_r, SEXP q_r,
		   SEXP samples_r, SEXP nSamples_r, 
		   SEXP betaPrior_r, SEXP betaNorm_r, SEXP beta_r, SEXP sigmaSqIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		   SEXP obsD_r, SEXP obsPredD_r, SEXP covModel_r, SEXP nugget_r, SEXP verbose_r, SEXP nReport_r){
    
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
    double *Z = REAL(Z_r);
    int q = INTEGER(q_r)[0];
    int nq = n*q;
    int qp = q*p;
    
    int nSamples = INTEGER(nSamples_r)[0];
    
    double *beta = REAL(beta_r);
    int sigmaSqIndx = INTEGER(sigmaSqIndx_r)[0]; 
    int tauSqIndx  = INTEGER(tauSqIndx_r)[0]; 
    int phiIndx = INTEGER(phiIndx_r)[0]; 
    int nuIndx  = INTEGER(nuIndx_r)[0];
    
    double *obsD = REAL(obsD_r);
    double *obsPredD = REAL(obsPredD_r);
    
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
    
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Prediction at %i locations.\n\n", q);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      if(!nugget){
    	Rprintf("tau.sq not included in the model (i.e., no nugget model).\n\n");
      }
      
    } 
    
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/  
    SEXP predSamples_r;
    PROTECT(predSamples_r = allocMatrix(REALSXP, q, nSamples)); nProtect++; 
    
    int status=1;
    double a, b;

    double *theta = (double *) R_alloc(3, sizeof(double)); //phi, nu, and perhaps more in the future
    double *C = (double *) R_alloc(nn, sizeof(double)); 
    double *c = (double *) R_alloc(nq, sizeof(double)); 
    double *z = (double *) R_alloc(n, sizeof(double)); 
    double *u = (double *) R_alloc(n, sizeof(double)); 
   
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

    double *tmp_p = (double *) R_alloc(p, sizeof(double)); 
    double *tmp_np = NULL;
    double *tmp_qp = NULL;
    double *tmp_nn = NULL;
    double *tmp_nq = NULL;
     
    if(betaPrior == "normal"){
      tmp_np = (double *) R_alloc(np, sizeof(double));
      tmp_qp = (double *) R_alloc(qp, sizeof(double));
      tmp_nn = (double *) R_alloc(nn, sizeof(double));
      tmp_nq = (double *) R_alloc(nq, sizeof(double));

      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, betaMu, &incOne, &zero, z, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, z, &incOne);
       
      F77_NAME(dsymm)(rside, lower, &n, &p, &one, betaC, &p, X, &n, &zero, tmp_np, &n);
      F77_NAME(dgemm)(ntran, ytran, &n, &n, &p, &one, tmp_np, &n, X, &n, &zero, tmp_nn, &n);
 
      F77_NAME(dsymm)(rside, lower, &q, &p, &one, betaC, &p, Z, &q, &zero, tmp_qp, &q);
      F77_NAME(dgemm)(ntran, ytran, &q, &n, &p, &one, tmp_qp, &q, X, &n, &zero, tmp_nq, &q);
    }

    //check if any prediction locations are observed
    int *fitted = (int *) R_alloc(q, sizeof(int));
    for(i = 0; i < q; i++){
      fitted[i] = 0;
      for(j = 0; j < n; j++){
    	if(obsPredD[i*n+j] == 0){
    	  fitted[i] = 1;
    	  break;
    	}
      }
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
      
      theta[0] = sigmaSq[s];
      theta[1] = phi[s];
      
      if(covModel == "matern"){
      	theta[2] = nu[s];
      }
      
      spCovLT(obsD, n, theta, covModel, C);
      spCov(obsPredD, nq, theta, covModel, c);
   
      if(nugget){
      	for(k = 0; k < n; k++){
      	  C[k*n+k] += tauSq[s];
      	}
      }

      if(betaPrior == "normal"){
      	for(k = 0; k < n; k++){
      	  for(l = k; l < n; l++){
      	    C[k*n+l] += tmp_nn[k*n+l];
      	  }
      	}
      
      	for(k = 0; k < n; k++){
      	  for(l = 0; l < q; l++){
      	    c[l*n+k] += tmp_nq[k*q+l];
      	  }
      	}
      }
      
      for(k = 0; k < nq; k++){
      	if(obsPredD[k] == 0 && nugget){
      	  c[k] += tauSq[s];
      	}
      }
	
      F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}//L_1
      
      if(betaPrior == "normal"){
      	F77_NAME(dcopy)(&n, z, &incOne, u, &incOne);
      }else{
      	F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s], &nSamples, &zero, u, &incOne);
       	F77_NAME(daxpy)(&n, &one, Y, &incOne, u, &incOne);
      }

      F77_NAME(dtrsv)(lower, ntran, nUnit, &n, C, &n, u, &incOne);//L_1u = (y-X\mu_beta) or (y-X\beta)
     
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &q, &one, C, &n, c, &n);//L_1v = c_0


      //for each location
      for(j = 0; j < q; j++){
	
	if(betaPrior == "normal"){
	  a = F77_NAME(ddot)(&p, &Z[j], &q, betaMu, &incOne) + F77_NAME(ddot)(&n, u, &incOne, &c[n*j], &incOne);
	  
	  F77_NAME(dsymv)(lower, &p, &one, betaC, &p, &Z[j], &q, &zero, tmp_p, &incOne);
	  if(nugget){
	    b = sqrt(F77_NAME(ddot)(&p, tmp_p, &incOne, &Z[j], &q) + sigmaSq[s] + tauSq[s] - F77_NAME(ddot)(&n, &c[n*j], &incOne, &c[n*j], &incOne));
	  }else{
	    b = sqrt(F77_NAME(ddot)(&p, tmp_p, &incOne, &Z[j], &q) + sigmaSq[s] - F77_NAME(ddot)(&n, &c[n*j], &incOne, &c[n*j], &incOne));	    
	  }
	}else{
	  a = F77_NAME(ddot)(&p, &Z[j], &q, &beta[s], &nSamples) + F77_NAME(ddot)(&n, u, &incOne, &c[n*j], &incOne);
	  
	  if(nugget){
	    b = sqrt(sigmaSq[s] + tauSq[s] - F77_NAME(ddot)(&n, &c[n*j], &incOne, &c[n*j], &incOne));
	  }else{
	    b = sqrt(sigmaSq[s] - F77_NAME(ddot)(&n, &c[n*j], &incOne, &c[n*j], &incOne));
	  }
	  
	}

	if(fitted[j] == 0){
	  REAL(predSamples_r)[s*q+j] = rnorm(a, b);
	}else{
	  REAL(predSamples_r)[s*q+j] = a;
	}	

	R_CheckUserInterrupt();
      	}//end prediction location loop
     
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
