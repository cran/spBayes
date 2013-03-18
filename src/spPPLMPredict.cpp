#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP spPPLMPredict(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP p_r, SEXP m_r, SEXP Z_r, SEXP q_r,
		     SEXP samples_r, SEXP nSamples_r, 
		     SEXP beta_r, SEXP sigmaSqIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		     SEXP nugget_r, SEXP knotsD_r, SEXP knotsObsD_r, SEXP knotsPredD_r, SEXP covModel_r, SEXP modPP_r, SEXP verbose_r, SEXP nReport_r){
    
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
    double *Z = REAL(Z_r);
    int q = INTEGER(q_r)[0];
    int qm = q*m;
    
    int nSamples = INTEGER(nSamples_r)[0];

    double *beta = REAL(beta_r);
    int sigmaSqIndx = INTEGER(sigmaSqIndx_r)[0]; 
    int tauSqIndx  = INTEGER(tauSqIndx_r)[0]; 
    int phiIndx = INTEGER(phiIndx_r)[0]; 
    int nuIndx  = INTEGER(nuIndx_r)[0]; 

    double *knotsD = REAL(knotsD_r);
    double *knotsObsD = REAL(knotsObsD_r);
    double *knotsPredD = REAL(knotsPredD_r);
  
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    bool modPP = static_cast<bool>(INTEGER(modPP_r)[0]);
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
      
      if(modPP){
	Rprintf("Using modified predictive process with %i knots.\n\n", m);
      }else{
	Rprintf("Using non-modified predictive process with %i knots.\n\n", m);
      }
   
    } 

    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/ 
    SEXP predSamples_r;
    PROTECT(predSamples_r = allocMatrix(REALSXP, q, nSamples)); nProtect++; 
  
    int status=1;
    double a, b;

    double *theta = (double *) R_alloc(3, sizeof(double)); //phi, nu, and perhaps more in the future
    double *P = (double *) R_alloc(nm, sizeof(double)); 
    double *K = (double *) R_alloc(mm, sizeof(double)); 
    double *D = (double *) R_alloc(n, sizeof(double)); 
    double *O = (double *) R_alloc(qm, sizeof(double)); 
    double *H = (double *) R_alloc(nm, sizeof(double)); 
    double *z = (double *) R_alloc(n, sizeof(double)); 
    double *L = (double *) R_alloc(mm, sizeof(double));

    double *sigmaSq = &REAL(samples_r)[sigmaSqIndx*nSamples];

    double *tauSq = NULL;
    if(nugget){
      tauSq =&REAL(samples_r)[tauSqIndx*nSamples];
    }
    double *phi = &REAL(samples_r)[phiIndx*nSamples];
    double *nu = NULL;
    if(covModel == "matern"){
      nu = &REAL(samples_r)[nuIndx*nSamples];
    }

    double *tmp_n = (double *) R_alloc(n, sizeof(double)); 
    double *tmp_n2 = (double *) R_alloc(n, sizeof(double)); 
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double *tmp_m2 = (double *) R_alloc(m, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));

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
      
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s], &nSamples, &zero, z, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, z, &incOne);
      
      theta[0] = sigmaSq[s];
      theta[1] = phi[s];
      
      if(covModel == "matern"){
	theta[2] = nu[s];
      }
      
      spCovLT(knotsD, m, theta, covModel, K);
      spCov(knotsObsD, nm, theta, covModel, P);
      spCov(knotsPredD, qm, theta, covModel, O);
      
      F77_NAME(dcopy)(&nm, P, &incOne, H, &incOne);

      //get a copy of K
      for(k = 0; k < m; k++){
	for(l = k; l < m; l++){
	  L[k*m+l] = K[k*m+l];
	}
      }

      F77_NAME(dpotrf)(lower, &m, K, &m, &info);//L_K
    
      //get D
      if(modPP){
	F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &m, &n, &one, K, &m, H, &m);
	
	for(k = 0; k < n; k++){
	  D[k] = sigmaSq[s] - F77_NAME(ddot)(&m, &H[k*m], &incOne, &H[k*m], &incOne);
	  if(nugget){
	    D[k] += tauSq[s];
	  }
	}
	
      }else{
	for(k = 0; k < n; k++){
	  D[k] = tauSq[s];
	}
      }	 
      
      //get H
      for(k = 0; k < n; k++){
	for(l = 0; l < m; l++){
	  H[k*m+l] = P[k*m+l]/sqrt(D[k]);//W'
	}
      }
      
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &n, &one, H, &m, H, &m, &zero, tmp_mm, &m);//W'W
      
      for(k = 0; k < m; k++){
      	for(l = k; l < m; l++){
      	  L[k*m+l] += tmp_mm[k*m+l];
      	}
      }
      
      F77_NAME(dpotrf)(lower, &m, L, &m, &info); if(info != 0){error("c++ error: dpotrf failed\n");}//L
      
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &m, &n, &one, L, &m, H, &m);//LH = W'
	  
      //get [v_1:v_2: ... : v_n]
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &m, &n, &one, K, &m, P, &m);//L_K [v_1:v_2: ... : v_n] = [c_1:c_2: ... : c_n]
      
      //get v
      for(k = 0; k < n; k++){
	tmp_n2[k] = z[k]/sqrt(D[k]);//D^{-1/2}(y-XB)
      }

      //for each location
      for(j = 0; j < q; j++){
	
      	F77_NAME(dtrsv)(lower, ntran, nUnit, &m, K, &m, &O[j*m], &incOne);//L_K u = c_0
	
	//get u
      	for(k = 0; k < n; k++){
      	  tmp_n[k] = F77_NAME(ddot)(&m, &P[k*m], &incOne, &O[j*m], &incOne)/sqrt(D[k]);//u = D^{-1/2}\tidle{c_0}
      	}
	
	//get w
	F77_NAME(dgemv)(ntran, &m, &n, &one, H, &m, tmp_n, &incOne, &zero, tmp_m, &incOne);//w = Hu

	a = F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n, &incOne) - F77_NAME(ddot)(&m, tmp_m, &incOne, tmp_m, &incOne);

	//get z
	F77_NAME(dgemv)(ntran, &m, &n, &one, H, &m, tmp_n2, &incOne, &zero, tmp_m2, &incOne); //z = HD^{-1/2}(y-XB)

	a = sigmaSq[s] - a;
	if(nugget){
	  a += tauSq[s];
	}

	b = F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n2, &incOne) - F77_NAME(ddot)(&m, tmp_m, &incOne, tmp_m2, &incOne);

	REAL(predSamples_r)[s*q+j] = rnorm(F77_NAME(ddot)(&p, &Z[j], &q, &beta[s], &nSamples) + b, sqrt(a));
	
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
