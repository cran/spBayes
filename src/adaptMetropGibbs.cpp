#include <algorithm>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

double ltd(SEXP ltd_r, SEXP theta_r, SEXP rho_r){
  
  SEXP call_r, tmp_r, result_r;
  double result;
  
  PROTECT(call_r = lang2(ltd_r, theta_r));
  PROTECT(tmp_r = eval(call_r, rho_r));
  
  if(!isNumeric(tmp_r))
    error("ltd: result of function call must be numeric");
  
  if (LENGTH(tmp_r) != 1)
    error("ltd: result of function call must be scalar");
  
  PROTECT(result_r = coerceVector(tmp_r, REALSXP));
  
  result = REAL(result_r)[0];
  
  if (result == R_PosInf)
    error("ltd: function returned +Inf");
  
  if (R_IsNaN(result) || R_IsNA(result))
    error("ltd: function returned NA or NaN");
  
  
  UNPROTECT(3);
  
  return result;
} 

extern "C" {

  SEXP adaptMetropGibbs(SEXP ltd_r, SEXP starting_r, SEXP tuning_r,
		    SEXP acceptRate_r,
		    SEXP nBatch_r, SEXP batchLength_r,
		    SEXP verbose_r, SEXP nTheta_r, 
		    SEXP reportBatch, SEXP rho_r){ 
 
    int verbose = INTEGER(verbose_r)[0];
    int nTheta = INTEGER(nTheta_r)[0];
    int nBatch = INTEGER(nBatch_r)[0];
    int batchLength = INTEGER(batchLength_r)[0];

    int b, i, j, s, k, nProtect = 0, inc = 1, info;
    int nSamples = nBatch*batchLength;

    double *accept = (double *) R_alloc(nTheta, sizeof(double));
    double *tuning = REAL(tuning_r);
    
    for(i = 0; i < nTheta; i++){
      accept[i] = 0.0;
      tuning[i] = log(tuning[i]);
    }

    SEXP acceptance_r, samples_r, theta_r, thetaCand_r;
    PROTECT(samples_r = allocMatrix(REALSXP, nTheta, nSamples)); nProtect++; 
    PROTECT(theta_r = allocVector(REALSXP, nTheta)); nProtect++; 
    PROTECT(acceptance_r = allocMatrix(REALSXP, nTheta, nBatch)); nProtect++; 
     
    double ltdCand, ltdCurrent, thetajCurrent, status = 0;
   
    F77_NAME(dcopy)(&nTheta, REAL(starting_r), &inc, REAL(theta_r), &inc);
    
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    ltdCurrent = ltd(ltd_r, theta_r, rho_r);

    GetRNGstate();
    
    for(b = 0, s = 0; b < nBatch; b++){

      for(i = 0; i < batchLength; i++, s++){
	
	for(j = 0; j < nTheta; j++){
	  
	  thetajCurrent = REAL(theta_r)[j];

	  REAL(theta_r)[j] = rnorm(thetajCurrent, exp(tuning[j]));

	  ltdCand = ltd(ltd_r, theta_r, rho_r);

	  if(runif(0.0, 1.0) <= exp(ltdCand - ltdCurrent)){
	    ltdCurrent = ltdCand;
	    accept[j]++;
	  }else{
	    REAL(theta_r)[j] = thetajCurrent;
	  }

	}
	
	F77_NAME(dcopy)(&nTheta, REAL(theta_r), &inc, &REAL(samples_r)[nTheta*s], &inc);
      }//end batch
      
      //adjust tuning
      for(j = 0; j < nTheta; j++){
	REAL(acceptance_r)[b*nTheta+j] = accept[j]/batchLength;
	
	if(accept[j]/batchLength > REAL(acceptRate_r)[j]){
	  tuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}else{
	  tuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}
	accept[j] = 0.0;
      }
      
      //report
      if(verbose){
	if(status == INTEGER(reportBatch)[0]){
	  Rprintf("Batch: %i of %i\n", b, nBatch);
	  Rprintf("Metropolis batch acceptance rate:\n");
	  for(j = 0, i = 0; j < nTheta; j++, i++){
	    Rprintf("%1.3f\t", REAL(acceptance_r)[b*nTheta+j]);
	    if(i == 5){
	      Rprintf("\n");
	      i = 0;
	    }
	  }
	  Rprintf("\n-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      R_CheckUserInterrupt();

    }//end batches
 
    PutRNGstate();    

    //return stuff
    SEXP result_r, resultNames_r;
            
    PROTECT(result_r = allocVector(VECSXP, 2)); nProtect++;
    PROTECT(resultNames_r = allocVector(VECSXP, 2)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, samples_r);
    SET_VECTOR_ELT(resultNames_r, 0, mkChar("p.theta.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, acceptance_r);
    SET_VECTOR_ELT(resultNames_r, 1, mkChar("acceptance"));

    namesgets(result_r, resultNames_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  } 
  
}
