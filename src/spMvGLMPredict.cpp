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

  SEXP spMvGLMPredict(SEXP family_r, SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP p_r, SEXP Z_r, SEXP q_r, SEXP obsD_r, SEXP obsPredD_r,
		      SEXP samples_r, SEXP wSamples_r, SEXP nSamples_r, SEXP covModel_r, SEXP verbose_r, SEXP nReport_r){
    

    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, b, s, ii, jj, info, nProtect= 0;
    const char *lower = "L";
    const char *upper = "U";
    const char *ntran = "N";
    const char *ytran = "T";
    const char *rside = "R";
    const char *lside = "L";
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
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int p = INTEGER(p_r)[0];
    double *Z = REAL(Z_r);
    int q = INTEGER(q_r)[0];
    int nLTr = m*(m-1)/2+m;

    double *obsD = REAL(obsD_r);
    double *obsPredD = REAL(obsPredD_r);
    double *samples = REAL(samples_r);
    double *wSamples = REAL(wSamples_r);
    
    int nSamples = INTEGER(nSamples_r)[0];
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    bool verbose = static_cast<bool>(INTEGER(verbose_r)[0]);
    int nReport = INTEGER(nReport_r)[0];
 
    /*****************************************
         Set-up sample matrices etc.
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
    
    int mm = m*m;
    int nm = n*m;
    int nmnm = nm*nm;
    int qm = q*m;
    int qmnm = qm*nm;
    int qmqm = qm*qm;
    int nmm = nm*m;

    SEXP wPred_r, yPred_r;

    PROTECT(wPred_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 
    double *wPred = REAL(wPred_r);

    PROTECT(yPred_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 
    double *yPred = REAL(yPred_r);
    
    double *S_obs = (double *) R_alloc(nmnm, sizeof(double));
    double *S_obsPred = (double *) R_alloc(qmnm, sizeof(double));

    double *beta = (double *) R_alloc(p, sizeof(double));
    double *A = (double *) R_alloc(mm, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
    
    double *tmp_nmm = (double *) R_alloc(nmm, sizeof(double));
    double *tmp_qm = (double *) R_alloc(qm, sizeof(double));
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *K = (double *) R_alloc(mm, sizeof(double));
    
    int status = 0;
    
    //check if any prediction locations are observed
    int *fitted = (int *) R_alloc(q, sizeof(int));
    for(i = 0; i < q; i++){
      fitted[i] = -1;
      for(j = 0; j < n; j++){
    	if(obsPredD[i*n+j] == 0){
    	  fitted[i] = j;
    	  break;
    	}
      }
    }
    
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tStarting prediction\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }
    
    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){
      
      F77_NAME(dcopy)(&p, &samples[s*nParams+betaIndx], &incOne, beta, &incOne);
      covExpand(&samples[s*nParams+AIndx], A, m);
      F77_NAME(dpotrf)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotrf failed\n");}   
      clearUT(A, m);
      
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, A, &m, A, &m, &zero, K, &m);
      
      for(i = 0; i < m; i++){
	phi[i] = samples[s*nParams+phiIndx+i];
	
	if(covModel == "matern"){
	  nu[i] = samples[s*nParams+nuIndx+i];
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
      	      S_obs[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
      	      for(h = 0; h < m; h++){
		 S_obs[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(obsD[jj*n+ii], phi[h], nu[h], covModel);
      	      }
      	    }
      	  }
      	}
      }
      // } //parallel for

      // #pragma omp parallel 
      // {
      // #pragma omp for private(ii, k, l, h)
      for(jj = 0; jj < q; jj++){
      	for(ii = 0; ii < n; ii++){	
      	  for(k = 0; k < m; k++){
      	    for(l = 0; l < m; l++){
      	      S_obsPred[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
      	      for(h = 0; h < m; h++){
		 S_obsPred[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(obsPredD[jj*n+ii], phi[h], nu[h], covModel);
      	      }
      	    }
      	  }
      	}
      }
      // } //parallel for
    
      F77_NAME(dpotrf)(lower, &nm, S_obs, &nm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &nm, S_obs, &nm, &info); if(info != 0){error("c++ error: dpotri failed\n");}	 

      F77_NAME(dgemv)(ntran, &qm, &p, &one, Z, &qm, &samples[s*nParams+betaIndx], &incOne, &zero, tmp_qm, &incOne);
    
      for(j = 0; j < q; j++){

	if(fitted[j] == -1){
	  
	  //get Mu
	  F77_NAME(dsymm)(lside, lower, &nm, &m, &one, S_obs, &nm, &S_obsPred[(j*m)*nm], &nm, &zero, tmp_nmm, &nm);
	  F77_NAME(dgemv)(ytran, &nm, &m, &one, tmp_nmm, &nm, &wSamples[s*nm], &incOne, &zero, tmp_m, &incOne);
	  
	  //get Sigma
	  F77_NAME(dgemm)(ytran, ntran, &m, &m, &nm, &one, tmp_nmm, &nm, &S_obsPred[(j*m)*nm], &nm, &zero, tmp_mm, &m);
	  
	  for(i = 0; i < mm; i++){
	    tmp_mm[i] = K[i] - tmp_mm[i];
	  }
	  
	  F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	  mvrnorm(&wPred[s*qm+j*m], tmp_m, tmp_mm, m, false);
	  
	}else{
	  F77_NAME(dcopy)(&m, &wSamples[s*nm+m*fitted[j]], &incOne, &wPred[s*qm+j*m], &incOne);
	}

	if(family == "binomial"){
	  for(i = 0; i < m; i++){
	    yPred[s*qm+j*m+i] = 1.0/(1.0+exp(-1.0*(tmp_qm[j*m+i]+wPred[s*qm+j*m+i])));//rbinom(1, 1.0/(1.0+exp(-1.0*(tmp_qm[i]+wPred[s*qm+i]))));
	  }
	}else if(family == "poisson"){
	  for(i = 0; i < m; i++){
	    yPred[s*qm+j*m+i] =  exp(tmp_qm[j*m+i]+wPred[s*qm+j*m+i]);//rpois(exp(tmp_qm[i]+wPred[s*qm+i]));
	  }	   
	}else{
	  error("c++ error: family misspecification\n");
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
      
     } //end sample loop
         
     
     PutRNGstate();

     //make return object
     SEXP result, resultNames;
     int nResultListObjs = 0;
     nResultListObjs = 2;
     
     PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
     PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
          
     SET_VECTOR_ELT(result, 0, wPred_r);
     SET_VECTOR_ELT(resultNames, 0, mkChar("p.w.predictive.samples"));

     SET_VECTOR_ELT(result, 1, yPred_r);
     SET_VECTOR_ELT(resultNames, 1, mkChar("p.y.predictive.samples"));
     
     namesgets(result, resultNames);
     
     //unprotect
     UNPROTECT(nProtect);
     
     return(result);

  }
}
