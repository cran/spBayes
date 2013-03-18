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

  SEXP spPPMvGLMPredict(SEXP family_r, SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP g_r, SEXP p_r, SEXP Z_r, SEXP q_r, SEXP knotsD_r, SEXP predKnotsD_r, 
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
    int g = INTEGER(g_r)[0];
    int p = INTEGER(p_r)[0];
    double *Z = REAL(Z_r);
    int q = INTEGER(q_r)[0];
    int nLTr = m*(m-1)/2+m;

    double *knotsD = REAL(knotsD_r);
    double *predKnotsD = REAL(predKnotsD_r);
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
    int gm = g*m;
    int gmgm = gm*gm;
    int qm = q*m;
    int qmgm = qm*gm;
    int qmqm = qm*qm;
  
    SEXP wPred_r, yPred_r;

    PROTECT(wPred_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 
    double *wPred = REAL(wPred_r);

    PROTECT(yPred_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 
    double *yPred = REAL(yPred_r);
    
    double *S_knots = (double *) R_alloc(gmgm, sizeof(double));
    double *S_predKnots = (double *) R_alloc(qmgm, sizeof(double));

    double *beta = (double *) R_alloc(p, sizeof(double));
    double *A = (double *) R_alloc(mm, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
    
    double *tmp_gm = (double *) R_alloc(gm, sizeof(double));
    double *tmp_qm = (double *) R_alloc(qm, sizeof(double));
     
    int status = 0;
 
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
      for(jj = 0; jj < g; jj++){
      	for(ii = jj; ii < g; ii++){	
      	  for(k = 0; k < m; k++){
      	    for(l = 0; l < m; l++){
      	      S_knots[(k+jj*m)*gm+(ii*m+l)] = 0.0; 
      	      for(h = 0; h < m; h++){
		 S_knots[(k+jj*m)*gm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(knotsD[jj*g+ii], phi[h], nu[h], covModel);
      	      }
      	    }
      	  }
      	}
      }
      // } //parallel for

      // #pragma omp parallel 
      // {
      // #pragma omp for private(ii, k, l, h)
      for(jj = 0; jj < g; jj++){
      	for(ii = 0; ii < q; ii++){	
      	  for(k = 0; k < m; k++){
      	    for(l = 0; l < m; l++){
      	      S_predKnots[(k+jj*m)*qm+(ii*m+l)] = 0.0; 
      	      for(h = 0; h < m; h++){
		 S_predKnots[(k+jj*m)*qm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(predKnotsD[jj*q+ii], phi[h], nu[h], covModel);
      	      }
      	    }
      	  }
      	}
      }
      // } //parallel for
    
      F77_NAME(dpotrf)(lower, &gm, S_knots, &gm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &gm, S_knots, &gm, &info); if(info != 0){error("c++ error: dpotri failed\n");}	 

      F77_NAME(dgemv)(ntran, &qm, &p, &one, Z, &qm, &samples[s*nParams+betaIndx], &incOne, &zero, tmp_qm, &incOne);
    	
      F77_NAME(dsymv)(lower, &gm, &one, S_knots, &gm, &wSamples[s*gm], &incOne, &zero, tmp_gm, &incOne);
      F77_NAME(dgemv)(ntran, &qm, &gm, &one, S_predKnots, &qm, tmp_gm, &incOne, &zero, &wPred[s*qm], &incOne);
            
      if(family == "binomial"){
	for(i = 0; i < qm; i++){
	  yPred[s*qm+i] = 1.0/(1.0+exp(-1.0*(tmp_qm[i]+wPred[s*qm+i])));//rbinom(1, 1.0/(1.0+exp(-1.0*(tmp_qm[i]+wPred[s*qm+i]))));
	}
      }else if(family == "poisson"){
	for(i = 0; i < qm; i++){
	  yPred[s*qm+i] =  exp(tmp_qm[i]+wPred[s*qm+i]);//rpois(exp(tmp_qm[i]+wPred[s*qm+i]));
	}	   
      }else{
	error("c++ error: family misspecification\n");
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
