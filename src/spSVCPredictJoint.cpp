#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {

  SEXP spSVCPredictJoint(SEXP m_r, SEXP n_r, SEXP KDiag_r, SEXP obsD_r, SEXP predObsD_r, SEXP predD_r, SEXP q_r,
			 SEXP samples_r, SEXP wSamples_r, SEXP nSamples_r, 
			 SEXP AIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 	   
			 SEXP covModel_r, 
			 SEXP verbose_r, SEXP nReport_r, SEXP nThreads_r){

    /*****************************************
                Common variables
    *****************************************/
    int i, j, k, l, b, s, h, info, nProtect=0;
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
    double *obsD = REAL(obsD_r);
    double *predObsD = REAL(predObsD_r);
    double *predD = REAL(predD_r);
    int m = INTEGER(m_r)[0];
    int mm = m*m;
    int n = INTEGER(n_r)[0];
    int nn = n*n;
    int nm = n*m;
    int nmnm = nm*nm;
    int q = INTEGER(q_r)[0];//number of prediction locations
    int qm = q*m;
    int qmnm = qm*nm;
    int qmqm = qm*qm;
    bool KDiag = static_cast<bool>(INTEGER(KDiag_r)[0]);
    int nLTr = m*(m-1)/2+m;
	
    double *samples = REAL(samples_r);
    double *wSamples = REAL(wSamples_r);
    int nSamples = INTEGER(nSamples_r)[0];

    int AIndx = INTEGER(AIndx_r)[0]; 
    int phiIndx = INTEGER(phiIndx_r)[0]; 
    int nuIndx  = INTEGER(nuIndx_r)[0]; 

    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    SEXP wPredSamples_r;
    PROTECT(wPredSamples_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 

    int status=1;
    double *A = (double *) R_alloc(mm, sizeof(double)); zeros(A, mm); //to simplify a future move to the more general cross-cov model
    double *K = (double *) R_alloc(nmnm, sizeof(double)); 
    double *B = (double *) R_alloc(qmnm, sizeof(double)); 
    double *C = (double *) R_alloc(qmqm, sizeof(double));
    double *tmp_nltr = (double *) R_alloc(nLTr, sizeof(double)); 
    double *tmp_qmnm = (double *) R_alloc(qmnm, sizeof(double)); 
    double *tmp_qm = (double *) R_alloc(qm, sizeof(double));
    double *tmp_qmqm = (double *) R_alloc(qmqm, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double)); zeros(nu, m); //this just remains empty of not matern

    double maxNu = 0; //needed for thread safe bessel
    
    if(covModel == "matern"){
      for(s = 0; s < nSamples; s++){
	for(i = 0; i < m; i++){
	  if(samples[(nuIndx+i)*nSamples+s] > maxNu){
	    maxNu = samples[(nuIndx+i)*nSamples+s];
	  }
	}
      }
    }

    int threadID = 0;
    int bessel_ws_inc = static_cast<int>(1.0+maxNu);
    double *bessel_ws = (double *) R_alloc(nThreads*bessel_ws_inc, sizeof(double));

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
      if(verbose){
    Rprintf("Source compiled with OpenMP, posterior sampling is using %i thread(s).\n", nThreads);
      }
#else
    if(nThreads > 1){
      warning("n.omp.threads = %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif  
    
    if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\tJoint sampling of predicted w\n");
	Rprintf("-------------------------------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }
    
    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){
      
      if(KDiag == false){
	dcopy_(&nLTr, &samples[AIndx*nSamples+s], &nSamples, tmp_nltr, &incOne);
      	covExpand(tmp_nltr, A, m);//note this is K, so we need chol
      	F77_NAME(dpotrf)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");} 
      	clearUT(A, m); //make sure upper tri is clear
      }

      for(k = 0; k < m; k++){

	if(KDiag){
	  A[k*m+k] = sqrt(samples[(AIndx+k)*nSamples+s]);
	}
	
	phi[k] = samples[(phiIndx+k)*nSamples+s]; 
	
	if(covModel == "matern"){
	  nu[k] = samples[(nuIndx+k)*nSamples+s]; 
	}
	
      }
      
      //construct covariance matrix
#ifdef _OPENMP
#pragma omp parallel for private(i, k, l, h, threadID)
#endif
      for(j = 0; j < n; j++){
#ifdef _OPENMP
	threadID = omp_get_thread_num();
#endif
	for(i = 0; i < n; i++){	
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      K[(k+j*m)*nm+(i*m+l)] = 0.0; 
	      for(h = 0; h < m; h++){
		K[(k+j*m)*nm+(i*m+l)] += A[k+m*h]*A[l+m*h]*spCorTS(obsD[j*n+i], phi[h], nu[h], covModel, &bessel_ws[threadID*bessel_ws_inc]);
	      }
	    }
	  }
	}
      }
      
#ifdef _OPENMP
#pragma omp parallel for private(i, k, l, h, threadID)
#endif
      for(j = 0; j < n; j++){
#ifdef _OPENMP
	threadID = omp_get_thread_num();
#endif	
	for(i = 0; i < q; i++){	
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      B[(k+j*m)*qm+(i*m+l)] = 0.0; 
	      for(h = 0; h < m; h++){
		B[(k+j*m)*qm+(i*m+l)] += A[k+m*h]*A[l+m*h]*spCorTS(predObsD[j*q+i], phi[h], nu[h], covModel, &bessel_ws[threadID*bessel_ws_inc]);
	      }
	    }
	  }
	}
      }
      
      //printMtrx(B, qm, nm);
      
#ifdef _OPENMP
#pragma omp parallel for private(i, k, l, h, threadID)
#endif
      for(j = 0; j < q; j++){
#ifdef _OPENMP
	threadID = omp_get_thread_num();
#endif
	for(i = 0; i < q; i++){	
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      C[(k+j*m)*qm+(i*m+l)] = 0.0; 
	      for(h = 0; h < m; h++){
		C[(k+j*m)*qm+(i*m+l)] += A[k+m*h]*A[l+m*h]*spCorTS(predD[j*q+i], phi[h], nu[h], covModel, &bessel_ws[threadID*bessel_ws_inc]);
	      }
	    }
	  }
	}
      }
         
      F77_NAME(dpotrf)(lower, &nm, K, &nm, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");}
      F77_NAME(dpotri)(lower, &nm, K, &nm, &info); if(info != 0){error("c++ error: dpotri failed\n");}     
      F77_NAME(dsymm)(rside, lower, &qm, &nm, &one, K, &nm, B, &qm, &zero, tmp_qmnm, &qm);
      
      //mu
      F77_NAME(dgemv)(ntran, &qm, &nm, &one, tmp_qmnm, &qm, &wSamples[s*nm], &incOne, &zero, tmp_qm, &incOne);

      //var
      F77_NAME(dgemm)(ntran, ytran, &qm, &qm, &nm, &one, tmp_qmnm, &qm, B, &qm, &zero, tmp_qmqm, &qm);

      for(i = 0; i < qmqm; i++){
	C[i] = C[i] - tmp_qmqm[i];
      }
            
      F77_NAME(dpotrf)(lower, &qm, C, &qm, &info); if(info != 0){error("c++ error: dpotrf failed 2\n");}
      
      mvrnorm(&REAL(wPredSamples_r)[s*qm], tmp_qm, C, qm, false);
      
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

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, wPredSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.w.predictive.samples")); 

    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);

    }
}
