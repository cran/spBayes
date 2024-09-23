#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

#define USE_FC_LEN_T
#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include "util.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

extern "C" {

  SEXP spSVCPredictMarg(SEXP m_r, SEXP n_r, SEXP KDiag_r, SEXP obsD_r, SEXP predObsD_r, SEXP q_r,
			SEXP samples_r, SEXP wSamples_r, SEXP nSamples_r, 
			SEXP AIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 	   
			SEXP covModel_r, 
			SEXP verbose_r, SEXP nReport_r, SEXP nThreads_r){

    /*****************************************
                Common variables
    *****************************************/
    int i, j, k, l, s, h, info, nProtect=0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *lside = "L";
    const double one = 1.0;
    const double zero = 0.0;
    const int incOne = 1;
    
    /*****************************************
                     Set-up
    *****************************************/
    double *obsD = REAL(obsD_r);
    double *predObsD = REAL(predObsD_r);
    int m = INTEGER(m_r)[0];
    int mm = m*m;
    int n = INTEGER(n_r)[0];
    int nm = n*m;
    int nmnm = nm*nm;
    int q = INTEGER(q_r)[0];//number of prediction locations
    int qm = q*m;
    int qmnm = qm*nm;
    //int qmqm = qm*qm;
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
    PROTECT(wPredSamples_r = Rf_allocMatrix(REALSXP, qm, nSamples)); nProtect++; 

    int status=1;
    double *A = (double *) R_alloc(mm, sizeof(double)); zeros(A, mm); //to simplify a future move to the more general cross-cov model
    double *K = (double *) R_alloc(nmnm, sizeof(double)); 
    double *B = (double *) R_alloc(qmnm, sizeof(double)); //if q is large then there will likley be mem issues, I can revise later, e.g., don't store predObsD or B
    double *C = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nltr = (double *) R_alloc(nLTr, sizeof(double)); 
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nmm = (double *) R_alloc(nm*m, sizeof(double));
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
      Rf_warning("n.omp.threads = %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif  
    
    if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\tPoint-wise sampling of predicted w\n");
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
      	F77_NAME(dpotrf)(lower, &m, A, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed 1\n");} 
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

      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, A, &m, A, &m, &zero, C, &m FCONE FCONE);
            
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
	      B[(i*m+l)*nm+(k+j*m)] = 0.0; 
	      for(h = 0; h < m; h++){
		B[(i*m+l)*nm+(k+j*m)] += A[k+m*h]*A[l+m*h]*spCorTS(predObsD[j*q+i], phi[h], nu[h], covModel, &bessel_ws[threadID*bessel_ws_inc]);//note the transpose for B differs from the joint code
	      }
	    }
	  }
	}
      }

      //printMtrx(B, nm, qm);
           
      F77_NAME(dpotrf)(lower, &nm, K, &nm, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed 1\n");}
      F77_NAME(dpotri)(lower, &nm, K, &nm, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotri failed\n");}     
 
      for(i = 0; i < q; i++){

	//mu
	F77_NAME(dsymm)(lside, lower, &nm, &m, &one, K, &nm, &B[(i*m)*nm], &nm, &zero, tmp_nmm, &nm FCONE FCONE);
	F77_NAME(dgemv)(ytran, &nm, &m, &one, tmp_nmm, &nm, &wSamples[s*nm], &incOne, &zero, tmp_m, &incOne FCONE);
	
	//var
	F77_NAME(dgemm)(ytran, ntran, &m, &m, &nm, &one, tmp_nmm, &nm, &B[(i*m)*nm], &nm, &zero, tmp_mm, &m FCONE FCONE);		

	for(j = 0; j < mm; j++){
	  tmp_mm[j] = C[j] - tmp_mm[j];
	}
	
      	F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed 2\n");}
      	mvrnorm(&REAL(wPredSamples_r)[s*qm+i*m], tmp_m, tmp_mm, m, false);
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

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, wPredSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("p.w.predictive.samples")); 

    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);

    }
}
