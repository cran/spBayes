#include <string>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

extern "C" {

  SEXP spSVCRecover(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP Z_r, SEXP m_r, SEXP n_r, SEXP coords_r, SEXP q_r, SEXP KDiag_r,
		    SEXP samples_r, SEXP nSamples_r, 
		    SEXP AIndx_r, SEXP tauSqIndx_r, SEXP phiIndx_r, SEXP nuIndx_r, 
		    SEXP betaPrior_r, SEXP betaNorm_r, 	   
		    SEXP nugget_r, SEXP covModel_r, 
		    SEXP beta_r, SEXP w_r,
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
    double *Y = REAL(Y_r);
    double *X = REAL(X_r);
    double *Z = REAL(Z_r);
    int p = INTEGER(p_r)[0];
    int pp = p*p;
    int m = INTEGER(m_r)[0];
    int mm = m*m;
    int n = INTEGER(n_r)[0];
    int nn = n*n;
    int np = n*p;
    int nm = n*m;
    int nmnm = nm*nm;
    int q = INTEGER(q_r)[0];
    bool KDiag = static_cast<bool>(INTEGER(KDiag_r)[0]);
    int nLTr = m*(m-1)/2+m;
	
    double *coords = REAL(coords_r);

    double *samples = REAL(samples_r);
    int nSamples = INTEGER(nSamples_r)[0];

    int AIndx = INTEGER(AIndx_r)[0]; 
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
    int nThreads = INTEGER(nThreads_r)[0];
    
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    SEXP betaSamples_r, wSamples_r;
    if(getBeta){
      PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++;
    }
    
    if(getW){ 
      PROTECT(wSamples_r = allocMatrix(REALSXP, nm, nSamples)); nProtect++; 
    }

    int status=1;
    double *A = (double *) R_alloc(mm, sizeof(double)); zeros(A, mm); //to simplify a future move to the more general cross-cov model
    double *K = (double *) R_alloc(nmnm, sizeof(double)); zeros(K, nmnm);
    double *C = (double *) R_alloc(nn, sizeof(double)); zeros(C, nn);    
    double *B = (double *) R_alloc(pp, sizeof(double));
    double *bb = (double *) R_alloc(p, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_m = (double *) R_alloc(nThreads*m, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *tmp_nm2 = (double *) R_alloc(nm, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nltr = (double *) R_alloc(nLTr, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double)); zeros(nu, m); //this just remains empty of not matern
    double tauSq;
  
    double *D = (double *) R_alloc(nn, sizeof(double));
    distN(coords, n, coords, n, q, D);
    
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

    double *XSbetaX = NULL;
    double *tmp_np = NULL;
    
    if((betaPrior == "normal") & getW){
      XSbetaX = (double *) R_alloc(nn, sizeof(double));
      tmp_np = (double *) R_alloc(np, sizeof(double));

      F77_NAME(dgemm)(ntran, ntran, &n, &p, &p, &one, X, &n, betaC, &p, &zero, tmp_np, &n);
      F77_NAME(dgemm)(ntran, ytran, &n, &n, &p, &one, tmp_np, &n, X, &n, &zero, XSbetaX, &n);
    }

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
      if(getW){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\tRecovering beta and w\n");
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

      if(nugget){
	tauSq = samples[tauSqIndx*nSamples+s];
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
		K[(k+j*m)*nm+(i*m+l)] += A[k+m*h]*A[l+m*h]*spCorTS(D[j*n+i], phi[h], nu[h], covModel, &bessel_ws[threadID*bessel_ws_inc]);
	      }
	    }
	  }
	  
	  //mk\tild{X} C \tild{X}^T
	  for(k = 0; k < m; k++){
	    tmp_m[k+threadID*m] = F77_NAME(ddot)(&m, &Z[j], &n, &K[(k+j*m)*nm+(i*m)], &incOne); 
	  }
	  
	  C[j*n+i] = F77_NAME(ddot)(&m, &tmp_m[threadID*m], &incOne, &Z[i], &n);	      
	}
      }
 
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
      	  F77_NAME(dpotrf)(lower, &nm, K, &nm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      	  F77_NAME(dpotri)(lower, &nm, K, &nm, &info); if(info != 0){error("c++ error: dpotri failed\n");}

      	  for(i = 0; i < n; i++){
	    
      	    F77_NAME(dgemm)(ytran, ntran, &m, &m, &incOne, &one, &Z[i], &n, &Z[i], &n, &zero, tmp_mm, &m);
	    
      	    for(j = 0; j < mm; j++){
      	      tmp_mm[j] *= 1.0/tauSq;
      	    }
	    
      	    for(k = 0; k < m; k++){
      	      F77_NAME(daxpy)(&m, &one, &tmp_mm[k*m], &incOne, &K[i*m*nm+k*nm+i*m], &incOne);
      	    }
      	  }
	   	  
      	  F77_NAME(dpotrf)(lower, &nm, K, &nm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      	  F77_NAME(dpotri)(lower, &nm, K, &nm, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	  
      	  F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &REAL(betaSamples_r)[s*p], &incOne, &zero, vU, &incOne);
      	  F77_NAME(daxpy)(&n, &one, Y, &incOne, vU, &incOne);
	
      	  for(i = 0; i < n; i++){
      	    vU[i] *= 1.0/tauSq;
      	  }

      	  zeros(tmp_nm, nm);
      	  for(k = 0, i = 0; k < n; k++){
      	    for(l = 0; l < m; l++, i++){
      	      tmp_nm[i] = Z[l*n+k]*vU[k];
      	    }
      	  }
	  
      	  F77_NAME(dsymv)(lower, &nm, &one, K, &nm, tmp_nm, &incOne, &zero, tmp_nm2, &incOne);
	
      	  F77_NAME(dpotrf)(lower, &nm, K, &nm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      	  mvrnorm(&REAL(wSamples_r)[s*nm], tmp_nm2, K, nm, false);
	      
      	}else{
      	  //for now no w for no nugget model
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
