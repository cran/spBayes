#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  
  SEXP spMisalignRecover(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,
			 SEXP samples_r, SEXP nSamples_r, 
			 SEXP betaPrior_r, SEXP betaNorm_r, 	   
			 SEXP nugget_r, SEXP covModel_r, 
			 SEXP beta_r, SEXP w_r,
			 SEXP verbose_r, SEXP nReport_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, s, ii, jj, kk, info, nProtect=0;
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
    int *p = INTEGER(p_r);
    int *n = INTEGER(n_r);
    int m = INTEGER(m_r)[0];
    int nLTr = m*(m-1)/2+m;

    int N = 0;
    int P = 0;
    for(i = 0; i < m; i++){
      N += n[i];
      P += p[i];
    }
    
    int mm = m*m;
    int NN = N*N;
    int NP = N*P;
    int PP = P*P;

    double *coordsD = REAL(coordsD_r);

    double *samples = REAL(samples_r);
    int nSamples = INTEGER(nSamples_r)[0];
  
    //priors
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));
    double *betaMu = NULL;
    double *betaC = NULL;
    
    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(P, sizeof(double));
      F77_NAME(dcopy)(&P, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);
      
      betaC = (double *) R_alloc(PP, sizeof(double)); 
      F77_NAME(dcopy)(&PP, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaC, &incOne);
    }
  
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    bool getBeta = static_cast<bool>(INTEGER(beta_r)[0]);
    bool getW = static_cast<bool>(INTEGER(w_r)[0]);
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    int nParams, AIndx, PsiIndx, phiIndx, nuIndx;
    
    if(!nugget && covModel != "matern"){
      nParams = nLTr+m;//A, phi
      AIndx = 0; phiIndx = nLTr;
    }else if(nugget && covModel != "matern"){
      nParams = nLTr+m+m;//A, diag(Psi), phi
      AIndx = 0; PsiIndx = nLTr; phiIndx = PsiIndx+m;
    }else if(!nugget && covModel == "matern"){
      nParams = nLTr+2*m;//A, phi, nu
      AIndx = 0; phiIndx = nLTr, nuIndx = phiIndx+m;
    }else{
      nParams = nLTr+3*m;//A, diag(Psi), phi, nu
      AIndx = 0; PsiIndx = nLTr, phiIndx = PsiIndx+m, nuIndx = phiIndx+m;
     }
    
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    SEXP betaSamples_r, wSamples_r;
    if(getBeta){
      PROTECT(betaSamples_r = allocMatrix(REALSXP, P, nSamples)); nProtect++;
    }
    
    if(getW){ 
      PROTECT(wSamples_r = allocMatrix(REALSXP, N, nSamples)); nProtect++; 
    }
    
    int status=1;
    
    double *C = (double *) R_alloc(NN, sizeof(double));
    double *A = (double *) R_alloc(mm, sizeof(double)); 
    double *Psi = (double *) R_alloc(m, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));

    double *B = (double *) R_alloc(PP, sizeof(double));
    double *bb = (double *) R_alloc(P, sizeof(double));
    double *tmp_P = (double *) R_alloc(P, sizeof(double));
    double *tmp_PP = (double *) R_alloc(PP, sizeof(double));
    
    int P1 = P+1;
    double *vU = (double *) R_alloc(N*P1, sizeof(double));
    
    double *betaCInv = NULL;
    double *betaCInvMu = NULL;

    int sl, sk;

    if(betaPrior == "normal"){
      betaCInv = (double *) R_alloc(PP, sizeof(double));
      betaCInvMu = (double *) R_alloc(P, sizeof(double));
      
      F77_NAME(dcopy)(&PP, betaC, &incOne, betaCInv, &incOne);
      F77_NAME(dpotrf)(lower, &P, betaCInv, &P, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &P, betaCInv, &P, &info); if(info != 0){error("c++ error: dpotri failed\n");}
      
      F77_NAME(dsymv)(lower, &P, &one, betaCInv, &P, betaMu, &incOne, &zero, betaCInvMu, &incOne);      
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
      
      covExpand(&samples[s*nParams+AIndx], A, m);//note this is K, so we need chol
      F77_NAME(dpotrf)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");} 
      clearUT(A, m); //make sure upper tri is clear
   
      for(k = 0; k < m; k++){
	phi[k] = samples[s*nParams+(phiIndx+k)];
	
	if(covModel == "matern"){
	  nu[k] = samples[s*nParams+(nuIndx+k)];
	}

      }

      if(nugget){
	for(k = 0; k < m; k++){
	  Psi[k] = samples[s*nParams+(PsiIndx+k)];
	}
      }
      
      //construct covariance matrix
      sl = sk = 0;
      
      for(k = 0; k < m; k++){
	sl = 0;
	for(l = 0; l < m; l++){
	  for(kk = 0; kk < n[k]; kk++){
	    for(jj = 0; jj < n[l]; jj++){
	      C[(sl+jj)*N+(sk+kk)] = 0.0;
	      for(ii = 0; ii < m; ii++){
		C[(sl+jj)*N+(sk+kk)] += A[k+m*ii]*A[l+m*ii]*spCor(coordsD[(sl+jj)*N+(sk+kk)], phi[ii], nu[ii], covModel);
	      }
	    }
	  }
	  sl += n[l];
	}
	sk += n[k];
      }
      
      if(nugget){
	sl = 0;
	for(l = 0; l < m; l++){
	  for(k = 0; k < n[l]; k++){
	    C[(sl+k)*N+(sl+k)] += Psi[l];
	  }
	  sl += n[l];
	}
      }
            
      F77_NAME(dpotrf)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotrf failed 2\n");}
      
      F77_NAME(dcopy)(&N, Y, &incOne, vU, &incOne);
      F77_NAME(dcopy)(&NP, X, &incOne, &vU[N], &incOne);
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &N, &P1, &one, C, &N, vU, &N);//L[v:U] = [y:X]
      
      //B
      F77_NAME(dgemm)(ytran, ntran, &P, &P, &N, &one, &vU[N], &N, &vU[N], &N, &zero, B, &P); //U'U
      
      if(betaPrior == "normal"){
	for(k = 0; k < P; k++){
	  for(l = k; l < P; l++){
	    B[k*P+l] += betaCInv[k*P+l];
	  }
	}
      }
      
      F77_NAME(dpotrf)(lower, &P, B, &P, &info); if(info != 0){error("c++ error: dpotrf failed 3\n");}
      F77_NAME(dpotri)(lower, &P, B, &P, &info); if(info != 0){error("c++ error: dpotri failed 4\n");}
      
      //bb
      F77_NAME(dgemv)(ytran, &N, &P, &one, &vU[N], &N, vU, &incOne, &zero, tmp_P, &incOne); //U'v
      
      if(betaPrior == "normal"){
	for(k = 0; k < P; k++){
	  tmp_P[k] += betaCInvMu[k];
	}
      }
      
      F77_NAME(dsymv)(lower, &P, &one, B, &P, tmp_P, &incOne, &zero, bb, &incOne); 
      F77_NAME(dpotrf)(lower, &P, B, &P, &info); if(info != 0){error("c++ error: dpotrf failed 5\n");}
      
      mvrnorm(&REAL(betaSamples_r)[s*P], bb, B, P, false);
      
      //get w
      if(getW){
	
	if(nugget){
	  
	  //construct covariance matrix
	  sl = sk = 0;
	  
	  for(k = 0; k < m; k++){
	    sl = 0;
	    for(l = 0; l < m; l++){
	      for(kk = 0; kk < n[k]; kk++){
		for(jj = 0; jj < n[l]; jj++){
		  C[(sl+jj)*N+(sk+kk)] = 0.0;
		  for(ii = 0; ii < m; ii++){
		    C[(sl+jj)*N+(sk+kk)] += A[k+m*ii]*A[l+m*ii]*spCor(coordsD[(sl+jj)*N+(sk+kk)], phi[ii], nu[ii], covModel);
		  }
		}
	      }
	      sl += n[l];
	    }
	    sk += n[k];
	  }
	  
	  F77_NAME(dpotrf)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotrf failed 6\n");}
	  F77_NAME(dpotri)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotri failed 7\n");}
	  
	  

	  sl = 0;
	  for(l = 0; l < m; l++){
	    for(k = 0; k < n[l]; k++){
	      C[(sl+k)*N+(sl+k)] += 1.0/Psi[l];
	    }
	    sl += n[l];
	  }
	  
	  F77_NAME(dpotrf)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotrf failed 10\n");}
	  F77_NAME(dpotri)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotri failed 11\n");}

	  F77_NAME(dgemv)(ntran, &N, &P, &negOne, X, &N, &REAL(betaSamples_r)[s*P], &incOne, &zero, vU, &incOne);
	  F77_NAME(daxpy)(&N, &one, Y, &incOne, vU, &incOne);

	  sl = 0;
	  for(l = 0; l < m; l++){
	    for(k = 0; k < n[l]; k++){
	     vU[N+sl+k] = vU[sl+k]*1.0/Psi[l];
	    }
	    sl += n[l];
	  }
	  
	  F77_NAME(dsymv)(lower, &N, &one, C, &N, &vU[N], &incOne, &zero, vU, &incOne);
	  
	  F77_NAME(dpotrf)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotrf failed 12\n");}
	  
	  mvrnorm(&REAL(wSamples_r)[s*N], vU, C, N, false);
	  
	}else{
	  F77_NAME(dgemv)(ntran, &N, &P, &negOne, X, &N, &REAL(betaSamples_r)[s*P], &incOne, &zero, &REAL(wSamples_r)[s*N], &incOne);
	  F77_NAME(daxpy)(&N, &one, Y, &incOne, &REAL(wSamples_r)[s*N], &incOne);
	}
	
	
      }

      R_CheckUserInterrupt();

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
    }
    
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
