#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  
  SEXP spMvLMRecover(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,
		     SEXP samples_r, SEXP nSamples_r, 
		     SEXP betaPrior_r, SEXP betaNorm_r, 	   
		     SEXP nugget_r, SEXP PsiDiag_r, SEXP covModel_r, 
		     SEXP beta_r, SEXP w_r,
		     SEXP verbose_r, SEXP nReport_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, s, ii, jj, info, nProtect=0;
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
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int nLTr = m*(m-1)/2+m;
    int nn = n*n;
    int mm = m*m;
    int nm = n*m;
    int nmnm = nm*nm;
    int nmp = nm*p;
    int pp = p*p;
    
    double *coordsD = REAL(coordsD_r);

    double *samples = REAL(samples_r);
    int nSamples = INTEGER(nSamples_r)[0];
  
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
    bool PsiDiag = static_cast<bool>(INTEGER(PsiDiag_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
    bool getBeta = static_cast<bool>(INTEGER(beta_r)[0]);
    bool getW = static_cast<bool>(INTEGER(w_r)[0]);
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    int nParams, AIndx, LIndx, phiIndx, nuIndx;
    
    if(!nugget && covModel != "matern"){
      nParams = nLTr+m;//A, phi
      AIndx = 0; phiIndx = nLTr;
    }else if(nugget && covModel != "matern"){
      if(PsiDiag){
	nParams = nLTr+m+m;//A, diag(Psi), phi
	AIndx = 0; LIndx = nLTr; phiIndx = LIndx+m;
      }else{
	nParams = 2*nLTr+m;//A, L, phi
	AIndx = 0; LIndx = nLTr; phiIndx = LIndx+nLTr;
      }
    }else if(!nugget && covModel == "matern"){
      nParams = nLTr+2*m;//A, phi, nu
      AIndx = 0; phiIndx = nLTr, nuIndx = phiIndx+m;
    }else{
      if(PsiDiag){
	nParams = nLTr+3*m;//A, diag(Psi), phi, nu
	AIndx = 0; LIndx = nLTr, phiIndx = LIndx+m, nuIndx = phiIndx+m;
      }else{
	nParams = 2*nLTr+2*m;//A, Psi, phi, nu
	AIndx = 0; LIndx = nLTr, phiIndx = LIndx+nLTr, nuIndx = phiIndx+m;
      }
    }
    
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
    
    double *C = (double *) R_alloc(nmnm, sizeof(double));
    double *W = (double *) R_alloc(nmnm, sizeof(double));
    double *A = (double *) R_alloc(mm, sizeof(double)); 
    double *L = (double *) R_alloc(mm, sizeof(double));
    double *Psi = (double *) R_alloc(mm, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
    double *theta = (double *) R_alloc(2, sizeof(double)); //phi, nu, and perhaps more in the future

    double *B = (double *) R_alloc(pp, sizeof(double));
    double *b = (double *) R_alloc(p, sizeof(double));
    double *bb = (double *) R_alloc(p, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm2 = (double *) R_alloc(mm, sizeof(double));
    //double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    
    int p1 = p+1;
    double *vU = (double *) R_alloc(nm*p1, sizeof(double));
    
    double *betaCInv = NULL;
    double *betaCInvMu = NULL;

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
	if(PsiDiag){
	  for(k = 0; k < m; k++){
	    Psi[k] = samples[s*nParams+(LIndx+k)];//first column of Psi holds the m tau.sq's
	  }
	}else{
	  covExpand(&samples[s*nParams+LIndx], Psi, m);
	}
      }

      
      //construct covariance matrix
      for(jj = 0; jj < n; jj++){
	for(ii = jj; ii < n; ii++){	
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      C[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
	      for(h = 0; h < m; h++){
		theta[0] = phi[h];
		if(covModel == "matern"){
		  theta[1] = nu[h];
		}
		C[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(coordsD[jj*n+ii], theta, covModel);
	      }
	    }
	  }
	}
      }

      if(nugget){
	if(PsiDiag){
	  for(l = 0; l < n; l++){
	    for(k = 0; k < m; k++){
	      C[(l*m+k)*nm+(l*m+k)] += Psi[k];
	    }
	  }
	}else{
	  for(i = 0; i < n; i++){
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
		C[(i*m+l)*nm+(i*m+k)] += Psi[l*m+k];
	      }
	    }
	  }
	}
      }
      
      F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotrf failed 2\n");}
      
      F77_NAME(dcopy)(&nm, Y, &incOne, vU, &incOne);
      F77_NAME(dcopy)(&nmp, X, &incOne, &vU[nm], &incOne);
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &nm, &p1, &one, C, &nm, vU, &nm);//L[v:U] = [y:X]
      
      //B
      F77_NAME(dgemm)(ytran, ntran, &p, &p, &nm, &one, &vU[nm], &nm, &vU[nm], &nm, &zero, B, &p); //U'U
      
      if(betaPrior == "normal"){
	for(k = 0; k < p; k++){
	  for(l = k; l < p; l++){
	    B[k*p+l] += betaCInv[k*p+l];
	  }
	}
      }
      
      F77_NAME(dpotrf)(lower, &p, B, &p, &info); if(info != 0){error("c++ error: dpotrf failed 3\n");}
      F77_NAME(dpotri)(lower, &p, B, &p, &info); if(info != 0){error("c++ error: dpotri failed 4\n");}
      
      //bb
      F77_NAME(dgemv)(ytran, &nm, &p, &one, &vU[nm], &nm, vU, &incOne, &zero, tmp_p, &incOne); //U'v
      
      if(betaPrior == "normal"){
	for(k = 0; k < p; k++){
	  tmp_p[k] += betaCInvMu[k];
	}
      }
      
      F77_NAME(dsymv)(lower, &p, &one, B, &p, tmp_p, &incOne, &zero, bb, &incOne); 
      F77_NAME(dpotrf)(lower, &p, B, &p, &info); if(info != 0){error("c++ error: dpotrf failed 5\n");}
      
      mvrnorm(&REAL(betaSamples_r)[s*p], bb, B, p, false);
      
      //get w
      if(getW){
	
	if(nugget){

	  //v2
	  //construct covariance matrix
	  for(jj = 0; jj < n; jj++){
	    for(ii = jj; ii < n; ii++){	
	      for(k = 0; k < m; k++){
	  	for(l = 0; l < m; l++){
	  	  C[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
	  	  for(h = 0; h < m; h++){
	  	    theta[0] = phi[h];
	  	    if(covModel == "matern"){
	  	      theta[1] = nu[h];
	  	    }
	  	    C[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(coordsD[jj*n+ii], theta, covModel);
	  	  }
	  	}
	      }
	    }
	  }

	  zeros(W, nmnm);

	  if(PsiDiag){
	    for(l = 0; l < n; l++){
	      for(k = 0; k < m; k++){
	  	C[(l*m+k)*nm+(l*m+k)] += Psi[k];
	  	W[(l*m+k)*nm+(l*m+k)] = Psi[k];
	      }
	    }
	  }else{
	    //fill in the upper tri of Psi for W
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
	  	Psi[l*m+k] =  Psi[k*m+l];
	      }
	    }
 
	    for(i = 0; i < n; i++){
	      for(k = 0; k < m; k++){
	  	for(l = 0; l < m; l++){
	  	  C[(i*m+l)*nm+(i*m+k)] += Psi[l*m+k];
	  	  W[(i*m+l)*nm+(i*m+k)] = Psi[l*m+k];
	  	}
	      }
	    }
	  }

	  //L
	  F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: 1 dpotrf failed\n");}

	  //W
	  F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &nm, &nm, &one, C, &nm, W, &nm);

	  //L_B
	  F77_NAME(dgemm)(ytran, ntran, &nm, &nm, &nm, &negOne, W, &nm, W, &nm, &zero, C, &nm);

	  if(PsiDiag){
	    for(l = 0; l < n; l++){
	      for(k = 0; k < m; k++){
	  	C[(l*m+k)*nm+(l*m+k)] += Psi[k];
	      }
	    }
	  }else{	    
	    for(i = 0; i < n; i++){
	      for(k = 0; k < m; k++){
	  	for(l = 0; l < m; l++){
	  	  C[(i*m+l)*nm+(i*m+k)] += Psi[l*m+k];
	  	}
	      }
	    }
	  }

	  F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: 2 dpotrf failed\n");}

	  F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, &REAL(betaSamples_r)[s*p], &incOne, &zero, vU, &incOne);
	  F77_NAME(daxpy)(&nm, &one, Y, &incOne, vU, &incOne);

	  if(!PsiDiag){
	    F77_NAME(dpotrf)(lower, &m, Psi, &m, &info); if(info != 0){error("c++ error: 3 dpotrf failed 8\n");}
	    F77_NAME(dpotri)(lower, &m, Psi, &m, &info); if(info != 0){error("c++ error: 4 dpotri failed 9\n");}
	  }

	  for(k = 0; k < n; k++){
	    if(PsiDiag){
	      for(l = 0; l < m; l++){
	  	vU[nm+k*m+l] = vU[k*m+l]/Psi[l];
	      }
	    }else{
	      F77_NAME(dsymv)(lower, &m, &one, Psi, &m, &vU[k*m], &incOne, &zero, &vU[nm+k*m], &incOne); 
	    }
	  }

	  F77_NAME(dtrmv)(lower, ytran, nUnit, &nm, C, &nm, &vU[nm], &incOne);
	  F77_NAME(dtrmv)(lower, ntran, nUnit, &nm, C, &nm, &vU[nm], &incOne);

	  for(k = 0; k < nm; k++){
	    vU[k] = rnorm(0, 1);
	  }

	  F77_NAME(dtrmv)(lower, ntran, nUnit, &nm, C, &nm, vU, &incOne);

	  for(k = 0; k < nm; k++){
	    REAL(wSamples_r)[s*nm+k] = vU[nm+k] + vU[k];
	  }

	  // //v1	  
	  // //construct covariance matrix
	  // for(jj = 0; jj < n; jj++){
	  //   for(ii = jj; ii < n; ii++){	
	  //     for(k = 0; k < m; k++){
	  // 	for(l = 0; l < m; l++){
	  // 	  C[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
	  // 	  for(h = 0; h < m; h++){
	  // 	    theta[0] = phi[h];
	  // 	    if(covModel == "matern"){
	  // 	      theta[1] = nu[h];
	  // 	    }
	  // 	    C[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(coordsD[jj*n+ii], theta, covModel);
	  // 	  }
	  // 	}
	  //     }
	  //   }
	  // }
	  
	  // F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotrf failed 6\n");}
	  // F77_NAME(dpotri)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotri failed 7\n");}
	  
	  // if(PsiDiag){
	  //   for(l = 0; l < n; l++){
	  //     for(k = 0; k < m; k++){
	  // 	C[(l*m+k)*nm+(l*m+k)] += 1.0/Psi[k];
	  //     }
	  //   }
	  // }else{
	  //   F77_NAME(dpotrf)(lower, &m, Psi, &m, &info); if(info != 0){error("c++ error: dpotrf failed 8\n");}
	  //   F77_NAME(dpotri)(lower, &m, Psi, &m, &info); if(info != 0){error("c++ error: dpotri failed 9\n");}
	    
	  //   for(i = 0; i < n; i++){
	  //     for(k = 0; k < m; k++){
	  // 	for(l = 0; l < m; l++){
	  // 	  C[(i*m+l)*nm+(i*m+k)] += Psi[l*m+k];
	  // 	}
	  //     }
	  //   }
	  // }
	  
	  // F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotrf failed 10\n");}
	  // F77_NAME(dpotri)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotri failed 11\n");}
	  
	  // F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, &REAL(betaSamples_r)[s*p], &incOne, &zero, vU, &incOne);
	  // F77_NAME(daxpy)(&nm, &one, Y, &incOne, vU, &incOne);
	  
	  // for(k = 0; k < n; k++){
	  //   if(PsiDiag){
	  //     for(l = 0; l < m; l++){
	  // 	vU[nm+k*m+l] = vU[k*m+l]*1.0/Psi[l];
	  //     }
	  //   }else{
	  //     F77_NAME(dsymv)(lower, &m, &one, Psi, &m, &vU[k*m], &incOne, &zero, &vU[nm+k*m], &incOne); 
	  //   }
	  // }
	  
	  // F77_NAME(dsymv)(lower, &nm, &one, C, &nm, &vU[nm], &incOne, &zero, vU, &incOne);
	  
	  // F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotrf failed 12\n");}
	  
	  // mvrnorm(&REAL(wSamples_r)[s*nm], vU, C, nm, false);
	  
	}else{
	  F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, &REAL(betaSamples_r)[s*p], &incOne, &zero, &REAL(wSamples_r)[s*nm], &incOne);
	  F77_NAME(daxpy)(&nm, &one, Y, &incOne, &REAL(wSamples_r)[s*nm], &incOne);
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
