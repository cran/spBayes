#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  SEXP spPPMvLMRecover(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP g_r, SEXP p_r, 
		       SEXP knotsD_r, SEXP knotsObsD_r, 
		       SEXP samples_r, SEXP beta_r,  SEXP nSamples_r, 
		       SEXP nugget_r, SEXP PsiDiag_r, SEXP covModel_r, 
		       SEXP modPP_r, SEXP verbose_r, SEXP nReport_r){

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
    int pp = p*p;
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int g = INTEGER(g_r)[0];
    int nLTr = m*(m-1)/2+m;
    int nm = n*m;
    int nmnm = nm*nm;
    int nmp = nm*p;
    int gm = g*m;
    int mm = m*m;
    int nmm = n*mm;
    int gmm = g*mm;
    int nmgm = nm*gm;
    int gmgm = gm*gm;
       
    double *knotsD = REAL(knotsD_r);
    double *knotsObsD = REAL(knotsObsD_r);
     
    double *samples = REAL(samples_r);
    double *beta = REAL(beta_r);
    int nSamples = INTEGER(nSamples_r)[0];
    
    bool modPP = static_cast<bool>(INTEGER(modPP_r)[0]);
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    bool PsiDiag = static_cast<bool>(INTEGER(PsiDiag_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));
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
    int status=1;
    SEXP wStrSamples_r; 
    SEXP wSamples_r;

    PROTECT(wStrSamples_r = allocMatrix(REALSXP, gm, nSamples)); nProtect++; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, nm, nSamples)); nProtect++; 
      
    double *P = (double *) R_alloc(nmgm, sizeof(double)); 
    double *K = (double *) R_alloc(gmgm, sizeof(double)); 
    double *D = (double *) R_alloc(nmm, sizeof(double)); 
    double *A = (double *) R_alloc(mm, sizeof(double)); 
    double *V = (double *) R_alloc(mm, sizeof(double)); 
    double *Psi = (double *) R_alloc(mm, sizeof(double)); zeros(Psi, mm); //must be cleared for diag Psi
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
    double *theta = (double *) R_alloc(2, sizeof(double)); //phi, nu, and perhaps more in the future
    
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nmgm  = (double *) R_alloc(nmgm, sizeof(double));
    double *tmp_nmgm2  = (double *) R_alloc(nmgm, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double)); 
    double *tmp_gm = (double *) R_alloc(gm, sizeof(double)); 
    double *tmp_gm2 = (double *) R_alloc(gm, sizeof(double)); 
    double *tmp_gmgm = (double *) R_alloc(gmgm, sizeof(double)); 
    
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tRecovering w\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){    
  
     covExpand(&samples[s*nParams+AIndx], A, m);//note this is K
      
      F77_NAME(dcopy)(&mm, A, &incOne, V, &incOne); //keep a copy of K

      F77_NAME(dpotrf)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");} //get A
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
     	    Psi[k*m+k] = samples[s*nParams+LIndx+k];//note, I use only the first column of Psi in other routines
     	  }
     	}else{
     	  covExpand(&samples[s*nParams+LIndx], Psi, m);//note, lower tri is Psi
     	}
      }
      
      for(jj = 0; jj < g; jj++){
     	for(ii = jj; ii < g; ii++){	
     	  for(k = 0; k < m; k++){
     	    for(l = 0; l < m; l++){
     	      K[(k+jj*m)*gm+(ii*m+l)] = 0.0; 
     	      for(h = 0; h < m; h++){
     		theta[0] = phi[h];
     		if(covModel == "matern"){
     		  theta[1] = nu[h];
     		}
     		K[(k+jj*m)*gm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(knotsD[jj*g+ii], theta, covModel);
     	      }
     	    }
     	  }
     	}
      }

      for(jj = 0; jj < n; jj++){
     	for(ii = 0; ii < g; ii++){	
     	  for(k = 0; k < m; k++){
     	    for(l = 0; l < m; l++){
     	      P[(k+jj*m)*gm+(ii*m+l)] = 0.0; 
     	      for(h = 0; h < m; h++){
     		theta[0] = phi[h];
     		if(covModel == "matern"){
     		  theta[1] = nu[h];
     		}
     		P[(k+jj*m)*gm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(knotsObsD[jj*g+ii], theta, covModel);
     	      }
     	    }
     	  }
     	}
      }

      F77_NAME(dpotrf)(lower, &gm, K, &gm, &info); if(info != 0){error("c++ error: dpotrf failed 2\n");}
      F77_NAME(dpotri)(lower, &gm, K, &gm, &info); if(info != 0){error("c++ error: dpotri failed\n");}
     
      //P K^{-1}
      F77_NAME(dsymm)(lside, lower, &gm, &nm, &one, K, &gm, P, &gm, &zero, tmp_nmgm, &gm);
      
      if(modPP){
	
     	for(k = 0; k < n; k++){
     	  F77_NAME(dgemm)(ytran, ntran, &m, &m, &gm, &one, &P[k*gm*m], &gm, &tmp_nmgm[k*gm*m], &gm, &zero, tmp_mm, &m);
	  
     	  for(l = 0; l < mm; l++){
     	    tmp_mm[l] = V[l] - tmp_mm[l];
     	    if(nugget){
     	      tmp_mm[l] += Psi[l];
     	    }
     	  }
	  
	  
	  F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed 3\n");}
	  F77_NAME(dpotri)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	  
	  F77_NAME(dcopy)(&mm, tmp_mm, &incOne, &D[k*mm], &incOne); //D^{-1}
	}
      }else{
	
	F77_NAME(dcopy)(&mm, Psi, &incOne, tmp_mm, &incOne);
	F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed 4\n");}
	F77_NAME(dpotri)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	
	for(k = 0; k < n; k++){
	  F77_NAME(dcopy)(&mm, tmp_mm, &incOne, &D[k*mm], &incOne); //D^{-1}
	}
      }
      
      for(k = 0; k < n; k++){
	F77_NAME(dsymm)(rside, lower, &gm, &m, &one, &D[k*mm], &m, &tmp_nmgm[(k*m)*gm], &gm, &zero, &tmp_nmgm2[(k*m)*gm], &gm);
      }
      
      //(C^{*-1} c) (1/E ct C^{*-1})
      F77_NAME(dgemm)(ntran, ytran, &gm, &gm, &nm, &one, tmp_nmgm2, &gm, tmp_nmgm, &gm, &zero, tmp_gmgm, &gm);
      
      for(i = 0; i < gmgm; i++){
	K[i] += tmp_gmgm[i];
      }
      
      //invert C_str
      F77_NAME(dpotrf)(lower, &gm, K, &gm, &info); if(info != 0){error("c++ error: dpotrf failed 5\n");}
      F77_NAME(dpotri)(lower, &gm, K, &gm, &info); if(info != 0){error("c++ error: dpotri failed\n");}
      
      F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, &beta[s], &nSamples, &zero, tmp_nm, &incOne);
      F77_NAME(daxpy)(&nm, &one, Y, &incOne, tmp_nm, &incOne);
      
      //(1/E ct C^{*-1})'(Y-XB)
      F77_NAME(dgemv)(ntran, &gm, &nm, &one, tmp_nmgm2, &gm, tmp_nm, &incOne, &zero, tmp_gm, &incOne);     
      F77_NAME(dsymv)(lower, &gm, &one, K, &gm, tmp_gm, &incOne, &zero, tmp_gm2, &incOne);
      
      //chol for the mvnorm and draw
      F77_NAME(dpotrf)(lower, &gm, K, &gm, &info); if(info != 0){error("c++ error: dpotrf failed 6\n");}
      mvrnorm(&REAL(wStrSamples_r)[s*gm], tmp_gm2, K, gm, false);
      
      //make \tild{w}
      F77_NAME(dgemv)(ytran, &gm, &nm, &one, tmp_nmgm, &gm, &REAL(wStrSamples_r)[s*gm], &incOne, &zero, &REAL(wSamples_r)[s*nm], &incOne);
      
      R_CheckUserInterrupt();
      
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
    int nResultListObjs = 2;
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, wSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.w.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, wStrSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.wStr.samples")); 
    
    namesgets(result_r, resultName_r);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
