#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  SEXP spPPMvLMPredict(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP g_r, SEXP p_r, SEXP Z_r, SEXP q_r, 
		     SEXP knotsD_r, SEXP knotsObsD_r, SEXP knotsPredD_r, 
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
    
    double *Z = REAL(Z_r);
    int q = INTEGER(q_r)[0];
    int qm = q*m;
    int qmqm = qm*qm;
    int qmp = qm*p;
    int nmqm = nm*qm;
    int mm = m*m;
    int nmm = n*mm;
    int gmm = g*mm;
    int qmgm = qm*gm;
    int nmgm = nm*gm;
    int gmgm = gm*gm;
       

    double *knotsD = REAL(knotsD_r);
    double *knotsObsD = REAL(knotsObsD_r);
    double *knotsPredD = REAL(knotsPredD_r);
    
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
  
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Prediction at %i locations.\n\n", q);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      if(!nugget){
	Rprintf("Psi not included in the model (i.e., no nugget model).\n\n");
      }
      
    } 
    
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status=1;
    SEXP predSamples_r;
    PROTECT(predSamples_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 
    
    double *P = (double *) R_alloc(nmgm, sizeof(double)); 
    double *K = (double *) R_alloc(gmgm, sizeof(double)); 
    double *D = (double *) R_alloc(nmm, sizeof(double)); 
    double *O = (double *) R_alloc(qmgm, sizeof(double));
    double *H = (double *) R_alloc(nmgm, sizeof(double)); 
    double *v = (double *) R_alloc(nm, sizeof(double)); 
    double *L = (double *) R_alloc(gmgm, sizeof(double));
    double *u = (double *) R_alloc(nmm, sizeof(double)); 
    double *z = (double *) R_alloc(qm, sizeof(double)); 
 
    double *A = (double *) R_alloc(mm, sizeof(double)); 
    double *V = (double *) R_alloc(mm, sizeof(double)); 
    double *Psi = (double *) R_alloc(mm, sizeof(double)); zeros(Psi, mm); //must be cleared for diag Psi
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
    double *theta = (double *) R_alloc(2, sizeof(double)); //phi, nu, and perhaps more in the future
    
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double)); 
    double *tmp_mm2 = (double *) R_alloc(mm, sizeof(double)); 
    double *tmp_gmm = (double *) R_alloc(gmm, sizeof(double)); 
    double *tmp_gm = (double *) R_alloc(gm, sizeof(double)); 
    double *tmp_gmgm = (double *) R_alloc(gmgm, sizeof(double)); 
    double *tmp_m = (double *) R_alloc(m, sizeof(double)); 
    double *tmp_m2 = (double *) R_alloc(m, sizeof(double)); 
    
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
      
      F77_NAME(dgemv)(ntran, &qm, &p, &one, Z, &qm, &beta[s], &nSamples, &zero, z, &incOne);

      F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, &beta[s], &nSamples, &zero, v, &incOne);
      F77_NAME(daxpy)(&nm, &one, Y, &incOne, v, &incOne);
   
      covExpand(&samples[s*nParams+AIndx], A, m);//note this is K
      
      F77_NAME(dcopy)(&mm, A, &incOne, V, &incOne); //keep a copy of K

      F77_NAME(dpotrf)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotrf failed\n");} //get A
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

      for(jj = 0; jj < q; jj++){
     	for(ii = 0; ii < g; ii++){	
     	  for(k = 0; k < m; k++){
     	    for(l = 0; l < m; l++){
     	      O[(k+jj*m)*gm+(ii*m+l)] = 0.0; 
     	      for(h = 0; h < m; h++){
     		theta[0] = phi[h];
     		if(covModel == "matern"){
     		  theta[1] = nu[h];
     		}
     		O[(k+jj*m)*gm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(knotsPredD[jj*g+ii], theta, covModel);
     	      }
     	    }
     	  }
     	}
      }

     F77_NAME(dcopy)(&nmgm, P, &incOne, H, &incOne);

     //get a copy of K
     for(k = 0; k < gm; k++){
       for(l = k; l < gm; l++){
     	 L[k*gm+l] = K[k*gm+l];
       }
     }

      F77_NAME(dpotrf)(lower, &gm, K, &gm, &info);//L_K
      
      if(modPP){
     	F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &gm, &nm, &one, K, &gm, H, &gm);
	
     	for(k = 0; k < n; k++){
     	  F77_NAME(dgemm)(ytran, ntran, &m, &m, &gm, &one, &H[k*gm*m], &gm, &H[k*gm*m], &gm, &zero, tmp_mm, &m);
	  
     	  for(l = 0; l < mm; l++){
     	    tmp_mm[l] = V[l] - tmp_mm[l];
     	    if(nugget){
     	      tmp_mm[l] += Psi[l];
     	    }
     	  }
	  
     	  F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");}
     	  F77_NAME(dtrtri)(lower, nUnit, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dtrtri failed\n");}

     	  F77_NAME(dcopy)(&mm, tmp_mm, &incOne, &D[k*mm], &incOne); //D^{-1/2}
     	}
	
      }else{

     	F77_NAME(dcopy)(&mm, Psi, &incOne, tmp_mm, &incOne);
     	F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed 2\n");}
     	F77_NAME(dtrtri)(lower, nUnit, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dtrtri failed\n");}
	
     	for(k = 0; k < n; k++){
     	  F77_NAME(dcopy)(&mm, tmp_mm, &incOne, &D[k*mm], &incOne); //D^{-1/2}
     	}
      }

      F77_NAME(dcopy)(&nmgm, P, &incOne, H, &incOne);//get clean copy of P
      for(k = 0; k < n; k++){
     	F77_NAME(dtrmm)(rside, lower, ytran, nUnit, &gm, &m, &one, &D[k*mm], &m, &H[k*m*gm], &gm);//W'
      }

      F77_NAME(dgemm)(ntran, ytran, &gm, &gm, &nm, &one, H, &gm, H, &gm, &zero, tmp_gmgm, &gm);//W'W
      
      for(k = 0; k < gm; k++){
      	for(l = k; l < gm; l++){
      	  L[k*gm+l] += tmp_gmgm[k*gm+l];
      	}
      }

      F77_NAME(dpotrf)(lower, &gm, L, &gm, &info); if(info != 0){error("c++ error: dpotrf failed 3\n");}//L

      //get H      
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &gm, &nm, &one, L, &gm, H, &gm);//LH = W'

      //get [v_1:v_2: ... : v_n]
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &gm, &nm, &one, K, &gm, P, &gm);//L_K [v_1:v_2: ... : v_n] = [c_1:c_2: ... : c_n]

      //get v = D^{-1/2}(y-XB)      
      for(k = 0; k < n; k++){
	F77_NAME(dtrmv)(lower, ntran, nUnit, &m, &D[k*mm], &m, &v[k*m], &incOne);
      }
      
      //for each location
      for(j = 0; j < q; j++){
	
	F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &gm, &m, &one, K, &gm, &O[(j*m)*gm], &gm);//L_K u = c_0

	//get u = D^{-1/2}\tidle{c_0}
      	for(k = 0; k < n; k++){
      	  F77_NAME(dgemm)(ytran, ntran, &m, &m, &gm, &one, &O[(j*m)*gm], &gm, &P[(k*m)*gm], &gm, &zero, &u[k*mm], &m);//u'v_1 ... u'v_n
	  F77_NAME(dtrmm)(rside, lower, ytran, nUnit, &m, &m, &one, &D[k*mm], &m, &u[k*mm], &m);
	}
   
	//get w = Hu
	F77_NAME(dgemm)(ntran, ytran, &gm, &m, &nm, &one, H, &gm, u, &m, &zero, tmp_gmm, &gm);

	//get \tilde(c_0)'\Sigma^{-1}\tilde(c_0)
	F77_NAME(dgemm)(ntran, ytran, &m, &m, &nm, &one, u, &m, u, &m, &zero, tmp_mm, &m);

	F77_NAME(dgemm)(ytran, ntran, &m, &m, &gm, &one, tmp_gmm, &gm, tmp_gmm, &gm, &zero, tmp_mm2, &m);

	for(i = 0; i < mm; i++){
	  tmp_mm[i] = V[i] - (tmp_mm[i] - tmp_mm2[i]);
	  if(nugget){
	    tmp_mm[i] += Psi[i];
	  }
	}

	F77_NAME(dgemv)(ntran, &gm, &nm, &one, H, &gm, v, &incOne, &zero, tmp_gm, &incOne); //z = Hv

	F77_NAME(dgemv)(ntran, &m, &nm, &one, u, &m, v, &incOne, &zero, tmp_m, &incOne); //u'v

	F77_NAME(dgemv)(ytran, &gm, &m, &one, tmp_gmm, &gm, tmp_gm, &incOne, &zero, tmp_m2, &incOne); //w'z

	for(i = 0; i < m; i++){
	  tmp_m[i] =  z[j*m+i] + (tmp_m[i] - tmp_m2[i]);
	}

	F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed 4\n");}   
	mvrnorm(&REAL(predSamples_r)[s*qm+j*m], tmp_m, tmp_mm, m, false);

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
