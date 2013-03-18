#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  SEXP spMvLMPredict(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP p_r, SEXP Z_r, SEXP q_r, SEXP obsD_r, SEXP obsPredD_r, 
		     SEXP samples_r, SEXP beta_r,  SEXP nSamples_r, 
		     SEXP betaPrior_r, SEXP betaNorm_r, 
		     SEXP nugget_r, SEXP PsiDiag_r, SEXP covModel_r, 
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
    int pp = p*p;
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int nLTr = m*(m-1)/2+m;
    int nm = n*m;
    int nmnm = nm*nm;
    int nmp = nm*p;

    double *Z = REAL(Z_r);
    int q = INTEGER(q_r)[0];
    int qm = q*m;
    int qmqm = qm*qm;
    int qmp = qm*p;
    int nmqm = nm*qm;
    int mm = m*m;
    int qmm = q*mm; 
    int mp = m*p; 
     
    double *obsD = REAL(obsD_r);
    double *obsPredD = REAL(obsPredD_r);
    
    double *samples = REAL(samples_r);
    double *beta = REAL(beta_r);
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
    SEXP predSamples_r;
    PROTECT(predSamples_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 
    
    int status=1;
  
    double *C = (double *) R_alloc(nmnm, sizeof(double)); 
    double *c = (double *) R_alloc(nmqm, sizeof(double)); 
    double *z = (double *) R_alloc(nm, sizeof(double)); 
    double *u = (double *) R_alloc(nm, sizeof(double)); 
   
    double *A = (double *) R_alloc(mm, sizeof(double)); 
    double *K = (double *) R_alloc(mm, sizeof(double)); 
    double *L = (double *) R_alloc(mm, sizeof(double));
    double *Psi = (double *) R_alloc(mm, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
    double *theta = (double *) R_alloc(2, sizeof(double)); //phi, nu, and perhaps more in the future

    double *tmp_p = (double *) R_alloc(p, sizeof(double)); 
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double)); 
    double *tmp_qm = (double *) R_alloc(qm, sizeof(double)); 
    double *tmp_qmm = NULL;
    double *tmp_nmp = NULL;
    double *tmp_qmp = NULL;
    double *tmp_nmnm = NULL;
    double *tmp_nmqm = NULL;
    double *tmp_m = (double *) R_alloc(m, sizeof(double)); 
     
    if(betaPrior == "normal"){
      tmp_qmm = (double *) R_alloc(qmm, sizeof(double));
      tmp_nmp = (double *) R_alloc(nmp, sizeof(double));
      tmp_qmp = (double *) R_alloc(qmp, sizeof(double));
      tmp_nmnm = (double *) R_alloc(nmnm, sizeof(double));
      tmp_nmqm = (double *) R_alloc(nmqm, sizeof(double));

      //I really need to do away with the zeros in Z and X
      F77_NAME(dgemv)(ytran, &p, &qm, &one, Z, &p, betaMu, &incOne, &zero, tmp_qm, &incOne);

      F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, betaMu, &incOne, &zero, z, &incOne);
      F77_NAME(daxpy)(&nm, &one, Y, &incOne, z, &incOne);
       
      F77_NAME(dsymm)(rside, lower, &nm, &p, &one, betaC, &p, X, &nm, &zero, tmp_nmp, &nm);
      F77_NAME(dgemm)(ntran, ytran, &nm, &nm, &p, &one, tmp_nmp, &nm, X, &nm, &zero, tmp_nmnm, &nm);
 
      F77_NAME(dsymm)(lside, lower, &p, &qm, &one, betaC, &p, Z, &p, &zero, tmp_qmp, &p);

      for(i = 0; i < q; i++){
	F77_NAME(dgemm)(ytran, ntran, &m, &m, &p, &one, &Z[i*mp], &p, &tmp_qmp[i*mp], &p, &zero, &tmp_qmm[i*mm], &m);
      }

      F77_NAME(dgemm)(ytran, ytran, &qm, &nm, &p, &one, tmp_qmp, &p, X, &nm, &zero, tmp_nmqm, &qm); 
    }

    //check if any prediction locations are observed
    int *fitted = (int *) R_alloc(q, sizeof(int));
    for(i = 0; i < q; i++){
      fitted[i] = 0;
      for(j = 0; j < n; j++){
	if(obsPredD[i*n+j] == 0){
	  fitted[i] = 1;
	  break;
	}
      }
    }

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
      
      covExpand(&samples[s*nParams+AIndx], A, m);//note this is K

      F77_NAME(dcopy)(&mm, A, &incOne, K, &incOne); //keep a copy of K

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
		C[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(obsD[jj*n+ii], theta, covModel);
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
      	  for(l = 0; l < n; l++){
      	    for(k = 0; k < m; k++){
      	      F77_NAME(daxpy)(&m, &one, &Psi[k*m], &incOne, &C[l*m*nm+k*nm+l*m], &incOne);
      	    }
      	  }
      	}
      }
      
      //c nmxqm
      for(l = 0; l < q; l++){
      	for(h = 0; h < n; h++){
	  
      	  zeros(tmp_mm, mm);
	  
      	  for(k = 0; k < m; k++){
      	    theta[0] = phi[k];
	    
      	    if(covModel == "matern"){
      	      theta[1] = nu[k];
      	    }
	    
      	    spCor(&obsPredD[l*n+h], 1, theta, covModel, &tmp_mm[k*m+k]);
      	  }
	  
      	  F77_NAME(dtrmm)(lside, lower, ntran, nUnit, &m, &m, &one, A, &m, tmp_mm, &m);
      	  F77_NAME(dtrmm)(rside, lower, ytran, nUnit, &m, &m, &one, A, &m, tmp_mm, &m);
	  
      	  if(nugget && obsPredD[l*n+h] == 0){
      	    for(k = 0; k < m; k++){
      	      if(PsiDiag){
      		tmp_mm[k*m+k] += Psi[k];
      	      }else{
      		F77_NAME(daxpy)(&m, &one, &Psi[k*m], &incOne, &tmp_mm[k*m], &incOne);
      	      }
      	    }
      	  }
	  
      	  for(k = 0; k < m; k++){
      	    F77_NAME(dcopy)(&m, &tmp_mm[k*m], &incOne, &c[l*m*nm+k*nm+h*m], &incOne);
      	  }
	  
      	}
      }

      if(betaPrior == "normal"){
	
      	for(k = 0; k < nm; k++){
      	  for(l = k; l < nm; l++){
      	    C[k*nm+l] += tmp_nmnm[k*nm+l];
      	  }
      	}
	
      	for(k = 0; k < nm; k++){
      	  for(l = 0; l < qm; l++){
      	    c[l*nm+k] += tmp_nmqm[k*qm+l];
      	  }
      	}
      }

      F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}//L_1
      
      if(betaPrior == "normal"){
      	F77_NAME(dcopy)(&nm, z, &incOne, u, &incOne);
      }else{
      	F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, &beta[s], &nSamples, &zero, u, &incOne);
      	F77_NAME(daxpy)(&nm, &one, Y, &incOne, u, &incOne);

	F77_NAME(dgemv)(ytran, &p, &qm, &one, Z, &p, &beta[s], &nSamples, &zero, tmp_qm, &incOne);
      }

      F77_NAME(dtrsv)(lower, ntran, nUnit, &nm, C, &nm, u, &incOne);//L_1u = (y-X\mu_beta) or (y-X\beta)
     
      F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &nm, &qm, &one, C, &nm, c, &nm);//L_1v = c_0

      //for each location
      for(j = 0; j < q; j++){

	//mu
	F77_NAME(dgemv)(ytran, &nm, &m, &one, &c[(j*m)*nm], &nm, u, &incOne, &zero, tmp_m, &incOne);
	F77_NAME(daxpy)(&m, &one, &tmp_qm[j*m], &incOne, tmp_m, &incOne);

	//var
	F77_NAME(dgemm)(ytran, ntran, &m, &m, &nm, &negOne, &c[(j*m)*nm], &nm, &c[(j*m)*nm], &nm, &zero, tmp_mm, &m);
	F77_NAME(daxpy)(&mm, &one, K, &incOne, tmp_mm, &incOne);
	
	if(nugget){
	  for(k = 0; k < m; k++){
	    if(PsiDiag){
	      tmp_mm[k*m+k] += Psi[k];
	    }else{
	      F77_NAME(daxpy)(&mm, &one, Psi, &incOne, tmp_mm, &incOne);
	    }
	  }
	}

      	if(betaPrior == "normal"){
	  F77_NAME(daxpy)(&mm, &one, &tmp_qmm[j*mm], &incOne, tmp_mm, &incOne);
      	}
	
	if(fitted[j] == 0){
	  F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");}   
	  mvrnorm(&REAL(predSamples_r)[s*qm+j*m], tmp_m, tmp_mm, m, false);
	}else{
	  F77_NAME(dcopy)(&m, tmp_m, &incOne, &REAL(predSamples_r)[s*qm+j*m], &incOne);
	}
	
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
