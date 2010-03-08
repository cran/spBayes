#include <iostream>
#include <string>
using namespace std;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"
#include "covmodel.h"
#include "covInvDet.h"


extern "C" {

  SEXP mvLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP PsiIW_r, SEXP LStarting_r, SEXP betaStarting_r,
	      SEXP LTuning_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r){


    /*****************************************
                Common variables
    *****************************************/
    int i,j,k,l,info,nProtect= 0;
    char const *lower = "L";
    char const *upper = "U";
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
    
    double *betaStarting = REAL(betaStarting_r);

    double *LStarting = REAL(LStarting_r);
    double PsiIW_df = REAL(VECTOR_ELT(PsiIW_r, 0))[0]; 
    double *PsiIW_S = REAL(VECTOR_ELT(PsiIW_r, 1));
    
    //tuning
    double *LTuning = REAL(LTuning_r);

    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    double *L = (double *) R_alloc(mm, sizeof(double));

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
       
      
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);

      Rprintf("\tbeta flat.\n\n");   
      Rprintf("\n"); 
      
      Rprintf("\tPsi IW hyperpriors df=%.5f, S=\n", PsiIW_df);
      printMtrx(PsiIW_S, m, m);
      Rprintf("\n"); 
      
      Rprintf("Metropolis tuning values:\n");
      Rprintf("\tL tuning:\n");
      Rprintf("\t"); printVec(LTuning, nLTr);
      Rprintf("\n"); 
   
      Rprintf("Metropolis starting values:\n");
      covExpand(LStarting, L, m);
      Rprintf("\tL starting\n");
      printMtrx(L, m, m);
      Rprintf("\n");

    }
    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/
    //parameters
    int nMParams = nLTr; 
    int LIndx = 0; 

    double *params = (double *) R_alloc(nMParams, sizeof(double));
    
    //set starting
    covTrans(LStarting, &params[LIndx], m);
  
    //Beta parameter and set starting
    double *beta = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, betaStarting, &incOne, beta, &incOne);

    //samples and random effects
    int nParams = p+nMParams;

    SEXP samples_r, accept_r;

    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    double *samples = REAL(samples_r);

    PROTECT(accept_r = allocMatrix(REALSXP, 1, 1)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, rtnStatus=0, accept=0, batchAccept = 0;
    double logPostCurrent = 0, logPostCand = 0, logDetK, SKtrace, logDetPsi;

    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *tmp_nm1 = (double *) R_alloc(nm, sizeof(double));

    double *candParams = (double *) R_alloc(nMParams, sizeof(double));
    double *Psi = (double *) R_alloc(mm, sizeof(double));
    double *C = (double *) R_alloc(nm*nm, sizeof(double)); zeros(C, nm*nm);
    double *I_n = (double *) R_alloc(nn, sizeof(double)); identity(I_n, n);
    double logMHRatio;
    
    double *S_beta = (double *) R_alloc(p*p, sizeof(double));
    double *Mu_beta = (double *) R_alloc(p, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double)); 
    double *tmp_nmp = (double *) R_alloc(nm*p, sizeof(double)); 

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
      
      //
      //Current
      //
      covTransInvExpand(&params[LIndx], L, m);
 
      //Psi = L'L
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, L, &m, L, &m, &zero, Psi, &m);
      
      logDetPsi = 0.0;
      F77_NAME(dpotrf)(lower, &m, Psi, &m, &info); if(info != 0){error("Cholesky failed\n");}
      for(i = 0; i < m; i++) logDetPsi += 2.0*log(Psi[i*m+i]);
      F77_NAME(dpotri)(lower, &m, Psi, &m, &info); if(info != 0){error("Cholesky failed\n");}

      kron(I_n, n, n, Psi, m, m, C, nm, nm);

      //
      //Update Beta
      //
      //finish the Gibbs
      F77_NAME(dsymm)(lside, lower, &nm, &p, &one, C, &nm, X, &nm, &zero, tmp_nmp, &nm);
      F77_NAME(dgemm)(ytran, ntran, &p, &p, &nm, &one, X, &nm, tmp_nmp, &nm, &zero, S_beta, &p);
      
      F77_NAME(dpotrf)(lower, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
      F77_NAME(dpotri)(lower, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky inverse failed\n" << endl;}
      
      F77_NAME(dsymv)(lower, &nm, &one, C, &nm, Y, &incOne, &zero, tmp_nm, &incOne);
      F77_NAME(dgemv)(ytran, &nm, &p, &one, X, &nm, tmp_nm, &incOne, &zero, tmp_p, &incOne);
      F77_NAME(dsymv)(lower, &p, &one, S_beta, &p, tmp_p, &incOne, &zero, Mu_beta, &incOne);
      
      //Gibbs draw
      //take lower for the chol for the mv draw
      F77_NAME(dpotrf)(lower, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
      mvrnorm(beta, Mu_beta, S_beta, p, false);

      //
      //Likelihood with Jacobian   
      //
      logPostCurrent = 0.0;
      
      //
      //Jacobian and IW priors for K = A'A and Psi = L'L
      //
      
      //L'L prior with jacob.
      logDetK = 0.0;
      SKtrace = 0.0; 
      
      for(i = 0; i < m; i++) logDetK += 2*log(L[i*m+i]);
      
      //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
      for(i = 0; i < m; i++) logPostCurrent += (m-i)*log(L[i*m+i])+log(L[i*m+i]);
      
      //get S*K^-1, already have the chol of Psi (i.e., L)
      F77_NAME(dpotri)(lower, &m, L, &m, &info); if(info != 0){cout << "c++ error: L Cholesky inverse failed\n" << endl;}
      F77_NAME(dsymm)(rside, lower, &m, &m, &one, L, &m, PsiIW_S, &m, &zero, tmp_mm, &m);
      for(i = 0; i < m; i++){SKtrace += tmp_mm[i*m+i];}
      logPostCurrent += -0.5*(PsiIW_df+m+1)*logDetK - 0.5*SKtrace;
      
      //Y-XB
      F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
      F77_NAME(daxpy)(&nm, &one, Y, &incOne, tmp_nm, &incOne);
      
      for(i = 0; i < n; i++){
	F77_NAME(dsymv)(lower, &m, &one, Psi, &m, &tmp_nm[i*m], &incOne, &zero, &tmp_nm1[i*m], &incOne);
      }
      
      logPostCurrent += -(n/2.0)*logDetPsi-0.5*F77_NAME(ddot)(&nm, tmp_nm, &incOne, tmp_nm1, &incOne);

      //
      //Candidate
      //

      //propose   
      for(i = 0; i < nLTr; i++){
	candParams[LIndx+i] = rnorm(params[LIndx+i], LTuning[i]);
      }

      covTransInvExpand(&candParams[LIndx], L, m);
      
      //Psi = L'L
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, L, &m, L, &m, &zero, Psi, &m);
      
      logDetPsi = 0.0;
      F77_NAME(dpotrf)(lower, &m, Psi, &m, &info); if(info != 0){error("Cholesky failed\n");}
      for(i = 0; i < m; i++) logDetPsi += 2.0*log(Psi[i*m+i]);
      F77_NAME(dpotri)(lower, &m, Psi, &m, &info); if(info != 0){error("Cholesky failed\n");}

      kron(I_n, n, n, Psi, m, m, C, nm, nm);

      //
      //Likelihood with Jacobian   
      //
      logPostCand = 0.0;
      
      //
      //Jacobian and IW priors for K = A'A and Psi = L'L
      //
      
      //L'L prior with jacob.
      logDetK = 0.0;
      SKtrace = 0.0; 
      
      for(i = 0; i < m; i++) logDetK += 2*log(L[i*m+i]);
      
      //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
      for(i = 0; i < m; i++) logPostCand += (m-i)*log(L[i*m+i])+log(L[i*m+i]);
      
      //get S*K^-1, already have the chol of Psi (i.e., L)
      F77_NAME(dpotri)(lower, &m, L, &m, &info); if(info != 0){cout << "c++ error: L Cholesky inverse failed\n" << endl;}
      F77_NAME(dsymm)(rside, lower, &m, &m, &one, L, &m, PsiIW_S, &m, &zero, tmp_mm, &m);
      for(i = 0; i < m; i++){SKtrace += tmp_mm[i*m+i];}
      logPostCand += -0.5*(PsiIW_df+m+1)*logDetK - 0.5*SKtrace;
      
      //Y-XB
      F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
      F77_NAME(daxpy)(&nm, &one, Y, &incOne, tmp_nm, &incOne);
      
      for(i = 0; i < n; i++){
	F77_NAME(dsymv)(lower, &m, &one, Psi, &m, &tmp_nm[i*m], &incOne, &zero, &tmp_nm1[i*m], &incOne);
      }
      
      logPostCand += -(n/2.0)*logDetPsi-0.5*F77_NAME(ddot)(&nm, tmp_nm, &incOne, tmp_nm1, &incOne);

      //
      //MH accept/reject	
      //      
  
      //MH ratio with adjustment
      logMHRatio = logPostCand - logPostCurrent;
      
      if(runif(0.0,1.0) <= exp(logMHRatio)){
	F77_NAME(dcopy)(&nMParams, candParams, &incOne, params, &incOne);
	accept++;
	batchAccept++;
      }
    
      /******************************
          Save samples and report
      *******************************/
      F77_NAME(dcopy)(&p, beta, &incOne, &samples[s*nParams], &incOne);
      F77_NAME(dcopy)(&nMParams, params, &incOne, &samples[s*nParams+p], &incOne);
      
      //report
      if(verbose){
	if(status == nReport){
	  Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
	  Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
	  Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept/s);
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	  batchAccept = 0;
	}
      }
      status++;
      
      
      R_CheckUserInterrupt();
    }//end sample loop
    PutRNGstate();
    
    //final status report
    if(verbose){
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
    }
    Rprintf("-------------------------------------------------\n");
    #ifdef Win32
    R_FlushConsole();
    #endif

    //untransform variance variables
    for(s = 0; s < nSamples; s++){
      covTransInv(&samples[s*nParams+p+LIndx], &samples[s*nParams+p+LIndx], m);
    }   


    //calculate acceptance rate
    REAL(accept_r)[0] = 100.0*accept/s;

    //make return object
    SEXP result, resultNames;
    
    int nResultListObjs = 2;

    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;

   //samples
    SET_VECTOR_ELT(result, 0, samples_r);
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));
  
    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    
  }
}
