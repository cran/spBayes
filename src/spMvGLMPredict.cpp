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


extern"C" {
  
  void dbdimm_(int *transa, int *mb, int *n, int *kb, const double *alpha, int *descra,
	       double *val, int *blda, int *ibdiag, int *nbdiag, int *lb,
	       double *b, int *ldb, const double *beta, double *c, int *ldc, double *work, int *lwork);
}


extern "C" {

  SEXP spMvGLMPredict(SEXP family_r, SEXP X_r, SEXP Y_r, SEXP isPp_r, SEXP isModPp_r, SEXP n_r, SEXP m_r, SEXP p_r, SEXP q_r,  
		      SEXP beta_r, SEXP A_r, SEXP phi_r, SEXP nu_r,
		      SEXP nPred_r, SEXP predX_r, SEXP obsD_r, SEXP predD_r, SEXP predObsD_r, SEXP obsKnotsD_r, SEXP knotsD_r, SEXP predKnotsD_r,
		      SEXP covModel_r, SEXP nSamples_r, SEXP w_r, SEXP w_str_r, SEXP spEffects_r, SEXP verbose_r){
    

    /*****************************************
                Common variables
    *****************************************/
    int i,j,k,l,info,nProtect= 0;
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

    string family = CHAR(STRING_ELT(family_r,0));
    double *Y = REAL(Y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int q = INTEGER(q_r)[0];
    bool verbose = static_cast<bool>(INTEGER(verbose_r)[0]);

    int g = INTEGER(nPred_r)[0];
    double *predX = REAL(predX_r);

    int nSamples = INTEGER(nSamples_r)[0];

    //covariance model
    string covModel = CHAR(STRING_ELT(covModel_r,0));

    //pre-computed effects
    bool spEffects = static_cast<bool>(INTEGER(spEffects_r)[0]);

    //if predictive process
    bool isPp = static_cast<bool>(INTEGER(isPp_r)[0]);
    bool isModPp = static_cast<bool>(INTEGER(isModPp_r)[0]);

    double *knotsD = NULL;
    double *obsKnotsD = NULL;
    double *predKnotsD = NULL;
    double *obsD = NULL;
    double *predObsD = NULL;
    double *predD = REAL(predD_r);

    if(isPp){
      knotsD = REAL(knotsD_r);
      obsKnotsD = REAL(obsKnotsD_r);
      predKnotsD = REAL(predKnotsD_r);
    }else{
      obsD = REAL(obsD_r);
      predObsD = REAL(predObsD_r);
    }

    double *beta = REAL(beta_r);
    double *A = REAL(A_r);
    double *phi = REAL(phi_r);
   
    double *nu = NULL;
    if(covModel == "matern"){
      nu = REAL(nu_r);
    }

    double *w = REAL(w_r);
    double *w_str = NULL;

    if(isPp){
      w_str = REAL(w_str_r); 
    }

    /*****************************************
        Set-up cov. model function pointer
    *****************************************/
    int nPramPtr = 1;
    
    void (covmodel::*cov1ParamPtr)(double, double &, double &) = NULL; 
    void (covmodel::*cov2ParamPtr)(double, double, double &, double&) = NULL;
    
    if(covModel == "exponential"){
      cov1ParamPtr = &covmodel::exponential;
    }else if(covModel == "spherical"){
      cov1ParamPtr = &covmodel::spherical;
    }else if(covModel == "gaussian"){
      cov1ParamPtr = &covmodel::gaussian;
    }else if(covModel == "matern"){
      cov2ParamPtr = &covmodel::matern;
      nPramPtr = 2;
    }else{
      error("c++ error: cov.model is not correctly specified");
    }
   
    //my covmodel object for calling cov function
    covmodel *covModelObj = new covmodel;

    /*****************************************
         Set-up sample matrices etc.
    *****************************************/
    int nn = n*n;
    int mm = m*m;
    int nm = n*m;
    int nq = n*q;
    int qq = q*q;
    int qm = q*m;
    int qmqm = qm*qm;
    int gm = g*m;
    int gqm = g*qm;
    int nLTr = m*(m-1)/2+m;

    SEXP wPred_r, yPred_r;

    PROTECT(wPred_r = allocMatrix(REALSXP, gm, nSamples)); nProtect++; 
    double *wPred = REAL(wPred_r);

    PROTECT(yPred_r = allocMatrix(REALSXP, gm, nSamples)); nProtect++; 
    double *yPred = REAL(yPred_r);
  
    double *S_obs = NULL;
    double *S_predObs = NULL;
    double *S_pred = NULL;
    double *S_predKnots = NULL;
    double *S_knots = NULL;   
    double *tmp_gmnm = NULL;
    double *tmp_gmqm = NULL;
    double *E = NULL;
    double *K = NULL;
 
    if(!isPp){
      S_obs = (double *) R_alloc(nm*nm, sizeof(double));
      S_predObs = (double *) R_alloc(gm*nm, sizeof(double));
      S_pred = (double *) R_alloc(gm*gm, sizeof(double));
      tmp_gmnm = (double *) R_alloc(gm*nm, sizeof(double));
    }else{
      S_predKnots = (double *) R_alloc(gm*qm, sizeof(double));
      S_knots = (double *) R_alloc(qm*qm, sizeof(double));
      K = (double *) R_alloc(mm, sizeof(double));
      tmp_gmqm = (double *) R_alloc(gm*qm, sizeof(double));

      if(isModPp){
	E = (double *) R_alloc(m*m*g, sizeof(double));
      }
    }
    
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm1 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm2 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_gm = (double *) R_alloc(gm, sizeof(double));
    double *tmp_gmgm = (double *) R_alloc(gm*gm, sizeof(double));
    double *AA = (double *) R_alloc(mm, sizeof(double));
  
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, nReport = 100;

     if(verbose){
       Rprintf("-------------------------------------------------\n");
       Rprintf("\t\tStarting prediction\n");
       Rprintf("-------------------------------------------------\n");
       #ifdef Win32
       R_FlushConsole();
       #endif
     }

     GetRNGstate();

     //
     //Non predictive process
     //
     if(!isPp){

       for(s = 0; s < nSamples; s++){
	 
	 zeros(AA, m*m);
	 for(i = 0, k = 0; i < m; i++){
	   for(j = i; j < m; j++, k++){
	     AA[i*m+j] = A[s*nLTr+k];
	   }
	 }

	 //make S_obs
	 for(i = 0; i < n; i++){
	   for(j = 0; j < n; j++){
	     
	     zeros(tmp_mm, mm);
	     
	     for(k = 0; k < m; k++){
	       if(nPramPtr == 1)
		 (covModelObj->*cov1ParamPtr)(phi[s*m+k], tmp_mm[k*m+k], obsD[j*n+i]);
	       else //i.e., 2 parameter matern
		 (covModelObj->*cov2ParamPtr)(phi[s*m+k], nu[s*m+k], tmp_mm[k*m+k], obsD[j*n+i]);
	     }
	     
	     F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, AA, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	     F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, AA, &m, &zero, tmp_mm, &m);
	     
	     for(k = 0; k < m; k++){
	       for(l = 0; l < m; l++){
		 S_obs[((j*m+l)*nm)+(i*m+k)] = tmp_mm[l*m+k];
		 tmp_mm[l*m+k] = 0.0; //zero out
	       }
	     }
	   }
	 }
	 
	 //make S_PredObs
	 for(i = 0; i < g; i++){
	   for(j = 0; j < n; j++){
	     
	     zeros(tmp_mm, mm);
	     
	     for(k = 0; k < m; k++){
	       if(nPramPtr == 1)
		 (covModelObj->*cov1ParamPtr)(phi[s*m+k], tmp_mm[k*m+k], predObsD[j*g+i]);
	       else //i.e., 2 parameter matern
		 (covModelObj->*cov2ParamPtr)(phi[s*m+k], nu[s*m+k], tmp_mm[k*m+k], predObsD[j*g+i]);
	     }
	     
	     F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, AA, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	     F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, AA, &m, &zero, tmp_mm, &m);
	     
	     for(k = 0; k < m; k++){
	       for(l = 0; l < m; l++){
		 S_predObs[((j*m+l)*gm)+(i*m+k)] = tmp_mm[l*m+k];
		 tmp_mm[l*m+k] = 0.0; //zero out
	       }
	     }
	   }
	 }

	 //
	 //make S_pred
	 //
	 for(i = 0; i < g; i++){
	   for(j = 0; j < g; j++){
	     
	     zeros(tmp_mm, mm);
	     
	     for(k = 0; k < m; k++){
	       if(nPramPtr == 1)
		 (covModelObj->*cov1ParamPtr)(phi[s*m+k], tmp_mm[k*m+k], predD[j*g+i]);
	       else //i.e., 2 parameter matern
		 (covModelObj->*cov2ParamPtr)(phi[s*m+k], nu[s*m+k], tmp_mm[k*m+k], predD[j*g+i]);
	     }
	     
	     F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, AA, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	     F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, AA, &m, &zero, tmp_mm, &m);
	     
	     for(k = 0; k < m; k++){
	       for(l = 0; l < m; l++){
		 S_pred[((j*m+l)*gm)+(i*m+k)] = tmp_mm[l*m+k];
		 tmp_mm[l*m+k] = 0.0; //zero out
	       }
	     }
	   }
	 }

	 F77_NAME(dpotrf)(lower, &nm, S_obs, &nm, &info); if(info != 0){error("c++ error: Cholesky failed in spMvLMPredict\n");}
	 F77_NAME(dpotri)(lower, &nm, S_obs, &nm, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spMvLMPredict\n");}	 

	 //get Mu
	 F77_NAME(dsymm)(rside, lower, &gm, &nm, &one, S_obs, &nm, S_predObs, &gm, &zero, tmp_gmnm, &gm);
	 F77_NAME(dgemv)(ntran, &gm, &nm, &one, tmp_gmnm, &gm, &w[s*nm], &incOne, &zero, tmp_gm, &incOne);

	 //get Sigma
	 F77_NAME(dgemm)(ntran, ytran, &gm, &gm, &nm, &one, tmp_gmnm, &gm, S_predObs, &gm, &zero, tmp_gmgm, &gm);
	 for(i = 0; i < gm*gm; i++) S_pred[i] -= tmp_gmgm[i];
	 
	 F77_NAME(dpotrf)(lower, &gm, S_pred, &gm, &info); if(info != 0){error("c++ error: Cholesky failed in spMvLMPredict\n");}
	 mvrnorm(&wPred[s*gm], tmp_gm, S_pred, gm, false);
	 
	 F77_NAME(dgemv)(ntran, &gm, &p, &one, predX, &gm, &beta[s*p], &incOne, &zero, tmp_gm, &incOne);
	
	 if(family == "binomial"){
	   for(i = 0; i < gm; i++){
	     yPred[s*gm+i] =  rbinom(1, 1.0/(1.0+exp(-1.0*(tmp_gm[i]+wPred[s*gm+i]))));
	   }
	 }else if(family == "poisson"){
	   for(i = 0; i < gm; i++){
	     yPred[s*gm+i] =  rpois(exp(tmp_gm[i]+wPred[s*gm+i]));	   
	   }
	 }else{
	   error("c++ error: family misspecification in spMvGLMPredict\n");
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
       
     }else{//predictive process prediction
       
       for(s = 0; s < nSamples; s++){	 
	 
	 zeros(AA, m*m);
	 for(i = 0, k = 0; i < m; i++){
	   for(j = i; j < m; j++, k++){
	     AA[i*m+j] = A[s*nLTr+k];
	   }
	 }

	 //make S_predKnots
	 for(i = 0; i < g; i++){
	   for(j = 0; j < q; j++){
	     
	     zeros(tmp_mm, mm);
	     
	     for(k = 0; k < m; k++){
	       if(nPramPtr == 1)
		 (covModelObj->*cov1ParamPtr)(phi[s*m+k], tmp_mm[k*m+k], predKnotsD[j*g+i]);
	       else //i.e., 2 parameter matern
		 (covModelObj->*cov2ParamPtr)(phi[s*m+k], nu[s*m+k], tmp_mm[k*m+k], predKnotsD[j*g+i]);
	     }
	     
	     F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, AA, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	     F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, AA, &m, &zero, tmp_mm, &m);
	     
	     for(k = 0; k < m; k++){
	       for(l = 0; l < m; l++){
		 S_predKnots[((j*m+l)*gm)+(i*m+k)] = tmp_mm[l*m+k];
		 tmp_mm[l*m+k] = 0.0; //zero out
	       }
	     }
	   }
	 }

	 //
	 //make S_knots
	 //
	 for(i = 0; i < q; i++){
	   for(j = 0; j < q; j++){
	     
	     zeros(tmp_mm, mm);
	     
	     for(k = 0; k < m; k++){
	       if(nPramPtr == 1)
		 (covModelObj->*cov1ParamPtr)(phi[s*m+k], tmp_mm[k*m+k], knotsD[j*q+i]);
	       else //i.e., 2 parameter matern
		 (covModelObj->*cov2ParamPtr)(phi[s*m+k], nu[s*m+k], tmp_mm[k*m+k], knotsD[j*q+i]);
	     }
	     
	     F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, AA, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	     F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, AA, &m, &zero, tmp_mm, &m);
	     
	     for(k = 0; k < m; k++){
	       for(l = 0; l < m; l++){
		 S_knots[((j*m+l)*qm)+(i*m+k)] = tmp_mm[l*m+k];
		 tmp_mm[l*m+k] = 0.0; //zero out
	       }
	     }
	   }
	 }

	 F77_NAME(dpotrf)(lower, &qm, S_knots, &qm, &info); if(info != 0){error("c++ error: Cholesky failed in spMvLMPredict\n");}
	 F77_NAME(dpotri)(lower, &qm, S_knots, &qm, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spMvLMPredict\n");}	 
	 
	 //S_predKnots S_knots^{-1}
	 F77_NAME(dsymm)(rside, lower, &gm, &qm, &one, S_knots, &qm, S_predKnots, &gm, &zero, tmp_gmqm, &gm);
	 
// 	 if(isModPp){
	   
// 	   //K = A'A and Psi = L'L
// 	   F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, AA, &m, AA, &m, &zero, K, &m);
	  
// 	   //C = ct C_str^{-1} t(ct)
// 	   F77_NAME(dgemm)(ntran, ytran, &gm, &gm, &qm, &one, tmp_gmqm, &gm, S_predKnots, &gm, &zero, tmp_gmgm, &gm);
	   
// 	   //make E = (Psi + V - blk[C])^{-1}
// 	   for(i = 0, j = 0; i < g; i++){
// 	     for(k = 0; k < m; k++){
// 	       for(l = 0; l < m; l++){
// 		 E[j] = K[l*m+k]-tmp_gmgm[(i*m+l)*gm+(i*m+k)];
// 		 j++;
// 	       }
// 	     } 
// 	   }
	   
// 	 }

	 F77_NAME(dgemv)(ntran, &gm, &p, &one, predX, &gm, &beta[s*p], &incOne, &zero, tmp_gm, &incOne);
	 F77_NAME(dgemv)(ntran, &gm, &qm, &one, tmp_gmqm, &gm, &w_str[s*qm], &incOne, &zero, &wPred[s*gm], &incOne);

	 if(family == "binomial"){
	   for(i = 0; i < gm; i++){
	     yPred[s*gm+i] =  rbinom(1, 1.0/(1.0+exp(-1.0*(tmp_gm[i]+wPred[s*gm+i]))));
	   }
	 }else if(family == "poisson"){
	   for(i = 0; i < gm; i++){
	     yPred[s*gm+i] =  rpois(exp(tmp_gm[i]+wPred[s*gm+i]));	   
	   }
	 }else{
	   error("c++ error: family misspecification in spMvGLMPredict\n");
	 }

	 //normal
// 	 for(i = 0; i < g; i++){  
   
// 	   if(isModPp){
// 	     F77_NAME(dpotrf)(lower, &m, &E[i*mm], &m, &info); if(info != 0){error("c++ error: Cholesky failed in spMvLMPredict\n");}
// 	     mvrnorm(&yPred[s*gm+i*m], &tmp_gm[i*m], &E[i*mm], m, false);
// 	   }else{
// 	     mvrnorm(&yPred[s*gm+i*m], &tmp_gm[i*m], LL, m, false);
// 	   }

// 	 }

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
     }

     PutRNGstate();

     //make return object
     SEXP result, resultNames;
     
     int nResultListObjs = 0;
     
     nResultListObjs = 2;
     
     PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
     PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
     
     
     SET_VECTOR_ELT(result, 0, wPred_r);
     SET_VECTOR_ELT(resultNames, 0, mkChar("w.pred"));

     SET_VECTOR_ELT(result, 1, yPred_r);
     SET_VECTOR_ELT(resultNames, 1, mkChar("y.pred"));
     
     namesgets(result, resultNames);
     
     //unprotect
     UNPROTECT(nProtect);
     
     return(result);

  }
}
