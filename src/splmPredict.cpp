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



extern "C" {

  SEXP spLMPredict(SEXP X_r, SEXP Y_r, SEXP isPp_r, SEXP isModPp_r, SEXP n_r, SEXP m_r, SEXP p_r, SEXP nugget_r, SEXP beta_r, SEXP sigmaSq_r, 
		   SEXP tauSq_r, SEXP phi_r, SEXP nu_r,
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

    double *Y = REAL(Y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    bool verbose = static_cast<bool>(INTEGER(verbose_r)[0]);

    int q = INTEGER(nPred_r)[0];
    double *predX = REAL(predX_r);

    int nSamples = INTEGER(nSamples_r)[0];

    //covariance model
    string covModel = CHAR(STRING_ELT(covModel_r,0));

    //pre-computed effects
    bool spEffects = static_cast<bool>(INTEGER(spEffects_r)[0]);

    //if predictive process
    bool isPp = static_cast<bool>(INTEGER(isPp_r)[0]);
    bool isModPp = static_cast<bool>(INTEGER(isModPp_r)[0]);

    int m = 0;
    double *knotsD = NULL;
    double *obsKnotsD = NULL;
    double *predKnotsD = NULL;
    double *obsD = NULL;
    double *predObsD = NULL;
    double *predD = REAL(predD_r);

    if(isPp){
      m = INTEGER(m_r)[0];
      knotsD = REAL(knotsD_r);
      obsKnotsD = REAL(obsKnotsD_r);
      predKnotsD = REAL(predKnotsD_r);
    }else{
      obsD = REAL(obsD_r);
      predObsD = REAL(predObsD_r);
    }

    int nn = n*n;
    int nm = n*m;
    int mm = m*m;
    int qn = q*n;
    int qq = q*q;
    int qm = q*m;
   
    //if nugget
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    double *tauSq = NULL;

    if(nugget){
      tauSq = REAL(tauSq_r);
    }

    //if matern
    double *nu = NULL;
    if(covModel == "matern"){
      nu = REAL(nu_r);
    }

    double *phi = NULL;
    double *sigmaSq = NULL;
    double *beta = NULL;

    phi = REAL(phi_r);
    sigmaSq = REAL(sigmaSq_r);
    beta = REAL(beta_r);

    /*****************************************
        Set-up cov. model function pointer
    *****************************************/
    bool onePramPtr = true;
    
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
      onePramPtr = false;
    }else{
      error("c++ error: cov.model is not correctly specified");
    }
   
    //my covmodel object for calling cov function
    covmodel *covModelObj = new covmodel;

    /*****************************************
         Set-up sample matrices etc.
    *****************************************/
    SEXP w_pred_r, y_pred_r;

    PROTECT(w_pred_r = allocMatrix(REALSXP, q, nSamples)); nProtect++; 
    double *w_pred = REAL(w_pred_r);

    PROTECT(y_pred_r = allocMatrix(REALSXP, q, nSamples)); nProtect++; 
    double *y_pred = REAL(y_pred_r);
  
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, nReport = 100;

    double *C = (double *) R_alloc(nn, sizeof(double)); 
    double *ct = NULL;
    double *C_str = NULL;
    double *wMu = NULL;
    double *w_strMu = NULL;
    double *E = NULL;
    double *Einv = NULL;
    double *tmp_m = NULL;
    double *tmp_nm = NULL;
    double *tmp_nm1 = NULL;
    double *tmp_mm = NULL;
    double *tmp_qm = NULL;
    double *tmp_qm1 = NULL;

    double *w = REAL(w_r);
    double *w_str = NULL;
    
    if(isPp){
      w_str = REAL(w_str_r);
    }  

    if(isPp){
      ct = (double *) R_alloc(nm, sizeof(double));
      C_str = (double *) R_alloc(mm, sizeof(double));
      w_strMu = (double *) R_alloc(m, sizeof(double));
      E = (double *) R_alloc(n, sizeof(double));
      Einv = (double *) R_alloc(n, sizeof(double)); 
      tmp_m = (double *) R_alloc(m, sizeof(double));
      tmp_nm = (double *) R_alloc(nm, sizeof(double));
      tmp_nm1 = (double *) R_alloc(nm, sizeof(double));
      tmp_mm = (double *) R_alloc(mm, sizeof(double));
      tmp_qm = (double *) R_alloc(qm, sizeof(double));
      tmp_qm1 = (double *) R_alloc(qm, sizeof(double));
    }else{
      wMu = (double *) R_alloc(n, sizeof(double)); zeros(wMu, n);
    }
     
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_nn = (double *) R_alloc(nn, sizeof(double));     
    double *tmp_qn = (double *) R_alloc(qn, sizeof(double));
    double *tmp_qn1 = (double *) R_alloc(qn, sizeof(double)); 
    double *tmp_qq = (double *) R_alloc(qq, sizeof(double));
    double *tmp_qq1 = (double *) R_alloc(qq, sizeof(double));
    double *tmp_q = (double *) R_alloc(q, sizeof(double));
    double *tmp_q1 = (double *) R_alloc(q, sizeof(double));
    double *tmp_q2 = (double *) R_alloc(q, sizeof(double));
    double sigmaSqTmp, tauSqTmp;

     GetRNGstate();
     
     if(verbose){
       Rprintf("-------------------------------------------------\n");
       Rprintf("\t\tStarting prediction\n");
       Rprintf("-------------------------------------------------\n");
       #ifdef Win32
       R_FlushConsole();
       #endif
     }

     status = 0; 

     if(!isPp){

       for(s = 0; s < nSamples; s++){
	 
	 //make the correlation matrix
	 for(i = 0; i < nn; i++){
	   if(onePramPtr)
	     (covModelObj->*cov1ParamPtr)(phi[s], tmp_nn[i], obsD[i]);
	   else //i.e., 2 parameter matern
	     (covModelObj->*cov2ParamPtr)(phi[s], nu[s], tmp_nn[i], obsD[i]);
	 }
	 
    	 for(i = 0; i < qn; i++){
	   if(onePramPtr)
	     (covModelObj->*cov1ParamPtr)(phi[s], tmp_qn[i], predObsD[i]);
	   else //i.e., 2 parameter matern
	     (covModelObj->*cov2ParamPtr)(phi[s], nu[s], tmp_qn[i], predObsD[i]);
	 }

    	 for(i = 0; i < qq; i++){
	   if(onePramPtr)
	     (covModelObj->*cov1ParamPtr)(phi[s], tmp_qq[i], predD[i]);
	   else //i.e., 2 parameter matern
	     (covModelObj->*cov2ParamPtr)(phi[s], nu[s], tmp_qq[i], predD[i]);
	 }

	 //scale by sigma^2
	 F77_NAME(dscal)(&nn, &sigmaSq[s], tmp_nn, &incOne);	
	 F77_NAME(dscal)(&qn, &sigmaSq[s], tmp_qn, &incOne);
	 F77_NAME(dscal)(&qq, &sigmaSq[s], tmp_qq, &incOne);
	 
	 //invert C_str
	 F77_NAME(dpotrf)(upper, &n, tmp_nn, &n, &info); if(info != 0){error("c++ error: Cholesky failed in spLMPredict\n");}
	 F77_NAME(dpotri)(upper, &n, tmp_nn, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spLMPredict\n");}
	 
	 //get Mu
	 F77_NAME(dsymm)(rside, upper, &q, &n, &one, tmp_nn, &n, tmp_qn, &q, &zero, tmp_qn1, &q);
	 F77_NAME(dgemv)(ntran, &q, &n, &one, tmp_qn1, &q, &w[s*n], &incOne, &zero, tmp_q, &incOne);

	 //get Sigma
	 F77_NAME(dgemm)(ntran, ytran, &q, &q, &n, &one, tmp_qn1, &q, tmp_qn, &q, &zero, tmp_qq1, &q);
	 for(i = 0; i < qq; i++) tmp_qq[i] -= tmp_qq1[i];

	 F77_NAME(dpotrf)(upper, &q, tmp_qq, &q, &info); if(info != 0){error("c++ error: Cholesky failed in spLMPredict\n");}
	 mvrnorm(&w_pred[s*q], tmp_q, tmp_qq, q, true);

	 F77_NAME(dgemv)(ntran, &q, &p, &one, predX, &q, &beta[s*p], &incOne, &zero, tmp_q, &incOne);

	 if(nugget){
	   for(i = 0; i < q; i++) y_pred[s*q+i] = rnorm(tmp_q[i]+w_pred[s*q+i], sqrt(tauSq[s]));
	 }else{
	   for(i = 0; i < q; i++) y_pred[s*q+i] = tmp_q[i]+w_pred[s*q+i];
	 }

	 report(s, nSamples, status, nReport, verbose);
       } //end sample loop

     }else{//predictive process prediction
       
       for(s = 0; s < nSamples; s++){	 
	 
	 //got w* above now get the mean components MVN(XB + ct C^{*-1} w* + \tild{\eps}, (sigma^2 - ct C^{*-1} c)^{-1})
	 
	 //make the correlation matrix
	 for(i = 0; i < mm; i++){
	   if(onePramPtr)
	     (covModelObj->*cov1ParamPtr)(phi[s], C_str[i], knotsD[i]);
	   else //i.e., 2 parameter matern
	     (covModelObj->*cov2ParamPtr)(phi[s], nu[s], C_str[i], knotsD[i]);
	 }
	 
	 for(i = 0; i < qm; i++){
	   if(onePramPtr)
	     (covModelObj->*cov1ParamPtr)(phi[s], tmp_qm[i], predKnotsD[i]);
	   else //i.e., 2 parameter matern
	     (covModelObj->*cov2ParamPtr)(phi[s], nu[s], tmp_qm[i], predKnotsD[i]);
	 }
	 
	 //scale by sigma^2
	 F77_NAME(dscal)(&mm, &sigmaSq[s], C_str, &incOne);	
	 F77_NAME(dscal)(&qm, &sigmaSq[s], tmp_qm, &incOne);
	 
	 //invert C_str
	 F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in spLMPredict\n");}
	 F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spLMPredict\n");}
	 
	 //ct C^{*-1}, where ct is now qxm 
	 F77_NAME(dsymm)(rside, upper, &q, &m, &one, C_str, &m, tmp_qm, &q, &zero, tmp_qm1, &q);
	 
	 //ct C^{*-1} c
	 F77_NAME(dgemm)(ntran, ytran, &q, &q, &m, &one, tmp_qm1, &q, tmp_qm, &q, &zero, tmp_qq, &q);

	 //ct C^{*-1} w*, tmp_q will be the mean
	 F77_NAME(dgemv)(ntran, &q, &m, &one, tmp_qm1, &q, &w_str[s*m], &incOne, &zero, &w_pred[s*q], &incOne);

	 //\tild{\eps}
	 if(isModPp){
	   for(i = 0; i < q; i++) tmp_q[i] = rnorm(0.0, sqrt(sigmaSq[s]-tmp_qq[i*q+i]));
	 }

	 //XB
	 F77_NAME(dgemv)(ntran, &q, &p, &one, predX, &q, &beta[s*p], &incOne, &zero, tmp_q2, &incOne);
	
	 if(isModPp){
	   for(i = 0; i < q; i++) y_pred[s*q+i] = rnorm(tmp_q2[i]+w_pred[s*q+i]+tmp_q[i], sqrt(tauSq[s]));
	 }else{
	   for(i = 0; i < q; i++) y_pred[s*q+i] = rnorm(tmp_q2[i]+w_pred[s*q+i], sqrt(tauSq[s]));
	 }

	 report(s, nSamples, status, nReport, verbose);
       } //end sample loop
     }

     PutRNGstate();

     //make return object
     SEXP result, resultNames;
     
     int nResultListObjs = 0;
     
     nResultListObjs = 2;
     
     PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
     PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
     
     
     SET_VECTOR_ELT(result, 0, w_pred_r);
     SET_VECTOR_ELT(resultNames, 0, mkChar("w.pred"));

     SET_VECTOR_ELT(result, 1, y_pred_r);
     SET_VECTOR_ELT(resultNames, 1, mkChar("y.pred"));
     
     namesgets(result, resultNames);
     
     //unprotect
     UNPROTECT(nProtect);
     
     return(result);

  }
}


//   //recover w or w_str if not pre-computed
//      if(!spEffects){
       
       
//        if(verbose){
// 	 Rprintf("-------------------------------------------------\n");
// 	 Rprintf("Recovering random spatial effects for prediction\n");
// 	 Rprintf("-------------------------------------------------\n");
//          #ifdef Win32
// 	 R_FlushConsole();
//          #endif
//        }

//        if(!isPp){
	 
// 	 for(s = 0; s < nSamples; s++){
	   
// 	   if(nugget){
// 	     //make the correlation matrix
// 	     for(i = 0; i < nn; i++){
// 	       if(onePramPtr)
// 		 (covModelObj->*cov1ParamPtr)(phi[s], C[i], obsD[i]);
// 	       else //i.e., 2 parameter matern
// 		 (covModelObj->*cov2ParamPtr)(phi[s], nu[s], C[i], obsD[i]);
// 	     }
	     
// 	     //invert C
// 	     F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in spLMPredict\n");}
// 	     F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spLMPredict\n");}
	     
// 	     //scale correlation matrix with 1/sigmasq and add 1/nugget to diag
// 	     sigmaSqTmp = 1.0/sigmaSq[s];
// 	     F77_NAME(dscal)(&nn, &sigmaSqTmp, C, &incOne);
	     
// 	     for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+1.0/tauSq[s];
	     
// 	     //invert C
// 	     F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in spLMPredict\n");}
// 	     F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spLMPredict\n");}
	     
// 	     //make w mu
// 	     F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s*p], &incOne, &zero, tmp_n, &incOne);
// 	     F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	     
// 	     tauSqTmp = 1.0/tauSq[s];
// 	     F77_NAME(dscal)(&n, &tauSqTmp, tmp_n, &incOne);
	     
// 	     F77_NAME(dsymv)(upper, &n, &one, C, &n, tmp_n, &incOne, &zero, wMu, &incOne);
	     
// 	     //chol for the mvnorm and draw
// 	     F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in spLMPredict\n");}
	     
// 	     mvrnorm(&w[s*n], wMu, C, n, true);
	     
// 	   }else{//no nugget so w is just resids
// 	     F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s*p], &incOne, &zero, &w[s*n], &incOne);
// 	     F77_NAME(daxpy)(&n, &one, Y, &incOne, &w[s*n], &incOne);
// 	   }

// 	   report(s, nSamples, status, nReport, verbose);
// 	 } //end sample loop
	 
//        }else{//is pp
	 
// 	 for(s = 0; s < nSamples; s++){	 
	   
// 	   //make the correlation matrix
// 	   for(i = 0; i < mm; i++){
// 	     if(onePramPtr)
// 	       (covModelObj->*cov1ParamPtr)(phi[s], C_str[i], knotsD[i]);
// 	     else //i.e., 2 parameter matern
// 	       (covModelObj->*cov2ParamPtr)(phi[s], nu[s], C_str[i], knotsD[i]);
// 	   }
	   
// 	   for(i = 0; i < nm; i++){
// 	     if(onePramPtr)
// 	       (covModelObj->*cov1ParamPtr)(phi[s], ct[i], obsKnotsD[i]);
// 	     else //i.e., 2 parameter matern
// 	       (covModelObj->*cov2ParamPtr)(phi[s], nu[s], ct[i], obsKnotsD[i]);
// 	   }
	   
// 	   //scale by sigma^2
// 	   F77_NAME(dscal)(&mm, &sigmaSq[s], C_str, &incOne);	
// 	   F77_NAME(dscal)(&nm, &sigmaSq[s], ct, &incOne);
	   
// 	   //invert C_str
// 	   F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in spLMPredict\n");}
// 	   F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spLMPredict\n");}
	   
// 	   //make w* Sigma
// 	   //ct C^{*-1}
// 	   F77_NAME(dsymm)(rside, upper, &n, &m, &one, C_str, &m, ct, &n, &zero, tmp_nm, &n);
	   	   
// 	   if(!isModPp){
// 	     for(i = 0; i < n; i++) Einv[i] = 1.0/(tauSq[s]);
// 	   }else{
// 	     //ct C^{*-1} c
// 	     F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, tmp_nn, &n);

// 	     for(i = 0; i < n; i++) Einv[i] = 1.0/(tauSq[s]+sigmaSq[s]-tmp_nn[i*n+i]);
// 	   }

// 	   diagmm(n, m, Einv, tmp_nm, tmp_nm1);
	   
// 	   //(C^{*-1} c) (1/E ct C^{*-1})
// 	   F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &one, tmp_nm, &n, tmp_nm1, &n, &zero, tmp_mm, &m);
	   
// 	   for(i = 0; i < mm; i++) C_str[i] = C_str[i] + tmp_mm[i];
	   
// 	   //invert C_str
// 	   F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
// 	   F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	   
// 	   //make w* mu
// 	   F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s*p], &incOne, &zero, tmp_n, &incOne);
// 	   F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	   
// 	   //(1/E ct C^{*-1})'(Y-XB)
// 	   F77_NAME(dgemv)(ytran, &n, &m, &one, tmp_nm1, &n, tmp_n, &incOne, &zero, tmp_m, &incOne);
// 	   F77_NAME(dsymv)(upper, &m, &one, C_str, &m, tmp_m, &incOne, &zero, w_strMu, &incOne);
	  
// 	   //chol for the mvnorm and draw
// 	   F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
// 	   mvrnorm(&w_str[s*m], w_strMu, C_str, m, true);
	  
// 	   //make \tild{w}
// 	   F77_NAME(dgemv)(ntran, &n, &m, &one, tmp_nm, &n, &w_str[s*m], &incOne, &zero, &w[s*n], &incOne);

// 	   report(s, nSamples, status, nReport, verbose);
// 	 } //end sample loop
//        } 
//      }
