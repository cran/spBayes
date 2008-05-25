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


  void diff2(int &length, double *a, double *b, double *c){
    for(int i = 0; i < length; i++) 
      c[i] = a[i]-b[i];
  }

  void diff3(int &length, double *a, double *b, double *c, double *d){
    for(int i = 0; i < length; i++) 
      d[i] = a[i]-b[i]-c[i];
  }

  
  SEXP splmDIC(SEXP X_r, SEXP Y_r, SEXP isPp_r, SEXP isModPp_r, SEXP n_r, SEXP m_r, SEXP p_r, SEXP nugget_r, SEXP beta_r, SEXP sigmaSq_r, 
	       SEXP tauSq_r, SEXP phi_r, SEXP nu_r,
	       SEXP obsD_r, SEXP obsKnotsD_r, SEXP knotsD_r,
	       SEXP covModel_r, SEXP nSamples_r, SEXP w_r, SEXP w_str_r, SEXP spEffects_r, SEXP DICMarg_r, SEXP DICUnmarg_r, SEXP verbose_r){
     
     
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
    bool verbose = static_cast<bool>(INTEGER(verbose_r)[0]);
    int nReport = 100;

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
    double *obsD = NULL;
    double *w = REAL(w_r);
    double *w_str = NULL;

    bool DICMarg = static_cast<bool>(INTEGER(DICMarg_r)[0]);
    bool DICUnmarg = static_cast<bool>(INTEGER(DICUnmarg_r)[0]);

    if(isPp){
      m = INTEGER(m_r)[0];
      knotsD = REAL(knotsD_r);
      obsKnotsD = REAL(obsKnotsD_r);
      w_str = REAL(w_str_r);
    }else{
      obsD = REAL(obsD_r);
    }

    int nn = n*n;
    int nm = n*m;
    int mm = m*m;
   
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
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, rtnStatus=0;
    double logPost = 0, detCand = 0, detCandC_str = 0, detCandE = 0;

    double *C = (double *) R_alloc(nn, sizeof(double)); 
    double *ct = NULL;
    double *C_str = NULL;
    double *wMu = NULL;
    double *w_strMu = NULL;
    double *E = NULL;
    double *Einv = NULL;
    double *tmp_nn = NULL;

    if(isPp){
      ct = (double *) R_alloc(nm, sizeof(double));
      C_str = (double *) R_alloc(mm, sizeof(double));
      w_strMu = (double *) R_alloc(m, sizeof(double));
      tmp_nn = (double *) R_alloc(nn, sizeof(double));
      E = (double *) R_alloc(n, sizeof(double));
      Einv = (double *) R_alloc(n, sizeof(double));
    }else{
      wMu = (double *) R_alloc(n, sizeof(double));
    }
 

    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_n1 = (double *) R_alloc(n, sizeof(double));
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *tmp_nm1 = (double *) R_alloc(nm, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));

    double tmp, sigmaSqTmp, tauSqTmp, tauSqInv, negTauSqInv;
   
    double *tmp_p = (double *) R_alloc(p, sizeof(double)); 
    double *tmp_np = (double *) R_alloc(n*p, sizeof(double)); 
 
    double *wMeans = (double *) R_alloc(n, sizeof(double)); zeros(wMeans, n);
    double *tildEpsMean = (double *) R_alloc(n, sizeof(double)); zeros(tildEpsMean, n);
    double *betaMean = (double *) R_alloc(p, sizeof(double)); zeros(betaMean, p);
    double nuMean, phiMean, tauSqMean, sigmaSqMean;

    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tCalculating DIC\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }
    
    GetRNGstate();

    
    if(!spEffects && DICUnmarg){
      if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\tRecovering spatial effects\n");
	Rprintf("-------------------------------------------------\n");
        #ifdef Win32
	R_FlushConsole();
        #endif
      } 
      
      
      for(s = 0; s < nSamples; s++){
	
	/******************************
         Recover w and w* if needed
	*******************************/
	if(!isPp){
	  if(nugget){
	    
	    //make the correlation matrix
	    for(i = 0; i < nn; i++){
	      if(onePramPtr)
		(covModelObj->*cov1ParamPtr)(phi[s], C[i], obsD[i]);
	      else //i.e., 2 parameter matern
		(covModelObj->*cov2ParamPtr)(phi[s], nu[s], C[i], obsD[i]);
	    }
	    
	    //invert C
	    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	    
	    //scale correlation matrix with 1/sigmasq and add 1/nugget to diag
	    sigmaSqTmp = 1.0/sigmaSq[s];
	    F77_NAME(dscal)(&nn, &sigmaSqTmp, C, &incOne);
	    
	    for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+1.0/tauSq[s];
	    
	    //invert C
	    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	    
	    //make w mu
	    F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s*p], &incOne, &zero, tmp_n, &incOne);
	    F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	    
	    tauSqTmp = 1.0/tauSq[s];
	    F77_NAME(dscal)(&n, &tauSqTmp, tmp_n, &incOne);
	    
	    F77_NAME(dsymv)(upper, &n, &one, C, &n, tmp_n, &incOne, &zero, wMu, &incOne);
	    
	    //chol for the mvnorm and draw
	    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    
	    mvrnorm(&w[s*n], wMu, C, n, true);
	    
	  }else{//no nugget so w is just resids
	    F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s*p], &incOne, &zero, &w[s*n], &incOne);
	    F77_NAME(daxpy)(&n, &one, Y, &incOne, &w[s*n], &incOne);
	  }
	}else{//using pp
	  
	  //w* ~ MVN(mu_w, Sigma_w)
	  //Sigma_w = [C^{*-1} + C^{*-1} C (1/E) C' C^{*-1}]^{-1}
	  //mu_w = Sigma_w [C^{*-1} C (1/E) (Y-XB)]
	  //where E = I \otimes (Psi + A'A) - Diag[C'(s_i) C^{*-1} C(s_i)]^n_{i=1}
	  //then w = C' C^{*-1} w*
	  
	  //make the correlation matrix
	  for(i = 0; i < mm; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi[s], C_str[i], knotsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi[s], nu[s], C_str[i], knotsD[i]);
	  }
	  
	  for(i = 0; i < nm; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi[s], ct[i], obsKnotsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi[s], nu[s], ct[i], obsKnotsD[i]);
	  }
	  
	  //scale by sigma^2
	  F77_NAME(dscal)(&mm, &sigmaSq[s], C_str, &incOne);	
	  F77_NAME(dscal)(&nm, &sigmaSq[s], ct, &incOne);
	  
	  //invert C_str
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  //make w* Sigma
	  //ct C^{*-1}
	  F77_NAME(dsymm)(rside, upper, &n, &m, &one, C_str, &m, ct, &n, &zero, tmp_nm, &n);
	  	  
	  if(!isModPp){
	    for(i = 0; i < n; i++) Einv[i] = 1.0/(tauSq[s]);
	  }else{
	    //ct C^{*-1} c
	    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, tmp_nn, &n);
	    
	    for(i = 0; i < n; i++) Einv[i] = 1.0/(tauSq[s]+sigmaSq[s]-tmp_nn[i*n+i]);
	   }

	  diagmm(n, m, Einv, tmp_nm, tmp_nm1);
	  
	  //(C^{*-1} c) (1/E ct C^{*-1})
	  F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &one, tmp_nm, &n, tmp_nm1, &n, &zero, tmp_mm, &m);
	  
	  for(i = 0; i < mm; i++) C_str[i] = C_str[i] + tmp_mm[i];
	  
	  //invert C_str
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  //make w* mu
	  F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s*p], &incOne, &zero, tmp_n, &incOne);
	  F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	  	  
	  //(1/E ct C^{*-1})'(Y-XB)
	  F77_NAME(dgemv)(ytran, &n, &m, &one, tmp_nm1, &n, tmp_n, &incOne, &zero, tmp_m, &incOne);
	  
	  F77_NAME(dsymv)(upper, &m, &one, C_str, &m, tmp_m, &incOne, &zero, w_strMu, &incOne);
	  
	  //chol for the mvnorm and draw
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  mvrnorm(&w_str[s*m], w_strMu, C_str, m, true);
	  
	  //make \tild{w}
	  F77_NAME(dgemv)(ntran, &n, &m, &one, tmp_nm, &n, &w_str[s*m], &incOne, &zero, &w[s*n], &incOne);
	}
	
	if(verbose){
	  if(status == 100){
	    Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
            #ifdef Win32
	    R_FlushConsole();
            #endif
	    status = 0;
	  }
	  status++;
	}
	R_CheckUserInterrupt();
	
      }//end sample loop
      
      if(verbose){
	Rprintf("Sampled: %i of %i, %3.2f%%\n", nSamples, nSamples, 100.0);
        #ifdef Win32
	R_FlushConsole();
        #endif
	status = 0;
	R_CheckUserInterrupt();
      }
          
    }//fi get w
    PutRNGstate();

 
    /******************************
              DIC setup
    *******************************/
    if(DICUnmarg){
      //get w row means      
      for(s = 0; s < nSamples; s++){
	for(j = 0; j < n; j++){
	  wMeans[j] += w[s*n+j];
	}
      }
      
      for(j = 0; j < n; j++)
	wMeans[j] = wMeans[j]/nSamples;
    }
    
    //get means for other parameters (if pp then we need to get mean \tild{\eps} below)
    nuMean = phiMean = tauSqMean = sigmaSqMean = 0;
    for(s = 0; s < nSamples; s++){

      for(i = 0; i < p; i++){
	betaMean[i] += beta[s*p+i];
      }

      sigmaSqMean += sigmaSq[s];

      if(nugget)
	tauSqMean += tauSq[s];

      phiMean += phi[s];

      if(covModel == "matern"){
	nuMean += nu[s];
      }

    }

    for(i = 0; i < p; i++){
      betaMean[i] = betaMean[i]/nSamples;
    }
 

    phiMean = phiMean/nSamples;
    tauSqMean = tauSqMean/nSamples;
    sigmaSqMean = sigmaSqMean/nSamples;
   
    if(covModel == "matern"){
      nuMean = nuMean/nSamples;
    }
 

    

    /*********************
            DICs
    **********************/
    SEXP DICMargResult;
    if(DICMarg){
      PROTECT(DICMargResult = allocMatrix(REALSXP, 4, 1)); nProtect++; //for Dbar, DbarOmega, pD, and DIC
    }
    
    SEXP DICUnmargResult;
    if(DICUnmarg){  
      PROTECT(DICUnmargResult = allocMatrix(REALSXP, 4, 1)); nProtect++; //for Dbar, DbarOmega, pD, and DIC
    }
    
    double *DMarg;
    double *DUnmarg;
    if(DICMarg){DMarg = (double *) R_alloc(nSamples, sizeof(double));}
    if(DICUnmarg){DUnmarg = (double *) R_alloc(nSamples, sizeof(double));}
  
    double DBarMarg, DBarMargOmega, pDMarg, DBarUnmarg, DBarUnmargOmega, pDUnmarg;
    DBarMarg=DBarMargOmega=pDMarg=DBarUnmarg=DBarUnmargOmega=pDUnmarg = 0.0;



    /*********************
        DICUnmarg DBar
    **********************/
    
    if(DICUnmarg){
      if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\tCalculating unmarginalized DIC\n");
	Rprintf("-------------------------------------------------\n");
        #ifdef Win32
	R_FlushConsole();
        #endif
      }
      
      status = 0;
      for(s = 0; s < nSamples; s++){
	
	if(!isPp){

	  if(nugget){

	    F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, &beta[s*p], &incOne, &zero, tmp_n, &incOne);
	    diff3(n, Y, tmp_n, &w[s*n], tmp_n1);
	    
	    DUnmarg[s] = n*log(tauSq[s])+1.0/tauSq[s]*F77_NAME(ddot)(&n, tmp_n1, &incOne, tmp_n1, &incOne); 
	
	  }else{//no such thing as unmarg so just give them the use marginalized

	    //make the correlation matrix
	    for(i = 0; i < nn; i++){
	      if(onePramPtr)
		(covModelObj->*cov1ParamPtr)(phi[s], C[i], obsD[i]);
	      else //i.e., 2 parameter matern
		(covModelObj->*cov2ParamPtr)(phi[s], nu[s], C[i], obsD[i]);
	    }
	    
	    F77_NAME(dscal)(&nn, &sigmaSq[s], C, &incOne);

	    //invert C and log det cov
	    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    
	    detCand = 0;
	    for(i = 0; i < n; i++) detCand += 2*log(C[i*n+i]);
	    
	    F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	    
	    //Y-XB
	    F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s*p], &incOne, &zero, tmp_n, &incOne);
	    F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);

	    F77_NAME(dsymv)(upper, &n, &one,  C, &n, tmp_n, &incOne, &zero, tmp_n1, &incOne);
	    DUnmarg[s] = detCand+F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n1, &incOne);
      
	  }


	}else{

	  //make the correlation matrix
	  for(i = 0; i < mm; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi[s], C_str[i], knotsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi[s], nu[s], C_str[i], knotsD[i]);
	  }
	  
	  for(i = 0; i < nm; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi[s], ct[i], obsKnotsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi[s], nu[s], ct[i], obsKnotsD[i]);
	  }
	  
	  //scale by sigma^2
	  F77_NAME(dscal)(&mm, &sigmaSq[s], C_str, &incOne);	
	  F77_NAME(dscal)(&nm, &sigmaSq[s], ct, &incOne);
	  
	  //invert C_str
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  //ct C^{*-1}
	  F77_NAME(dsymm)(rside, upper, &n, &m, &one, C_str, &m, ct, &n, &zero, tmp_nm, &n);
	  
	  F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, &beta[s*p], &incOne, &zero, tmp_n, &incOne);
	  diff3(n, Y, tmp_n, &w[s*n], tmp_n1);
	  
	  if(isModPp){

	    //ct C^{*-1} c
	    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, tmp_nn, &n);
	    
	  //\tild{\eps}
	    for(i = 0; i < n; i++){
	      
	      tmp = rnorm(1.0/(1.0/(sigmaSq[s]-tmp_nn[i*n+i]) + 1.0/tauSq[s])*1.0/tauSq[s]*tmp_n1[i], 
			  sqrt(1.0/(1.0/(sigmaSq[s]-tmp_nn[i*n+i]) + 1.0/tauSq[s])));
	      
	      
	      tildEpsMean[i] += tmp;
	      tmp_n1[i] = tmp_n1[i]-tmp;
	    }
	  }
	  
	  DUnmarg[s] = n*log(tauSq[s])+1.0/tauSq[s]*F77_NAME(ddot)(&n, tmp_n1, &incOne, tmp_n1, &incOne);

	}

	if(verbose){
	  if(status == 100){
	    Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
            #ifdef Win32
	    R_FlushConsole();
            #endif
	    status = 0;
	  }
	  status++;
	}
	R_CheckUserInterrupt();
      }
      

      DBarUnmarg = 0;
      for(i = 0; i < nSamples; i++)
	DBarUnmarg += DUnmarg[i];
      
      DBarUnmarg = DBarUnmarg/nSamples;
      
      if(isPp){
	if(isModPp){
	  for(i = 0; i < n; i++){
	    tildEpsMean[i] = tildEpsMean[i]/nSamples;
	  }
	}
      }
      
      if(!isPp){
	

	if(nugget){
	F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, betaMean, &incOne, &zero, tmp_n, &incOne);
	diff3(n, Y, tmp_n, wMeans, tmp_n1);
	
	DBarUnmargOmega = n*log(tauSqMean)+1.0/tauSqMean*F77_NAME(ddot)(&n, tmp_n1, &incOne, tmp_n1, &incOne); 

	}else{
	  
	  //make the correlation matrix
	  for(i = 0; i < nn; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phiMean, C[i], obsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phiMean, nuMean, C[i], obsD[i]);
	  }
	  
	  F77_NAME(dscal)(&nn, &sigmaSqMean, C, &incOne);
	  
	  //invert C and log det cov
	  F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  
	  detCand = 0;
	  for(i = 0; i < n; i++) detCand += 2*log(C[i*n+i]);
	  
	  F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  //Y-XB
	  F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, betaMean, &incOne, &zero, tmp_n, &incOne);
	  F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	  
	  F77_NAME(dsymv)(upper, &n, &one,  C, &n, tmp_n, &incOne, &zero, tmp_n1, &incOne);
	  DBarUnmargOmega = detCand+F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n1, &incOne);
	  
	}

      }else{
		
	for(i = 0; i < mm; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiMean, C_str[i], knotsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiMean, nuMean, C_str[i], knotsD[i]);
	}
	
	for(i = 0; i < nm; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiMean, ct[i], obsKnotsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiMean, nuMean, ct[i], obsKnotsD[i]);
	}
		
	//scale by sigma^2
	F77_NAME(dscal)(&mm, &sigmaSqMean, C_str, &incOne);	
	F77_NAME(dscal)(&nm, &sigmaSqMean, ct, &incOne);
	
	//invert C_str
	F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	
	//ct C^{*-1}
	F77_NAME(dsymm)(rside, upper, &n, &m, &one, C_str, &m, ct, &n, &zero, tmp_nm, &n);
	
	//ct C^{*-1} c
	F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, tmp_nn, &n);

	F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, betaMean, &incOne, &zero, tmp_n, &incOne);
	diff3(n, Y, tmp_n, wMeans, tmp_n1);

	if(isModPp){
	  for(i = 0; i < n; i++){
	    tmp_n1[i] = tmp_n1[i]-tildEpsMean[i];
	  }
	}

	DBarUnmargOmega = n*log(tauSqMean)+1.0/tauSqMean*F77_NAME(ddot)(&n, tmp_n1, &incOne, tmp_n1, &incOne); 
      }
      
    }//end Unmarg




    /*********************
        DICMarg DBar
    **********************/
    
    if(DICMarg){
      if(verbose){
	Rprintf("-------------------------------------------------\n");
	Rprintf("\tCalculating marginalized DIC\n");
	Rprintf("-------------------------------------------------\n");
        #ifdef Win32
	R_FlushConsole();
        #endif
      }
      
      status = 0;
      for(s = 0; s < nSamples; s++){
	
	if(!isPp){
	  
	  //make the correlation matrix
	  for(i = 0; i < nn; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi[s], C[i], obsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi[s], nu[s], C[i], obsD[i]);
	  }
	  
	  F77_NAME(dscal)(&nn, &sigmaSq[s], C, &incOne);
	  
	  if(nugget){
	    for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+tauSq[s];
	  }
	  
	  //invert C and log det cov
	  F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  
	  detCand = 0;
	  for(i = 0; i < n; i++) detCand += 2*log(C[i*n+i]);
	  
	  F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	}else{
	  
	  /******************************
           Predictive process
	  *******************************/
	  //make the correlation matrix
	  for(i = 0; i < mm; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi[s], C_str[i], knotsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi[s], nu[s], C_str[i], knotsD[i]);
	  }
	  
	  for(i = 0; i < nm; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi[s], ct[i], obsKnotsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi[s], nu[s], ct[i], obsKnotsD[i]);
	  }
	  
	  //scale by sigma^2
	  F77_NAME(dscal)(&mm, &sigmaSq[s], C_str, &incOne);	
	  F77_NAME(dscal)(&nm, &sigmaSq[s], ct, &incOne);
	  
	  /******************************
           Sherman-Woodbury-Morrison
	  *******************************/
	  if(!nugget) tauSq[s] = 1e-10;//ridge the matrix if no nugget model
	  //Unmodified predictive process
	  
	  if(!isModPp){  
	    tauSqInv = 1.0/tauSq[s];
	    negTauSqInv = -1.0*tauSqInv;
	    
	    //t(ct) 1/tau^2  ct = tmp_mm
	    F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &tauSqInv, ct, &n, ct, &n, &zero, tmp_mm, &m);
	    
	    //[C* + t(ct) 1/tau^2  ct]^{-1} = tmp_mm, and get the det on the way
	    F77_NAME(daxpy)(&mm, &one, C_str, &incOne, tmp_mm, &incOne);
	    
	    F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    
	    //get log det cov
	    detCand = 0;
	    for(i = 0; i < m; i++) detCand += 2.0*log(tmp_mm[i*m+i]);
	    
	    F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	    
	    //-1/tau^2 ct  tmp_mm = tmp_nm
	    F77_NAME(dsymm)(rside, upper, &n, &m, &negTauSqInv, tmp_mm, &m, ct, &n, &zero, tmp_nm, &n);
	    
	    //diag(1/tau^2) + tmp_nm ct 1/tau^2 = C
	    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &tauSqInv, tmp_nm, &n, ct, &n, &zero, C, &n);
	    for(i = 0; i < n; i++) C[i*n+i] = tauSqInv+C[i*n+i];
	    
	    //finish getting the log det cov 
	    detCand += n*log(tauSq[s]);
	    F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
	    for(i = 0; i < m; i++) detCand -= 2.0*log(C_str[i*m+i]);
	    
	  }else{//Modified predictive process
	    //make ct C
	    F77_NAME(dcopy)(&mm, C_str, &incOne, tmp_mm, &incOne);
	    F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    
	    detCandC_str = 0;
	    for(i = 0; i < m; i++) 
	      detCandC_str += 2.0*log(tmp_mm[i*m+i]);
	    
	    F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	    
	    F77_NAME(dsymm)(rside, upper, &n, &m, &one, tmp_mm, &m, ct, &n, &zero, tmp_nm, &n);
	    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, tmp_nn, &n);
	    
	    detCandE = 0;
	    for(i = 0; i < n; i++){ 
	      E[i] = tauSq[s]+sigmaSq[s]-tmp_nn[i*n+i];
	      Einv[i] = 1.0/E[i];
	      detCandE += 2.0*log(sqrt(E[i]));
	    }
	    
	    //make Einv Ct
	    //F77_NAME(dsymm)(lside, upper, &n, &m, &one, Einv, &n, ct, &n, &zero, tmp_nm, &n);
	    diagmm(n, m, Einv, ct, tmp_nm);
	    
	    //make C* + t(Ct) E.inv Ct
	    F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &one, ct, &n, tmp_nm, &n, &zero, tmp_mm, &m);
	    F77_NAME(daxpy)(&mm, &one, C_str, &incOne, tmp_mm, &incOne);
	    
	    //get log(|tmp_mm|) then tmp_mm^{-1}
	    detCand = 0.0;
	    F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    for(i = 0; i < m; i++) detCand += 2.0*log(tmp_mm[i*m+i]);
	    F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	    
	    detCand = detCandE+detCand-detCandC_str;
	    
	    //C = Einv - Einv ct (C* + t(ct) Einv ct)^{-1} t(ct) Einv
	    F77_NAME(dsymm)(rside, upper, &n, &m, &one, tmp_mm, &m, tmp_nm, &n, &zero, tmp_nm1, &n);
	    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm1, &n, tmp_nm, &n, &zero, C, &n);
	    
	    F77_NAME(dscal)(&nn, &negOne, C, &incOne);
	    for(i = 0; i < n; i++) 
	      C[i*n+i] = Einv[i]+C[i*n+i];
	  }
	}
	
	//Y-XB
	F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, &beta[s*p], &incOne, &zero, tmp_n, &incOne);
	F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	
	F77_NAME(dsymv)(upper, &n, &one,  C, &n, tmp_n, &incOne, &zero, tmp_n1, &incOne);
	DMarg[s] = detCand+F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n1, &incOne);
	
	if(verbose){
	  if(status == 100){
	    Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
            #ifdef Win32
	    R_FlushConsole();
            #endif
	    status = 0;
	  }
	  status++;
	}
	R_CheckUserInterrupt();
      }


      DBarMarg = 0;
      for(i = 0; i < nSamples; i++)
	DBarMarg += DMarg[i];
      
      DBarMarg = DBarMarg/nSamples;
      
      
      if(!isPp){
	
	//make the correlation matrix
	for(i = 0; i < nn; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiMean, C[i], obsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiMean, nuMean, C[i], obsD[i]);
	}
	
	F77_NAME(dscal)(&nn, &sigmaSqMean, C, &incOne);
	
	if(nugget){
	  for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+tauSqMean;
	}
	
	//invert C and log det cov
	F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	
	detCand = 0;
	for(i = 0; i < n; i++) detCand += 2*log(C[i*n+i]);
	
	F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	
      }else{
	
	/******************************
           Predictive process
	*******************************/
	//make the correlation matrix
	for(i = 0; i < mm; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiMean, C_str[i], knotsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiMean, nuMean, C_str[i], knotsD[i]);
	}
	
	for(i = 0; i < nm; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phiMean, ct[i], obsKnotsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phiMean, nuMean, ct[i], obsKnotsD[i]);
	}
	
	//scale by sigma^2
	F77_NAME(dscal)(&mm, &sigmaSqMean, C_str, &incOne);	
	F77_NAME(dscal)(&nm, &sigmaSqMean, ct, &incOne);
	
	/******************************
           Sherman-Woodbury-Morrison
	*******************************/
	if(!isModPp){  
	  tauSqInv = 1.0/tauSqMean;
	  negTauSqInv = -1.0*tauSqInv;
	  
	  //t(ct) 1/tau^2  ct = tmp_mm
	  F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &tauSqInv, ct, &n, ct, &n, &zero, tmp_mm, &m);
	  
	  //[C* + t(ct) 1/tau^2  ct]^{-1} = tmp_mm, and get the det on the way
	  F77_NAME(daxpy)(&mm, &one, C_str, &incOne, tmp_mm, &incOne);
	  
	  F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  
	  //get log det cov
	  detCand = 0;
	  for(i = 0; i < m; i++) detCand += 2.0*log(tmp_mm[i*m+i]);
	  
	  F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  //-1/tau^2 ct  tmp_mm = tmp_nm
	  F77_NAME(dsymm)(rside, upper, &n, &m, &negTauSqInv, tmp_mm, &m, ct, &n, &zero, tmp_nm, &n);
	  
	  //diag(1/tau^2) + tmp_nm ct 1/tau^2 = C
	  F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &tauSqInv, tmp_nm, &n, ct, &n, &zero, C, &n);
	  for(i = 0; i < n; i++) C[i*n+i] = tauSqInv+C[i*n+i];
	  
	  //finish getting the log det cov 
	  detCand += n*log(tauSqMean);
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
	  for(i = 0; i < m; i++) detCand -= 2.0*log(C_str[i*m+i]);
	  
	}else{//Modified predictive process
	  
	  //make ct C
	  F77_NAME(dcopy)(&mm, C_str, &incOne, tmp_mm, &incOne);
	  F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  
	  detCandC_str = 0;
	  for(i = 0; i < m; i++) 
	    detCandC_str += 2.0*log(tmp_mm[i*m+i]);
	  
	  F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  F77_NAME(dsymm)(rside, upper, &n, &m, &one, tmp_mm, &m, ct, &n, &zero, tmp_nm, &n);
	  F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, tmp_nn, &n);
	  
	  detCandE = 0;
	  for(i = 0; i < n; i++){ 
	    E[i] = tauSqMean+sigmaSqMean-tmp_nn[i*n+i];
	    Einv[i] = 1.0/E[i];
	    detCandE += 2.0*log(sqrt(E[i]));
	  }
	  
	  //make Einv Ct
	  //F77_NAME(dsymm)(lside, upper, &n, &m, &one, Einv, &n, ct, &n, &zero, tmp_nm, &n);
	  diagmm(n, m, Einv, ct, tmp_nm);
	  
	  //make C* + t(Ct) E.inv Ct
	  F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &one, ct, &n, tmp_nm, &n, &zero, tmp_mm, &m);
	  F77_NAME(daxpy)(&mm, &one, C_str, &incOne, tmp_mm, &incOne);
	  
	  //get log(|tmp_mm|) then tmp_mm^{-1}
	  detCand = 0.0;
	  F77_NAME(dpotrf)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  for(i = 0; i < m; i++) detCand += 2.0*log(tmp_mm[i*m+i]);
	  F77_NAME(dpotri)(upper, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  detCand = detCandE+detCand-detCandC_str;
	  
	  //C = Einv - Einv ct (C* + t(ct) Einv ct)^{-1} t(ct) Einv
	  F77_NAME(dsymm)(rside, upper, &n, &m, &one, tmp_mm, &m, tmp_nm, &n, &zero, tmp_nm1, &n);
	  F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm1, &n, tmp_nm, &n, &zero, C, &n);
	  
	  F77_NAME(dscal)(&nn, &negOne, C, &incOne);
	  for(i = 0; i < n; i++) 
	    C[i*n+i] = Einv[i]+C[i*n+i];
	}
      }
      
      //Y-XB
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, betaMean, &incOne, &zero, tmp_n, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
      
      F77_NAME(dsymv)(upper, &n, &one,  C, &n, tmp_n, &incOne, &zero, tmp_n1, &incOne);
      DBarMargOmega = detCand+F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n1, &incOne);
	
      
    }//end Marg
   

    /*********************
         return
    **********************/
    
    //make the result list object
    int nResultListObjs = 0;
    if(DICMarg) nResultListObjs++;
    if(DICUnmarg) nResultListObjs++;
    if(!spEffects  && DICUnmarg) nResultListObjs++;
    
    SEXP dicResults, dicResultsRowNames, dicResultsColNames, dicResultDimName;
    PROTECT(dicResultsRowNames = allocVector(VECSXP, 4)); nProtect++;
    PROTECT(dicResultsColNames = allocVector(VECSXP, 1)); nProtect++;
    SET_VECTOR_ELT(dicResultsRowNames, 0, mkChar("bar.D"));
    SET_VECTOR_ELT(dicResultsRowNames, 1, mkChar("D.bar.Omega"));
    SET_VECTOR_ELT(dicResultsRowNames, 2, mkChar("pD"));
    SET_VECTOR_ELT(dicResultsRowNames, 3, mkChar("DIC"));
    SET_VECTOR_ELT(dicResultsColNames, 0, mkChar("value"));
    PROTECT(dicResultDimName = allocVector(VECSXP, 2)); nProtect++;
    SET_VECTOR_ELT(dicResultDimName, 0, dicResultsRowNames);
    SET_VECTOR_ELT(dicResultDimName, 1, dicResultsColNames);
    
    
    SEXP result, resultNames;
    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //set result list elements
    int resultListIndx = 0;
    if(DICMarg){
      SET_VECTOR_ELT(result, resultListIndx, DICMargResult);
      REAL(DICMargResult)[0] = DBarMarg;
      REAL(DICMargResult)[1] = DBarMargOmega;
      REAL(DICMargResult)[2] = DBarMarg - DBarMargOmega;
      REAL(DICMargResult)[3] = DBarMarg + DBarMarg - DBarMargOmega; 
      setAttrib(DICMargResult, R_DimNamesSymbol, dicResultDimName);
      
      SET_VECTOR_ELT(resultNames, resultListIndx, mkChar("DIC.marg")); 
      resultListIndx++;
    }
    
    if(DICUnmarg){
      SET_VECTOR_ELT(result, resultListIndx, DICUnmargResult);
      REAL(DICUnmargResult)[0] = DBarUnmarg;
      REAL(DICUnmargResult)[1] = DBarUnmargOmega;
      REAL(DICUnmargResult)[2] = DBarUnmarg - DBarUnmargOmega;
      REAL(DICUnmargResult)[3] = DBarUnmarg + DBarUnmarg - DBarUnmargOmega; 
      setAttrib(DICUnmargResult, R_DimNamesSymbol, dicResultDimName);

      SET_VECTOR_ELT(resultNames, resultListIndx, mkChar("DIC.unmarg")); 
      resultListIndx++;
    }
    
    if(!spEffects  && DICUnmarg){//calculated it here
      //w observed
      SET_VECTOR_ELT(result, resultListIndx, w_r);
      SET_VECTOR_ELT(resultNames, resultListIndx, mkChar("sp.effects"));    
    }
    
    namesgets(result, resultNames);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    


























    }
}
