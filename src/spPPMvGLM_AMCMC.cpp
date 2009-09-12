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

  SEXP spPPMvGLM_AMCMC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,SEXP family_r,
		       SEXP isModPp_r, SEXP q_r, SEXP knotsD_r, SEXP coordsKnotsD_r,
		       SEXP betaPrior_r, SEXP betaNorm_r, SEXP KIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
		       SEXP phiStarting_r, SEXP AStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP w_strStarting_r,
		       SEXP phiTuning_r, SEXP ATuning_r, SEXP nuTuning_r , SEXP betaTuning_r, SEXP w_strTuning_r,
		       SEXP covModel_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r){


    /*****************************************
                Common variables
    *****************************************/
    int i,j,k,l,b,r,t,info,nProtect= 0;
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
    int pp = p*p;
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int q = INTEGER(q_r)[0];
    int nn = n*n;
    int mm = m*m;
    int nm = n*m;
    int nq = n*q;
    int qq = q*q;
    int qm = q*m;
    int qmqm = qm*qm;
    int nLTr = m*(m-1)/2+m;

    string family = CHAR(STRING_ELT(family_r,0));

    //covariance model
    string covModel = CHAR(STRING_ELT(covModel_r,0));

    //if predictive process
    bool isModPp = static_cast<bool>(INTEGER(isModPp_r)[0]);

    double *knotsD = REAL(knotsD_r);
    double *coordsKnotsD = REAL(coordsKnotsD_r);

    //priors and starting
    string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));

    double *betaMu = NULL;
    double *betaSd = NULL;
    
    if(betaPrior == "normal"){
      betaMu = REAL(VECTOR_ELT(betaNorm_r, 0)); 
      betaSd = REAL(VECTOR_ELT(betaNorm_r, 1));
    }
    
    double *phiUnif = REAL(phiUnif_r);
    double KIW_df = REAL(VECTOR_ELT(KIW_r, 0))[0]; double *KIW_S = REAL(VECTOR_ELT(KIW_r, 1));

    double *phiStarting = REAL(phiStarting_r);
    double *AStarting = REAL(AStarting_r);
    double *betaStarting = REAL(betaStarting_r);
    double *w_strStarting = REAL(w_strStarting_r);

    //if matern
    double *nuUnif = NULL;
    double *nuStarting = NULL;

    if(covModel == "matern"){
      nuUnif = REAL(nuUnif_r);
      nuStarting = REAL(nuStarting_r);
    }

    int nBatch = INTEGER(nBatch_r)[0];
    int batchLength = INTEGER(batchLength_r)[0];
    double acceptRate  = REAL(acceptRate_r)[0];
    int nSamples = nBatch*batchLength;

    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    double *A = (double *) R_alloc(mm, sizeof(double));

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      if(isModPp)
	Rprintf("Using modified predictive process with %i knots.\n\n", q);
      else
	Rprintf("Using non-modified predictive process with %i knots.\n\n", q);
      
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);

      Rprintf("Priors and hyperpriors:\n");
  
      if(betaPrior == "flat"){
	Rprintf("\tbeta flat.\n");
      }else{
	Rprintf("\tbeta normal:\n");
	Rprintf("\t\tmu:"); printVec(betaMu, p);
	Rprintf("\t\tsd:"); printVec(betaSd, p);Rprintf("\n");
      }
      Rprintf("\n");
 
      Rprintf("\tK IW hyperpriors df=%.5f\n", KIW_df);
      printMtrx(KIW_S, m, m);
      Rprintf("\n"); 

      Rprintf("\tphi Unif hyperpriors\n");
      Rprintf("\t");   
      for(i = 0; i < m; i++){
	Rprintf("(%.5f, %.5f) ", phiUnif[i*2], phiUnif[i*2+1]);
      }
      Rprintf("\n\n");   

      if(covModel == "matern"){
	Rprintf("\tnu Unif hyperpriors\n");
	Rprintf("\t");
	for(i = 0; i < m; i++){
	  Rprintf("(%.5f, %.5f) ", nuUnif[i*2], nuUnif[i*2+1]);
	}
	Rprintf("\n\n");   
      }

      Rprintf("Adaptive Metropolis with target acceptance rate: %.1f%\n", 100*acceptRate);
      Rprintf("\n");   

      Rprintf("Metropolis starting values:\n");
  
      Rprintf("\tbeta starting:\n");
      Rprintf("\t"); printVec(betaStarting, p);
      Rprintf("\n"); 

      covExpand(AStarting, A, m);
      Rprintf("\tA starting\n");
      printMtrx(A, m, m);
      Rprintf("\n"); 

      Rprintf("\tphi starting\n");
      Rprintf("\t"); printVec(phiStarting, m);
      Rprintf("\n");   

      if(covModel == "matern"){
	Rprintf("\t"); Rprintf("\tnu starting\n");
	printVec(nuStarting, m);
      }
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
         Set-up MCMC sample matrices etc.
    *****************************************/
    //spatial parameters
    int nParams, betaIndx, AIndx, LIndx, phiIndx, nuIndx;

    if(covModel != "matern"){
      nParams = p+nLTr+m;//A, phi
      betaIndx = 0; AIndx = betaIndx+p; phiIndx = AIndx+nLTr;
    }else {
      nParams = p+nLTr+2*m;//A, phi, nu
      betaIndx = 0; AIndx = betaIndx+p; phiIndx = AIndx+nLTr, nuIndx = phiIndx+m;
    }
    
    double *spParams = (double *) R_alloc(nParams, sizeof(double));
    
    //set starting
    F77_NAME(dcopy)(&p, betaStarting, &incOne, &spParams[betaIndx], &incOne);

    covTrans(AStarting, &spParams[AIndx], m);

    for(i = 0; i < m; i++){
      spParams[phiIndx+i] = logit(phiStarting[i], phiUnif[i*2], phiUnif[i*2+1]);
    }

    if(covModel == "matern"){
      for(i = 0; i < m; i++){
	spParams[nuIndx+i] = logit(nuStarting[i], nuUnif[i*2], nuUnif[i*2+1]);
      }
    }

    double *w_str = (double *) R_alloc(qm, sizeof(double));
    F77_NAME(dcopy)(&qm, w_strStarting, &incOne, w_str, &incOne);
   
    //set tuning
    double *spTuning = (double *) R_alloc(nParams, sizeof(double));
    
    F77_NAME(dcopy)(&p, REAL(betaTuning_r), &incOne, &spTuning[betaIndx], &incOne);

    for(i = 0; i < nLTr; i++)
      spTuning[AIndx+i] = REAL(ATuning_r)[i];
    
    for(i = 0; i < m; i++){
      spTuning[phiIndx+i] = REAL(phiTuning_r)[i];
      
      if(covModel == "matern"){
	spTuning[nuIndx+i] = REAL(nuTuning_r)[i];
      }
    }
    
    for(i = 0; i < nParams; i++)
      spTuning[i] = log(spTuning[i]);

    double *w_strTuning = (double *) R_alloc(qm, sizeof(double));
    F77_NAME(dcopy)(&qm, REAL(w_strTuning_r), &incOne, w_strTuning, &incOne);

    for(i = 0; i < qm; i++)
      w_strTuning[i] = log(w_strTuning[i]);

    //samples and random effects
    SEXP w_r, w_str_r, samples_r, accept_r, accept_w_str_r, tuning_r, tuning_w_str_r;
    PROTECT(w_r = allocMatrix(REALSXP, nm, nSamples)); nProtect++; 
    PROTECT(w_str_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 
    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    PROTECT(accept_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++; //just to monitor acceptance rate
    PROTECT(accept_w_str_r = allocMatrix(REALSXP, qm, nBatch)); nProtect++; //just to monitor acceptance rate
    PROTECT(tuning_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++; //just to monitor acceptance rate
    PROTECT(tuning_w_str_r = allocMatrix(REALSXP, qm, nBatch)); nProtect++; //just to monitor acceptance rate

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0, logDetK, SKtrace, spParamsjCurrent, w_strjCurrent;
 
    double *accept = (double *) R_alloc(nParams, sizeof(double));
    double *accept_w_str = (double *) R_alloc(qm, sizeof(double));

    double *ct = (double *) R_alloc(nm*qm, sizeof(double));
    double *C_str = (double *) R_alloc(qmqm, sizeof(double));

    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm1 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *tmp_qm = (double *) R_alloc(qm, sizeof(double));
    double *tmp_nmqm = (double *) R_alloc(nm*qm, sizeof(double));

    double *beta = (double *) R_alloc(p, sizeof(double));
    double *K = (double *) R_alloc(mm, sizeof(double));
    double *theta = (double *) R_alloc(mm, sizeof(double));
    double *w = (double *) R_alloc(nm, sizeof(double)); zeros(w, nm);

    double logMHRatio;

    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    logPostCurrent = R_NegInf;
    
    GetRNGstate();
    for(b = 0, s = 0; b < nBatch; b++){
      
      for(i = 0; i < batchLength; i++, s++){
	
	for(j = 0; j < nParams; j++){
	  
	  //propose
	  spParamsjCurrent = spParams[j];
	  spParams[j] = rnorm(spParamsjCurrent, exp(spTuning[j]));
	  
	  //extract and transform
	  F77_NAME(dcopy)(&p, &spParams[betaIndx], &incOne, beta, &incOne);
	  
	  covTransInvExpand(&spParams[AIndx], A, m);
	  
	  for(k = 0; k < m; k++){
	    theta[k] = logitInv(spParams[phiIndx+k], phiUnif[k*2], phiUnif[k*2+1]);
	    
	    if(covModel == "matern"){
	      theta[m+k] = logitInv(spParams[nuIndx+k], nuUnif[k*2], nuUnif[k*2+1]);
	    }
	  }
	  
	  //K = A'A
	  F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, A, &m, A, &m, &zero, K, &m);
	  
	  //make ct
	  for(r = 0; r < n; r++){
	    for(t = 0; t < q; t++){
	      
	      zeros(tmp_mm, mm);
	      
	      for(k = 0; k < m; k++){
		if(nPramPtr == 1)
		  (covModelObj->*cov1ParamPtr)(theta[k], tmp_mm[k*m+k], coordsKnotsD[t*n+r]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(theta[k], theta[m+k], tmp_mm[k*m+k], coordsKnotsD[t*n+r]);
	      }
	      
	      F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, A, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, A, &m, &zero, tmp_mm, &m);
	      
	      for(k = 0; k < m; k++){
		for(l = 0; l < m; l++){
		  ct[((t*m+l)*nm)+(r*m+k)] = tmp_mm[l*m+k];
		  tmp_mm[l*m+k] = 0.0; //zero out
		}
	      }
	    }
	  }
	  
	  //
	  //make C_str
	  //
	  for(r = 0; r < q; r++){
	    for(t = 0; t < q; t++){
	      
	      zeros(tmp_mm, mm);
	      
	      for(k = 0; k < m; k++){
		if(nPramPtr == 1)
		  (covModelObj->*cov1ParamPtr)(theta[k], tmp_mm[k*m+k], knotsD[t*q+r]);
		else //i.e., 2 parameter matern
		  (covModelObj->*cov2ParamPtr)(theta[k], theta[m+k], tmp_mm[k*m+k], knotsD[t*q+r]);
	      }
	      
	      F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, A, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, A, &m, &zero, tmp_mm, &m);
	      
	      for(k = 0; k < m; k++){
		for(l = 0; l < m; l++){
		  C_str[((t*m+l)*qm)+(r*m+k)] = tmp_mm[l*m+k];
		  tmp_mm[l*m+k] = 0.0; //zero out
		}
	      }
	    }
	  }
	  
	  detCand = 0.0;
	  F77_NAME(dpotrf)(lower, &qm, C_str, &qm, &info); if(info != 0){error("mvPPCovInvDet: Cholesky failed (2)\n");}
	  for(k = 0; k < qm; k++) detCand += 2.0*log(C_str[k*qm+k]);
	  F77_NAME(dpotri)(lower, &qm, C_str, &qm, &info); if(info != 0){error("mvPPCovInvDet: Cholesky failed (3)\n");}
	  
	  F77_NAME(dsymm)(rside, lower, &nm, &qm, &one, C_str, &qm, ct, &nm, &zero, tmp_nmqm, &nm);
	  
	  //make \tild{w}
	  F77_NAME(dgemv)(ntran, &nm, &qm, &one, tmp_nmqm, &nm, w_str, &incOne, &zero, w, &incOne);
	  
	  //
	  //Likelihood with Jacobian   
	  //
	  logPostCand = 0.0;
	  
	  if(betaPrior == "normal"){
	    for(k = 0; k < p; k++){
	      logPostCand += dnorm(beta[k], betaMu[k], betaSd[k], 1);
	    }
	  }      
	  
	  //
	  //Jacobian and IW priors for K = A'A
	  //
	  
	  //AtA prior with jacob.
	  logDetK = 0.0;
	  SKtrace = 0.0;
	  
	  for(k = 0; k < m; k++) logDetK += 2*log(A[k*m+k]);
	  
	  //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	  for(k = 0; k < m; k++) logPostCand += (m-k)*log(A[k*m+k])+log(A[k*m+k]);
	  
	  //get S*K^-1, already have the chol of K (i.e., A)
	  F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){cout << "c++ error: Cand A Cholesky inverse failed\n" << endl;}
	  F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
	  for(k = 0; k < m; k++){SKtrace += tmp_mm[k*m+k];}
	  logPostCand += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;
	  
	  for(k = 0; k < m; k++){
	    logPostCand += log(theta[k] - phiUnif[k*2]) + log(phiUnif[k*2+1] - theta[k]); 
	    
	    if(covModel == "matern"){
	      logPostCand += log(theta[m+k] - nuUnif[k*2]) + log(nuUnif[k*2+1] - theta[m+k]);  
	    }
	  }
	  
	  F77_NAME(dgemv)(ntran, &nm, &p, &one, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
	  
	  if(family == "binomial"){
	    logPostCand += logit_logpost(nm, Y, tmp_nm, w);
	  }else if(family == "poisson"){
	    logPostCand += poisson_logpost(nm, Y, tmp_nm, w);
	  }else{
	    error("c++ error: family misspecification in spGLM\n");
	  }
	  
	  //(-1/2) * tmp_n` *  C^-1 * tmp_n
	  F77_NAME(dsymv)(lower, &qm, &one,  C_str, &qm, w_str, &incOne, &zero, tmp_qm, &incOne);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&qm, w_str, &incOne, tmp_qm, &incOne);
	  
	  //
	  //MH accept/reject	
	  //      
	  
	  //MH ratio with adjustment
	  logMHRatio = logPostCand - logPostCurrent;

	  if(runif(0.0,1.0) <= exp(logMHRatio)){
	    logPostCurrent = logPostCand;
	    accept[j]++;
	  }else{
	    spParams[j] = spParamsjCurrent;
	  }

	}//end spParams


	/************************/
	//update the sp. matrices
	/************************/
	//extract and transform
	F77_NAME(dcopy)(&p, &spParams[betaIndx], &incOne, beta, &incOne);
	
	covTransInvExpand(&spParams[AIndx], A, m);
	
	for(k = 0; k < m; k++){
	  theta[k] = logitInv(spParams[phiIndx+k], phiUnif[k*2], phiUnif[k*2+1]);
	  
	  if(covModel == "matern"){
	    theta[m+k] = logitInv(spParams[nuIndx+k], nuUnif[k*2], nuUnif[k*2+1]);
	  }
	}
	
	//make ct
	for(r = 0; r < n; r++){
	  for(t = 0; t < q; t++){
	    
	    zeros(tmp_mm, mm);
	    
	    for(k = 0; k < m; k++){
	      if(nPramPtr == 1)
		(covModelObj->*cov1ParamPtr)(theta[k], tmp_mm[k*m+k], coordsKnotsD[t*n+r]);
	      else //i.e., 2 parameter matern
		(covModelObj->*cov2ParamPtr)(theta[k], theta[m+k], tmp_mm[k*m+k], coordsKnotsD[t*n+r]);
	    }
	    
	    F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, A, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	    F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, A, &m, &zero, tmp_mm, &m);
	    
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
		ct[((t*m+l)*nm)+(r*m+k)] = tmp_mm[l*m+k];
		tmp_mm[l*m+k] = 0.0; //zero out
	      }
	    }
	  }
	}
	
	//
	//make C_str
	//
	for(r = 0; r < q; r++){
	  for(t = 0; t < q; t++){
	    
	    zeros(tmp_mm, mm);
	    
	    for(k = 0; k < m; k++){
	      if(nPramPtr == 1)
		(covModelObj->*cov1ParamPtr)(theta[k], tmp_mm[k*m+k], knotsD[t*q+r]);
	      else //i.e., 2 parameter matern
		(covModelObj->*cov2ParamPtr)(theta[k], theta[m+k], tmp_mm[k*m+k], knotsD[t*q+r]);
	    }
	    
	    F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, A, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	    F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, A, &m, &zero, tmp_mm, &m);
	    
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
		C_str[((t*m+l)*qm)+(r*m+k)] = tmp_mm[l*m+k];
		tmp_mm[l*m+k] = 0.0; //zero out
	      }
	    }
	  }
	}
	  
	detCand = 0.0;
	F77_NAME(dpotrf)(lower, &qm, C_str, &qm, &info); if(info != 0){error("mvPPCovInvDet: Cholesky failed (2)\n");}
	for(k = 0; k < qm; k++) detCand += 2.0*log(C_str[k*qm+k]);
	F77_NAME(dpotri)(lower, &qm, C_str, &qm, &info); if(info != 0){error("mvPPCovInvDet: Cholesky failed (3)\n");}
	
	F77_NAME(dsymm)(rside, lower, &nm, &qm, &one, C_str, &qm, ct, &nm, &zero, tmp_nmqm, &nm);
	
	for(j = 0; j < qm; j++){
	  
	  //propose
	  w_strjCurrent = w_str[j];
	  w_str[j] = rnorm(w_strjCurrent, exp(w_strTuning[j]));
	  
	  //make \tild{w}
	  F77_NAME(dgemv)(ntran, &nm, &qm, &one, tmp_nmqm, &nm, w_str, &incOne, &zero, w, &incOne);
	  
	  //
	  //Likelihood with Jacobian   
	  //
	  logPostCand = 0.0;
	  
	  if(betaPrior == "normal"){
	    for(k = 0; k < p; k++){
	      logPostCand += dnorm(beta[k], betaMu[k], betaSd[k], 1);
	    }
	  }      
	  
	  //
	  //Jacobian and IW priors for K = A'A
	  //

	  covTransInvExpand(&spParams[AIndx], A, m);	  

	  //AtA prior with jacob.
	  logDetK = 0.0;
	  SKtrace = 0.0;
	  
	  for(k = 0; k < m; k++) logDetK += 2*log(A[k*m+k]);
	  
	  //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	  for(k = 0; k < m; k++) logPostCand += (m-k)*log(A[k*m+k])+log(A[k*m+k]);
	  
	  //get S*K^-1, already have the chol of K (i.e., A)
	  F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){cout << "c++ error: Cand A Cholesky inverse failed\n" << endl;}
	  F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
	  for(k = 0; k < m; k++){SKtrace += tmp_mm[k*m+k];}
	  logPostCand += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;
	  
	  for(k = 0; k < m; k++){
	    logPostCand += log(theta[k] - phiUnif[k*2]) + log(phiUnif[k*2+1] - theta[k]); 
	    
	    if(covModel == "matern"){
	      logPostCand += log(theta[m+k] - nuUnif[k*2]) + log(nuUnif[k*2+1] - theta[m+k]);  
	    }
	  }
	  
	  F77_NAME(dgemv)(ntran, &nm, &p, &one, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
	  
	  if(family == "binomial"){
	    logPostCand += logit_logpost(nm, Y, tmp_nm, w);
	  }else if(family == "poisson"){
	    logPostCand += poisson_logpost(nm, Y, tmp_nm, w);
	  }else{
	    error("c++ error: family misspecification in spGLM\n");
	  }
	  
	  //(-1/2) * tmp_n` *  C^-1 * tmp_n
	  F77_NAME(dsymv)(lower, &qm, &one,  C_str, &qm, w_str, &incOne, &zero, tmp_qm, &incOne);
	  logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&qm, w_str, &incOne, tmp_qm, &incOne);

	  //
	  //MH accept/reject	
	  //      
	  
	  //MH ratio with adjustment
	  logMHRatio = logPostCand - logPostCurrent;
	  
	  if(runif(0.0,1.0) <= exp(logMHRatio)){
	    logPostCurrent = logPostCand;
	    accept_w_str[j]++;
	  }else{
	    w_str[j] = w_strjCurrent;
	  }
	  
	}//end w
	

	/******************************
               Save samples
	*******************************/
	F77_NAME(dcopy)(&nParams, spParams, &incOne, &REAL(samples_r)[s*nParams], &incOne);
	F77_NAME(dcopy)(&nm, w, &incOne, &REAL(w_r)[s*nm], &incOne);
	F77_NAME(dcopy)(&qm, w_str, &incOne, &REAL(w_str_r)[s*qm], &incOne);
	
	R_CheckUserInterrupt();
      }//end batch
      
      //adjust tuning
      for(j = 0; j < nParams; j++){
	REAL(accept_r)[b*nParams+j] = accept[j]/batchLength;
	REAL(tuning_r)[b*nParams+j] = spTuning[j];
	
	if(accept[j]/batchLength > acceptRate){
	  spTuning[j] += min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}else{
	  spTuning[j] -= min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}
	accept[j] = 0.0;
      }
      
      for(j = 0; j < qm; j++){
	REAL(accept_w_str_r)[b*qm+j] = accept_w_str[j]/batchLength;
	REAL(tuning_w_str_r)[b*qm+j] = w_strTuning[j];
	
	if(accept_w_str[j]/batchLength > acceptRate){
	  w_strTuning[j] += min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}else{
	  w_strTuning[j] -= min(0.01, 1.0/sqrt(static_cast<double>(b)));
	}
	accept_w_str[j] = 0.0;
      }
    
      //report
      if(verbose){
	if(status == nReport){
	  Rprintf("Batch: %i of %i, %3.2f%%\n", b, nBatch, 100.0*b/nBatch);
	  Rprintf("\tparameter\tacceptance\ttuning\n");	  
	  for(j = 0; j < p; j++){
	    Rprintf("\tbeta[%i]\t\t%3.1f%\t\t%1.5f\n", j, 100.0*REAL(accept_r)[b*nParams+betaIndx+j], exp(spTuning[betaIndx+j]));
	  }
	  
	  covTransInvExpand(&spParams[AIndx], A, m);
	  
	  for(k = 0, i = 0; k < m; k++){
	    for(j = k; j < m; j++, i++){
	      Rprintf("\tA[%i,%i]\t\t%3.1f%\t\t%1.5f\n",j, k, 100.0*REAL(accept_r)[b*nParams+AIndx+j], exp(spTuning[AIndx+i]));
	    }
	  }
	  
	  for(j = 0; j < m; j++){
	    Rprintf("\tphi[%i]\t\t%3.1f%\t\t%1.5f\n", j, 100.0*REAL(accept_r)[b*nParams+phiIndx+j], exp(spTuning[phiIndx+j]));	
	    if(covModel == "matern")
	      Rprintf("\tnu[%i]\t\t%3.1f%\t\t%1.5f\n", j, 100.0*REAL(accept_r)[b*nParams+nuIndx+j], exp(spTuning[nuIndx+j]));
	  } 
	  
	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      
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
      
      covTransInv(&REAL(samples_r)[s*nParams+AIndx], &REAL(samples_r)[s*nParams+AIndx], m);
      
      for(i = 0; i < m; i++){
    	REAL(samples_r)[s*nParams+phiIndx+i] = logitInv(REAL(samples_r)[s*nParams+phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
    	if(covModel == "matern"){
    	  REAL(samples_r)[s*nParams+nuIndx+i] = logitInv(REAL(samples_r)[s*nParams+nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
    	}
      }
    }   
    
    //make return object
    SEXP result, resultNames;
    
    int nResultListObjs = 7;
    
    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result, 0, samples_r);
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));

    SET_VECTOR_ELT(result, 2, accept_w_str_r);
    SET_VECTOR_ELT(resultNames, 2, mkChar("acceptance.w.str"));
    
    SET_VECTOR_ELT(result, 3, w_r);
    SET_VECTOR_ELT(resultNames, 3, mkChar("sp.effects"));

    SET_VECTOR_ELT(result, 4, w_str_r);
    SET_VECTOR_ELT(resultNames, 4, mkChar("sp.effects.knots"));

    SET_VECTOR_ELT(result, 5, tuning_r);
    SET_VECTOR_ELT(resultNames, 5, mkChar("tuning"));

    SET_VECTOR_ELT(result, 6, tuning_w_str_r);
    SET_VECTOR_ELT(resultNames, 6, mkChar("tuning.w.str"));
  
    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    
  }
}
