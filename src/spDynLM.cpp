#include <algorithm>
#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP spDynLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP Nt_r, SEXP coordsD_r,
	       SEXP beta0Norm_r, SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r, SEXP sigmaEtaIW_r, 
	       SEXP betaStarting_r, SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r, SEXP sigmaEtaStarting_r, 
	       SEXP phiTuning_r, SEXP nuTuning_r, SEXP covModel_r, SEXP nSamples_r, SEXP missing_r, SEXP getFitted_r, SEXP verbose_r, SEXP nReport_r){
    
    /*****************************************
                Common variables
    *****************************************/
    int i, j, s, t, info, nProtect=0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    const double one = 1.0;
    const double two = 2.0;
    const double zero = 0.0;
    const int incOne = 1;

    /*****************************************
                     Set-up
    *****************************************/
    double *Y = REAL(Y_r); //n x Nt
    double *X = REAL(X_r); //p x Nt*n
    int p = INTEGER(p_r)[0];
    int pp = p*p;
    int n = INTEGER(n_r)[0];
    int Nt = INTEGER(Nt_r)[0];
    int Ntp = Nt*p;
    int nn = n*n;
    int Ntn = Nt*n;
         
    double *coordsD = REAL(coordsD_r);

    std::string covModel = CHAR(STRING_ELT(covModel_r,0));

    //priors  
    double *m0 = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(beta0Norm_r, 0)), &incOne, m0, &incOne);
    
    double *sigma0 = (double *) R_alloc(pp, sizeof(double)); 
    F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(beta0Norm_r, 1)), &incOne, sigma0, &incOne);
    
    double *sigmaSqIG = REAL(sigmaSqIG_r); 
    double *phiUnif = REAL(phiUnif_r); 
    double *tauSqIG = REAL(tauSqIG_r); 
 
    //matern
    double *nuUnif = NULL;
    if(covModel == "matern"){
      nuUnif = REAL(nuUnif_r);
    }
  
    double *SigmaEtaIW_S = (double *) R_alloc(pp, sizeof(double));
    double SigmaEtaIW_df = REAL(VECTOR_ELT(sigmaEtaIW_r, 0))[0]; SigmaEtaIW_S = REAL(VECTOR_ELT(sigmaEtaIW_r, 1));

    int nSamples = INTEGER(nSamples_r)[0];
    int *missing = INTEGER(missing_r);
    bool getFitted = static_cast<bool>(INTEGER(getFitted_r)[0]);
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    int nMissing = 0;
    bool anyMissing = false;
    for(i = 0; i < Ntn; i++){
      if(missing[i] == 1){
	nMissing++;
	anyMissing = true;
      }
    }
  
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations in %i time steps.\n\n", n, Nt);
      Rprintf("Number of missing observations %i.\n\n", nMissing);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
   
      Rprintf("Priors and hyperpriors:\n");
      
      Rprintf("\tbeta normal:\n");
      Rprintf("\tm_0:"); printVec(m0, p);
      Rprintf("\tSigma_0:\n"); printMtrx(sigma0, p, p);
      Rprintf("\n");
      
      for(i = 0, j=1; i < Nt; i++, j++){
	Rprintf("\tsigma.sq_t=%i IG hyperpriors shape=%.5f and scale=%.5f\n", j, sigmaSqIG[i*2], sigmaSqIG[i*2+1]);
	
	Rprintf("\ttau.sq_t=%i IG hyperpriors shape=%.5f and scale=%.5f\n", j, tauSqIG[i*2], tauSqIG[i*2+1]); 
	
	Rprintf("\tphi_t=%i Unif hyperpriors a=%.5f and b=%.5f\n", j, phiUnif[i*2], phiUnif[i*2+1]);
	if(covModel == "matern"){
	  Rprintf("\tnu_t=%i Unif hyperpriors a=%.5f and b=%.5f\n", j, nuUnif[i*2], nuUnif[i*2+1]);	  
	}
	Rprintf("\t---\n");
      }
    } 

    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/ 
    //parameters
    int nTheta, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;

    if(covModel != "matern"){
      nTheta = 3;//sigma^2, tau^2, phi
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; 
    }else{
      nTheta = 4;//sigma^2, tau^2, phi, nu
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;//sigma^2, tau^2, phi, nu
    }
    
    int NtnTheta = Nt*nTheta;

    double *beta0 = (double *) R_alloc(p, sizeof(double)); zeros(beta0, p);
    double *theta = (double *) R_alloc(NtnTheta, sizeof(double)); 
    double *beta = (double *) R_alloc(Ntp, sizeof(double));
    double *u = (double *) R_alloc(Ntn, sizeof(double)); zeros(u, Ntn);
    double *sigmaEta = (double *) R_alloc(pp, sizeof(double)); 
  
    //starting
    for(i = 0; i < Nt; i++){

      F77_NAME(dcopy)(&p, &REAL(betaStarting_r)[p*i], &incOne, &beta[p*i], &incOne); //p x Nt

      theta[nTheta*i+sigmaSqIndx] = REAL(sigmaSqStarting_r)[i]; //nTheta x Nt
      
      theta[nTheta*i+tauSqIndx] = REAL(tauSqStarting_r)[i];
           
      theta[nTheta*i+phiIndx] = logit(REAL(phiStarting_r)[i], phiUnif[i*2], phiUnif[i*2+1]);
      
      if(covModel == "matern"){
	theta[nTheta*i+nuIndx] = logit(REAL(nuStarting_r)[i], nuUnif[i*2], nuUnif[i*2+1]);
      }
    }

    F77_NAME(dcopy)(&pp, REAL(sigmaEtaStarting_r), &incOne, sigmaEta, &incOne);

    //tuning 
    double *phiTuning = (double *) R_alloc(Nt, sizeof(double));
    double *nuTuning = NULL;
    if(covModel == "matern"){
      nuTuning = (double *) R_alloc(Nt, sizeof(double));
    }

    for(i = 0; i < Nt; i++){
      phiTuning[i] = REAL(phiTuning_r)[i];
      
      if(covModel == "matern"){
	nuTuning[i] = REAL(nuTuning_r)[i];
      }
    }

    //return stuff  
    SEXP thetaSamples_r, beta0Samples_r, betaSamples_r, uSamples_r, sigmaEtaSamples_r, ySamples_r;
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta*Nt, nSamples)); nProtect++; 
    PROTECT(beta0Samples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++; 
    PROTECT(betaSamples_r = allocMatrix(REALSXP, p*Nt, nSamples)); nProtect++; 
    PROTECT(uSamples_r = allocMatrix(REALSXP, n*Nt, nSamples)); nProtect++; 
    PROTECT(sigmaEtaSamples_r = allocMatrix(REALSXP, pp, nSamples)); nProtect++; 

    double *thetaSamples = REAL(thetaSamples_r);
    double *beta0Samples = REAL(beta0Samples_r);
    double *betaSamples = REAL(betaSamples_r);
    double *uSamples = REAL(uSamples_r);
    double *sigmaEtaSamples = REAL(sigmaEtaSamples_r);
    double *ySamples = NULL;

    if(getFitted){
      PROTECT(ySamples_r = allocMatrix(REALSXP, n*Nt, nSamples)); nProtect++; 
      ySamples = REAL(ySamples_r);
    }
    
    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status = 1, accept = 0, batchAccept=0;
    double logMHRatio, logPostCand, logPost, logDetCand, logDet;
  
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_pp2 = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double tmp;
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_n2 = (double *) R_alloc(n, sizeof(double));
    double *tmp_n3 = (double *) R_alloc(n, sizeof(double));
    double *C = (double *) R_alloc(nn, sizeof(double));
    double *C2 = (double *) R_alloc(nn, sizeof(double));
    //double *C3 = (double *) R_alloc(nn, sizeof(double));
    double *gamma = (double *) R_alloc(3, sizeof(double)); //sigma^2, phi, nu

    //Sigma_0^{-1}
    F77_NAME(dpotrf)(lower, &p, sigma0, &p, &info); if(info != 0){error("c++ error: dpotrf failed 0\n");}
    F77_NAME(dpotri)(lower, &p, sigma0, &p, &info); if(info != 0){error("c++ error: dpotri failed 0.5\n");}

    //X'_tX_t
    double *XX = (double *) R_alloc(Nt*pp, sizeof(double));
    for(t = 0; t < Nt; t++){
      F77_NAME(dgemm)(ntran, ytran, &p, &p, &n, &one, &X[t*n*p], &p, &X[t*n*p], &p, &zero, &XX[t*pp], &p); 
    }
    
    //SigmaEtaIW_S is Upsilon^{-1}
    F77_NAME(dpotrf)(lower, &p, SigmaEtaIW_S, &p, &info); if(info != 0){error("c++ error: dpotrf failed 1.1\n");}
    F77_NAME(dpotri)(lower, &p, SigmaEtaIW_S, &p, &info); if(info != 0){error("c++ error: dpotri failed 2.01\n");}

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

      //Sigma_eta^{-1}
      F77_NAME(dpotrf)(lower, &p, sigmaEta, &p, &info); if(info != 0){error("c++ error: dpotrf1 failed\n");}
      F77_NAME(dpotri)(lower, &p, sigmaEta, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}

      //update Beta_0
      F77_NAME(dcopy)(&pp, sigmaEta, &incOne, tmp_pp, &incOne);
      F77_NAME(daxpy)(&pp, &one, sigma0, &incOne, tmp_pp, &incOne);

      //Sigma_Beta0
      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf2 failed\n");}
      F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}

      //mu_Beta0
      F77_NAME(dsymv)(lower, &p, &one, sigma0, &p, m0, &incOne, &zero, tmp_p, &incOne);
      F77_NAME(dsymv)(lower, &p, &one, sigmaEta, &p, beta, &incOne, &one, tmp_p, &incOne);//first column of beta is Beta_1
      F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &incOne, &zero, tmp_p2, &incOne); 

     F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf3 failed\n");}
     mvrnorm(beta0, tmp_p2, tmp_pp, p);

     F77_NAME(dcopy)(&p, beta0, &incOne, &beta0Samples[s*p], &incOne);
     
     ////////////////////
     //update y_t
     ////////////////////
     if(anyMissing || getFitted){
       for(t = 0, j = 0; t < Nt; t++){
	 F77_NAME(dgemv)(ytran, &p, &n, &one, &X[t*n*p], &p, &beta[t*p], &incOne, &zero, tmp_n3, &incOne);
	 for(i = 0; i < n; i++, j++){
	   if(missing[j] == 1){
	     Y[n*t+i] = rnorm(tmp_n3[i] + u[n*t+i], sqrt(theta[t*nTheta+tauSqIndx]));
	     if(getFitted){
	       ySamples[s*Ntn+t*n+i] = Y[n*t+i];
	     }
	   }else{
	     if(getFitted){
	       ySamples[s*Ntn+t*n+i] = rnorm(tmp_n3[i] + u[n*t+i], sqrt(theta[t*nTheta+tauSqIndx]));
	     }
	   }
	 }
       }
     }
     
     for(t = 0; t < Nt; t++){
       
       ////////////////////
       //update B_t
       ////////////////////
       //Sigma_Beta_t
       F77_NAME(dcopy)(&pp, sigmaEta, &incOne, tmp_pp, &incOne);
       if(t < (Nt-1)){
	 F77_NAME(dscal)(&pp, &two, tmp_pp, &incOne);
       }
       
       tmp = 1.0/theta[t*nTheta+tauSqIndx]; 
       F77_NAME(daxpy)(&pp, &tmp, &XX[t*pp], &incOne, tmp_pp, &incOne);
   
       F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf4 failed\n");}
       F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
       
       //mu_Beta_t
       for(i = 0; i < n; i++){
	 tmp_n[i] = (Y[n*t+i] - u[n*t+i])*tmp;
       }
       F77_NAME(dgemv)(ntran, &p, &n, &one, &X[t*n*p], &p, tmp_n, &incOne, &zero, tmp_p, &incOne);
       
       if(t == 0){
	 for(i = 0; i < p; i++){
	   tmp_p2[i] = beta0[i]+beta[p*(t+1)+i];
	 }
       }if(t == (Nt-1)){
	 F77_NAME(dcopy)(&p, &beta[p*(t-1)], &incOne, tmp_p2, &incOne);
       }else{
	 for(i = 0; i < p; i++){
	   tmp_p2[i] = beta[p*(t-1)+i]+beta[p*(t+1)+i];
	 }
       }

       F77_NAME(dsymv)(lower, &p, &one, sigmaEta, &p, tmp_p2, &incOne, &one, tmp_p, &incOne); 
       F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &incOne, &zero, tmp_p2, &incOne); 

       F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf5 failed\n");}
       mvrnorm(&beta[t*p], tmp_p2, tmp_pp, p);

       ////////////////////
       //update u_t
       ////////////////////
       gamma[0] = theta[t*nTheta+sigmaSqIndx];
       gamma[1] = logitInv(theta[t*nTheta+phiIndx], phiUnif[t*2], phiUnif[t*2+1]);
       if(covModel == "matern"){
       	 gamma[2] = logitInv(theta[t*nTheta+nuIndx], nuUnif[t*2], nuUnif[t*2+1]);
       }
	 
       spCovLT(coordsD, n, gamma, covModel, C);
    
       F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf6 failed\n");}
       F77_NAME(dpotri)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotri failed\n");}

       F77_NAME(dgemv)(ytran, &p, &n, &one, &X[t*n*p], &p, &beta[t*p], &incOne, &zero, tmp_n3, &incOne);

       if(t < (Nt-1)){

       	 //t+1
       	 gamma[0] = theta[(t+1)*nTheta+sigmaSqIndx];
       	 gamma[1] = logitInv(theta[(t+1)*nTheta+phiIndx], phiUnif[(t+1)*2], phiUnif[(t+1)*2+1]);
         if(covModel == "matern"){
       	   gamma[2] = logitInv(theta[(t+1)*nTheta+nuIndx], nuUnif[(t+1)*2], nuUnif[(t+1)*2+1]);
       	 }
	 
       	 spCovLT(coordsD, n, gamma, covModel, C2);
	 
       	 F77_NAME(dpotrf)(lower, &n, C2, &n, &info); if(info != 0){error("c++ error: dpotrf7 failed\n");}
       	 F77_NAME(dpotri)(lower, &n, C2, &n, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	 
       	 //mu
       	 F77_NAME(dsymv)(lower, &n, &one, C2, &n, &u[(t+1)*n], &incOne, &zero, tmp_n, &incOne);

       	 if(t > 0){
       	   F77_NAME(dsymv)(lower, &n, &one, C, &n, &u[(t-1)*n], &incOne, &one, tmp_n, &incOne);
       	 }
	 
       	 for(i = 0; i < n; i++){
       	   tmp_n[i] += (Y[t*n+i] - tmp_n3[i])*tmp;
       	 }
	 
       	 //sigma
       	 F77_NAME(daxpy)(&nn, &one, C, &incOne, C2, &incOne);
	 
       	 for(i = 0; i < n; i++){
       	   C2[i*n+i] += tmp;
       	 }
	 
       	 F77_NAME(dpotrf)(lower, &n, C2, &n, &info); if(info != 0){error("c++ error: dpotrf8 failed\n");}
       	 F77_NAME(dpotri)(lower, &n, C2, &n, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	 	 
       	 //finish mu
       	 F77_NAME(dsymv)(lower, &n, &one, C2, &n, tmp_n, &incOne, &zero, tmp_n2, &incOne);

       	 F77_NAME(dpotrf)(lower, &n, C2, &n, &info); if(info != 0){error("c++ error: dpotrf9 failed\n");}
       	 mvrnorm(&u[t*n], tmp_n2, C2, n);

       }else{//t==Nt
	 
       	 F77_NAME(dcopy)(&nn, C, &incOne, C2, &incOne);

       	 //mu
       	 F77_NAME(dsymv)(lower, &n, &one, C2, &n, &u[(t-1)*n], &incOne, &zero, tmp_n, &incOne);

       	 for(i = 0; i < n; i++){
       	   tmp_n[i] += (Y[t*n+i] - tmp_n3[i])*tmp;
       	 }

       	 //sigma
       	 for(i = 0; i < n; i++){
       	   C2[i*n+i] += tmp;
       	 }
	 
       	 F77_NAME(dpotrf)(lower, &n, C2, &n, &info); if(info != 0){error("c++ error: dpotrf10 failed\n");}
       	 F77_NAME(dpotri)(lower, &n, C2, &n, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	 
       	 //finish mu	 
       	 F77_NAME(dsymv)(lower, &n, &one, C2, &n, tmp_n, &incOne, &zero, tmp_n2, &incOne);
	 
       	 F77_NAME(dpotrf)(lower, &n, C2, &n, &info); if(info != 0){error("c++ error: dpotrf11 failed\n");}
       	 mvrnorm(&u[t*n], tmp_n2, C2, n);
       }
 
       ////////////////////
       //update tau^2
       ////////////////////
       for(i = 0; i < n; i++){
       	 tmp_n[i] = Y[t*n+i] - tmp_n3[i] - u[t*n+i];
       }
       theta[t*nTheta+tauSqIndx] = 1.0/rgamma(tauSqIG[t*2]+n/2.0, 
       					      1.0/(tauSqIG[t*2+1]+0.5*F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n, &incOne)));
       
       ////////////////////
       //update sigma^2
       ////////////////////
       F77_NAME(dcopy)(&n, &u[n*t], &incOne, tmp_n, &incOne);
       if(t > 0){
       	 for(i = 0; i < n; i++){
       	   tmp_n[i] = tmp_n[i] - u[n*(t-1)+i];
       	 }
       }

       F77_NAME(dsymv)(lower, &n, &one, C, &n, tmp_n, &incOne, &zero, tmp_n2, &incOne);
       
       theta[t*nTheta+sigmaSqIndx] = 1.0/rgamma(sigmaSqIG[t*2]+n/2.0, 
       						1.0/(sigmaSqIG[t*2+1]+0.5*F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n2, &incOne)*theta[t*nTheta+sigmaSqIndx]));

       ////////////////////
       //update phi
       ////////////////////
       if(t > 0){
       	 for(i = 0; i < n; i++){
       	   tmp_n[i] = u[n*t+i]-u[n*(t-1)+i];
       	 }
       }else{
       	 F77_NAME(dcopy)(&n, &u[n*t], &incOne, tmp_n, &incOne);
       }

       //current
       gamma[0] = theta[t*nTheta+sigmaSqIndx];
       gamma[1] = logitInv(theta[t*nTheta+phiIndx], phiUnif[t*2], phiUnif[t*2+1]);
       if(covModel == "matern"){
       	 gamma[2] = logitInv(theta[t*nTheta+nuIndx], nuUnif[t*2], nuUnif[t*2+1]);
       }

       spCovLT(coordsD, n, gamma, covModel, C);
       
       logDet = 0;
       F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf12 failed\n");}
       for(i = 0; i < n; i++) logDet += 2*log(C[i*n+i]);
       F77_NAME(dpotri)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	 
       F77_NAME(dsymv)(lower, &n, &one, C, &n, tmp_n, &incOne, &zero, tmp_n2, &incOne);
       logPost = -0.5*logDet-0.5*F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n2, &incOne);

       logPost += log(gamma[1] - phiUnif[t*2]) + log(phiUnif[t*2+1] - gamma[1]);
       if(covModel == "matern"){
       	 logPost += log(gamma[2] - nuUnif[t*2]) + log(nuUnif[t*2+1] - gamma[2]);
       }
       
       //cand
       gamma[0] = theta[t*nTheta+sigmaSqIndx];
       gamma[1] = logitInv(rnorm(theta[t*nTheta+phiIndx], phiTuning[t]), phiUnif[t*2], phiUnif[t*2+1]);
       if(covModel == "matern"){
       	 gamma[2] = logitInv(rnorm(theta[t*nTheta+nuIndx], nuTuning[t]), nuUnif[t*2], nuUnif[t*2+1]);
       }
       
       spCovLT(coordsD, n, gamma, covModel, C);
       
       logDetCand = 0;
       F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf13 failed\n");}
       for(i = 0; i < n; i++) logDetCand += 2*log(C[i*n+i]);
       F77_NAME(dpotri)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotri failed\n");}
       
       F77_NAME(dsymv)(lower, &n, &one, C, &n, tmp_n, &incOne, &zero, tmp_n2, &incOne);
       logPostCand = -0.5*logDetCand-0.5*F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n2, &incOne);
       
       logPostCand += log(gamma[1] - phiUnif[t*2]) + log(phiUnif[t*2+1] - gamma[1]);
       if(covModel == "matern"){
       	 logPostCand += log(gamma[2] - nuUnif[t*2]) + log(nuUnif[t*2+1] - gamma[2]);
       }

       logMHRatio = logPostCand - logPost;
       
       if(runif(0.0,1.0) <= exp(logMHRatio)){
       	 theta[t*nTheta+phiIndx] = logit(gamma[1], phiUnif[t*2], phiUnif[t*2+1]);
       	 if(covModel == "matern"){
       	   theta[t*nTheta+nuIndx] = logit(gamma[2], nuUnif[t*2], nuUnif[t*2+1]);
       	 }
       	 accept++;
       	 batchAccept++;
       }

       R_CheckUserInterrupt();
     }//end Nt
     
     F77_NAME(dcopy)(&Ntp, beta, &incOne, &betaSamples[s*Ntp], &incOne);
     F77_NAME(dcopy)(&Ntn, u, &incOne, &uSamples[s*Ntn], &incOne);
     F77_NAME(dcopy)(&NtnTheta, theta, &incOne, &thetaSamples[s*NtnTheta], &incOne);
     
     //update Sigma_eta
     for(i = 0; i < p; i++){
       tmp_p[i] = beta[i]-beta0[i];
     }
     F77_NAME(dgemm)(ntran, ytran, &p, &p, &incOne, &one, tmp_p, &p, tmp_p, &p, &zero, tmp_pp, &p); 
     
     for(t = 1; t < Nt; t++){
       for(i = 0; i < p; i++){
     	 tmp_p[i] = beta[p*t+i]-beta[p*(t-1)+i];
       }
       F77_NAME(dgemm)(ntran, ytran, &p, &p, &incOne, &one, tmp_p, &p, tmp_p, &p, &one, tmp_pp, &p); 
     }

     F77_NAME(daxpy)(&pp, &one, SigmaEtaIW_S, &incOne, tmp_pp, &incOne);
     
     rwish(tmp_pp, SigmaEtaIW_df+Nt, p, sigmaEta, tmp_pp2, 1);

     // F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
     // F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
     
     // for(i = 1; i < p; i++){
     //   for(j = 0; j < i; j++){
     // 	 tmp_pp[i*p+j] = tmp_pp[j*p+i];
     //   }
     // }
     // riwishart(tmp_pp, SigmaEtaIW_df+Nt, p, C, C2, C3, sigmaEta);

     F77_NAME(dcopy)(&pp, sigmaEta, &incOne, &sigmaEtaSamples[s*pp], &incOne);
     
     //report
     if(verbose){
       if(status == nReport){
	 Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
	 Rprintf("Report interval Mean Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/(nReport*Nt));
	 Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept/(s*Nt));
	 Rprintf("-------------------------------------------------\n");
         #ifdef Win32
	 R_FlushConsole();
         #endif
	 status = 0;
	 batchAccept = 0;
       }
     }
     status++;

    }//end sample loop
    
    PutRNGstate();

    //untransform variance variables
    for(s = 0; s < nSamples; s++){
      for(t = 0; t < Nt; t++){
    	thetaSamples[s*NtnTheta+t*nTheta+phiIndx] = logitInv(thetaSamples[s*NtnTheta+t*nTheta+phiIndx], phiUnif[t*2], phiUnif[t*2+1]);
	
    	if(covModel == "matern"){
	  thetaSamples[s*NtnTheta+t*nTheta+nuIndx] = logitInv(thetaSamples[s*NtnTheta+t*nTheta+nuIndx], nuUnif[t*2], nuUnif[t*2+1]);
    	}
      }
    }
    
    // //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 5;
    if(getFitted){
      nResultListObjs++;
    }

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    //samples
    SET_VECTOR_ELT(result_r, 0, beta0Samples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.0.samples")); 

    SET_VECTOR_ELT(result_r, 1, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.beta.samples")); 

    SET_VECTOR_ELT(result_r, 2, sigmaEtaSamples_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("p.sigma.eta.samples"));

    SET_VECTOR_ELT(result_r, 3, uSamples_r);
    SET_VECTOR_ELT(resultName_r, 3, mkChar("p.u.samples"));

    SET_VECTOR_ELT(result_r, 4, thetaSamples_r);
    SET_VECTOR_ELT(resultName_r, 4, mkChar("p.theta.samples"));

    if(getFitted){
      SET_VECTOR_ELT(result_r, 5, ySamples_r);
      SET_VECTOR_ELT(resultName_r, 5, mkChar("p.y.samples"));
    }

    namesgets(result_r, resultName_r);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
