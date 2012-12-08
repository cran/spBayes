#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP spPPGLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r, SEXP family_r, SEXP weights_r,
	       SEXP m_r, SEXP knotsD_r, SEXP knotsCoordsD_r, 
	       SEXP betaPrior_r, SEXP betaNorm_r, SEXP sigmaSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	       SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r, SEXP w_strStarting_r,
	       SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP nuTuning_r, SEXP betaTuning_r, SEXP w_strTuning_r,
	       SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r){
    
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
    int pp = p*p;
    int n = INTEGER(n_r)[0];

    std::string family = CHAR(STRING_ELT(family_r,0));

    int *weights = INTEGER(weights_r);

    //covariance model
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));

    int m = INTEGER(m_r)[0];
    double *knotsD = REAL(knotsD_r);
    double *knotsCoordsD = REAL(knotsCoordsD_r);

    //priors and starting
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));

    double *betaMu = NULL;
    double *betaSd = NULL;
    
    if(betaPrior == "normal"){
      betaMu = REAL(VECTOR_ELT(betaNorm_r, 0)); 
      betaSd = REAL(VECTOR_ELT(betaNorm_r, 1));
    }
    
    double *sigmaSqIG = REAL(sigmaSqIG_r);
    double *phiUnif = REAL(phiUnif_r);

    double phiStarting = REAL(phiStarting_r)[0];
    double sigmaSqStarting = REAL(sigmaSqStarting_r)[0];
    double *betaStarting = REAL(betaStarting_r);
    double *w_strStarting = REAL(w_strStarting_r);

    double sigmaSqIGa = sigmaSqIG[0]; double sigmaSqIGb = sigmaSqIG[1];
    double phiUnifa = phiUnif[0]; double phiUnifb = phiUnif[1];

    //if matern
    double *nuUnif = NULL;
    double nuStarting = 0;
    double nuUnifa = 0, nuUnifb = 0;
    if(covModel == "matern"){
      nuUnif = REAL(nuUnif_r);
      nuStarting = REAL(nuStarting_r)[0];
      nuUnifa = nuUnif[0]; nuUnifb = nuUnif[1]; 
    }

    //tuning
    double *betaTuning = (double *) R_alloc(p*p, sizeof(double)); 
    F77_NAME(dcopy)(&pp, REAL(betaTuning_r), &incOne, betaTuning, &incOne);
    double phiTuning = sqrt(REAL(phiTuning_r)[0]);
    double sigmaSqTuning = sqrt(REAL(sigmaSqTuning_r)[0]);
    double *w_strTuning = REAL(w_strTuning_r);
   
    double nuTuning = 0;
    if(covModel == "matern")
      nuTuning = sqrt(REAL(nuTuning_r)[0]);

    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      Rprintf("Using non-modified predictive process with %i knots.\n\n", m);
    
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
  
      Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);
      Rprintf("\n");
      
      Rprintf("\tphi Unif hyperpriors a=%.5f and b=%.5f\n", phiUnifa, phiUnifb);
      Rprintf("\n");
      
      if(covModel == "matern"){
	Rprintf("\tnu Unif hyperpriors a=%.5f and b=%.5f\n", nuUnifa, nuUnifb);	  
      }

      Rprintf("Metropolis tuning values:\n");
  
      Rprintf("\tbeta tuning:\n");
      printMtrx(betaTuning, p, p);
      Rprintf("\n"); 

      Rprintf("\tsigma.sq tuning: %.5f\n", sigmaSqTuning);
      Rprintf("\n");

      Rprintf("\tphi tuning: %.5f\n", phiTuning);
      Rprintf("\n");

      if(covModel == "matern"){
	Rprintf("\tnu tuning: %.5f\n", nuTuning);
	Rprintf("\n");
      }

      Rprintf("Metropolis starting values:\n");
  
      Rprintf("\tbeta starting:\n");
      Rprintf("\t"); printVec(betaStarting, p);
      Rprintf("\n"); 

      Rprintf("\tsigma.sq starting: %.5f\n", sigmaSqStarting);
      Rprintf("\n");

      Rprintf("\tphi starting: %.5f\n", phiStarting);
      Rprintf("\n");

      if(covModel == "matern"){
	Rprintf("\tnu starting: %.5f\n", nuStarting);
	Rprintf("\n");
      }

    } 

    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/
    int nn = n*n, nm = n*m, mm = m*m;

    //spatial parameters
    int nParams, betaIndx, sigmaSqIndx, phiIndx, nuIndx;

    if(covModel != "matern"){
      nParams = p+2;//sigma^2, phi
      betaIndx = 0; sigmaSqIndx = betaIndx+p; phiIndx = sigmaSqIndx+1;
    }else{
      nParams = p+3;//sigma^2, phi, nu
      betaIndx = 0; sigmaSqIndx = betaIndx+p; phiIndx = sigmaSqIndx+1; nuIndx = phiIndx+1;
    }

    double *spParams = (double *) R_alloc(nParams, sizeof(double));
    
    //set starting
    F77_NAME(dcopy)(&p, betaStarting, &incOne, &spParams[betaIndx], &incOne);

    spParams[sigmaSqIndx] = log(sigmaSqStarting);

    spParams[phiIndx] = logit(phiStarting, phiUnifa, phiUnifb);

    if(covModel == "matern") 
      spParams[nuIndx] = logit(nuStarting, nuUnifa, nuUnifb);

    double *wCurrent = (double *) R_alloc(n, sizeof(double));
    double *w_strCurrent = (double *) R_alloc(m, sizeof(double));
    F77_NAME(dcopy)(&m, w_strStarting, &incOne, w_strCurrent, &incOne);

    //samples and random effects
    SEXP w_r, w_str_r, samples_r, accept_r;

    PROTECT(w_r = allocMatrix(REALSXP, n, nSamples)); nProtect++; 
    double *w = REAL(w_r); zeros(w, n*nSamples);

    PROTECT(w_str_r = allocMatrix(REALSXP, m, nSamples)); nProtect++; 
    double *w_str = REAL(w_str_r); zeros(w_str, m*nSamples);

    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    double *samples = REAL(samples_r);

    PROTECT(accept_r = allocMatrix(REALSXP, 1, 1)); nProtect++;


    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, rtnStatus=0, accept=0, batchAccept = 0;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0;
  
    double *P = (double *) R_alloc(nm, sizeof(double));
    double *K = (double *) R_alloc(mm, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *theta = (double *) R_alloc(3, sizeof(double)); //phi, nu, and perhaps more in the future

    double *candSpParams = (double *) R_alloc(nParams, sizeof(double));
    double *w_strCand = (double *) R_alloc(m, sizeof(double));
    double *wCand = (double *) R_alloc(n, sizeof(double));
    double sigmaSq, phi, nu;
    double *beta = (double *) R_alloc(p, sizeof(double));

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
    for(s = 0; s < nSamples; s++){
 
      //propose   
      mvrnorm(&candSpParams[betaIndx], &spParams[betaIndx], betaTuning, p, false);
      F77_NAME(dcopy)(&p, &candSpParams[betaIndx], &incOne, beta, &incOne);

      candSpParams[sigmaSqIndx] = rnorm(spParams[sigmaSqIndx], sigmaSqTuning);
      sigmaSq = theta[0] = exp(candSpParams[sigmaSqIndx]);

      candSpParams[phiIndx] = rnorm(spParams[phiIndx], phiTuning);
      phi = theta[1] = logitInv(candSpParams[phiIndx], phiUnifa, phiUnifb);

      if(covModel == "matern"){
	candSpParams[nuIndx] = rnorm(spParams[nuIndx], nuTuning);
	nu = theta[2] = logitInv(candSpParams[nuIndx], nuUnifa, nuUnifb);
      }

      for(i = 0; i < m; i++){
	w_strCand[i] = rnorm(w_strCurrent[i], sqrt(w_strTuning[i]));
      }
      
      //construct covariance matrices 
      spCovLT(knotsD, m, theta, covModel, K);
      spCov(knotsCoordsD, nm, theta, covModel, P);
    
      //invert C and log det cov
      detCand = 0;
      F77_NAME(dpotrf)(lower, &m, K, &m, &info); if(info != 0){error("c++ error: Cholesky failed in spGLM\n");}
      for(i = 0; i < m; i++) detCand += 2*log(K[i*m+i]);
      F77_NAME(dpotri)(lower, &m, K, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in spGLM\n");}
      
      //make \tild{w}
      F77_NAME(dsymv)(lower, &m, &one, K, &m, w_strCand, &incOne, &zero, tmp_m, &incOne);     
      F77_NAME(dgemv)(ytran, &m, &n, &one, P, &m, tmp_m, &incOne, &zero, wCand, &incOne);
      
      //Likelihood with Jacobian  
      logPostCand = 0.0;

      if(betaPrior == "normal"){
	for(i = 0; i < p; i++){
	  logPostCand += dnorm(beta[i], betaMu[i], betaSd[i], 1);
	}
      }

      logPostCand += -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq); 
       
      logPostCand += log(phi - phiUnifa) + log(phiUnifb - phi); 

      if(covModel == "matern"){
	logPostCand += log(nu - nuUnifa) + log(nuUnifb - nu);   
      }

      F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
      
      if(family == "binomial"){
	logPostCand += binomial_logpost(n, Y, tmp_n, wCand, weights);
      }else if(family == "poisson"){
	logPostCand += poisson_logpost(n, Y, tmp_n, wCand, weights);
      }else{
	error("c++ error: family misspecification in spGLM\n");
      }

      //(-1/2) * tmp_n` *  C^-1 * tmp_n
      logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&m, w_strCand, &incOne, tmp_m, &incOne);

      //
      //MH accept/reject	
      //      
  
      //MH ratio with adjustment
      logMHRatio = logPostCand - logPostCurrent;
      
      if(runif(0.0,1.0) <= exp(logMHRatio)){
	F77_NAME(dcopy)(&nParams, candSpParams, &incOne, spParams, &incOne);
	F77_NAME(dcopy)(&n, wCand, &incOne, wCurrent, &incOne);
	F77_NAME(dcopy)(&m, w_strCand, &incOne, w_strCurrent, &incOne);
	logPostCurrent = logPostCand;
	accept++;
	batchAccept++;
      }
      
      /******************************
          Save samples and report
      *******************************/
      F77_NAME(dcopy)(&nParams, spParams, &incOne, &samples[s*nParams], &incOne);
      F77_NAME(dcopy)(&n, wCurrent, &incOne, &w[s*n], &incOne);
      F77_NAME(dcopy)(&m, w_strCurrent, &incOne, &w_str[s*m], &incOne);
      
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
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    //untransform variance variables
    for(s = 0; s < nSamples; s++){
      samples[s*nParams+sigmaSqIndx] = exp(samples[s*nParams+sigmaSqIndx]);

      samples[s*nParams+phiIndx] = logitInv(samples[s*nParams+phiIndx], phiUnifa, phiUnifb);

      if(covModel == "matern")
	samples[s*nParams+nuIndx] = logitInv(samples[s*nParams+nuIndx], nuUnifa, nuUnifb);
    }
   
    //calculate acceptance rate
    REAL(accept_r)[0] = 100.0*accept/s;

    //make return object
    SEXP result, resultNames;
    
    int nResultListObjs = 4;

    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;

   //samples
    SET_VECTOR_ELT(result, 0, samples_r);
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.beta.theta.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));
    
    SET_VECTOR_ELT(result, 2, w_r);
    SET_VECTOR_ELT(resultNames, 2, mkChar("p.w.samples"));

    SET_VECTOR_ELT(result, 3, w_str_r);
    SET_VECTOR_ELT(resultNames, 3, mkChar("p.w.knots.samples"));
  
    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    
  }
}
