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

void proposeUnif(double *paramVec, double &tuning, int &indx, double &candParam, double *candParamVec,
		double &a, double &b){
  
  bool pass = true;
  int alarm = 0;
  
  while(pass){
    candParamVec[indx] = rnorm(paramVec[indx], tuning);
    candParam = candParamVec[indx];
    if(candParam > a && candParam < b) pass = false;
    
    if(alarm == 1000){
      Rprintf("warning: likely a starting value well outside the a uniform prior's support.\n");
      alarm = 0;
    }
    alarm++;
    
    R_CheckUserInterrupt();
  }
}


void proposeIG(double *paramVec, double &tuning, int &indx, double &candParam, double *candParamVec){
  candParamVec[indx] = rnorm(paramVec[indx], tuning);
  candParam = exp(candParamVec[indx]);
}

extern "C" {

  SEXP splm(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coordsD_r,
	    SEXP isPp_r, SEXP isModPp_r, SEXP m_r, SEXP knotsD_r, SEXP coordsKnotsD_r, SEXP nugget_r, 
	    SEXP betaFixed_r, SEXP sigmaSqFixed_r, SEXP tauSqFixed_r, SEXP phiFixed_r, SEXP nuFixed_r,
	    SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	    SEXP phiStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r, SEXP betaStarting_r,
	    SEXP phiTuning_r, SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP nuTuning_r, 
	    SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r, SEXP spEffects_r){


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

    double *coordsD = REAL(coordsD_r);

    //covariance model
    string covModel = CHAR(STRING_ELT(covModel_r,0));

    //if predictive process
    bool isPp = static_cast<bool>(INTEGER(isPp_r)[0]);
    bool isModPp = static_cast<bool>(INTEGER(isModPp_r)[0]);
    int m = 0;
    double *knotsD = NULL;
    double *coordsKnotsD = NULL;
    if(isPp){
      m = INTEGER(m_r)[0];
      knotsD = REAL(knotsD_r);
      coordsKnotsD = REAL(coordsKnotsD_r);
    }

    //priors and starting
    double *sigmaSqIG = REAL(sigmaSqIG_r);
    double *phiUnif = REAL(phiUnif_r);
    double phiStarting = REAL(phiStarting_r)[0];
    double sigmaSqStarting = REAL(sigmaSqStarting_r)[0];
    double *betaStarting = REAL(betaStarting_r);

    double sigmaSqIGa = sigmaSqIG[0]; double sigmaSqIGb = sigmaSqIG[1];
    double phiUnifa = phiUnif[0]; double phiUnifb = phiUnif[1];

    //if nugget
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    double *tauSqIG = NULL;
    double tauSqStarting = 0;
    double tauSqIGa = 0, tauSqIGb = 0;
    if(nugget){
      tauSqIG = REAL(tauSqIG_r);
      tauSqStarting = REAL(tauSqStarting_r)[0];
      tauSqIGa = tauSqIG[0]; tauSqIGb = tauSqIG[1]; 
    }

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
    double phiTuning = sqrt(REAL(phiTuning_r)[0]);
    double sigmaSqTuning = sqrt(REAL(sigmaSqTuning_r)[0]);
    double tauSqTuning = 0;
    double nuTuning = 0;
    if(nugget)
      tauSqTuning = sqrt(REAL(tauSqTuning_r)[0]);

    if(covModel == "matern")
      nuTuning = sqrt(REAL(nuTuning_r)[0]);


    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    bool spEffects = static_cast<bool>(INTEGER(spEffects_r)[0]);

    //if fixed
    bool betaFixed = static_cast<bool>(INTEGER(betaFixed_r)[0]);
    bool sigmaSqFixed = static_cast<bool>(INTEGER(sigmaSqFixed_r)[0]);
    bool tauSqFixed = static_cast<bool>(INTEGER(tauSqFixed_r)[0]);
    bool phiFixed = static_cast<bool>(INTEGER(phiFixed_r)[0]);
    bool nuFixed = static_cast<bool>(INTEGER(nuFixed_r)[0]);
    

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      if(isPp){
	if(isModPp)
	  Rprintf("Using modified predictive process with %i knots.\n\n", m);
	else
	  Rprintf("Using non-modified predictive process with %i knots.\n\n", m);
      }
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);

      Rprintf("Priors and hyperpriors:\n");
      if(!betaFixed)
	Rprintf("\tbeta flat.\n");   
      
      if(!sigmaSqFixed)
	Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);

      if(nugget){
	if(!tauSqFixed)
	  Rprintf("\ttau.sq IG hyperpriors shape=%.5f and scale=%.5f\n", tauSqIGa, tauSqIGb); 
      }

      if(!phiFixed)
	Rprintf("\tphi Unif hyperpriors a=%.5f and b=%.5f\n", phiUnifa, phiUnifb);
      
      if(covModel == "matern"){
	if(!nuFixed)
	  Rprintf("\tnu Unif hyperpriors a=%.5f and b=%.5f\n", nuUnifa, nuUnifb);	  
      }



      if(!betaFixed || !sigmaSqFixed || !tauSqFixed || !phiFixed){
	Rprintf("\n");	
	Rprintf("Starting values:\n");
      }

      if(!betaFixed){
	Rprintf("\t");
	for(i = 0; i < p; i++)
	  Rprintf("beta_%i=%.5f ", i, betaStarting[i]);
	Rprintf("\n");
      }
      if(!sigmaSqFixed){
	Rprintf("\tsigma.sq=%.5f\n", sigmaSqStarting);
      }
      if(!tauSqFixed){
	if(nugget)
	  Rprintf("\ttau.sq=%.5f\n", tauSqStarting);
      }
      if(!phiFixed){
	Rprintf("\tphi=%.5f\n", phiStarting);
      }
      if(!nuFixed){
	if(covModel == "matern")
	  Rprintf("\tnu=%.5f\n", nuStarting);
      }



      if(!sigmaSqFixed || !tauSqFixed || !phiFixed){
	Rprintf("\n");  
	Rprintf("MH tuning values:\n");
      }
      
      if(!sigmaSqFixed)
	Rprintf("\tsigma.sq=%.5f\n", pow(sigmaSqTuning,2));
      
      if(!tauSqFixed){
	if(nugget)
	Rprintf("\ttau.sq=%.5f\n", pow(tauSqTuning,2));
      }

      if(!phiFixed){
      Rprintf("\tphi=%.5f\n", pow(phiTuning,2));
      }

      if(!nuFixed){
	if(covModel == "matern")
	  Rprintf("\tnu=%.5f\n", pow(nuTuning,2));
      }
      


      if(sigmaSqFixed || tauSqFixed || phiFixed){
	Rprintf("\n");
	Rprintf("Fixed values:\n");
      }

      if(betaFixed){
	Rprintf("\t");
	for(i = 0; i < p; i++)
	  Rprintf("beta_%i=%.5f ", i, betaStarting[i]);
	Rprintf("\n");
      }
      if(sigmaSqFixed){
	Rprintf("\tsigma.sq=%.5f\n", sigmaSqStarting);
      }
      if(tauSqFixed){
	if(nugget)
	  Rprintf("\ttau.sq=%.5f\n", tauSqStarting);
      }
      if(phiFixed){
	Rprintf("\tphi=%.5f\n", phiStarting);
      }
      if(nuFixed){
	if(covModel == "matern")
	  Rprintf("\tnu=%.5f\n", nuStarting);
      }

    } 

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
         Set-up MCMC sample matrices etc.
    *****************************************/
    int nn = n*n, nm = n*m, mm = m*m;

    //spatial parameters
    int nSpParams, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;

    if(!nugget && covModel != "matern"){
      nSpParams = 2;//sigma^2, phi
      sigmaSqIndx = 0; phiIndx = 1;
    }else if(nugget && covModel != "matern"){
      nSpParams = 3;//sigma^2, tau^2, phi
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
    }else if(!nugget && covModel == "matern"){
      nSpParams = 3;//sigma^2, phi, nu
      sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
    }else{
      nSpParams = 4;//sigma^2, tau^2, phi, nu
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;//sigma^2, tau^2, phi, nu
    }
    
    double *spParams = (double *) R_alloc(nSpParams, sizeof(double));
    
    //set starting
    spParams[sigmaSqIndx] = log(sigmaSqStarting);
    if(nugget) spParams[tauSqIndx] = log(tauSqStarting);
    spParams[phiIndx] = phiStarting;
    if(covModel == "matern") spParams[nuIndx] = nuStarting;

    //Beta parameter and set starting
    double *beta = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, betaStarting, &incOne, beta, &incOne);

    //samples, random effects, predictions
    int nParams = p+nSpParams;
    double *w = NULL;
    double *w_str = NULL;

    double *samples = NULL;
    SEXP w_r, w_str_r, samples_r, accept_r;

    PROTECT(w_r = allocMatrix(REALSXP, n, nSamples)); nProtect++; 
    w = REAL(w_r); zeros(w, n*nSamples);

    if(isPp){
      PROTECT(w_str_r = allocMatrix(REALSXP, m, nSamples)); nProtect++; 
      w_str = REAL(w_str_r); zeros(w_str, m*nSamples);
    }

    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    samples = REAL(samples_r);

    PROTECT(accept_r = allocMatrix(REALSXP, 1, 1)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, rtnStatus=0, accept=-1;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0, detCandC_str = 0, detCandE = 0;
    bool first = true, accepted = true;

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
      if(isModPp){//not needed for the non-modified pp
	tmp_nn = (double *) R_alloc(nn, sizeof(double));
      }
      E = (double *) R_alloc(n, sizeof(double));
      Einv = (double *) R_alloc(n, sizeof(double));
    }else{
      wMu = (double *) R_alloc(n, sizeof(double));
    }
 
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_n1 = (double *) R_alloc(n, sizeof(double));
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *tmp_nm1 = (double *) R_alloc(nm, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));

    double *candSpParams = (double *) R_alloc(nSpParams, sizeof(double));
    double sigmaSq, tauSq, phi, nu, sigmaSqTmp, tauSqTmp, tauSqInv, negTauSqInv;
    double logMHRatio, HastAdj;
    
    double *S_beta = (double *) R_alloc(p*p, sizeof(double));
    double *Mu_beta = (double *) R_alloc(p, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double)); 
    double *tmp_np = (double *) R_alloc(n*p, sizeof(double)); 

    //set fixed values
    if(sigmaSqFixed){
      sigmaSq = exp(spParams[sigmaSqIndx]);
      candSpParams[sigmaSqIndx] = spParams[sigmaSqIndx];
    }
    
    if(nugget){
      if(tauSqFixed){
	tauSq = exp(spParams[tauSqIndx]);
	candSpParams[tauSqIndx] = spParams[tauSqIndx];
      }
    }    

    if(phiFixed){
      phi = spParams[phiIndx];
      candSpParams[phiIndx] = spParams[phiIndx];
    }

    if(covModel == "matern"){
      if(nuFixed){
	nu = spParams[nuIndx];
	candSpParams[nuIndx] = spParams[nuIndx];
      }
    }

    //check if all MH parameters are fixed
    bool allFixed = false;
    if(sigmaSqFixed && tauSqFixed && phiFixed){
      if(covModel == "matern"){
	if(nuFixed)
	  allFixed = true;
      }else{
	allFixed = true;
      }
    }


    GetRNGstate();
    for(s = 0; s < nSamples; s++){
      
      //propose
      if(!sigmaSqFixed)
	proposeIG(spParams, sigmaSqTuning, sigmaSqIndx, sigmaSq, candSpParams);
      
      if(nugget){
	if(!tauSqFixed)
	  proposeIG(spParams, tauSqTuning, tauSqIndx, tauSq, candSpParams);
      }      

      if(!phiFixed)
	proposeUnif(spParams, phiTuning, phiIndx, phi, candSpParams, phiUnifa, phiUnifb);

     if(covModel == "matern"){
	if(!nuFixed)
	  proposeUnif(spParams, nuTuning, nuIndx, nu, candSpParams, nuUnifa, nuUnifb);
      }

    
      /******************************
          Non predictive process
      *******************************/
      if(!isPp){
	
	//make the correlation matrix
	for(i = 0; i < nn; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phi, C[i], coordsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phi, nu, C[i], coordsD[i]);
	}

	F77_NAME(dscal)(&nn, &sigmaSq, C, &incOne);
		
	if(nugget){
	  for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+tauSq;
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
	    (covModelObj->*cov1ParamPtr)(phi, C_str[i], knotsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phi, nu, C_str[i], knotsD[i]);
	}
	
	for(i = 0; i < nm; i++){
	  if(onePramPtr)
	    (covModelObj->*cov1ParamPtr)(phi, ct[i], coordsKnotsD[i]);
	  else //i.e., 2 parameter matern
	    (covModelObj->*cov2ParamPtr)(phi, nu, ct[i], coordsKnotsD[i]);
	}

	//scale by sigma^2
	F77_NAME(dscal)(&mm, &sigmaSq, C_str, &incOne);	
	F77_NAME(dscal)(&nm, &sigmaSq, ct, &incOne);

	/******************************
           Sherman-Woodbury-Morrison
	*******************************/
	if(!nugget) tauSq = 1e-10;//ridge the matrix if no nugget model
	
 	//Unmodified predictive process
	if(!isModPp){  
	  tauSqInv = 1.0/tauSq;
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
	  detCand += n*log(tauSq);
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
	  for(i = 0; i < m; i++) detCand -= 2.0*log(C_str[i*m+i]);
	  
	}else{//Modified predictive process
	  
	  if(!nugget) 
	    tauSq = 1e-10;//ridge the matrix if no nugget model
	  
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
	    E[i] = tauSq+sigmaSq-tmp_nn[i*n+i];
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
      }//finish PP
      
      /******************************
           Likelihood and MH
      *******************************/
      //Likelihood with Jacobian
      logPostCand = 0.0;
      
      if(!sigmaSqFixed)
	logPostCand = -1.0*(1.0+sigmaSqIGa)*log(sigmaSq)-sigmaSqIGb/sigmaSq+log(sigmaSq);
      
      if(nugget){
	if(!tauSqFixed)
	  logPostCand += -1.0*(1.0+tauSqIGa)*log(tauSq)-tauSqIGb/tauSq+log(tauSq);
      }
      
      //Hastings adjustment for truncated normal
      HastAdj = 0.0;
      //for phi
      if(!phiFixed){
	HastAdj = log(dTNorm(spParams[phiIndx], phi, phiTuning, phiUnifa, phiUnifb))-
	  log(dTNorm(phi, spParams[phiIndx], phiTuning, phiUnifa, phiUnifb));
      }
      
      //for nu
      if(covModel == "matern"){
	if(!nuFixed){
	  HastAdj += log(dTNorm(spParams[nuIndx], nu, nuTuning, nuUnifa, nuUnifb))-
	    log(dTNorm(nu, spParams[nuIndx], nuTuning, nuUnifa, nuUnifb));
	}
      }
      
      //Y-XB
      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
      
      //(-1/2) * tmp_n` *  C^-1 * tmp_n
      F77_NAME(dsymv)(upper, &n, &one,  C, &n, tmp_n, &incOne, &zero, tmp_n1, &incOne);
      logPostCand += -0.5*detCand-0.5*F77_NAME(ddot)(&n, tmp_n, &incOne, tmp_n1, &incOne);
      
      //
      //MH accept/reject	
      //      
      accepted = false;
      
      if(first){
	first = false;
	accepted = true; //force the first iteration to accept
      }

      //force accept if all MH parameters are fixed
      if(allFixed) accepted = true;

      //MH ratio with adjustment
      logMHRatio = logPostCand - logPostCurrent + HastAdj;
      
      if(runif(0.0,1.0) <= exp(logMHRatio)){
	accepted = true;
      }

      if(accepted){
	F77_NAME(dcopy)(&nSpParams, candSpParams, &incOne, spParams, &incOne);
	logPostCurrent = logPostCand;
	accept++;
      }
        
      /******************************
                Update Beta
      *******************************/
      if(!betaFixed){
	
	if(accepted){
	  //finish the Gibbs
	  //C^-1 X = tmp_np
	  F77_NAME(dsymm)(lside, upper, &n, &p, &one, C, &n, X, &n, &zero, tmp_np, &n);
	  
	  //(t(X) tmp_np)^{-1} = S_beta
	  F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, tmp_np, &n, &zero, S_beta, &p);
	  
	  F77_NAME(dpotrf)(upper, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
	  F77_NAME(dpotri)(upper, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky inverse failed\n" << endl;}
	  
	  //S_beta tmp_np Y = Mu_beta
	  F77_NAME(dsymv)(upper, &n, &one, C, &n, Y, &incOne, &zero, tmp_n, &incOne);
	  F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, tmp_n, &incOne, &zero, tmp_p, &incOne);
	  F77_NAME(dsymv)(upper, &p, &one, S_beta, &p, tmp_p, &incOne, &zero, Mu_beta, &incOne);
	  
	  //Gibbs draw
	  //take upper for the chol for the mv draw
	  F77_NAME(dpotrf)(upper, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
	}
       
	mvrnorm(beta, Mu_beta, S_beta, p, true);
      }
      
      /******************************
         Recover w and w* if needed
      *******************************/
      if(spEffects){
	
	//Get previously accepted sample
	sigmaSq = exp(spParams[sigmaSqIndx]);
	
	if(nugget)
	  tauSq = exp(spParams[tauSqIndx]);
	
        phi = spParams[phiIndx];
	
	if(covModel == "matern")
	  nu = spParams[nuIndx];
	

	if(!isPp){
	  if(nugget){
	    
	    //make the correlation matrix
	    for(i = 0; i < nn; i++){
	      if(onePramPtr)
		(covModelObj->*cov1ParamPtr)(phi, C[i], coordsD[i]);
	      else //i.e., 2 parameter matern
		(covModelObj->*cov2ParamPtr)(phi, nu, C[i], coordsD[i]);
	    }
	    
	    //invert C
	    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	    
	    //scale correlation matrix with 1/sigmasq and add 1/nugget to diag
	    sigmaSqTmp = 1.0/sigmaSq;
	    F77_NAME(dscal)(&nn, &sigmaSqTmp, C, &incOne);
	    
	    for(i = 0; i < n; i++) C[i*n+i] = C[i*n+i]+1.0/tauSq;
	    
	    //invert C
	    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    F77_NAME(dpotri)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	    
	    //make w mu
	    F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
	    F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	    
	    tauSqTmp = 1.0/tauSq;
	    F77_NAME(dscal)(&n, &tauSqTmp, tmp_n, &incOne);
	    
	    F77_NAME(dsymv)(upper, &n, &one, C, &n, tmp_n, &incOne, &zero, wMu, &incOne);
	    
	    //chol for the mvnorm and draw
	    F77_NAME(dpotrf)(upper, &n, C, &n, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	    
	    mvrnorm(&w[s*n], wMu, C, n, true);
	    
	  }else{//no nugget so w is just resids
	    F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, &w[s*n], &incOne);
	    F77_NAME(daxpy)(&n, &one, Y, &incOne, &w[s*n], &incOne);
	  }
	}else{//using pp
	  
	  if(!nugget) tauSq = 1e-10;//ridge the matrix if no nugget model
	  
	  //make the correlation matrix
	  for(i = 0; i < mm; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi, C_str[i], knotsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi, nu, C_str[i], knotsD[i]);
	  }
	  
	  for(i = 0; i < nm; i++){
	    if(onePramPtr)
	      (covModelObj->*cov1ParamPtr)(phi, ct[i], coordsKnotsD[i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(phi, nu, ct[i], coordsKnotsD[i]);
	  }
	  
	  //scale by sigma^2
	  F77_NAME(dscal)(&mm, &sigmaSq, C_str, &incOne);	
	  F77_NAME(dscal)(&nm, &sigmaSq, ct, &incOne);
	  
	  //invert C_str
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  //make w* Sigma
	  //ct C^{*-1}
	  F77_NAME(dsymm)(rside, upper, &n, &m, &one, C_str, &m, ct, &n, &zero, tmp_nm, &n);
	  
	  
	  //w* ~ MVN(mu_w, Sigma_w)
	  //unmodified
	  //Sigma_w = [C^{*-1} + C^{*-1} C (1/tau^2 I_n) C' C^{*-1}]^{-1}
	  //mu_w = Sigma_w [C^{*-1} C (1/tau^2 I_n) (Y-XB)]
	  
	  //modified			      
	  //Sigma_w = [C^{*-1} + C^{*-1} C (1/E) C' C^{*-1}]^{-1}
	  //mu_w = Sigma_w [C^{*-1} C (1/E) (Y-XB)]
	  //where E = I \otimes (Psi + A'A) - Diag[C'(s_i) C^{*-1} C(s_i)]^n_{i=1}
	
	  //then w = C' C^{*-1} w*

	  //Unmodified predictive process
	  if(!isModPp){
	    for(i = 0; i < n; i++) Einv[i] = 1.0/(tauSq);
	  }else{//Modified predictive process
	    //ct C^{*-1} c
	    F77_NAME(dgemm)(ntran, ytran, &n, &n, &m, &one, tmp_nm, &n, ct, &n, &zero, tmp_nn, &n);

	    for(i = 0; i < n; i++) Einv[i] = 1.0/(tauSq+sigmaSq-tmp_nn[i*n+i]);
	  }
	  
	  diagmm(n, m, Einv, tmp_nm, tmp_nm1); //(1/E ct C^{*-1})

	  //(C^{*-1} c) (1/E ct C^{*-1})
	  F77_NAME(dgemm)(ytran, ntran, &m, &m, &n, &one, tmp_nm, &n, tmp_nm1, &n, &zero, tmp_mm, &m);

	  for(i = 0; i < mm; i++) C_str[i] = C_str[i] + tmp_mm[i];
	  
	  //invert C_str
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  F77_NAME(dpotri)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed in sp.lm\n");}
	  
	  //make w* mu
	  F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, beta, &incOne, &zero, tmp_n, &incOne);
	  F77_NAME(daxpy)(&n, &one, Y, &incOne, tmp_n, &incOne);
	  	  
	  //(1/E ct C^{*-1})'(Y-XB)
	  F77_NAME(dgemv)(ytran, &n, &m, &one, tmp_nm1, &n, tmp_n, &incOne, &zero, tmp_m, &incOne);
	  
	  F77_NAME(dsymv)(upper, &m, &one, C_str, &m, tmp_m, &incOne, &zero, w_strMu, &incOne);
	  
	  //chol for the mvnorm and draw
	  F77_NAME(dpotrf)(upper, &m, C_str, &m, &info); if(info != 0){error("c++ error: Cholesky failed in sp.lm\n");}
	  mvrnorm(&w_str[s*m], w_strMu, C_str, m, true);
	  
	  //make \tild{w}
	  F77_NAME(dgemv)(ntran, &n, &m, &one, tmp_nm, &n, &w_str[s*m], &incOne, &zero, &w[s*n], &incOne);
	}//fi pp
      }//end spEffects
      
          
      /******************************
          Save samples and report
      *******************************/
   
      F77_NAME(dcopy)(&p, beta, &incOne, &samples[s*nParams], &incOne);
      F77_NAME(dcopy)(&nSpParams, spParams, &incOne, &samples[s*nParams+p], &incOne);

     //report
     if(verbose){
       if(status == nReport){
	 Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
	 Rprintf("MH Acceptance rate: %3.2f%%\n", 100.0*accept/s);
	 Rprintf("-------------------------------------------------\n");
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
   
    //untransform variance variables
    for(s = 0; s < nSamples; s++){
      samples[s*nParams+p+sigmaSqIndx] = exp(samples[s*nParams+p+sigmaSqIndx]);
      if(nugget)
	samples[s*nParams+p+tauSqIndx] = exp(samples[s*nParams+p+tauSqIndx]);
    }
    
    //final status report
    if(verbose){
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
    }
    Rprintf("-------------------------------------------------\n");
    #ifdef Win32
    R_FlushConsole();
    #endif

    //calculate acceptance rate
    REAL(accept_r)[0] = 100.0*accept/s;


    //make return object
    SEXP result, resultNames;
    
    int nResultListObjs = 0;
    
    if(spEffects){
      if(!isPp){
	nResultListObjs = 3;
      }else{
	nResultListObjs = 4;
      }
    }else{
      nResultListObjs = 2;
    }

    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, nResultListObjs)); nProtect++;

   //samples
    SET_VECTOR_ELT(result, 0, samples_r);
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));

    if(nResultListObjs == 3){
      SET_VECTOR_ELT(result, 2, w_r);
      SET_VECTOR_ELT(resultNames, 2, mkChar("sp.effects"));
    }

    if(nResultListObjs == 4){
      SET_VECTOR_ELT(result, 2, w_r);
      SET_VECTOR_ELT(resultNames, 2, mkChar("sp.effects"));

      SET_VECTOR_ELT(result, 3, w_str_r);
      SET_VECTOR_ELT(resultNames, 3, mkChar("sp.effects.knots"));
    }

    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    
  }
}
