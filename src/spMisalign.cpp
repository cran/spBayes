#include <algorithm>
#include <string>
// #ifdef _OPENMP
// #include <omp.h>
// #endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {

  SEXP spMisalign(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,
		  SEXP betaPrior_r, SEXP betaNorm_r, 
		  SEXP KPrior_r, SEXP KPriorName_r, 
		  SEXP PsiPrior_r, 
		  SEXP nuUnif_r, SEXP phiUnif_r,
		  SEXP phiStarting_r, SEXP AStarting_r, SEXP PsiStarting_r, SEXP nuStarting_r, 
		  SEXP phiTuning_r, SEXP ATuning_r, SEXP PsiTuning_r, SEXP nuTuning_r, 
		  SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r){

    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, b, s, ii, jj, kk, info, nProtect= 0;
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
    int *p = INTEGER(p_r);
    int *n = INTEGER(n_r);
    int m = INTEGER(m_r)[0];
    int nLTr = m*(m-1)/2+m;

    int N = 0;
    int P = 0;
    for(i = 0; i < m; i++){
      N += n[i];
      P += p[i];
    }

    int mm = m*m;
    int NN = N*N;
    int NP = N*P;
    int PP = P*P;

    double *coordsD = REAL(coordsD_r);

    std::string covModel = CHAR(STRING_ELT(covModel_r,0));

    //priors
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));
    double *betaMu = NULL;
    double *betaC = NULL;
    
    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(P, sizeof(double));
      F77_NAME(dcopy)(&P, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);
      
      betaC = (double *) R_alloc(PP, sizeof(double)); 
      F77_NAME(dcopy)(&PP, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaC, &incOne);
    }

    double *phiUnif = REAL(phiUnif_r);

    std::string KPriorName = CHAR(STRING_ELT(KPriorName_r,0));
    double KIW_df = 0; double *KIW_S = NULL;
    double *ANormMu = NULL; double *ANormC = NULL;

    if(KPriorName == "IW"){
      KIW_S = (double *) R_alloc(mm, sizeof(double));
      KIW_df = REAL(VECTOR_ELT(KPrior_r, 0))[0]; KIW_S = REAL(VECTOR_ELT(KPrior_r, 1));
    }else{//assume A normal (can add more specifications later)
      ANormMu = (double *) R_alloc(nLTr, sizeof(double));
      ANormC = (double *) R_alloc(nLTr, sizeof(double));
      
      for(i = 0; i < nLTr; i++){
	ANormMu[i] = REAL(VECTOR_ELT(KPrior_r, 0))[i];
	ANormC[i] = REAL(VECTOR_ELT(KPrior_r, 1))[i];
      }
    }

    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    double *PsiIGa = NULL; double *PsiIGb = NULL;

    if(nugget){
      PsiIGa = (double *) R_alloc(m, sizeof(double));
      PsiIGb = (double *) R_alloc(m, sizeof(double));
      
      for(i = 0; i < m; i++){
	PsiIGa[i] = REAL(VECTOR_ELT(PsiPrior_r, 0))[i];
	PsiIGb[i] = REAL(VECTOR_ELT(PsiPrior_r, 1))[i];
      }
    }
 
    //matern
    double *nuUnif = NULL;
    if(covModel == "matern"){
      nuUnif = REAL(nuUnif_r);
    }

    bool amcmc = static_cast<bool>(INTEGER(amcmc_r)[0]);
    int nBatch = INTEGER(nBatch_r)[0];
    int batchLength = INTEGER(batchLength_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    int nSamples = nBatch*batchLength;
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
 
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i outcome variables.\n\n", m);
      Rprintf("Number of observations within each outcome:"); printVec(n, m);
      Rprintf("\nNumber of covariates for each outcome (including intercept if specified):"); printVec(p, m);
      Rprintf("\nTotal number of observations: %i\n\n", N);
      Rprintf("Total number of covariates (including intercept if specified): %i\n\n", P);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      if(amcmc){
	Rprintf("Using adaptive MCMC.\n\n");
	Rprintf("\tNumber of batches %i.\n", nBatch);
	Rprintf("\tBatch length %i.\n", batchLength);
	Rprintf("\ttarget acceptance rate %.5f.\n", acceptRate);
	Rprintf("\n");
      }else{
	Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      }
      
      if(!nugget){
	Rprintf("Psi not included in the model (i.e., no nugget model).\n\n");
      }

      Rprintf("Priors and hyperpriors:\n");
      
      if(betaPrior == "flat"){
	Rprintf("\tbeta flat.\n");
      }else{
	Rprintf("\tbeta normal:\n");
	Rprintf("\tmu:"); printVec(betaMu, P);
	Rprintf("\tcov:\n"); printMtrx(betaC, P, P);
      }
      Rprintf("\n");
      
      if(KPriorName == "IW"){
	Rprintf("\tK IW hyperpriors df=%.5f, S=\n", KIW_df);
	printMtrx(KIW_S, m, m);
      }else{
	Rprintf("\tA Normal hyperpriors\n");
	Rprintf("\t\tparameter\tmean\tvar\n");
	for(j = 0, i = 0; j < m; j++){
	  for(k = j; k < m; k++, i++){
	    Rprintf("\t\tA[%i,%i]\t\t%3.1f\t%1.2f\n", j+1, k+1, ANormMu[i], ANormC[i]);
	  }
	}
      }
      Rprintf("\n"); 
      
      if(nugget){
	Rprintf("\tDiag(Psi) IG hyperpriors\n");
	Rprintf("\t\tparameter\tshape\tscale\n");
	for(j = 0; j < m; j++){
	  Rprintf("\t\tPsi[%i,%i]\t%3.1f\t%1.2f\n", j+1, j+1, PsiIGa[j], PsiIGb[j]);
	}
      }
      Rprintf("\n");  

      Rprintf("\tphi Unif hyperpriors\n");
      Rprintf("\t\tparameter\ta\tb\n");
      for(j = 0; j < m; j++){
	Rprintf("\t\tphi[%i]\t\t%0.5f\t%0.5f\n", j+1, phiUnif[j*2], phiUnif[j*2+1]);
      }
      Rprintf("\n");   
      
      if(covModel == "matern"){
	Rprintf("\tnu Unif hyperpriors\n");
	for(j = 0; j < m; j++){
	  Rprintf("\t\tnu[%i]\t\t%0.5f\t%0.5f\n", j+1, nuUnif[j*2], nuUnif[j*2+1]);
	}
	Rprintf("\n");   
      }
      
    }
 
    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/
    //spatial parameters
    int nParams, AIndx, PsiIndx, phiIndx, nuIndx;

    if(!nugget && covModel != "matern"){
      nParams = nLTr+m;//A, phi
      AIndx = 0; phiIndx = nLTr;
    }else if(nugget && covModel != "matern"){
      nParams = nLTr+m+m;//A, diag(Psi), phi
      AIndx = 0; PsiIndx = nLTr; phiIndx = PsiIndx+m;
    }else if(!nugget && covModel == "matern"){
      nParams = nLTr+2*m;//A, phi, nu
      AIndx = 0; phiIndx = nLTr, nuIndx = phiIndx+m;
    }else{
      nParams = nLTr+3*m;//A, diag(Psi), phi, nu
      AIndx = 0; PsiIndx = nLTr, phiIndx = PsiIndx+m, nuIndx = phiIndx+m;
     }
    
    double *params = (double *) R_alloc(nParams, sizeof(double));

    //starting
    covTrans(REAL(AStarting_r), &params[AIndx], m);

    if(nugget){
      for(i = 0; i < m; i++){
	params[PsiIndx+i] = log(REAL(PsiStarting_r)[i]);
      }   
    }

    for(i = 0; i < m; i++){
      params[phiIndx+i] = logit(REAL(phiStarting_r)[i], phiUnif[i*2], phiUnif[i*2+1]);
      
      if(covModel == "matern"){
    	params[nuIndx+i] = logit(REAL(nuStarting_r)[i], nuUnif[i*2], nuUnif[i*2+1]);
      }
    }

    //tuning and fixed
    double *tuning = (double *) R_alloc(nParams, sizeof(double));
    int *fixed = (int *) R_alloc(nParams, sizeof(int)); zeros(fixed, nParams);

    for(i = 0; i < nLTr; i++){
      tuning[AIndx+i] = REAL(ATuning_r)[i];
      if(tuning[AIndx+i] == 0){
    	fixed[AIndx+i] = 1;
      }
    }
    
    if(nugget){
      for(i = 0; i < m; i++){
	tuning[PsiIndx+i] = REAL(PsiTuning_r)[i];
	if(tuning[PsiIndx+i] == 0){
	  fixed[PsiIndx+i] = 1;
	}
      }	
    }

    for(i = 0; i < m; i++){
      tuning[phiIndx+i] = REAL(phiTuning_r)[i];
      if(tuning[phiIndx+i] == 0){
    	fixed[phiIndx+i] = 1;
      }
      
      if(covModel == "matern"){
    	tuning[nuIndx+i] = REAL(nuTuning_r)[i];
    	if(tuning[nuIndx+i] == 0){
    	  fixed[nuIndx+i] = 1;
    	}
      }
    }

    for(i = 0; i < nParams; i++){
      tuning[i] = log(sqrt(tuning[i]));
    }

    //return stuff  
    SEXP samples_r, accept_r, tuning_r;
    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++;

    if(amcmc){
      PROTECT(accept_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++; 
      PROTECT(tuning_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++;  
    }else{
      PROTECT(accept_r = allocMatrix(REALSXP, 1, nSamples/nReport)); nProtect++; 
    }

    // /*****************************************
    //    Set-up MCMC alg. vars. matrices etc.
    // *****************************************/
    int status=1, batchAccept=0, reportCnt=0;
    double logMHRatio =0, logPostCurrent = R_NegInf, logPostCand = 0, det = 0, paramsjCurrent = 0;
    double Q, logDetK, SKtrace;
    
    double *paramsCurrent = (double *) R_alloc(nParams, sizeof(double));
    double *accept = (double *) R_alloc(nParams, sizeof(double)); zeros(accept, nParams);
    
    double *C = (double *) R_alloc(NN, sizeof(double)); 
    double *K = (double *) R_alloc(mm, sizeof(double));
    double *Psi = (double *) R_alloc(m, sizeof(double));
    double *A = (double *) R_alloc(mm, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));

    int P1 = P+1;
    double *vU = (double *) R_alloc(N*P1, sizeof(double));
    double *z = (double *) R_alloc(N, sizeof(double));
    double *tmp_N = (double *) R_alloc(N, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_PP = (double *) R_alloc(PP, sizeof(double));
    double *tmp_P = (double *) R_alloc(P, sizeof(double));
    double *tmp_NN = NULL;
    double *Cbeta = NULL;

    if(betaPrior == "normal"){
      tmp_NN = (double *) R_alloc(NN, sizeof(double));
      Cbeta = (double *) R_alloc(NN, sizeof(double));
      
      F77_NAME(dgemv)(ntran, &N, &P, &negOne, X, &N, betaMu, &incOne, &zero, z, &incOne);
      F77_NAME(daxpy)(&N, &one, Y, &incOne, z, &incOne);

      F77_NAME(dsymm)(rside, lower, &N, &P, &one, betaC, &P, X, &N, &zero, vU, &N);
      F77_NAME(dgemm)(ntran, ytran, &N, &N, &P, &one, vU, &N, X, &N, &zero, tmp_NN, &N);
    }
     
    int sl, sk;

    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    GetRNGstate();
    
    for(b = 0, s = 0; b < nBatch; b++){
      for(i = 0; i < batchLength; i++, s++){
    	for(j = 0; j < nParams; j++){
	  
    	  //propose
    	  if(amcmc){
    	    if(fixed[j] == 1){
    	      paramsjCurrent = params[j];
    	    }else{
    	      paramsjCurrent = params[j];
    	      params[j] = rnorm(paramsjCurrent, exp(tuning[j]));
    	    }
    	  }else{
    	    F77_NAME(dcopy)(&nParams, params, &incOne, paramsCurrent, &incOne);
	    
    	    for(j = 0; j < nParams; j++){
    	      if(fixed[j] == 1){
    		params[j] = params[j];
    	      }else{
    		params[j] = rnorm(params[j], exp(tuning[j]));
    	      }
    	    }
    	  }
	  
    	  //extract and transform
    	  covTransInvExpand(&params[AIndx], A, m);
	  
    	  for(k = 0; k < m; k++){
    	    phi[k] = logitInv(params[phiIndx+k], phiUnif[k*2], phiUnif[k*2+1]);
	    
    	    if(covModel == "matern"){
    	      nu[k] = logitInv(params[nuIndx+k], nuUnif[k*2], nuUnif[k*2+1]);
    	    }	  
    	  }
	  
    	  if(nugget){
	    for(k = 0; k < m; k++){
	      Psi[k] = exp(params[PsiIndx+k]);
	    }
	  }
	  
	  //construct covariance matrix
	  sl = sk = 0;
	  
	  for(k = 0; k < m; k++){
	    sl = 0;
	    for(l = 0; l < m; l++){
	      for(kk = 0; kk < n[k]; kk++){
		for(jj = 0; jj < n[l]; jj++){
		  C[(sl+jj)*N+(sk+kk)] = 0.0;
		  for(ii = 0; ii < m; ii++){
		    C[(sl+jj)*N+(sk+kk)] += A[k+m*ii]*A[l+m*ii]*spCor(coordsD[(sl+jj)*N+(sk+kk)], phi[ii], nu[ii], covModel);
		  }
		}
	      }
	      sl += n[l];
	    }
	    sk += n[k];
	  }
	  
    	  if(nugget){
    	    sl = 0;
	    for(l = 0; l < m; l++){
	      for(k = 0; k < n[l]; k++){
	    	C[(sl+k)*N+(sl+k)] += Psi[l];
	      }
	      sl += n[l];
	    }
    	  }

    	  if(betaPrior == "normal"){    
    	    for(k = 0; k < N; k++){
    	      for(l = k; l < N; l++){
    	    	Cbeta[k*N+l] = C[k*N+l]+tmp_NN[k*N+l];
    	      }
    	    }
	    
    	    det = 0;
    	    F77_NAME(dpotrf)(lower, &N, Cbeta, &N, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
    	    for(k = 0; k < N; k++) det += 2*log(Cbeta[k*N+k]);
	    
    	    F77_NAME(dcopy)(&N, z, &incOne, tmp_N, &incOne);
    	    F77_NAME(dtrsv)(lower, ntran, nUnit, &N, Cbeta, &N, tmp_N, &incOne);//u = L^{-1}(y-X'beta)
	    
    	    Q = pow(F77_NAME(dnrm2)(&N, tmp_N, &incOne),2);
    	  }else{//beta flat
    	    det = 0;
    	    F77_NAME(dpotrf)(lower, &N, C, &N, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
    	    for(k = 0; k < N; k++) det += 2*log(C[k*N+k]);
	    
    	    F77_NAME(dcopy)(&N, Y, &incOne, vU, &incOne);
    	    F77_NAME(dcopy)(&NP, X, &incOne, &vU[N], &incOne);

    	    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &N, &P1, &one, C, &N, vU, &N);//L^{-1}[v:U] = [y:X]
	    
    	    F77_NAME(dgemm)(ytran, ntran, &P, &P, &N, &one, &vU[N], &N, &vU[N], &N, &zero, tmp_PP, &P); //U'U
    	    F77_NAME(dpotrf)(lower, &P, tmp_PP, &P, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
    	    for(k = 0; k < P; k++) det += 2*log(tmp_PP[k*P+k]);
	    
    	    F77_NAME(dgemv)(ytran, &N, &P, &one, &vU[N], &N, vU, &incOne, &zero, tmp_P, &incOne); //U'v
    	    F77_NAME(dtrsv)(lower, ntran, nUnit, &P, tmp_PP, &P, tmp_P, &incOne);

    	    Q = pow(F77_NAME(dnrm2)(&N, vU, &incOne),2) - pow(F77_NAME(dnrm2)(&P, tmp_P, &incOne),2) ;
    	  }
	  
    	  //
    	  //priors, jacobian adjustments, and likelihood
    	  //
    	  logPostCand = 0.0;
	  
    	  if(KPriorName == "IW"){
    	    logDetK = 0.0;
    	    SKtrace = 0.0;
	    
    	    for(k = 0; k < m; k++){logDetK += 2*log(A[k*m+k]);}
	    
    	    //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
    	    for(k = 0; k < m; k++){logPostCand += (m-k)*log(A[k*m+k])+log(A[k*m+k]);}
	    
    	    //S*K^-1
    	    F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
    	    F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
    	    for(k = 0; k < m; k++){SKtrace += tmp_mm[k*m+k];}
    	    logPostCand += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;
    	  }else{	     
    	    for(k = 0; k < nLTr; k++){
    	      logPostCand += dnorm(params[AIndx+k], ANormMu[k], sqrt(ANormC[k]), 1);
    	    }
    	  }
	  
    	  if(nugget){
	    for(k = 0; k < m; k++){
	      logPostCand += -1.0*(1.0+PsiIGa[k])*log(Psi[k])-PsiIGb[k]/Psi[k]+log(Psi[k]);
	    }
	  }
	  
    	  for(k = 0; k < m; k++){
    	    logPostCand += log(phi[k] - phiUnif[k*2]) + log(phiUnif[k*2+1] - phi[k]); 
	    
    	    if(covModel == "matern"){
    	      logPostCand += log(nu[k] - nuUnif[k*2]) + log(nuUnif[k*2+1] - nu[k]);  
    	    }
    	  }
	  
    	  logPostCand += -0.5*det-0.5*Q;
	  
    	  //
    	  //MH accept/reject	
    	  //      
    	  logMHRatio = logPostCand - logPostCurrent;
	  
    	  if(runif(0.0,1.0) <= exp(logMHRatio)){
    	    logPostCurrent = logPostCand;
	    
    	    if(amcmc){
    	      accept[j]++;
    	    }else{
    	      accept[0]++;
    	      batchAccept++;
    	    }
	    
    	  }else{
    	    if(amcmc){
    	      params[j] = paramsjCurrent;
    	    }else{
    	      F77_NAME(dcopy)(&nParams, paramsCurrent, &incOne, params, &incOne);
    	    }
    	  }
	  
    	  if(!amcmc){
    	    break;
    	  }
	}//end params
	
    	/******************************
               Save samples
    	*******************************/
    	F77_NAME(dcopy)(&nParams, params, &incOne, &REAL(samples_r)[s*nParams], &incOne);
	
    	R_CheckUserInterrupt();
      }//end batch
      
      //adjust tuning
      if(amcmc){
    	for(j = 0; j < nParams; j++){
    	  REAL(accept_r)[b*nParams+j] = accept[j]/batchLength;
    	  REAL(tuning_r)[b*nParams+j] = tuning[j];
	  
    	  if(accept[j]/batchLength > acceptRate){
    	    tuning[j] += std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
    	  }else{
    	    tuning[j] -= std::min(0.01, 1.0/sqrt(static_cast<double>(b)));
    	  }
    	  accept[j] = 0.0;
    	}
      }
      
      //report
      if(status == nReport){
	
    	if(verbose){
    	  if(amcmc){
    	    Rprintf("Batch: %i of %i, %3.2f%%\n", b+1, nBatch, 100.0*(b+1)/nBatch);
    	    Rprintf("\tparameter\tacceptance\ttuning\n");
    	    for(j = 0, i = 0; j < m; j++){
    	      for(k = j; k < m; k++, i++){
    		Rprintf("\tA[%i,%i]\t\t%3.1f\t\t%1.2f\n", j+1, k+1, 100.0*REAL(accept_r)[b*nParams+AIndx+i], exp(tuning[AIndx+i]));
    	      }
    	    }
    	    if(nugget){
	      for(j = 0; j < m; j++){
		Rprintf("\tPsi[%i,%i]\t%3.1f\t\t%1.2f\n", j+1, j+1, 100.0*REAL(accept_r)[b*nParams+PsiIndx+j], exp(tuning[PsiIndx+j]));
	      }
	    }
    	    for(j = 0; j < m; j++){
    	      Rprintf("\tphi[%i]\t\t%3.1f\t\t%1.2f\n", j+1, 100.0*REAL(accept_r)[b*nParams+phiIndx+j], exp(tuning[phiIndx+j]));
    	    }
    	    if(covModel == "matern"){
    	      Rprintf("\n");
    	      for(j = 0; j < m; j++){
    		Rprintf("\tnu[%i]\t\t%3.1f\t\t%1.2f\n", j+1, 100.0*REAL(accept_r)[b*nParams+nuIndx+j], exp(tuning[nuIndx+j]));
    	      } 
    	    }
    	  }else{
    	    Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
    	    Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
    	    Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept[0]/s);
    	  }
    	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
    	  R_FlushConsole();
          #endif
    	}

    	if(!amcmc){
    	  REAL(accept_r)[reportCnt] = 100.0*batchAccept/nReport;
    	  reportCnt++;
    	}
	
    	status = 0;
    	batchAccept = 0;
      }
      status++;
      
    }//end sample loop
    
    PutRNGstate();
    
    //untransform variance variables
    for(s = 0; s < nSamples; s++){
      
      covTransInv(&REAL(samples_r)[s*nParams+AIndx], &REAL(samples_r)[s*nParams+AIndx], m);
      
      if(nugget){
	for(i = 0; i < m; i++){
	  REAL(samples_r)[s*nParams+PsiIndx+i] = exp(REAL(samples_r)[s*nParams+PsiIndx+i]);
	}
      }
      
      for(i = 0; i < m; i++){
    	REAL(samples_r)[s*nParams+phiIndx+i] = logitInv(REAL(samples_r)[s*nParams+phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
    	if(covModel == "matern"){
    	  REAL(samples_r)[s*nParams+nuIndx+i] = logitInv(REAL(samples_r)[s*nParams+nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
    	}
      }
    }
    
    //make return object
    SEXP result_r, resultName_r;  
    int nResultListObjs = 2;

    if(amcmc){
      nResultListObjs++;
    }
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, samples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.theta.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, accept_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("acceptance"));
    
    if(amcmc){
      SET_VECTOR_ELT(result_r, 2, tuning_r);
      SET_VECTOR_ELT(resultName_r, 2, mkChar("tuning"));
    }
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
   
    return(result_r);
  }
}


  
