#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <string>
#include "util.h"

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

extern "C" {

  SEXP spSVC(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP Z_r, SEXP m_r, SEXP n_r, SEXP coords_r, SEXP q_r, 
	     SEXP KDiag_r, SEXP betaPrior_r, SEXP betaNorm_r, SEXP KPrior_r, SEXP tauSqIG_r, SEXP nuUnif_r, SEXP phiUnif_r,
	     SEXP phiStarting_r, SEXP AStarting_r, SEXP tauSqStarting_r, SEXP nuStarting_r,
	     SEXP phiTuning_r, SEXP ATuning_r, SEXP tauSqTuning_r, SEXP nuTuning_r, 
	     SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r,
	     SEXP nThreads_r){

    /*****************************************
                Common variables
    *****************************************/
    int i, j, k, l, b, s, ii, jj, h, info, nProtect=0;
    char const *lower = "L";
    char const *nUnit = "N";
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
    double *Z = REAL(Z_r);
    int p = INTEGER(p_r)[0];
    int pp = p*p;
    int m = INTEGER(m_r)[0];
    int nLTr = m*(m-1)/2+m;
    int mm = m*m;
    int n = INTEGER(n_r)[0];
    int nn = n*n;
    int np = n*p;
    int nm = n*m;
    int nmnm = nm*nm;
    int q = INTEGER(q_r)[0];
    
    double *coords = REAL(coords_r);

    std::string covModel = CHAR(STRING_ELT(covModel_r,0));

    //priors
    bool KDiag = static_cast<bool>(INTEGER(KDiag_r)[0]);
    
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));
    double *betaMu = NULL;
    double *betaC = NULL;
    
    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(p, sizeof(double));
      F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);

      betaC = (double *) R_alloc(pp, sizeof(double)); 
      F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaC, &incOne);
    }

    double *Ka = NULL;
    double *Kb = NULL;
    double *phiUnifa = (double *) R_alloc(m, sizeof(double));
    double *phiUnifb = (double *) R_alloc(m, sizeof(double));
    double *nuUnifa = NULL;
    double *nuUnifb = NULL;
    double maxNuUnifb = 0; //needed for thread safe bessel
    if(covModel == "matern"){
      nuUnifa = (double *) R_alloc(m, sizeof(double));
      nuUnifb = (double *) R_alloc(m, sizeof(double));
    }

    if(KDiag){

      Ka = (double *) R_alloc(m, sizeof(double));
      Kb = (double *) R_alloc(m, sizeof(double));
      for(i = 0; i < m; i++){
	Ka[i] = REAL(VECTOR_ELT(KPrior_r, 0))[i];
	Kb[i] = REAL(VECTOR_ELT(KPrior_r, 1))[i];
      }
      
    }else{

      Ka = (double *) R_alloc(1, sizeof(double));
      Ka[0] = REAL(VECTOR_ELT(KPrior_r, 0))[0];
      Kb = (double *) R_alloc(mm, sizeof(double));
      F77_NAME(dcopy)(&mm, REAL(VECTOR_ELT(KPrior_r, 1)), &incOne, Kb, &incOne);
    }

    
    for(i = 0; i < m; i++){
      phiUnifa[i] = REAL(phiUnif_r)[i*2];
      phiUnifb[i] = REAL(phiUnif_r)[i*2+1];
      
      if(covModel == "matern"){
	nuUnifa[i] = REAL(nuUnif_r)[i*2];
	nuUnifb[i] = REAL(nuUnif_r)[i*2+1];
	if(nuUnifb[i] > maxNuUnifb){
	  maxNuUnifb = nuUnifb[i];
	}
      }
    
    }

    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);
    double tauSqIGa = 0, tauSqIGb = 0;
    if(nugget){
      tauSqIGa = REAL(tauSqIG_r)[0]; tauSqIGb = REAL(tauSqIG_r)[1]; 
    }

    bool amcmc = static_cast<bool>(INTEGER(amcmc_r)[0]);
    int nBatch = INTEGER(nBatch_r)[0];
    int batchLength = INTEGER(batchLength_r)[0];
    double acceptRate = REAL(acceptRate_r)[0];
    int nSamples = nBatch*batchLength;
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];

#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads = %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif    
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i.\n\n", p);
      Rprintf("Number of space varying covariates %i.\n\n", m);      
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      if(amcmc){
	Rprintf("Using adaptive MCMC.\n\n");
	Rprintf("\tNumber of batches %i.\n", nBatch);
	Rprintf("\tBatch length %i.\n", batchLength);
	Rprintf("\tTarget acceptance rate %.5f.\n", acceptRate);
	Rprintf("\n");
      }else{
	Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      }

      if(!nugget){
	Rprintf("tau.sq not included in the model (i.e., no nugget model).\n\n");
      }

      Rprintf("Priors and hyperpriors:\n");
      
      if(betaPrior == "flat"){
	Rprintf("\tbeta flat.\n");
      }else{
	Rprintf("\tbeta normal:\n");
	Rprintf("\tmu:"); printVec(betaMu, p);
	Rprintf("\tcov:\n"); printMtrx(betaC, p, p);
	Rprintf("\n");
      }
 
      if(KDiag){
	Rprintf("\tDiag(K) IG hyperpriors\n");
	Rprintf("\t\tparameter\tshape\t\tscale\n");
	for(j = 0; j < m; j++){
	  Rprintf("\t\tK[%i,%i]\t\t%5f\t%5f\n", j+1, j+1, Ka[j], Kb[j]);
	}
	Rprintf("\n"); 
      }else{
	Rprintf("\tK IW hyperpriors:\n\tdf: %.5f\n\tS:\n", Ka[0]);
	printMtrx(Kb, m, m);
	Rprintf("\n"); 
      }
      
      Rprintf("\tphi Unif lower bound hyperpriors:"); printVec(phiUnifa, m);
      Rprintf("\tphi Unif upper bound hyperpriors:"); printVec(phiUnifb, m);
      
      if(covModel == "matern"){
	Rprintf("\n"); 
	Rprintf("\tnu Unif lower bound hyperpriors:"); printVec(nuUnifa, m);
	Rprintf("\tnu Unif upper bound hyperpriors:"); printVec(nuUnifb, m);
      }
      
      if(nugget){
	Rprintf("\n"); 
	Rprintf("\ttau.sq IG hyperpriors shape=%.5f and scale=%.5f\n", tauSqIGa, tauSqIGb); 
      }
      
#ifdef _OPENMP
   Rprintf("\nSource compiled with OpenMP, posterior sampling is using %i thread(s).\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 

    /*****************************************
         Set-up MCMC sample matrices etc.
    *****************************************/ 
    //parameters
    int nParams, AIndx, tauSqIndx = 0, phiIndx, nuIndx = 0;

    if(KDiag){
      j = m;
    }else{
      j = nLTr;
    }
    
    if(!nugget && covModel != "matern"){
      nParams = j+m;//A, phi
      AIndx = 0; phiIndx = AIndx+j;
    }else if(nugget && covModel != "matern"){
      nParams = j + m + 1;//A, tau^2, phi
      AIndx = 0; tauSqIndx = AIndx+j; phiIndx = tauSqIndx+1;  
    }else if(!nugget && covModel == "matern"){
      nParams = j + 2*m;//A, phi, nu
      AIndx = 0; phiIndx = AIndx+j; nuIndx = phiIndx+m;
    }else{
      nParams = j + 2*m + 1;//A, tau^2,  phi, nu
      AIndx = 0; tauSqIndx = AIndx+j; phiIndx = tauSqIndx+1; nuIndx = phiIndx+m;
    }
    
    double *params = (double *) R_alloc(nParams, sizeof(double));
    
    //starting
    if(KDiag){
      for(i = 0; i < m; i++){
	params[AIndx+i] = log(REAL(AStarting_r)[i]);
      }
    }else{
      covTrans(REAL(AStarting_r), &params[AIndx], m);
    }
    
    for(i = 0; i < m; i++){  
      params[phiIndx+i] = logit(REAL(phiStarting_r)[i], phiUnifa[i], phiUnifb[i]);
      
      if(covModel == "matern"){
	params[nuIndx+i] = logit(REAL(nuStarting_r)[i], nuUnifa[i], nuUnifb[i]);
      }
      
    }
    
    if(nugget){
      params[tauSqIndx] = log(REAL(tauSqStarting_r)[0]);
    }
    
    //tuning and fixed
    double *tuning = (double *) R_alloc(nParams, sizeof(double));
    int *fixed = (int *) R_alloc(nParams, sizeof(int)); zeros(fixed, nParams);
    
    for(i = 0; i < j; i++){//j was defined above as the number of A parameters
      tuning[AIndx+i] = REAL(ATuning_r)[i];
      if(tuning[AIndx+i] == 0){
	fixed[AIndx+i] = 1;
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
    
    if(nugget){
      tuning[tauSqIndx] = REAL(tauSqTuning_r)[0];
      if(tuning[tauSqIndx] == 0){
	fixed[tauSqIndx] = 1;
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
      PROTECT(accept_r = allocMatrix(REALSXP, 1, floor(static_cast<double>(nSamples/nReport)))); nProtect++; 
    }

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status=1, batchAccept=0, reportCnt=0;
    double logMHRatio =0, logPostCurrent = R_NegInf, logPostCand = 0, det = 0, paramsjCurrent = 0;
   
    double *paramsCurrent = (double *) R_alloc(nParams, sizeof(double));
    double *accept = (double *) R_alloc(nParams, sizeof(double)); zeros(accept, nParams);

    double Q, logDetK, SKtrace;
    double *A = (double *) R_alloc(mm, sizeof(double)); zeros(A, mm); 
    double *K = (double *) R_alloc(nmnm, sizeof(double)); zeros(K, nmnm);
    double *C = (double *) R_alloc(nn, sizeof(double)); zeros(C, nn);
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double)); zeros(nu, m); //this just remains empty of not matern
    double tauSq = 0;
    
    double *D = (double *) R_alloc(nn, sizeof(double));
    distN(coords, n, coords, n, q, D);
    
    int p1 = p+1;
    double *vU = (double *) R_alloc(n*p1, sizeof(double));

    double *z = (double *) R_alloc(n, sizeof(double)); 
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    
    int threadID = 0;
    double *tmp_m = (double *) R_alloc(nThreads*m, sizeof(double));
    int bessel_ws_inc = static_cast<int>(1.0+maxNuUnifb);
    double *bessel_ws = (double *) R_alloc(nThreads*bessel_ws_inc, sizeof(double));

    double *tmp_nn = NULL;
    double *Cbeta = NULL;
    
    if(betaPrior == "normal"){
      tmp_nn = (double *) R_alloc(nn, sizeof(double));
      Cbeta = (double *) R_alloc(nn, sizeof(double));

      F77_NAME(dgemv)(ntran, &n, &p, &negOne, X, &n, betaMu, &incOne, &zero, z, &incOne);
      F77_NAME(daxpy)(&n, &one, Y, &incOne, z, &incOne);

      F77_NAME(dsymm)(rside, lower, &n, &p, &one, betaC, &p, X, &n, &zero, vU, &n);
      F77_NAME(dgemm)(ntran, ytran, &n, &n, &p, &one, vU, &n, X, &n, &zero, tmp_nn, &n);
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
	  if(KDiag == false){
	    covTransInvExpand(&params[AIndx], A, m);
	  }
	  
	  for(k = 0; k < m; k++){

	    if(KDiag){
	      A[k*m+k] = sqrt(exp(params[AIndx+k]));
	    }
	    
	    phi[k] = logitInv(params[phiIndx+k], phiUnifa[k], phiUnifb[k]);
	    
	    if(covModel == "matern"){
	      nu[k] = logitInv(params[nuIndx+k], nuUnifa[k], nuUnifb[k]);
	    }
	  
	  }

	  if(nugget){
	    tauSq = exp(params[tauSqIndx]);
	  }

	  //construct covariance matrix

#ifdef _OPENMP
#pragma omp parallel for private(ii, k, l, h, threadID)
#endif
	  for(jj = 0; jj < n; jj++){
#ifdef _OPENMP
	    threadID = omp_get_thread_num();
#endif
	    for(ii = 0; ii < n; ii++){	
	      for(k = 0; k < m; k++){
		for(l = 0; l < m; l++){
		  K[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
		  for(h = 0; h < m; h++){
		    K[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCorTS(D[jj*n+ii], phi[h], nu[h], covModel, &bessel_ws[threadID*bessel_ws_inc]);
		  }
		}
	      }
	      
	      //mk\tild{X} K \tild{X}^T
	      for(k = 0; k < m; k++){
		tmp_m[k+threadID*m] = F77_NAME(ddot)(&m, &Z[jj], &n, &K[(k+jj*m)*nm+(ii*m)], &incOne); 
	      }
	      
	      C[jj*n+ii] = F77_NAME(ddot)(&m, &tmp_m[threadID*m], &incOne, &Z[ii], &n);	      
	    }
	  }

	  if(nugget){
	    for(k = 0; k < n; k++){
	      C[k*n+k] += tauSq;
	    }
	  }

	  if(betaPrior == "normal"){

	    for(k = 0; k < n; k++){
	      for(l = k; l < n; l++){
	    	Cbeta[k*n+l] = C[k*n+l]+tmp_nn[k*n+l];
	      }
	    }

	    det = 0;
	    F77_NAME(dpotrf)(lower, &n, Cbeta, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	    for(k = 0; k < n; k++) det += log(Cbeta[k*n+k]);
	    
	    F77_NAME(dcopy)(&n, z, &incOne, tmp_n, &incOne);
	    F77_NAME(dtrsv)(lower, ntran, nUnit, &n, Cbeta, &n, tmp_n, &incOne);//L[u] = (y-X'beta)

	    Q = pow(F77_NAME(dnrm2)(&n, tmp_n, &incOne),2);
	  }else{//beta flat
	    det = 0;

	    F77_NAME(dpotrf)(lower, &n, C, &n, &info); if(info != 0){error("c++ error: dpotrf failed\n");}//L
	    for(k = 0; k < n; k++) det += log(C[k*n+k]);
	    
	    F77_NAME(dcopy)(&n, Y, &incOne, vU, &incOne);
	    F77_NAME(dcopy)(&np, X, &incOne, &vU[n], &incOne);
	    F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &n, &p1, &one, C, &n, vU, &n);//L[v:U] = [y:X]

	    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, &vU[n], &n, &vU[n], &n, &zero, tmp_pp, &p); //U'U
	    F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	    for(k = 0; k < p; k++) det += log(tmp_pp[k*p+k]);
	    
	    F77_NAME(dgemv)(ytran, &n, &p, &one, &vU[n], &n, vU, &incOne, &zero, tmp_p, &incOne); //U'v
	    F77_NAME(dtrsv)(lower, ntran, nUnit, &p, tmp_pp, &p, tmp_p, &incOne);

	    Q = pow(F77_NAME(dnrm2)(&n, vU, &incOne),2) - pow(F77_NAME(dnrm2)(&p, tmp_p, &incOne),2) ;
	  }
	  
	  //
	  //priors, jacobian adjustments, and likelihood
	  //
	  logPostCand = 0;

	  if(KDiag == false){
	    
	    logDetK = 0.0;
	    SKtrace = 0.0;
	    
	    for(k = 0; k < m; k++){logDetK += 2*log(A[k*m+k]);}
	    
	    //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	    for(k = 0; k < m; k++){logPostCand += (m-k)*log(A[k*m+k])+log(A[k*m+k]);}
	    
	    //S*K^-1
	    F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	    F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, Kb, &m, &zero, tmp_mm, &m);
	    for(k = 0; k < m; k++){SKtrace += tmp_mm[k*m+k];}
	    logPostCand += -0.5*(Ka[0]+m+1)*logDetK - 0.5*SKtrace;
	    
	  }
	  
	  for(k = 0; k < m; k++){

	    if(KDiag){
	      logPostCand += -1.0*(1.0+Ka[k])*log(pow(A[k*m+k],2))-Kb[k]/pow(A[k*m+k],2)+log(pow(A[k*m+k],2));
	    }
	    
	    logPostCand += log(phi[k] - phiUnifa[k]) + log(phiUnifb[k] - phi[k]); 
	    
	    if(covModel == "matern"){
	      logPostCand += log(nu[k] - nuUnifa[k]) + log(nuUnifb[k] - nu[k]);   
	    }
	  }
	  
	  if(nugget){
	    logPostCand += -1.0*(1.0+tauSqIGa)*log(tauSq)-tauSqIGb/tauSq+log(tauSq);
	  }
	  
	  
	  logPostCand += -det-0.5*Q;
	  
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
      	  REAL(tuning_r)[b*nParams+j] = exp(tuning[j]);
	  
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

	    for(k = 0; k < m; k++){
	   
	      if(KDiag){
		Rprintf("\tsigma.sq[%i]\t%3.1f\t\t%1.5f\n", k+1, 100.0*REAL(accept_r)[b*nParams+AIndx+k], exp(tuning[AIndx+k]));
	      }
	      
	      Rprintf("\tphi[%i]\t\t%3.1f\t\t%1.5f\n", k+1, 100.0*REAL(accept_r)[b*nParams+phiIndx+k], exp(tuning[phiIndx+k]));
	      if(covModel == "matern"){
		Rprintf("\tnu[%i]\t\t%3.1f\t\t%1.5f\n", k+1, 100.0*REAL(accept_r)[b*nParams+nuIndx+k], exp(tuning[nuIndx+k]));
	      }

	    }

	    for(j = 0, i = 0; j < m; j++){
	      for(k = j; k < m; k++, i++){
		Rprintf("\tA[%i,%i]\t\t%3.1f\t\t%1.5f\n", j+1, k+1, 100.0*REAL(accept_r)[b*nParams+AIndx+i], exp(tuning[AIndx+i]));
	      }
	    } 
	   
      	    if(nugget){
	      Rprintf("\ttau.sq\t\t%3.1f\t\t%1.5f\n", 100.0*REAL(accept_r)[b*nParams+tauSqIndx], exp(tuning[tauSqIndx]));
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

	batchAccept = 0;
	status = 0;

      }
      status++;
      
    }//end sample loop
    
    PutRNGstate();
 
    //untransform variance variables
    for(s = 0; s < nSamples; s++){

      for(k = 0; k < m; k++){

	if(KDiag){
	  REAL(samples_r)[s*nParams+AIndx+k] = exp(REAL(samples_r)[s*nParams+AIndx+k]);
	}
	
	REAL(samples_r)[s*nParams+phiIndx+k] = logitInv(REAL(samples_r)[s*nParams+phiIndx+k], phiUnifa[k], phiUnifb[k]);
	
	if(covModel == "matern"){
	  REAL(samples_r)[s*nParams+nuIndx+k] = logitInv(REAL(samples_r)[s*nParams+nuIndx+k], nuUnifa[k], nuUnifb[k]);
	}
      }

      if(KDiag == false){
	covTransInv(&REAL(samples_r)[s*nParams+AIndx], &REAL(samples_r)[s*nParams+AIndx], m);
      }
      
      if(nugget){
	REAL(samples_r)[s*nParams+tauSqIndx] = exp(REAL(samples_r)[s*nParams+tauSqIndx]);
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
