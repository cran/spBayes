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

  SEXP spPPMvLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP g_r, SEXP coordsD_r, SEXP knotsD_r, SEXP knotsCoordsD_r, SEXP modPP_r, 
	      SEXP betaPrior_r, SEXP betaNorm_r, 
	      SEXP KPrior_r, SEXP KPriorName_r, 
	      SEXP PsiPrior_r, SEXP PsiPriorName_r, SEXP PsiDiag_r, 
	      SEXP nuUnif_r, SEXP phiUnif_r,
	      SEXP betaStarting_r, SEXP phiStarting_r, SEXP AStarting_r, SEXP LStarting_r, SEXP nuStarting_r, 
	      SEXP phiTuning_r, SEXP ATuning_r, SEXP LTuning_r, SEXP nuTuning_r, 
	      SEXP nugget_r, SEXP covModel_r, SEXP amcmc_r, SEXP nBatch_r, SEXP batchLength_r, SEXP acceptRate_r, SEXP verbose_r, SEXP nReport_r){

    /*****************************************
                Common variables
    *****************************************/
    int h, i, j, k, l, b, s, ii, jj, info, nProtect= 0;
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
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int g = INTEGER(g_r)[0];
    int nLTr = m*(m-1)/2+m;
    int nn = n*n;
    int mm = m*m;
    int nm = n*m;
    int nmnm = nm*nm;
    int nmp = nm*p;
    int pp = p*p;
    int nmm = nm*m;

    int gm = g*m;
    int gmgm = gm*gm;
    int ng = n*g;
    int nmgm = nm*gm;
    int gmp = gm*p;
    
    double *coordsD = REAL(coordsD_r);
    double *knotsD = REAL(knotsD_r);
    double *knotsCoordsD = REAL(knotsCoordsD_r);

    bool modPP = static_cast<bool>(INTEGER(modPP_r)[0]);
    std::string covModel = CHAR(STRING_ELT(covModel_r,0));

    //priors
    std::string betaPrior = CHAR(STRING_ELT(betaPrior_r,0));
    double *betaMu = NULL;
    double *betaC = NULL;
    
    if(betaPrior == "normal"){
      betaMu = (double *) R_alloc(p, sizeof(double));
      F77_NAME(dcopy)(&p, REAL(VECTOR_ELT(betaNorm_r, 0)), &incOne, betaMu, &incOne);
      
      betaC = (double *) R_alloc(pp, sizeof(double)); 
      F77_NAME(dcopy)(&pp, REAL(VECTOR_ELT(betaNorm_r, 1)), &incOne, betaC, &incOne);
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
    std::string PsiPriorName;
    bool PsiDiag = static_cast<bool>(INTEGER(PsiDiag_r)[0]);
    double PsiIW_df = 0; double *PsiIW_S = NULL;
    double *LNormMu = NULL; double *LNormC = NULL;
    double *PsiIGa = NULL; double *PsiIGb = NULL;

    if(nugget){
      PsiPriorName = CHAR(STRING_ELT(PsiPriorName_r,0));

      if(PsiDiag){
	PsiIGa = (double *) R_alloc(m, sizeof(double));
	PsiIGb = (double *) R_alloc(m, sizeof(double));
	
	for(i = 0; i < m; i++){
	  PsiIGa[i] = REAL(VECTOR_ELT(PsiPrior_r, 0))[i];
	  PsiIGb[i] = REAL(VECTOR_ELT(PsiPrior_r, 1))[i];
	}
      }else{
	if(PsiPriorName == "IW"){
	  PsiIW_S = (double *) R_alloc(mm, sizeof(double));
	  PsiIW_df = REAL(VECTOR_ELT(PsiPrior_r, 0))[0]; PsiIW_S = REAL(VECTOR_ELT(PsiPrior_r, 1));
	}else{//assume A normal (can add more specifications later)
	  LNormMu = (double *) R_alloc(nLTr, sizeof(double));
	  LNormC = (double *) R_alloc(nLTr, sizeof(double));
	  
	  for(i = 0; i < nLTr; i++){
	    LNormMu[i] = REAL(VECTOR_ELT(PsiPrior_r, 0))[i];
	    LNormC[i] = REAL(VECTOR_ELT(PsiPrior_r, 1))[i];
	  }
	}
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

    //other stuff
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
 
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tGeneral model description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());
      
      if(modPP){
	Rprintf("Using modified predictive process with %i knots.\n\n", g);
      }else{
	Rprintf("Using non-modified predictive process with %i knots.\n\n", g);
      }

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
	Rprintf("tau.sq not included in the model (i.e., no nugget model).\n\n");
      }

      Rprintf("Priors and hyperpriors:\n");
      
      if(betaPrior == "flat"){
	Rprintf("\tbeta flat.\n");
      }else{
	Rprintf("\tbeta normal:\n");
	Rprintf("\tmu:"); printVec(betaMu, p);
	Rprintf("\tcov:\n"); printMtrx(betaC, p, p);
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
	if(PsiPriorName == "IW"){
	  Rprintf("\tPsi IW hyperpriors df=%.5f, S=\n", PsiIW_df);
	  printMtrx(PsiIW_S, m, m);
	  Rprintf("\n"); 
	}else{
	  if(PsiDiag){
	    Rprintf("\tDiag(Psi) IG hyperpriors\n");
	    Rprintf("\t\tparameter\tshape\tscale\n");
	    for(j = 0; j < m; j++){
	      Rprintf("\t\tPsi[%i,%i]\t%3.1f\t%1.2f\n", j+1, j+1, PsiIGa[j], PsiIGb[j]);
	    }
	  }else{
	    Rprintf("\tL Normal hyperpriors\n");
	    Rprintf("\t\tparameter\tmean\tvar\n");
	    for(j = 0, i = 0; j < m; j++){
	      for(k = j; k < m; k++, i++){
		Rprintf("\t\tL[%i,%i]\t\t%3.1f\t%1.2f\n", j+1, k+1, LNormMu[i], LNormC[i]);
	      }
	    }
	  }
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
    int nParams, AIndx, LIndx, phiIndx, nuIndx;

    if(!nugget && covModel != "matern"){
      nParams = nLTr+m;//A, phi
      AIndx = 0; phiIndx = nLTr;
    }else if(nugget && covModel != "matern"){
      if(PsiDiag){
	nParams = nLTr+m+m;//A, diag(Psi), phi
	AIndx = 0; LIndx = nLTr; phiIndx = LIndx+m;
      }else{
	nParams = 2*nLTr+m;//A, L, phi
	AIndx = 0; LIndx = nLTr; phiIndx = LIndx+nLTr;
      }
    }else if(!nugget && covModel == "matern"){
      nParams = nLTr+2*m;//A, phi, nu
      AIndx = 0; phiIndx = nLTr, nuIndx = phiIndx+m;
    }else{
      if(PsiDiag){
	nParams = nLTr+3*m;//A, diag(Psi), phi, nu
	AIndx = 0; LIndx = nLTr, phiIndx = LIndx+m, nuIndx = phiIndx+m;
      }else{
	nParams = 2*nLTr+2*m;//A, Psi, phi, nu
	AIndx = 0; LIndx = nLTr, phiIndx = LIndx+nLTr, nuIndx = phiIndx+m;
      }
    }
    
    double *params = (double *) R_alloc(nParams, sizeof(double));

    //starting
    double *beta = (double *) R_alloc(p, sizeof(double)); 
    F77_NAME(dcopy)(&p, REAL(betaStarting_r), &incOne, beta, &incOne);

    covTrans(REAL(AStarting_r), &params[AIndx], m);

    if(nugget){
      if(PsiDiag){
	for(i = 0; i < m; i++){
	  params[LIndx+i] = log(REAL(LStarting_r)[i]);
	}
      }else{
	covTrans(REAL(LStarting_r), &params[LIndx], m);
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
      if(PsiDiag){
	for(i = 0; i < m; i++){
	  tuning[LIndx+i] = REAL(LTuning_r)[i];
	  if(tuning[LIndx+i] == 0){
	    fixed[LIndx+i] = 1;
	  }
	}	
      }else{
	for(i = 0; i < nLTr; i++){
	  tuning[LIndx+i] = REAL(LTuning_r)[i];
	  if(tuning[LIndx+i] == 0){
	    fixed[LIndx+i] = 1;
	  }
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
    SEXP samples_r, accept_r, tuning_r, betaSamples_r;
    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++;
    PROTECT(accept_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++; 
    PROTECT(tuning_r = allocMatrix(REALSXP, nParams, nBatch)); nProtect++;  
    PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++; 

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int status=0, batchAccept=0;
    double logMHRatio =0, logPostCurrent = R_NegInf, logPostCand = 0, paramsjCurrent = 0;
    double det = 0, detCurrent = 0, SKtrace, logDetK, Q;
    double priors = 0, priorsCurrent = 0;
    bool updateBeta = true;

    double *paramsCurrent = (double *) R_alloc(nParams, sizeof(double));
    double *accept = (double *) R_alloc(nParams, sizeof(double)); zeros(accept, nParams);

    double *C = (double *) R_alloc(nmnm, sizeof(double)); 
    double *P = (double *) R_alloc(nmgm, sizeof(double)); 
    double *K = (double *) R_alloc(gmgm, sizeof(double)); 
    double *D = (double *) R_alloc(nmm, sizeof(double)); 
    double *H = (double *) R_alloc(nmgm, sizeof(double));
    

    double *V = (double *) R_alloc(mm, sizeof(double));//AA'
    double *Psi = (double *) R_alloc(mm, sizeof(double));  zeros(Psi, mm); //must be cleared for diag Psi
    double *A = (double *) R_alloc(mm, sizeof(double));
    double *L = (double *) R_alloc(mm, sizeof(double));
    double *phi = (double *) R_alloc(m, sizeof(double));
    double *nu = (double *) R_alloc(m, sizeof(double));
   
    double *u = (double *) R_alloc(nm, sizeof(double)); 

    //X is transposed on the R side to simplify recovery of beta
    F77_NAME(dgemv)(ytran, &p, &nm, &negOne, X, &p, beta, &incOne, &zero, u, &incOne);
    F77_NAME(daxpy)(&nm, &one, Y, &incOne, u, &incOne);

    double *DCurrent = (double *) R_alloc(nmm, sizeof(double)); 
    double *HCurrent = (double *) R_alloc(nmgm, sizeof(double));
       
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double)); 
    double *tmp_nmp = (double *) R_alloc(nmp, sizeof(double));
    double *tmp_gmgm = (double *) R_alloc(gmgm, sizeof(double));
    double *tmp_gmp = (double *) R_alloc(gmp, sizeof(double));
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_pp2 = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_gm = (double *) R_alloc(gm, sizeof(double));
    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
   
    double *betaCInv = NULL;
    double *betaCInvMu = NULL;
    
    if(betaPrior == "normal"){
      betaCInv = (double *) R_alloc(pp, sizeof(double));
      betaCInvMu = (double *) R_alloc(p, sizeof(double));
      
      F77_NAME(dcopy)(&pp, betaC, &incOne, betaCInv, &incOne);
      F77_NAME(dpotrf)(lower, &p, betaCInv, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &p, betaCInv, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
      
      F77_NAME(dsymv)(lower, &p, &one, betaCInv, &p, betaMu, &incOne, &zero, betaCInvMu, &incOne);      
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
	  covTransInvExpand(&params[AIndx], A, m);
	  F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, A, &m, A, &m, &zero, V, &m);
	  
	  for(k = 0; k < m; k++){
	    phi[k] = logitInv(params[phiIndx+k], phiUnif[k*2], phiUnif[k*2+1]);
	    
	    if(covModel == "matern"){
	      nu[k] = logitInv(params[nuIndx+k], nuUnif[k*2], nuUnif[k*2+1]);
	    }
	  }
	  
	  if(nugget){
	    if(PsiDiag){
	      for(k = 0; k < m; k++){
		Psi[k*m+k] = exp(params[LIndx+k]);//note, I use only the first column of Psi in other routines
	      }
	    }else{
	      covTransInvExpand(&params[LIndx], L, m);
	      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, L, &m, L, &m, &zero, Psi, &m);
	    }
	  }
	  
	  //construct covariance matrix
          // #pragma omp parallel 
          // {
          // #pragma omp for private(ii, k, l, h)
	  for(jj = 0; jj < g; jj++){
	    for(ii = jj; ii < g; ii++){	
	      for(k = 0; k < m; k++){
		for(l = 0; l < m; l++){
		  K[(k+jj*m)*gm+(ii*m+l)] = 0.0; 
		  for(h = 0; h < m; h++){
		    K[(k+jj*m)*gm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(knotsD[jj*g+ii], phi[h], nu[h], covModel);
		  }
		}
	      }
	    }
	  }
	  // } //parallel for
	  
        // #pragma omp parallel 
        // {
        // #pragma omp for private(ii, k, l, h) 
	for(jj = 0; jj < n; jj++){
	  for(ii = 0; ii < g; ii++){	
	    for(k = 0; k < m; k++){
	      for(l = 0; l < m; l++){
		P[(k+jj*m)*gm+(ii*m+l)] = 0.0; 
		for(h = 0; h < m; h++){
		  P[(k+jj*m)*gm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(knotsCoordsD[jj*g+ii], phi[h], nu[h], covModel);
		}
	      }
	    }
	  }
	}
	// } //parallel for
	
	//get D, det(D), and D^{-1/2}
	det = 0;
	
	if(modPP){
	  
	  for(k = 0; k < gm; k++){
	    for(l = k; l < gm; l++){
	      tmp_gmgm[k*gm+l] = K[k*gm+l];
	    }
	  }
	  
	  F77_NAME(dpotrf)(lower, &gm, tmp_gmgm, &gm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}//L_K
	  
	  F77_NAME(dcopy)(&nmgm, P, &incOne, H, &incOne);
	  
	  F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &gm, &nm, &one, tmp_gmgm, &gm, H, &gm);
	  
	  for(k = 0; k < n; k++){
	    F77_NAME(dgemm)(ytran, ntran, &m, &m, &gm, &one, &H[k*gm*m], &gm, &H[k*gm*m], &gm, &zero, tmp_mm, &m);
	    
	    for(l = 0; l < mm; l++){
	      tmp_mm[l] = Psi[l] + V[l] - tmp_mm[l]; 
	    }
	    
	    F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	    
	    for(l = 0; l < m; l++){
	  	det += 2*log(tmp_mm[l*m+l]); //det(D)
	      }

	      F77_NAME(dtrtri)(lower, nUnit, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dtrtri failed\n");}

	      F77_NAME(dcopy)(&mm, tmp_mm, &incOne, &D[k*mm], &incOne); //D^{-1/2}
	    }

	  }else{
	    F77_NAME(dcopy)(&mm, Psi, &incOne, tmp_mm, &incOne);

	    if(PsiDiag){
	      for(k = 0; k < m; k++){
	  	det += log(tmp_mm[k*m+k]);
	  	tmp_mm[k*m+k] = 1.0/sqrt(tmp_mm[k*m+k]);
	      }
	      det *= n;
	    }else{
	      F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	      for(k = 0; k < m; k++){
	  	det += 2*log(tmp_mm[k*m+k]);
	      }
	      det *= n;
	      F77_NAME(dtrtri)(lower, nUnit, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dtrtri failed\n");}
	    }

	    for(k = 0; k < n; k++){
	      F77_NAME(dcopy)(&mm, tmp_mm, &incOne, &D[k*mm], &incOne); //D^{-1/2}
	    }
	  }

	  F77_NAME(dcopy)(&nmgm, P, &incOne, H, &incOne);//note, copy needed only in amcmc case, this could be dropped later
	  for(k = 0; k < n; k++){
	    F77_NAME(dtrmm)(rside, lower, ytran, nUnit, &gm, &m, &one, &D[k*mm], &m, &H[k*m*gm], &gm);//W'
	  }
	  
	  F77_NAME(dgemm)(ntran, ytran, &gm, &gm, &nm, &one, H, &gm, H, &gm, &zero, tmp_gmgm, &gm);//W'W

	  for(k = 0; k < gm; k++){
	    for(l = k; l < gm; l++){
	      K[k*gm+l] += tmp_gmgm[k*gm+l];
	    }
	  }

	  F77_NAME(dpotrf)(lower, &gm, K, &gm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}//L

	  F77_NAME(dtrsm)(lside, lower, ntran, nUnit, &gm, &nm, &one, K, &gm, H, &gm);//LH = W'

	  F77_NAME(dcopy)(&nm, u, &incOne, tmp_nm, &incOne);//note, copy needed only in amcmc case, this could be dropped later
	  for(k = 0; k < n; k++){
	    F77_NAME(dtrmv)(lower, ntran, nUnit, &m, &D[k*mm], &m, &tmp_nm[k*m], &incOne);//v
	  }

	  F77_NAME(dgemv)(ntran, &gm, &nm, &one, H, &gm, tmp_nm, &incOne, &zero, tmp_gm, &incOne); //w = Hv
	  
	  Q = F77_NAME(ddot)(&nm, tmp_nm, &incOne, tmp_nm, &incOne) - F77_NAME(ddot)(&gm, tmp_gm, &incOne, tmp_gm, &incOne);
	  
	  F77_NAME(dgemm)(ntran, ytran, &gm, &gm, &nm, &negOne, H, &gm, H, &gm, &zero, tmp_gmgm, &gm);//-HH'
	  
	  for(k = 0; k < gm; k++){
	    tmp_gmgm[k*gm+k] += 1.0; //J check
	  }

	  F77_NAME(dpotrf)(lower, &gm, tmp_gmgm, &gm, &info); if(info != 0){error("c++ error: dpotrf failed\n");}//L_J
	  
	  for(k = 0; k < gm; k++){
	    det -= 2*log(tmp_gmgm[k*gm+k]);
	  }

	  //
	  //priors, jacobian adjustments, and likelihood
	  //
	  priors = 0.0;
	  
	  if(KPriorName == "IW"){
	    logDetK = 0.0;
	    SKtrace = 0.0;
	    
	    for(k = 0; k < m; k++){logDetK += 2*log(A[k*m+k]);}
	    
	    //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	    for(k = 0; k < m; k++){priors += (m-k)*log(A[k*m+k])+log(A[k*m+k]);}
	    
	    //S*K^-1
	    F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	    F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
	    for(k = 0; k < m; k++){SKtrace += tmp_mm[k*m+k];}
	    priors += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;
	  }else{	     
	    for(k = 0; k < nLTr; k++){
	      priors += dnorm(params[AIndx+k], ANormMu[k], sqrt(ANormC[k]), 1);
	    }
	  }
	  
	  if(nugget){
	    if(PsiDiag){
	      for(k = 0; k < m; k++){
	  	priors += -1.0*(1.0+PsiIGa[k])*log(Psi[k*m+k])-PsiIGb[k]/Psi[k*m+k]+log(Psi[k*m+k]);
	      }
	    }else{
	      if(PsiPriorName == "IW"){
	  	logDetK = 0.0;
	  	SKtrace = 0.0; 
	
	  	for(k = 0; k < m; k++){logDetK += 2*log(L[k*m+k]);}
		
	  	//jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	  	for(k = 0; k < m; k++){priors += (m-k)*log(L[k*m+k])+log(L[k*m+k]);}
		
	  	//get S*K^-1
	  	F77_NAME(dpotri)(lower, &m, L, &m, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	  	F77_NAME(dsymm)(rside, lower, &m, &m, &one, L, &m, PsiIW_S, &m, &zero, tmp_mm, &m);
	  	for(k = 0; k < m; k++){SKtrace += tmp_mm[k*m+k];}
	  	priors += -0.5*(PsiIW_df+m+1)*logDetK - 0.5*SKtrace;
	      }else{
	  	for(k = 0; k < nLTr; k++){
	  	  priors += dnorm(params[LIndx+k], LNormMu[k], sqrt(LNormC[k]), 1);
	  	}
	      }
	    }
	  }
	  
	  for(k = 0; k < m; k++){
	    priors += log(phi[k] - phiUnif[k*2]) + log(phiUnif[k*2+1] - phi[k]); 
	    
	    if(covModel == "matern"){
	      priors += log(nu[k] - nuUnif[k*2]) + log(nuUnif[k*2+1] - nu[k]);  
	    }
	  }
	  
	  logPostCand = priors-0.5*det-0.5*Q;

	  //
	  //MH accept/reject	
	  //      
	  logMHRatio = logPostCand - logPostCurrent;
	  
	  if(runif(0.0,1.0) <= exp(logMHRatio)){
	    logPostCurrent = logPostCand;

	    F77_NAME(dcopy)(&nmm, D, &incOne, DCurrent, &incOne);
	    F77_NAME(dcopy)(&nmgm, H, &incOne, HCurrent, &incOne);
	    detCurrent = det;
	    priorsCurrent = priors;

	    //set to true so beta's mu and var are updated
	    updateBeta = true;

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
	
	//only update beta's mu and var if theta has changed
	if(updateBeta){
	  
	  F77_NAME(dcopy)(&nm, Y, &incOne, tmp_nm, &incOne);
	  F77_NAME(dcopy)(&nmp, X, &incOne, tmp_nmp, &incOne);

	  for(k = 0; k < n; k++){
	    F77_NAME(dtrmv)(lower, ntran, nUnit, &m, &DCurrent[k*mm], &m, &tmp_nm[k*m], &incOne);//\tilda{y}
	    F77_NAME(dtrmm)(rside, lower, ytran, nUnit, &p, &m, &one, &DCurrent[k*mm], &m, &tmp_nmp[k*m*p], &p);//V 
	  }
  
	  F77_NAME(dgemm)(ntran, ytran, &gm, &p, &nm, &one, HCurrent, &gm, tmp_nmp, &p, &zero, tmp_gmp, &gm);//\tilda{V}
	  F77_NAME(dgemm)(ntran, ytran, &p, &p, &nm, &one, tmp_nmp, &p, tmp_nmp, &p, &zero, tmp_pp, &p);//V'V
	  F77_NAME(dgemm)(ytran, ntran, &p, &p, &gm, &one, tmp_gmp, &gm, tmp_gmp, &gm, &zero, tmp_pp2, &p);//\tilda{V}'\tilda{V}
	  
	  for(k = 0; k < p; k++){
	    for(l = k; l < p; l++){
	      tmp_pp[k*p+l] -= tmp_pp2[k*p+l];//Z
	      
	      if(betaPrior == "normal"){
		tmp_pp[k*p+l] += betaCInv[k*p+l];
	      }
	      
	    }
	  }
 	  
	  //B
	  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	  F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
   
	  //b
	  F77_NAME(dgemv)(ntran, &p, &nm, &one, tmp_nmp, &p, tmp_nm, &incOne, &zero, tmp_p, &incOne); //V'\tilda{y}
	  F77_NAME(dgemv)(ntran, &gm, &nm, &one, HCurrent, &gm, tmp_nm, &incOne, &zero, tmp_gm, &incOne); //H\tilda{y}
	  F77_NAME(dgemv)(ytran, &gm, &p, &one, tmp_gmp, &gm, tmp_gm, &incOne, &zero, tmp_pp2, &incOne); //\tilda{V}'(H\tilda{y})
	  
	  for(k = 0; k < p; k++){
	    tmp_p[k] -= tmp_pp2[k]; 
	    
	    if(betaPrior == "normal"){
	      tmp_p[k] += betaCInvMu[k];
	    }
	  }
	  
	  F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &incOne, &zero, tmp_p2, &incOne); //Bb
	  F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}

	}//end updateBeta

	//set to false so beta's mu and var are only updated when theta has changed
       	updateBeta = false;
	
	//draw beta
	mvrnorm(beta, tmp_p2, tmp_pp, p, false);//note, tmp_p2 and tmp_pp will carry over when theta is not updated
	
	//update logPostCurrent (beta changes on every iteration so logPostCurrent must be updated)
	F77_NAME(dgemv)(ytran, &p, &nm, &negOne, X, &p, beta, &incOne, &zero, u, &incOne);
	F77_NAME(daxpy)(&nm, &one, Y, &incOne, u, &incOne);//u
	  
	F77_NAME(dcopy)(&nm, u, &incOne, tmp_nm, &incOne);
	for(k = 0; k < n; k++){
	  F77_NAME(dtrmv)(lower, ntran, nUnit, &m, &DCurrent[k*mm], &m, &tmp_nm[k*m], &incOne);//v
	}
	  
	F77_NAME(dgemv)(ntran, &gm, &nm, &one, HCurrent, &gm, tmp_nm, &incOne, &zero, tmp_gm, &incOne); //w = Hv
	  
	Q = F77_NAME(ddot)(&nm, tmp_nm, &incOne, tmp_nm, &incOne) - F77_NAME(ddot)(&gm, tmp_gm, &incOne, tmp_gm, &incOne);

	logPostCurrent = priorsCurrent-0.5*detCurrent-0.5*Q;

	/******************************
               Save samples
	*******************************/
	F77_NAME(dcopy)(&p, beta, &incOne, &REAL(betaSamples_r)[s*p], &incOne);
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
      if(verbose){
      	if(status == nReport){
      	  if(amcmc){
      	    Rprintf("Batch: %i of %i, %3.2f%%\n", b, nBatch, 100.0*b/nBatch);
      	    Rprintf("\tparameter\tacceptance\ttuning\n");
      	    for(j = 0, i = 0; j < m; j++){
      	      for(k = j; k < m; k++, i++){
      		Rprintf("\tA[%i,%i]\t\t%3.1f%\t\t%1.2f\n", j+1, k+1, 100.0*REAL(accept_r)[b*nParams+AIndx+i], exp(tuning[AIndx+i]));
      	      }
      	    }
      	    if(nugget){
      	      if(PsiDiag){
      		for(j = 0; j < m; j++){
      		  Rprintf("\tPsi[%i,%i]\t%3.1f%\t\t%1.2f\n", j+1, j+1, 100.0*REAL(accept_r)[b*nParams+LIndx+j], exp(tuning[LIndx+j]));
      		}
      	      }else{
      		for(j = 0, i = 0; j < m; j++){
      		  for(k = j; k < m; k++, i++){
      		    Rprintf("\tL[%i,%i]\t\t%3.1f%\t\t%1.2f\n", j+1, k+1, 100.0*REAL(accept_r)[b*nParams+LIndx+i], exp(tuning[LIndx+i]));
      		  }
      		}
      	      }
      	    }
      	    for(j = 0; j < m; j++){
      	      Rprintf("\tphi[%i]\t\t%3.1f%\t\t%1.2f\n", j+1, 100.0*REAL(accept_r)[b*nParams+phiIndx+j], exp(tuning[phiIndx+j]));
      	    }
      	    if(covModel == "matern"){
      	      Rprintf("\n");
      	      for(j = 0; j < m; j++){
      		Rprintf("\tnu[%i]\t\t%3.1f%\t\t%1.2f\n", j+1, 100.0*REAL(accept_r)[b*nParams+nuIndx+j], exp(tuning[nuIndx+j]));
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
      	  status = 0;
      	  batchAccept = 0;
      	}
      }
      status++;
      
    }//end samples
    
    PutRNGstate();
    
    //untransform variance variables
    for(s = 0; s < nSamples; s++){
      
      covTransInv(&REAL(samples_r)[s*nParams+AIndx], &REAL(samples_r)[s*nParams+AIndx], m);
      
      if(nugget){
	if(PsiDiag){
	  for(i = 0; i < m; i++){
	    REAL(samples_r)[s*nParams+LIndx+i] = exp(REAL(samples_r)[s*nParams+LIndx+i]);
	  }
	}else{
	  covTransInv(&REAL(samples_r)[s*nParams+LIndx], &REAL(samples_r)[s*nParams+LIndx], m);
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
    int nResultListObjs = 4;
        
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    //samples
    SET_VECTOR_ELT(result_r, 0, samples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.theta.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, accept_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("acceptance"));
    
    SET_VECTOR_ELT(result_r, 2, tuning_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("tuning"));
    
    SET_VECTOR_ELT(result_r, 3, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 3, mkChar("p.beta.samples"));
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}


  
