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
#include "covInvDet.h"

extern"C" {
  
  void dbdimm_(int *transa, int *mb, int *n, int *kb, const double *alpha, int *descra,
	       double *val, int *blda, int *ibdiag, int *nbdiag, int *lb,
	       double *b, int *ldb, const double *beta, double *c, int *ldc, double *work, int *lwork);
}

extern "C" {

  SEXP spPPMvLM(SEXP Y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coordsD_r,
		SEXP isModPp_r, SEXP q_r, SEXP knotsD_r, SEXP coordsKnotsD_r,
		SEXP KIW_r, SEXP PsiIW_r, SEXP nuUnif_r, SEXP phiUnif_r,
		SEXP phiStarting_r, SEXP AStarting_r, SEXP LStarting_r, SEXP nuStarting_r, SEXP betaStarting_r,
		SEXP phiTuning_r, SEXP ATuning_r, SEXP LTuning_r, SEXP nuTuning_r, 
		SEXP nugget_r, SEXP covModel_r, SEXP nSamples_r, SEXP verbose_r, SEXP nReport_r){


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
    int m = INTEGER(m_r)[0];
    int nLTr = m*(m-1)/2+m;
    int q = INTEGER(q_r)[0];
    int nn = n*n;
    int mm = m*m;
    int nm = n*m;
    int nq = n*q;
    int qq = q*q;
    int qm = q*m;
    int qmqm = qm*qm;

    double *coordsD = REAL(coordsD_r);

    //covariance model
    string covModel = CHAR(STRING_ELT(covModel_r,0));

    //if predictive process
    bool isModPp = static_cast<bool>(INTEGER(isModPp_r)[0]);

    double *knotsD = REAL(knotsD_r);
    double *coordsKnotsD = REAL(coordsKnotsD_r);

    //priors and starting
    double *phiUnif = REAL(phiUnif_r);
    double KIW_df = REAL(VECTOR_ELT(KIW_r, 0))[0]; double *KIW_S = REAL(VECTOR_ELT(KIW_r, 1));

    double *phiStarting = REAL(phiStarting_r);
    double *AStarting = REAL(AStarting_r);
    double *betaStarting = REAL(betaStarting_r);

    //if nugget
    bool nugget = static_cast<bool>(INTEGER(nugget_r)[0]);

    double *LStarting = NULL;
    double PsiIW_df = 0; double *PsiIW_S = NULL;
    if(nugget){
      LStarting = REAL(LStarting_r);
      PsiIW_df = REAL(VECTOR_ELT(PsiIW_r, 0))[0]; PsiIW_S = REAL(VECTOR_ELT(PsiIW_r, 1));
    }

    //if matern
    double *nuUnif = NULL;
    double *nuStarting = NULL;

    if(covModel == "matern"){
      nuUnif = REAL(nuUnif_r);
      nuStarting = REAL(nuStarting_r);
    }

    //tuning
    double *phiTuning = REAL(phiTuning_r);
    double *ATuning = REAL(ATuning_r);
    double *LTuning = NULL;
    double *nuTuning = NULL;

    if(nugget)
      LTuning = REAL(LTuning_r);

    if(covModel == "matern")
      nuTuning = REAL(nuTuning_r);

    int nSamples = INTEGER(nSamples_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];

    double *A = (double *) R_alloc(mm, sizeof(double));
    double *L = (double *) R_alloc(mm, sizeof(double));

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

      if(!nugget && !isModPp){
	Rprintf("Warning: setting Psi = Diag[0.01].\n\n");
      }

      Rprintf("Priors and hyperpriors:\n");
      Rprintf("\tbeta flat.\n\n");   
      Rprintf("\tK IW hyperpriors df=%.5f\n", KIW_df);
      Rprintf("\t"); printMtrx(KIW_S, m, m);
      Rprintf("\n"); 

      if(nugget){
	Rprintf("\tPsi IW hyperpriors df=%.5f\n", PsiIW_df);
	Rprintf("\t"); printMtrx(PsiIW_S, m, m);
	Rprintf("\n"); 
      }
 
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

      Rprintf("Metropolis tuning values:\n");
  
      Rprintf("\tA tuning:\n");
      Rprintf("\t"); printVec(ATuning, nLTr);
      Rprintf("\n"); 

      if(nugget){
	Rprintf("\tL tuning:\n");
	Rprintf("\t"); printVec(LTuning, nLTr);
	Rprintf("\n"); 
      }

      Rprintf("\tphi tuning\n");
      Rprintf("\t"); printVec(phiTuning, m);
      Rprintf("\n");   

      if(covModel == "matern"){
	Rprintf("\tnu tuning\n");
	Rprintf("\t"); printVec(nuTuning, m);
	Rprintf("\n");
      }

      Rprintf("Metropolis starting values:\n");
      
      covExpand(AStarting, A, m);
      Rprintf("\tA starting\n");
      printMtrx(A, m, m);
      Rprintf("\n");

      if(nugget){
	covExpand(LStarting, L, m);
	Rprintf("\tL starting\n");
	printMtrx(L, m, m);
	Rprintf("\n");
      }

      Rprintf("\tphi starting\n");
      Rprintf("\t"); printVec(phiStarting, m);
      Rprintf("\n");   

      if(covModel == "matern"){
	Rprintf("\tnu starting\n");
	Rprintf("\t"); printVec(nuStarting, m);
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
    int nSpParams, AIndx, LIndx, phiIndx, nuIndx;

    if(!nugget && covModel != "matern"){
      nSpParams = nLTr+m;//A, phi
      AIndx = 0; phiIndx = nLTr;
    }else if(nugget && covModel != "matern"){
      nSpParams = 2*nLTr+m;//A, L, phi
      AIndx = 0; LIndx = nLTr; phiIndx = LIndx+nLTr;
    }else if(!nugget && covModel == "matern"){
      nSpParams = nLTr+2*m;//A, phi, nu
      AIndx = 0; phiIndx = nLTr, nuIndx = phiIndx+m;
    }else{
      nSpParams = 2*nLTr+2*m;//A, Phi, phi, nu
      AIndx = 0; LIndx = nLTr, phiIndx = LIndx+nLTr, nuIndx = phiIndx+m;
    }
    
    double *spParams = (double *) R_alloc(nSpParams, sizeof(double));
    
    //set starting
    covTrans(AStarting, &spParams[AIndx], m);

    if(nugget){
      covTrans(LStarting, &spParams[LIndx], m);
    }

    for(i = 0; i < m; i++){
      spParams[phiIndx+i] = logit(phiStarting[i], phiUnif[i*2], phiUnif[i*2+1]);
    }

    if(covModel == "matern"){
      for(i = 0; i < m; i++){
	spParams[nuIndx+i] = logit(nuStarting[i], nuUnif[i*2], nuUnif[i*2+1]);
      }
    }

    //Beta parameter and set starting
    double *beta = (double *) R_alloc(p, sizeof(double));
    F77_NAME(dcopy)(&p, betaStarting, &incOne, beta, &incOne);

    //samples and random effects
    int nParams = p+nSpParams;

    SEXP w_r, w_str_r, samples_r, accept_r;

    PROTECT(w_r = allocMatrix(REALSXP, nm, nSamples)); nProtect++; 
    double *w = REAL(w_r); zeros(w, nm*nSamples);

    PROTECT(w_str_r = allocMatrix(REALSXP, qm, nSamples)); nProtect++; 
    double *w_str = REAL(w_str_r); zeros(w_str, qm*nSamples);

    PROTECT(samples_r = allocMatrix(REALSXP, nParams, nSamples)); nProtect++; 
    double *samples = REAL(samples_r);

    PROTECT(accept_r = allocMatrix(REALSXP, 1, 1)); nProtect++;

    /*****************************************
       Set-up MCMC alg. vars. matrices etc.
    *****************************************/
    int s=0, status=0, rtnStatus=0, accept=0, batchAccept = 0;
    double logPostCurrent = 0, logPostCand = 0, detCand = 0, det = 0,
      logDetK, SKtrace;

    double *C = (double *) R_alloc(nm*nm, sizeof(double)); 
    double *ct = (double *) R_alloc(nm*qm, sizeof(double));
    double *C_str = (double *) R_alloc(qm*qm, sizeof(double));
    double *E = (double *) R_alloc(m*m*n, sizeof(double));

    //for sparse mult
    int ytranInt = 1;
    int ntranInt = 0;
    
    int mb = n; 
    int kb = n; //Number of block columns in matrix A
    int blda = n; //leading block dimension of val(), i.e., block leading dimension
    int nbdiag = 1; //the number of non-zero block diagonals in A.
    int lb = m;//dimension of dense blocks composing A
    
    int lwork = mb*lb*lb; if(mb*lb*lb < qm) lwork = qm;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    
    //for main diagonal block matrix only!
    int ibdiag = 0;
    int descra = 0;
    
    //Other stuff
    double * tmp_nmnm = NULL;
    if(isModPp){//not needed for the non-modified pp
      tmp_nmnm = (double *) R_alloc(nm*nm, sizeof(double));
    }

    double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm1 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_mm2 = (double *) R_alloc(mm, sizeof(double));
    double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
    double *tmp_nm1 = (double *) R_alloc(nm, sizeof(double));

    double *tmp_qm = (double *) R_alloc(qm, sizeof(double));
    double *tmp_qm1 = (double *) R_alloc(qm, sizeof(double));

    double *tmp_nmqm = (double *) R_alloc(nm*qm, sizeof(double));
    double *tmp_nmqm1 = (double *) R_alloc(nm*qm, sizeof(double));
    double *tmp_qmqm = (double *) R_alloc(qm*qm, sizeof(double));
    double *tmp_qmqm1 = (double *) R_alloc(qm*qm, sizeof(double));

    double *candSpParams = (double *) R_alloc(nSpParams, sizeof(double));
    double *K = (double *) R_alloc(mm, sizeof(double));
    double *Psi = (double *) R_alloc(mm, sizeof(double)); 
    double *theta = (double *) R_alloc(mm, sizeof(double));
    double logMHRatio;
    
    double *S_beta = (double *) R_alloc(p*p, sizeof(double));
    double *Mu_beta = (double *) R_alloc(p, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double)); 
    double *tmp_nmp = (double *) R_alloc(nm*p, sizeof(double)); 
    int brute = 0;

    //fix nugget if needed
    if(!nugget){
      zeros(Psi, mm);

      if(!isModPp)
	for(i = 0; i < m; i++) Psi[i*m+i] = 0.01;
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
    for(s = 0; s < nSamples; s++){
      
      //
      //Current
      //
      covTransInvExpand(&spParams[AIndx], A, m);

      if(nugget){
       covTransInvExpand(&spParams[LIndx], L, m);
      }      

      for(i = 0; i < m; i++){
	theta[i] = logitInv(spParams[phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
	if(covModel == "matern"){
	  theta[m+i] = logitInv(spParams[nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
	}
      }

      //K = A'A and Psi = L'L
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, A, &m, A, &m, &zero, K, &m);
      
      if(nugget){
	F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, L, &m, L, &m, &zero, Psi, &m);
      }     

      if(isModPp){  
	det = mvMPPCovInvDet(knotsD, coordsKnotsD, C, C_str, ct, E, n, m, q, Psi, K, theta, 
			     tmp_mm, tmp_mm1, tmp_mm2, tmp_nmqm, tmp_nmqm1, tmp_qmqm, tmp_qmqm1, lwork, work,
			     covModel, nPramPtr, covModelObj, cov1ParamPtr, cov2ParamPtr, brute);

      }else{
	det = mvPPCovInvDet(knotsD, coordsKnotsD, C, C_str, ct, E, n, m, q, Psi, K, theta, 
			    tmp_mm, tmp_mm1, tmp_mm2, tmp_nmqm, tmp_nmqm1, tmp_qmqm, tmp_qmqm1, lwork, work,
			    covModel, nPramPtr, covModelObj, cov1ParamPtr, cov2ParamPtr, brute);

      }

      //
      //Update Beta
      //
      //finish the Gibbs
      F77_NAME(dsymm)(lside, lower, &nm, &p, &one, C, &nm, X, &nm, &zero, tmp_nmp, &nm);
      F77_NAME(dgemm)(ytran, ntran, &p, &p, &nm, &one, X, &nm, tmp_nmp, &nm, &zero, S_beta, &p);
      
      F77_NAME(dpotrf)(lower, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
      F77_NAME(dpotri)(lower, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky inverse failed\n" << endl;}
      
      F77_NAME(dsymv)(lower, &nm, &one, C, &nm, Y, &incOne, &zero, tmp_nm, &incOne);
      F77_NAME(dgemv)(ytran, &nm, &p, &one, X, &nm, tmp_nm, &incOne, &zero, tmp_p, &incOne);
      F77_NAME(dsymv)(lower, &p, &one, S_beta, &p, tmp_p, &incOne, &zero, Mu_beta, &incOne);
      
      //Gibbs draw
      //take lower for the chol for the mv draw
      F77_NAME(dpotrf)(lower, &p, S_beta, &p, &info); if(info != 0){cout << "c++ error: Cholesky failed\n" << endl;}
      mvrnorm(beta, Mu_beta, S_beta, p, false);

      //
      //Likelihood with Jacobian   
      //
      logPostCurrent = 0.0;
      
      //
      //Jacobian and IW priors for K = A'A and Psi = L'L
      //
      //A'A prior with jacob.
      logDetK = 0.0;
      SKtrace = 0.0;

      for(i = 0; i < m; i++) logDetK += 2*log(A[i*m+i]);
      
      //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
      for(i = 0; i < m; i++) logPostCurrent += (m-i)*log(A[i*m+i])+log(A[i*m+i]);
      
      //get S*K^-1, already have the chol of K (i.e., A)
      F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){cout << "c++ error: A Cholesky inverse failed\n" << endl;}
      F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
      for(i = 0; i < m; i++){SKtrace += tmp_mm[i*m+i];}
      logPostCurrent += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;

      if(nugget){
	
	//L'L prior with jacob.
	logDetK = 0.0;
	SKtrace = 0.0; 
	
	for(i = 0; i < m; i++) logDetK += 2*log(L[i*m+i]);
	
	//jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	for(i = 0; i < m; i++) logPostCurrent += (m-i)*log(L[i*m+i])+log(L[i*m+i]);
	
	//get S*K^-1, already have the chol of Psi (i.e., L)
	F77_NAME(dpotri)(lower, &m, L, &m, &info); if(info != 0){cout << "c++ error: L Cholesky inverse failed\n" << endl;}
	F77_NAME(dsymm)(rside, lower, &m, &m, &one, L, &m, PsiIW_S, &m, &zero, tmp_mm, &m);
	for(i = 0; i < m; i++){SKtrace += tmp_mm[i*m+i];}
	logPostCurrent += -0.5*(PsiIW_df+m+1)*logDetK - 0.5*SKtrace;

      }

      for(i = 0; i < m; i++){
	logPostCurrent += log(theta[i] - phiUnif[i*2]) + log(phiUnif[i*2+1] - theta[i]); 
      
	if(covModel == "matern"){
	  logPostCurrent += log(theta[m+i] - nuUnif[i*2]) + log(nuUnif[i*2+1] - theta[m+i]);  
	}
      }

      //Y-XB
      F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
      F77_NAME(daxpy)(&nm, &one, Y, &incOne, tmp_nm, &incOne);
      
      //(-1/2) * (Y-XB)` * C * (Y-XB)
      F77_NAME(dsymv)(lower, &nm, &one, C, &nm, tmp_nm, &incOne, &zero, tmp_nm1, &incOne);
      logPostCurrent += -0.5*det-0.5*F77_NAME(ddot)(&nm, tmp_nm, &incOne, tmp_nm1, &incOne);

      //
      //Candidate
      //

      //propose   
      for(i = 0; i < nLTr; i++){
	candSpParams[AIndx+i] = rnorm(spParams[AIndx+i], ATuning[i]);

	if(nugget){
	  candSpParams[LIndx+i] = rnorm(spParams[LIndx+i], LTuning[i]);
	}
      }

      covTransInvExpand(&candSpParams[AIndx], A, m);

      if(nugget){
	covTransInvExpand(&candSpParams[LIndx], L, m);
      }

      for(i = 0; i < m; i++){
	candSpParams[phiIndx+i] = rnorm(spParams[phiIndx+i], phiTuning[i]);
	theta[i] = logitInv(candSpParams[phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
	if(covModel == "matern"){
	  candSpParams[nuIndx+i] = rnorm(spParams[nuIndx+i], nuTuning[i]);
	  theta[m+i] = logitInv(candSpParams[nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
	}
      }
      
      //K = A'A and Psi = L'L
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, A, &m, A, &m, &zero, K, &m);

      if(nugget){
	F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, L, &m, L, &m, &zero, Psi, &m);
      }

      if(isModPp){  
	det = mvMPPCovInvDet(knotsD, coordsKnotsD, C, C_str, ct, E, n, m, q, Psi, K, theta, 
			     tmp_mm, tmp_mm1, tmp_mm2, tmp_nmqm, tmp_nmqm1, tmp_qmqm, tmp_qmqm1, lwork, work,
			     covModel, nPramPtr, covModelObj, cov1ParamPtr, cov2ParamPtr, brute);
      }else{
	det = mvPPCovInvDet(knotsD, coordsKnotsD, C, C_str, ct, E, n, m, q, Psi, K, theta, 
			    tmp_mm, tmp_mm1, tmp_mm2, tmp_nmqm, tmp_nmqm1, tmp_qmqm, tmp_qmqm1, lwork, work,
			    covModel, nPramPtr, covModelObj, cov1ParamPtr, cov2ParamPtr, brute);
      }
      
      //
      //Likelihood with Jacobian   
      //
      logPostCand = 0.0;
      
      //
      //Jacobian and IW priors for K = A'A and Psi = L'L
      //
      //AtA prior with jacob.
      logDetK = 0.0;
      SKtrace = 0.0;

      for(i = 0; i < m; i++) logDetK += 2*log(A[i*m+i]);
      
      //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
      for(i = 0; i < m; i++) logPostCand += (m-i)*log(A[i*m+i])+log(A[i*m+i]);
      
      //get S*K^-1, already have the chol of K (i.e., A)
      F77_NAME(dpotri)(lower, &m, A, &m, &info); if(info != 0){cout << "c++ error: Cand A Cholesky inverse failed\n" << endl;}
      F77_NAME(dsymm)(rside, lower, &m, &m, &one, A, &m, KIW_S, &m, &zero, tmp_mm, &m);
      for(i = 0; i < m; i++){SKtrace += tmp_mm[i*m+i];}
      logPostCand += -0.5*(KIW_df+m+1)*logDetK - 0.5*SKtrace;

      if(nugget){
	
	//L'L prior with jacob.
	logDetK = 0.0;
	SKtrace = 0.0; 

	for(i = 0; i < m; i++) logDetK += 2*log(L[i*m+i]);
	
	//jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
	for(i = 0; i < m; i++) logPostCand += (m-i)*log(L[i*m+i])+log(L[i*m+i]);
	
	//get S*K^-1, already have the chol of Psi (i.e., L)
	F77_NAME(dpotri)(lower, &m, L, &m, &info); if(info != 0){cout << "c++ error: Cand L Cholesky inverse failed\n" << endl;}
	F77_NAME(dsymm)(rside, lower, &m, &m, &one, L, &m, PsiIW_S, &m, &zero, tmp_mm, &m);
	for(i = 0; i < m; i++){SKtrace += tmp_mm[i*m+i];}
	logPostCand += -0.5*(PsiIW_df+m+1)*logDetK - 0.5*SKtrace;

      }

      for(i = 0; i < m; i++){
	logPostCand += log(theta[i] - phiUnif[i*2]) + log(phiUnif[i*2+1] - theta[i]); 
      
	if(covModel == "matern"){
	  logPostCand += log(theta[m+i] - nuUnif[i*2]) + log(nuUnif[i*2+1] - theta[m+i]);  
	}
      }

      //Y-XB
      F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
      F77_NAME(daxpy)(&nm, &one, Y, &incOne, tmp_nm, &incOne);
      
      //(-1/2) * (Y-XB)` * C * (Y-XB)
      F77_NAME(dsymv)(lower, &nm, &one, C, &nm, tmp_nm, &incOne, &zero, tmp_nm1, &incOne);
      logPostCand += -0.5*det-0.5*F77_NAME(ddot)(&nm, tmp_nm, &incOne, tmp_nm1, &incOne);

      //
      //MH accept/reject	
      //      
  
      //MH ratio with adjustment
      logMHRatio = logPostCand - logPostCurrent;
      
      if(runif(0.0,1.0) <= exp(logMHRatio)){
	F77_NAME(dcopy)(&nSpParams, candSpParams, &incOne, spParams, &incOne);
	accept++;
	batchAccept++;
      }
    
      /******************************
         Recover w and w* if needed
      *******************************/
      //
      //Current
      //
      covTransInvExpand(&spParams[AIndx], A, m);

      if(nugget){
       covTransInvExpand(&spParams[LIndx], L, m);
      }      

      for(i = 0; i < m; i++){
	theta[i] = logitInv(spParams[phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
	if(covModel == "matern"){
	  theta[m+i] = logitInv(spParams[nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
	}
      }

      //K = A'A and Psi = L'L
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, A, &m, A, &m, &zero, K, &m);
      
      if(nugget){
	F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, L, &m, L, &m, &zero, Psi, &m);
      }
   
      //make ct
      for(i = 0; i < n; i++){
	for(j = 0; j < q; j++){
	  
	  zeros(tmp_mm, mm);
	  
	  for(k = 0; k < m; k++){
	    if(nPramPtr == 1)
	      (covModelObj->*cov1ParamPtr)(theta[k], tmp_mm[k*m+k], coordsKnotsD[j*n+i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(theta[k], theta[m+k], tmp_mm[k*m+k], coordsKnotsD[j*n+i]);
	  }
	  
	  F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, A, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	  F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, A, &m, &zero, tmp_mm, &m);
	  
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      ct[((j*m+l)*nm)+(i*m+k)] = tmp_mm[l*m+k];
	      tmp_mm[l*m+k] = 0.0; //zero out
	    }
	  }
	}
      }
      
      
      //
      //make C_str
      //
      for(i = 0; i < q; i++){
	for(j = 0; j < q; j++){
	  
	  zeros(tmp_mm, mm);
	  
	  for(k = 0; k < m; k++){
	    if(nPramPtr == 1)
	      (covModelObj->*cov1ParamPtr)(theta[k], tmp_mm[k*m+k], knotsD[j*q+i]);
	    else //i.e., 2 parameter matern
	      (covModelObj->*cov2ParamPtr)(theta[k], theta[m+k], tmp_mm[k*m+k], knotsD[j*q+i]);
	  }
	  
	  F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, A, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
	  F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, A, &m, &zero, tmp_mm, &m);
	  
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      C_str[((j*m+l)*qm)+(i*m+k)] = tmp_mm[l*m+k];
	      tmp_mm[l*m+k] = 0.0; //zero out
	    }
	  }
	}
      }
      

      F77_NAME(dpotrf)(lower, &qm, C_str, &qm, &info); if(info != 0){error("c++ error: recover w Cholesky failed\n");}
      F77_NAME(dpotri)(lower, &qm, C_str, &qm, &info); if(info != 0){error("c++ error: recover w Cholesky failed\n");}

      //ct C_str^{-1}
      F77_NAME(dsymm)(rside, lower, &nm, &qm, &one, C_str, &qm, ct, &nm, &zero, tmp_nmqm, &nm);
      
      if(isModPp){
	
	//C = ct C_str^{-1} t(ct)
	F77_NAME(dgemm)(ntran, ytran, &nm, &nm, &qm, &one, tmp_nmqm, &nm, ct, &nm, &zero, C, &nm);
	
	//make E = (Psi + V - blk[C])^{-1}
	for(i = 0, j = 0; i < n; i++){
	  
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      tmp_mm[l*m+k] = Psi[l*m+k]+K[l*m+k]-C[(i*m+l)*nm+(i*m+k)];
	    }
	  }
	  
	  //must must be a full upper/lower for later.
	  mtrxInv(tmp_mm, m, info);
	  
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      E[j] = tmp_mm[l*m+k];
	      j++;
	    }
	  } 
	}
	
      }else{
	
	//must must be a full upper/lower for later.
	F77_NAME(dcopy)(&mm, Psi, &incOne, tmp_mm, &incOne);
	mtrxInv(tmp_mm, m, info);
	
	//this is a bit wasteful but makes life easier below.
	for(i = 0, j = 0; i < n; i++){
	  for(k = 0; k < m; k++){
	    for(l = 0; l < m; l++){
	      E[j] = tmp_mm[l*m+k];
	      j++;
	    }
	  } 
	}
	
      }
      
      //make tmp_nmqm1 = E ct C_str^{-1}
      F77_NAME(dbdimm)(&ntranInt, &mb, &qm, &kb, &one, &descra,
		       E, &blda, &ibdiag, &nbdiag, &lb,
		       tmp_nmqm, &nm, &zero, tmp_nmqm1, &nm, work, &lwork);
      
      //tmp_qmqm = (C_str^{-1} + (Cstr^{-1}ct)' E ct C_str^{-1})^{-1}
      F77_NAME(dgemm)(ytran, ntran, &qm, &qm, &nm, &one, tmp_nmqm, &nm, tmp_nmqm1, &nm, &zero, tmp_qmqm, &qm);
      F77_NAME(daxpy)(&qmqm, &one, C_str, &incOne, tmp_qmqm, &incOne);
      F77_NAME(dpotrf)(lower, &qm, tmp_qmqm, &qm, &info); if(info != 0){error("c++ error: recover w Cholesky failed\n");}
      F77_NAME(dpotri)(lower, &qm, tmp_qmqm, &qm, &info); if(info != 0){error("c++ error: recover w Cholesky failed\n");}
      
      //C_str^{-1} ct' E (Y-XB)
      F77_NAME(dgemv)(ntran, &nm, &p, &negOne, X, &nm, beta, &incOne, &zero, tmp_nm, &incOne);
      F77_NAME(daxpy)(&nm, &one, Y, &incOne, tmp_nm, &incOne);
      F77_NAME(dgemv)(ytran, &nm, &qm, &one, tmp_nmqm1, &nm, tmp_nm, &incOne, &zero, tmp_qm, &incOne);
      
      //(C_str^{-1} + [Cstr^{-1}ct]' E ct C_str^{-1}) (C_str^{-1} ct' E (Y-XB))
      F77_NAME(dsymv)(lower, &qm, &one, tmp_qmqm, &qm, tmp_qm, &incOne, &zero, tmp_qm1, &incOne);	
      
      
      F77_NAME(dpotrf)(lower, &qm, tmp_qmqm, &qm, &info); if(info != 0){error("c++ error: recover w Cholesky failed\n");}
      mvrnorm(&w_str[s*qm], tmp_qm1, tmp_qmqm, qm, false);
      
      //get w
      F77_NAME(dgemv)(ntran, &nm, &qm, &one, tmp_nmqm, &nm, &w_str[s*qm], &incOne, &zero, &w[s*nm], &incOne);

      /******************************
          Save samples and report
      *******************************/
      F77_NAME(dcopy)(&p, beta, &incOne, &samples[s*nParams], &incOne);
      F77_NAME(dcopy)(&nSpParams, spParams, &incOne, &samples[s*nParams+p], &incOne);
      
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
    }
    Rprintf("-------------------------------------------------\n");
    #ifdef Win32
    R_FlushConsole();
    #endif

    //untransform variance variables
    for(s = 0; s < nSamples; s++){
 
      covTransInv(&samples[s*nParams+p+AIndx], &samples[s*nParams+p+AIndx], m);
     
      if(nugget){
	covTransInv(&samples[s*nParams+p+LIndx], &samples[s*nParams+p+LIndx], m);
      }
	
      for(i = 0; i < m; i++){
	samples[s*nParams+p+phiIndx+i] = logitInv(samples[s*nParams+p+phiIndx+i], phiUnif[i*2], phiUnif[i*2+1]);
	
	if(covModel == "matern"){
	  samples[s*nParams+p+nuIndx+i] = logitInv(samples[s*nParams+p+nuIndx+i], nuUnif[i*2], nuUnif[i*2+1]);
	}
      }
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
    SET_VECTOR_ELT(resultNames, 0, mkChar("p.samples")); 

    SET_VECTOR_ELT(result, 1, accept_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("acceptance"));
    
    SET_VECTOR_ELT(result, 2, w_r);
    SET_VECTOR_ELT(resultNames, 2, mkChar("sp.effects"));
    
    SET_VECTOR_ELT(result, 3, w_str_r);
    SET_VECTOR_ELT(resultNames, 3, mkChar("sp.effects.knots"));
    
    namesgets(result, resultNames);
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
    
    
  }
}
