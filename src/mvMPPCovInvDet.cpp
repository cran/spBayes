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

extern"C" {
  
  void dbdimm_(int *transa, int *mb, int *n, int *kb, double *alpha, int *descra,
	       double *val, int *blda, int *ibdiag, int *nbdiag, int *lb,
	       double *b, int *ldb, double *beta, double *c, int *ldc, double *work, int *lwork);
}

double mvMPPCovInvDet(double *knotsD, double *coordsKnotsD, double *C, double *C_str, double *ct, 
			       double *E, int n, int m, int q, double *Psi, double *V, double *theta, 
			       double *tmp_mm, double *tmp_mm1, double *tmp_mm2, double *tmp_nmqm, 
			       double *tmp_nmqm1, double *tmp_qmqm, double *tmp_qmqm1, int lwork, double *work,
			       string covModel, int nPramPtr, covmodel *covModelObj, 
			       void (covmodel::*cov1ParamPtr)(double, double &, double &),
			       void (covmodel::*cov2ParamPtr)(double, double, double &, double&),
			       int brute){
  
  //n is the number of coords
  //m is the number of response variables
  //q is the number of knots
  //for sparse E is m*m*n
  
  int nm = n*m;
  int mm = m*m;
  int nn = n*n;
  int qm = q*m;
  int nmnm = nm*nm;
  int qmqm = qm*qm;
  int i, j, k, l, info;
  double logDet = 0.0, logDetC_str = 0.0, logDetE = 0.0;

  //for sparse mult
  int ytranInt = 1;
  int ntranInt = 0;
  
  int mb = n; 
  int kb = n; //Number of block columns in matrix A
  int blda = n; //leading block dimension of val(), i.e., block leading dimension
  int nbdiag = 1; //the number of non-zero block diagonals in A.
  int lb = m;//dimension of dense blocks composing A
  
//   int lwork = mb*lb*lb; if(mb*lb*lb < qm) lwork = qm;
//   double *work = (double *) R_alloc(lwork, sizeof(double));
  
  //for main diagonal block matrix only!
  int ibdiag = 0;
  int descra = 0;


  //BLAS and LAPACK vars
  char const *lower = "L";
  char const *upper = "U";
  char const *ntran = "N";
  char const *ytran = "T";
  char const *rside = "R";
  char const *lside = "L";
  double one = 1.0;
  double negOne = -1.0;
  double zero = 0.0;
  int incOne = 1;
  
  double *A = tmp_mm2;

  //Get A
  F77_NAME(dcopy)(&mm, V, &incOne, A, &incOne);
  F77_NAME(dpotrf)(lower, &m, A, &m, &info); 
  if(info != 0){error("mvMPPCovInvDet: Cholesky failed (1)\n");}

  //clear upper tri
  for(i = 0; i < m-1; i++){
    for(j = i+1; j < m; j++){
      A[j*m+i] = 0.0;
    }
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
  
  if(!brute){
 
    //////////////////
    //SWM
    ///////////////////
    //log(|C_str|) and tmp_qmqm = C_str^{-1|
    logDetC_str = 0.0;
    F77_NAME(dcopy)(&qmqm, C_str, &incOne, tmp_qmqm, &incOne);
    F77_NAME(dpotrf)(lower, &qm, tmp_qmqm, &qm, &info); if(info != 0){error("mvMPPCovInvDet: Cholesky failed (2)\n");}
    for(i = 0; i < qm; i++) logDetC_str += 2.0*log(tmp_qmqm[i*qm+i]);
    F77_NAME(dpotri)(lower, &qm, tmp_qmqm, &qm, &info); if(info != 0){error("mvMPPCovInvDet: Cholesky failed (3)\n");}
    
    //note, this is a total waist since I'm only going to take the block diag, but for now...
    //C = ct C_str^{-1} t(ct)
    F77_NAME(dsymm)(rside, lower, &nm, &qm, &one, tmp_qmqm, &qm, ct, &nm, &zero, tmp_nmqm, &nm);
    F77_NAME(dgemm)(ntran, ytran, &nm, &nm, &qm, &one, tmp_nmqm, &nm, ct, &nm, &zero, C, &nm);
    
    //make E = (Psi + V - blk[C])^{-1}
    logDetE = 0.0;
    for(i = 0, j = 0; i < n; i++){
      
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  tmp_mm[l*m+k] = Psi[l*m+k]+V[l*m+k]-C[(i*m+l)*nm+(i*m+k)];
	}
      }
      
      //     F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){cout << "c++ error: Cholesky failed - tmp_mm\n" << endl;}
      //     for(j = 0; j < m; j++) logDetE += 2.0*log(tmp_mm[j*m+j]);
      //     F77_NAME(dpotri)(lower, &m, tmp_mm, &m, &info); if(info != 0){cout << "c++ error: Cholesky inverse failed - tmp_mm\n" << endl;}
      
      //must must be a full upper/lower for later.
      logDetE += mtrxInvLogDet(tmp_mm, m, info);
      
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  E[j] = tmp_mm[l*m+k];
	  j++;
	}
      }
      
    }
    
    //make tmp_nmqm = E ct
    //F77_NAME(dsymm)(lside, lower, &nm, &qm, &one, E, &nm, ct, &nm, &zero, tmp_nmqm, &nm);
    F77_NAME(dbdimm)(&ntranInt, &mb, &qm, &kb, &one, &descra,
		     E, &blda, &ibdiag, &nbdiag, &lb,
		     ct, &nm, &zero, tmp_nmqm, &nm, work, &lwork);
    
    //make tmp_qmqm = C_str + t(ct) E ct
    F77_NAME(dgemm)(ytran, ntran, &qm, &qm, &nm, &one, ct, &nm, tmp_nmqm, &nm, &zero, tmp_qmqm1, &qm);
    F77_NAME(daxpy)(&qmqm, &one, C_str, &incOne, tmp_qmqm1, &incOne);
    
    //get log(|tmp_qmqm|) and log(|C|) then tmp_qmqm^{-1}
    logDet = 0.0;
    F77_NAME(dpotrf)(lower, &qm, tmp_qmqm1, &qm, &info); if(info != 0){error("mvMPPCovInvDet: Cholesky failed (4)\n");}
    for(j = 0; j < qm; j++) logDet += 2.0*log(tmp_qmqm1[j*qm+j]);
    F77_NAME(dpotri)(lower, &qm, tmp_qmqm1, &qm, &info); if(info != 0){error("mvMPPCovInvDet: Cholesky failed (5)\n");}
    
    logDet = logDetE+logDet-logDetC_str;
    
    //finally C = E - E ct (C_str + t(ct) E ct)^{-1} t(ct) E
    F77_NAME(dsymm)(rside, lower, &nm, &qm, &one, tmp_qmqm1, &qm, tmp_nmqm, &nm, &zero, tmp_nmqm1, &nm);
    
    //F77_NAME(dgemm)(ntran, ytran, &nm, &nm, &qm, &one, tmp_nmqm1, &nm, tmp_nmqm, &nm, &zero, C, &nm);
    //for(i = 0; i < nmnm; i++) 
    // C[i] = E[i]-C[i];
    
    F77_NAME(dgemm)(ntran, ytran, &nm, &nm, &qm, &negOne, tmp_nmqm1, &nm, tmp_nmqm, &nm, &zero, C, &nm);
    
    for(i = 0, j = 0; i < n; i++){
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  C[(i*m+l)*nm+(i*m+k)] = E[j]+C[(i*m+l)*nm+(i*m+k)];
	  j++;
	}
      }
    }
    
  }else{
    
    //////////////////
    //Brute force
    ///////////////////
    //C_str^{-1}
    F77_NAME(dpotrf)(lower, &qm, C_str, &qm, &info); if(info != 0){error("mvMPPCovInvDet: Cholesky failed (6)\n");}
    F77_NAME(dpotri)(lower, &qm, C_str, &qm, &info); if(info != 0){error("mvMPPCovInvDet: Cholesky failed (7)\n");}
    
    //C = ct C_str^{-1} t(ct)
    F77_NAME(dsymm)(rside, lower, &nm, &qm, &one, C_str, &qm, ct, &nm, &zero, tmp_nmqm, &nm);
    F77_NAME(dgemm)(ntran, ytran, &nm, &nm, &qm, &one, tmp_nmqm, &nm, ct, &nm, &zero, C, &nm);
    
    //make E = (Psi + V - blk[C])
    for(i = 0; i < n; i++){
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  C[(i*m+l)*nm+(i*m+k)] = Psi[l*m+k]+V[l*m+k];
	}
      }
    }
    
    //C_str^{-1}
    logDet = 0.0;
    F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("mvMPPCovInvDet: Cholesky failed (8)\n");}
    for(i = 0; i < nm; i++) logDet += 2.0*log(C[i*nm+i]);
    F77_NAME(dpotri)(lower, &nm, C, &nm, &info); if(info != 0){error("mvMPPCovInvDet: Cholesky failed (9)\n");}
    
}
  
    return logDet;
}
