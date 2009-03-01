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

double mvCovInvDet(double *coordsD, double *C, int n, int m, double *Psi, double *V, double *theta, 
		   double *tmp_mm, double *tmp_mm1, double *tmp_mm2, 
		   string covModel, int nPramPtr, covmodel *covModelObj, 
		   void (covmodel::*cov1ParamPtr)(double, double &, double &),
		   void (covmodel::*cov2ParamPtr)(double, double, double &, double&)){
  
  //n is the number of coords
  //m is the number of response variables
  int nm = n*m;
  int mm = m*m;
  int nn = n*n;
  int nmnm = nm*nm;

  int i, j, k, l, info;
  double logDet = 0.0;

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
  if(info != 0){error("mvCovInvDet: Cholesky failed (1)\n");}

  //clear upper tri
  for(i = 0; i < m-1; i++){
    for(j = i+1; j < m; j++){
      A[j*m+i] = 0.0;
    }
  }
  
  //
  //make C
  //
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      
      zeros(tmp_mm, mm);
      
      for(k = 0; k < m; k++){
	if(nPramPtr == 1)
	  (covModelObj->*cov1ParamPtr)(theta[k], tmp_mm[k*m+k], coordsD[j*n+i]);
	else //i.e., 2 parameter matern
	  (covModelObj->*cov2ParamPtr)(theta[k], theta[m+k], tmp_mm[k*m+k], coordsD[j*n+i]);
      }
      
      F77_NAME(dgemm)(ntran, ntran, &m, &m, &m, &one, A, &m, tmp_mm, &m, &zero, tmp_mm1, &m);
      F77_NAME(dgemm)(ntran, ytran, &m, &m, &m, &one, tmp_mm1, &m, A, &m, &zero, tmp_mm, &m);
      
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  C[((j*m+l)*nm)+(i*m+k)] = tmp_mm[l*m+k];
	  tmp_mm[l*m+k] = 0.0; //zero out
	}
      }
    }
  }
  
  for(i = 0; i < n; i++){
    for(k = 0; k < m; k++){
      for(l = 0; l < m; l++){
	C[(i*m+l)*nm+(i*m+k)] += Psi[l*m+k];
      }
    }
  }
  
  //C^{-1}
  logDet = 0.0;
  F77_NAME(dpotrf)(lower, &nm, C, &nm, &info); if(info != 0){error("mvCovInvDet: Cholesky failed (8)\n");}
  for(i = 0; i < nm; i++) logDet += 2.0*log(C[i*nm+i]);
  F77_NAME(dpotri)(lower, &nm, C, &nm, &info); if(info != 0){error("mvCovInvDet: Cholesky failed (9)\n");}
  
  return logDet;
}
