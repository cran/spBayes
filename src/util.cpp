// Andrew O. Finley
// Dept. of Forest Resources
// University of Minnesota
// afinley@stat.umn.edu 
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2004 Andrew O. Finley

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
using namespace std;

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include "util.h"

void updateThetaGibbs(double *x, double *y, int &xnrow, int &xncol, double *fixedEffectSamples, 
		 double *covInv, double *tmpXRowCol, double *tmpXColCol, double *tmpXRow, 
		 double *tmpXCol, double *tmpXCol1, string &thetaPrior, double *thetaPriorMu, double *thetaPriorV){
  int info;
  const int incOne = 1;
  const double one = 1.0;
  const double zero = 0.0;
  const char upper = 'U';
  const char lower = 'L';
  const char ntran = 'N';
  const char ytran = 'T';
  const char rside = 'R';
  const char lside = 'L';
  
  //B ~ N(chol2inv(chol(t(x)%*%inv(cov)%*%x))%*%t(x)%*%s%*%y, chol2inv(chol(t(x)%*%inv(cov)%*%x)))
  //assumed the upper chol2inv(chol(cov)) was sent in as covInv

  //(t(x)%*%s%*%x)^{-1}
  F77_NAME(dsymm)(&lside, &upper, &xnrow, &xncol, &one, covInv, &xnrow, x, &xnrow, &zero, tmpXRowCol, &xnrow);
  F77_NAME(dgemm)(&ytran, &ntran, &xncol, &xncol, &xnrow, &one, x, &xnrow, tmpXRowCol, &xnrow, &zero, tmpXColCol, &xncol);

  if(thetaPrior == "NORMAL"){
    F77_NAME(daxpy)(&xncol, &one, thetaPriorV, &incOne, tmpXColCol, &incOne);
  }

  F77_CALL(dpotrf)(&lower, &xncol, tmpXColCol, &xncol, &info); if(info != 0){error("Cholesky failed\n");}
  F77_CALL(dpotri)(&lower, &xncol, tmpXColCol, &xncol, &info); if(info != 0){error("Cholesky inverse failed\n");}

  //mvrnorm mean
  //chol2inv(chol(t(x)%*%s%*%x))%*%t(x)%*%s%*%y
  F77_NAME(dsymv)(&upper, &xnrow, &one, covInv, &xnrow, y, &incOne, &zero, tmpXRow, &incOne);
  F77_NAME(dgemv)(&ytran, &xnrow, &xncol, &one, x, &xnrow, tmpXRow, &incOne, &zero, tmpXCol, &incOne);

  if(thetaPrior == "NORMAL"){
    F77_NAME(dgemv)(&ntran, &xncol, &xncol, &one, thetaPriorV, &xncol, thetaPriorMu, &incOne, &one, tmpXCol, &incOne);
  }

  F77_NAME(dsymv)(&lower, &xncol, &one, tmpXColCol, &xncol, tmpXCol, &incOne, &zero, tmpXCol1, &incOne);

  //my mvrnorm wants a lower Cholesky so
  F77_CALL(dpotrf)(&lower, &xncol, tmpXColCol, &xncol, &info); if(info != 0){error("Cholesky failed third\n");}
  mvrnorm(fixedEffectSamples, tmpXCol1, tmpXColCol, xncol); 
}


void mvrnorm(double *des, double *mu, double *cholCov, int dim){
  
  int i;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;

  //make some std norm draws
  for(i = 0; i < dim; i++)
    des[i] = rnorm(0.0, 1.0);

  //mult this vector by the lower triangle of the cholCov
  F77_NAME(dtrmv)("L", "N", "N", &dim, cholCov, &dim, des, &inc);

  //add the mean to the result
  F77_NAME(daxpy)(&dim, &one, mu, &inc, des, &inc);

}

void mvrnorm(double *des, double *mu, double *cholCov, int dim, bool upper){
  
  int i;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  
  //make some std norm draws
  for(i = 0; i < dim; i++)
    des[i] = rnorm(0.0, 1.0);

  //mult this vector by the lower triangle of the cholCov
  if(upper)
    F77_NAME(dtrmv)("U", "T", "N", &dim, cholCov, &dim, des, &inc);
  else
    F77_NAME(dtrmv)("L", "N", "N", &dim, cholCov, &dim, des, &inc);

  //add the mean to the result
  F77_NAME(daxpy)(&dim, &one, mu, &inc, des, &inc);

}

void showMatrix(double *x, int xnrow, int xncol){
  int i,j;
  for(i = 0; i < xnrow; i++){
    for(j = 0; j < xncol; j++){
      cout << x[j*xnrow+i] << "\t";
    }
    cout << endl;
  }      
}

SEXP getListElement (SEXP list, char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  
  if(elmt == R_NilValue){
    Rprintf("\nlist element %s not found\n", str);
  }
  return elmt;
}

void myCholInv(double *a, double *p, int &n, double &logdet){

  int i,j,k;
  double sum;
  logdet = 0;

  for(i = 0; i < n; i++){
    for(j = i; j < n; j++){
      for(sum = a[(j*n)+i], k = i-1; k >= 0; k--)
	sum -= a[(k*n)+i]*a[(k*n)+j];
      if(i == j){
	if(sum <=0.0)
	  error("matrix not pd");
	p[i]=sqrt(sum);
      }else{
	a[(i*n)+j] = sum/p[i];
      }
    }
  }

  for(i = 0; i < n; i++)
    logdet += log(p[i]);
  logdet = 2*logdet;

  for(i = 0; i < n; i++){
    a[(i*n)+i] = 1.0/p[i];
    for(j = i+1; j <  n; j++){
      sum=0.0;
      for(k = i; k < j; k++)
	sum -= a[(k*n)+j]*a[(i*n)+k];
      a[(i*n)+j]=sum/p[j];
    }
  }
  
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      a[(j*n)+i] = 0.0;
    }
  }
  for (i = 0; i < n; i++) {
    a[(i*n)+i] *= a[(i*n)+i];
    for (k = i + 1; k < n; k++) {
      a[(i*n)+i] += a[(i*n)+k] * a[(i*n)+k];
    }
    for (j = i + 1; j < n; j++) {
      for (k = j; k < n; k++) {
	a[(j*n)+i] += a[(i*n)+k] * a[(j*n)+k];
      }
    }
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      a[(j*n)+i] = a[(i*n)+j];
    }
  }
}

void myChol(double *a, double *p, int &n){

  int i,j,k;
  double sum;

  for(i = 0; i < n; i++){
    for(j = i; j < n; j++){
      for(sum = a[(j*n)+i], k = i-1; k >= 0; k--)
	sum -= a[(k*n)+i]*a[(k*n)+j];
      if(i == j){
	if(sum <=0.0)
	  error("matrix not pd");
	p[i]=sqrt(sum);
      }else{
	a[(i*n)+j] = sum/p[i];
      }
    }
  }

  for(i = 0; i < n; i++)
    a[(i*n)+i] = p[i];

}
