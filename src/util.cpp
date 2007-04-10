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
#include <fstream>
#include <string>
#include <sstream>
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


void writeRMatrix(string outfile, double * a, int nrow, int ncol){
    ofstream file(outfile.c_str());
    if ( !file ) {
      cerr << "Data file could not be opened." << endl;
      exit(1);
    }
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol-1; j++){
      file << fixed << a[j*nrow+i] << " ";
    }
    file << fixed << a[(ncol-1)*nrow+i] << endl;    

  }
  file.close();
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


void zeros(double *x, int length){
  for(int i = 0; i < length; i++)
    x[i] = 0.0;
}

void identity(double *x, int &nrow){

  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < nrow; j++){
      if(i != j)
	x[j*nrow+i] = 0.0;
      else
	x[j*nrow+i] = 1.0;
    }
  }

}

void kron(double *a, int &dima1, int &dima2, 
	  double *b, int &dimb1, int &dimb2, 
	  double *c, int &dimc1, int &dimc2){
  
  int i, j, k, l;
  
  for (k = 0; k < dima1; k++) {
    for (l = 0; l < dima2; l++) {
      for (i = 0; i < dimb1; i++) {
	for (j = 0; j < dimb2; j++) {
	  c[(l*dimb2+j)*dimc1+(k*dimb1+i)] = a[l*dima1+k] * b[j*dimb1+i];
	}
      }
    }
  }
}

void setLowerChol(double *A, double *S, int dim){
  int i, j, k;
  
  zeros(A, dim*dim);
  for(i = 0, k = 0; i < dim; i++){
    for(j = i; j < dim; j++, k++){
      A[i*dim+j] = S[k];
    }
  }
}


string toString(int &x) {
  ostringstream oss;
  oss << x;
  return oss.str();
}
 

double dTNorm(double x, double mu, double sd, double a, double b){
  if(x < a || x > b)
    return 0.0;
  else
    return dnorm(x, mu, sd, false)/(pnorm(b, mu, sd, true, false) - pnorm(a, mu, sd, true, false));
}
