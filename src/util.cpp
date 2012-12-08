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

// void showMatrix(double *x, int xnrow, int xncol){
//   int i,j;
//   for(i = 0; i < xnrow; i++){
//     for(j = 0; j < xncol; j++){
//       cout << x[j*xnrow+i] << "\t";
//     }
//     cout << endl;
//   }      
// }



// void writeRMatrix(string outfile, double * a, int nrow, int ncol){
//     ofstream file(outfile.c_str());
//     if ( !file ) {
//       cout << "Data file could not be opened." << endl;
//       //exit(1);
//     }
  

//   for(int i = 0; i < nrow; i++){
//     for(int j = 0; j < ncol-1; j++){
//       file << setprecision(20) << fixed << a[j*nrow+i] << "\t";
//     }
//     file << setprecision(20) << fixed << a[(ncol-1)*nrow+i] << endl;    

//   }
//   file.close();
// }

SEXP getList(SEXP list, const char *str){
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

SEXP getGetList(SEXP list, const char *str1, const char *str2){
  SEXP list2 = getList(list, str1);
  return getList(list2, str2);  
}


void zeros(double *x, int length){
  for(int i = 0; i < length; i++)
    x[i] = 0.0;
}

void zeros(int *x, int length){
  for(int i = 0; i < length; i++)
    x[i] = 0;
}

void iden(double *x, int &nrow){

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

void diagmm(int &nrow_b, int &ncol_b, double *a, double *b, double *c){
  for(int i = 0; i < nrow_b; i++){
    for(int j = 0; j < ncol_b; j++){
      c[j*nrow_b+i] = b[j*nrow_b+i] * a[i];
    }
  }
}



void subsetCovRow(double *x, int n, int p, int begin, int end, double *cov, double *means){
  
  int nSubset = end-begin+1;
  int i,j,k;
  
  for(i = 0; i < p; i++){
    means[i] = 0.0;
    for(j = 0; j < p; j++){
      cov[j*p+i] = 0.0;
    }
  }
  
  
  for(i = 0; i < p; i++){
    for(j = 0; j < nSubset; j++){
      means[i] += x[(i*n)+(begin+j)];
    }
    means[i] = means[i]/nSubset;
  }
  
  
  for(i = 0; i < p; i++){
    for(j = i; j < p; j++){
      for(k = 0; k < nSubset; k++){
	cov[i*p+j] += (x[(i*n)+(begin+k)]-means[i])*(x[(j*n)+(begin+k)]-means[j]);
      }
      cov[i*p+j] = cov[i*p+j]/(nSubset-1);
    }
  }
   
}


void subsetCovCol(double *x, int p, int begin, int end, double *cov, double *means){
  
  int nSubset = end-begin+1;
  int i,j,k;
  
  for(i = 0; i < p; i++){
    means[i] = 0.0;
    for(j = 0; j < p; j++){
      cov[j*p+i] = 0.0;
    }
  }
  
  for(i = 0; i < p; i++){
    for(j = 0; j < nSubset; j++){
      means[i] += x[(begin+j)*p+i];
    }
    means[i] = means[i]/nSubset;
  }
  
  
  for(i = 0; i < p; i++){
    for(j = i; j < p; j++){
      for(k = 0; k < nSubset; k++){
	cov[i*p+j] += (x[(begin+k)*p+i]-means[i])*(x[(begin+k)*p+j]-means[j]);
      }
      cov[i*p+j] = cov[i*p+j]/(nSubset-1);
    }
  }

}

  

// double mtrxInvLogDet(double *m, int dim, int info){

//   double logDet = 0.0;
//   int i,j;

//   F77_NAME(dpotrf)("L", &dim, m, &dim, &info); if(info != 0){cout << "c++ error: mtrxInvLogDet Cholesky failed\n" << endl;}
//   for(j = 0; j < dim; j++) logDet += 2.0*log(m[j*dim+j]);
//   F77_NAME(dpotri)("L", &dim, m, &dim, &info); if(info != 0){cout << "c++ error: mtrxInvLogDet Cholesky inverse failed\n" << endl;}
  
//   for(i = 1; i < dim; i++){
//     for(j = 0; j < i; j++){
//       m[i*dim+j] = m[j*dim+i];
//     }
//   }

//   return(logDet);
// }


// void mtrxInv(double *m, int dim, int info){

//   int i,j;

//   F77_NAME(dpotrf)("L", &dim, m, &dim, &info); if(info != 0){cout << "c++ error: mtrxInv Cholesky failed\n" << endl;}
//   F77_NAME(dpotri)("L", &dim, m, &dim, &info); if(info != 0){cout << "c++ error: mtrxInv Cholesky inverse failed\n" << endl;}
  
//   for(i = 1; i < dim; i++){
//     for(j = 0; j < i; j++){
//       m[i*dim+j] = m[j*dim+i];
//     }
//   }

// }

double logit(double theta, double a, double b){
  return log((theta-a)/(b-theta));
}

double logitInv(double z, double a, double b){
  return b-(b-a)/(1+exp(z));
}


void covTransInv(double *z, double *v, int m){
  int i, j, k;

  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      v[k] = z[k];
      if(i == j)
	v[k] = exp(z[k]);
    }
  }

}

void covTrans(double *v, double *z, int m){
  int i, j, k;

  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      z[k] = v[k];
      if(i == j)
	z[k] = log(v[k]);
    }
  }

}

void covTransInvExpand(double *v, double *z, int m){
  int i, j, k;
  
  zeros(z, m*m);
  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      z[i*m+j] = v[k];
      if(i == j)
	z[i*m+j] = exp(z[i*m+j]);
    }
  }
  
}

void covExpand(double *v, double *z, int m){
  int i, j, k;
  
  zeros(z, m*m);
  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      z[i*m+j] = v[k];
    }
  }
  
}

void printMtrx(double *m, int nRow, int nCol){

  int i, j;

  for(i = 0; i < nRow; i++){
    Rprintf("\t");
    for(j = 0; j < nCol; j++){
      Rprintf("%.3f\t", m[j*nRow+i]);
    }
    Rprintf("\n");
  }
}

void printVec(double *m, int n){

  Rprintf("\t");
    for(int j = 0; j < n; j++){
      Rprintf("%.3f\t", m[j]);
    }
    Rprintf("\n");
}



double logit_logpost(int &n, double *Y, double *eta, double *w){
  double loglike = 0.0;
  int i;  

  for(i = 0; i < n; i++)
    loglike += Y[i]*(eta[i]+w[i]);
  
  for(i = 0; i < n; i++)
    loglike -= log(1.0+exp(eta[i]+w[i]));

  return loglike;
}

double binomial_logpost(int &n, double *Y, double *eta, double *w, int *weights){
  double p, loglike = 0.0;
  int i;  

  for(i = 0; i < n; i++){
    p = 1.0/(1+exp(-eta[i]-w[i]));
    loglike += Y[i]*log(p)+(static_cast<double>(weights[i])-Y[i])*log(1.0-p);
  }

  return loglike;
}

double poisson_logpost(int &n, double *Y, double *eta, double *w, int *r){
  double loglike = 0.0;
  int i;  

  for(i = 0; i < n; i++){
    loglike += Y[i]*(eta[i]+w[i]+log(r[i]))-exp(eta[i]+w[i]+log(r[i]));
  }

  return loglike;
}

double binomial_logpost(int &n, double *Y, double *eta, int *weights){
  double p, loglike = 0.0;
  int i;  

  for(i = 0; i < n; i++){
    p = 1.0/(1+exp(-eta[i]));
    loglike += Y[i]*log(p)+(static_cast<double>(weights[i])-Y[i])*log(1.0-p);
  }

  return loglike;
}

double poisson_logpost(int &n, double *Y, double *eta, int *r){
  double loglike = 0.0;
  int i;  

  for(i = 0; i < n; i++){
    loglike += Y[i]*(eta[i]+log(r[i]))-exp(eta[i]+log(r[i]));
  }

  return loglike;
}

void report(int &s, int &nSamples, int &status, int &nReport, bool &verbose){

  if(verbose){
    if(status == nReport){
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
      //Rprintf("---------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
      status = 0;
    }
  }
  status++;
  
  R_CheckUserInterrupt();  
}

// void dnscsr(double *A, int nrowA, int ncolA, double *a, int *indxa, int *pntrb, int *pntre){
//   //untested!!!
//   //int nzmaxa = nrowA*ncolA;
//   //int *indxa = (int *) R_alloc(nzmaxa, sizeof(int));
//   //int *pntrb = (int *) R_alloc(nrowA, sizeof(int)); 
//   //int *pntre = (int *) R_alloc(nrowA, sizeof(int)); 
//   //double *a = (double *) R_alloc(nzmaxa, sizeof(double));
  
//   //Convert dense matrix to NIST one-based indexing (following intel MKL docs)
//   //todo: add eps if needed
  
//   int i, j, k, l, ll;
//   for(i = 0; i < nrowA; i++){
//     pntrb[i]=pntre[i]=indxa[i]=0;
//   }
  
//   bool first;
//   for(i = 0, l = 0, ll = 0; i < nrowA; i++, ll = 0){
//     first = true;
//     for(j = 0; j < ncolA; j++){
//       if(A[j*nrowA+i] != 0.0){
// 	a[k] = A[j*nrowA+i];
// 	indxa[k] = j+1;
// 	if(first){
// 	  pntrb[l] = k+1;
// 	  first = false;
// 	  l++;
// 	}
// 	ll++;
// 	k++;
//       }
//     }
//     if(!first){
//       pntre[l-1]=pntrb[l-1]+ll;
//     }else{
//       pntrb[l] = 0;
//       pntre[l] = 0;
//       l++;
//     }
//   }
// }

void spCor(double *D, int n, double *theta, string &covModel, double *C){
  int i;
  
  if(covModel == "exponential"){
    
    for(i = 0; i < n; i++){
      C[i] = exp(-1.0*theta[0]*D[i]);
    }
    
  }else if(covModel == "spherical"){
    
    for(i = 0; i < n; i++){
      if(D[i] > 0 && D[i] <= 1.0/theta[0]){
	C[i] = 1.0 - 1.5*theta[0]*D[i] + 0.5*pow(theta[0]*D[i],3);
      }else if(D[i] >= 1.0/theta[0]){
	C[i] = 0.0;
      }else{
	C[i] = 1.0;
      }
    }
    
  }else if(covModel == "gaussian"){
    
    for(i = 0; i < n; i++){
      C[i] = exp(-1.0*(pow(theta[0]*D[i],2)));
    }
    
  }else if(covModel == "matern"){
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    for(i = 0; i < n; i++){
      if(D[i]*theta[0] > 0.0){
	C[i] = pow(D[i]*theta[0], theta[1])/(pow(2, theta[1]-1)*gammafn(theta[1]))*bessel_k(D[i]*theta[0], theta[1], 1.0);
      }else{
	C[i] = 1.0;
      }
    }
    
 }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

double spCor(double D, double *theta, string &covModel){
  
  if(covModel == "exponential"){
    
    return exp(-1.0*theta[0]*D);
    
  }else if(covModel == "spherical"){
    
    if(D > 0 && D <= 1.0/theta[0]){
      return 1.0 - 1.5*theta[0]*D + 0.5*pow(theta[0]*D,3);
    }else if(D >= 1.0/theta[0]){
      return 0.0;
    }else{
      return 1.0;
    }
    
  }else if(covModel == "gaussian"){
    
    return exp(-1.0*(pow(theta[0]*D,2)));
    
  }else if(covModel == "matern"){
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    if(D*theta[0] > 0.0){
      return pow(D*theta[0], theta[1])/(pow(2, theta[1]-1)*gammafn(theta[1]))*bessel_k(D*theta[0], theta[1], 1.0);
    }else{
      return 1.0;
    }    
    
  }else{
    error("c++ error: cov.model is not correctly specified");
    return 0;
  }
}

double spCor(double D, double phi, double nu, string &covModel){
  
  if(covModel == "exponential"){
    
    return exp(-1.0*phi*D);
    
  }else if(covModel == "spherical"){
    
    if(D > 0 && D <= 1.0/phi){
      return 1.0 - 1.5*phi*D + 0.5*pow(phi*D,3);
    }else if(D >= 1.0/phi){
      return 0.0;
    }else{
      return 1.0;
    }
  }else if(covModel == "matern"){
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    if(D*phi > 0.0){
      return pow(D*phi, nu)/(pow(2, nu-1)*gammafn(nu))*bessel_k(D*phi, nu, 1.0);
    }else{
      return 1.0;
    } 
  }else if(covModel == "gaussian"){
    
    return exp(-1.0*(pow(phi*D,2)));
      
  }else{
    error("c++ error: cov.model is not correctly specified");
    return 0;
  }
}



void spCov(double *D, int n, double *theta, string &covModel, double *C){
  int i;
  
  if(covModel == "exponential"){
    
    for(i = 0; i < n; i++){
      C[i] = theta[0]*exp(-1.0*theta[1]*D[i]);
    }
    
  }else if(covModel == "spherical"){
    
    for(i = 0; i < n; i++){
      if(D[i] > 0 && D[i] <= 1.0/theta[1]){
	C[i] = theta[0]*(1.0 - 1.5*theta[1]*D[i] + 0.5*pow(theta[1]*D[i],3));
      }else if(D[i] >= 1.0/theta[1]){
	C[i] = 0.0;
      }else{
	C[i] = theta[0];
      }
    }
    
  }else if(covModel == "gaussian"){
    
    for(i = 0; i < n; i++){
      C[i] = theta[0]*exp(-1.0*(pow(theta[1]*D[i],2)));
    }
    
  }else if(covModel == "matern"){
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    for(i = 0; i < n; i++){
      if(D[i]*theta[1] > 0.0){
	C[i] = theta[0]*pow(D[i]*theta[1], theta[2])/(pow(2, theta[2]-1)*gammafn(theta[2]))*bessel_k(D[i]*theta[1], theta[2], 1.0);
      }else{
	C[i] = theta[0];
      }
    }
    
 }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

void spCovLT(double *D, int n, double *theta, string &covModel, double *C){
  int i,j;
  
  if(covModel == "exponential"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	C[i*n+j] = theta[0]*exp(-1.0*theta[1]*D[i*n+j]);
      }
    }
    
  }else if(covModel == "spherical"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	if(D[i*n+j] > 0 && D[i*n+j] <= 1.0/theta[1]){
	  C[i*n+j] = theta[0]*(1.0 - 1.5*theta[1]*D[i*n+j] + 0.5*pow(theta[1]*D[i*n+j],3));
	}else if(D[i*n+j] >= 1.0/theta[1]){
	  C[i*n+j] = 0.0;
	}else{
	  C[i*n+j] = theta[0];
	}
      }
    }
    
  }else if(covModel == "gaussian"){
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	C[i*n+j] = theta[0]*exp(-1.0*(pow(theta[1]*D[i*n+j],2)));
      }
    }

  }else if(covModel == "matern"){
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    for(i = 0; i < n; i++){
      for(j = i; j < n; j++){ 
	if(D[i*n+j]*theta[1] > 0.0){
	  C[i*n+j] = theta[0]*pow(D[i*n+j]*theta[1], theta[2])/(pow(2, theta[2]-1)*gammafn(theta[2]))*bessel_k(D[i*n+j]*theta[1], theta[2], 1.0);
	}else{
	  C[i*n+j] = theta[0];
	}
      }
    }
    
 }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

void transpose(double *m, int w, int h){
  int start, next, i;
  double tmp;
  
  for(start = 0; start <= w * h - 1; start++) {
    next = start;
    i = 0;
    do{	i++;
      next = (next % h) * w + next / h;
    }while (next > start);
    if(next < start || i == 1) continue;
    
    tmp = m[next = start];
    do{
      i = (next % h) * w + next / h;
      m[next] = (i == start) ? tmp : m[i];
      next = i;
    }while (next > start);
  }
}

void clearUT(double *m, int n){
  for(int i = 1; i < n; i++){
    for(int j = 0; j < i; j++){
      m[i*n+j] = 0;
    }
  }
}
