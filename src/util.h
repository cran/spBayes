#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

void mvrnorm(double *des, double *mu, double * cholCov, int dim);

void mvrnorm(double *des, double *mu, double * cholCov, int dim, bool upper);

SEXP getList(SEXP list, const char *str);

SEXP getGetList(SEXP list, const char *str1, const char *str2);

void zeros(double *x, int length);

void zeros(int *x, int length);

void iden(double *x, int &nrow);

void kron(double *a, int &dima1, int &dima2, 
	  double *b, int &dimb1, int &dimb2, 
	  double *c, int &dimc1, int &dimc2);

void setLowerChol(double *A, double *S, int dim);

double dTNorm(double x, double mu, double sd, double a, double b);

void diagmm(int &nrow_b, int &ncol_b, double *a, double *b, double *c);

void subsetCovRow(double *x, int n, int p, int begin, int end, double *cov, double *means);

void subsetCovCol(double *x, int p, int begin, int end, double *cov, double *means);

double logit(double theta, double a, double b);

double logitInv(double z, double a, double b);

void covTransInv(double *z, double *v, int m);

void covTrans(double *v, double *z, int m);

void covTransInvExpand(double *v, double *z, int m);

void covExpand(double *v, double *z, int m);

void printMtrx(double *m, int nRow, int nCol);

void printVec(double *m, int n);

void printVec(int *m, int n);

double logit_logpost(int &n, double *Y, double *eta, double *w);

double binomial_logpost(int &n, double *Y, double *eta, double *w, int *r);

double poisson_logpost(int &n, double *Y, double *eta, double *w, int *r);

double binomial_logpost(int &n, double *Y, double *eta, int *r);

double poisson_logpost(int &n, double *Y, double *eta, int *r);

void report(int &s, int &nSamples, int &status, int &nReport, bool &verbose);

void spCor(double *D, int n, double *theta, std::string &covModel, double *C);

double spCor(double D, double *theta, std::string &covModel);

double spCor(double D, double phi, double nu, std::string &covModel);

void spCov(double *D, int n, double *theta, std::string &covModel, double *C); 

void spCovLT(double *D, int n, double *theta, std::string &covModel, double *C); 

void transpose(double *m, int w, int h);

void clearUT(double *m, int n);

void rwish(double *S, int v, int p, double *Z, double *tmp_pp, int iwish);
