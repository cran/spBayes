// Andrew O. Finley
// Dept. of Forest Resources
// University of Minnesota
// finleya@msu.edu 
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2004 Andrew O. Finley

#include <iostream>
#include <iomanip>
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

void updateThetaGibbs(double *x, double *y, int &nxrow, int &nxcol, double *fixedEffectSamples, 
		      double *cov, double *tmpXRowCol, double *tmpXCowCol, double *tmpXRow, 
		      double *tmpXCol, double *tmpXCol1, string &thetaPrior, double *thetaPriorMu, double *thetaPriorV);

void mvrnorm(double *des, double *mu, double * cholCov, int dim);
void mvrnorm(double *des, double *mu, double * cholCov, int dim, bool upper);

void showMatrix(double *x, int xnrow, int xncol);
void writeRMatrix(string outfile, double * a, int nrow, int ncol);

SEXP getListElement (SEXP list, const char *str);

void zeros(double *x, int length);

void identity(double *x, int &nrow);

void kron(double *a, int &dima1, int &dima2, 
	  double *b, int &dimb1, int &dimb2, 
	  double *c, int &dimc1, int &dimc2);

void setLowerChol(double *A, double *S, int dim);

string toString(int &x);

double dTNorm(double x, double mu, double sd, double a, double b);

void diagmm(int &nrow_b, int &ncol_b, double *a, double *b, double *c);


void subsetCovRow(double *x, int n, int p, int begin, int end, double *cov, double *means);
void subsetCovCol(double *x, int p, int begin, int end, double *cov, double *means);

double mtrxInvLogDet(double *m, int dim, int info);

void mtrxInv(double *m, int dim, int info);


double logit(double theta, double a, double b);

double logitInv(double z, double a, double b);

void covTransInv(double *z, double *v, int m);

void covTrans(double *v, double *z, int m);

void covTransInvExpand(double *v, double *z, int m);

void zeroUpperTri(double *v, int m);

void printMtrx(double *m, int nRow, int nCol);

void printVec(double *m, int n);

double logit_logpost(int &n, double *Y, double *eta, double *w);

double poisson_logpost(int &n, double *Y, double *eta, double *w);

void report(int &s, int &nSamples, int &status, int &nReport, bool &verbose);
