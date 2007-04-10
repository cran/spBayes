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

SEXP getListElement (SEXP list, char *str);

void zeros(double *x, int length);

void identity(double *x, int &nrow);

void kron(double *a, int &dima1, int &dima2, 
	  double *b, int &dimb1, int &dimb2, 
	  double *c, int &dimc1, int &dimc2);

void setLowerChol(double *A, double *S, int dim);

string toString(int &x);

double dTNorm(double x, double mu, double sd, double a, double b);
