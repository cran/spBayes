#include <iostream>
#include <string>
using namespace std;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern "C" {
  
  SEXP spMvDIC(SEXP n_r, SEXP m_r, SEXP Q_r, SEXP Psi_r, SEXP tmp_mm_r, SEXP tmp_nm_r){

    //BLAS and LAPACK vars
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

    double logDet = 0.0;
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    int nm = n*m;
    int mm = m*m;
    int i, j, k, l, info=0;

    SEXP D_r;
    PROTECT(D_r = allocVector(REALSXP, 1)); 
    REAL(D_r)[0] = 0.0;
    
    F77_NAME(dcopy)(&mm, REAL(Psi_r), &incOne, REAL(tmp_mm_r), &incOne);
    
    logDet = n*mtrxInvLogDet(REAL(tmp_mm_r), m, info);
    
    for(i = 0; i < n; i++){
      F77_NAME(dsymv)(lower, &m, &one, REAL(tmp_mm_r), &m, &REAL(Q_r)[i*m], &incOne, &zero, &REAL(tmp_nm_r)[i*m], &incOne);
    }
    
    REAL(D_r)[0] = logDet+F77_NAME(ddot)(&nm, REAL(tmp_nm_r), &incOne, REAL(Q_r), &incOne);
    
    
    UNPROTECT(1);

    return(D_r);
  }
  
}
