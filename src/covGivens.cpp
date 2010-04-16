#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"


extern "C" {
  
  SEXP covGivens(SEXP p_R, SEXP theta, SEXP lambda, SEXP cov, SEXP tmp1, SEXP tmp2){
    int i,j,k,p,pp;
    char ntran = 'N';
    char ytran = 'T';
    double zero = 0.0;
    double one = 1.0;
    int inc = 1;
    p = INTEGER(p_R)[0];
    pp = p*p;
    iden(REAL(tmp1), p);
    iden(REAL(tmp2), p);
    iden(REAL(cov), p);
    

    //D = P \Lambda P^T
     for(i = 0, k=0; i < p-1; i++){
      for(j = i+1; j < p; j++, k++){

	iden(REAL(cov), p);
 	REAL(cov)[i*p+i] = REAL(cov)[j*p+j] = cos(REAL(theta)[k]);	
	REAL(cov)[j*p+i] = REAL(cov)[i*p+j] = sin(REAL(theta)[k]);
	REAL(cov)[i*p+j] *= -1.0;

	F77_NAME(dgemm)(&ntran, &ntran, &p, &p, &p, &one, REAL(tmp2), &p, REAL(cov), &p, &zero, REAL(tmp1), &p);
	F77_NAME(dcopy)(&pp, REAL(tmp1), &inc, REAL(tmp2), &inc);
   
      }
     }
     
     zeros(REAL(cov), pp);
     for(i = 0; i < p; i++) REAL(cov)[i*p+i] = REAL(lambda)[i];
     F77_NAME(dgemm)(&ntran, &ntran, &p, &p, &p, &one, REAL(tmp2), &p, REAL(cov), &p, &zero, REAL(tmp1), &p);
     F77_NAME(dgemm)(&ntran, &ytran, &p, &p, &p, &one, REAL(tmp1), &p, REAL(tmp2), &p, &zero, REAL(cov), &p);
     
     return(R_NilValue);
  }
}


      

