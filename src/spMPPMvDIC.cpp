#include <string>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

extern"C" {

  SEXP spMPPMvDIC(SEXP Q_r, SEXP knotsD_r, SEXP coordsKnotsD_r, SEXP n_r, SEXP m_r, SEXP q_r, 
		  SEXP Psi_r, SEXP V_r, SEXP phi_r, SEXP nu_r, SEXP covModel_r, SEXP CEps_r){
    
  int n = INTEGER(n_r)[0];
  int m = INTEGER(m_r)[0];
  int q = INTEGER(q_r)[0];

  int nm = n*m;
  int mm = m*m;
  int nn = n*n;
  int qm = q*m;
  int nmnm = nm*nm;
  int qmqm = qm*qm;
  int nmqm = nm*qm;
  int h, i, j, k, l, ii, jj, info;
 
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

  double *ct = (double *) R_alloc(nmqm, sizeof(double));
  double *C_str = (double *) R_alloc(qmqm, sizeof(double));
  double *C = (double *) R_alloc(nmnm, sizeof(double));
  double *tmp_nmqm = (double *) R_alloc(nmqm, sizeof(double));
  double *tmp_mm = (double *) R_alloc(mm, sizeof(double));
  double *tmp_mm1 = (double *) R_alloc(mm, sizeof(double));
  double *tmp_nm = (double *) R_alloc(nm, sizeof(double));
  double *A = (double *) R_alloc(mm, sizeof(double));
  double *theta = (double *) R_alloc(2, sizeof(double));
 
  double logDet = 0.0;

  std::string covModel = CHAR(STRING_ELT(covModel_r,0));
  
  //Get A
  F77_NAME(dcopy)(&mm, REAL(V_r), &incOne, A, &incOne);
  F77_NAME(dpotrf)(lower, &m, A, &m, &info); if(info != 0){error("c++ error: dpotrf failed 1\n");}

  SEXP D_r;
  PROTECT(D_r = allocVector(REALSXP, 1)); 
  REAL(D_r)[0] = 0.0;

  //clear upper tri
  clearUT(A, m);
  
  //make ct
  for(jj = 0; jj < q; jj++){
    for(ii = 0; ii < n; ii++){	
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  ct[(k+jj*m)*nm+(ii*m+l)] = 0.0; 
	  for(h = 0; h < m; h++){
	    theta[0] = REAL(phi_r)[h];
	    if(covModel == "matern"){
	      theta[1] = REAL(nu_r)[h];
	    }
	    ct[(k+jj*m)*nm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(REAL(coordsKnotsD_r)[jj*n+ii], theta, covModel);
	  }
	}
      }
    }
  }

  //make C_str
  for(jj = 0; jj < q; jj++){
    for(ii = jj; ii < q; ii++){	
      for(k = 0; k < m; k++){
	for(l = 0; l < m; l++){
	  C_str[(k+jj*m)*qm+(ii*m+l)] = 0.0; 
	  for(h = 0; h < m; h++){
	    theta[0] = REAL(phi_r)[h];
	    if(covModel == "matern"){
	      theta[1] = REAL(nu_r)[h];
	    }
	    C_str[(k+jj*m)*qm+(ii*m+l)] += A[k+m*h]*A[l+m*h]*spCor(REAL(knotsD_r)[jj*q+ii], theta, covModel);
	  }
	}
      }
    }
  }

  F77_NAME(dpotrf)(lower, &qm, C_str, &qm, &info); if(info != 0){error("c++ error: dpotrf failed 2\n");}
  F77_NAME(dpotri)(lower, &qm, C_str, &qm, &info); if(info != 0){error("c++ error: dpotri failed 3\n");}
  
  //C = ct C_str^{-1} t(ct)
  F77_NAME(dsymm)(rside, lower, &nm, &qm, &one, C_str, &qm, ct, &nm, &zero, tmp_nmqm, &nm);
  F77_NAME(dgemm)(ntran, ytran, &nm, &nm, &qm, &one, tmp_nmqm, &nm, ct, &nm, &zero, C, &nm);
    
  for(i = 0; i < n; i++){
    
    for(k = 0; k < m; k++){
      for(l = 0; l < m; l++){
	tmp_mm[l*m+k] = REAL(Psi_r)[l*m+k]+REAL(V_r)[l*m+k]-C[(i*m+l)*nm+(i*m+k)];
      }
    }
    
    F77_NAME(dcopy)(&mm, tmp_mm, &incOne, &REAL(CEps_r)[i*mm], &incOne);

    F77_NAME(dpotrf)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotrf failed 4\n");}
    for(j = 0; j < m; j++) logDet += 2.0*log(tmp_mm[j*m+j]);
    F77_NAME(dpotri)(lower, &m, tmp_mm, &m, &info); if(info != 0){error("c++ error: dpotri failed 5\n");}
     
    F77_NAME(dsymv)(lower, &m, &one, tmp_mm, &m, &REAL(Q_r)[i*m], &incOne, &zero, &tmp_nm[i*m], &incOne);
  }
  
  REAL(D_r)[0] = logDet+F77_NAME(ddot)(&nm, tmp_nm, &incOne, REAL(Q_r), &incOne);
  
  UNPROTECT(1);
  
  return(D_r);
  }

}
