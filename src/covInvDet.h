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

/* knotsD       mxm knot distance matrix */
/* coordsKnotsD nxm coordinate distance matrix */
/* C            nxn empty */
/* C_str        mxm empty */
/* ct           nxm empty */
/* n            # of coordinates */
/* m            # of knots */
/* tau_sq       residual variance*/
/* sigma.sq     spatial variance */
/* theta        1 if nPramPtr is 1 or 2 if nPramPtr is 2 these are the spatial range parameters */
/* ... */

double pPCovInvDet(double *knotsD, double *coordsKnotsD, double *C, double *C_str, double *ct, 
		   int n, int m, double tauSq, double sigmaSq, double* theta, 
		   double *tmp_mm, double *tmp_nm,
		   string covModel, int nPramPtr, covmodel *covModelObj, 
		   void (covmodel::*cov1ParamPtr)(double, double &, double &),
		   void (covmodel::*cov2ParamPtr)(double, double, double &, double&), int brute);

/* knotsD       mxm knot distance matrix */
/* coordsKnotsD nxm coordinate distance matrix */
/* C            nxn empty */
/* C_str        mxm empty */
/* ct           nxm empty */
/* E            n empty */
/* n            # of coordinates */
/* m            # of knots */
/* tau_sq       residual variance*/
/* sigma.sq     spatial variance */
/* theta        1 if nPramPtr is 1 or 2 if nPramPtr is 2 these are the spatial range parameters */
/* ... */

/*note: tau_sq > 0*/
double mPPCovInvDet(double *knotsD, double *coordsKnotsD, double *C, double *C_str, double *ct, 
		    double *E, double *Einv, int n, int m, double tauSq, double sigmaSq, double* theta, 
		    double *tmp_mm, double *tmp_nm, double *tmp_nm1, double *tmp_nn,
		    string covModel, int nPramPtr, covmodel *covModelObj, 
		    void (covmodel::*cov1ParamPtr)(double, double &, double &),
		    void (covmodel::*cov2ParamPtr)(double, double, double &, double&), int brute);

/* knotsD       qxq knot distance matrix */
/* coordsKnotsD nxq coordinate distance matrix */
/* C            nmxnm empty */
/* C_str        qmxqm empty */
/* ct           nmxqm empty */
/* E            nmm empty */
/* n            # of coordinates */
/* m            # of response variables */
/* q            # of knots */
/* Psi          mxm residual covariance matrix */
/* V            mxm spatial covariance matrix */
/* theta        m if nPramPtr is 1 or 2*m if nPramPtr is 2 these are the spatial range parameters */
/* ... */

double mvMPPCovInvDet(double *knotsD, double *coordsKnotsD, double *C, double *C_str, double *ct, 
		      double *E, int n, int m, int q, double *Psi, double *V, double *theta, 
		      double *tmp_mm, double *tmp_mm1, double *tmp_mm2, double *tmp_nmqm, 
		      double *tmp_nmqm1, double *tmp_qmqm, double *tmp_qmqm1, int lwork, double *work,
		      string covModel, int nPramPtr, covmodel *covModelObj, 
		      void (covmodel::*cov1ParamPtr)(double, double &, double &),
		      void (covmodel::*cov2ParamPtr)(double, double, double &, double&), int brute);

/* knotsD       qxq knot distance matrix */
/* coordsKnotsD nxq coordinate distance matrix */
/* C            nmxnm empty */
/* C_str        qmxqm empty */
/* ct           nmxqm empty */
/* E            nmm empty */
/* n            # of coordinates */
/* m            # of response variables */
/* q            # of knots */
/* Psi          mxm residual covariance matrix */
/* V            mxm spatial covariance matrix */
/* theta        m if nPramPtr is 1 or 2*m if nPramPtr is 2 these are the spatial range parameters */
/* ... */

/*note: Psi must be pd*/
double mvPPCovInvDet(double *knotsD, double *coordsKnotsD, double *C, double *C_str, double *ct, 
		     double *E, int n, int m, int q, double *Psi, double *V, double *theta, 
		     double *tmp_mm, double *tmp_mm1, double *tmp_mm2, double *tmp_nmqm, 
		     double *tmp_nmqm1, double *tmp_qmqm, double *tmp_qmqm1, int lwork, double *work,
		     string covModel, int nPramPtr, covmodel *covModelObj, 
		     void (covmodel::*cov1ParamPtr)(double, double &, double &),
		     void (covmodel::*cov2ParamPtr)(double, double, double &, double&), int brute);


/* coordsD      nxn coordinate distance matrix */
/* C            nmxnm empty */
/* n            # of coordinates */
/* m            # of response variables */
/* Psi          mxm residual covariance matrix */
/* V            mxm spatial covariance matrix */
/* theta        m if nPramPtr is 1 or 2*m if nPramPtr is 2 these are the spatial range parameters */
/* ... */
/* nugget       if Psi is to be used */

double mvCovInvDet(double *coordsD, double *C, int n, int m, double *Psi, double *V, double *theta, 
		   double *tmp_mm, double *tmp_mm1, double *tmp_mm2,
		   string covModel, int nPramPtr, covmodel *covModelObj, 
		   void (covmodel::*cov1ParamPtr)(double, double &, double &),
		   void (covmodel::*cov2ParamPtr)(double, double, double &, double&));
