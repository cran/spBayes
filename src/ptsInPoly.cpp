#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

extern "C" {

  SEXP ptsInPoly(SEXP verts_r, SEXP nVerts_r, SEXP pts_r, SEXP nPts_r, SEXP inPtIndx_r, SEXP nInPts_r){
    
    int nVerts = INTEGER(nVerts_r)[0];
    double *verts = REAL(verts_r);

    int nPts = INTEGER(nPts_r)[0];
    double *pts = REAL(pts_r);

    int *inPtsIndx = INTEGER(inPtIndx_r);
    int *nInPts = INTEGER(nInPts_r);
    nInPts[0] = 0;

    int i, j, crossings;
    static const double eps = 1e-7;
    double px, py, x1, x2, dx, dy, k, m, y2;
    
    for(j = 0; j < nPts; j++){
      
      crossings = 0;
      
      px = pts[j];
      py = pts[nPts+j];
      
      for (i = 0; i < nVerts; i++ ){
	
	if ( verts[i] < verts[(i+1)%nVerts] ){
	  x1 = verts[i];
	  x2 = verts[(i+1)%nVerts];
	} else {
	  x1 = verts[(i+1)%nVerts];
	  x2 = verts[i];
	}
	
	if ( px > x1 && px <= x2 && (py < verts[nVerts+i] || py <= verts[nVerts+(i+1)%nVerts])){
	  
	  dx = verts[(i+1)%nVerts] - verts[i];
	  dy = verts[nVerts+(i+1)%nVerts] - verts[nVerts+i];
	  
	  if ( fabs(dx) < eps ){
	    k = R_PosInf;
	  } else {
	    k = dy/dx;
	  }
	  
	  m = verts[nVerts+i] - k * verts[i];
	  
	  y2 = k * px + m;
	  if ( py <= y2 ){
	    crossings++;
	  }
	}
      }//end vert loop
      
      if(crossings % 2 == 1 ){
	inPtsIndx[nInPts[0]] = j; //correct for R ob1 index
	nInPts[0]++;
      }
      
    }//end point loop

    return(nInPts_r);
  }  
}
