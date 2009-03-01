#include <iostream>
using namespace std;

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "hex_utils.h"



extern "C" {

  SEXP hexGrid(SEXP nodes_per_layer_r, SEXP box_r){

    int nProtect = 0;
    int info, n, layers;
    int nodes_per_layer = INTEGER(nodes_per_layer_r)[0];
    double *box = REAL(box_r);

    SEXP hx_hy_r;
    PROTECT(hx_hy_r = allocMatrix(REALSXP, 2, 1)); nProtect++;
    double *hx = &REAL(hx_hy_r)[0];
    double *hy = &REAL(hx_hy_r)[1];

    SEXP layers_r;
    PROTECT(layers_r = allocMatrix(INTSXP, 1, 1)); nProtect++;
 
    layers = hex_grid_layers(nodes_per_layer, box);
    INTEGER(layers_r)[0] = layers;

    hex_grid_h(nodes_per_layer, box, hx, hy);
    n = hex_grid_n(nodes_per_layer, box);

    SEXP p_r; PROTECT(p_r = allocMatrix(REALSXP, 2*n, 1)); nProtect++;
    double *p = REAL(p_r);

    info = hex_grid_points(nodes_per_layer, layers, n, box, p);
    if(info != 0) error("nodes per layer < 1");
    
    //return obj
    SEXP result, resultNames;
    PROTECT(result = allocVector(VECSXP, 4)); nProtect++;
    PROTECT(resultNames = allocVector(VECSXP, 4)); nProtect++;

    //samples
    SET_VECTOR_ELT(result, 0, hx_hy_r);
    SET_VECTOR_ELT(resultNames, 0, mkChar("hx.hy")); 

    SET_VECTOR_ELT(result, 1, layers_r);
    SET_VECTOR_ELT(resultNames, 1, mkChar("layers"));

    SET_VECTOR_ELT(result, 2, nodes_per_layer_r);
    SET_VECTOR_ELT(resultNames, 2, mkChar("nodes.per.layer")); 

    SET_VECTOR_ELT(result, 3, p_r);
    SET_VECTOR_ELT(resultNames, 3, mkChar("hex.centroids"));

    namesgets(result, resultNames);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
  }  
}
