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
#include <string>
using namespace std;

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include "covmodel.h"

covmodel::covmodel(){}

void covmodel::spherical(double phi, double *r, double *d, int &length){
  int i;
  for(i = 0; i < length; i++){
    if(d[i] > 0 && d[i] <= 1.0/phi){
      r[i] = 1.0 - 1.5*phi*d[i] + 0.5*pow(phi*d[i],3);
    }else if(d[i] >= 1.0/phi){
      r[i] = 0.0;
    }else{
      r[i] = 1.0;
    }
    
  }
}

void covmodel::exponential(double phi, double *r, double *d, int &length){
  int i;
  for(i = 0; i < length; i++) r[i] = exp(-1.0*phi*d[i]);
}

void covmodel::gaussian(double phi, double *r, double *d, int &length){
  int i;
  for(i = 0; i < length; i++) r[i] = exp(-1.0*(pow(phi*d[i],2)));
}

void covmodel::matern(double phi, double nu, double *r, double *d, int &length){
  //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
  //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)

  //may be underflow problems here, but perhaps the bessel_k is smart.
  int i;
  for(i = 0; i < length; i++){
    if(d[i]*phi > 0.0)
      r[i] = pow(d[i]*phi, nu)/(pow(2, nu-1)*gammafn(nu))*bessel_k(d[i]*phi, nu, 1.0);
    else
      r[i] = 1.0;
  }
}


//
//for single values
//


void covmodel::spherical(double phi, double &r, double &d){

    if(d > 0 && d <= 1.0/phi){
      r = 1.0 - 1.5*phi*d + 0.5*pow(phi*d,3);
    }else if(d >= 1.0/phi){
      r = 0.0;
    }else{
      r = 1.0;
    }
    
}


void covmodel::exponential(double phi, double &r, double &d){
  r = exp(-1.0*phi*d);
}

void covmodel::gaussian(double phi, double &r, double &d){
    r = exp(-1.0*(pow(phi*d,2)));
}

void covmodel::matern(double phi, double nu, double &r, double &d){
  //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
  //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)

  //may be underflow problems here, but perhaps the bessel_k is smart.
    if(d*phi > 0.0)
      r = pow(d*phi, nu)/(pow(2, nu-1)*gammafn(nu))*bessel_k(d*phi, nu, 1.0);
    else
      r = 1.0;
}

