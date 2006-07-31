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

#ifndef COVMODEL_H
#define COVMODEL_H

#include <iostream>
#include <string>
using namespace std;


class covmodel
{
 public:
  covmodel();
  void spherical(double phi, double *r, double *d, int &length);
  void exponential(double phi, double *r, double *d, int &length);
  void gaussian(double phi, double *r, double *d, int &length);
  void matern(double phi, double nu, double *r, double *d, int &length);


  void spherical(double phi, double &r, double &d);
  void exponential(double phi, double &r, double &d);
  void gaussian(double phi, double &r, double &d);
  void matern(double phi, double nu, double &r, double &d);

 private:
};
#endif

