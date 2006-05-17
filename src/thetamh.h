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

#ifndef THETAMH_H
#define THETAMH_H

#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "vparammh.h"

class thetamh : public vparammh
{
 public:
  thetamh(string fname, int ns);
  
  /************************
      virtual  functions
  ************************/
  virtual bool propose(double p);
  virtual double logPrior();
  
 private:
};
#endif

