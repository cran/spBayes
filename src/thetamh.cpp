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
#include <vector>
#include <cmath>
using namespace std;

#include "vparammh.h"
#include "thetamh.h"


thetamh::thetamh(string fname, int ns): vparammh(fname, ns)
{
  prior = "THETA";
  logged = false; 
  samples.clear();
  priorParams.clear();
}

bool thetamh::propose(double p){
  return false;
}

double thetamh::logPrior(){
  return 0.0;
}
