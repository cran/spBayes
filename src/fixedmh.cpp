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
#include "fixedmh.h"


fixedmh::fixedmh(string fname, int ns): vparammh(fname, ns)
{
  prior = "FIXED";
  logged = false; 
  samples.resize(nsamples, 0.0);
  priorParams.resize(2, 0.0);

}

bool fixedmh::propose(double p){
  //just call this once for the fixed
  for(int i = 0; i < nsamples; i++)
    samples[i] = starting;
  return false;
}

double fixedmh::logPrior(){
  return 0.0;
}
