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

#include <Rinternals.h>
#include "vparammh.h"
#include "igmh.h"


igmh::igmh(string fname, int ns): vparammh(fname, ns)
{
  prior = "IG";
  logged = true;
  samples.resize(nsamples, 0.0);
  priorParams.resize(2, 0.0);
}
 
bool igmh::propose(double p){
  currentSample++; //needs to come before the new proposal
  accpeted++; //this is accpeted-- if proposal is rejected
  samples[currentSample] = p;
  return false;
}

double igmh::logPrior(){
  //the samples of the inverse gamma are logged for efficiency in the metrop so exp 
  return -1.0*(priorParams[0]+1.0)*samples[currentSample]-priorParams[1]/exp(samples[currentSample]);
}
