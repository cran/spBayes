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
#include "hc.h"


hc::hc(string fname, int ns): vparammh(fname, ns)
{
  prior = "HC";
  logged = true;
  samples.resize(nsamples, 0.0);
  priorParams.resize(1, 0.0);
}
 
bool hc::propose(double p){
  currentSample++; //needs to come before the new proposal
  accpeted++; //this is accpeted-- if proposal is rejected
  samples[currentSample] = p;
  return false;
}

double hc::logPrior(){
  //the samples of the half-cauchy are logged for efficiency in the metrop so exp 
  return -1.0*log(pow(exp(samples[currentSample]), 2) + pow(priorParams[0],2));
}
