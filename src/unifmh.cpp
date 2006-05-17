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
#include "unifmh.h"


unifmh::unifmh(string fname, int ns): vparammh(fname, ns)
{
  prior = "UNIF";
  logged = false; 
  samples.resize(nsamples, 0.0);
  priorParams.resize(2, 0.0);
}

bool unifmh::propose(double p){
  if(p > priorParams[0] && p < priorParams[1]){
    currentSample++; //needs to come before the new proposal
    accpeted++; //this is accpeted-- if proposal is rejected
    samples[currentSample] = p;
    return false;
  }else{
    return true;//hit me again
  }
}

double unifmh::logPrior(){
  return 0.0;
}
