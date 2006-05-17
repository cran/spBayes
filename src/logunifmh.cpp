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
using namespace std;

#include "vparammh.h"
#include "logunifmh.h"


logunifmh::logunifmh(string fname, int ns): vparammh(fname, ns)
{
  prior = "LOGUNIF";
  logged = true; 
  samples.resize(nsamples, 0.0);
  priorParams.resize(2, 0.0);
}

bool logunifmh::propose(double p){
  if(exp(p) > priorParams[0] && exp(p) < priorParams[1]){
    currentSample++; //needs to come before the new proposal
    accpeted++; //this is accpeted-- if proposal is rejected
    samples[currentSample] = p;
    return false;
  }else{
    return true;//hit me again
  }
}

double logunifmh::logPrior(){
  return 0.0;
}
