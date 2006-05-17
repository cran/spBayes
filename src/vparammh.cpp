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
#include <cmath>
#include <vector>
#include <cmath>
using namespace std;

#include <Rinternals.h>
#include "vparammh.h"


vparammh::vparammh(string fname, int ns): formalName(fname), nsamples(ns)
{
  starting = 0.0;
  tuning = 0.0;
  currentSample = 0;
  accpeted = 0;
  sampleOrder = 0;
  logged = false;
}

vparammh::~vparammh(){
  samples.clear();
  priorParams.clear();
}

void vparammh::setStarting(double s){
  if(logged){
    if(s < 0){
      error("c++ error: a starting value is less than 0 for a logged parameter, see vparammh and setStarting()");
    }
    starting = samples[0] = log(s);
  }else{
    starting = samples[0] = s;
  }
}

double vparammh::getCurrentSampleTrans(){
  if(logged)
    return exp(samples[currentSample]);
  else
    return samples[currentSample]; 
}

double vparammh::getCurrentSampleNoTrans(){
    return samples[currentSample]; 
}

double vparammh::getSampleTrans(int indx){
  if(indx >= nsamples)
      error("c++ error: index too large in getCurrentSamplsTrans");

  if(logged)
    return exp(samples[indx]);
  else
    return samples[indx]; 
}

double vparammh::getSampleNoTrans(int indx){
  if(indx >= nsamples)
      error("c++ error: index too large in getCurrentSampleNoTrans");

    return samples[indx]; 
}

//make this virtual sometime
void vparammh::setPriorParams(int indx, double p){
  priorParams[indx] = p;
}

void vparammh::rejectProposal(){
  samples[currentSample] = samples[currentSample-1];
  accpeted--;
}

double vparammh::getSampleTransMean(){
  double mean = 0;
  int i;
  for(i = 1; i < nsamples; i++){//start at 1 to avoid including starting values
    if(logged)
      mean += exp(samples[i]);
    else
      mean += samples[i];
  }

  mean = mean/(nsamples-1);
  return(mean);
}

double vparammh::getSampleTransMean(int start){
  double mean = 0;
  int i;
  for(i = start; i < nsamples; i++){
    if(logged)
      mean += exp(samples[i]);
    else
      mean += samples[i];
  }

  mean = mean/(nsamples-start);
  return(mean);
}


void vparammh::show()
{
  Rprintf("\tparameter's name:\t\t\t%s\n", formalName.c_str());
  if(logged)
    Rprintf("\tparameter's starting:\t\t\t%f\n", exp(starting));
  else
    Rprintf("\tparameter's starting:\t\t\t%f\n", starting);
  Rprintf("\tparameter's tuning:\t\t\t%f\n", tuning);
  Rprintf("\tparameter's prior:\t\t\t%s\n", prior.c_str());
  Rprintf("\tprior's parameter(s):\n");
  Rprintf("\t\t\t");
  for(int i = 0; i < priorParams.size(); i++)
    Rprintf("\t%f", priorParams[i]);
  
  //Rprintf("\n\tparameter's current sample index:\t%d\n", currentSample);
  //Rprintf("\tparameter's accepted:\t\t\t%d\n", accpeted);
  Rprintf("\n\tparameter's relative sample order:\t%d\n", sampleOrder);
  Rprintf("\n-------------\n");
}
