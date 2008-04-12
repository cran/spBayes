// Andrew O. Finley
// Dept. of Forest Resources
// University of Minnesota
// finleya@msu.edu 
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2004 Andrew O. Finley

#ifndef VPRIOR_H
#define VPRIOR_H

#include <iostream>
#include <string>
#include <vector>
using namespace std;
#include "util.h"

class vprior
{
 public:
  vprior(string fn, int dimH1, int dimH2, int nP, int nS);
  virtual ~vprior(){/*cout<<"Destructor: Base"<<endl;*/}
  /************************
      common functions
  ************************/
  string getFormalName(){return formalName;}
  void setSubParIndx(int sIndx){subParIndx = sIndx;}
  int getSubParIndx(){return subParIndx;}
  int getNPar(){return nPar;}
  int getnSamples(){return nSamples;}
  int getSampleOrder(){return sampleOrder;}
  void setSampleOrder(int o){sampleOrder = o;}
  int getAcceptance(){return acceptance;}
  string getParName();

  virtual void setHyper1(double *hyper1) = 0;
  virtual void setHyper2(double *hyper2) = 0;
  virtual double* getHyper1() = 0;
  virtual double* getHyper2() = 0;
  virtual void setTuning(double *t) = 0;
  virtual void setStarting(double *s) = 0;

  virtual double logLikelihood() = 0;
  virtual double logHastingsAdj() = 0;
  virtual void propose() = 0;
  virtual void reject() = 0;
  virtual double getCurrentSampleTrans(int indx) = 0;
  virtual double* getSamples() = 0;
  virtual void setCurrentSample(double*) = 0;
  virtual double* getCurrentSample() = 0;
  virtual void transSamples() = 0;
  virtual void sampleMeansFromStart(int s) = 0;
  virtual double* getSampleMeans() = 0;
  virtual void show() = 0;
  
 
  /************************
      common variables
  ************************/
  string prior;
  string formalName;
  int subParIndx;
  int dimHyper1;
  int dimHyper2;
  int nPar;
  int nSamples;
  int sampleOrder;
  int current;
  int acceptance;
  double nil;
};

#endif

