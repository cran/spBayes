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

#ifndef VPARAMMH_H
#define VPARAMMH_H

#include <iostream>
#include <string>
#include <vector>
using namespace std;

class vparammh
{
 public:
  vparammh(string fname, int ns);
  ~vparammh();
  /************************
      common functions
  ************************/
  void setFormalName(string n){formalName = n;}
  string getFormalName(){return formalName;}
  string & getPrior(){return prior;}
  void setPriorParams(int indx, double p);
  double & getPriorParam(int indx){return priorParams[indx];}
  void setTuning(double t){tuning = t;}
  double getTuning(){return tuning;}
  void setStarting(double s);
  double getStarting(){return starting;}
  double getCurrentSampleTrans(); //returns the sample untransfprmed (i.e., exp if logged)
  double getCurrentSampleNoTrans(); //returns the sample as is (i.e., logged or unlogged)
  double getSampleTrans(int indx); //returns the sample untransfprmed (i.e., exp if logged)
  double getSampleNoTrans(int indx); //returns the sample as is (i.e., logged or unlogged)
  void setSampleOrder(int o){sampleOrder = o;}
  int getSampleOrder(){return sampleOrder;}
  double getAcceptRate(){return 100.0*accpeted/(1.0*nsamples);}
  void rejectProposal();
  double getSampleTransMean();
  double getSampleTransMean(int start);
  void show();



  /************************
      virtual  functions
  ************************/
  virtual bool propose(double p) = 0;
  virtual double logPrior() = 0;

 protected:

  /************************
      common variables
  ************************/
  string formalName;
  string prior;
  vector<double> priorParams;
  double starting;
  double tuning;
  int nsamples;
  vector<double> samples;
  int currentSample;
  int accpeted;
  int sampleOrder;
  bool logged;
};

#endif

