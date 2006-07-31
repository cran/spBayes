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
#include "vprior.h"
#include "util.h"


vprior::vprior(string fn, int dimH1, int dimH2, int nP, int nS): 
  formalName(fn), dimHyper1(dimH1), dimHyper2(dimH2), nPar(nP), nSamples(nS)
{
  nSamples++; //add an extra sample to hold the starting values
  current = 0;
  sampleOrder = -1;
  subParIndx = -1;
  acceptance = 0;
  nil = 0.0;
}

string vprior::getParName(){
  string name = formalName;
  if(subParIndx != -1)
    name += "_"+toString(subParIndx);

  return name;
}
