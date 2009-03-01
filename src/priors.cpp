#include <iostream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "util.h"
#include "vprior.h"
#include "priors.h"



/************************
     Inverse-gamma 
************************/
ig::ig(string fn, int dimH1, int dimH2, int nP, int nS): 
  vprior(fn, dimH1, dimH2, nP, nS)
{
  prior = "IG";
  samples = (double *) R_alloc(nPar*nSamples, sizeof(double)); zeros(samples, nPar*nSamples);
  shape = (double *) R_alloc(nPar, sizeof(double)); 
  scale = (double *) R_alloc(nPar, sizeof(double));
  tuning = (double *) R_alloc(nPar, sizeof(double));
  sampleMeans = (double *) R_alloc(nPar, sizeof(double));
}

void ig::setHyper1(double *hyper1){
  for(int i = 0; i < nPar; i++)
    shape[i] = hyper1[i];
}

void ig::setHyper2(double *hyper2){
  for(int i = 0; i < nPar; i++)
    scale[i] = hyper2[i];
}

double* ig::getHyper1(){return shape;}

double* ig::getHyper2(){return scale;}

void ig::setTuning(double *t){
  for(int i = 0; i < nPar; i++)
    tuning[i] = t[i];
}

void ig::setStarting(double *s){
  current = 0;
  for(int i = 0; i < nPar; i++){
    if(s[i] <= 0)
      error("c++ error: ig::setStarting given starting values is <= 0");
    samples[i] = log(s[i]);
  }
}

double ig::logLikelihood(){
  double logLike = 0;
  
  //with Jacobian
  for(int i = 0; i < nPar; i++)
    logLike += -1.0*(shape[i]+1.0)*samples[current*nPar+i] - scale[i]/exp(samples[current*nPar+i]) + samples[current*nPar+i];

  return logLike;
}


double ig::logHastingsAdj(){
  return 0.0;
}


void ig::propose(){
  current++;
  acceptance++;
  for(int i = 0; i < nPar; i++)
    samples[current*nPar+i] = rnorm(samples[(current-1)*nPar+i], tuning[i]);
}

void ig::reject(){
  int incOne = 1;
  acceptance--;
  F77_NAME(dcopy)(&nPar, &samples[(current-1)*nPar], &incOne, &samples[current*nPar], &incOne);
}

double ig::getCurrentSampleTrans(int indx){
  if(indx > nPar-1)
    error("c++ error: ig::getCurrentSampleTrans");
  return exp(samples[current*nPar+indx]);
}

double* ig::getSamples(){return samples;}

void ig::setCurrentSample(double*){
  error("c++ error: ig::setCurrentSample not implemented");
}

double* ig::getCurrentSample(){
  error("c++ error: ig::getCurrentSample not implemented");
  return &nil;
}

void ig::transSamples(){
  for(int i = 0; i < nSamples; i++)
    for(int j = 0; j < nPar; j++)
      samples[i*nPar+j] = exp(samples[i*nPar+j]);
}

void ig::sampleMeansFromStart(int s){
  int n = 0;
  int i,j;
  for(j = 0; j < nPar; j++)
    sampleMeans[j] = 0.0;

  for(i = s; i < nSamples; i++){
    n++;
    for(j = 0; j < nPar; j++){
      sampleMeans[j] += exp(samples[i*nPar+j]);
    }
  }

  for(j = 0; j < nPar; j++)
    sampleMeans[j] = sampleMeans[j]/n;
}

double* ig::getSampleMeans(){return sampleMeans;}

void ig::show(){

  if(subParIndx == -1)
    Rprintf("name: %s\n", formalName.c_str());
  else
    Rprintf("name: %s_%i\n", formalName.c_str(), subParIndx);

  Rprintf("\tnumber of samples: %i\n", nSamples-1);//ob1 for starting
  Rprintf("\tnumber of parameters: %i\n", nPar);
  Rprintf("\tsample order: %i\n", sampleOrder);
  Rprintf("\tprior: %s\n", prior.c_str());
  Rprintf("\t\tparam. number and shape, scale hyperparameter:\n");
  for(int i = 0; i < nPar; i++)
    Rprintf("\t\t\t%i %f %f\n", i+1, shape[i], scale[i]);
  
  Rprintf("\tMetropolis-Hastings param. number, tuning, and starting value:\n");
  for(int i = 0; i < nPar; i++)
    Rprintf("\t\t%i %f %f\n", i+1, tuning[i], exp(samples[i]));
  Rprintf("-------------------------------------------------\n");

}

/************************
     Inverse-Wishart 
************************/
iwish::iwish(string fn, int dimH1, int dimH2, int nP, int nS): 
  vprior(fn, dimH1, dimH2, nP, nS)
{
  prior = "IWISH";
  S = (double *) R_alloc(dimHyper2*dimHyper2, sizeof(double));
  tuning = (double *) R_alloc(nPar*nPar, sizeof(double));
  samples = (double *) R_alloc(nPar*nSamples, sizeof(double)); zeros(samples, nPar*nSamples);
  A = (double *) R_alloc(dimHyper2*dimHyper2, sizeof(double));
  Atmp = (double *) R_alloc(dimHyper2*dimHyper2, sizeof(double));
  sampleMeans = (double *) R_alloc(nPar, sizeof(double));
}

void iwish::setHyper1(double *hyper1){
  df = *hyper1;
}

void iwish::setHyper2(double *hyper2){
  for(int i = 0; i < dimHyper2*dimHyper2; i++) S[i] = hyper2[i];
}

double* iwish::getHyper1(){return &df;}

double* iwish::getHyper2(){return S;}

void iwish::setTuning(double *t){
  for(int i = 0; i < nPar*nPar; i++) tuning[i] = t[i];
}

void iwish::setStarting(double *s){
  int i, j, k;
  current = 0;

  for(i = 0, k = 0; i < dimHyper2; i++){
    for(j = i; j < dimHyper2; j++, k++){
      samples[k] = s[k];
      if(i == j){
	if(samples[k] <= 0.0)
	  error("c++ error: iwish::setStarting given starting value on diag <= 0");
	samples[k] = log(samples[k]);
      }
    }
  }
}

double iwish::logLikelihood(){

  double logLike = 0;
  double logDetK = 0;
  double SKtrace = 0;
  int m = dimHyper2;
  char lower = 'L';
  char right = 'R';
  double zero = 0.0;
  double one = 1.0;
  int info, i;
  
  //get current samples into A and trans the diag
  setA();
  
  //get logDetK
  for(i = 0; i < m; i++)
    logDetK += log(A[i*m+i]);
  logDetK = 2*logDetK;
  
  //jacobian \sum_{i=1}^{m} (m-i+1)*log(a_ii)+log(a_ii)
  for(i = 0; i < m; i++) 
    logLike += (m-i)*log(A[i*m+i])+log(A[i*m+i]);
  
  //get S*K^-1, already have the chol of K (i.e., A)
  F77_CALL(dpotri)(&lower, &m, A, &m, &info); if(info != 0){error("c++ error: Cholesky inverse failed\n");}
  
  F77_CALL(dsymm)(&right, &lower, &m, &m, &one, A, &m, S, &m, &zero, Atmp, &m);
  
  for(i = 0; i < m; i++){SKtrace += Atmp[i*m+i];}
  
  logLike += -0.5*(df+m+1)*logDetK - 0.5*SKtrace;
  
  return logLike;
}

double iwish::logHastingsAdj(){
  return 0.0;
}

void iwish::propose(){
  current++;
  acceptance++;
  mvrnorm(&samples[current*nPar], &samples[(current-1)*nPar], tuning, nPar);
}

void iwish::reject(){
  int incOne = 1;
  acceptance--;
  F77_NAME(dcopy)(&nPar, &samples[(current-1)*nPar], &incOne, &samples[current*nPar], &incOne);
}

double iwish::getCurrentSampleTrans(int indx){
  return nil;
}

double* iwish::getSamples(){return samples;}

void iwish::setCurrentSample(double*){
  error("c++ error: iwish::setCurrentSample not implemented");
}

double* iwish::getCurrentSample(){
  setA();
  return A;
}

void iwish::setA(){
  int i, j , k;

  //lower tri A
  zeros(A, dimHyper2*dimHyper2);
  for(i = 0, k = 0; i < dimHyper2; i++){
    for(j = i; j < dimHyper2; j++, k++){
      A[i*dimHyper2+j] = samples[current*nPar+k];
      if(i == j)
	A[i*dimHyper2+j] = exp(A[i*dimHyper2+j]);
    }
  }

}

void iwish::transSamples(){
  int indx = 0;
  for(int i = 0; i < nSamples; i++){
    indx = 0;
    for(int j = 0; j < dimHyper2; j++){
      samples[i*nPar+indx] =  exp(samples[i*nPar+indx]);
      indx = indx+dimHyper2-j;
    }
  }
}


void iwish::sampleMeansFromStart(int s){
  int n = 0;
  int diag,i,j,k;

  for(j = 0; j < nPar; j++)
    sampleMeans[j] = 0.0;

  //first all sum all elements
  for(i = s; i < nSamples; i++)
    for(j = 0; j < nPar; j++)
      sampleMeans[j] += samples[i*nPar+j];
  
  //zero the diag elements
  diag = 0;
  for(j = 0; j < dimHyper2; j++){
    sampleMeans[diag] = 0.0;
    diag = diag+dimHyper2-j;
  }

  //now get the diag's trans sum
  diag = 0;
  for(i = s; i < nSamples; i++){
    n++;
    diag = 0;
    for(j = 0; j < dimHyper2; j++){
       sampleMeans[diag] += exp(samples[i*nPar+diag]);
      diag = diag+dimHyper2-j;
    }
  }

  //finally mean
  for(j = 0; j < nPar; j++)
    sampleMeans[j] = sampleMeans[j]/n;
}

double* iwish::getSampleMeans(){
  int i,j,k;

  zeros(A, dimHyper2*dimHyper2);
  for(i = 0, k = 0; i < dimHyper2; i++)
    for(j = i; j < dimHyper2; j++, k++)
      A[i*dimHyper2+j] = sampleMeans[k];
  
  return A;
}

void iwish::show(){
  int i, j, k;
  if(subParIndx == -1)
    Rprintf("name: %s\n", formalName.c_str());
  else
    Rprintf("name: %s_%i\n", formalName.c_str(), subParIndx);
  Rprintf("\tnumber of samples: %i\n", nSamples-1);//ob1 for starting
  Rprintf("\tnumber of parameters: %i\n", nPar);
  Rprintf("\tsample order: %i\n", sampleOrder);
  Rprintf("\tprior: %s\n", prior.c_str());
  Rprintf("\t\tdf hyperparameter: %f\n", df);
  Rprintf("\t\tscale hyperparameter:\n");
  Rprintf("\t\t\t");  
  for(i = 0; i < dimHyper2; i++){
    for(j = 0; j < dimHyper2; j++){
      Rprintf("%f ", S[j*dimHyper2+i]);
    }
    Rprintf("\n\t\t\t");     
  }
  Rprintf("\n\t"); 
  Rprintf("Metropolis-Hastings tuning:\n");
  Rprintf("\t\t");
  //show them what they entered not the lower tri.
  for(i = 0; i < nPar; i++){
    for(j = 0; j < nPar; j++){
      Rprintf("%f ", tuning[i*nPar+j]);
    }
    Rprintf("\n\t\t");     
  }
  Rprintf("\n\t"); 

  setA();
  Rprintf("parameter starting value:\n");
  Rprintf("\t\t");
  for(i = 0; i < dimHyper2; i++){
    for(j = 0; j < dimHyper2; j++){
      Rprintf("%f ", A[j*dimHyper2+i]);
    }
    Rprintf("\n\t\t"); 
  }    
  Rprintf("\n"); 
  Rprintf("-------------------------------------------------\n");
}
 

/************************
         Uniform 
************************/
unif::unif(string fn, int dimH1, int dimH2, int nP, int nS): 
  vprior(fn, dimH1, dimH2, nP, nS)
{
  prior = "UNIF";
  samples = (double *) R_alloc(nPar*nSamples, sizeof(double)); zeros(samples, nPar*nSamples);
  a = (double *) R_alloc(nPar, sizeof(double)); 
  b = (double *) R_alloc(nPar, sizeof(double));
  tuning = (double *) R_alloc(nPar, sizeof(double));
  sampleMeans = (double *) R_alloc(nPar, sizeof(double));
}

void unif::setHyper1(double *hyper1){
  for(int i = 0; i < nPar; i++)
    a[i] = hyper1[i];
}

void unif::setHyper2(double *hyper2){
  for(int i = 0; i < nPar; i++)
    b[i] = hyper2[i];
}

double* unif::getHyper1(){return a;}

double* unif::getHyper2(){return b;}

void unif::setTuning(double *t){
  for(int i = 0; i < nPar; i++)
    tuning[i] = t[i];
}

void unif::setStarting(double *s){
  current = 0;
  for(int i = 0; i < nPar; i++){
    samples[i] = s[i];
  }
}

double unif::logLikelihood(){
  return 0;
}


double unif::logHastingsAdj(){
  double adj = 0.0;
  
  for(int i = 0; i < nPar; i++){
    adj += log(dTNorm(samples[(current-1)*nPar+i], samples[current*nPar+i], tuning[i], a[i], b[i]))-
      log(dTNorm(samples[current*nPar+i], samples[(current-1)*nPar+i], tuning[i], a[i], b[i]));
  }

  return adj;
}


void unif::propose(){
  current++;
  acceptance++;
  bool again = true;

  for(int i = 0; i < nPar; i++){
    again = true;
    while(again){
      samples[current*nPar+i] = rnorm(samples[(current-1)*nPar+i], tuning[i]);
      if(samples[current*nPar+i] > a[i] && samples[current*nPar+i] < b[i])
	again = false;
      
      R_CheckUserInterrupt();
    }
  }

}

void unif::reject(){
  int incOne = 1;
  acceptance--;
  F77_NAME(dcopy)(&nPar, &samples[(current-1)*nPar], &incOne, &samples[current*nPar], &incOne);
}

double unif::getCurrentSampleTrans(int indx){
  if(indx > nPar-1)
    error("c++ error: unif::getCurrentSampleTrans");
  
  return samples[current*nPar+indx];
}

double* unif::getSamples(){return samples;}

void unif::setCurrentSample(double*){
  error("c++ error: unif::setCurrentSample not implemented");
}

void unif::transSamples(){}

double* unif::getCurrentSample(){
  error("c++ error: iwish::getCurrentSample not implemented");
  return &nil;
}

void unif::sampleMeansFromStart(int s){
  int n = 0;
  int i,j;
  for(j = 0; j < nPar; j++)
    sampleMeans[j] = 0.0;

  for(i = s; i < nSamples; i++){
    n++;
    for(j = 0; j < nPar; j++){
      sampleMeans[j] += samples[i*nPar+j];
    }
  }

  for(j = 0; j < nPar; j++)
    sampleMeans[j] = sampleMeans[j]/n;
}

double* unif::getSampleMeans(){return sampleMeans;}

void unif::show(){
  string tmp;

  if(subParIndx == -1)
    Rprintf("name: %s\n", formalName.c_str());
  else
    Rprintf("name: %s_%i\n", formalName.c_str(), subParIndx);
  Rprintf("\tnumber of samples: %i\n", nSamples-1);//ob1 for starting
  Rprintf("\tnumber of parameters: %i\n", nPar);
  Rprintf("\tsample order: %i\n", sampleOrder);
  Rprintf("\tprior: %s\n", prior.c_str());
  Rprintf("\t\tparam. number and a and b hyperparameter:\n");
  for(int i = 0; i < nPar; i++)
    Rprintf("\t\t\t%i %f %f\n", i+1, a[i], b[i]);
  
  Rprintf("\tMetropolis-Hastings param. number, tuning, and starting value:\n");
  for(int i = 0; i < nPar; i++)
    Rprintf("\t\t%i %f %f\n", i+1, tuning[i], samples[i]);
  Rprintf("-------------------------------------------------\n");

}



/************************
       Half-Cauchy
************************/
hc::hc(string fn, int dimH1, int dimH2, int nP, int nS): 
  vprior(fn, dimH1, dimH2, nP, nS)
{
  prior = "HC";
  samples = (double *) R_alloc(nPar*nSamples, sizeof(double)); zeros(samples, nPar*nSamples);
  a = (double *) R_alloc(nPar, sizeof(double)); 
  tuning = (double *) R_alloc(nPar, sizeof(double));
  sampleMeans = (double *) R_alloc(nPar, sizeof(double));
}

void hc::setHyper1(double *hyper1){
  for(int i = 0; i < nPar; i++)
    a[i] = hyper1[i];
}

void hc::setHyper2(double *hyper2){
  error("c++ error: hc::seyHyper2 not implemented");
}

double* hc::getHyper1(){return a;}

double* hc::getHyper2(){
  error("c++ error: hc::getHyper2 not implemented");
  return &nil;
}

void hc::setTuning(double *t){
  for(int i = 0; i < nPar; i++)
    tuning[i] = t[i];
}

void hc::setStarting(double *s){
  current = 0;
  for(int i = 0; i < nPar; i++){
    if(s[i] <= 0)
      error("c++ error: hc::setStarting given starting values is <= 0");
    samples[i] = log(s[i]);
  }
}

double hc::logLikelihood(){
  double logLike = 0;
  
  //with Jacobian
  for(int i = 0; i < nPar; i++)
    logLike += -1.0*log(exp(samples[current*nPar+i]) + pow(a[i],2)) + samples[current*nPar+i];

  return logLike;
}


double hc::logHastingsAdj(){
  return 0.0;
}


void hc::propose(){
  current++;
  acceptance++;
  for(int i = 0; i < nPar; i++)
    samples[current*nPar+i] = rnorm(samples[(current-1)*nPar+i], tuning[i]);
}

void hc::reject(){
  int incOne = 1;
  acceptance--;
  F77_NAME(dcopy)(&nPar, &samples[(current-1)*nPar], &incOne, &samples[current*nPar], &incOne);
}

double hc::getCurrentSampleTrans(int indx){
  if(indx > nPar-1)
    error("c++ error: hc::getCurrentSampleTrans");
  
  return exp(samples[current*nPar+indx]);
}

double* hc::getSamples(){
  return samples;
}

void hc::setCurrentSample(double*){
  error("c++ error: hc::setCurrentSample not implemented");
}

double* hc::getCurrentSample(){
  error("c++ error: hc::getCurrentSample not implemented");
  return &nil;
}

void hc::transSamples(){
  for(int i = 0; i < nSamples; i++)
    for(int j = 0; j < nPar; j++)
      samples[i*nPar+j] = exp(samples[i*nPar+j]);
}

void hc::sampleMeansFromStart(int s){
  int n = 0;
  int i,j;

  for(j = 0; j < nPar; j++)
    sampleMeans[j] = 0.0;

  for(i = s; i < nSamples; i++){
    n++;
    for(j = 0; j < nPar; j++){
      sampleMeans[j] += exp(samples[i*nPar+j]);
    }
  }

  for(j = 0; j < nPar; j++)
    sampleMeans[j] = sampleMeans[j]/n;
}

double* hc::getSampleMeans(){return sampleMeans;}

void hc::show(){

  if(subParIndx == -1)
    Rprintf("name: %s\n", formalName.c_str());
  else
    Rprintf("name: %s_%i\n", formalName.c_str(), subParIndx);
  Rprintf("\tnumber of samples: %i\n", nSamples-1);//ob1 for starting
  Rprintf("\tnumber of parameters: %i\n", nPar);
  Rprintf("\tsample order: %i\n", sampleOrder);
  Rprintf("\tprior: %s\n", prior.c_str());
  Rprintf("\t\tparam. number and a hyperparameter:\n");
  for(int i = 0; i < nPar; i++)
    Rprintf("\t\t\t%i %f\n", i+1, a[i]);
  
  Rprintf("\tMetropolis-Hastings param. number, tuning, and starting value:\n");
  for(int i = 0; i < nPar; i++)
    Rprintf("\t\t%i %f %f\n", i+1, tuning[i], exp(samples[i]));
  Rprintf("-------------------------------------------------\n");

}


/************************
       Fixed
************************/
fixedpar::fixedpar(string fn, int dimH1, int dimH2, int nP, int nS): 
  vprior(fn, dimH1, dimH2, nP, nS)
{
  prior = "FIXED";
  samples = (double *) R_alloc(nPar*nSamples, sizeof(double)); zeros(samples, nPar*nSamples);
  sampleMeans = (double *) R_alloc(nPar, sizeof(double));
}

void fixedpar::setHyper1(double *hyper1){
  error("c++ error: fixedpar::setHyper1 not implemented");
}

void fixedpar::setHyper2(double *hyper2){
  error("c++ error: fixedpar::setHyper2 not implemented");
}

double* fixedpar::getHyper1(){
  error("c++ error: fixedpar::getHyper1 not implemented");
  return &nil;
}

double* fixedpar::getHyper2(){
  error("c++ error: fixedpar::getHyper2 not implemented");
  return &nil;
}

void fixedpar::setTuning(double *t){}

void fixedpar::setStarting(double *s){
  int i;
  int incOne = 1;

  for(i = 0; i < nPar; i++){
    samples[i] = s[i];
  }

  //just fill up the samples with fixed
  for(i = 1; i < nSamples; i++)
    F77_NAME(dcopy)(&nPar, &samples[(i-1)*nPar], &incOne, &samples[i*nPar], &incOne);
}

double fixedpar::logLikelihood(){
  return 0;
}

double fixedpar::logHastingsAdj(){
  return 0.0;
}


void fixedpar::propose(){}

void fixedpar::reject(){}

double fixedpar::getCurrentSampleTrans(int indx){
  if(indx > nPar-1)
    error("c++ error: fixedpar::getCurrentSampleTrans");
  
  return samples[current*nPar+indx];
}

double* fixedpar::getSamples(){
  return samples;
}

void fixedpar::setCurrentSample(double*){
  error("c++ error: fixedpar::setCurrentSample not implemented");
}

double* fixedpar::getCurrentSample(){
  error("c++ error: fixedpar::getCurrentSample not implemented");
  return &nil;
}

void fixedpar::transSamples(){}

void fixedpar::sampleMeansFromStart(int s){
  int n = 0;
  int j;

  for(j = 0; j < nPar; j++)
    sampleMeans[j] = samples[j];
}

double* fixedpar::getSampleMeans(){return sampleMeans;}

void fixedpar::show(){

  if(subParIndx == -1)
    Rprintf("name: %s\n", formalName.c_str());
  else
    Rprintf("name: %s_%i\n", formalName.c_str(), subParIndx);
  Rprintf("\tnumber of samples: %i\n", nSamples-1);//ob1 for starting
  Rprintf("\tnumber of fixed parameters: %i\n", nPar);
  Rprintf("\tsample order (null): %i\n", sampleOrder);
  Rprintf("\tparam. number and fixed value:\n");
  for(int i = 0; i < nPar; i++)
    Rprintf("\t\t%i %f\n", i+1, samples[i]);
  Rprintf("-------------------------------------------------\n");
}

/************************
     Fixed matrix
************************/

fixedmtrx::fixedmtrx(string fn, int dimH1, int dimH2, int nP, int nS): 
  vprior(fn, dimH1, dimH2, nP, nS)
{
  prior = "FIXEDMTRX";
  samples = (double *) R_alloc(nPar*nSamples, sizeof(double)); zeros(samples, nPar*nSamples);
  A = (double *) R_alloc(dimHyper2*dimHyper2, sizeof(double));
  sampleMeans = (double *) R_alloc(nPar, sizeof(double));
}

void fixedmtrx::setHyper1(double *hyper1){
  error("c++ error: fixedmtrx::setHyper1 not implemented");
}

void fixedmtrx::setHyper2(double *hyper2){
  error("c++ error: fixedmtrx::setHyper2 not implemented");
}

double* fixedmtrx::getHyper1(){
  error("c++ error: fixedmtrx::getHyper1 not implemented");
  return &nil;
}

double* fixedmtrx::getHyper2(){
  error("c++ error: fixedmtrx::getHyper2 not implemented");
  return &nil;
}

void fixedmtrx::setTuning(double *t){}

void fixedmtrx::setStarting(double *s){
  int i, j , k;
  int incOne = 1;

  for(i = 0, k = 0; i < dimHyper2; i++){
    for(j = i; j < dimHyper2; j++, k++){
      samples[k] = s[k];
    }
  }

  //just fill up the samples with fixed
  for(i = 1; i < nSamples; i++)
    F77_NAME(dcopy)(&nPar, &samples[(i-1)*nPar], &incOne, &samples[i*nPar], &incOne);

}

double fixedmtrx::logLikelihood(){
  return 0;
}

double fixedmtrx::logHastingsAdj(){
  return 0.0;
}

void fixedmtrx::propose(){}

void fixedmtrx::reject(){}

double fixedmtrx::getCurrentSampleTrans(int indx){
  error("c++ error: fixedmtrx::getCurrentSampleTrans not implemented");
  return nil;
}

double* fixedmtrx::getSamples(){
  return samples;
}

void fixedmtrx::setCurrentSample(double*){
  error("c++ error: fixedmtrx::setCurrentSample not implemented");
}

double* fixedmtrx::getCurrentSample(){
  setLowerTriA();
  return A;
}

void fixedmtrx::setLowerTriA(){
  int i, j , k;

  //lower tri A
  zeros(A, dimHyper2*dimHyper2);
  for(i = 0, k = 0; i < dimHyper2; i++){
    for(j = i; j < dimHyper2; j++, k++){
      A[i*dimHyper2+j] = samples[current*nPar+k];
    }
  }
}

void fixedmtrx::transSamples(){}

void fixedmtrx::sampleMeansFromStart(int s){
  int n = 0;
  int j;
  for(j = 0; j < nPar; j++)
    sampleMeans[j] = samples[j];
}

double* fixedmtrx::getSampleMeans(){
    int i, j, k;
  //lower tri A
  zeros(A, dimHyper2*dimHyper2);
  for(i = 0, k = 0; i < dimHyper2; i++){
    for(j = i; j < dimHyper2; j++, k++){
      A[i*dimHyper2+j] = sampleMeans[k];
    }
  }

  return A;
}

void fixedmtrx::show(){
  int i, j, k;

  if(subParIndx == -1)
    Rprintf("name: %s\n", formalName.c_str());
  else
    Rprintf("name: %s_%i\n", formalName.c_str(), subParIndx);
  Rprintf("\tnumber of samples: %i\n", nSamples-1);//ob1 for starting
  Rprintf("\tnumber of parameters: %i\n", nPar);
  Rprintf("\tsample order (null): %i\n", sampleOrder);

  setLowerTriA(); //currently I only use lower tri
  Rprintf("parameter fixed value:\n");
  Rprintf("\t\t");
  for(i = 0; i < dimHyper2; i++){
    for(j = 0; j < dimHyper2; j++){
      Rprintf("%f ", A[j*dimHyper2+i]);
    }
    Rprintf("\n\t\t"); 
  }    
  Rprintf("\n"); 
  Rprintf("-------------------------------------------------\n");
}
 

/************************
     Beta parameter
************************/
betapar::betapar(string fn, int dimH1, int dimH2, int nP, int nS, bool n, bool g): 
  vprior(fn, dimH1, dimH2, nP, nS)
{
  prior = "Beta";
  normalPrior = n;
  gibbs = g;
  
  if(normalPrior){
    mu =  (double *) R_alloc(nPar, sizeof(double)); zeros(mu, nPar);
    sigma = (double *) R_alloc(nPar*nPar, sizeof(double)); zeros(sigma, nPar*nPar);
    tmpNParVec =  (double *) R_alloc(nPar, sizeof(double)); zeros(tmpNParVec, nPar);
    tmpNParVec1 =  (double *) R_alloc(nPar, sizeof(double)); zeros(tmpNParVec1, nPar);
  }else{
    mu = NULL;
    sigma = NULL;
  }

  if(!gibbs){
    tuning = (double *) R_alloc(nPar*nPar, sizeof(double)); zeros(tuning, nPar*nPar);
  }else{
    tuning = NULL;
  }

  samples = (double *) R_alloc(nPar*nSamples, sizeof(double)); zeros(samples, nPar*nSamples);
  sampleMeans = (double *) R_alloc(nPar, sizeof(double));
}

void betapar::setHyper1(double *hyper1){
  if(normalPrior){
    for(int i = 0; i < nPar; i++) mu[i] = hyper1[i];
  }else{
    error("c++ error: betapar::setHyper1 normal prior not in use");
  }

}

void betapar::setHyper2(double *hyper2){
  if(normalPrior){
    for(int i = 0; i < nPar*nPar; i++) sigma[i] = hyper2[i];
  }else{
    error("c++ error: betapar::setHyper2 normal prior not in use");
  }
}

double* betapar::getHyper1(){return mu;}

double* betapar::getHyper2(){return sigma;}

void betapar::setTuning(double *t){
  if(!gibbs){
    for(int i = 0; i < nPar*nPar; i++) tuning[i] = t[i];
  }else{
    error("c++ error: betapar::setTuning MH not in use");
  }
}

void betapar::setStarting(double *s){
  current = 0;
  for(int i = 0; i < nPar; i++)
    samples[i] = s[i];
}

double betapar::logLikelihood(){

  if(normalPrior){
    int incOne = 1;
    double zero = 0.0;
    double one = 1.0;
    double negOne = -1.0;
    char ntran = 'N';

    F77_NAME(dcopy)(&nPar, mu, &incOne, tmpNParVec, &incOne);
    F77_NAME(daxpy)(&nPar, &negOne, &samples[current*nPar], &incOne, tmpNParVec, &incOne);
    F77_NAME(dgemv)(&ntran, &nPar, &nPar, &one, sigma, &nPar, tmpNParVec, &incOne, &zero, tmpNParVec1, &incOne);
    return -0.5*F77_NAME(ddot)(&nPar, tmpNParVec, &incOne, tmpNParVec1, &incOne);

  }else{//flat
    return 0.0;
  }
}

double betapar::logHastingsAdj(){
  return 0.0;
}

void betapar::propose(){
  //no effect if Gibbs
  if(!gibbs){//only used for MH
    current++;
    acceptance++;
    mvrnorm(&samples[current*nPar], &samples[(current-1)*nPar], tuning, nPar);
  }
}

void betapar::reject(){
  //no effect if Gibbs
  if(!gibbs){
    acceptance--;
    int incOne = 1;
    F77_NAME(dcopy)(&nPar, &samples[(current-1)*nPar], &incOne, &samples[current*nPar], &incOne);
  }
}

double betapar::getCurrentSampleTrans(int indx){
  error("c++ error: betapar::getCurrentSampleTrans not implemented");
  return nil;
}

double* betapar::getSamples(){
  return samples;
}

void betapar::setCurrentSample(double *beta){
  if(gibbs){
    int incOne = 1;
    current++; 
    F77_NAME(dcopy)(&nPar, beta, &incOne, &samples[current*nPar], &incOne);
  }else{
    error("c++ error: betapar::setCurrentBeta not Gibbs");
  }
}

double* betapar::getCurrentSample(){
   return &samples[current*nPar];
}

void betapar::transSamples(){}

void betapar::sampleMeansFromStart(int s){
  int n = 0;
  int i, j;

  for(j = 0; j < nPar; j++)
    sampleMeans[j] = 0.0;

  for(i = s; i < nSamples; i++){
    n++;
    for(j = 0; j < nPar; j++){
      sampleMeans[j] += samples[i*nPar+j];
    }
  }

  for(j = 0; j < nPar; j++)
    sampleMeans[j] = sampleMeans[j]/n;
}

double* betapar::getSampleMeans(){return sampleMeans;}

void betapar::show(){
  int i, j, k;

  Rprintf("parameter vector name: %s\n", formalName.c_str());
  if(gibbs)
    Rprintf("update method: Gibbs\n");
  else
    Rprintf("\tupdate method: Metropolis-Hastings\n");    
  Rprintf("\tnumber of samples: %i\n", nSamples-1);//ob1 for starting
  Rprintf("\tnumber of parameters: %i\n", nPar);
  Rprintf("\tsample order: %i\n", sampleOrder);
  if(normalPrior)
    Rprintf("\tprior: NORMAL\n");
  else
    Rprintf("\tprior: FLAT\n");

  if(normalPrior){
    Rprintf("\t\tmu hyperparameter:\n");
    Rprintf("\t\t\t");
    for(i = 0; i < nPar; i++)
      Rprintf("%f ", mu[i]);
    Rprintf("\n");
      
    Rprintf("\t\tSigma hyperparameter:\n");
    Rprintf("\t\t\t");
    for(i = 0; i < nPar; i++){
      for(j = 0; j < nPar; j++){
	Rprintf("%f ", sigma[j*nPar+i]);
      }
      Rprintf("\n\t\t\t");     
    }
  }
  Rprintf("\n\t");

  if(!gibbs){
    Rprintf("Metropolis-Hastings tuning:\n");
    Rprintf("\t\t");
    //show them what they entered not the lower tri.
    for(i = 0; i < nPar; i++){
      for(j = 0; j < nPar; j++){
	Rprintf("%f ", tuning[i*nPar+j]);
      }
      Rprintf("\n\t\t");     
    }
    Rprintf("\n\t"); 
  }
  
  Rprintf("parameter starting value:\n");
  Rprintf("\t\t");
  for(i = 0; i < nPar; i++)
    Rprintf("%f ", samples[i]);
  
  Rprintf("\n"); 
  Rprintf("-------------------------------------------------\n");
}
 
 
