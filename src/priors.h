#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "vprior.h"

/************************
     Inverse-gamma
************************/
#ifndef IG_H
#define IG_H

class ig : public vprior
{
 public:
  ig(string fn, int dimH1, int dimH2, int nP, int nS);
  ~ig(){/*cout<<"Destructor: igmh"<<endl;*/}

  virtual void setHyper1(double *hyper1);
  virtual void setHyper2(double *hyper2);
  virtual double* getHyper1();
  virtual double* getHyper2();
  virtual void setTuning(double *t);
  virtual void setStarting(double *s);
  virtual double logLikelihood();
  virtual double logHastingsAdj();
  virtual void propose();
  virtual void reject();
  virtual double getCurrentSampleTrans(int indx);
  virtual double* getSamples();
  virtual void setCurrentSample(double*);
  virtual double* getCurrentSample();
  virtual void transSamples();
  virtual void sampleMeansFromStart(int s);
  virtual double* getSampleMeans();
  virtual void show();

 private:
  double *shape;
  double *scale;
  double *tuning;
  double *samples;
  double *sampleMeans;
};
#endif

/************************
     Inverse-Wishart
************************/
#ifndef IWISH_H
#define IWISH_H

class iwish : public vprior
{
 public:
  iwish(string fn, int dimH1, int dimH2, int nP, int nS);
  ~iwish(){/*cout<<"Destructor: iwish"<<endl;*/}

  virtual void setHyper1(double *hyper1);
  virtual void setHyper2(double *hyper2);
  virtual double* getHyper1();
  virtual double* getHyper2();
  virtual void setTuning(double *t);
  virtual void setStarting(double *s);
  virtual double logLikelihood();
  virtual double logHastingsAdj();
  virtual void propose();
  virtual void reject();
  virtual double getCurrentSampleTrans(int indx);
  virtual double* getSamples();
  virtual void setCurrentSample(double*);
  virtual double* getCurrentSample();
  virtual void transSamples();
  virtual void sampleMeansFromStart(int s);
  virtual double* getSampleMeans();
  virtual void show();

  void setA();

 private:
  double df;
  double *S;
  double *tuning;
  double *samples;
  double *A;
  double *Atmp;
  double *sampleMeans;
};
#endif


/************************
         Uniform
************************/
#ifndef UNIFMH_H
#define UNIFMH_H

class unif : public vprior
{
 public:
  unif(string fn, int dimH1, int dimH2, int nP, int nS);
  ~unif(){/*cout<<"Destructor: unif"<<endl;*/}

  virtual void setHyper1(double *hyper1);
  virtual void setHyper2(double *hyper2);
  virtual double* getHyper1();
  virtual double* getHyper2();
  virtual void setTuning(double *t);
  virtual void setStarting(double *s);
  virtual double logLikelihood();
  virtual double logHastingsAdj();
  virtual void propose();
  virtual void reject();
  virtual double getCurrentSampleTrans(int indx);
  virtual double* getSamples();
  virtual void setCurrentSample(double*);
  virtual double* getCurrentSample();
  virtual void transSamples();
  virtual void sampleMeansFromStart(int s);
  virtual double* getSampleMeans();
  virtual void show();

 private:
  double *a;
  double *b;
  double *tuning;
  double *samples;
  double *sampleMeans;
};

#endif


/************************
       Half-Cauchy
************************/
#ifndef HC_H
#define HC_H

class hc : public vprior
{
 public:
  hc(string fn, int dimH1, int dimH2, int nP, int nS);
  ~hc(){/*cout<<"Destructor: hc"<<endl;*/}

  virtual void setHyper1(double *hyper1);
  virtual void setHyper2(double *hyper2);
  virtual double* getHyper1();
  virtual double* getHyper2();
  virtual void setTuning(double *t);
  virtual void setStarting(double *s);
  virtual double logLikelihood();
  virtual double logHastingsAdj();
  virtual void propose();
  virtual void reject();
  virtual double getCurrentSampleTrans(int indx);
  virtual double* getSamples();
  virtual void setCurrentSample(double*);
  virtual double* getCurrentSample();
  virtual void transSamples();
  virtual void sampleMeansFromStart(int s);
  virtual double* getSampleMeans();
  virtual void show();

 private:
  double *a;
  double *tuning;
  double *samples;
  double *sampleMeans;
};

#endif


/************************
   Fixed parameter
************************/
#ifndef FIXEDPAR_H
#define FIXEDPAR_H

class fixedpar : public vprior
{
 public:
  fixedpar(string fn, int dimH1, int dimH2, int nP, int nS);
  ~fixedpar(){/*cout<<"Destructor: fixed"<<endl;*/}

  virtual void setHyper1(double *hyper1);
  virtual void setHyper2(double *hyper2);
  virtual double* getHyper1();
  virtual double* getHyper2();
  virtual void setTuning(double *t);
  virtual void setStarting(double *s);
  virtual double logLikelihood();
  virtual double logHastingsAdj();
  virtual void propose();
  virtual void reject();
  virtual double getCurrentSampleTrans(int indx);
  virtual double* getSamples();
  virtual void setCurrentSample(double*);
  virtual double* getCurrentSample();
  virtual void transSamples();
  virtual void sampleMeansFromStart(int s);
  virtual double* getSampleMeans();
  virtual void show();

 private:
  double *samples;
  double *sampleMeans;
};

#endif

/************************
 Fixed matrix parameter
************************/
#ifndef FIXEDMTRX_H
#define FIXEDMTRX_H

class fixedmtrx : public vprior
{
 public:
  fixedmtrx(string fn, int dimH1, int dimH2, int nP, int nS);
  ~fixedmtrx(){/*cout<<"Destructor: fixedmtrx"<<endl;*/}

  virtual void setHyper1(double *hyper1);
  virtual void setHyper2(double *hyper2);
  virtual double* getHyper1();
  virtual double* getHyper2();
  virtual void setTuning(double *t);
  virtual void setStarting(double *s);
  virtual double logLikelihood();
  virtual double logHastingsAdj();
  virtual void propose();
  virtual void reject();
  virtual double getCurrentSampleTrans(int indx);
  virtual double* getSamples();
  virtual void setCurrentSample(double*);
  virtual double* getCurrentSample();
  virtual void transSamples();
  virtual void sampleMeansFromStart(int s);
  virtual double* getSampleMeans();
  virtual void show();

  void setLowerTriA();

 private:
  double *samples;
  double *A;
  double *sampleMeans;
};
#endif


/************************
   Beta parameter
***********************/
#ifndef BETAPAR_H
#define BETAPAR_H

class betapar : public vprior
{
 public:
  betapar(string fn, int dimH1, int dimH2, int nP, int nS, bool n, bool g);
  ~betapar(){/*cout<<"Destructor: betapar"<<endl;*/}

  virtual void setHyper1(double *hyper1);
  virtual void setHyper2(double *hyper2);
  virtual double* getHyper1();
  virtual double* getHyper2();
  virtual void setTuning(double *t);
  virtual void setStarting(double *s);
  virtual double logLikelihood();
  virtual double logHastingsAdj();
  virtual void propose();//for MH
  virtual void reject();
  virtual double getCurrentSampleTrans(int indx);
  virtual double* getSamples();
  virtual void setCurrentSample(double*);
  virtual double* getCurrentSample();
  virtual void transSamples();
  virtual void sampleMeansFromStart(int s);
  virtual double* getSampleMeans();
  virtual void show();

 private:
  double *samples;
  double *tuning;
  double *mu;
  double *sigma;
  bool normalPrior;
  bool gibbs;
  double *tmpNParVec;
  double *tmpNParVec1;
  double *sampleMeans;

};
#endif
