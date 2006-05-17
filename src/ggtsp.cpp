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
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <math.h>
using namespace std;

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include "vparammh.h"
#include "igmh.h"
#include "unifmh.h"
#include "logunifmh.h"
#include "fixedmh.h"
#include "thetamh.h"
#include "hc.h"
#include "util.h"
#include "covmodel.h"

extern "C" {
  
  SEXP ggtsp(SEXP args)
  {
    SEXP varList, coords, xdims, fixedList, paramControlList, paramControlListNames, tmpSEX1, tmpSEX2, tmpSEX3;
    int i,j, verbose, nProtect= 0;

    /*****************************************************************
       get the args and setup
    ******************************************************************/
    /*********************
       verbose
    **********************/
    if(!isInteger(getListElement(args, "verbose")))
      error("c++ error: verbose must be of type integer");
    
    verbose = INTEGER(getListElement(args, "verbose"))[0];

    /*********************
       DIC
    **********************/
    int dic;
    int dicStart;
    if(!isInteger(getListElement(args, "DIC")))
      error("c++ error: DIC must be of type integer");
    
    dic = INTEGER(getListElement(args, "DIC"))[0];

    if(!isInteger(getListElement(args, "DIC.start")))
      error("c++ error: DIC.start must be of type integer");
    
    dicStart = INTEGER(getListElement(args, "DIC.start"))[0];

    /*********************
       run
    **********************/
    SEXP runList;
    int nSamples = 0;

    PROTECT(runList = getListElement(args, "run")); nProtect++;
    if(!isNewList(runList))
      error("c++ error: run must be a list");
    
    if(!isInteger(getListElement(runList, "n.samples")))
      error("c++ error: n.samples must be an int");
    
    nSamples = INTEGER(getListElement(runList, "n.samples"))[0];
    
    //The starting values of the fixed and var parameters are  
    //stored in the first row of the return sample matrix.  nSamples is used to control the
    //main sampling loop etc, so ++ to accommodate the starting values.
    nSamples++;

    /*********************
       covariance model
    **********************/
    if(!isString(getListElement(args, "cov.model")))
      error("c++ error: cov.model needs to be a string");
    
    string covModel = CHAR(STRING_ELT(getListElement(args, "cov.model"), 0));
    bool onePramPtr = true;

    void (covmodel::*cov1ParamPtr)(double, double *r, double *d, int &length) = NULL; 
    void (covmodel::*cov2ParamPtr)(double, double, double *r, double *d, int &length) = NULL;

    if(covModel == "exponential"){
      cov1ParamPtr = &covmodel::exponential;
    }else if(covModel == "spherical"){
      cov1ParamPtr = &covmodel::spherical;
    }else if(covModel == "gaussian"){
      cov1ParamPtr = &covmodel::gaussian;
    }else if(covModel == "matern"){
      cov2ParamPtr = &covmodel::matern;
      onePramPtr = false;
    }else{
      error("c++ error: cov.model is not correctly specified");
    }
   
    //my covmodel object for calling cov function
    covmodel *covModelObj = new covmodel;

    /*********************
       var
    **********************/
    //get var list
    PROTECT(varList = getListElement(args, "var")); nProtect++;
    if(!isNewList(varList))
      error("c++ error: var must be a list");

    //
    //for the free parameters
    //

    //get the free parameter control list
    PROTECT(paramControlList = getListElement(varList, "var.update.control")); nProtect++;
    if(!isNewList(paramControlList))
      error("c++ error: var.update.control must be a list");

    //make sure there is at least one free parameter
    if(length(paramControlList) == 0)
      error("c++ error: there must be at least one free parameter in the var.update.control");

    //get the formal names from the paramControlList
    PROTECT(paramControlListNames = getAttrib(paramControlList, R_NamesSymbol)); nProtect++;

    //
    //set up the paramvec
    //

    //make vector virtual of parameter pointers
    int nFreeParams = length(paramControlListNames);
    vector<vparammh *>paramVec(nFreeParams);
    string tmpPrior;
    
 
    for(i = 0; i < nFreeParams; i++){
      
      //get each parameter and add to the paramVec
      //tmpSEX1 is the given parameter list element (and is a list itself)
      PROTECT(tmpSEX1 = getListElement(paramControlList, CHAR(STRING_ELT(paramControlListNames, i))));
           
      //tmpSEX2 is the prior object (i.e., list) for the given parameter tmpSEX1
      PROTECT(tmpSEX2 = getListElement(tmpSEX1, "prior"));
      
      tmpPrior = CHAR(STRING_ELT(getListElement(tmpSEX2, "dist"), 0));
      
      if(tmpPrior == "IG"){ //inverse gamma
	paramVec[i] = new igmh(CHAR(STRING_ELT(paramControlListNames, i)), nSamples);
	
	//tmpSEX3 is the list of parameters associated with the specific prior (e.g., scale, shape, a, b, etc.)
 	PROTECT(tmpSEX3 = getListElement(tmpSEX2, "params"));
	paramVec[i]->setPriorParams(0, REAL(getListElement(tmpSEX3, "shape"))[0]);
	paramVec[i]->setPriorParams(1, REAL(getListElement(tmpSEX3, "scale"))[0]);
	
      }else if(tmpPrior == "UNIF"){ //uniform
	
	paramVec[i] = new unifmh(CHAR(STRING_ELT(paramControlListNames, i)), nSamples);
	
	//tmpSEX3 is the list of parameters associated with the specific prior (e.g., scale, shape, a, b, etc.)
 	PROTECT(tmpSEX3 = getListElement(tmpSEX2, "params"));
	paramVec[i]->setPriorParams(0, REAL(getListElement(tmpSEX3, "a"))[0]);
	paramVec[i]->setPriorParams(1, REAL(getListElement(tmpSEX3, "b"))[0]);
	
      }else if(tmpPrior == "LOGUNIF"){ //log uniform
	paramVec[i] = new logunifmh(CHAR(STRING_ELT(paramControlListNames, i)), nSamples);
	
	//tmpSEX3 is the list of parameters associated with the specific prior (e.g., scale, shape, a, b, etc.)
	PROTECT(tmpSEX3 = getListElement(tmpSEX2, "params"));
	paramVec[i]->setPriorParams(0, REAL(getListElement(tmpSEX3, "a"))[0]);
	paramVec[i]->setPriorParams(1, REAL(getListElement(tmpSEX3, "b"))[0]);
	
      }else if(tmpPrior == "HC"){ //half-cauchy
	
	paramVec[i] = new hc(CHAR(STRING_ELT(paramControlListNames, i)), nSamples);
	
	//tmpSEX3 is the list of parameters associated with the specific prior (e.g., scale, shape, a, b, etc.)
	PROTECT(tmpSEX3 = getListElement(tmpSEX2, "params"));
	paramVec[i]->setPriorParams(0, REAL(getListElement(tmpSEX3, "a"))[0]);
	
      }else{
	error("c++ error: parameter type %s not found for parameter %s", tmpPrior.c_str(), CHAR(STRING_ELT(paramControlListNames, i)));
      }
      
      //set the common stuff in the base class (recall even for logged parameters igmh and logunif the virtual class
      //handles logging of the starting value, so its cool to set starting etc for all parameters at once)

      paramVec[i]->setStarting(REAL(getListElement(tmpSEX1, "starting"))[0]);
      paramVec[i]->setTuning(REAL(getListElement(tmpSEX1, "tuning"))[0]);
      paramVec[i]->setSampleOrder(static_cast<int>(REAL(getListElement(tmpSEX1, "sample.order"))[0]));
      UNPROTECT(3); //tmpSEX1, tmpSEX2, tmpSEX3
    }

    //
    //set up the sample order map with key being the order and value portion a vector of vparammh pointers to paramMap->second
    //
    map<int, vector<vparammh*> > sampleOrderMap;
    map<int, vector<vparammh*> >::iterator sampleOrderMapIter;
    vector<vparammh*> tmpVec1(0);//just a temporary vector for init, will be copied into the sampleOrderMap

    //use the fact that maps only have unique keys, to make a unique key list of sample order
    for(i = 0; i < nFreeParams; i++){
      sampleOrderMap.insert(pair<int, vector<vparammh*> >(paramVec[i]->getSampleOrder(), tmpVec1));
    } 

    //now add vparammh* pointers to their respective order vectors
    for(i = 0; i < nFreeParams; i++){
      for(sampleOrderMapIter = sampleOrderMap.begin(); sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++){
	if(paramVec[i]->getSampleOrder() == sampleOrderMapIter->first)
	  sampleOrderMapIter->second.push_back(paramVec[i]);
      }
    }
    
//     //show sample order
//     if(verbose){
//       Rprintf("Metropolis-Hastings sampling order:\n");
//       for(sampleOrderMapIter = sampleOrderMap.begin(), j=1; sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++, j++){
// 	Rprintf("\tSample set %i: ", j);
// 	for(i = 0; i < sampleOrderMapIter->second.size(); i++){
// 	  Rprintf("%s ",sampleOrderMapIter->second[i]->getFormalName().c_str());
// 	}
// 	Rprintf("\n");
//       }
//     }
    
    //
    //for the fixed parameters, so I'm setting these up using the virtual parameter class
    //This works fine but I'm not sure it is the best approach.  So the way it works is the propose() in fixedmh class
    //is called once to fill the entire sample vector with  the starting value (i.e., fixed value), 
    //a bit wasteful but it allows the fixed parameters to be in the same vector of virtual 
    //pointers as the other parameters, and all is well.
    //
    
    //get the fixed parameter control list (recycle paramControlList)
    PROTECT(paramControlList = getListElement(varList, "fixed.var.params")); nProtect++;
    if(!isNewList(paramControlList))
      error("c++ error: fixed.var.params must be a list");

    //make sure there are fixed parameters
    if(length(paramControlList) != 0){
      
      //get their names
      PROTECT(paramControlListNames = getAttrib(paramControlList, R_NamesSymbol)); nProtect++;

      int nFixedParams = length(paramControlListNames);

      //resize the parameter vector to accommodate the fixed parameters
      for(i = 0; i < nFixedParams; i++){
	
	paramVec.resize(paramVec.size()+1);

	paramVec[paramVec.size()-1] = new fixedmh(CHAR(STRING_ELT(paramControlListNames, i)), nSamples);
	
	//get each parameter and add to the previous paramVec
	//tmpSEX1 is the given parameter list element (and is a list itself)
	PROTECT(tmpSEX1 = getListElement(paramControlList, CHAR(STRING_ELT(paramControlListNames, i))));
	
	//only need to set starting for fixed parameters
	paramVec[paramVec.size()-1]->setStarting(REAL(getListElement(tmpSEX1, "fixed"))[0]);
	
	//now just call its proposal to fill the sample vector with fixed values, the 
	//propose values is not used but is needed for the virtual function call so just stick something, e.g.,0.0
	//just make sure the starting values is called before
	paramVec[paramVec.size()-1]->propose(REAL(getListElement(tmpSEX1, "fixed"))[0]);
	
	UNPROTECT(1); //tmpSEX1
      }
      
    }

//     //show the error structure
//     if(verbose){
//     Rprintf("Random effects and associated specifications:\n");
//     for(i = 0; i < paramVec.size(); i++){
//     paramVec[i]->show();
//     }
//     }

    //
    //finally
    //
    //Ok now set parameter specific pointers to the paramVec elements.
    //I want to make all of my correlation functions generic (i.e., they will all take the same parameters 
    //phi and nu), so if matern then nu is in the paramVec as either free or fixed, otherwise nu is null.

    vparammh *phi = NULL, *sigmasq = NULL, *tausq = NULL, *nu = NULL;
    bool useTausq = false;

    for(i = 0; i < paramVec.size(); i++){
      if("phi" == paramVec[i]->getFormalName()){
	phi = paramVec[i];
      }else if("sigmasq" == paramVec[i]->getFormalName()){
	sigmasq = paramVec[i];
      }else if("tausq" == paramVec[i]->getFormalName()){
	tausq  = paramVec[i];
	useTausq = true;
      }else if("nu" == paramVec[i]->getFormalName()){
	nu = paramVec[i];
      }else{
	error("c++ error: somethings wrong, covariance model abbreviation not found, should have been caught on the R side");;
      }
    }

    /*********************
       fixed
    **********************/
    int xnrow, xncol, length, incOne = 1;
    double *x, *y;
    SEXP xColNames, priorMu, priorPrecision;

    //get error list
    PROTECT(fixedList = getListElement(args, "fixed")); nProtect++;
    if(!isNewList(fixedList))
      error("c++ error: fixed must be a list");

    //
    //get the X and Y matrices
    //

    //get the X matrix
    if(!isMatrix(getListElement(fixedList, "X"))){
      error("c++ error: fixed list element 'X' is not a matrix \n");
    }
    
    PROTECT(xdims = getAttrib(getListElement(fixedList, "X"), R_DimSymbol)); nProtect++;      
    xnrow = INTEGER(xdims)[0];
    xncol = INTEGER(xdims)[1];

    x = (double *) R_alloc(xnrow*xncol, sizeof(double));
    length = xnrow*xncol;

    F77_NAME(dcopy)(&length, REAL(getListElement(fixedList, "X")), &incOne, x, &incOne);

    //get the X matrix names for later
    PROTECT(xColNames = VECTOR_ELT(getAttrib(getListElement(fixedList, "X"), R_DimNamesSymbol), 1)); nProtect++;

    //get the Y vector
    if(!isMatrix(getListElement(fixedList, "Y"))){
      error("c++ error: fixed list element 'Y' is not a matrix \n");
    }
 
    y = (double *) R_alloc(xnrow, sizeof(double));

    F77_NAME(dcopy)(&xnrow, REAL(getListElement(fixedList, "Y")), &incOne, y, &incOne);

    //
    //fixed effects samples
    //
    SEXP fixedEffectSamples, samplesColNames, samplesRowNames, dimNames;

    PROTECT(fixedEffectSamples = allocMatrix(REALSXP, xncol, nSamples)); nProtect++;

    //add column names, first the fixed effects then the random effects
    PROTECT(samplesColNames = allocVector(STRSXP, xncol)); nProtect++;
    PROTECT(samplesRowNames = allocVector(STRSXP, nSamples)); nProtect++;
    PROTECT(dimNames = allocVector(VECSXP, 2)); nProtect++;

    //make column names
    for(i = 0; i < xncol; i++)
      SET_VECTOR_ELT(samplesColNames, i, STRING_ELT(xColNames, i));

    SET_VECTOR_ELT(dimNames, 1, samplesRowNames);
    SET_VECTOR_ELT(dimNames, 0, samplesColNames);
    setAttrib(fixedEffectSamples, R_DimNamesSymbol, dimNames);

    //set to zero for fun
    for (i = 0; i < nSamples*xncol; i++)
      REAL(fixedEffectSamples)[i] = 0.0;

    //
    //fixed effects starting values
    //

    //add the starting values directly to the fixedEffectSample matrix
    if(!isVector(getListElement(fixedList, "fixed.effects.starting"))){
      error("c++ error: fixed list element 'fixed.effects.starting' is not a vector\n");
    }

    F77_NAME(dcopy)(&xncol, &REAL(getListElement(fixedList, "fixed.effects.starting"))[0], &incOne, &REAL(fixedEffectSamples)[0], &incOne);
    
    //
    //how to update fixed effects theta and tuning matrix 
    //
    if(!isString(getListElement(fixedList, "update.theta.method")))
      error("c++ error: update.theta.method needs to be a string");
    
    string updateThetaMethod = CHAR(STRING_ELT(getListElement(fixedList, "update.theta.method"), 0));
    bool useGibbsThetaUpdate;
    int thetaTuningLength = xncol*xncol;
    double *thetaTuning;
    
    if(updateThetaMethod == "MH"){
      useGibbsThetaUpdate = false;

      //allocate then copy the tuning matrix
      thetaTuning = (double *) R_alloc(thetaTuningLength, sizeof(double));

      if(!isMatrix(getListElement(fixedList, "tuning.matrix"))){
	error("c++ error: fixed list element 'tuning.matrix' is not a matrix \n");
      }
      
      F77_NAME(dcopy)(&thetaTuningLength, REAL(getListElement(fixedList, "tuning.matrix")), &incOne, thetaTuning, &incOne);
      
      //So we would like to use sample order with the mh option for updating the fixed effects
      //the solution is to include a dummy parameter in the paramVec then subsiquently in the 
      //sampleOrderMap.  This is a bit ad hoc but, but should serve just fine.


      paramVec.resize(paramVec.size()+1);

      paramVec[paramVec.size()-1] = new thetamh("Beta", 0);
      paramVec[paramVec.size()-1]->setSampleOrder(static_cast<int>(REAL(getListElement(fixedList, "sample.order"))[0]));
      
      //As the sampleOrderMap has already been established above, check if the sample order of the theta is any
      //of the pre-existing order and add it to the corresponding second, otherwise add a new key to the 
      //sampleOrderMap

      bool found = false;
      for(sampleOrderMapIter = sampleOrderMap.begin(); sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++){
	if(paramVec[paramVec.size()-1]->getSampleOrder() == sampleOrderMapIter->first){
	  sampleOrderMapIter->second.push_back(paramVec[paramVec.size()-1]);
	  found = true;
	}
      }
      
      if(!found){//if sample order not in there then add it for the theta
	//tmpVec1 is just a temporary vector for init, will be copied into the sampleOrderMap
	sampleOrderMap.insert(pair<int, vector<vparammh*> >(paramVec[paramVec.size()-1]->getSampleOrder(), tmpVec1));
	sampleOrderMap[paramVec[paramVec.size()-1]->getSampleOrder()].push_back(paramVec[paramVec.size()-1]);
      }
      
    }else if(updateThetaMethod == "GIBBS"){
      useGibbsThetaUpdate = true;
    }else{
      error("c++ error: update.theta.method directive in the fixed list must be either MH or GIBBS");
    }

    //show sample order
    if(verbose){
      Rprintf("Metropolis-Hastings sampling order:\n");
      for(sampleOrderMapIter = sampleOrderMap.begin(), j=1; sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++, j++){
	Rprintf("\tSample set %i: ", j);
	for(i = 0; i < sampleOrderMapIter->second.size(); i++){
	  Rprintf("%s ",sampleOrderMapIter->second[i]->getFormalName().c_str());
	}
	Rprintf("\n");
      }
    }

    //
    //use flat or normal prior on the fixed
    //
    if(!isString(getListElement(fixedList, "prior")))
      error("c++ error: theta prior needs to be a string");
    
    string thetaPrior = CHAR(STRING_ELT(getListElement(fixedList, "prior"), 0));

    double *thetaPriorMu;
    double *thetaPriorV;
    int thetaPriorVLength = xncol*xncol;

    if(thetaPrior == "NORMAL"){

      //for the mu
      if(!isMatrix(getListElement(fixedList, "prior.mu"))){
	error("c++ error: prior.mu is not a matrix \n");
      }

      thetaPriorMu = (double *) R_alloc(xncol, sizeof(double));
      F77_NAME(dcopy)(&xncol, REAL(getListElement(fixedList, "prior.mu")), &incOne, thetaPriorMu, &incOne);


      //for the precision
      if(!isMatrix(getListElement(fixedList, "prior.precision"))){
	error("c++ error: prior.precision is not a matrix \n");
      }

      thetaPriorV = (double *) R_alloc(thetaPriorVLength, sizeof(double));
      F77_NAME(dcopy)(&thetaPriorVLength, REAL(getListElement(fixedList, "prior.precision")), &incOne, thetaPriorV, &incOne);
    }

    /*********************
       distance matrix
    **********************/
    //get error list
    PROTECT(coords = getListElement(args, "coords")); nProtect++;
    if(!isMatrix(coords))
      error("c++ error: coords is not a matrix \n");
    
    
    PROTECT(xdims = getAttrib(coords, R_DimSymbol)); nProtect++;      
    if(INTEGER(xdims)[0] != xnrow)
      error("c++ error: rows of coords != nxrows \n");
    if(INTEGER(xdims)[1] != 2)
      error("c++ error: cols of coords != 2 \n");
      
    int D_length = xnrow*xnrow;
    double* D = (double *) R_alloc(D_length, sizeof(double));
    
    double dmax = 0;
    for(i = 0; i < xnrow; i++){
      for(j = 0; j < xnrow; j++){
	D[i*xnrow+j] = sqrt(pow(REAL(coords)[i]-REAL(coords)[j],2)+pow(REAL(coords)[xnrow+i]-REAL(coords)[xnrow+j],2));
	if(D[i*xnrow+j] > dmax)
	  dmax = D[i*xnrow+j];
      }
    }
    //cout << "max distance: " << dmax << endl << endl;
    //showMatrix(D, xnrow, xnrow);

    /*********************
       DIC set-up
    **********************/
    int nDIC = nSamples - dicStart;
    double *DD;
    if(dic){
      DD = (double *) R_alloc(nDIC, sizeof(double));
    }

    /*****************************************************************
       sampling
    ******************************************************************/
    if(verbose){
      Rprintf("\n------------------------------------------------------------------------\n");
      Rprintf("\t\t\tStart sampling");
      Rprintf("\n------------------------------------------------------------------------\n");
      Rprintf("Samples collected: ");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    int s=0, info=0, status=0, rtnStatus=0, thetaAccept = 0;
    double logPostCurrent = 0;
    double logPostCand = 0;
    double phiTmp, sigmasqTmp, tausqTmp;
    bool accept = true;
    bool first = true;
    const char lower = 'L';
    const char upper = 'U';
    const char ntran = 'N';
    const char ytran = 'Y';
    const char rside = 'R';
    const char lside = 'L';
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    double logDetCov = 0;
    double logDetCovCurrent = 0;
    double r = 0;
    double rUnif = 0;
    double *R; R = (double *) R_alloc(D_length, sizeof(double));
    double *RInvCurrent; RInvCurrent = (double *) R_alloc(D_length, sizeof(double));
    double *tmpXRow; tmpXRow = (double *) R_alloc(xnrow, sizeof(double));
    double *tmpXRow1; tmpXRow1 = (double *) R_alloc(xnrow, sizeof(double));
    double *tmpXCol; tmpXCol = (double *) R_alloc(xncol, sizeof(double));
    double *tmpXCol1; tmpXCol1 = (double *) R_alloc(xncol, sizeof(double));
    double *tmpXRowCol; tmpXRowCol = (double *) R_alloc(xnrow*xncol, sizeof(double));
    double *tmpXColCol; tmpXColCol = (double *) R_alloc(xncol*xncol, sizeof(double));
    bool redraw = true;
    double proposal = 0;
    double junk;
    int job = 01;
    int dicIndx = 0;
    bool candTheta = false;

   GetRNGstate();
    for(s = 0; s < nSamples-1; s++){//off-by-one because of the starting values
      
      //for each subset of parameters 
      for(sampleOrderMapIter = sampleOrderMap.begin(); sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++){

	candTheta = false;//should the a proposal theta be considered in this mh iter

	//get proposal for each parameter in the subset
	for(i = 0; i < sampleOrderMapIter->second.size(); i++){

	  redraw = true;
	  while(redraw){
	    if(sampleOrderMapIter->second[i]->getFormalName() != "Beta"){
	      
	      proposal = rnorm(sampleOrderMapIter->second[i]->getCurrentSampleNoTrans(), 
			       sampleOrderMapIter->second[i]->getTuning());
	      
	      //returns falus if proposal is in the prior's support, currently for the uniforms
	      redraw = sampleOrderMapIter->second[i]->propose(proposal);
	      
	    }else{//for mh theta
	      mvrnorm(&REAL(fixedEffectSamples)[(s+1)*xncol], &REAL(fixedEffectSamples)[s*xncol], thetaTuning, xncol);
	      candTheta = true;
	      redraw = false;
	    }
	    R_CheckUserInterrupt();
	  }
	}

 	//make the correlation matrix
 	if(onePramPtr)
 	  (covModelObj->*cov1ParamPtr)(phi->getCurrentSampleTrans(), R, D, D_length);
 	else //i.e., 2 parameter matern
 	  (covModelObj->*cov2ParamPtr)(phi->getCurrentSampleTrans(), nu->getCurrentSampleTrans(), R, D, D_length);
	
 	//scale correlation matrix with sigmasq
 	sigmasqTmp = sigmasq->getCurrentSampleTrans();
 	F77_NAME(dscal)(&D_length, &sigmasqTmp, R, &incOne);
	
 	//use tausq
 	if(useTausq)
 	  for(i = 0; i < xnrow; i++) R[i*xnrow+i] = R[i*xnrow+i]+tausq->getCurrentSampleTrans();
	
	//get the sum of the logged posterior with their current parameter values (some of which are proposals)
	logPostCand = 0;
	for(i = 0; i < paramVec.size(); i++){
	  logPostCand += paramVec[i]->logPrior();
	}

	//chol decom
	F77_NAME(dpofa)(R, &xnrow, &xnrow, &info);
	if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}

	//get logDet
	logDetCov = 0;

	for(i = 0; i < xnrow; i++)
	  logDetCov += log(R[i*xnrow+i]);
	logDetCov = 2*logDetCov;
	logPostCand += -0.5*logDetCov;

	//finish the invert
	F77_NAME(dpodi)(R, &xnrow, &xnrow, &junk, &job);
	if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}

	
	if(candTheta){//using proposed theta
	  //tmp = y - x*theta
	  F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &one, x, &xnrow, &REAL(fixedEffectSamples)[(s+1)*xncol], &incOne, &zero, tmpXRow, &incOne);
	  F77_NAME(daxpy)(&xnrow, &negOne, y, &incOne, tmpXRow, &incOne);
	  
	  //(-1/2) * tmp` * R^{-1} * tmp
	  F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
	  logPostCand += -0.5*F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne);
	  
	  if(thetaPrior == "NORMAL"){
	    F77_NAME(dcopy)(&xncol, thetaPriorMu, &incOne, tmpXCol, &incOne);
	    F77_NAME(daxpy)(&xncol, &negOne, &REAL(fixedEffectSamples)[(s+1)*xncol], &incOne, tmpXCol, &incOne);
	    F77_NAME(dgemv)(&ntran, &xncol, &xncol, &one, thetaPriorV, &xncol, tmpXCol, &incOne, &zero, tmpXCol1, &incOne);
	    logPostCand += -0.5*F77_NAME(ddot)(&xncol, tmpXCol, &incOne, tmpXCol1, &incOne);
	  }
	}
	else{//theta is not being considered so use currently accepted theta
	  //tmp = y - x*theta
	  F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &one, x, &xnrow, &REAL(fixedEffectSamples)[s*xncol], &incOne, &zero, tmpXRow, &incOne);
	  F77_NAME(daxpy)(&xnrow, &negOne, y, &incOne, tmpXRow, &incOne);
	  
	  //(-1/2) * tmp` * R^{-1} * tmp
	  F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
	  logPostCand += -0.5*F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne);
	  
	  if(thetaPrior == "NORMAL"){
	    F77_NAME(dcopy)(&xncol, thetaPriorMu, &incOne, tmpXCol, &incOne);
	    F77_NAME(daxpy)(&xncol, &negOne, &REAL(fixedEffectSamples)[s*xncol], &incOne, tmpXCol, &incOne);
	    F77_NAME(dgemv)(&ntran, &xncol, &xncol, &one, thetaPriorV, &xncol, tmpXCol, &incOne, &zero, tmpXCol1, &incOne);
	    logPostCand += -0.5*F77_NAME(ddot)(&xncol, tmpXCol, &incOne, tmpXCol1, &incOne);
	  } 
	}

	if(first){
	  //if this is the very first sample loop
	  logPostCurrent = logPostCand;
	  F77_NAME(dcopy)(&D_length, R, &incOne, RInvCurrent, &incOne);
 	  logDetCovCurrent = logDetCov;
	  first = false;
	}
	
	r = logPostCand - logPostCurrent;
	r = exp(r);

	rUnif = runif(0.0,1.0);
	
	if(r >= 1){
	  logPostCurrent = logPostCand;
	  F77_NAME(dcopy)(&D_length, R, &incOne, RInvCurrent, &incOne);
 	  logDetCovCurrent = logDetCov;

	  if(candTheta)
	    thetaAccept++;
	}
	else if(rUnif < r){
	  logPostCurrent = logPostCand;
	  F77_NAME(dcopy)(&D_length, R, &incOne, RInvCurrent, &incOne);
 	  logDetCovCurrent = logDetCov;

	  if(candTheta)
	    thetaAccept++;
	}
	else{//reject
	  if(candTheta){//if theta is being considered then reject accordingly
	    
	    for(i = 0; i < sampleOrderMapIter->second.size(); i++){
	      if(sampleOrderMapIter->second[i]->getFormalName() != "Beta"){
		sampleOrderMapIter->second[i]->rejectProposal();
	      }else{
		F77_NAME(dcopy)(&xncol, &REAL(fixedEffectSamples)[s*xncol], &incOne, &REAL(fixedEffectSamples)[(s+1)*xncol], &incOne); 
	      }
	    }

	  }else{//if theta is not being considered then simply reject all

	    for(i = 0; i < sampleOrderMapIter->second.size(); i++)
	      sampleOrderMapIter->second[i]->rejectProposal();

	  }
 	}
      }
      
      //update theta if using gibbs
      if(useGibbsThetaUpdate){
	updateThetaGibbs(x, y, xnrow, xncol, &REAL(fixedEffectSamples)[(s+1)*xncol], RInvCurrent, tmpXRowCol, tmpXColCol, tmpXRow, tmpXCol, tmpXCol1, thetaPrior, thetaPriorMu, thetaPriorV);
      }

      //
      //DIC
      //
      if(dic && s+1 >= dicStart){
	//tmp = y - x*theta
	F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &one, x, &xnrow, &REAL(fixedEffectSamples)[(s+1)*xncol], &incOne, &zero, tmpXRow, &incOne);
	F77_NAME(daxpy)(&xnrow, &negOne, y, &incOne, tmpXRow, &incOne);
	
	//(-1/2) * tmp` * R^{-1} * tmp
	F77_NAME(dsymv)(&upper, &xnrow, &one, RInvCurrent, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
	DD[dicIndx] = logDetCovCurrent+F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne); //i.e., DD = -2.0*(-0.5*logDetCov - 0.5*dotResult);
	dicIndx++;
      }

      //report
      if(verbose){
	if(status == 100){
	  Rprintf("%i...", s);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
	if(rtnStatus == 800){
	  Rprintf("\n\t");
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  rtnStatus = 0;
	}
	rtnStatus++;
	status++;
      }
      
      R_CheckUserInterrupt();
    }
    
    PutRNGstate(); 
    
    if(verbose){
      Rprintf("\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    /*********************
        finish DIC
    **********************/
    double DDBar, DDBarOmega, pD, DIC;
    DDBar = DDBarOmega = pD = DIC = 0.0;

    if(dic){
      for(i = 0; i < nDIC; i++){
	DDBar += DD[i];
      }
      DDBar = DDBar/nDIC;
    }
  
    //now for DDBarOmega

    //get the means of the params
    double *thetaMean; 
    double sigmasqMean, tausqMean, phiMean, nuMean;
    sigmasqMean = tausqMean = phiMean = nuMean = 0.0;

    if(dic){
      thetaMean = (double *) R_alloc(xncol, sizeof(double));
      
      //zero
      for(i = 0; i < xncol; i++){thetaMean[i] = 0.0;}
      
      for(i = 0; i < xncol; i++){
	for(s = dicStart; s < nSamples; s++){
	  thetaMean[i] += REAL(fixedEffectSamples)[s*xncol+i];
	}
	thetaMean[i] = thetaMean[i]/nDIC;
      }
   
      //get the means of the variance parameters
      sigmasqMean = sigmasq->getSampleTransMean(dicStart);
      if(useTausq)
	tausqMean = tausq->getSampleTransMean(dicStart);
      if(onePramPtr){
	phiMean = phi->getSampleTransMean(dicStart);
      }else{ //i.e., 2 parameter matern
	phiMean = phi->getSampleTransMean(dicStart);
	nuMean = nu->getSampleTransMean(dicStart);
      }

      //now calculate DIC stuff for the sample means
      
      //make the correlation matrix
      if(onePramPtr)
	(covModelObj->*cov1ParamPtr)(phiMean, R, D, D_length);
      else //i.e., 2 parameter matern
	(covModelObj->*cov2ParamPtr)(phiMean, nuMean, R, D, D_length);
      
      //scale correlation matrix with sigmasq
      F77_NAME(dscal)(&D_length, &sigmasqMean, R, &incOne);
      
      //use tausq
      if(useTausq)
	for(i = 0; i < xnrow; i++) R[i*xnrow+i] = R[i*xnrow+i]+tausqMean;
	
      //chol decom
      F77_NAME(dpofa)(R, &xnrow, &xnrow, &info);
      if(info != 0){error("c++ error: Cholesky failed (2), see ggt.sp documentation \n");}

      //get logDet
      logDetCov = 0;
      for(i = 0; i < xnrow; i++)
	logDetCov += log(R[i*xnrow+i]);
      logDetCov = 2*logDetCov;
      
      //finish the invert
      F77_NAME(dpodi)(R, &xnrow, &xnrow, &junk, &job);
      if(info != 0){error("c++ error: Cholesky inverse failed (2), see ggt.sp documentation \n");}

      //NaN check 
      for(i = 0; i < xnrow*xnrow; i++){
	if(isnan(R[i]))
	  error("NaN discovered in the covariance matrix (2), see ggt.sp documentation");
      }

      //tmp = y - x*theta
      F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &one, x, &xnrow, &thetaMean[0], &incOne, &zero, tmpXRow, &incOne);
      F77_NAME(daxpy)(&xnrow, &negOne, y, &incOne, tmpXRow, &incOne);
      
      //(-1/2) * tmp` * R^{-1} * tmp
      F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
      DDBarOmega = logDetCov + F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne); //i.e., DDBarOmega = -2.0*(-0.5*logDetCov - 0.5*dotResult);
    }

    /*********************
     make return object
    **********************/
    //get the parameter samples out of their respective vectors and place then in a single matrix
    //then add column headers
    SEXP varSamples, varSamplesColNames, varSamplesRowNames, varDimNames, dicResults, dicResultsRowNames, dicResultsColNames, dicResultDimName;
    
    //remove theta from the paramVec if needed
    for(i = 0; i < paramVec.size(); i++){
      if(paramVec[i]->getFormalName() == "Beta"){
	paramVec.erase(paramVec.begin()+i);
	break;
      }
    }
    
    PROTECT(varSamples = allocMatrix(REALSXP, nSamples, paramVec.size())); nProtect++;

    for(i = 0; i < paramVec.size(); i++){
      for(j = 0; j < nSamples; j++){
	REAL(varSamples)[(i*nSamples)+j] = paramVec[i]->getSampleTrans(j);
      }
    }
    
    PROTECT(varSamplesColNames = allocVector(STRSXP, paramVec.size())); nProtect++;
    PROTECT(varSamplesRowNames = allocVector(STRSXP, nSamples)); nProtect++;
    PROTECT(varDimNames = allocVector(VECSXP, 2)); nProtect++;

    for(i = 0; i < paramVec.size(); i++)
      SET_VECTOR_ELT(varSamplesColNames, i, mkChar(paramVec[i]->getFormalName().c_str()));

    SET_VECTOR_ELT(varDimNames, 0, varSamplesRowNames); //should be default null
    SET_VECTOR_ELT(varDimNames, 1, varSamplesColNames);
    setAttrib(varSamples, R_DimNamesSymbol, varDimNames);

    //put acceptance in a matrix with parameters as row names
    SEXP acceptance, acceptColNames, acceptRowNames, acceptDimNames;

    PROTECT(acceptColNames = allocVector(STRSXP, 1)); nProtect++;
    SET_VECTOR_ELT(acceptColNames, 0, mkChar("acceptance.rate"));

    if(useGibbsThetaUpdate){
      PROTECT(acceptance = allocMatrix(REALSXP, paramVec.size(), 1)); nProtect++; 
 
      for(i = 0; i < paramVec.size(); i++)
	REAL(acceptance)[i] = paramVec[i]->getAcceptRate();
      
      PROTECT(acceptRowNames = allocVector(STRSXP, paramVec.size())); nProtect++;

      for(i = 0; i < paramVec.size(); i++)
	SET_VECTOR_ELT(acceptRowNames, i, mkChar(paramVec[i]->getFormalName().c_str()));
    }else{
      PROTECT(acceptance = allocMatrix(REALSXP, paramVec.size()+1, 1)); nProtect++;  

      for(i = 0; i < paramVec.size(); i++)
	REAL(acceptance)[i] = paramVec[i]->getAcceptRate();
      REAL(acceptance)[paramVec.size()] = 100.0*thetaAccept/nSamples;

      PROTECT(acceptRowNames = allocVector(STRSXP, paramVec.size()+1)); nProtect++;

      for(i = 0; i < paramVec.size(); i++)
	SET_VECTOR_ELT(acceptRowNames, i, mkChar(paramVec[i]->getFormalName().c_str()));
      SET_VECTOR_ELT(acceptRowNames, paramVec.size(), mkChar("Beta"));
    }

    PROTECT(acceptDimNames = allocVector(VECSXP, 2)); nProtect++;
    SET_VECTOR_ELT(acceptDimNames, 0, acceptRowNames);
    SET_VECTOR_ELT(acceptDimNames, 1, acceptColNames);
    setAttrib(acceptance, R_DimNamesSymbol, acceptDimNames);

    //make DIC return matrix
    if(dic){
      PROTECT(dicResults = allocMatrix(REALSXP, 4, 1)); nProtect++; //for Dbar, DbarOmega, pD, and DIC
      REAL(dicResults)[0] = DDBar;
      REAL(dicResults)[1] = DDBarOmega;
      REAL(dicResults)[2] = DDBar - DDBarOmega;
      REAL(dicResults)[3] = DDBar + DDBar - DDBarOmega;   
    
      PROTECT(dicResultsRowNames = allocVector(STRSXP, 4)); nProtect++;
      SET_VECTOR_ELT(dicResultsRowNames, 0, mkChar("bar.D"));
      SET_VECTOR_ELT(dicResultsRowNames, 1, mkChar("D.bar.Theta"));
      SET_VECTOR_ELT(dicResultsRowNames, 2, mkChar("pD"));
      SET_VECTOR_ELT(dicResultsRowNames, 3, mkChar("DIC"));

      PROTECT(dicResultsColNames = allocVector(STRSXP, 1)); nProtect++;
      SET_VECTOR_ELT(dicResultsColNames, 0, mkChar("value"));

      PROTECT(dicResultDimName = allocVector(VECSXP, 2)); nProtect++;
      SET_VECTOR_ELT(dicResultDimName, 0, dicResultsRowNames);
      SET_VECTOR_ELT(dicResultDimName, 1, dicResultsColNames);
      setAttrib(dicResults, R_DimNamesSymbol, dicResultDimName);
    }

    //make the result list object
    int nResultListObjs = 3;
    
    if(dic){
      nResultListObjs = 4;
    }

    SEXP result, resultNames;
    
    PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultNames = allocVector(STRSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result, 0, fixedEffectSamples);
    SET_VECTOR_ELT(result, 1, varSamples);
    SET_VECTOR_ELT(result, 2, acceptance);
    
    SET_STRING_ELT(resultNames, 0, mkChar("fixedEffectSamples"));
    SET_STRING_ELT(resultNames, 1, mkChar("varParameterSamples"));
    SET_STRING_ELT(resultNames, 2, mkChar("acceptance"));

    if(dic){
      SET_VECTOR_ELT(result, 3, dicResults);
      SET_STRING_ELT(resultNames, 3, mkChar("DIC"));
    }

    namesgets(result, resultNames);
    
    /*********************
     clean-up and return
    **********************/

    delete covModelObj;

    //first explicitly delete the sample vectors and then the param class in the paramVec
    for(i = 0; i < paramVec.size(); i++) delete paramVec[i];
    paramVec.clear();
    
    //no need to delete the objects again because the sampleOrderMap just has pointers to the paramVec
    for(sampleOrderMapIter = sampleOrderMap.begin(); sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++){
      sampleOrderMapIter->second.clear();
    }
    sampleOrderMap.clear();
   
    //unprotect
    UNPROTECT(nProtect);
    
    return(result);
  }
}
  
