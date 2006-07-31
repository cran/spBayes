#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "vprior.h"
#include "priors.h"
#include "covmodel.h"
#include "util.h"

extern "C" {
  
  SEXP ggtSp(SEXP args){
    int i,j,k,l,nProtect= 0;


    /*****************************************
      run control etc. and some basic checks
    *****************************************/

    if(!isInteger(getListElement(args, "n.samples")))
      error("c++ error: n.samples must be of type integer");

    int nSamples = INTEGER(getListElement(args, "n.samples"))[0];

    if(!isInteger(getListElement(args, "m")))
      error("c++ error: m must be of type integer");
    
    int m = INTEGER(getListElement(args, "m"))[0];
    
    if(!isInteger(getListElement(args, "verbose")))
      error("c++ error: verbose must be of type integer");

    if(!isInteger(getListElement(args, "linpack")))
      error("c++ error: linpack must be of type integer");

    bool linpack = INTEGER(getListElement(args, "linpack"))[0];
    
    bool verbose = INTEGER(getListElement(args, "verbose"))[0];

    bool spEffects = INTEGER(getListElement(args, "sp.effects"))[0];

    bool DIC = INTEGER(getListElement(args, "DIC"))[0];

    int DICStart = INTEGER(getListElement(args, "DIC.start"))[0];


    /*****************************************
          Get the covariance model
    *****************************************/
    if(!isString(getListElement(args, "cov.model")))
      error("c++ error: cov.model needs to be a string");
    
    string covModel = CHAR(STRING_ELT(getListElement(args, "cov.model"), 0));
    
    /*****************************************
           Parse variance parameters
    *****************************************/
    SEXP VarList;
    PROTECT(VarList = getListElement(args, "var.control")); nProtect++;
    bool noPsi = INTEGER(getListElement(VarList, "no.Psi"))[0];

    /*********************
       Parse the K list
    **********************/
    SEXP AList, APriorList, APrior, APriorHyper;
    string PriorDist;
    int ACase = 0;
    vector<vprior*> AParams(0);
 
    PROTECT(AList = getListElement(VarList, "K")); nProtect++;
    if(!isNewList(AList))
      error("c++ error: var must be a list");
    
    ACase = INTEGER(getListElement(AList, "case"))[0];
    
    if(ACase == 1){//all diag elements have share the same parameter
      
      PROTECT(APrior = getListElement(AList, "prior")); 
      PriorDist = CHAR(STRING_ELT(getListElement(APrior, "dist"), 0));
      PROTECT(APriorHyper = getListElement(APrior, "params"));
      
      if(PriorDist == "IG"){
	AParams.push_back(new ig("K", 1, 1, 1, nSamples));
	AParams[0]->setHyper1(&REAL(getListElement(APriorHyper, "shape"))[0]);
	AParams[0]->setHyper2(&REAL(getListElement(APriorHyper, "scale"))[0]);
	AParams[0]->setSampleOrder(INTEGER(getListElement(AList, "sample.order"))[0]);
      }else if(PriorDist == "HC"){
	AParams.push_back(new hc("K", 1, 0, 1, nSamples));
	AParams[0]->setHyper1(&REAL(getListElement(APriorHyper, "a"))[0]);
	AParams[0]->setSampleOrder(INTEGER(getListElement(AList, "sample.order"))[0]);
      }else if(PriorDist == "LOGUNIF"){
	AParams.push_back(new logunif("K", 1, 1, 1, nSamples));
	AParams[0]->setHyper1(&REAL(getListElement(APriorHyper, "a"))[0]);
	AParams[0]->setHyper2(&REAL(getListElement(APriorHyper, "b"))[0]);
	AParams[0]->setSampleOrder(INTEGER(getListElement(AList, "sample.order"))[0]);
      }else if(PriorDist == "FIXED"){
	AParams.push_back(new fixedpar("K", 0, 0, 1, nSamples));
	AParams[0]->setSampleOrder(-1);
      }else{
	error("c++ error: prior misspecified for parameter A, case 2");
      }
      //attributes common to all priors (ignored for fixed)
      AParams[0]->setTuning(&REAL(getListElement(AList, "tuning"))[0]);
      AParams[0]->setStarting(&REAL(getListElement(AList, "starting"))[0]);	
      UNPROTECT(2);
     
    }else if(ACase == 2){//each diag element has its own parameter
  
      //APriorList is a list of length m in this case
      PROTECT(APriorList = getListElement(AList, "prior")); nProtect++; 

      for(i = 0; i < length(APriorList); i++){
	
	//get the next prior objects
	PROTECT(APrior = VECTOR_ELT(APriorList, i));
	PriorDist = CHAR(STRING_ELT(getListElement(APrior, "dist"), 0));
	PROTECT(APriorHyper = getListElement(APrior, "params"));
	
	if(PriorDist == "IG"){
	  AParams.push_back(new ig("K", 1, 1, 1, nSamples));
	  AParams[i]->setHyper1(&REAL(getListElement(APriorHyper, "shape"))[0]);
	  AParams[i]->setHyper2(&REAL(getListElement(APriorHyper, "scale"))[0]);
	  AParams[i]->setSampleOrder(INTEGER(getListElement(AList, "sample.order"))[i]);
	}else if(PriorDist == "HC"){
	  AParams.push_back(new hc("K", 1, 0, 1, nSamples));
	  AParams[i]->setHyper1(&REAL(getListElement(APriorHyper, "a"))[0]);
	  AParams[i]->setSampleOrder(INTEGER(getListElement(AList, "sample.order"))[i]);
	}else if(PriorDist == "LOGUNIF"){
	  AParams.push_back(new logunif("K", 1, 1, 1, nSamples));
	  AParams[i]->setHyper1(&REAL(getListElement(APriorHyper, "a"))[0]);
	  AParams[i]->setHyper2(&REAL(getListElement(APriorHyper, "b"))[0]);
	  AParams[i]->setSampleOrder(INTEGER(getListElement(AList, "sample.order"))[i]);
	}else if(PriorDist == "FIXED"){
	  AParams.push_back(new fixedpar("K", 0, 0, 1, nSamples));
	  AParams[i]->setSampleOrder(-1);
	}else{
	  error("c++ error: prior misspecified for parameter A, case 2");
	}
	//attributes common to all priors (ignored for fixed)
	AParams[i]->setTuning(&REAL(getListElement(AList, "tuning"))[i]);
	AParams[i]->setStarting(&REAL(getListElement(AList, "starting"))[i]);	
	AParams[i]->setSubParIndx(i);
	UNPROTECT(2);
      }
    }else if(ACase == 3){//full covariance matrix (A`A = K)
      if(INTEGER(getListElement(AList, "fixed"))[0]){
	AParams.push_back(new fixedmtrx("K", 1, m, (m*m-m)/2+m, nSamples));
	AParams[0]->setStarting(REAL(getListElement(AList, "starting")));
	AParams[0]->setSampleOrder(-1);// this is the default in vprior
      }else{
	PROTECT(APrior = getListElement(AList, "prior")); 
	PriorDist = CHAR(STRING_ELT(getListElement(APrior, "dist"), 0));
	PROTECT(APriorHyper = getListElement(APrior, "params"));
	
	if(PriorDist != "IWISH")
	  error("c++ error: prior misspecified for parameter A, case 3");
	
	AParams.push_back(new iwish("K", 1, m, (m*m-m)/2+m, nSamples));
	AParams[0]->setHyper1(REAL(getListElement(APriorHyper, "df")));
	AParams[0]->setHyper2(REAL(getListElement(APriorHyper, "S")));
	AParams[0]->setTuning(REAL(getListElement(AList, "tuning")));
	AParams[0]->setStarting(REAL(getListElement(AList, "starting")));
	AParams[0]->setSampleOrder(INTEGER(getListElement(AList, "sample.order"))[0]);
	UNPROTECT(2);
      }
    }else{
      error("c++ error: ACase is misspecified");
    }
    
    
    /*********************
     Parse the Psi list
    **********************/
    vector<vprior*> PsiParams(0);
    int PsiCase = 0;

    if(!noPsi){
      
      SEXP PsiList, PsiPriorList, PsiPrior, PsiPriorHyper;
      
      PROTECT(PsiList = getListElement(VarList, "Psi")); nProtect++;
      if(!isNewList(PsiList))
	error("c++ error: var must be a list");
      
      PsiCase = INTEGER(getListElement(PsiList, "case"))[0];
      
      if(PsiCase == 1){//all diag elements have share the same parameter
	
	PROTECT(PsiPrior = getListElement(PsiList, "prior")); 
	PriorDist = CHAR(STRING_ELT(getListElement(PsiPrior, "dist"), 0));
	PROTECT(PsiPriorHyper = getListElement(PsiPrior, "params"));
	
	if(PriorDist == "IG"){
	  PsiParams.push_back(new ig("Psi", 1, 1, 1, nSamples));
	  PsiParams[0]->setHyper1(&REAL(getListElement(PsiPriorHyper, "shape"))[0]);
	  PsiParams[0]->setHyper2(&REAL(getListElement(PsiPriorHyper, "scale"))[0]);
	  PsiParams[0]->setSampleOrder(INTEGER(getListElement(PsiList, "sample.order"))[0]);
	}else if(PriorDist == "HC"){
	  PsiParams.push_back(new hc("Psi", 1, 0, 1, nSamples));
	  PsiParams[0]->setHyper1(&REAL(getListElement(PsiPriorHyper, "a"))[0]);
	  PsiParams[0]->setSampleOrder(INTEGER(getListElement(PsiList, "sample.order"))[0]);
	}else if(PriorDist == "LOGUNIF"){
	  PsiParams.push_back(new logunif("Psi", 1, 1, 1, nSamples));
	  PsiParams[0]->setHyper1(&REAL(getListElement(PsiPriorHyper, "a"))[0]);
	  PsiParams[0]->setHyper2(&REAL(getListElement(PsiPriorHyper, "b"))[0]);
	  PsiParams[0]->setSampleOrder(INTEGER(getListElement(PsiList, "sample.order"))[0]);
	}else if(PriorDist == "FIXED"){
	  PsiParams.push_back(new fixedpar("Psi", 0, 0, 1, nSamples));
	  PsiParams[0]->setSampleOrder(-1);
	}else{
	  error("c++ error: prior misspecified for parameter Psi, case 2");
	}
	//attributes common to all priors (ignored for fixed)
	PsiParams[0]->setTuning(&REAL(getListElement(PsiList, "tuning"))[0]);
	PsiParams[0]->setStarting(&REAL(getListElement(PsiList, "starting"))[0]);	
	UNPROTECT(2);
	
      }else if(PsiCase == 2){//each diag element has its own parameter
	
	//PsiPriorList is a list of length m in this case
	PROTECT(PsiPriorList = getListElement(PsiList, "prior")); nProtect++; 
	
	for(i = 0; i < length(PsiPriorList); i++){
	  
	  //get the next prior objects
	  PROTECT(PsiPrior = VECTOR_ELT(PsiPriorList, i));
	  PriorDist = CHAR(STRING_ELT(getListElement(PsiPrior, "dist"), 0));
	  PROTECT(PsiPriorHyper = getListElement(PsiPrior, "params"));
	  
	  if(PriorDist == "IG"){
	    PsiParams.push_back(new ig("Psi", 1, 1, 1, nSamples));
	    PsiParams[i]->setHyper1(&REAL(getListElement(PsiPriorHyper, "shape"))[0]);
	    PsiParams[i]->setHyper2(&REAL(getListElement(PsiPriorHyper, "scale"))[0]);
	    PsiParams[i]->setSampleOrder(INTEGER(getListElement(PsiList, "sample.order"))[i]);
	  }else if(PriorDist == "HC"){
	    PsiParams.push_back(new hc("Psi", 1, 0, 1, nSamples));
	    PsiParams[i]->setHyper1(&REAL(getListElement(PsiPriorHyper, "a"))[0]);
	    PsiParams[i]->setSampleOrder(INTEGER(getListElement(PsiList, "sample.order"))[i]);
	  }else if(PriorDist == "LOGUNIF"){
	    PsiParams.push_back(new logunif("Psi", 1, 1, 1, nSamples));
	    PsiParams[i]->setHyper1(&REAL(getListElement(PsiPriorHyper, "a"))[0]);
	    PsiParams[i]->setHyper2(&REAL(getListElement(PsiPriorHyper, "b"))[0]);
	    PsiParams[i]->setSampleOrder(INTEGER(getListElement(PsiList, "sample.order"))[i]);
	  }else if(PriorDist == "FIXED"){
	    PsiParams.push_back(new fixedpar("Psi", 0, 0, 1, nSamples));
	    PsiParams[i]->setSampleOrder(-1);
	  }else{
	    error("c++ error: prior misspecified for parameter Psi, case 2");
	  }
	  //attributes common to all priors (ignored for fixed)
	  PsiParams[i]->setTuning(&REAL(getListElement(PsiList, "tuning"))[i]);
	  PsiParams[i]->setStarting(&REAL(getListElement(PsiList, "starting"))[i]);
	  PsiParams[i]->setSubParIndx(i);
	  UNPROTECT(2);
	}
      }else if(PsiCase == 3){//full covariance matrix (A`A = K)
	if(INTEGER(getListElement(PsiList, "fixed"))[0]){
	  PsiParams.push_back(new fixedmtrx("Psi", 1, m, (m*m-m)/2+m, nSamples));
	  PsiParams[0]->setStarting(REAL(getListElement(PsiList, "starting")));
	  PsiParams[0]->setSampleOrder(-1);// this is the default in vprior
	}else{
	  PROTECT(PsiPrior = getListElement(PsiList, "prior")); 
	  PriorDist = CHAR(STRING_ELT(getListElement(PsiPrior, "dist"), 0));
	  PROTECT(PsiPriorHyper = getListElement(PsiPrior, "params"));
	  
	  if(PriorDist != "IWISH")
	    error("c++ error: prior misspecified for parameter Psi, case 3");
	  
	  PsiParams.push_back(new iwish("Psi", 1, m, (m*m-m)/2+m, nSamples));
	  PsiParams[0]->setHyper1(REAL(getListElement(PsiPriorHyper, "df")));
	  PsiParams[0]->setHyper2(REAL(getListElement(PsiPriorHyper, "S")));
	  PsiParams[0]->setTuning(REAL(getListElement(PsiList, "tuning")));
	  PsiParams[0]->setStarting(REAL(getListElement(PsiList, "starting")));
	  PsiParams[0]->setSampleOrder(INTEGER(getListElement(PsiList, "sample.order"))[0]);
	  UNPROTECT(2);
	}
      }else{
	error("c++ error: PsiCase is misspecified");
      }
      
    }
  
    /*********************
     Parse the Phi list
    **********************/
    vector<vprior*> PhiParams(0);
    
    SEXP PhiList, PhiPriorList, PhiPrior, PhiPriorHyper;
    int PhiCase = 0;
    
    PROTECT(PhiList = getListElement(VarList, "phi")); nProtect++;
    if(!isNewList(PhiList))
      error("c++ error: var must be a list");
    
    PhiCase = INTEGER(getListElement(PhiList, "case"))[0];
    
    if(PhiCase == 1){//all diag elements have share the same parameter (i.e., separable)
      
      PROTECT(PhiPrior = getListElement(PhiList, "prior")); 
      PriorDist = CHAR(STRING_ELT(getListElement(PhiPrior, "dist"), 0));
      PROTECT(PhiPriorHyper = getListElement(PhiPrior, "params"));
      
      if(PriorDist == "IG"){
	PhiParams.push_back(new ig("Phi", 1, 1, 1, nSamples));
	PhiParams[0]->setHyper1(&REAL(getListElement(PhiPriorHyper, "shape"))[0]);
	PhiParams[0]->setHyper2(&REAL(getListElement(PhiPriorHyper, "scale"))[0]);
	PhiParams[0]->setSampleOrder(INTEGER(getListElement(PhiList, "sample.order"))[0]);
      }else if(PriorDist == "HC"){
	PhiParams.push_back(new hc("Phi", 1, 0, 1, nSamples));
	PhiParams[0]->setHyper1(&REAL(getListElement(PhiPriorHyper, "a"))[0]);
	PhiParams[0]->setSampleOrder(INTEGER(getListElement(PhiList, "sample.order"))[0]);
      }else if(PriorDist == "LOGUNIF"){
	PhiParams.push_back(new logunif("Phi", 1, 1, 1, nSamples));
	PhiParams[0]->setHyper1(&REAL(getListElement(PhiPriorHyper, "a"))[0]);
	PhiParams[0]->setHyper2(&REAL(getListElement(PhiPriorHyper, "b"))[0]);
	PhiParams[0]->setSampleOrder(INTEGER(getListElement(PhiList, "sample.order"))[0]);
      }else if(PriorDist == "UNIF"){
	PhiParams.push_back(new logunif("Phi", 1, 1, 1, nSamples));
	PhiParams[0]->setHyper1(&REAL(getListElement(PhiPriorHyper, "a"))[0]);
	PhiParams[0]->setHyper2(&REAL(getListElement(PhiPriorHyper, "b"))[0]);
	PhiParams[0]->setSampleOrder(INTEGER(getListElement(PhiList, "sample.order"))[0]);
      }else if(PriorDist == "FIXED"){
	PhiParams.push_back(new fixedpar("Phi", 0, 0, 1, nSamples));
	PhiParams[0]->setSampleOrder(-1);
      }else{
	error("c++ error: prior misspecified for parameter Phi, case 2");
      }
      //attributes common to all priors (ignored for fixed)
      PhiParams[0]->setTuning(&REAL(getListElement(PhiList, "tuning"))[0]);
      PhiParams[0]->setStarting(&REAL(getListElement(PhiList, "starting"))[0]);	
      UNPROTECT(2);
      
    }else if(PhiCase == 2){//each diag element has its own parameter (i.e., non-separable)
      
      //PhiPriorList is a list of length m in this case
      PROTECT(PhiPriorList = getListElement(PhiList, "prior")); nProtect++; 
      
      for(i = 0; i < length(PhiPriorList); i++){
	
	//get the next prior objects
	PROTECT(PhiPrior = VECTOR_ELT(PhiPriorList, i));
	PriorDist = CHAR(STRING_ELT(getListElement(PhiPrior, "dist"), 0));
	PROTECT(PhiPriorHyper = getListElement(PhiPrior, "params"));
	
	if(PriorDist == "IG"){
	  PhiParams.push_back(new ig("Phi", 1, 1, 1, nSamples));
	  PhiParams[i]->setHyper1(&REAL(getListElement(PhiPriorHyper, "shape"))[0]);
	  PhiParams[i]->setHyper2(&REAL(getListElement(PhiPriorHyper, "scale"))[0]);
	  PhiParams[i]->setSampleOrder(INTEGER(getListElement(PhiList, "sample.order"))[i]);
	}else if(PriorDist == "HC"){
	  PhiParams.push_back(new hc("Phi", 1, 0, 1, nSamples));
	  PhiParams[i]->setHyper1(&REAL(getListElement(PhiPriorHyper, "a"))[0]);
	  PhiParams[i]->setSampleOrder(INTEGER(getListElement(PhiList, "sample.order"))[i]);
	}else if(PriorDist == "LOGUNIF"){
	  PhiParams.push_back(new logunif("Phi", 1, 1, 1, nSamples));
	  PhiParams[i]->setHyper1(&REAL(getListElement(PhiPriorHyper, "a"))[0]);
	  PhiParams[i]->setHyper2(&REAL(getListElement(PhiPriorHyper, "b"))[0]);
	  PhiParams[i]->setSampleOrder(INTEGER(getListElement(PhiList, "sample.order"))[i]);
	}else if(PriorDist == "UNIF"){
	  PhiParams.push_back(new logunif("Phi", 1, 1, 1, nSamples));
	  PhiParams[i]->setHyper1(&REAL(getListElement(PhiPriorHyper, "a"))[0]);
	  PhiParams[i]->setHyper2(&REAL(getListElement(PhiPriorHyper, "b"))[0]);
	  PhiParams[i]->setSampleOrder(INTEGER(getListElement(PhiList, "sample.order"))[i]);
	}else if(PriorDist == "FIXED"){
	  PhiParams.push_back(new fixedpar("Phi", 0, 0, 1, nSamples));
	  PhiParams[i]->setSampleOrder(-1);
	}else{
	  error("c++ error: prior misspecified for parameter Phi, case 2");
	}
	//attributes common to all priors (ignored for fixed)
	PhiParams[i]->setTuning(&REAL(getListElement(PhiList, "tuning"))[i]);
	PhiParams[i]->setStarting(&REAL(getListElement(PhiList, "starting"))[i]);
	PhiParams[i]->setSubParIndx(i);
	UNPROTECT(2);
      }
      
    }else{
      error("c++ error: PhiCase is misspecified");
    }
    
  
    /*********************
     Parse the Nu list
    **********************/
    vector<vprior*> NuParams(0);
    
    if(covModel == "matern"){
      SEXP NuList, NuPriorList, NuPrior, NuPriorHyper;
      int NuCase = 0;
      
      PROTECT(NuList = getListElement(VarList, "nu")); nProtect++;
      if(!isNewList(NuList))
	error("c++ error: var must be a list");
      
      NuCase = INTEGER(getListElement(NuList, "case"))[0];
      
      if(NuCase == 1){//all diag elements have share the same parameter (i.e., separable)
	
	PROTECT(NuPrior = getListElement(NuList, "prior")); 
	PriorDist = CHAR(STRING_ELT(getListElement(NuPrior, "dist"), 0));
	PROTECT(NuPriorHyper = getListElement(NuPrior, "params"));
	
	if(PriorDist == "IG"){
	  NuParams.push_back(new ig("Nu", 1, 1, 1, nSamples));
	  NuParams[0]->setHyper1(&REAL(getListElement(NuPriorHyper, "shape"))[0]);
	  NuParams[0]->setHyper2(&REAL(getListElement(NuPriorHyper, "scale"))[0]);
	  NuParams[0]->setSampleOrder(INTEGER(getListElement(NuList, "sample.order"))[0]);
	}else if(PriorDist == "HC"){
	  NuParams.push_back(new hc("Nu", 1, 0, 1, nSamples));
	  NuParams[0]->setHyper1(&REAL(getListElement(NuPriorHyper, "a"))[0]);
	  NuParams[0]->setSampleOrder(INTEGER(getListElement(NuList, "sample.order"))[0]);
	}else if(PriorDist == "LOGUNIF"){
	  NuParams.push_back(new logunif("Nu", 1, 1, 1, nSamples));
	  NuParams[0]->setHyper1(&REAL(getListElement(NuPriorHyper, "a"))[0]);
	  NuParams[0]->setHyper2(&REAL(getListElement(NuPriorHyper, "b"))[0]);
	  NuParams[0]->setSampleOrder(INTEGER(getListElement(NuList, "sample.order"))[0]);
	}else if(PriorDist == "UNIF"){
	  NuParams.push_back(new logunif("Nu", 1, 1, 1, nSamples));
	  NuParams[0]->setHyper1(&REAL(getListElement(NuPriorHyper, "a"))[0]);
	  NuParams[0]->setHyper2(&REAL(getListElement(NuPriorHyper, "b"))[0]);
	  NuParams[0]->setSampleOrder(INTEGER(getListElement(NuList, "sample.order"))[0]);
	}else if(PriorDist == "FIXED"){
	  NuParams.push_back(new fixedpar("Nu", 0, 0, 1, nSamples));
	  NuParams[0]->setSampleOrder(-1);
	}else{
	  error("c++ error: prior misspecified for parameter Nu, case 2");
	}
	//attributes common to all priors (ignored for fixed)
	NuParams[0]->setTuning(&REAL(getListElement(NuList, "tuning"))[0]);
	NuParams[0]->setStarting(&REAL(getListElement(NuList, "starting"))[0]);	
	UNPROTECT(2);
	
      }else if(NuCase == 2){//each diag element has its own parameter (i.e., non-separable)
	
	//NuPriorList is a list of length m in this case
	PROTECT(NuPriorList = getListElement(NuList, "prior")); nProtect++; 
	
	for(i = 0; i < length(NuPriorList); i++){
	  
	  //get the next prior objects
	  PROTECT(NuPrior = VECTOR_ELT(NuPriorList, i));
	  PriorDist = CHAR(STRING_ELT(getListElement(NuPrior, "dist"), 0));
	  PROTECT(NuPriorHyper = getListElement(NuPrior, "params"));
	  
	  if(PriorDist == "IG"){
	    NuParams.push_back(new ig("Nu", 1, 1, 1, nSamples));
	    NuParams[i]->setHyper1(&REAL(getListElement(NuPriorHyper, "shape"))[0]);
	    NuParams[i]->setHyper2(&REAL(getListElement(NuPriorHyper, "scale"))[0]);
	    NuParams[i]->setSampleOrder(INTEGER(getListElement(NuList, "sample.order"))[i]);
	  }else if(PriorDist == "HC"){
	    NuParams.push_back(new hc("Nu", 1, 0, 1, nSamples));
	    NuParams[i]->setHyper1(&REAL(getListElement(NuPriorHyper, "a"))[0]);
	    NuParams[i]->setSampleOrder(INTEGER(getListElement(NuList, "sample.order"))[i]);
	  }else if(PriorDist == "LOGUNIF"){
	    NuParams.push_back(new logunif("Nu", 1, 1, 1, nSamples));
	    NuParams[i]->setHyper1(&REAL(getListElement(NuPriorHyper, "a"))[0]);
	    NuParams[i]->setHyper2(&REAL(getListElement(NuPriorHyper, "b"))[0]);
	    NuParams[i]->setSampleOrder(INTEGER(getListElement(NuList, "sample.order"))[i]);
	  }else if(PriorDist == "UNIF"){
	    NuParams.push_back(new logunif("Nu", 1, 1, 1, nSamples));
	    NuParams[i]->setHyper1(&REAL(getListElement(NuPriorHyper, "a"))[0]);
	    NuParams[i]->setHyper2(&REAL(getListElement(NuPriorHyper, "b"))[0]);
	    NuParams[i]->setSampleOrder(INTEGER(getListElement(NuList, "sample.order"))[i]);
	  }else if(PriorDist == "FIXED"){
	    NuParams.push_back(new fixedpar("Nu", 0, 0, 1, nSamples));
	    NuParams[i]->setSampleOrder(-1);
	  }else{
	    error("c++ error: prior misspecified for parameter Nu, case 2");
	  }
	  //attributes common to all priors (ignored for fixed)
	  NuParams[i]->setTuning(&REAL(getListElement(NuList, "tuning"))[i]);
	  NuParams[i]->setStarting(&REAL(getListElement(NuList, "starting"))[i]);
	  NuParams[i]->setSubParIndx(i);
	  UNPROTECT(2);
	}
	
      }else{
	error("c++ error: NuCase is misspecified");
      }
      
    }


    /*****************************************
           Parse Beta parameters
    *****************************************/
    SEXP BetaList, BetaPrior, BetaPriorHyper;
    PROTECT(BetaList = getListElement(args, "beta.control")); nProtect++;
    PROTECT(BetaPrior = getListElement(BetaList, "prior"));  nProtect++;

    string BetaPriorDist = CHAR(STRING_ELT(getListElement(BetaPrior, "dist"), 0));
    string BetaUpdate = CHAR(STRING_ELT(getListElement(BetaList, "update"), 0));
    int nBeta = INTEGER(getListElement(BetaList, "n.beta"))[0];

    bool BetaNormal = false;
    bool BetaGibbs = true;
    
    if(BetaPriorDist == "NORMAL")
      BetaNormal = true;
    
    if(BetaUpdate == "MH")
      BetaGibbs = false;

    //allocates hypers interally based on prior and update bool
    vprior* Beta = new betapar("Beta", 0, 0, nBeta, nSamples, BetaNormal, BetaGibbs);
    if(BetaGibbs)
      Beta->setSampleOrder(-1);
    else
      Beta->setSampleOrder(INTEGER(getListElement(BetaList, "sample.order"))[0]);

    Beta->setStarting(&REAL(getListElement(BetaList, "starting"))[0]);
    
    if(BetaPriorDist == "NORMAL"){
      PROTECT(BetaPriorHyper = getListElement(BetaPrior, "params")); nProtect++;
      Beta->setHyper1(REAL(getListElement(BetaPriorHyper, "mu")));
      Beta->setHyper2(REAL(getListElement(BetaPriorHyper, "precision")));
    }


    //set tuning if needed
    if(!BetaGibbs){
      Beta->setTuning(&REAL(getListElement(BetaList, "tuning"))[0]);
    }

    /*********************
    mk common param vector
    **********************/
    vector<vprior*> params(0);
    for(i = 0; i < AParams.size(); i++)
      params.push_back(AParams[i]);

    if(!noPsi)
      for(i = 0; i < PsiParams.size(); i++)
	params.push_back(PsiParams[i]);
    
    for(i = 0; i < PhiParams.size(); i++)
      params.push_back(PhiParams[i]);

    if(covModel == "matern")
      for(i = 0; i < NuParams.size(); i++)
	params.push_back(NuParams[i]);
    
    params.push_back(Beta);
    
    
    /*********************
     describe model to fit
    **********************/
    //just for the report
    //check D matrix
    SEXP dims;
    if(!isMatrix(getListElement(args, "D")))
      error("c++ error: fixed list element 'D' is not a matrix\n");
    dims = getAttrib(getListElement(args, "D"), R_DimSymbol);

    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tGeneral model description\n");
      Rprintf("-------------------------------------------------\n");
      Rprintf("Model fit with %i observations.\n\n", INTEGER(dims)[0]);

      if(m == 1){
	Rprintf("Univariate response model.\n");
      }else{
	Rprintf("Multivariate model with %i response variables.\n\n", m);
	if(ACase == 1)
	  Rprintf("Diagonal cross-covariance matrix constructed with\na common diagonal parameter K.\n\n");	  
  	if(ACase == 2)
	  Rprintf("Diagonal cross-covariance matrix constructed with\nindependent diagonal parameters K_0:K_%i.\n\n", m-1);	  
  	if(ACase == 3)
	  Rprintf("Full cross-covariance matrix constructed with lower\ntriangle chol(K).\n\n");	  
	
	if(!noPsi){
	  if(PsiCase == 1)
	    Rprintf("Diagonal cross-covariance non-spatial error matrix\nconstructed with a common diagonal parameter Psi.\n\n");	  
	  if(PsiCase == 2)
	    Rprintf("Diagonal cross-covariance non-spatial error matrix\nconstructed with independent diagonal parameters Psi_0:Psi_%i.\n\n", m-1);	  
	  if(PsiCase == 3)
	    Rprintf("Full cross-covariance non-spatial error matrix\nconstructed with lower triangle chol(Psi).\n\n");	  
	}

	Rprintf("Using the %s spatial correlation model.\n\n", covModel.c_str());	

	if(covModel != "matern"){
	  if(PhiCase == 1)
	    Rprintf("Common phi parameter, separable model.\n");
	  if(PhiCase == 2)	 
	    Rprintf("Independent phi parameters, non-separable model.\n");	  
	}else{
	  if(PhiCase == 1)
	    Rprintf("Common phi and nu parameters, separable model.\n");
	  if(PhiCase == 2)	 
	    Rprintf("Independent phi and nu parameters, non-separable model.\n");
	}
      }
    }
    
    /*********************
       show params
    **********************/
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tModel parameters\n");
      Rprintf("-------------------------------------------------\n");
      for(j = 0; j < params.size(); j++)
	params[j]->show();
    }
    
    /*********************
      mk sample order map
    **********************/
    map<int, vector<vprior*> > sampleOrderMap;
    map<int, vector<vprior*> >::iterator sampleOrderMapIter;
    vector<vprior*> tmpVec1(0);

    //exclude non-MH params (i.e., getSampleOrder() == -1
    for(i = 0; i < params.size(); i++){
      if(params[i]->getSampleOrder() != -1){
	sampleOrderMap.insert(pair<int, vector<vprior*> >(params[i]->getSampleOrder(), tmpVec1));
      }
    } 

    //now add vprior* pointers to their respective order vectors.
    for(i = 0; i < params.size(); i++){
      for(sampleOrderMapIter = sampleOrderMap.begin(); sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++){
	if(params[i]->getSampleOrder() != -1){
	  if(params[i]->getSampleOrder() == sampleOrderMapIter->first)
	    sampleOrderMapIter->second.push_back(params[i]);
	}
      }
    }

    //show sample order
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling order\n");
      Rprintf("-------------------------------------------------\n");
      for(sampleOrderMapIter = sampleOrderMap.begin(), j=1; sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++, j++){
	Rprintf("Sample set %i: ", j);
	for(i = 0; i < sampleOrderMapIter->second.size(); i++){
	  if(sampleOrderMapIter->second[i]->getSubParIndx() != -1)
	    Rprintf("%s_%i ",sampleOrderMapIter->second[i]->getFormalName().c_str(),sampleOrderMapIter->second[i]->getSubParIndx());
	  else
	    Rprintf("%s ",sampleOrderMapIter->second[i]->getFormalName().c_str()); 
	}
	Rprintf("\n");
      }
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    /*****************************************
        set-up cov. model function pointer
    *****************************************/
    bool onePramPtr = true;
    
    void (covmodel::*cov1ParamPtr)(double, double &, double &) = NULL; 
    void (covmodel::*cov2ParamPtr)(double, double, double &, double&) = NULL;
    
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
      get X, Y, and dist
    **********************/
    int xnrow, xncol, xlength, dnrow, dlength;
    int incOne = 1;
    double *X, *Y;

    //check X matrix
    if(!isMatrix(getListElement(args, "X")))
      error("c++ error: fixed list element 'X' is not a matrix\n");
        
    PROTECT(dims = getAttrib(getListElement(args, "X"), R_DimSymbol)); nProtect++;      
    xnrow = INTEGER(dims)[0];
    xncol = INTEGER(dims)[1];

    xlength = xnrow*xncol;
    X = (double *) R_alloc(xlength, sizeof(double));

    F77_NAME(dcopy)(&xlength, REAL(getListElement(args, "X")), &incOne, X, &incOne);

    //check Y matrix
    if(!isMatrix(getListElement(args, "Y")))
      error("c++ error: fixed list element 'Y' is not a matrix\n");
     
    Y = (double *) R_alloc(xnrow, sizeof(double));

    F77_NAME(dcopy)(&xnrow, REAL(getListElement(args, "Y")), &incOne, Y, &incOne);

    //check D matrix
    if(!isMatrix(getListElement(args, "D")))
      error("c++ error: fixed list element 'D' is not a matrix\n");

    //just check dims
    dims = getAttrib(getListElement(args, "D"), R_DimSymbol);
    if(xnrow/m != INTEGER(dims)[0])
      error("c++ error: dim of 'D' is wrong\n");
    
    dnrow = xnrow/m;
    dlength = dnrow*dnrow;
    double* D = (double *) R_alloc(dlength, sizeof(double));

    F77_NAME(dcopy)(&dlength, REAL(getListElement(args, "D")), &incOne, D, &incOne);

    /*********************
         make w matrix 
           if needed
    **********************/
    SEXP w;
    if(spEffects){
      PROTECT(w = allocMatrix(REALSXP, xnrow, nSamples)); nProtect++;
      zeros(REAL(w), xnrow*nSamples);
    }


    /*********************
       DIC set-up
    **********************/
    int nDIC = nSamples - DICStart;
    double *DD;
    if(DIC){
      DD = (double *) R_alloc(nDIC, sizeof(double));
    }


    /*********************
      sampling
    **********************/
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }
    
    int s=0, info=0, status=0, rtnStatus=0;
    double logPostCurrent = 0;
    double logPostCand = 0;
    bool accept = true;
    bool first = true;
    const char lower = 'L';
    const char upper = 'U';
    const char ntran = 'N';
    const char ytran = 'T';
    const char rside = 'R';
    const char lside = 'L';
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    double junk;
    int job = 01;
    double sigmasqTmp;
    double tausqTmp;
    double logDetR = 0;
    double logDetCovCurrent = 0;
    double r = 0;
    double rUnif = 0;
    int rnrow = dnrow*m;
    int rlength = rnrow*rnrow;
    double *R = (double *) R_alloc(rlength, sizeof(double));
    double *RCurrent = (double *) R_alloc(rlength, sizeof(double));
    double *IKPsi = (double *) R_alloc(rlength, sizeof(double)); zeros(IKPsi, rlength);
    //double *wCov = (double *) R_alloc(rlength, sizeof(double)); identity(wCov, rnrow);
    double *wMu = (double *) R_alloc(rnrow, sizeof(double));
    int Alength = m*m;
    double *mmblk = (double *) R_alloc(Alength, sizeof(double));
    double *tmp = (double *) R_alloc(Alength, sizeof(double));
    double *A = (double *) R_alloc(Alength, sizeof(double));
    double *Psi = (double *) R_alloc(Alength, sizeof(double));
    double *tmpXRow = (double *) R_alloc(xnrow, sizeof(double));
    double *tmpXRow1 = (double *) R_alloc(xnrow, sizeof(double));
    double *tmpXCol = (double *) R_alloc(xncol, sizeof(double));
    double *tmpXCol1 = (double *) R_alloc(xncol, sizeof(double));
    double *tmpXRowCol = (double *) R_alloc(xnrow*xncol, sizeof(double));
    double *tmpXColCol = (double *) R_alloc(xncol*xncol, sizeof(double));
    int dicIndx = 0;

   GetRNGstate();
   for(s = 0; s < nSamples; s++){
     
     //for each sample order subset
     for(sampleOrderMapIter = sampleOrderMap.begin(); sampleOrderMapIter != sampleOrderMap.end(); sampleOrderMapIter++){
       
       //propose within the subset
       for(i = 0; i < sampleOrderMapIter->second.size(); i++)
	 sampleOrderMapIter->second[i]->propose();
       
       //form covariance matrix
       if(m == 1){
	 
	 //make the correlation matrix
	 for(i = 0; i < dlength; i++){
	   if(onePramPtr)
	     (covModelObj->*cov1ParamPtr)(PhiParams[0]->getCurrentSampleTrans(0), R[i], D[i]);
	   else //i.e., 2 parameter matern
	     (covModelObj->*cov2ParamPtr)(PhiParams[0]->getCurrentSampleTrans(0), NuParams[0]->getCurrentSampleTrans(0), R[i], D[i]);
	 }

	 //scale correlation matrix with sigmasq
	 sigmasqTmp = AParams[0]->getCurrentSampleTrans(0);
	 F77_NAME(dscal)(&dlength, &sigmasqTmp, R, &incOne);
	 
	 //use Psi
	 if(!noPsi){
	   for(i = 0; i < rnrow; i++) 
	     R[i*rnrow+i] = R[i*rnrow+i]+PsiParams[0]->getCurrentSampleTrans(0);
	 }	 

       }else{ //m > 1
	 
	 //get A sample
	 zeros(A, Alength);
	 if(ACase == 1){
	   for(i = 0; i < m; i++){A[i*m+i] = sqrt(AParams[0]->getCurrentSampleTrans(0));}
	 }else if(ACase == 2){
	   for(i = 0; i < m; i++){A[i*m+i] = sqrt(AParams[i]->getCurrentSampleTrans(0));}
	 }else{
	   F77_CALL(dcopy)(&Alength, AParams[0]->getCurrentSample(), &incOne, A, &incOne);
	 }

	 //make the correlation matrix
	 zeros(mmblk, Alength);
	 for(i = 0; i < dnrow; i++){
	   for(j = 0; j < dnrow; j++){
	     
	     for(k = 0; k < m; k++){
	       
	       if(PhiCase == 1){//separable
		 if(onePramPtr)
		   (covModelObj->*cov1ParamPtr)(PhiParams[0]->getCurrentSampleTrans(0), mmblk[k*m+k], D[j*dnrow+i]);
		 else //i.e., 2 parameter matern
		   (covModelObj->*cov2ParamPtr)(PhiParams[0]->getCurrentSampleTrans(0), NuParams[0]->getCurrentSampleTrans(0), mmblk[k*m+k], D[j*dnrow+i]);
	       }else{//require NuCase == PhiCase == 2 separable
		 if(onePramPtr)
		   (covModelObj->*cov1ParamPtr)(PhiParams[k]->getCurrentSampleTrans(0), mmblk[k*m+k], D[j*dnrow+i]);
		 else //i.e., 2 parameter matern
		   (covModelObj->*cov2ParamPtr)(PhiParams[k]->getCurrentSampleTrans(0), NuParams[k]->getCurrentSampleTrans(0), mmblk[k*m+k], D[j*dnrow+i]);
	       }
	     }
	      
	     F77_NAME(dgemm)(&ntran, &ntran, &m, &m, &m, &one, A, &m, mmblk, &m, &zero, tmp, &m);
	     F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, tmp, &m, A, &m, &zero, mmblk, &m);
	     
	     for(k = 0; k < m; k++){
	       for(l = 0; l < m; l++){
		 R[((j*m+l)*rnrow)+(i*m+k)] = mmblk[l*m+k];
		 mmblk[l*m+k] = 0.0; //zero out
	       }
	     }
	   }
	 }
	 
	 //use Psi
	 if(!noPsi){
	   
	   //get Psi sample
	   zeros(Psi, Alength);
	   if(PsiCase == 1){
	     for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiParams[0]->getCurrentSampleTrans(0));}
	   }else if(PsiCase == 2){
	     for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiParams[i]->getCurrentSampleTrans(0));}
	   }else{
	     F77_CALL(dcopy)(&Alength, PsiParams[0]->getCurrentSample(), &incOne, Psi, &incOne);
	   }
	   
	   F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, Psi, &m, Psi, &m, &zero, tmp, &m);	
	   
	   for(i = 0; i < dnrow; i++){
	     for(k = 0; k < m; k++){
	       for(l = 0; l < m; l++){
		 R[((i*m+l)*rnrow)+(i*m+k)] = R[((i*m+l)*rnrow)+(i*m+k)]+tmp[l*m+k];
	       }
	     }
	   }
	 }
       }

       //chol decom
       if(!linpack){
	 F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	 if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
       }else{
	 F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	 if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
       }
       
       //get logDet
       logDetR = 0;
       for(i = 0; i < rnrow; i++)
	 logDetR += log(R[i*rnrow+i]);
       logDetR = 2*logDetR;
       
       //finish the invert
       if(!linpack){
	 F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	 if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
       }else{
	 F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	 if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
       }
       
       for(i = 0; i < rlength; i++){
	 if(isnan(R[i]))
	   error("c++ error: Lapack dpotri routine produced an undiscovered NAN in the R^{-1} computation. This is a rogue bug we have yet to track down.  Rerun with linpack=TRUE in the run.control list.  If this error persists when linpack=TRUE please report to package maintainer.\n");
       }

       //evaluate likelihood
       logPostCand = 0;
       for(i = 0; i < params.size(); i++){ 
	 logPostCand += params[i]->logLikelihood();
       }

       F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, Beta->getCurrentSample(), 
		       &incOne, &zero, tmpXRow, &incOne);
       F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
       
       //(-1/2) * tmp` * R^{-1} * tmp
       F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
       logPostCand += -0.5*logDetR-0.5*F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne);
      
       if(first){
	 F77_CALL(dcopy)(&rlength, R, &incOne, RCurrent, &incOne);
	 logPostCurrent = logPostCand;
	 logDetCovCurrent = logDetR;
	 first = false;
       }
       
       if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){
	 F77_CALL(dcopy)(&rlength, R, &incOne, RCurrent, &incOne);
	 logPostCurrent = logPostCand;
	 logDetCovCurrent = logDetR;
       }else{
	 //reject within the subset
	 for(i = 0; i < sampleOrderMapIter->second.size(); i++)
	   sampleOrderMapIter->second[i]->reject();
       }
 
     }

     //update beta if Gibbs
     if(BetaGibbs){
       updateThetaGibbs(X, Y, xnrow, xncol, tmpXCol, 
			RCurrent, tmpXRowCol, tmpXColCol, tmpXRow, tmpXCol, tmpXCol1, 
			BetaPriorDist, Beta->getHyper1(), Beta->getHyper2());
       Beta->setCurrentSample(tmpXCol);
     }

     //
     //DIC
     //
     if(DIC && s >= DICStart){
       
       F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, Beta->getCurrentSample(), 
		       &incOne, &zero, tmpXRow, &incOne);
       F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
       
       //(-1/2) * tmp` * R^{-1} * tmp
       F77_NAME(dsymv)(&upper, &xnrow, &one, RCurrent, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
       DD[dicIndx] = logDetCovCurrent+F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne); //i.e., DD = -2.0*(-0.5*logDetCov - 0.5*dotResult);
       dicIndx++;
      }
     
     //recover w
     if(spEffects){
       
       if(m == 1){
	 if(!noPsi){
	   
	   //make the correlation matrix
	   for(i = 0; i < dlength; i++){
	     if(onePramPtr)
	       (covModelObj->*cov1ParamPtr)(PhiParams[0]->getCurrentSampleTrans(0), R[i], D[i]);
	     else //i.e., 2 parameter matern
	       (covModelObj->*cov2ParamPtr)(PhiParams[0]->getCurrentSampleTrans(0), NuParams[0]->getCurrentSampleTrans(0), R[i], D[i]);
	   }
	   
	   //chol decom
	   if(!linpack){
	     F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }
	   
	   //finish the invert
	   if(!linpack){
	     F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	     if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	   }
	   
	   //scale correlation matrix with sigmasq
	   sigmasqTmp = 1.0/AParams[0]->getCurrentSampleTrans(0);
	   F77_NAME(dscal)(&dlength, &sigmasqTmp, R, &incOne);
	   
	   //use Psi
	   for(i = 0; i < rnrow; i++) 
	     R[i*rnrow+i] = R[i*rnrow+i]+1.0/PsiParams[0]->getCurrentSampleTrans(0);
	   
	   //chol decom
	   if(!linpack){
	     F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }
	   
	   //finish the invert
	   if(!linpack){
	     F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	     if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	   }
	   
	   //make w mu
	   F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, Beta->getCurrentSample(), 
			   &incOne, &zero, tmpXRow, &incOne);
	   
	   F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
	   
	   tausqTmp = 1.0/PsiParams[0]->getCurrentSampleTrans(0);
	   F77_NAME(dscal)(&xnrow, &tausqTmp, tmpXRow, &incOne);
	   
	   F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, wMu, &incOne);
	   
	   //chol decom for the mvnorm
	   if(!linpack){
	     F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }
	   
	   mvrnorm(&REAL(w)[s*xnrow], wMu, R, xnrow, true);
	   
	 }else{
	   F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, Beta->getCurrentSample(), 
			   &incOne, &zero, &REAL(w)[s*xnrow], &incOne);
	   
	   F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, &REAL(w)[s*xnrow], &incOne);
	 }
       }else{ //m > 1
	 
	 if(!noPsi){
	   //get A sample
	   zeros(A, Alength);
	   if(ACase == 1){
	     for(i = 0; i < m; i++){A[i*m+i] = sqrt(AParams[0]->getCurrentSampleTrans(0));}
	   }else if(ACase == 2){
	     for(i = 0; i < m; i++){A[i*m+i] = sqrt(AParams[i]->getCurrentSampleTrans(0));}
	   }else{
	     F77_CALL(dcopy)(&Alength, AParams[0]->getCurrentSample(), &incOne, A, &incOne);
	   }
	   
	   //make the correlation matrix
	   zeros(mmblk, Alength);
	   for(i = 0; i < dnrow; i++){
	     for(j = 0; j < dnrow; j++){
	       
	       for(k = 0; k < m; k++){
		 
		 if(PhiCase == 1){//separable
		   if(onePramPtr)
		     (covModelObj->*cov1ParamPtr)(PhiParams[0]->getCurrentSampleTrans(0), mmblk[k*m+k], D[j*dnrow+i]);
		   else //i.e., 2 parameter matern
		     (covModelObj->*cov2ParamPtr)(PhiParams[0]->getCurrentSampleTrans(0), NuParams[0]->getCurrentSampleTrans(0), mmblk[k*m+k], D[j*dnrow+i]);
		 }else{//require NuCase == PhiCase == 2 separable
		   if(onePramPtr)
		     (covModelObj->*cov1ParamPtr)(PhiParams[k]->getCurrentSampleTrans(0), mmblk[k*m+k], D[j*dnrow+i]);
		   else //i.e., 2 parameter matern
		     (covModelObj->*cov2ParamPtr)(PhiParams[k]->getCurrentSampleTrans(0), NuParams[k]->getCurrentSampleTrans(0), mmblk[k*m+k], D[j*dnrow+i]);
		 }
	       }
	       
	       F77_NAME(dgemm)(&ntran, &ntran, &m, &m, &m, &one, A, &m, mmblk, &m, &zero, tmp, &m);
	       F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, tmp, &m, A, &m, &zero, mmblk, &m);
	       
	       for(k = 0; k < m; k++){
		 for(l = 0; l < m; l++){
		   R[((j*m+l)*rnrow)+(i*m+k)] = mmblk[l*m+k];
		   mmblk[l*m+k] = 0.0; //zero out
		 }
	       }
	     }
	   }
	   
	   //chol decom
	   if(!linpack){
	     F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }
	   
	   //finish the invert
	   if(!linpack){
	     F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	     if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	   }
	   
	   //get Psi sample
	   zeros(Psi, Alength);
	   if(PsiCase == 1){
	     for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiParams[0]->getCurrentSampleTrans(0));}
	   }else if(PsiCase == 2){
	     for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiParams[i]->getCurrentSampleTrans(0));}
	   }else{
	     F77_CALL(dcopy)(&Alength, PsiParams[0]->getCurrentSample(), &incOne, Psi, &incOne);
	   }
	   
	   //already chol in lower tri
	   F77_NAME(dpotri)(&lower, &m, Psi, &m, &info);
	   
	   //transpose fill upper tri
	   for(i = 0; i < m; i++)
	     for(j = i; j < m; j++)
	       Psi[j*m+i] = Psi[i*m+j];
	   
	   //block diag (aka I_n \otimes Psi^{-1}
 	   for(i = 0; i < dnrow; i++){
 	     for(k = 0; k < m; k++){
 	       for(l = 0; l < m; l++){
 		 IKPsi[((i*m+l)*rnrow)+(i*m+k)] = Psi[l*m+k];
 	       }
 	     }
 	   }
	   
	   F77_NAME(daxpy)(&rlength, &one, IKPsi, &incOne, R, &incOne);
	   
	   //chol decom
	   if(!linpack){
	     F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }
	   
	   //finish the invert
	   if(!linpack){
	     F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
	     if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
	   }
	   
	   //make w mu
	   F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, Beta->getCurrentSample(), 
			   &incOne, &zero, tmpXRow, &incOne);
	   F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
	   
	   F77_NAME(dsymv)(&upper, &xnrow, &one, IKPsi, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
	   F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow1, &incOne, &zero, wMu, &incOne);
	   
	   
	   //chol decom for the mvnorm
	   if(!linpack){
	     F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }else{
	     F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
	     if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
	   }
	   
	   mvrnorm(&REAL(w)[s*xnrow], wMu, R, xnrow, true);
	 }else{
	   F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, Beta->getCurrentSample(), 
			   &incOne, &zero, &REAL(w)[s*xnrow], &incOne);
	   F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, &REAL(w)[s*xnrow], &incOne);
	 }
	 
       }
     }
     
     //report
     if(verbose){
       if(status == 100){
	 Rprintf("Sampled: %i of %i, %3.2f\%\n", s, nSamples, 100.0*s/nSamples);
	 Rprintf("Acceptance rates\n");
	 for(i = 0; i < params.size(); i++){
	   if(params[i]->getSubParIndx() == -1)
	     Rprintf("\tparam: %s\t\tacc. rate: %3.2f\n", params[i]->formalName.c_str(), 100.0*params[i]->getAcceptance()/s);
	   else
	     Rprintf("\tparam: %s_%i\t\tacc. rate: %3.2f\n", params[i]->formalName.c_str(), params[i]->getSubParIndx(), 100.0*params[i]->getAcceptance()/s);
	 }
	 Rprintf("-------------------------------------------------\n");
         #ifdef Win32
	 R_FlushConsole();
        #endif
	 status = 0;
       }
       status++;
     }
     
     R_CheckUserInterrupt();
   }//end sample loop
   
   //final status report
   if(verbose){
     Rprintf("Sampled: %i of %i, %3.2f\%\n", s, nSamples, 100.0*s/nSamples);
     Rprintf("Acceptance rates\n");
     for(i = 0; i < params.size(); i++){
       if(params[i]->getSubParIndx() == -1)
	 Rprintf("\tparam: %s\t\tacc. rate: %3.2f\n", params[i]->formalName.c_str(), 100.0*params[i]->getAcceptance()/nSamples);
       else
	 Rprintf("\tparam: %s_%i\t\tacc. rate: %3.2f\n", params[i]->formalName.c_str(), params[i]->getSubParIndx(), 100.0*params[i]->getAcceptance()/nSamples);
     }
     Rprintf("-------------------------------------------------\n");
     #ifdef Win32
     R_FlushConsole();
     #endif
   }
   
   PutRNGstate();
   
   /*********************
         finish DIC
   **********************/
    double DDBar, DDBarOmega, pD;
    DDBar = DDBarOmega = pD = 0.0;

   if(DIC){

     //always inc DICStart
     DICStart++;
     
     for(i = 0; i < nDIC; i++)
       DDBar += DD[i];

     DDBar = DDBar/nDIC;

     //get means of transformed samples
     for(i = 0; i < params.size(); i++)
       params[i]->sampleMeansFromStart(DICStart);
     
     //form covariance matrix
     if(m == 1){
       
       //make the correlation matrix
       for(i = 0; i < dlength; i++){
	 if(onePramPtr)
	   (covModelObj->*cov1ParamPtr)(PhiParams[0]->getSampleMeans()[0], R[i], D[i]);
	 else //i.e., 2 parameter matern
	   (covModelObj->*cov2ParamPtr)(PhiParams[0]->getSampleMeans()[0], NuParams[0]->getSampleMeans()[0], R[i], D[i]);
       }
       
       //scale correlation matrix with sigmasq
       sigmasqTmp = AParams[0]->getSampleMeans()[0];
       F77_NAME(dscal)(&dlength, &sigmasqTmp, R, &incOne);
       
       //use Psi
       if(!noPsi){
	 for(i = 0; i < rnrow; i++) 
	   R[i*rnrow+i] = R[i*rnrow+i]+PsiParams[0]->getSampleMeans()[0];
       }	 
       
     }else{ //m > 1
       
       //get A sample
       zeros(A, Alength);
       if(ACase == 1){
	 for(i = 0; i < m; i++){A[i*m+i] = sqrt(AParams[0]->getSampleMeans()[0]);}
       }else if(ACase == 2){
	 for(i = 0; i < m; i++){A[i*m+i] = sqrt(AParams[i]->getSampleMeans()[0]);}
       }else{
	 F77_CALL(dcopy)(&Alength, AParams[0]->getSampleMeans(), &incOne, A, &incOne);
       }
       
       //make the correlation matrix
       zeros(mmblk, Alength);
       for(i = 0; i < dnrow; i++){
	 for(j = 0; j < dnrow; j++){
	   
	   for(k = 0; k < m; k++){
	     
	     if(PhiCase == 1){//separable
	       if(onePramPtr)
		 (covModelObj->*cov1ParamPtr)(PhiParams[0]->getSampleMeans()[0], mmblk[k*m+k], D[j*dnrow+i]);
	       else //i.e., 2 parameter matern
		 (covModelObj->*cov2ParamPtr)(PhiParams[0]->getSampleMeans()[0], NuParams[0]->getSampleMeans()[0], mmblk[k*m+k], D[j*dnrow+i]);
	     }else{//require NuCase == PhiCase == 2 separable
	       if(onePramPtr)
		 (covModelObj->*cov1ParamPtr)(PhiParams[k]->getSampleMeans()[0], mmblk[k*m+k], D[j*dnrow+i]);
	       else //i.e., 2 parameter matern
		 (covModelObj->*cov2ParamPtr)(PhiParams[k]->getSampleMeans()[0], NuParams[k]->getSampleMeans()[0], mmblk[k*m+k], D[j*dnrow+i]);
	     }
	   }
	   
	   F77_NAME(dgemm)(&ntran, &ntran, &m, &m, &m, &one, A, &m, mmblk, &m, &zero, tmp, &m);
	   F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, tmp, &m, A, &m, &zero, mmblk, &m);
	   
	   for(k = 0; k < m; k++){
	     for(l = 0; l < m; l++){
	       R[((j*m+l)*rnrow)+(i*m+k)] = mmblk[l*m+k];
	       mmblk[l*m+k] = 0.0; //zero out
	     }
	   }
	 }
       }
       
       //use Psi
       if(!noPsi){
	 
	 //get Psi sample
	 zeros(Psi, Alength);
	 if(PsiCase == 1){
	   for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiParams[0]->getSampleMeans()[0]);}
	 }else if(PsiCase == 2){
	   for(i = 0; i < m; i++){Psi[i*m+i] = sqrt(PsiParams[i]->getSampleMeans()[0]);}
	 }else{
	   F77_CALL(dcopy)(&Alength, PsiParams[0]->getSampleMeans(), &incOne, Psi, &incOne);
	 }
	 
	 F77_NAME(dgemm)(&ntran, &ytran, &m, &m, &m, &one, Psi, &m, Psi, &m, &zero, tmp, &m);	
	 
	 for(i = 0; i < dnrow; i++){
	   for(k = 0; k < m; k++){
	     for(l = 0; l < m; l++){
	       R[((i*m+l)*rnrow)+(i*m+k)] = R[((i*m+l)*rnrow)+(i*m+k)]+tmp[l*m+k];
	     }
	   }
	 }
       }
     }
     
     //chol decom
     if(!linpack){
       F77_NAME(dpotrf)(&upper, &rnrow, R, &rnrow, &info);
       if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
     }else{
       F77_NAME(dpofa)(R, &rnrow, &rnrow, &info);
       if(info != 0){error("c++ error: Cholesky failed (1), see ggt.sp documentation\n");}
     }
     
     //get logDet
     logDetR = 0;
     for(i = 0; i < rnrow; i++)
       logDetR += log(R[i*rnrow+i]);
     logDetR = 2*logDetR;
     
     //finish the invert
     if(!linpack){
       F77_NAME(dpotri)(&upper, &rnrow, R, &rnrow, &info);
       if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
     }else{
       F77_NAME(dpodi)(R, &rnrow, &rnrow, &junk, &job);
       if(info != 0){error("c++ error: Cholesky inverse failed (1), see ggt.sp documentation\n");}
     }
     
     for(i = 0; i < rlength; i++){
       if(isnan(R[i]))
	 error("c++ error: Lapack dpotri routine produced an undiscovered NAN in the R^{-1} computation. This is a rogue bug we have yet to track down.  Rerun with linpack=TRUE in the run.control list.  If this error persists when linpack=TRUE please report to package maintainer.\n");
     }
     
     F77_NAME(dgemv)(&ntran, &xnrow, &xncol, &negOne, X, &xnrow, Beta->getSampleMeans(), 
		     &incOne, &zero, tmpXRow, &incOne);
     F77_NAME(daxpy)(&xnrow, &one, Y, &incOne, tmpXRow, &incOne);
     
     //(-1/2) * tmp` * R^{-1} * tmp
     F77_NAME(dsymv)(&upper, &xnrow, &one, R, &xnrow, tmpXRow, &incOne, &zero, tmpXRow1, &incOne);
     DDBarOmega = logDetR+F77_NAME(ddot)(&xnrow, tmpXRow, &incOne, tmpXRow1, &incOne); //i.e., DD = -2.0*(-0.5*logDetCov - 0.5*dotResult);
   }

   /*********************
         return
   **********************/

   //make full sample matrix
   int nSampleCols = 0;
   SEXP samples;

   //transpose those parameters that need it
   for(i = 0; i < params.size(); i++) params[i]->transSamples();

   for(i = 0; i < params.size(); i++) nSampleCols += params[i]->getNPar();
   
   PROTECT(samples = allocMatrix(REALSXP, nSampleCols, nSamples+1)); nProtect++;//nSamples ob1 for starting
   zeros(REAL(samples), (nSamples+1)*nSampleCols);

   int offset = 0;
   for(j = 0; j < params.size(); j++){
     if(j > 0) offset += params[j-1]->getNPar();
     for(i = 0; i < params[j]->getNPar(); i++){
       for(s = 0; s < nSamples+1; s++){
	 REAL(samples)[s*nSampleCols+i+offset] = params[j]->getSamples()[s*params[j]->getNPar()+i];
       }
     }
   }

   //add parameter names
   SEXP samplesRowNames, samplesDimNames;
   PROTECT(samplesRowNames = allocVector(STRSXP, nSampleCols)); nProtect++;

   for(i = 0, k = 0; i < params.size(); i++)
     for(j = 0; j < params[i]->getNPar(); j++, k++)
       SET_VECTOR_ELT(samplesRowNames, k, mkChar(params[i]->getParName().c_str()));
   
   PROTECT(samplesDimNames = allocVector(VECSXP, 2)); nProtect++;
   SET_VECTOR_ELT(samplesDimNames, 0, samplesRowNames);
   setAttrib(samples, R_DimNamesSymbol, samplesDimNames);   


   //make the result list object
   int nResultListObjs = 2;
   
   if(spEffects){
     nResultListObjs++;
   }
   
   if(DIC){
     nResultListObjs++;
   }

   
   SEXP result, resultNames;
   
   PROTECT(result = allocVector(VECSXP, nResultListObjs)); nProtect++;
   PROTECT(resultNames = allocVector(STRSXP, nResultListObjs)); nProtect++;

   //set result list elements   

   //samples
   SET_VECTOR_ELT(result, 0, samples);
   SET_STRING_ELT(resultNames, 0, mkChar("p.samples"));    
   
   //acceptance
   SEXP acceptance, acceptColNames, acceptRowNames, acceptDimNames;
   PROTECT(acceptColNames = allocVector(STRSXP, 1)); nProtect++;
   SET_VECTOR_ELT(acceptColNames, 0, mkChar("acceptance.rate"));

   PROTECT(acceptance = allocMatrix(REALSXP, params.size(), 1)); nProtect++; 
   for(i = 0; i < params.size(); i++)
     REAL(acceptance)[i] = 100.0*params[i]->getAcceptance()/nSamples;

   PROTECT(acceptRowNames = allocVector(STRSXP, params.size())); nProtect++;
   
   for(i = 0; i < params.size(); i++)
     SET_VECTOR_ELT(acceptRowNames, i, mkChar(params[i]->getParName().c_str()));
   
   PROTECT(acceptDimNames = allocVector(VECSXP, 2)); nProtect++;
   SET_VECTOR_ELT(acceptDimNames, 0, acceptRowNames);
   SET_VECTOR_ELT(acceptDimNames, 1, acceptColNames);
   setAttrib(acceptance, R_DimNamesSymbol, acceptDimNames);
   
   SET_VECTOR_ELT(result, 1, acceptance);
   SET_STRING_ELT(resultNames, 1, mkChar("acceptance"));
   
   //spatial effects
   if(spEffects){
     SET_VECTOR_ELT(result, 2, w);
     SET_STRING_ELT(resultNames, 2, mkChar("sp.effects"));
   }
   

   //DIC return
   SEXP dicResults, dicResultsRowNames, dicResultsColNames, dicResultDimName;
   if(DIC){
     PROTECT(dicResults = allocMatrix(REALSXP, 4, 1)); nProtect++; //for Dbar, DbarOmega, pD, and DIC
     REAL(dicResults)[0] = DDBar;
     REAL(dicResults)[1] = DDBarOmega;
     REAL(dicResults)[2] = DDBar - DDBarOmega;
     REAL(dicResults)[3] = DDBar + DDBar - DDBarOmega;   
     
     PROTECT(dicResultsRowNames = allocVector(STRSXP, 4)); nProtect++;
     SET_VECTOR_ELT(dicResultsRowNames, 0, mkChar("bar.D"));
     SET_VECTOR_ELT(dicResultsRowNames, 1, mkChar("D.bar.Omega"));
     SET_VECTOR_ELT(dicResultsRowNames, 2, mkChar("pD"));
     SET_VECTOR_ELT(dicResultsRowNames, 3, mkChar("DIC"));
     
     PROTECT(dicResultsColNames = allocVector(STRSXP, 1)); nProtect++;
     SET_VECTOR_ELT(dicResultsColNames, 0, mkChar("value"));
     
     PROTECT(dicResultDimName = allocVector(VECSXP, 2)); nProtect++;
     SET_VECTOR_ELT(dicResultDimName, 0, dicResultsRowNames);
     SET_VECTOR_ELT(dicResultDimName, 1, dicResultsColNames);
     setAttrib(dicResults, R_DimNamesSymbol, dicResultDimName);
   }

    if(DIC){
      if(spEffects){
	SET_VECTOR_ELT(result, 3, dicResults);
	SET_STRING_ELT(resultNames, 3, mkChar("DIC"));
      }else{
	SET_VECTOR_ELT(result, 2, dicResults);
	SET_STRING_ELT(resultNames, 2, mkChar("DIC"));
      }
    }


   namesgets(result, resultNames);
   
   
   //clean-up (note if the user stops this process then 
   //the prior objects will not be destroyed, but what to do?)
   for(j = 0; j < params.size(); j++)
     delete params[j];
   
   AParams.clear();
   if(!noPsi)
     PsiParams.clear();
   PhiParams.clear();
   NuParams.clear();
   
   //unprotect
   UNPROTECT(nProtect);
   
   return(result);
  }
}














