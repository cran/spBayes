sp.predict <- function(ggt.sp.obj, pred.coords, pred.covars, start=1, end, thin=1, verbose=TRUE, ...){

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  if(missing(ggt.sp.obj)){stop("error: sp.predict expects an output object of class ggt.sp\n")}  
  if(class(ggt.sp.obj) != "ggt.sp"){stop("error: sp.predict requires an output object of class ggt.sp\n")}
  if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
  if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
  if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}

  m <- ggt.sp.obj$m
  
  if(!missing(pred.covars)){
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(nrow(pred.coords)*m != nrow(pred.covars)){stop("error: nrow(pred.coords) must be the same number of rows as pred.covars/m\n")}
  }else{
    stop("error: pred.covars must be specified\n")
  }
  
  if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}
  
  ##subsample if specified
  samples <- as.matrix(ggt.sp.obj$p.samples)

  ##make sure there is not just 2 samples
  if(nrow(samples) == 2){
    samples <- t(samples[2:nrow(samples),])##get rid of the starting value row
  }else{
    samples <- samples[2:nrow(samples),]##get rid of the starting value row
  }
  n.samples <- nrow(samples)
  
  if(missing(end))
    end <- n.samples
  
  if(!is.numeric(start) || start >= n.samples)
    stop("error: invalid start")
  if(!is.numeric(end) || end > n.samples) 
    stop("error: invalid end")
  if(!is.numeric(thin) || thin >= n.samples) 
    stop("error: invalid thin")

  
  ##if spatial effect previously calculated
  if("sp.effects" %in% names(ggt.sp.obj))
    sp.effect <- TRUE
  else
    sp.effect <- FALSE

  sp.effect <- as.integer(sp.effect)

  w <- NULL
  if(sp.effect){ ##precalculated
    w <- ggt.sp.obj$sp.effects
    w <- w[,seq(start, end, by=as.integer(thin))]
    storage.mode(w) <- "double"
  }

  samples <- samples[seq(start, end, by=as.integer(thin)),]
  n.samples <- nrow(samples)

  ##get other stuff out of ggt.sp.obj
  X <- as.matrix(ggt.sp.obj$X)
  Y <- as.matrix(ggt.sp.obj$Y)
  obs.coords <- as.matrix(ggt.sp.obj$coords)
  
  ##get parameter samples
  ##K
  K.case <- ggt.sp.obj$K.case
  if(K.case == 1){
    K <- as.matrix(samples[,1])
  }else if(K.case == 2){
    K <- samples[,1:m]
  }else{
    K <- samples[,1:((m^2-m)/2+m)]
  }
  
  samples <- samples[,(ncol(K)+1):ncol(samples)]

  ##if K.case is 3 then take the chol of each sample matrix.
  A.chol <- function(x, m){
    A <- matrix(0, m, m)
    A[lower.tri(A, diag=TRUE)] <- x
    A[upper.tri(A, diag=FALSE)] <- t(A)[upper.tri(A, diag=FALSE)]
    t(chol(A))[lower.tri(A, diag=TRUE)]
  }

  if(K.case == 3){
    K <- t(apply(K, 1, A.chol, m))
  }
  
  K <- t(K) ##trans for easy BLAS

  
  ##Psi
  no.Psi <- ggt.sp.obj$no.Psi
  Psi <- NULL
  Psi.case <- NULL
  if(!no.Psi){
    Psi.case <- ggt.sp.obj$Psi.case
    if(Psi.case == 1){
      Psi <- as.matrix(samples[,1])
    }else if(Psi.case == 2){
      Psi <- samples[,1:m]
    }else{
      Psi <- samples[,1:((m^2-m)/2+m)]
    }
    
    samples <- samples[,(ncol(Psi)+1):ncol(samples)]

    ##if Psi.case is 3 then take the chol of each sample matrix.
    if(Psi.case == 3){
      Psi <- t(apply(Psi, 1, A.chol, m))
    }
    
    Psi <- t(Psi) ##trans for easy BLAS
  }


  ##phi
  phi.case <- ggt.sp.obj$phi.case
  if(phi.case == 1){
    phi <- as.matrix(samples[,1])
  }else{
    phi <- samples[,1:m]
  }

  samples <- samples[,(ncol(phi)+1):ncol(samples)] 
  phi <- t(phi) ##trans for easy BLAS

  
  cov.model <- ggt.sp.obj$cov.model
  nu <- NULL
  if(cov.model == "matern"){ ##recall phi case == nu case
    if(phi.case == 1){
      nu <- as.matrix(samples[,1])
    }else{
      nu <- samples[,1:m]
    }
    samples <- samples[,(ncol(nu)+1):ncol(samples)] 
    nu <- t(nu) ##trans for easy BLAS

  }

  ##beta
  beta <- t(samples)

  ##just double check the dim of beta against the X and pred.X
  if(nrow(beta) != ncol(X)){
    stop("error: the number of X columns does not equal the number of sampled beta parameters\n")
  }

  ##get the prediction covars
  if(ncol(pred.covars) != ncol(X))
    {stop("error: the number of columns in the matrix or data frame specified by pred.covars must be equal to the number of covariates specified in ggt.sp.obj, which is ",ncol(X),"\n")}
  
  pred.X <- as.matrix(pred.covars)
  
  ##just double check dims
  if(ncol(pred.X) != nrow(beta))
    stop("error: this should never happen, the number of prediction covariates != number of sampled beta\n")

  ##make distance matrices
  ##observed
  obs.coords <- ggt.sp.obj$coords
  obs.D <- as.matrix(dist(cbind(obs.coords[,1], obs.coords[,2])))
  
  ##predicted
  pred.D <- as.matrix(dist(cbind(pred.coords[,1], pred.coords[,2])))

  ##predicted and observed
  predObs.D <- matrix(0, nrow(pred.coords), nrow(obs.coords))
  for(i in 1:nrow(obs.coords)){
    predObs.D[,i] <- sqrt((pred.coords[,1]-obs.coords[i,1])^2 + (pred.coords[,2]-obs.coords[i,2])^2)
  }

  ##check for zeros
  if(any(round(predObs.D, 10) == 0.0))
    stop("error: one or more predicted points coordinates is in the observed set of coordinates\n")
  
  n.pred <- as.integer(nrow(pred.coords))
  n.obs <- as.integer(nrow(obs.coords))
  
  ##get the right storage mode
  storage.mode(predObs.D) <- "double"
  storage.mode(pred.D) <- "double"
  storage.mode(obs.D) <- "double"
  storage.mode(pred.X) <- "double"
  storage.mode(beta) <- "double"
  storage.mode(X) <- "double"
  storage.mode(Y) <- "double"
  
  storage.mode(K) <- "double"
  if(!no.Psi)
    storage.mode(Psi) <- "double"
  storage.mode(phi) <- "double"
  if(cov.model == "matern")
    storage.mode(nu) <- "double"

  if(sp.effect)
    storage.mode(w) <- "double"

  args <- list("X"=X, "xrows"=as.integer(nrow(X)), "xcols"=as.integer(ncol(X)), "Y"=Y, "obs.D"=obs.D, "n.obs"=n.obs,
               "pred.X"=pred.X, "pred.D"=pred.D, "n.pred"=n.pred,
               "predObs.D"=predObs.D,
               "K"=K, "K.case"=K.case,
               "Psi"=Psi, "Psi.case"=Psi.case, "no.Psi"=as.integer(no.Psi),
               "phi"=phi, "phi.case"=phi.case,
               "nu"=nu,
               "beta"=beta,
               "cov.model"=cov.model,
               "n.samples"=n.samples,
               "m"=m,
               "sp.effect"=sp.effect, "w"=w,
               "verbose"=verbose)
  
  out <- .Call("spPredict",args)

  if(sp.effect)
    out$sp.effect <- w
  
  out
}
