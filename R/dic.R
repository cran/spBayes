sp.DIC <- function(ggt.sp.obj, DIC.marg=TRUE, DIC.unmarg=TRUE, start=1, end, thin=1, verbose=TRUE, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  if(missing(ggt.sp.obj)){stop("error: DIC expects an output object of class ggt.sp\n")}  
  if(class(ggt.sp.obj) != "ggt.sp"){stop("error: DIC requires an output object of class ggt.sp\n")}
  if(!is.logical(DIC.marg)){stop("error: DIC.marg must be of type logical\n")}
  if(!is.logical(DIC.unmarg)){stop("error: DIC.unmarg must be of type logical\n")}
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
  m <- ggt.sp.obj$m
  X <- as.matrix(ggt.sp.obj$X)
  Y <- as.matrix(ggt.sp.obj$Y)
  
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
  
  K <- t(K) ##trans for easy BLAS
  K.means <- rowMeans(K)

  
  ##Psi
  no.Psi <- ggt.sp.obj$no.Psi
  Psi <- NULL
  Psi.case <- NULL
  Psi.means <- NULL
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
    Psi <- t(Psi) ##trans for easy BLAS
    Psi.means <- rowMeans(Psi)
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
  phi.means <- rowMeans(phi)
  
  
  cov.model <- ggt.sp.obj$cov.model
  nu <- NULL
  nu.means <- NULL
  if(cov.model == "matern"){ ##recall phi case == nu case
    if(phi.case == 1){
      nu <- as.matrix(samples[,1])
    }else{
      nu <- samples[,1:m]
    }
    samples <- samples[,(ncol(nu)+1):ncol(samples)] 
    nu <- t(nu) ##trans for easy BLAS
    nu.means <- rowMeans(nu)
  }

  ##beta
  beta <- t(samples)
  beta.means <- rowMeans(beta)
  
  ##just double check the dim of beta against the X and pred.X
  if(nrow(beta) != ncol(X)){
    stop("error: the number of X columns does not equal the number of sampled beta parameters\n")
  }
  
  ##make distance matrices
  ##observed
  obs.coords <- as.matrix(ggt.sp.obj$coords)
  D <- as.matrix(dist(cbind(obs.coords[,1], obs.coords[,2])))
  n <- as.integer(nrow(obs.coords))
  
  ##get the right storage mode
  storage.mode(D) <- "double"
  storage.mode(beta) <- "double"
  storage.mode(X) <- "double"
  storage.mode(Y) <- "double"
  storage.mode(K) <- "double"
  storage.mode(K.means) <- "double"
  if(!no.Psi){
    storage.mode(Psi) <- "double"
    storage.mode(Psi.means) <- "double"
  }
  storage.mode(phi) <- "double"
  storage.mode(phi.means) <- "double"
  if(cov.model == "matern"){
    storage.mode(nu) <- "double"
    storage.mode(nu.means) <- "double"
  }
  if(sp.effect)
    storage.mode(w) <- "double"

  args <- list("DIC.marg"=as.integer(DIC.marg), "DIC.unmarg"=as.integer(DIC.unmarg),
               "X"=X, "xrows"=as.integer(nrow(X)), "xcols"=as.integer(ncol(X)), "Y"=Y, "D"=D, "n"=n,
               "K"=K, "K.case"=K.case, "K.means"=K.means,
               "Psi"=Psi, "Psi.case"=Psi.case, "Psi.means"=Psi.means, "no.Psi"=as.integer(no.Psi),
               "phi"=phi, "phi.case"=phi.case, "phi.means"=phi.means,
               "nu"=nu, "nu.means"=nu.means,
               "beta"=beta, "beta.means"=beta.means,
               "cov.model"=cov.model,
               "n.samples"=n.samples,
               "m"=m,
               "sp.effect"=sp.effect, "w"=w,
               "verbose"=verbose)
  
  out <- .Call("dic",args)
 
  out
}
