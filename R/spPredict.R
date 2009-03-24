spPredict <- function(sp.obj, pred.coords, pred.covars, start=1, end, thin=1, verbose=TRUE, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  if(missing(sp.obj)){stop("error: spPredict expects sp.obj\n")}
  if(!class(sp.obj) %in% c("spGGT","bayesGeostatExact","spLM","spMvLM", "spGLM", "spMvGLM")){
    stop("error: requires an output object of class spGGT, bayesGeostatExact, spLM, spMvLM, spGLM, or spMvGLM\n")
  }
  if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
  if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
  if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}
  
  obj.class <- class(sp.obj)

  ##
  ##prediction
  ##
  if(obj.class == "spMvGLM"){
    
    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp
    
    if(is.pp)
      knot.coords <- sp.obj$knot.coords

    family <- sp.obj$family
    Y <- sp.obj$Y
    X <- sp.obj$X
    n <- sp.obj$n
    m <- sp.obj$m
    p <- sp.obj$p
    q <- sp.obj$q
    obs.coords <- sp.obj$coords
    knots.D <- sp.obj$knots.D
    obs.D <- sp.obj$coords.D
    obs.knots.D <- sp.obj$coords.knots.D
    cov.model <- sp.obj$cov.model
    nugget <- sp.obj$nugget
    n.samples <- sp.obj$n.samples
    samples <- sp.obj$p.samples
    sp.effects <- sp.obj$recovered.effects

    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    
    if(ncol(pred.covars) != ncol(X))
      stop(paste("error: pred.covars must have ",p," columns\n"))
    
    ##thin samples and spatial effects if pre-computed
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    samples <- samples[seq(start, end, by=as.integer(thin)),]
    n.samples <- nrow(samples)

    w <- NULL
    w.str <- NULL
    
   #if(sp.effects){##Currently I'm forcing the sp.effects in spMvLM
      w <- sp.obj$sp.effects[,seq(start, end, by=as.integer(thin))]
      if(is.pp){
        w.str <- sp.obj$sp.effects.knots[,seq(start, end, by=as.integer(thin))]
      }
   # }else{##allocate to recieve the computed values
   #   w <- matrix(0, n, n.samples) #no need to compute w for the pp, but might as well on the way.
   #   if(is.pp){
   #     w.str <- matrix(0, m, n.samples)
   #   }
   # }
    
    ##get samples
    beta <- t(samples[,1:p])##transpose to simply BLAS call
    A <- NULL
    L <- NULL
    phi <- NULL
    nu <- NULL

    nltr <- m*(m-1)/2+m
    
    A.chol <- function(x, m){
      A <- matrix(0, m, m)
      A[lower.tri(A, diag=TRUE)] <- x
      A[upper.tri(A, diag=FALSE)] <- t(A)[upper.tri(A, diag=FALSE)]
      t(chol(A))[lower.tri(A, diag=TRUE)]
    }
   
    if(cov.model != "matern"){
      A <- samples[,(p+1):(p+nltr)]; A <- t(apply(A, 1, A.chol, m)); A <- t(A)
      phi <- t(samples[,(p+nltr+1):(p+nltr+m)])
      
    }else{
      A <- samples[,(p+1):(p+nltr)]; A <- t(apply(A, 1, A.chol, m)); A <- t(A)
      phi <- t(samples[,(p+nltr+1):(p+nltr+m)])
      nu <- t(samples[,(p+nltr+m+1):(p+nltr+2*m)])
    }

    pred.D <- as.matrix(dist(pred.coords))
    n.pred <- nrow(pred.D)
    pred.knots.D <- NULL
    pred.obs.D <- NULL
    
    if(is.pp){
      pred.knots.D <- matrix(0, n.pred, q)
      
      for(i in 1:n.pred){
        pred.knots.D[i,] <- sqrt((pred.coords[i,1]-knot.coords[,1])^2 + (pred.coords[i,2]-knot.coords[,2])^2)
      } 
    }else{
      pred.obs.D <- matrix(0, n.pred, n)
      
      for(i in 1:n.pred){
        pred.obs.D[i,] <- sqrt((pred.coords[i,1]-obs.coords[,1])^2 + (pred.coords[i,2]-obs.coords[,2])^2)
      } 
    }

    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(is.pp) <- "integer"
    storage.mode(modified.pp) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(q) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(A) <- "double"
    storage.mode(phi) <- "double"
    storage.mode(nu) <- "double"
    storage.mode(n.pred) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(obs.D) <- "double"
    storage.mode(pred.D) <- "double"
    storage.mode(pred.obs.D) <- "double"
    storage.mode(obs.knots.D) <- "double"
    storage.mode(knots.D) <- "double"
    storage.mode(pred.knots.D) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(w) <- "double"
    storage.mode(w.str) <- "double"
    storage.mode(sp.effects) <- "integer"
    storage.mode(verbose) <- "integer"

    out <- .Call("spMvGLMPredict", family, X, Y, is.pp, modified.pp, n, m, p, q, beta, A, phi, nu,
                 n.pred, pred.covars, obs.D, pred.D, pred.obs.D, obs.knots.D, knots.D, pred.knots.D,
                 cov.model, n.samples, w, w.str, sp.effects, verbose)
    out

     
  }else if(obj.class == "spGLM"){

    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp
    
    if(is.pp)
      knot.coords <- sp.obj$knot.coords

    family <- sp.obj$family
    Y <- sp.obj$Y
    X <- sp.obj$X
    n <- sp.obj$n
    m <- sp.obj$m
    p <- sp.obj$p
    obs.coords <- sp.obj$coords
    knots.D <- sp.obj$knots.D
    obs.D <- sp.obj$coords.D
    obs.knots.D <- sp.obj$coords.knots.D
    cov.model <- sp.obj$cov.model
    n.samples <- sp.obj$n.samples
    samples <- sp.obj$p.samples
    sp.effects <- sp.obj$recovered.effects
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    
    if(ncol(pred.covars) != ncol(X))
      stop(paste("error: pred.covars must have ",p," columns\n"))
    
    ##thin samples and spatial effects if pre-computed
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    samples <- samples[seq(start, end, by=as.integer(thin)),]
    n.samples <- nrow(samples)

    w <- NULL
    w.str <- NULL
    
    if(sp.effects){
      w <- sp.obj$sp.effects[,seq(start, end, by=as.integer(thin))]
      if(is.pp){
        w.str <- sp.obj$sp.effects.knots[,seq(start, end, by=as.integer(thin))]
      }
    }else{##allocate to recieve the computed values
      w <- matrix(0, n, n.samples) #no need to compute w for the pp, but might as well on the way.
      if(is.pp){
        w.str <- matrix(0, m, n.samples)
      }
    }
    
    ##get samples
    beta <- t(samples[,1:p])##transpose to simply BLAS call
    sigma.sq <- NULL
    phi <- NULL
    nu <- NULL

    if(cov.model != "matern"){
      sigma.sq <- samples[,"sigma.sq"]
      phi <- samples[,"phi"]
    }else{
      sigma.sq <- samples[,"sigma.sq"]
      phi <- samples[,"phi"]
      nu <- samples[,"nu"]
    }

    pred.D <- as.matrix(dist(pred.coords))
    n.pred <- nrow(pred.D)
    pred.knots.D <- NULL
    pred.obs.D <- NULL
    
    if(is.pp){
      pred.knots.D <- matrix(0, n.pred, m)
      
      for(i in 1:n.pred){
        pred.knots.D[i,] <- sqrt((pred.coords[i,1]-knot.coords[,1])^2 + (pred.coords[i,2]-knot.coords[,2])^2)
      } 
    }else{
      pred.obs.D <- matrix(0, n.pred, n)
      
      for(i in 1:n.pred){
        pred.obs.D[i,] <- sqrt((pred.coords[i,1]-obs.coords[,1])^2 + (pred.coords[i,2]-obs.coords[,2])^2)
      } 
    }
    
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(is.pp) <- "integer"
    storage.mode(modified.pp) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(sigma.sq) <- "double"
    storage.mode(phi) <- "double"
    storage.mode(nu) <- "double"
    storage.mode(n.pred) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(obs.D) <- "double"
    storage.mode(pred.D) <- "double"
    storage.mode(pred.obs.D) <- "double"
    storage.mode(obs.knots.D) <- "double"
    storage.mode(knots.D) <- "double"
    storage.mode(pred.knots.D) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(w) <- "double"
    storage.mode(w.str) <- "double"
    storage.mode(sp.effects) <- "integer"
    storage.mode(verbose) <- "integer"

    out <- .Call("spGLMPredict", family, X, Y, is.pp, modified.pp, n, m, p, beta, sigma.sq, phi, nu,
                 n.pred, pred.covars, obs.D, pred.D, pred.obs.D, obs.knots.D, knots.D, pred.knots.D,
                 cov.model, n.samples, w, w.str, sp.effects, verbose)
    out

  }else if(obj.class == "spMvLM"){

    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp
    
    if(is.pp)
      knot.coords <- sp.obj$knot.coords
    
    Y <- sp.obj$Y
    X <- sp.obj$X
    n <- sp.obj$n
    m <- sp.obj$m
    p <- sp.obj$p
    q <- sp.obj$q
    obs.coords <- sp.obj$coords
    knots.D <- sp.obj$knots.D
    obs.D <- sp.obj$coords.D
    obs.knots.D <- sp.obj$coords.knots.D
    cov.model <- sp.obj$cov.model
    nugget <- sp.obj$nugget
    n.samples <- sp.obj$n.samples
    samples <- sp.obj$p.samples
    sp.effects <- sp.obj$recovered.effects

    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    
    if(ncol(pred.covars) != ncol(X))
      stop(paste("error: pred.covars must have ",p," columns\n"))
    
    ##thin samples and spatial effects if pre-computed
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    samples <- samples[seq(start, end, by=as.integer(thin)),]
    n.samples <- nrow(samples)

    w <- NULL
    w.str <- NULL
    
   #if(sp.effects){##Currently I'm forcing the sp.effects in spMvLM
      w <- sp.obj$sp.effects[,seq(start, end, by=as.integer(thin))]
      if(is.pp){
        w.str <- sp.obj$sp.effects.knots[,seq(start, end, by=as.integer(thin))]
      }
   # }else{##allocate to recieve the computed values
   #   w <- matrix(0, n, n.samples) #no need to compute w for the pp, but might as well on the way.
   #   if(is.pp){
   #     w.str <- matrix(0, m, n.samples)
   #   }
   # }
    
    ##get samples
    beta <- t(samples[,1:p])##transpose to simply BLAS call
    A <- NULL
    L <- NULL
    phi <- NULL
    nu <- NULL

    nltr <- m*(m-1)/2+m
    
    A.chol <- function(x, m){
      A <- matrix(0, m, m)
      A[lower.tri(A, diag=TRUE)] <- x
      A[upper.tri(A, diag=FALSE)] <- t(A)[upper.tri(A, diag=FALSE)]
      t(chol(A))[lower.tri(A, diag=TRUE)]
    }
   
    if(!nugget && cov.model != "matern"){
      A <- samples[,(p+1):(p+nltr)]; A <- t(apply(A, 1, A.chol, m)); A <- t(A)
      phi <- t(samples[,(p+nltr+1):(p+nltr+m)])
      
    }else if(nugget && cov.model != "matern"){
      A <- samples[,(p+1):(p+nltr)];A <- t(apply(A, 1, A.chol, m)); A <- t(A)
      L <- samples[,(p+nltr+1):(p+2*nltr)]; L <- t(apply(L, 1, A.chol, m)); L <- t(L)
      phi <- t(samples[,(p+2*nltr+1):(p+2*nltr+m)])
      
    }else if(!nugget && cov.model == "matern"){
      A <- samples[,(p+1):(p+nltr)]; A <- t(apply(A, 1, A.chol, m)); A <- t(A)
      phi <- t(samples[,(p+nltr+1):(p+nltr+m)])
      nu <- t(samples[,(p+nltr+m+1):(p+nltr+2*m)])
      
    }else{
      A <- samples[,(p+1):(p+nltr)]; A <- t(apply(A, 1, A.chol, m)); A <- t(A)
      L <- samples[,(p+nltr+1):(p+2*nltr)]; L <- t(apply(L, 1, A.chol, m)); L <- t(L)
      phi <- t(samples[,(p+2*nltr+1):(p+2*nltr+m)])
      nu <- t(samples[,(p+2*nltr+m+1):(p+2*nltr+2*m)])
      
    }

    pred.D <- as.matrix(dist(pred.coords))
    n.pred <- nrow(pred.D)
    pred.knots.D <- NULL
    pred.obs.D <- NULL
    
    if(is.pp){
      pred.knots.D <- matrix(0, n.pred, q)
      
      for(i in 1:n.pred){
        pred.knots.D[i,] <- sqrt((pred.coords[i,1]-knot.coords[,1])^2 + (pred.coords[i,2]-knot.coords[,2])^2)
      } 
    }else{
      pred.obs.D <- matrix(0, n.pred, n)
      
      for(i in 1:n.pred){
        pred.obs.D[i,] <- sqrt((pred.coords[i,1]-obs.coords[,1])^2 + (pred.coords[i,2]-obs.coords[,2])^2)
      } 
    }

    ##fix the nugget if needed
    if(is.pp && !nugget){
      
      if(modified.pp){
        L <- rep(rep(0,m), n.samples)
      }else{
        L <- rep(diag(sqrt(0.01))[lower.tri(diag(m), diag=TRUE)], n.samples)
      }
      nugget <- TRUE
    }
   
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(is.pp) <- "integer"
    storage.mode(modified.pp) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(q) <- "integer"
    storage.mode(nugget) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(A) <- "double"
    storage.mode(L) <- "double"
    storage.mode(phi) <- "double"
    storage.mode(nu) <- "double"
    storage.mode(n.pred) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(obs.D) <- "double"
    storage.mode(pred.D) <- "double"
    storage.mode(pred.obs.D) <- "double"
    storage.mode(obs.knots.D) <- "double"
    storage.mode(knots.D) <- "double"
    storage.mode(pred.knots.D) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(w) <- "double"
    storage.mode(w.str) <- "double"
    storage.mode(sp.effects) <- "integer"
    storage.mode(verbose) <- "integer"

    out <- .Call("spMvLMPredict",X, Y, is.pp, modified.pp, n, m, p, q, nugget, beta, A, L, phi, nu,
                 n.pred, pred.covars, obs.D, pred.D, pred.obs.D, obs.knots.D, knots.D, pred.knots.D,
                 cov.model, n.samples, w, w.str, sp.effects, verbose)
    out

    
  }else if(obj.class == "spLM"){

    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp
    
    if(is.pp)
      knot.coords <- sp.obj$knot.coords
    
    Y <- sp.obj$Y
    X <- sp.obj$X
    n <- sp.obj$n
    m <- sp.obj$m
    p <- sp.obj$p
    obs.coords <- sp.obj$coords
    knots.D <- sp.obj$knots.D
    obs.D <- sp.obj$coords.D
    obs.knots.D <- sp.obj$coords.knots.D
    cov.model <- sp.obj$cov.model
    nugget <- sp.obj$nugget
    n.samples <- sp.obj$n.samples
    samples <- sp.obj$p.samples
    sp.effects <- sp.obj$recovered.effects

    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    
    if(ncol(pred.covars) != ncol(X))
      stop(paste("error: pred.covars must have ",p," columns\n"))
    
    ##thin samples and spatial effects if pre-computed
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    samples <- samples[seq(start, end, by=as.integer(thin)),]
    n.samples <- nrow(samples)

    w <- NULL
    w.str <- NULL
    
    if(sp.effects){
      w <- sp.obj$sp.effects[,seq(start, end, by=as.integer(thin))]
      if(is.pp){
        w.str <- sp.obj$sp.effects.knots[,seq(start, end, by=as.integer(thin))]
      }
    }else{##allocate to recieve the computed values
      w <- matrix(0, n, n.samples) #no need to compute w for the pp, but might as well on the way.
      if(is.pp){
        w.str <- matrix(0, m, n.samples)
      }
    }
    
    ##get samples
    beta <- t(samples[,1:p])##transpose to simply BLAS call
    sigma.sq <- NULL
    tau.sq <- NULL
    phi <- NULL
    nu <- NULL

    if(!nugget && cov.model != "matern"){
      sigma.sq <- samples[,"sigma.sq"]
      phi <- samples[,"phi"]
      
    }else if(nugget && cov.model != "matern"){
      sigma.sq <- samples[,"sigma.sq"]
      tau.sq <- samples[,"tau.sq"]
      phi <- samples[,"phi"]
      
    }else if(!nugget && cov.model == "matern"){

      sigma.sq <- samples[,"sigma.sq"]
      phi <- samples[,"phi"]
      nu <- samples[,"nu"]

    }else{
      sigma.sq <- samples[,"sigma.sq"]
      tau.sq <- samples[,"tau.sq"]
      phi <- samples[,"phi"]
      nu <- samples[,"nu"]
    }

    pred.D <- as.matrix(dist(pred.coords))
    n.pred <- nrow(pred.D)
    pred.knots.D <- NULL
    pred.obs.D <- NULL
    
    if(is.pp){
      pred.knots.D <- matrix(0, n.pred, m)
      
      for(i in 1:n.pred){
        pred.knots.D[i,] <- sqrt((pred.coords[i,1]-knot.coords[,1])^2 + (pred.coords[i,2]-knot.coords[,2])^2)
      } 
    }else{
      pred.obs.D <- matrix(0, n.pred, n)
      
      for(i in 1:n.pred){
        pred.obs.D[i,] <- sqrt((pred.coords[i,1]-obs.coords[,1])^2 + (pred.coords[i,2]-obs.coords[,2])^2)
      } 
    }

    ##fix the nugget if needed
    if(is.pp && !nugget){
      
      if(modified.pp){
        tau.sq <- rep(0.0, n.samples)
      }else{
        tau.sq <- rep(0.01, n.samples)
      }
      nugget <- TRUE
    }
    
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(is.pp) <- "integer"
    storage.mode(modified.pp) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(nugget) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(sigma.sq) <- "double"
    storage.mode(tau.sq) <- "double"
    storage.mode(phi) <- "double"
    storage.mode(nu) <- "double"
    storage.mode(n.pred) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(obs.D) <- "double"
    storage.mode(pred.D) <- "double"
    storage.mode(pred.obs.D) <- "double"
    storage.mode(obs.knots.D) <- "double"
    storage.mode(knots.D) <- "double"
    storage.mode(pred.knots.D) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(w) <- "double"
    storage.mode(w.str) <- "double"
    storage.mode(sp.effects) <- "integer"
    storage.mode(verbose) <- "integer"

    out <- .Call("spLMPredict",X, Y, is.pp, modified.pp, n, m, p, nugget, beta, sigma.sq, tau.sq, phi, nu,
                 n.pred, pred.covars, obs.D, pred.D, pred.obs.D, obs.knots.D, knots.D, pred.knots.D,
                 cov.model, n.samples, w, w.str, sp.effects, verbose)
    out
    
  }else if(obj.class == "bayesGeostatExact"){
    
    X <- sp.obj$args$X
    n <- sp.obj$args$n
    p <- sp.obj$args$p
    Y <- sp.obj$args$Y
    coords <- sp.obj$args$coords
    cov.model <- sp.obj$args$cov.model
    samples <- sp.obj$p.samples
    phi <- sp.obj$args$phi
    n.samples <- sp.obj$args$n.samples
    
    if(cov.model == "matern")
      nu <- sp.obj$args$nu
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    
    if(ncol(pred.covars) != ncol(X))
      stop(paste("error: pred.covars must have ",p," columns\n"))
    
    ##thin
    if(missing(end))
      end <- n.samples
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    samples <- samples[seq(start, end, by=as.integer(thin)),]
    n.samples <- nrow(samples)
    
    ##get samples
    beta <- as.matrix(samples[,1:p])
    tau.sq <- samples[,"tau.sq"]
    sigma.sq <- samples[,"sigma.sq"]    
    
    ##make R
    D <- as.matrix(dist(coords))
    
    if(cov.model == "exponential"){
      R <- exp(-phi*D)
    }else if(cov.model == "matern"){
      R <- (D*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D*phi, nu=nu)
      diag(R) <- 1
    }else if(cov.model == "gaussian"){
      R <- exp(-1*((phi*D)^2))
    }else if(cov.model == "spherical"){
      R <- D
      R[TRUE] <- 1
      R[D > 0 & D < 1/phi] <- 1-1.5*phi*D[D > 0 & D <= 1/phi]+0.5*((phi*D[D > 0 & D <= 1/phi])^3)
      R[D >= 1/phi] <- 0   
    }else{
      stop("error: in spPredict, specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
    }
    
    
    n.pred <- nrow(pred.coords)
    
    y.pred <- matrix(0, n.pred, n.samples)
    
    R.eigen <- eigen(R)
    R.vals <- R.eigen$values
    R.vecs <- R.eigen$vectors
    R.vects.t <- t(R.vecs)


    if(verbose)
      cat("Predicting ...\n")

    report <- 1
    
    ##for each pred point by each sample
    for(i in 1:n.pred){
      
      D.pred <- sqrt((pred.coords[i,1]-coords[,1])^2 + (pred.coords[i,2]-coords[,2])^2)
      
      if(cov.model == "exponential"){
        gamma <- exp(-phi*D.pred)
      }else if(cov.model == "matern"){
        gamma <- (D.pred*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D.pred*phi, nu=nu)
        diag(gamma) <- 1
      }else if(cov.model == "gaussian"){
        gamma <- exp(-1*((phi*D.pred)^2))
      }else if(cov.model == "spherical"){
        gamma <- D.pred
        gamma[TRUE] <- 1
        gamma[D.pred > 0 & D.pred < 1/phi] <- 1-1.5*phi*D.pred[D.pred > 0 & D.pred <= 1/phi]+0.5*((phi*D.pred[D.pred > 0 & D.pred <= 1/phi])^3)
        gamma[D.pred >= 1/phi] <- 0   
      }else{
        stop("error: in spPredict, specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
      }
      
      gamma <- as.matrix(gamma)
      
      for(s in 1:n.samples){
        
        R.inv <- R.vecs%*%diag(1/(R.vals+tau.sq[s]/sigma.sq[s]))%*%t(R.vecs)
        
        mu <- pred.covars[i,]%*%beta[s,]+t(gamma)%*%R.inv%*%(Y-X%*%beta[s,])
        S <- sigma.sq[s]*(1-t(gamma)%*%R.inv%*%gamma)+tau.sq[s]
        
        y.pred[i,s] <- rnorm(1, mu, sqrt(S))
        
      }

      
      if(verbose){
        if(report == 10){
          cat(paste("Percent complete: ",100*i/n.pred,"\n",sep=""))
          report <- 0
        }
        report <- report+1
      }
    }
    
    y.pred
    
  }else if(obj.class == "spGGT"){
    
    m <- sp.obj$m
    
    if(!missing(pred.covars)){
      if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
      if(nrow(pred.coords)*m != nrow(pred.covars)){stop("error: nrow(pred.coords) must be the same number of rows as pred.covars/m\n")}
    }else{
      stop("error: pred.covars must be specified\n")
    }
    
    if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}
    
    ##subsample if specified
    samples <- as.matrix(sp.obj$p.samples)
    
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
    if("sp.effects" %in% names(sp.obj))
      sp.effect <- TRUE
    else
      sp.effect <- FALSE
    
    sp.effect <- as.integer(sp.effect)
    
    w <- NULL
    if(sp.effect){ ##precalculated
      w <- sp.obj$sp.effects
      w <- w[,seq(start, end, by=as.integer(thin))]
      storage.mode(w) <- "double"
    }
    
    samples <- samples[seq(start, end, by=as.integer(thin)),]
    n.samples <- nrow(samples)
    
    ##get other stuff out of sp.obj
    X <- as.matrix(sp.obj$X)
    Y <- as.matrix(sp.obj$Y)
    obs.coords <- as.matrix(sp.obj$coords)
    
    ##get parameter samples
    ##K
    K.case <- sp.obj$K.case
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
    no.Psi <- sp.obj$no.Psi
    Psi <- NULL
    Psi.case <- NULL
    if(!no.Psi){
      Psi.case <- sp.obj$Psi.case
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
    phi.case <- sp.obj$phi.case
    if(phi.case == 1){
      phi <- as.matrix(samples[,1])
    }else{
      phi <- samples[,1:m]
    }
    
    samples <- samples[,(ncol(phi)+1):ncol(samples)] 
    phi <- t(phi) ##trans for easy BLAS
    
    
    cov.model <- sp.obj$cov.model
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
      {stop("error: the number of columns in the matrix or data frame specified by pred.covars must be equal to the number of covariates specified in sp.obj, which is ",ncol(X),"\n")}
    
    pred.X <- as.matrix(pred.covars)
    
    ##just double check dims
    if(ncol(pred.X) != nrow(beta))
      stop("error: this should never happen, the number of prediction covariates != number of sampled beta\n")
    
    ##make distance matrices
    ##observed
    obs.coords <- sp.obj$coords
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
    
  }else{
    stop("error: requires an output object of class spGGT, bayesGeostatExact, or spLM\n")
  }
  
}
