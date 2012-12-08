spPredict <- function(sp.obj, pred.coords, pred.covars, start=1, end, thin=1, verbose=TRUE, n.report=100, ...){
  
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
  if(!class(sp.obj) %in% c("spLM", "spMvLM", "spGLM", "spMvGLM","bayesLMRef","nonSpGLM","nonSpMvGLM")){
    stop("error: requires an output object of class spLM, spMvLM, spGLM, spMvGLM, bayesLMRef, nonSpGLM, or nonSpMvGLM\n")
  }

  obj.class <- class(sp.obj)

  ##
  ##non spatial model prediction
  ##

  if(obj.class %in% c("nonSpGLM", "nonSpMvGLM")){

    X <- sp.obj$X
    Y <- sp.obj$Y
    family <- sp.obj$family
    p.beta.samples <- sp.obj$p.beta.samples
    n.samples <- nrow(p.beta.samples)
    
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
       
    s.indx <- seq(start, end, by=as.integer(thin))
    
    p.beta.samples <- p.beta.samples[s.indx,,drop=FALSE]
    n.samples <- nrow(p.beta.samples)

    out <- list()
    
    if(family == "binomial"){
      out$p.predictive.samples <- apply(p.beta.samples, 1,  function(s){1/(1+exp(-pred.covars%*%s))})          
    }else{##poisson
       out$p.predictive.samples <- apply(p.beta.samples, 1,  function(s){exp(pred.covars%*%s)}) 
     }

    return(out)
  }

  if(obj.class == "bayesLMRef"){

    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    n <- nrow(X)
    p.beta.tauSq.samples <- sp.obj$p.beta.tauSq.samples
    n.samples <- nrow(p.beta.tauSq.samples)

    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
       
    s.indx <- seq(start, end, by=as.integer(thin))
    
    p.beta.tauSq.samples <- p.beta.tauSq.samples[s.indx,]
    n.samples <- nrow(p.beta.tauSq.samples)
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}

    out <- list()
    out$p.predictive.samples <- apply(p.beta.tauSq.samples, 1, function(s){rnorm(n, pred.covars%*%s[1:p], sqrt(s[p+1]))})
    
    return(out)
  }

  
  ##
  ##spatial model prediction
  ##
  if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
  if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
  if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}
  
  if(obj.class %in% c("spGLM", "spMvGLM")){

    family <- sp.obj$family
    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    m <- 1 ##for spGLM
    if(obj.class == "spMvGLM"){m <- sp.obj$m}
    obs.coords <- sp.obj$coords
    n <- nrow(obs.coords)
    cov.model <- sp.obj$cov.model
    p.beta.theta.samples <- sp.obj$p.beta.theta.samples
    n.samples <- nrow(p.beta.theta.samples)
    is.pp <- sp.obj$is.pp
    q <- nrow(pred.coords)
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
    
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
       
    s.indx <- seq(start, end, by=as.integer(thin))

    p.beta.theta.samples <- t(p.beta.theta.samples[s.indx,,drop=FALSE])
    n.samples <- ncol(p.beta.theta.samples)

    knots.obs.D  <- NULL
    knots.D <- NULL
    knots.pred.D <- NULL
    obs.pred.D <- NULL
    obs.D <- NULL
    pred.D <- NULL
    p.w.samples <- NULL
    
    if(is.pp){
      knot.coords <- sp.obj$knot.coords
      g <- nrow(knot.coords)
      p.w.samples <- sp.obj$p.w.knots.samples[,s.indx,drop=FALSE]
      knots.D <- iDist(knot.coords)
      pred.knots.D <- iDist(pred.coords, knot.coords)
    }else{
      p.w.samples <- sp.obj$p.w.samples[,s.indx,drop=FALSE]
      obs.pred.D <- iDist(obs.coords, pred.coords)
      obs.D <- iDist(obs.coords)
    }
    
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(q) <- "integer"
    storage.mode(p.beta.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(p.w.samples) <- "double"   
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    if(is.pp){
      
      storage.mode(g) <- "integer"
      storage.mode(knots.D) <- "double"
      storage.mode(pred.knots.D) <- "double"

      out <- .Call("spPPMvGLMPredict", family, X, Y, n, m, g, p, pred.covars, q, knots.D, pred.knots.D, 
                   p.beta.theta.samples, p.w.samples, n.samples,
                   cov.model, verbose, n.report)
    }else{
      
      storage.mode(obs.pred.D) <- "double"
      storage.mode(obs.D) <- "double"
       
      out <- .Call("spMvGLMPredict", family, X, Y, n, m, p, pred.covars, q, obs.D, obs.pred.D, 
                   p.beta.theta.samples, p.w.samples, n.samples, cov.model,
                   verbose, n.report)
    }

    out
        
  }else if(obj.class == "spLM"){

    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    n <- nrow(X)
    obs.coords <- sp.obj$coords
    cov.model <- sp.obj$cov.model
    p.theta.samples <- sp.obj$p.theta.samples
    n.samples <- nrow(p.theta.samples)
    is.pp <- sp.obj$is.pp
    nugget <- sp.obj$nugget
    beta.prior <- sp.obj$beta.prior
    beta.Norm <- sp.obj$beta.Norm
    ##r.indx <- sp.obj$r.indx
    get.beta <- sp.obj$get.beta
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
    
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    s.indx <- seq(start, end, by=as.integer(thin))
    
    beta <- NULL;
    
    if(is.pp){
      beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
      knot.coords <- sp.obj$knot.coords
      m <- nrow(knot.coords)
      modified.pp <- sp.obj$modified.pp
    }   
   
    p.theta.samples <- p.theta.samples[s.indx,,drop=FALSE]##note, for flat we could use the p.theta.recover.samples
    
    n.samples <- nrow(p.theta.samples)
    
    ##recover beta if needed
    if(!is.pp && beta.prior == "flat"){
      ## if(all(s.indx %in% r.indx)){
      ##   beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
      ## }else{
        beta <- spRecover(sp.obj, get.beta=TRUE, get.w=FALSE, start=start, end=end, thin=thin)$p.beta.recover.samples
      ## }
    }
    
    sigma.sq.indx <- 0; tau.sq.indx <- 0; phi.indx <- 0; nu.indx <- 0
    
    if(!nugget && cov.model != "matern"){
      sigma.sq.indx <- 0; phi.indx <- 1
    }else if(nugget && cov.model != "matern"){
      sigma.sq.indx <- 0; tau.sq.indx <- 1; phi.indx <- 2
    }else if(!nugget && cov.model == "matern"){
      sigma.sq.indx <- 0; phi.indx <- 1; nu.indx <- 2
    }else{
      sigma.sq.indx <- 0; tau.sq.indx <- 1; phi.indx <- 2; nu.indx <- 3
    }
    
    knots.obs.D  <- NULL
    knots.D <- NULL
    knots.pred.D <- NULL
    obs.pred.D <- NULL
    obs.D <- NULL
    
    if(is.pp){
      knots.obs.D <- iDist(knot.coords, obs.coords)
      knots.D <- iDist(knot.coords)
      knots.pred.D <- iDist(knot.coords, pred.coords)
    }else{
      obs.pred.D <- iDist(obs.coords, pred.coords)
      obs.D <- iDist(obs.coords)
    }
    
    n.pred <- nrow(pred.coords)
    
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(n.pred) <- "integer"
    storage.mode(p.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(sigma.sq.indx) <- "integer"
    storage.mode(tau.sq.indx) <- "integer"
    storage.mode(phi.indx) <- "integer"
    storage.mode(nu.indx) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    if(is.pp){
      
      storage.mode(m) <- "integer"
      storage.mode(knots.D) <- "double"
      storage.mode(knots.obs.D) <- "double"
      storage.mode(knots.pred.D) <- "double"
      storage.mode(modified.pp) <- "integer"
      
      out <- .Call("spPPLMPredict", X, Y, n, p, m, pred.covars, n.pred,
                   p.theta.samples, n.samples,
                   beta, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                   nugget, knots.D, knots.obs.D, knots.pred.D, cov.model, modified.pp, verbose, n.report)
    }else{
      
     storage.mode(obs.pred.D) <- "double"
     storage.mode(obs.D) <- "double"
     
     out <- .Call("spLMPredict", X, Y, n, p, pred.covars, n.pred,
                  p.theta.samples, n.samples,
                  beta.prior, beta.Norm, beta, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                  obs.D, obs.pred.D, cov.model, nugget, verbose, n.report)
   }
    
    out
  
  }else if(obj.class == "spMvLM"){

    X <- sp.obj$X
    Y <- sp.obj$Y
    p <- ncol(X)
    m <- sp.obj$m
    obs.coords <- sp.obj$coords
    n <- nrow(obs.coords)
    cov.model <- sp.obj$cov.model
    p.theta.samples <- sp.obj$p.theta.samples
    n.samples <- nrow(p.theta.samples)
    is.pp <- sp.obj$is.pp
    nugget <- sp.obj$nugget
    Psi.diag <- sp.obj$Psi.diag
    beta.prior <- sp.obj$beta.prior
    beta.Norm <- sp.obj$beta.Norm
    ##r.indx <- sp.obj$r.indx
    get.beta <- sp.obj$get.beta
    
    q <- nrow(pred.coords)
    
    ##check covars
    if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
    
    ##subsamples
    if(missing(end)){end <- n.samples}
    
    if(!is.numeric(start) || start >= n.samples)
      stop("error: invalid start")
    if(!is.numeric(end) || end > n.samples) 
      stop("error: invalid end")
    if(!is.numeric(thin) || thin >= n.samples) 
      stop("error: invalid thin")
    
    s.indx <- seq(start, end, by=as.integer(thin))
    
    beta <- NULL;
    
    if(is.pp){
      beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
      knot.coords <- sp.obj$knot.coords
      g <- nrow(knot.coords)
      modified.pp <- sp.obj$modified.pp
    }  
       
    p.theta.samples <- t(p.theta.samples[s.indx,,drop=FALSE])##note, for flat we could use the p.theta.recover.samples
    n.samples <- ncol(p.theta.samples)
    
    ##recover beta if needed (note, beta samples not needed for beta normal)
    if(!is.pp && beta.prior == "flat"){
      ## if(all(s.indx %in% r.indx)){
      ##   beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
      ## }else{
        beta <- spRecover(sp.obj, get.beta=TRUE, get.w=FALSE, start=start, end=end, thin=thin)$p.beta.recover.samples
      ## }
    }
        
    knots.obs.D  <- NULL
    knots.D <- NULL
    knots.pred.D <- NULL
    obs.pred.D <- NULL
    obs.D <- NULL
    
    if(is.pp){
      knots.obs.D <- iDist(knot.coords, obs.coords)
      knots.D <- iDist(knot.coords)
      knots.pred.D <- iDist(knot.coords, pred.coords)
    }else{
      obs.pred.D <- iDist(obs.coords, pred.coords)
      obs.D <- iDist(obs.coords)
    }
    
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(pred.covars) <- "double"
    storage.mode(q) <- "integer"
    storage.mode(p.theta.samples) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(nugget) <- "integer"
    storage.mode(Psi.diag) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.report) <- "integer"
    
    if(is.pp){
      
      storage.mode(g) <- "integer"
      storage.mode(knots.D) <- "double"
      storage.mode(knots.obs.D) <- "double"
      storage.mode(knots.pred.D) <- "double"
      storage.mode(modified.pp) <- "integer"

      out <- .Call("spPPMvLMPredict", X, Y, n, m, g, p, pred.covars, q,
                   knots.D, knots.obs.D, knots.pred.D, 
                   p.theta.samples, beta, n.samples,
                   nugget, Psi.diag, cov.model,
                   modified.pp, verbose, n.report)
    }else{
      
     storage.mode(obs.pred.D) <- "double"
     storage.mode(obs.D) <- "double"


     Z <- t(pred.covars)
     storage.mode(Z) <- "double"
     
     out <- .Call("spMvLMPredict", X, Y, n, m, p, Z, q, obs.D, obs.pred.D, 
                  p.theta.samples, beta, n.samples,
                  beta.prior, beta.Norm,
                  nugget, Psi.diag, cov.model,
                  verbose, n.report)
   }
    
    out

    
  }else{
    stop("error: wrong class\n")
  }
  
}
