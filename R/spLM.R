spLM <- function(formula, data = parent.frame(), coords, knots, 
                 starting, tuning, priors, cov.model,
                 modified.pp=TRUE, amcmc, n.samples,
                 verbose=TRUE, n.report=100, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  ####################################################
  ##Formula
  ####################################################
  if(missing(formula)){stop("error: formula must be specified")}
  
  ##if(class(formula) == "formula"){
  if(inherits(formula, "formula")){
      
    holder <- parseFormula(formula, data)
    Y <- holder[[1]]
    X <- as.matrix(holder[[2]])
    x.names <- holder[[3]]

  }else{
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)
  
  ##make sure storage mode is correct
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

  ####################################################
  ##sampling method
  ####################################################
  n.batch <- 0
  batch.length <- 0
  accept.rate <- 0
  is.amcmc <- TRUE
  
  if(missing(amcmc)){
   
    if(missing(n.samples)){stop("error: n.samples needs to be specified")}
    
    n.batch <- n.samples
    batch.length <- 1
    is.amcmc <- FALSE
        
  }else{
    
    names(amcmc) <- tolower(names(amcmc))
    
    if(!"n.batch" %in% names(amcmc)){stop("error: n.batch must be specified in amcmc list")}
    n.batch <- amcmc[["n.batch"]]
    
    if(!"batch.length" %in% names(amcmc)){stop("error: batch.length must be specified in amcmc list")}
    batch.length <- amcmc[["batch.length"]]

    if(!"accept.rate" %in% names(amcmc)){
      warning("accept.rate was not specified in the amcmc list and was therefore set to the default 0.43")
      accept.rate <- 0.43
    }else{
      accept.rate <- amcmc[["accept.rate"]]
    }

  }

  storage.mode(is.amcmc) <- "integer"
  storage.mode(n.batch) <- "integer"
  storage.mode(batch.length) <- "integer"
  storage.mode(accept.rate) <- "double"

  ####################################################
  ##Fit non-spatial model if specified
  ####################################################
  if(missing(coords)){
    return(bayesLMRef(lm(Y~X-1), n.batch*batch.length))
  }
  
  ####################################################
  ##Distance matrices
  ####################################################
  
  ####################
  ##Coords
  #################### 
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or then number of rows is different than
          data used in the model formula")
  }
  
  ####################
  ##Knots
  ####################
  is.pp <- FALSE
  
  if(!missing(knots)){
    
    if(is.vector(knots) && length(knots) %in% c(2,3)){
      
      ##allow single knot dim
      if(knots[1] > 1){
        x.knots <- seq(min(coords[,1]), max(coords[,1]), length.out=knots[1])
      }else{
        x.knots <- (max(coords[,1])-min(coords[,1]))/2
      }
      
      if(knots[2] > 1){
        y.knots <- seq(min(coords[,2]), max(coords[,2]), length.out=knots[2])
      }else{
        y.knots <- (max(coords[,2])-min(coords[,2]))/2
      }
      
      ##if not single knot then adjust out half distance on all sides
      if(length(knots) == 2){
        if(knots[1] > 1){
          x.int <- (x.knots[2]-x.knots[1])/2
          x.knots <- seq(min(x.knots)-x.int, max(x.knots)+x.int, length.out=knots[1])
        }
        
        if(knots[2] > 1){
          y.int <- (y.knots[2]-y.knots[1])/2
          y.knots <- seq(min(y.knots)-y.int, max(y.knots)+y.int, length.out=knots[2])
        }
        
        knot.coords <- as.matrix(expand.grid(x.knots, y.knots))
        is.pp <- TRUE
      }else{   
        if(knots[1] > 1){
          x.int <- knots[3]
          x.knots <- seq(min(x.knots)-x.int, max(x.knots)+x.int, length.out=knots[1])
        }
        
        if(knots[2] > 1){
          y.int <- knots[3]
          y.knots <- seq(min(y.knots)-y.int, max(y.knots)+y.int, length.out=knots[2])
        }
        
        knot.coords <- as.matrix(expand.grid(x.knots, y.knots))
        is.pp <- TRUE
      }
      
    }else if(is.matrix(knots) && ncol(knots) == 2){
      knot.coords <- knots
      is.pp <- TRUE
    }else{
      stop("error: knots is misspecified")
    }
  }
 
  m <- 0
  coords.D <- 0
  knots.D <- 0
  knots.coords.D <- 0
  
  if(is.pp){
    knots.D <- iDist(knot.coords)
    m <- nrow(knots.D)
    knots.coords.D <- iDist(knot.coords, coords)
  }else{
    coords.D <- iDist(coords)
  }

  storage.mode(modified.pp) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(coords.D) <- "double"
  storage.mode(knots.D) <- "double"
  storage.mode(knots.coords.D) <- "double"

  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}
  
  ####################################################
  ##Priors
  ####################################################
  beta.Norm <- 0
  beta.prior <- "flat"
  sigma.sq.IG <- 0
  tau.sq.IG <- 0
  nu.Unif <- 0
  phi.Unif <- 0
  nugget <- FALSE

  if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
  names(priors) <- tolower(names(priors))

  if("beta.norm" %in% names(priors)){
    beta.Norm <- priors[["beta.norm"]]
    if(!is.list(beta.Norm) || length(beta.Norm) != 2){stop("error: beta.Norm must be a list of length 2")}
    if(length(beta.Norm[[1]]) != p ){stop(paste("error: beta.Norm[[1]] must be a vector of length, ",p, "",sep=""))}
    if(length(beta.Norm[[2]]) != p^2 ){stop(paste("error: beta.Norm[[2]] must be a ",p,"x",p," covariance matrix",sep=""))}
    beta.prior <- "normal"
  }
  
  if(!"sigma.sq.ig" %in% names(priors)){stop("error: sigma.sq.IG must be specified")}
  sigma.sq.IG <- priors[["sigma.sq.ig"]]
  
  if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){stop("error: sigma.sq.IG must be a vector of length 2")}
  if(any(sigma.sq.IG <= 0)){stop("error: sigma.sq.IG must be a positive vector of length 2")}
  
  if("tau.sq.ig" %in% names(priors)){
    tau.sq.IG <- priors[["tau.sq.ig"]]
    
    if(!is.vector(tau.sq.IG) || length(tau.sq.IG) != 2){stop("error: tau.sq.IG must be a vector of length 2")}
    if(any(tau.sq.IG <= 0)){stop("error: tau.sq.IG must be a positive vector of length 2")}
    nugget <- TRUE
  }
  
  if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified")}
  phi.Unif <- priors[["phi.unif"]]
  
  if(!is.vector(phi.Unif) || length(phi.Unif) != 2){stop("error: phi.Unif must be a vector of length 2")}
  if(any(phi.Unif <= 0, phi.Unif[1] >= phi.Unif[2])){stop("error: phi.Unif must be a positive vector of length 2 with element 1 < element 2")}
  
  if(cov.model == "matern"){
    
    if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified")}
    nu.Unif <- priors[["nu.unif"]]
    
    if(!is.vector(nu.Unif) || length(nu.Unif) != 2){stop("error: nu.Unif must be a vector of length 2")}
    if(any(nu.Unif <= 0, nu.Unif[1] >= nu.Unif[2])){stop("error: nu.Unif must be a positive vector of length 2 with element 1 < element 2")}
  }
  
  storage.mode(sigma.sq.IG) <- "double"
  storage.mode(tau.sq.IG) <- "double"
  storage.mode(nu.Unif) <- "double"
  storage.mode(phi.Unif) <- "double"
  storage.mode(nugget) <- "integer"

  ####################################################
  ##Starting values
  ####################################################
  beta.starting <- 0
  sigma.sq.starting <- 0
  tau.sq.starting <- 0
  phi.starting <- 0
  nu.starting <- 0

  if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
  
  names(starting) <- tolower(names(starting))   
  
  if(is.pp){
    if(!"beta" %in% names(starting)){
      beta.starting <- as.vector(coefficients(lm(Y~X-1)))
    }else{
      beta.starting <- starting[["beta"]]
    }
  }
 
  if(!"sigma.sq" %in% names(starting)){stop("error: sigma.sq must be specified in starting value list")}
  sigma.sq.starting <- starting[["sigma.sq"]][1]

  if(nugget){
    if(!"tau.sq" %in% names(starting)){stop("error: a prior was spcified for tau.sq therefore tau.sq must be specified in starting value list")}
    tau.sq.starting <- starting[["tau.sq"]][1]
  }
  
  if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting value list")}
  phi.starting <- starting[["phi"]][1]
  
  if(cov.model == "matern"){
    if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting value list")}
    nu.starting <- starting[["nu"]][1]
  }

  storage.mode(beta.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(sigma.sq.starting) <- "double"
  storage.mode(tau.sq.starting) <- "double"
  storage.mode(nu.starting) <- "double"

  ####################################################
  ##Tuning values
  ####################################################
  phi.tuning <- 0
  sigma.sq.tuning <- 0
  tau.sq.tuning <- 0
  nu.tuning <- 0
  
  if(missing(tuning)){stop("error: tuning value vector for the spatial parameters must be specified")}
  
  names(tuning) <- tolower(names(tuning))
  
  if(!"sigma.sq" %in% names(tuning)){stop("error: sigma.sq must be specified in tuning value list")}
  sigma.sq.tuning <- tuning[["sigma.sq"]][1]
  
  if(nugget){
    if(!"tau.sq" %in% names(tuning)){stop("error: tau.sq must be specified in tuning value list")}
    tau.sq.tuning <- tuning[["tau.sq"]][1]
  }
  
  if(!"phi" %in% names(tuning)){stop("error: phi must be specified in tuning value list")}
  phi.tuning <- tuning[["phi"]][1]
    
  if(cov.model == "matern"){
    if(!"nu" %in% names(tuning)){stop("error: nu must be specified in tuning value list")}
    nu.tuning <- tuning[["nu"]][1]
  }    
    
  storage.mode(phi.tuning) <- "double"
  storage.mode(sigma.sq.tuning) <- "double"
  storage.mode(tau.sq.tuning) <- "double"
  storage.mode(nu.tuning) <- "double"
  
  ####################################################
  ##Other stuff
  ####################################################
  storage.mode(n.report) <- "integer"
  storage.mode(verbose) <- "integer"

  ####################################################
  ##Fix for no nugget pp model
  ####################################################
  if(is.pp && !nugget){
    
    if(!modified.pp){
      modified.pp <- TRUE
      warning("Because the specified predictive process model has no nugget, a modified predictive process was used.")
    }

    if(any(knots.coords.D == 0)){

      jit <- 1000*.Machine$double.eps
      warning(paste("Because the specified predictive process model has no nugget, a shift of distance ",jit," was added to those ", sum(knots.coords.D == 0), " knot locations that coincide with observed locations.",sep=""))
      
      while(any(knots.coords.D == 0)){
        ind <- which(knots.coords.D == 0, arr.ind = TRUE)
        knot.coords[ind[1,]] <- knot.coords[ind[1,]]+jit       
        knots.coords.D <- iDist(knot.coords, coords)
        knots.D <- iDist(knot.coords)
      }
    }
     
  }

  ####################################################
  ##Pack it up and off it goes
  ####################################################
  ptm <- proc.time()
  
  if(is.pp){
    out <- .Call("spPPLM", Y, X, p, n, m, knots.D, knots.coords.D,
                 modified.pp, beta.prior, beta.Norm, sigma.sq.IG, tau.sq.IG, nu.Unif, phi.Unif,
                 beta.starting, phi.starting, sigma.sq.starting, tau.sq.starting, nu.starting,
                 phi.tuning, sigma.sq.tuning, tau.sq.tuning, nu.tuning,
                 nugget, cov.model, is.amcmc, n.batch, batch.length, accept.rate, verbose, n.report)
    
  }else{
    out <- .Call("spLM", Y, X, p, n, coords.D,
                 beta.prior, beta.Norm, sigma.sq.IG, tau.sq.IG, nu.Unif, phi.Unif,
                 phi.starting, sigma.sq.starting, tau.sq.starting, nu.starting,
                 phi.tuning, sigma.sq.tuning, tau.sq.tuning, nu.tuning,
                 nugget, cov.model, is.amcmc, n.batch, batch.length, accept.rate, verbose, n.report)   
  }
  
  run.time <- proc.time() - ptm
  
  ##parameter names
  if(is.pp){
    out$p.beta.samples <- mcmc(t(out$p.beta.samples))
    colnames(out$p.beta.samples) <- x.names
  }
  
  out$p.theta.samples <- mcmc(t(out$p.theta.samples))
  
  col.names <- rep("null",ncol(out$p.theta.samples))
  
  if(!nugget && cov.model != "matern"){
    col.names[1:2] <- c("sigma.sq", "phi")
  }else if(nugget && cov.model != "matern"){
    col.names[1:3] <- c("sigma.sq", "tau.sq", "phi")
  }else if(!nugget && cov.model == "matern"){
    col.names[1:3] <- c("sigma.sq", "phi", "nu")
  }else{
    col.names[1:4] <- c("sigma.sq", "tau.sq", "phi", "nu")
  }
  
  colnames(out$p.theta.samples) <- col.names
  
  out$Y <- Y
  out$X <- X
  out$coords <- coords
  out$is.pp <- is.pp
  if(is.pp){
    out$knot.coords <- knot.coords
    out$modified.pp <- modified.pp
  }  
  out$cov.model <- cov.model
  out$nugget <- nugget
  out$beta.prior <- beta.prior
  out$beta.Norm <- beta.Norm
  out$x.names <- x.names
  out$run.time <- run.time
  
  class(out) <- "spLM"
  out  
}

