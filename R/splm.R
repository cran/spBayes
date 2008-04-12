sp.lm <- function(formula, data = parent.frame(), coords, knots,
                  fixed, starting, sp.tuning, priors, cov.model, sp.effects=TRUE, n.samples, verbose=TRUE, n.report=100, ...){

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
  ##formula
  ####################################################
  if(missing(formula)){stop("error: formula must be specified")}
  
  if(class(formula) == "formula"){
    
    holder <- parse.formula(formula, data)
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
  storage.mode(sp.effects) <- "integer"

  ####################################################
  ##Distance matrices
  ####################################################
  
  ####################
  ##Coords
  ####################
  if(missing(coords)){stop("error: coords must be specified")}
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(ncol(coords) != 2 || nrow(coords) != n){
    stop("error: either the coords have more than two columns or then number of rows is different than
          data used in the model formula")
  }
  
  coords.D <- as.matrix(dist(coords))
  storage.mode(coords.D) <- "double"
  
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
  knots.D <- 0
  coords.knots.D <- 0
  
  if(is.pp){
    knots.D <- as.matrix(dist(knot.coords))
    m <- nrow(knots.D)
    coords.knots.D <- matrix(0, n, m) ##this is for c^t

    for(i in 1:n){
      coords.knots.D[i,] <- sqrt((knot.coords[,1]-coords[i,1])^2+
                                 (knot.coords[,2]-coords[i,2])^2)
    }

    storage.mode(is.pp) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(knots.D) <- "double"
    storage.mode(coords.knots.D) <- "double"
  }

  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}

  ####################################################
  ##Fixed parameters
  ####################################################
  beta.fixed <- FALSE
  sigma.sq.fixed <- FALSE
  tau.sq.fixed <- FALSE
  phi.fixed <- FALSE
  nu.fixed <- FALSE

  if(!missing(fixed)){
    names(fixed) <- tolower(names(fixed))
    
    if("beta" %in% names(fixed))
      beta.fixed <- TRUE
    
    if("sigma.sq" %in% names(fixed))
      sigma.sq.fixed <- TRUE

    if("tau.sq" %in% names(fixed))
      tau.sq.fixed <- TRUE

    if("phi" %in% names(fixed))
      phi.fixed <- TRUE

    if("nu" %in% names(fixed))
      nu.fixed <- TRUE
  }
  
  storage.mode(beta.fixed) <- "integer"
  storage.mode(sigma.sq.fixed) <- "integer"
  storage.mode(tau.sq.fixed) <- "integer"
  storage.mode(phi.fixed) <- "integer"
  storage.mode(nu.fixed) <- "integer"

  ####################################################
  ##Starting values
  ####################################################

  ####################
  ##Free parameters
  ####################
  all.fixed <- FALSE

  if(cov.model == "matern"){
    if(all(sigma.sq.fixed, tau.sq.fixed, phi.fixed, nu.fixed))
      all.fixed <- TRUE
  }
  
  if(cov.model != "matern"){
    if(all(sigma.sq.fixed, tau.sq.fixed, phi.fixed))
      all.fixed <- TRUE
  }
    
  beta.starting <- 0
  sigma.sq.starting <- 0
  tau.sq.starting <- 0
  phi.starting <- 0
  nu.starting <- 0

  nugget <- FALSE
  if(tau.sq.fixed)
    nugget <- TRUE
  
  if(!all.fixed){
    
    if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
    
    names(starting) <- tolower(names(starting))   

    if(!sigma.sq.fixed && !"sigma.sq" %in% names(starting)){stop("error: sigma.sq must be specified as fixed or free")}
    
    if("tau.sq" %in% names(starting))
      nugget <- TRUE
    
    if(!phi.fixed && !"phi" %in% names(starting)){stop("error: phi must be specified as fixed or free")}
    
    if(cov.model == "matern"){
      if(!nu.fixed && !"nu" %in% names(starting)){stop("error: nu must be specified as fixed or free")}
    }


    if("beta" %in% names(starting)) beta.starting <- starting[["beta"]]
    if("sigma.sq" %in% names(starting)) sigma.sq.starting <- starting[["sigma.sq"]][1]
    if(nugget && "tau.sq" %in% names(starting)) tau.sq.starting <- starting[["tau.sq"]][1]
    if("phi" %in% names(starting)) phi.starting <- starting[["phi"]][1]
    if(cov.model == "matern" && "nu" %in% names(starting)) nu.starting <- starting[["nu"]][1]
  }

  ##Special case for beta, if missing beta in fixed and free, then just assume free
  if(!beta.fixed){
    if(!missing(starting)){
      names(starting) <- tolower(names(starting))
      if("beta" %in% names(starting)){
        if(!all(is.vector(starting[["beta"]]), length(starting[["beta"]]) == p)){stop("error: beta must be a vector of length p")}
        beta.starting <- starting[["beta"]]
      }else{
        beta.starting <- unname(coefficients(lm(Y~X-1))) ##assumed free and no starting values given
      }
    }else{
      beta.starting <- unname(coefficients(lm(Y~X-1))) ##assumed free and no starting values given
      }
  }else{
    if(!all(is.vector(fixed[["beta"]]), length(fixed[["beta"]]) == p)){stop("error: beta must be a vector of length p")}
    beta.starting <- fixed[["beta"]]
  }

  if(sigma.sq.fixed) sigma.sq.starting <- fixed[["sigma.sq"]][1]
  if(tau.sq.fixed && nugget) tau.sq.starting <- fixed[["tau.sq"]][1]
  if(phi.fixed) phi.starting <- fixed[["phi"]][1]
  if(nu.fixed && cov.model == "matern") nu.starting <- fixed[["nu"]][1]

  
  storage.mode(nugget) <- "integer"
  storage.mode(beta.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(sigma.sq.starting) <- "double"
  storage.mode(tau.sq.starting) <- "double"
  storage.mode(nu.starting) <- "double"

  ####################################################
  ##Priors
  ####################################################

  ##Priors

  sigma.sq.IG <- 0
  tau.sq.IG <- 0
  nu.Unif <- 0
  phi.Unif <- 0

  if(!all.fixed){

    if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
    names(priors) <- tolower(names(priors))
    
    if(!sigma.sq.fixed){
      if(!"sigma.sq.ig" %in% names(priors)){stop("error: sigma.sq.IG must be specified")}
      sigma.sq.IG <- priors[["sigma.sq.ig"]]
      
      if(!is.vector(sigma.sq.IG) || length(sigma.sq.IG) != 2){stop("error: sigma.sq.IG must be a vector of length 2")}
      if(any(sigma.sq.IG <= 0)){stop("error: sigma.sq.IG must be a positive vector of length 2")}
    }
    
    if(nugget && !tau.sq.fixed){
      if(!"tau.sq.ig" %in% names(priors)){stop("error: tau.sq.IG must be specified")}
      tau.sq.IG <- priors[["tau.sq.ig"]]

      
      if(!is.vector(tau.sq.IG) || length(tau.sq.IG) != 2){stop("error: tau.sq.IG must be a vector of length 2")}
      if(any(tau.sq.IG <= 0)){stop("error: tau.sq.IG must be a positive vector of length 2")}
    }
    
    if(!phi.fixed){
      if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified")}
      phi.Unif <- priors[["phi.unif"]]
      
      if(!is.vector(phi.Unif) || length(phi.Unif) != 2){stop("error: phi.Unif must be a vector of length 2")}
      if(any(phi.Unif <= 0, phi.Unif[1] >= phi.Unif[2])){stop("error: phi.Unif must be a positive vector of length 2 with element 1 < element 2")}
    }
    
    if(!nu.fixed && cov.model == "matern"){

      if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified")}
      nu.Unif <- priors[["nu.unif"]]
      
      if(!is.vector(nu.Unif) || length(nu.Unif) != 2){stop("error: nu.Unif must be a vector of length 2")}
      if(any(nu.Unif <= 0, nu.Unif[1] >= nu.Unif[2])){stop("error: nu.Unif must be a positive vector of length 2 with element 1 < element 2")}
    }
  }
  
  storage.mode(sigma.sq.IG) <- "double"
  storage.mode(tau.sq.IG) <- "double"
  storage.mode(nu.Unif) <- "double"
  storage.mode(phi.Unif) <- "double"

  ####################################################
  ##Tuning values
  ####################################################

  ####################
  ##Tuning params
  ####################

  phi.tuning <- 0
  sigma.sq.tuning <- 0
  tau.sq.tuning <- 0
  nu.tuning <- 0
  
  if(!all.fixed){
    if(missing(sp.tuning)){stop("error: sp.tuning value vector for the spatial parameters must be specified")}
    
    names(sp.tuning) <- tolower(names(sp.tuning))

    if(!sigma.sq.fixed){
      if(!"sigma.sq" %in% names(sp.tuning)){stop("error: sigma.sq must be specified in tuning value list")}
      sigma.sq.tuning <- sp.tuning[["sigma.sq"]][1]
    }
    
    if(nugget && !tau.sq.fixed){
      if(!"tau.sq" %in% names(sp.tuning)){stop("error: tau.sq must be specified in tuning value list")}
      tau.sq.tuning <- sp.tuning[["tau.sq"]][1]
    }
    
    if(!phi.fixed){
      if(!"phi" %in% names(sp.tuning)){stop("error: phi must be specified in tuning value list")}
      phi.tuning <- sp.tuning[["phi"]][1]
    }
    
    if(!nu.fixed && cov.model == "matern"){
      if(!"nu" %in% names(sp.tuning)){stop("error: nu must be specified in tuning value list")}
      nu.tuning <- sp.tuning[["nu"]][1]
    }    
  }
  
  storage.mode(phi.tuning) <- "double"
  storage.mode(sigma.sq.tuning) <- "double"
  storage.mode(tau.sq.tuning) <- "double"
  storage.mode(nu.tuning) <- "double"

  
  ####################################################
  ##Other stuff
  ####################################################
  if(missing(n.samples)){stop("error: n.samples need to be specified")}
  storage.mode(n.samples) <- "integer"
  storage.mode(n.report) <- "integer"
  storage.mode(verbose) <- "integer"


  ####################################################
  ##Pack it up and off it goes
  ####################################################
  
  out <- .Call("splm", Y, X, p, n, coords.D,
               is.pp, m, knots.D, coords.knots.D, nugget, 
               beta.fixed,sigma.sq.fixed,tau.sq.fixed,phi.fixed,nu.fixed,
               sigma.sq.IG, tau.sq.IG, nu.Unif, phi.Unif,
               phi.starting, sigma.sq.starting, tau.sq.starting, nu.starting, beta.starting,
               phi.tuning, sigma.sq.tuning, tau.sq.tuning, nu.tuning, 
               cov.model, n.samples, verbose, n.report, sp.effects)

  out$coords <- coords
  out$is.pp <- is.pp
  
  if(is.pp)
    out$knot.coords <- knot.coords
  
  out$p.samples <- mcmc(t(out$p.samples))

  out$Y <- Y
  out$X <- X
  out$n <- n
  out$m <- m
  out$p <- p
  out$knots.D <- knots.D
  out$coords.D <- coords.D
  out$coords.knots.D <- coords.knots.D
  out$cov.model <- cov.model
  out$nugget <- nugget
  out$verbose <- verbose
  out$n.samples <- n.samples
  out$recovered.effects <- sp.effects
  
  col.names <- rep("null",ncol(out$p.samples))
  
  col.names[1:p] <- x.names
  if(!nugget && cov.model != "matern"){
    col.names[(p+1):(p+2)] <- c("sigma.sq", "phi")
  }else if(nugget && cov.model != "matern"){
    col.names[(p+1):(p+3)] <- c("sigma.sq", "tau.sq", "phi")
  }else if(!nugget && cov.model == "matern"){
    col.names[(p+1):(p+3)] <- c("sigma.sq", "phi", "nu")
  }else{
    col.names[(p+1):(p+4)] <- c("sigma.sq", "tau.sq", "phi", "nu")
  }
    
  colnames(out$p.samples) <- col.names
  

  class(out) <- "sp.lm"
  out  
}

