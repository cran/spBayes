spMvLM <- function(formula, data = parent.frame(), coords, knots,
                   starting, sp.tuning, priors, cov.model, 
                   modified.pp = TRUE, n.samples, sub.samples, verbose=TRUE, n.report=100, ...){
  
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
    
    holder <- parseFormula(formula, data)
    Y <- holder[[1]]
    X <- holder[[2]]
    x.names <- holder[[3]]
    m <- 1
    
  }else if(class(formula) == "list"){
    
    mv.mats <- mkMats(formula, data)
    Y <- mv.mats[[1]]
    X <- mv.mats[[2]]
    x.names <- mv.mats[[3]]
    m <- length(formula)
    
  }else{
    stop("error: formula is misspecified")
  }

  p <- ncol(X)
  n <- nrow(X)/m

  ##make sure storage mode is correct
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(m) <- "integer"
  storage.mode(p) <- "integer"
  storage.mode(n) <- "integer"

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

  q <- 0
  knots.D <- 0
  coords.knots.D <- 0
  
  if(is.pp){
    knots.D <- as.matrix(dist(knot.coords))
    q <- nrow(knots.D)
    coords.knots.D <- matrix(0, n, q) ##this is for c^t

    for(i in 1:n){
      coords.knots.D[i,] <- sqrt((knot.coords[,1]-coords[i,1])^2+
                                 (knot.coords[,2]-coords[i,2])^2)
    }

    storage.mode(modified.pp) <- "integer"
    storage.mode(q) <- "integer"
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
  ##Starting values
  ####################################################

  beta.starting <- 0
  A.starting <- 0
  L.starting <- 0
  phi.starting <- 0
  nu.starting <- 0

  nugget <- FALSE
  
  if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
  
  names(starting) <- tolower(names(starting))   

  if(!"beta" %in% names(starting)){
    beta.starting <- coefficients(lm(Y~X-1))
  }else{
    beta.starting <- starting[["beta"]]
  }
     
  if(!"a" %in% names(starting)){stop("error: A must be specified in starting")}
  A.starting <- starting[["a"]]
  if(length(A.starting) != m*(m-1)/2+m){stop(paste("error: A must be of length ",m*(m-1)/2+m," in starting value list",sep=""))}
  
  if("l" %in% names(starting)){
    L.starting <- starting[["l"]]
    if(length(L.starting) != m*(m-1)/2+m){stop(paste("error: L must be of length ",m*(m-1)/2+m," in starting value list",sep=""))}
     nugget <- TRUE
  }  

  if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting")}
  phi.starting <- starting[["phi"]]
  if(length(phi.starting) != m){stop(paste("error: phi must be of length ",m," in starting value list",sep=""))}
  
  if(cov.model == "matern"){
    if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting")}
    nu.starting <- starting[["nu"]]
    if(length(nu.starting) != m){stop(paste("error: nu must be of length ",m," in starting value list",sep=""))}
  }
              
  storage.mode(nugget) <- "integer"
  storage.mode(beta.starting) <- "double"
  storage.mode(phi.starting) <- "double"
  storage.mode(A.starting) <- "double"
  storage.mode(L.starting) <- "double"
  storage.mode(nu.starting) <- "double"

  ####################################################
  ##Priors
  ####################################################
  K.IW <- 0
  Psi.IW <- 0
  nu.Unif <- 0
  phi.Unif <- 0

  if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
  names(priors) <- tolower(names(priors))
  
  if(!"k.iw" %in% names(priors)){stop("error: K.IW must be specified")}
  K.IW <- priors[["k.iw"]]
  if(!is.list(K.IW) || length(K.IW) != 2){stop("error: K.IW must be a list of length 2")}
  if(length(K.IW[[1]]) != 1 ){stop("error: K.IW[[1]] must be of length 1 (i.e., the IW df hyperparameter)")}   
  if(length(K.IW[[2]]) != m^2 ){stop(paste("error: K.IW[[2]] must be a vector or matrix of length, ",m^2, ", (i.e., the IW scale matrix hyperparameter)",sep=""))}
  
  if(nugget){
    if(!"psi.iw" %in% names(priors)){stop("error: Psi.IW must be specified")}
    Psi.IW <- priors[["psi.iw"]]

    if(!is.list(Psi.IW) || length(Psi.IW) != 2){stop("error: Psi.IW must be a list of length 2")}
    if(length(Psi.IW[[1]]) != 1 ){stop("error: Psi.IW[[1]] must be of length 1 (i.e., the IW df hyperparameter)")}   
    if(length(Psi.IW[[2]]) != m^2 ){stop(paste("error: Psi.IW[[2]] must be a vector or matrix of length, ",m^2, ", (i.e., the IW scale matrix hyperparameter)",sep=""))}
  }
  
  if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified")}
  phi.Unif <- priors[["phi.unif"]]
  
  if(!is.vector(phi.Unif) || length(phi.Unif) != 2*m){stop(paste("error: phi.Unif must be a vector of length ",2*m,sep=""))}
  if(any(phi.Unif[seq(2,2*m,2)]-phi.Unif[seq(1,2*m,2)] <= 0)){stop("error: phi.Unif has zero support")}
  
  if(cov.model == "matern"){
    
    if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified")}
    nu.Unif <- priors[["nu.unif"]]
    
    if(!is.vector(nu.Unif) || length(nu.Unif) != 2*m){stop(paste("error: nu.Unif must be a vector of length ",2*m,sep=""))}
    if(any(nu.Unif[seq(2,2*m,2)]-nu.Unif[seq(1,2*m,2)] <= 0)){stop("error: nu.Unif has zero support")}
  }
  
  #storage.mode(K.IW) <- "double"
  #storage.mode(Psi.IW) <- "double"
  storage.mode(nu.Unif) <- "double"
  storage.mode(phi.Unif) <- "double"

  ####################################################
  ##Tuning values
  ####################################################

  phi.tuning <- 0
  A.tuning <- 0
  L.tuning <- 0
  nu.tuning <- 0
  
  if(missing(sp.tuning)){stop("error: sp.tuning value vector for the spatial parameters must be specified")}
  
  names(sp.tuning) <- tolower(names(sp.tuning))
  
  if(!"a" %in% names(sp.tuning)){stop("error: A must be specified in tuning value list")}
  A.tuning <- as.vector(sp.tuning[["a"]])
  if(length(A.tuning) != m*(m-1)/2+m){stop(paste("error: A must be of length ",m*(m-1)/2+m," in tuning value list",sep=""))}
  
  if(nugget){
    if(!"l" %in% names(sp.tuning)){stop("error: L must be specified in tuning value list")}
    L.tuning <- as.vector(sp.tuning[["l"]])
    if(length(L.tuning) != m*(m-1)/2+m){stop(paste("error: L must be of length ",m*(m-1)/2+m," in tuning value list",sep=""))}
  }
  
  if(!"phi" %in% names(sp.tuning)){stop("error: phi must be specified in tuning value list")}
  phi.tuning <- sp.tuning[["phi"]]
  if(length(phi.tuning) != m){stop(paste("error: phi must be of length ",m," in tuning value list",sep=""))}
    
  if(cov.model == "matern"){
    if(!"nu" %in% names(sp.tuning)){stop("error: nu must be specified in tuning value list")}
    nu.tuning <- sp.tuning[["nu"]]
    if(length(nu.tuning) != m){stop(paste("error: nu must be of length ",m," in tuning value list",sep=""))}
  }    
    
  storage.mode(phi.tuning) <- "double"
  storage.mode(A.tuning) <- "double"
  storage.mode(L.tuning) <- "double"
  storage.mode(nu.tuning) <- "double"
  
  
  ####################################################
  ##Other stuff
  ####################################################
  if(missing(n.samples)){stop("error: n.samples need to be specified")}

  if(missing(sub.samples)){sub.samples <- c(1, n.samples, 1)}
  if(length(sub.samples) != 3 || any(sub.samples > n.samples) ){stop("error: sub.samples misspecified")}
  
  storage.mode(n.samples) <- "integer"
  storage.mode(n.report) <- "integer"
  storage.mode(verbose) <- "integer"


  ####################################################
  ##Pack it up and off it goes
  ####################################################

  ##force modified predictive process if nugget == 0
  if(is.pp){
    if(!nugget && !modified.pp){
      warning("spMvLM: switching to modified predictive process for no nugget model.")
      modified.pp <- TRUE
      storage.mode(modified.pp) <- "integer"
    }
  }
  
  if(is.pp){
    
    out <- .Call("spPPMvLM", Y, X, p, n, m, coords.D,
                 modified.pp, q, knots.D, coords.knots.D,                
                 K.IW, Psi.IW, nu.Unif, phi.Unif,
                 phi.starting, A.starting, L.starting, nu.starting, beta.starting,
                 phi.tuning, A.tuning, L.tuning, nu.tuning, 
                 nugget, cov.model, n.samples, verbose, n.report)
    
  }else{
    
    out <- .Call("spMvLM", Y, X, p, n, m, coords.D,
                 K.IW, Psi.IW, nu.Unif, phi.Unif,
                 phi.starting, A.starting, L.starting, nu.starting, beta.starting,
                 phi.tuning, A.tuning, L.tuning, nu.tuning, 
                 nugget, cov.model, n.samples, verbose, n.report)
    
  }

  out$coords <- coords
  out$is.pp <- is.pp
  out$modified.pp <- modified.pp
  
  if(is.pp){out$knot.coords <- knot.coords}
  
  out$Y <- Y
  out$X <- X
  out$n <- n
  out$m <- m
  out$q <- q
  out$p <- p
  out$knots.D <- knots.D
  out$coords.D <- coords.D
  out$coords.knots.D <- coords.knots.D
  out$cov.model <- cov.model
  out$nugget <- nugget
  out$verbose <- verbose
  #out$n.samples <- n.samples
  out$sub.samples <- sub.samples
  out$recovered.effects <- TRUE ##forced recovery

  ##subsample 
  out$sp.effects <- out$sp.effects[,seq(sub.samples[1], sub.samples[2], by=as.integer(sub.samples[3]))]
  if(is.pp){out$sp.effects.knots <- out$sp.effects.knots[,seq(sub.samples[1], sub.samples[2], by=as.integer(sub.samples[3]))]}
  out$p.samples <- mcmc(t(out$p.samples[,seq(sub.samples[1], sub.samples[2], by=as.integer(sub.samples[3]))]))
  out$n.samples <- nrow(out$p.samples)##get adjusted n.samples
  
  col.names <- rep("null",ncol(out$p.samples))

  nltr <- m*(m-1)/2+m
  
  col.names[1:p] <- x.names
  if(!nugget && cov.model != "matern"){
    col.names[(p+1):ncol(out$p.samples)] <- c(rep("K",nltr), paste("phi_",1:m,sep=""))
  }else if(nugget && cov.model != "matern"){
    col.names[(p+1):ncol(out$p.samples)] <- c(rep("K",nltr), rep("Psi",nltr), paste("phi_",1:m,sep=""))
  }else if(!nugget && cov.model == "matern"){
    col.names[(p+1):ncol(out$p.samples)] <- c(rep("K",nltr), paste("phi.",1:m,sep=""), paste("nu_",1:m,sep=""))
  }else{
    col.names[(p+1):ncol(out$p.samples)] <- c(rep("K",nltr), rep("Psi",nltr), paste("phi_",1:m,sep=""), paste("nu_",1:m,sep=""))
  }
    
  colnames(out$p.samples) <- col.names

  AtA <- function(x, m){
    A <- matrix(0, m, m)
    A[lower.tri(A, diag=TRUE)] <- x
    (A%*%t(A))[lower.tri(A, diag=TRUE)]
  }

  ##K.names <- paste("K_",1:nltr,sep="")
  K.names <- paste("K[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")

  colnames(out$p.samples)[colnames(out$p.samples)%in%"K"] <- K.names
  out$p.samples[,K.names] <- t(apply(out$p.samples[,K.names], 1, AtA, m))
  
  if(nugget){
    ##Psi.names <- paste("Psi_",1:nltr,sep="")
    Psi.names <- paste("Psi[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
    colnames(out$p.samples)[colnames(out$p.samples)%in%"Psi"] <- Psi.names
    out$p.samples[,Psi.names] <- t(apply(out$p.samples[,Psi.names], 1, AtA, m)) 
  }

  class(out) <- "spMvLM"
  out  
}

