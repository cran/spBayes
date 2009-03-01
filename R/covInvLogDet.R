covInvLogDet <- function(coords, knots, cov.model, sigma.sq, tau.sq, theta, modified.pp=TRUE, SWM=TRUE, ...){
  
  
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
  ##Distance matrices
  ####################################################
  
  ####################
  ##Coords
  ####################
  if(missing(coords)){stop("error: coords must be specified")}
  if(!is.matrix(coords) || ncol(coords) != 2){
    stop("error: coords is misspecified")
  }

  ####################
  ##Knots
  ####################
  if(missing(knots)){stop("error: knots must be specified")}
  if(!is.matrix(knots) || ncol(knots) != 2){
    stop("error: knots is misspecified")
  }

  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}

  if(modified.pp){
    
    n <- nrow(coords)
    m <- nrow(knots)
    knots.D <- as.matrix(dist(knots))
    coords.D <- as.matrix(dist(coords))
    coords.knots.D <- matrix(0, n, m)
    for(i in 1:n){
      coords.knots.D[i,] <- sqrt((knots[,1]-coords[i,1])^2+
                                 (knots[,2]-coords[i,2])^2)
    }
    
    C.inv <- matrix(0,n,n)
    C.str <- matrix(0,m,m)
    ct <- matrix(0,n,m)
    E <- matrix(0,n,1)
    E.inv <- matrix(0,n,1)
    tmp.mm <- matrix(0,m,m)
    tmp.nm <- matrix(0,n,m)
    tmp.nm1 <-matrix(0,n,m)
    tmp.nn <- matrix(0,n,n)
    brute <- !SWM
    
    storage.mode(m) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(knots) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(knots.D) <- "double"
    storage.mode(coords.D) <- "double"
    storage.mode(coords.knots.D) <- "double"
    storage.mode(C.inv) <- "double"
    storage.mode(C.str) <- "double"
    storage.mode(ct) <- "double"
    storage.mode(E) <- "double"
    storage.mode(E.inv) <- "double"
    storage.mode(tau.sq) <- "double"
    storage.mode(sigma.sq) <- "double"
    storage.mode(theta) <- "double"
    storage.mode(tmp.mm) <- "double"
    storage.mode(tmp.nm) <- "double"
    storage.mode(tmp.nm1) <- "double"
    storage.mode(tmp.nn) <- "double"
    storage.mode(brute) <- "integer"

    log.det <- .Call("mPPCovInvDet_wrap", knots.D, coords.knots.D, C.inv, C.str, ct,
                     E, E.inv, n, m, tau.sq, sigma.sq, theta, tmp.mm, tmp.nm, tmp.nm1, tmp.nn,
                     cov.model, brute)

    out <- list()
    out$log.det <- log.det
    out$C.inv <- C.inv
    
    if(brute){
      out$C.inv[lower.tri(C.inv,F)] <- t(C.inv)[lower.tri(C.inv,F)]
    }

    out$C <- .Call("mPPCov", knots.D, coords.knots.D, n, m, tau.sq, sigma.sq, theta, cov.model)
    
    return(out);
    
  }else{
    
    n <- nrow(coords)
    m <- nrow(knots)
    knots.D <- as.matrix(dist(knots))
    coords.D <- as.matrix(dist(coords))
    coords.knots.D <- matrix(0, n, m)
    for(i in 1:n){
      coords.knots.D[i,] <- sqrt((knots[,1]-coords[i,1])^2+
                                 (knots[,2]-coords[i,2])^2)
    }
    
    C.inv <- matrix(0,n,n)
    C.str <- matrix(0,m,m)
    ct <- matrix(0,n,m)
    tmp.mm <- matrix(0,m,m)
    tmp.nm <- matrix(0,n,m)
    brute <- !SWM
    
    storage.mode(m) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(knots) <- "double"
    storage.mode(coords) <- "double"
    storage.mode(knots.D) <- "double"
    storage.mode(coords.D) <- "double"
    storage.mode(coords.knots.D) <- "double"
    storage.mode(C.inv) <- "double"
    storage.mode(C.str) <- "double"
    storage.mode(ct) <- "double"
    storage.mode(tau.sq) <- "double"
    storage.mode(sigma.sq) <- "double"
    storage.mode(theta) <- "double"
    storage.mode(tmp.mm) <- "double"
    storage.mode(tmp.nm) <- "double"
    storage.mode(brute) <- "integer"
    
    log.det <- .Call("pPCovInvDet_wrap", knots.D, coords.knots.D, C.inv, C.str, ct,
                     n, m, tau.sq, sigma.sq, theta, tmp.mm, tmp.nm, cov.model, brute)

    out <- list()
    out$log.det <- log.det
    out$C.inv <- C.inv
    
    if(brute){
      out$C.inv[lower.tri(C.inv,F)] <- t(C.inv)[lower.tri(C.inv,F)]
    }

    out$C <- .Call("pPCov", knots.D, coords.knots.D, n, m, tau.sq, sigma.sq, theta, cov.model)
    
    return(out);
    
  }

}
