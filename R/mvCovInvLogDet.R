mvCovInvLogDet <- function(coords, knots, cov.model, V, Psi, theta, modified.pp=TRUE, SWM=TRUE, ...){
  
  
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
  is.pp <- TRUE
  if(missing(knots)){
    is.pp <- FALSE
  }

  if(is.pp){
    if(!is.matrix(knots) || ncol(knots) != 2){
      stop("error: knots is misspecified")
    }
  }

  ####################################################
  ##Covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}
  
  if(is.pp){
    
    n <- nrow(coords)
    m <- nrow(V)
    q <- nrow(knots)
    nq <- n*q
    mn <- m*n
    qm <- m*q
    
    knots.D <- as.matrix(dist(knots))
    coords.D <- as.matrix(dist(coords))
    coords.knots.D <- matrix(0, n, q)
    for(i in 1:n){
      coords.knots.D[i,] <- sqrt((knots[,1]-coords[i,1])^2+
                                 (knots[,2]-coords[i,2])^2)
    }
    C.inv <- matrix(0,mn,mn)
    C.str <- matrix(0,qm,qm)
    ct <- matrix(0,mn,qm)
    E <- rep(0, n*m*m)
    
    tmp.mm <- matrix(0,m,m)
    tmp.mm1 <- matrix(0,m,m)
    tmp.mm2 <- matrix(0,m,m)
    tmp.nmqm <- matrix(0,mn,qm)
    tmp.nmqm1 <- matrix(0,mn,qm)
    tmp.qmqm <- matrix(0,qm,qm)
    tmp.qmqm1 <- matrix(0,qm,qm)
    brute <- !SWM
    
    lwork <- n*m*m;
    if(n*m*m < q*m)
      lwork <- q*m;
    
    work <- rep(0,lwork)
    
    storage.mode(m) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(q) <- "integer"
    storage.mode(knots.D) <- "double"
    storage.mode(coords.D) <- "double"
    storage.mode(coords.knots.D) <- "double"
    storage.mode(C.inv) <- "double"
    storage.mode(C.str) <- "double"
    storage.mode(ct) <- "double"
    storage.mode(E) <- "double"
    storage.mode(V) <- "double"
    storage.mode(Psi) <- "double"
    storage.mode(theta) <- "double"
    storage.mode(tmp.mm) <- "double"
    storage.mode(tmp.mm1) <- "double"
    storage.mode(tmp.mm2) <- "double"
    storage.mode(tmp.nmqm) <- "double"
    storage.mode(tmp.nmqm1) <- "double"
    storage.mode(tmp.qmqm) <- "double"
    storage.mode(tmp.qmqm1) <- "double"
    storage.mode(lwork) <- "integer"
    storage.mode(work) <- "double"
    storage.mode(brute) <- "integer"
    
    if(modified.pp){
      log.det <-  .Call("mvMPPCovInvDet_wrap", knots.D, coords.knots.D, C.inv, C.str, ct, 
                        E, n, m, q, Psi, V, theta, 
                        tmp.mm, tmp.mm1, tmp.mm2, tmp.nmqm, 
                        tmp.nmqm1, tmp.qmqm, tmp.qmqm1, lwork, work,
                        cov.model, brute)
      
      out <- list()
      out$log.det <- log.det
      out$C.inv <- C.inv
      
      if(brute){
        out$C.inv[upper.tri(C.inv,F)] <- t(C.inv)[upper.tri(C.inv,F)]
      }
      
      out$C <- .Call("mvMPPCov", knots.D, coords.knots.D, n, m, q, Psi, V, theta, cov.model)
      
    }else{   
      log.det <-  .Call("mvPPCovInvDet_wrap", knots.D, coords.knots.D, C.inv, C.str, ct, 
                        E, n, m, q, Psi, V, theta, 
                        tmp.mm, tmp.mm1, tmp.mm2, tmp.nmqm, 
                        tmp.nmqm1, tmp.qmqm, tmp.qmqm1, lwork, work,
                        cov.model, brute)
      
      out <- list()
      out$log.det <- log.det
      out$C.inv <- C.inv
      
      if(brute){
        out$C.inv[upper.tri(C.inv,F)] <- t(C.inv)[upper.tri(C.inv,F)]
      }
      
      out$C <- .Call("mvPPCov", knots.D, coords.knots.D, n, m, q, Psi, V, theta, cov.model)
    }
    
    return(out);
    
  }else{##non-pp
    
    n <- nrow(coords)
    m <- nrow(V)
    mn <- m*n
    coords.D <- as.matrix(dist(coords))
    C.inv <- matrix(0,mn,mn)

    tmp.mm <- matrix(0,m,m)
    tmp.mm1 <- matrix(0,m,m)
    tmp.mm2 <- matrix(0,m,m)
    
    storage.mode(m) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(coords.D) <- "double"
    storage.mode(C.inv) <- "double"
    storage.mode(V) <- "double"
    storage.mode(Psi) <- "double"
    storage.mode(theta) <- "double"
    storage.mode(tmp.mm) <- "double"
    storage.mode(tmp.mm1) <- "double"
    storage.mode(tmp.mm2) <- "double"

    log.det <-  .Call("mvCovInvDet_wrap", coords.D, C.inv, n, m, Psi, V, theta, 
                      tmp.mm, tmp.mm1, tmp.mm2, cov.model)
        
    out <- list()
    out$log.det <- log.det
    out$C.inv <- C.inv
    out$C.inv[upper.tri(C.inv,F)] <- t(C.inv)[upper.tri(C.inv,F)]

    out$C <- .Call("mvCov", coords.D, n, m, Psi, V, theta, cov.model)

    return(out);
  }
}
