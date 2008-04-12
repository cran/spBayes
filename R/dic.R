sp.DIC <- function(sp.obj, DIC.marg=TRUE, DIC.unmarg=TRUE, start=1, end, thin=1, verbose=TRUE, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  
  if(missing(sp.obj)){stop("error: sp.DIC expects sp.obj\n")}
  if(!class(sp.obj) %in% c("ggt.sp","sp.lm","bayes.geostat.exact")){
    stop("error: sp.DIC requires an output object of class ggt.sp, sp.lm, or bayes.geostat.exact\n")}
  if(!is.logical(DIC.marg)){stop("error: DIC.marg must be of type logical\n")}
  if(!is.logical(DIC.unmarg)){stop("error: DIC.unmarg must be of type logical\n")}
  if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}

  
  if(class(sp.obj) == "sp.lm"){

    is.pp <- sp.obj$is.pp
    
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
    beta <- t(sp.obj$p.samples[,1:p])##transpose to simply BLAS call
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
   
    if(is.pp && !nugget){ ##need to ridge
      tau.sq <- rep(1e-10, n.samples)
      nugget <- TRUE;
    }
  
    storage.mode(X) <- "double"
    storage.mode(Y) <- "double"
    storage.mode(is.pp) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(m) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(nugget) <- "integer"
    storage.mode(beta) <- "double"
    storage.mode(sigma.sq) <- "double"
    storage.mode(tau.sq) <- "double"
    storage.mode(phi) <- "double"
    storage.mode(nu) <- "double"
    storage.mode(obs.D) <- "double"
    storage.mode(obs.knots.D) <- "double"
    storage.mode(knots.D) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(w) <- "double"
    storage.mode(w.str) <- "double"
    storage.mode(sp.effects) <- "integer"
    storage.mode(DIC.marg) <- "integer"
    storage.mode(DIC.unmarg) <- "integer"
    
    out <- .Call("splmDIC",X, Y, is.pp, n, m, p, nugget, beta, sigma.sq, tau.sq, phi, nu,
                 obs.D, obs.knots.D, knots.D, cov.model, n.samples, w, w.str, sp.effects,
                 DIC.marg, DIC.unmarg, verbose)
    out
    
    
  }else if(class(sp.obj) == "bayes.geostat.exact"){
    
    X <- sp.obj$args$X
    n <- sp.obj$args$n
    p <- sp.obj$args$p
    Y <- sp.obj$args$Y
    coords <- sp.obj$args$coords
    cov.model <- sp.obj$args$cov.model
    samples <- sp.obj$p.samples
    phi <- sp.obj$args$phi
    n.samples <- sp.obj$args$n.samples
    alpha <- sp.obj$args$alpha
    
    if(cov.model == "matern")
      nu <- sp.obj$args$nu
    
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
      stop("error: in sp.predict, specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
    }

    R.eigen <- eigen(R)
    R.vals <- R.eigen$values
    R.vecs <- R.eigen$vectors
    R.vects.t <- t(R.vecs)

    d <- rep(0, n.samples)
    
    ##
    ##marginalized
    ##
    
    marg.DIC <- matrix(0,4,1)
    if(DIC.marg){
      
      for(s in 1:n.samples){
        R.inv <- R.vecs%*%diag(1/(sigma.sq[s]*R.vals+tau.sq[s]))%*%R.vects.t
        R.log.det <- sum(log(sigma.sq[s]*R.vals+tau.sq[s]))
        d[s] <- R.log.det+(t(Y-X%*%beta[s,])%*%R.inv%*%(Y-X%*%beta[s,]))
      }
      
      d.bar <- mean(d)
      
      sigma.sq.mu <- mean(sigma.sq)
      tau.sq.mu <- mean(tau.sq)
      beta.mu <- as.matrix(colMeans(beta))
      
      R.inv <- R.vecs%*%diag(1/(sigma.sq.mu*R.vals+tau.sq.mu))%*%R.vects.t
      R.log.det <- sum(log(sigma.sq.mu*R.vals+tau.sq.mu))
      
      d.bar.omega <- R.log.det+(t(Y-X%*%beta.mu)%*%R.inv%*%(Y-X%*%beta.mu))
      pd <- d.bar - d.bar.omega
      dic <- d.bar - pd
      
      
      rownames(marg.DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
      marg.DIC[1,1] <- d.bar
      marg.DIC[2,1] <- d.bar.omega
      marg.DIC[3,1] <- pd
      marg.DIC[4,1] <- dic
    }
    
    ##
    ##unmarginalized
    ##

    ##
    ##recover spatial effects
    ##
    unmarg.DIC <- matrix(0,4,1)
    
    if(DIC.unmarg){
      cat(paste("\nRecovering spatial effects ... \n", sep=""))
      
##      w <- matrix(0, n, n.samples)
##      
##      for(s in 1:n.samples){
##        
##        S.w <- R.vecs%*%diag(1/(sigma.sq[s]*R.vals+tau.sq[s]))%*%R.vects.t
##        
##        S.mu <- S.w%*%(Y - X%*%as.matrix(beta[s,]))/tau.sq[s]
##        
##        S.w.sq <- R.vecs%*%diag(sqrt(1/(sigma.sq[s]*R.vals+tau.sq[s])))
##        
##        w[,s] <- S.w.sq%*%as.matrix(rnorm(n))+S.mu
##      }
      
      w <- matrix(0, n, n.samples)
      
      R.inv <- chol2inv(chol(R))
      
      V.sp <- chol2inv(chol(R.inv + (1/alpha)*diag(nrow(R))))
      resid.posterior <- matrix(rep(Y, times=n.samples), nrow=length(Y), ncol=n.samples) -  X%*%t(beta) 
      sp.posterior.mean <- (1/alpha)*t(V.sp%*%resid.posterior)
      
      V.sp.root <- t(chol(V.sp)) ## chol returns "upper-triangular"; so t(); 
            
      for (s in 1:n.samples) {
        
        ## Using rmvnorm is slow - it calculates the matrix square-root each time
        ## sp.effects[s,] <- rmvnorm(1, sp..posterior.mean[s,], sigma.sq.posterior.samples[s]*V.sp)
        
        ## Instead use the pre-computed V.sp.root
        z <- rnorm(nrow(V.sp), 0, 1)
        w[,s] <- sp.posterior.mean[s,] + sqrt(sigma.sq[s])*V.sp.root%*%z
        
      }
         
      
      for(s in 1:n.samples){
        d[s] <- n*log(tau.sq[s])+1/(tau.sq[s])*(t(Y-X%*%beta[s,]-w[,s])%*%(Y-X%*%beta[s,]-w[,s]))
      }
      
      d.bar <- mean(d)
      
      sigma.sq.mu <- mean(sigma.sq)
      tau.sq.mu <- mean(tau.sq)
      beta.mu <- as.matrix(colMeans(beta))
      w.mu <- as.matrix(rowMeans(w))
      
      d.bar.omega <- n*log(tau.sq.mu)+1/(tau.sq.mu)*(t(Y-X%*%beta.mu-w.mu)%*%(Y-X%*%beta.mu-w.mu))
      pd <- d.bar - d.bar.omega
      dic <- d.bar - pd
      
      rownames(unmarg.DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
      unmarg.DIC[1,1] <- d.bar
      unmarg.DIC[2,1] <- d.bar.omega
      unmarg.DIC[3,1] <- pd
      unmarg.DIC[4,1] <- dic
    }

    out <- list()
    if(DIC.marg)
      out$DIC.marg <- marg.DIC

    if(DIC.unmarg){
      out$DIC.unmarg <- unmarg.DIC
      out$sp.effects <- w
    }

    out

  }else if(class(sp.obj) == "ggt.sp"){
    
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
    m <- sp.obj$m
    X <- as.matrix(sp.obj$X)
    Y <- as.matrix(sp.obj$Y)
    
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
    K.means <- rowMeans(K)
    
    ##Psi
    no.Psi <- sp.obj$no.Psi
    Psi <- NULL
    Psi.case <- NULL
    Psi.means <- NULL
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
      Psi.means <- rowMeans(Psi)
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
    phi.means <- rowMeans(phi)
    
    
    cov.model <- sp.obj$cov.model
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
    obs.coords <- as.matrix(sp.obj$coords)
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
  }else{
    stop("error: sp.DIC requires an output object of class ggt.sp or bayes.geostat.exact\n")
  }
}
