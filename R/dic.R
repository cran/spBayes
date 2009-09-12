spDIC <- function(sp.obj, start=1, end, thin=1, verbose=TRUE, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }

  if(missing(sp.obj)){stop("error: spDIC expects sp.obj\n")}
  if(!class(sp.obj) %in% c("spGGT","spLM","bayesGeostatExact","bayesLMConjugate", "bayesLMRef", "spMvLM", "spGLM")){
    stop("error: spDIC requires an output object of class bayesLMConjugate, bayesLMRef, spGGT, spLM, bayesGeostatExact, spMvLM, or spGLM\n")}
  if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}
  
  if(class(sp.obj) == "spGLM"){
    
    X <- sp.obj$X
    n <- sp.obj$n
    p <- sp.obj$p
    Y <- sp.obj$Y
    family <- sp.obj$family
    samples <- sp.obj$p.samples
    w <- sp.obj$sp.effects
    n.samples <- sp.obj$n.samples
    
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
    w <- w[,seq(start, end, by=as.integer(thin))]
    
    n.samples <- nrow(samples)
    
    ##get samples
    beta <- as.matrix(samples[,1:p])
    
    d <- rep(0, n.samples)

    cat("-------------------------------------------------\n\t\tCalculating DIC\n-------------------------------------------------\n")
    status <- 0
    
    for(s in 1:n.samples){
      
      if(family == "poisson"){
        d[s] <- -2*sum(-exp(X%*%beta[s,]+w[,s])+Y*(X%*%beta[s,]+w[,s]))
      }else if(family == "binomial"){
        d[s] <- -2*sum(Y*(X%*%beta[s,]+w[,s])-log(1+exp(X%*%beta[s,]+w[,s])))
      }else{
        stop("error: family is misspecified")
      }
      
      #if(verbose){
      #  if(status == 100){
      #    cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
      #    status <- 0
      #  }
      #  status <- status+1
      #}
    }
    
    d.bar <- mean(d)
    
    beta.mu <- as.matrix(colMeans(beta))
    w.mu <- as.matrix(rowMeans(w))
    
    if(family == "poisson"){
      d.bar.omega <- -2*sum(-exp(X%*%beta.mu+w.mu)+Y*(X%*%beta.mu+w.mu))
    }else if(family == "binomial"){
      d.bar.omega <- -2*sum(Y*(X%*%beta.mu+w.mu)-log(1+exp(X%*%beta.mu+w.mu)))
    }else{
      stop("error: family is misspecified")
    }
    
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd

    DIC <- matrix(0,4,1)
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
    
    out <- list()    
    out$DIC <- DIC
 
    out

  }else if(class(sp.obj) == "spMvLM"){
    
    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp
    
    if(is.pp)
      knot.coords <- sp.obj$knot.coords
    
    Y <- sp.obj$Y
    X <- sp.obj$X
    n <- sp.obj$n
    m <- sp.obj$m##number or response variables
    q <- sp.obj$q##number of knots
    p <- sp.obj$p
    nltr <- m*(m-1)/2+m

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
    w.mean <- NULL
    w.str <- NULL
    w.str.mean <- NULL
    
    if(sp.effects){
      w <- sp.obj$sp.effects[,seq(start, end, by=as.integer(thin))]
      w.mean <- rowMeans(w)
      if(is.pp){
        w.str <- sp.obj$sp.effects.knots[,seq(start, end, by=as.integer(thin))]
        w.str.mean <- rowMeans(w.str)
      }
    }else{##allocate to recieve the computed values
      w <- matrix(0, n*m, n.samples) #no need to compute w for the pp, but might as well on the way.
      w.mean <- rep(0, n*m)
      if(is.pp){
        w.str <- matrix(0, q*m, n.samples)
        w.str.mean <- rep(0, q*m)
      }
    }
    
    ##get samples
    samples.mean <- colMeans(samples)
    
    beta <- t(samples[,1:p])##transpose to simplify BLAS call
    beta.mean <- samples.mean[1:p]

    K <- NULL
    A <- NULL
    Psi <- NULL
    phi <- NULL
    nu <- NULL
    K.mean <- NULL
    A.mean <- NULL
    Psi.mean <- NULL
    phi.mean <- NULL
    nu.mean <- NULL

    A.chol <- function(x, m){
      A <- matrix(0, m, m)
      A[lower.tri(A, diag=TRUE)] <- x
      A[upper.tri(A, diag=FALSE)] <- t(A)[upper.tri(A, diag=FALSE)]
      t(chol(A))[lower.tri(A, diag=TRUE)]
    }

    if(!nugget && cov.model != "matern"){
      
      K <- t(samples[,paste("K_",1:nltr,sep="")]); A <- t(apply(K, 2, A.chol, m))
      phi <- t(samples[,paste("phi_",1:m,sep="")])
  
      K.mean <- samples.mean[paste("K_",1:nltr,sep="")]; A.mean <- A.chol(K.mean, m) 
      phi.mean <- samples.mean[paste("phi_",1:m,sep="")]
      
    }else if(nugget && cov.model != "matern"){
      K <- t(samples[,paste("K_",1:nltr,sep="")]); A <- t(apply(K, 2, A.chol, m))
      Psi <- t(samples[,paste("Psi_",1:nltr,sep="")])
      phi <- t(samples[,paste("phi_",1:m,sep="")])

      K.mean <- samples.mean[paste("K_",1:nltr,sep="")]; A.mean <- A.chol(K.mean, m)
      Psi.mean <- samples.mean[paste("Psi_",1:nltr,sep="")]
      phi.mean <- samples.mean[paste("phi_",1:m,sep="")]
      
    }else if(!nugget && cov.model == "matern"){
      K <- t(samples[,paste("K_",1:nltr,sep="")]); A <- t(apply(K, 2, A.chol, m))
      phi <- t(samples[,paste("phi_",1:m,sep="")])
      nu <- t(samples[,paste("nu_",1:m,sep="")])

      K.mean <- samples.mean[paste("K_",1:nltr,sep="")]; A.mean <- A.chol(K.mean, m)
      phi.mean <- samples.mean[paste("phi_",1:m,sep="")]
      nu.mean <- samples.mean[paste("nu_",1:m,sep="")]

    }else{
      K <- t(samples[,paste("K_",1:nltr,sep="")]); A <- t(apply(K, 2, A.chol, m))
      Psi <- t(samples[,paste("Psi_",1:nltr,sep="")])
      phi <- t(samples[,paste("phi_",1:m,sep="")])
      nu <- t(samples[,paste("nu_",1:m,sep="")])
      
      K.mean <- samples.mean[paste("K_",1:nltr,sep="")]; A.mean <- A.chol(K.mean, m)
      Psi.mean <- samples.mean[paste("Psi_",1:nltr,sep="")]
      phi.mean <- samples.mean[paste("phi_",1:m,sep="")]
      nu.mean <- samples.mean[paste("nu_",1:m,sep="")]
    }


    lower2full <- function(x, m){
      A <- matrix(0, m, m)
      A[lower.tri(A, diag=TRUE)] <- x
      A[upper.tri(A, diag=FALSE)] <- t(A)[upper.tri(A, diag=FALSE)]
      A
    }
    
    if(!nugget){
      stop("Unmarginalized DIC cannot be calculated for a model with Psi=0")
    }

    status <- 0

    d <- rep(0, n.samples)
    DIC <- matrix(0,4,1)
    
    cat("-------------------------------------------------\n\t\tCalculating DIC\n-------------------------------------------------\n")
  
    for(s in 1:n.samples){
      
      KK <- lower2full(K[,s], m)
      
      PP <- matrix(0, m, m)
      if(nugget){
        PP <- lower2full(Psi[,s], m)
      }
      
      theta <- phi[,s]
      
      if(cov.model=="matern"){
        theta <- c(phi[,s], nu[,s])
      }
      
      Q <- Y-X%*%beta[,s]-w[,s]
      
      tmp.mm <- matrix(0, m, m)
      tmp.nm <- matrix(0, n*m, 1)
      
      storage.mode(tmp.mm) <- "double"
      storage.mode(tmp.nm) <- "double"
      storage.mode(Q) <- "double"
      storage.mode(KK) <- "double"
      storage.mode(PP) <- "double"
      storage.mode(theta) <- "double"
      storage.mode(knots.D) <- "double"
      storage.mode(obs.knots.D) <- "double"
      storage.mode(q) <- "integer"
      storage.mode(n) <- "integer"
      storage.mode(m) <- "integer"
      
      if(!is.pp){
        
        d[s] <- .Call("spMvDIC", n, m, Q, PP, tmp.mm, tmp.nm)
        
      }else{##predictive process
        
        if(modified.pp){            
          
          d[s] <- .Call("spMPPMvDIC", Q, knots.D, obs.knots.D, n, m, q, PP, KK, theta, cov.model);
          
        }else{
          
          d[s] <- .Call("spMvDIC", n, m, Q, PP, tmp.mm, tmp.nm)
          
        }
      }
      
      if(verbose){
        if(status == 100){
          cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
          status <- 0
        }
        status <- status+1
      }
      
    }
    
    d.bar <- mean(d)
    
    ##
    ##Get d.bar.omega
    ##
    
    KK <- lower2full(K.mean, m)
    PP <- matrix(0, m, m)
    
    if(nugget){
      PP <- lower2full(Psi.mean, m)
    }
    
    theta <- phi.mean
    
    if(cov.model=="matern"){
      theta <- c(phi.mean, nu.mean)
    }
    
    Q <- Y-X%*%beta.mean-w.mean
    
    storage.mode(KK) <- "double"
    storage.mode(PP) <- "double"
    storage.mode(theta) <- "double"
    storage.mode(Q) <- "double"
    
    if(!is.pp){
      
      d.bar.omega <- .Call("spMvDIC", n, m, Q, PP, tmp.mm, tmp.nm)
      
    }else{##predictive process
      if(modified.pp){
        
        d.bar.omega <- .Call("spMPPMvDIC", Q, knots.D, obs.knots.D, n, m, q, PP, KK, theta, cov.model);
        
      }else{
        
        d.bar.omega <- .Call("spMvDIC", n, m, Q, PP, tmp.mm, tmp.nm)
      }
    }
    
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
  
    out <- list()
    out$DIC <- DIC
    
    out  
            
  }else if(class(sp.obj) == "bayesLMRef"){
    
    cat("-------------------------------------------------\n\t\tCalculating DIC\n-------------------------------------------------\n")
          
    X <- sp.obj$X
    n <- nrow(X)
    p <- ncol(X)
    Y <- sp.obj$Y
    samples <- sp.obj$p.samples
    n.samples <- sp.obj$n.samples

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
    sigma.sq <- samples[,"sigma.sq"]    
    
    d <- rep(0, n.samples)
    
    ##
    ##DIC
    ##
    
    DIC <- matrix(0,4,1)
    
    for(s in 1:n.samples){
      d[s] <- n*log(sigma.sq[s])+1/(sigma.sq[s])*(t(Y-X%*%beta[s,])%*%(Y-X%*%beta[s,]))
    }
    cat(paste("Sampled: ",n.samples," of ",n.samples,", ",100,"%\n", sep=""))
    
    d.bar <- mean(d)
    
    sigma.sq.mu <- mean(sigma.sq)
    beta.mu <- as.matrix(colMeans(beta))
    
    d.bar.omega <- n*log(sigma.sq.mu)+1/(sigma.sq.mu)*(t(Y-X%*%beta.mu)%*%(Y-X%*%beta.mu))
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic

    DIC
    
  }else if(class(sp.obj) == "bayesLMConjugate"){

    cat("-------------------------------------------------\n\t\tCalculating DIC\n-------------------------------------------------\n")
          
    X <- sp.obj$X
    n <- nrow(X)
    p <- ncol(X)
    Y <- sp.obj$Y
    samples <- sp.obj$p.samples
    n.samples <- sp.obj$n.samples

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
    sigma.sq <- samples[,"sigma.sq"]    
    
    d <- rep(0, n.samples)
    
    ##
    ##DIC
    ##
    
    DIC <- matrix(0,4,1)
    
    for(s in 1:n.samples){
      d[s] <- n*log(sigma.sq[s])+1/(sigma.sq[s])*(t(Y-X%*%beta[s,])%*%(Y-X%*%beta[s,]))
    }
    cat(paste("Sampled: ",n.samples," of ",n.samples,", ",100,"%\n", sep=""))
    
    d.bar <- mean(d)
    
    sigma.sq.mu <- mean(sigma.sq)
    beta.mu <- as.matrix(colMeans(beta))
    
    d.bar.omega <- n*log(sigma.sq.mu)+1/(sigma.sq.mu)*(t(Y-X%*%beta.mu)%*%(Y-X%*%beta.mu))
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic

    DIC
    
  }else if(class(sp.obj) == "spLM"){

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
   
    if(!nugget){
      stop("Unmarginalized DIC cannot be calculated for a model with tau.sq=0")
    }

    ##Removing DIC.marg
    DIC.marg <- FALSE
    DIC.unmarg <- TRUE
    
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
    storage.mode(obs.D) <- "double"
    storage.mode(obs.knots.D) <- "double"
    storage.mode(knots.D) <- "double"
    storage.mode(n.samples) <- "integer"
    storage.mode(w) <- "double"
    storage.mode(w.str) <- "double"
    storage.mode(sp.effects) <- "integer"
    storage.mode(DIC.marg) <- "integer"
    storage.mode(DIC.unmarg) <- "integer"
    
    out <- .Call("splmDIC",X, Y, is.pp, modified.pp, n, m, p, nugget, beta, sigma.sq, tau.sq, phi, nu,
                 obs.D, obs.knots.D, knots.D, cov.model, n.samples, w, w.str, sp.effects,
                 DIC.marg, DIC.unmarg, verbose)

    out$DIC <- out$DIC.unmarg
    out$DIC.unmarg <- NULL
    out
    
    
  }else if(class(sp.obj) == "bayesGeostatExact"){
    
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
      stop("error: in spPredict, specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
    }

    R.eigen <- eigen(R)
    R.vals <- R.eigen$values
    R.vecs <- R.eigen$vectors
    R.vects.t <- t(R.vecs)

    d <- rep(0, n.samples)

    DIC <- matrix(0,4,1)
    status <- 0
    
    cat("-------------------------------------------------\n\tRecovering spatial effects\n-------------------------------------------------\n")
    
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
      
      if(verbose){
        if(status == 10){
          cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
          status <- 0
        }
        status <- status+1
      }
    }
    cat(paste("Sampled: ",n.samples," of ",n.samples,", ",100,"%\n", sep=""))
    
    cat("-------------------------------------------------\n\tCalculating DIC\n-------------------------------------------------\n")
    
    for(s in 1:n.samples){
      d[s] <- n*log(tau.sq[s])+1/(tau.sq[s])*(t(Y-X%*%beta[s,]-w[,s])%*%(Y-X%*%beta[s,]-w[,s]))
    }
    cat(paste("Sampled: ",n.samples," of ",n.samples,", ",100,"%\n", sep=""))
    
    d.bar <- mean(d)
    
    sigma.sq.mu <- mean(sigma.sq)
    tau.sq.mu <- mean(tau.sq)
    beta.mu <- as.matrix(colMeans(beta))
    w.mu <- as.matrix(rowMeans(w))
    
    d.bar.omega <- n*log(tau.sq.mu)+1/(tau.sq.mu)*(t(Y-X%*%beta.mu-w.mu)%*%(Y-X%*%beta.mu-w.mu))
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
      
    out <- list()
    out$DIC <- DIC
    out$sp.effects <- w

    out

  }else if(class(sp.obj) == "spGGT"){

    cat("Computing DIC and associated statistics ...")
      
    
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

    if(no.Psi){
      stop("Unmarginalized DIC cannot be calculated for a model with Psi=0")
    }
    
    ##Removing DIC.marg
    DIC.marg <- FALSE
    DIC.unmarg <- FALSE
    
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

    out$DIC <- out$DIC.unmarg
    out$DIC.unmarg <- NULL
    out
    
  }else{
    stop("error: wrong class of input object.\n")
  }
}
