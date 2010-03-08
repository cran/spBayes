spDiag <- function(sp.obj, start=1, end, thin=1, n.report=100, verbose=TRUE, ...){
  
  ####################################################
  ##Check for unused args, thin, etc.
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  if(missing(sp.obj)){stop("error: spDIC expects sp.obj\n")}
  if(!class(sp.obj) %in% c("spGGT","spLM","bayesGeostatExact","bayesLMConjugate", "bayesLMRef", "spMvLM", "spGLM", "mvLM")){
    stop("error: spDIC requires an output object of class bayesLMConjugate, bayesLMRef, spGGT, spLM, bayesGeostatExact, spMvLM, spGLM or mvLM\n")}
  if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}

  lower2full <- function(x, m){
    A <- matrix(0, m, m)
    A[lower.tri(A, diag=TRUE)] <- x
    A[upper.tri(A, diag=FALSE)] <- t(A)[upper.tri(A, diag=FALSE)]
    A
  }

  GP <- function(Y.rep, y){
    mu.rep <- apply(Y.rep, 1, mean)
    var.rep <- apply(Y.rep, 1, var)
    G <- sum((y-mu.rep)^2)
    P <- sum(var.rep)
    D <- G+P

    GPD <- matrix(0,3,1)
    rownames(GPD) <- c("G","P","D")
    colnames(GPD) <- "value"
    GPD[1,1] <- G
    GPD[2,1] <- P
    GPD[3,1] <- D
    GPD
  }
  
  ##thin samples
  samples <- sp.obj$p.samples
  n.samples <- nrow(samples)

  if(missing(end))
    end <- n.samples
  
  if(!is.numeric(start) || start >= n.samples)
    stop("error: invalid start")
  if(!is.numeric(end) || end > n.samples) 
    stop("error: invalid end")
  if(!is.numeric(thin) || thin >= n.samples) 
    stop("error: invalid thin")
  
  samples <- samples[seq(start, end, by=as.integer(thin)),,drop=FALSE]
  n.samples <- nrow(samples)

  if(n.samples <= 1){stop("error: too few MCMC samples")}
  
  if("sp.effects" %in% names(sp.obj)){
    w <- sp.obj$sp.effects[,seq(start, end, by=as.integer(thin))]
  }
  
  ####################################################
  ##Class specific
  ####################################################  
  if(class(sp.obj) == "mvLM"){
    Y <- sp.obj$Y
    X <- sp.obj$X
    n <- sp.obj$n
    m <- sp.obj$m##number or response variables
    p <- sp.obj$p
    nltr <- m*(m-1)/2+m

    ##extract samples
    beta <- as.matrix(samples[,1:p])
    
    Psi.names <- paste("Psi[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
    Psi <- as.matrix(samples[,Psi.names])
  
    ##get sample means
    beta.mu <- apply(beta, 2, mean)
    Psi.mu <- apply(Psi, 2, mean)
       
    status <- 0

    d <- rep(0, n.samples)
    DIC <- matrix(0,4,1)

    Y.rep <- matrix(0, n*m, n.samples)

    if(verbose)
      cat("-------------------------------------------------\n\t\tCalculating DIC & GP\n-------------------------------------------------\n")
  
    for(s in 1:n.samples){
      
      PP <- lower2full(Psi[s,], m)

      Q <- Y-X%*%beta[s,]

      d[s] <- n*determinant(PP)$modulus+t(Q)%*%(diag(n)%x%chol2inv(chol(PP)))%*%Q

      ##GP
      mu <- X%*%beta[s,]
      for(i in 1:n){
        Y.rep[((i-1)*m+1):((i-1)*m+m),s] <- mvrnorm(1, mu[((i-1)*m+1):((i-1)*m+m)], PP)
      }
      
      if(verbose){
        if(status == n.report){
          cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
          status <- 0
        }
        status <- status+1
      }
      
    }
    
    d.bar <- mean(d)
    
    ##Get d.bar.omega
    PP <- lower2full(Psi.mu, m)
        
    Q <- Y-X%*%beta.mu
    
    d.bar.omega <- n*determinant(PP)$modulus+t(Q)%*%(diag(n)%x%chol2inv(chol(PP)))%*%Q
      
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    colnames(DIC) <- c("value")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
  
    out <- list()
    out$DIC <- DIC
    out$GP <- GP(Y.rep,Y)
    
    out      
  }else if(class(sp.obj) == "spGLM"){
    
    X <- sp.obj$X
    n <- sp.obj$n
    p <- sp.obj$p
    Y <- sp.obj$Y
    family <- sp.obj$family

    if(family == "binomial")
      weights <- sp.obj$weights
    
    ##extract samples
    beta <- as.matrix(samples[,1:p])
    
    ##get sample means
    beta.mu <- apply(beta, 2, mean)
    w.mu <- as.matrix(apply(w,1,mean))
       
    status <- 0

    d <- rep(0, n.samples)
    DIC <- matrix(0,4,1)

    ##Y.rep <- matrix(0, n, n.samples)

    if(verbose)
      cat("-------------------------------------------------\n\t\tCalculating DIC & GP\n-------------------------------------------------\n")
     
    for(s in 1:n.samples){
      
      if(family == "poisson"){
        d[s] <- -2*sum(-exp(X%*%beta[s,]+w[,s])+Y*(X%*%beta[s,]+w[,s]))

        ##GP
        ##Y.rep[,s] <- rpois(n, exp(X%*%beta[s,]+w[,s]))

      }else if(family == "binomial"){
        pp <- 1/(1+exp(-X%*%beta[s,]-w[,s]))
        d[s] <- -2*sum(Y*log(pp)+(weights-Y)*log(1-pp))
        
        ##GP
        ##Y.rep[,s] <- rbinom(n, weights, prob=1/(1+exp(-(X%*%beta[s,]+w[,s]))))

      }else{
        stop("error: family is misspecified")
      }
      
    }
    
    d.bar <- mean(d)
        
    if(family == "poisson"){
      d.bar.omega <- -2*sum(-exp(X%*%beta.mu+w.mu)+Y*(X%*%beta.mu+w.mu))
    }else if(family == "binomial"){
      pp <- 1/(1+exp(-X%*%beta.mu-w.mu))
      d.bar.omega <- -2*sum(Y*log(pp)+(weights-Y)*log(1-pp))
    }else{
      stop("error: family is misspecified")
    }
    
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd

    DIC <- matrix(0,4,1)
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    colnames(DIC) <- c("value")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
    
    out <- list()    
    out$DIC <- DIC
    ##out$GP <- GP(Y.rep,Y)
    
    out

  }else if(class(sp.obj) == "spMvLM"){

    if(!sp.obj$nugget){
      stop("DIC cannot be computed for a no nugget model.")
    }
    
    Y <- sp.obj$Y
    X <- sp.obj$X
    n <- sp.obj$n
    m <- sp.obj$m##number or response variables
    q <- sp.obj$q##number of knots
    p <- sp.obj$p
    nltr <- m*(m-1)/2+m

    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp
    
    ##extract samples
    beta <- as.matrix(samples[,1:p])
    
    Psi.names <- paste("Psi[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
    Psi <- as.matrix(samples[,Psi.names])

    ##get sample means
    beta.mu <- apply(beta, 2, mean)
    Psi.mu <- apply(Psi, 2, mean)
    w.mu <- as.matrix(apply(w,1,mean))

    status <- 0
    
    d <- rep(0, n.samples)
    DIC <- matrix(0,4,1)

    Y.rep <- matrix(0, n*m, n.samples)

    if(verbose)
      cat("-------------------------------------------------\n\t\tCalculating DIC & GP\n-------------------------------------------------\n")
      
    ##if non-pp or non-modified pp
    if(is.pp && modified.pp){
      
      obs.coords <- sp.obj$coords
      knots.D <- sp.obj$knots.D
      obs.D <- sp.obj$coords.D
      obs.knots.D <- sp.obj$coords.knots.D
      cov.model <- sp.obj$cov.model
      
      K.names <- paste("K[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
      phi.names <- paste("phi_",1:m,sep="")
      nu.names <- paste("nu_",1:m,sep="")
      
      K <- as.matrix(samples[,K.names])
      phi <- as.matrix(samples[,phi.names])
      
      K.mu <- apply(K, 2, mean) 
      phi.mu <- apply(phi, 2, mean);
      
      if(cov.model=="matern"){
        nu <- as.matrix(samples[,nu.names])
        nu.mu <- apply(nu, 2, mean);
      }

      ##note, C_e_r is the mnxmn mxm block diag covariance matrix marginalized over \tile{\eps} use to compute GP
      C.eps.r <- matrix(0, m, n*m)
      
      for(s in 1:n.samples){
        
        KK <- lower2full(K[s,], m)
        PP <- lower2full(Psi[s,], m)
        
        theta <- phi[s,]
        
        if(cov.model=="matern"){
          theta <- c(phi[s,], nu[s,])
        }
        
        Q <- Y-X%*%beta[s,]-w[,s]
        
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
           
        d[s] <- .Call("spMPPMvDIC", Q, knots.D, obs.knots.D, n, m, q, PP, KK, theta, cov.model, C.eps.r);

        ##GP
        mu <- X%*%beta[s,]+w[,s]
        for(i in 1:n){
          Y.rep[((i-1)*m+1):((i-1)*m+m),s] <- mvrnorm(1, mu[((i-1)*m+1):((i-1)*m+m)], C.eps.r[,((i-1)*m+1):((i-1)*m+m)])
        }
        
        if(verbose){
          if(status == n.report){
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
      
      KK <- lower2full(K.mu, m)
      PP <- lower2full(Psi.mu, m)
       
      theta <- phi.mu
      
      if(cov.model=="matern"){
        theta <- c(phi.mu, nu.mu)
      }
      
      Q <- Y-X%*%beta.mu-w.mu
      
      storage.mode(KK) <- "double"
      storage.mode(PP) <- "double"
      storage.mode(theta) <- "double"
      storage.mode(Q) <- "double"

      ##note C.eps.r is not really needed here
      d.bar.omega <- .Call("spMPPMvDIC", Q, knots.D, obs.knots.D, n, m, q, PP, KK, theta, cov.model, C.eps.r)
      
      pd <- d.bar - d.bar.omega
      dic <- d.bar + pd
      
      rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
      colnames(DIC) <- c("value")
      DIC[1,1] <- d.bar
      DIC[2,1] <- d.bar.omega
      DIC[3,1] <- pd
      DIC[4,1] <- dic
      
      out <- list()
      out$DIC <- DIC
      out$GP <- GP(Y.rep,Y)
      
      out
      
    }else{

      for(s in 1:n.samples){
        
        PP <- lower2full(Psi[s,], m)
        
        Q <- Y-X%*%beta[s,]-w[,s]
        
        d[s] <- n*determinant(PP)$modulus+t(Q)%*%(diag(n)%x%chol2inv(chol(PP)))%*%Q

        ##GP
        mu <- X%*%beta[s,]+w[,s]
        for(i in 1:n){
          Y.rep[((i-1)*m+1):((i-1)*m+m),s] <- mvrnorm(1, mu[((i-1)*m+1):((i-1)*m+m)], PP)
        }
        
        if(verbose){
          if(status == n.report){
            cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
            status <- 0
          }
          status <- status+1
        }
        
      }
      
      d.bar <- mean(d)
      
      ##Get d.bar.omega
      PP <- lower2full(Psi.mu, m)
      
      Q <- Y-X%*%beta.mu-w.mu
      
      d.bar.omega <- n*determinant(PP)$modulus+t(Q)%*%(diag(n)%x%chol2inv(chol(PP)))%*%Q
      
      pd <- d.bar - d.bar.omega
      dic <- d.bar + pd
      
      rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
      colnames(DIC) <- c("value")
      DIC[1,1] <- d.bar
      DIC[2,1] <- d.bar.omega
      DIC[3,1] <- pd
      DIC[4,1] <- dic
      
      out <- list()
      out$DIC <- DIC
      out$GP <- GP(Y.rep,Y)
      
      out      
    }
            
  }else if(class(sp.obj) == "bayesLMRef"){

    if(verbose)
      cat("-------------------------------------------------\n\t\tCalculating DIC & GP\n-------------------------------------------------\n")
          
    X <- sp.obj$X
    n <- nrow(X)
    p <- ncol(X)
    Y <- sp.obj$Y
      
    ##get samples
    beta <- as.matrix(samples[,1:p])
    sigma.sq <- samples[,"sigma.sq"]    
    
    d <- rep(0, n.samples)

    Y.rep <- matrix(0, n, n.samples)

    status <- 0
    
    ##
    ##DIC
    ##
    
    DIC <- matrix(0,4,1)
    
    for(s in 1:n.samples){
      d[s] <- n*log(sigma.sq[s])+1/(sigma.sq[s])*(t(Y-X%*%beta[s,])%*%(Y-X%*%beta[s,]))

      ##GP
      Y.rep[,s] <- rnorm(n, X%*%beta[s,], sqrt(sigma.sq[s]))

      if(verbose){
        if(status == n.report){
          cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
          status <- 0
        }
        status <- status+1
      }
    }
    
    d.bar <- mean(d)
    
    sigma.sq.mu <- mean(sigma.sq)
    beta.mu <- as.matrix(colMeans(beta))
    
    d.bar.omega <- n*log(sigma.sq.mu)+1/(sigma.sq.mu)*(t(Y-X%*%beta.mu)%*%(Y-X%*%beta.mu))
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    colnames(DIC) <- c("value")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic

    out <- list()
    out$DIC <- DIC
    out$GP <- GP(Y.rep,Y)

    out
     
  }else if(class(sp.obj) == "bayesLMConjugate"){

    if(verbose)
      cat("-------------------------------------------------\n\t\tCalculating DIC & GP\n-------------------------------------------------\n")
          
    X <- sp.obj$X
    n <- nrow(X)
    p <- ncol(X)
    Y <- sp.obj$Y
    
    ##get samples
    beta <- as.matrix(samples[,1:p])
    sigma.sq <- samples[,"sigma.sq"]    
    
    d <- rep(0, n.samples)

    Y.rep <- matrix(0, n, n.samples)

    status <- 0
    
    ##
    ##DIC
    ##
    
    DIC <- matrix(0,4,1)
    
    for(s in 1:n.samples){
      d[s] <- n*log(sigma.sq[s])+1/(sigma.sq[s])*(t(Y-X%*%beta[s,])%*%(Y-X%*%beta[s,]))

      ##GP
      Y.rep[,s] <- rnorm(n, X%*%beta[s,], sqrt(sigma.sq[s]))
      
      if(verbose){
        if(status == n.report){
          cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
          status <- 0
        }
        status <- status+1
      }
    }
    
    d.bar <- mean(d)
    
    sigma.sq.mu <- mean(sigma.sq)
    beta.mu <- as.matrix(colMeans(beta))
    
    d.bar.omega <- n*log(sigma.sq.mu)+1/(sigma.sq.mu)*(t(Y-X%*%beta.mu)%*%(Y-X%*%beta.mu))
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    colnames(DIC) <- c("value")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic

    out <- list()
    out$DIC <- DIC
    out$GP <- GP(Y.rep,Y)
    out
    
  }else if(class(sp.obj) == "spLM"){

    if(!sp.obj$nugget){
      stop("DIC cannot be computed for a no nugget model.")
    }
    
    is.pp <- sp.obj$is.pp
    modified.pp <- sp.obj$modified.pp
    
    Y <- sp.obj$Y
    X <- sp.obj$X
    n <- sp.obj$n
    m <- sp.obj$m
    p <- sp.obj$p

    ##get samples
    beta <- as.matrix(samples[,1:p])
    tau.sq <- samples[,"tau.sq"]
 
    ##get sample means
    beta.mu <- apply(beta, 2, mean)
    tau.sq.mu <- mean(tau.sq)
    w.mu <- as.matrix(apply(w,1,mean))

    status <- 0
    
    d <- rep(0, n.samples)
    DIC <- matrix(0,4,1)

    Y.rep <- matrix(0, n, n.samples)

    if(verbose)
      cat("-------------------------------------------------\n\t\tCalculating DIC & GP\n-------------------------------------------------\n")
      
    ##if non-pp or non-modified pp
    if(is.pp && modified.pp){

      knots.D <- sp.obj$knots.D
      obs.knots.D <- sp.obj$coords.knots.D
      cov.model <- sp.obj$cov.model
      
      sigma.sq <- as.matrix(samples[,"sigma.sq"])
      phi <- as.matrix(samples[,"phi"])
      
      sigma.sq.mu <- mean(sigma.sq) 
      phi.mu <- mean(phi);
      
      if(cov.model=="matern"){
        nu <- as.matrix(samples[,"nu"])
        nu.mu <- mean(nu);
      }
      
      for(s in 1:n.samples){

        Q <- Y-X%*%beta[s,]-w[,s]
        
        if(cov.model != "matern"){
          ct <- sigma.sq[s]*spCor(obs.knots.D, cov.model, phi[s])
          C.str <- sigma.sq[s]*spCor(knots.D, cov.model, phi[s])
        }else{
          ct <- sigma.sq[s]*spCor(obs.knots.D, cov.model, phi[s], nu[s])
          C.str <- sigma.sq[s]*spCor(knots.D, cov.model, phi[s], nu[s])
        }

        ##ct C^{*-1} c
        C <- ct%*%chol2inv(chol(C.str))%*%t(ct)
          
        e <- tau.sq[s]+sigma.sq[s]-diag(C)
        
        d[s] <- sum(log(e))+t(Q)%*%(Q/e)

        ##GP
        mu <- X%*%beta[s,]+w[,s]
        Y.rep[,s] <- rnorm(n, mu, sqrt(e))
        
        if(verbose){
          if(status == n.report){
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
      
      Q <- Y-X%*%beta.mu-w.mu
      
      if(cov.model != "matern"){
        ct <- sigma.sq.mu*spCor(obs.knots.D, cov.model, phi.mu)
        C.str <- sigma.sq.mu*spCor(knots.D, cov.model, phi.mu)
      }else{
        ct <- sigma.sq.mu*spCor(obs.knots.D, cov.model, phi.mu, nu.mu)
        C.str <- sigma.sq.mu*spCor(knots.D, cov.model, phi.mu, nu.mu)
      }
      
      ##ct C^{*-1} c
      C <- ct%*%chol2inv(chol(C.str))%*%t(ct)
      
      e <- tau.sq.mu+sigma.sq.mu-diag(C)
      
      d.bar.omega <-  sum(log(e))+t(Q)%*%(Q/e)
      
      pd <- d.bar - d.bar.omega
      dic <- d.bar + pd
      
      rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
      colnames(DIC) <- c("value")
      DIC[1,1] <- d.bar
      DIC[2,1] <- d.bar.omega
      DIC[3,1] <- pd
      DIC[4,1] <- dic
      
      out <- list()
      out$DIC <- DIC
      out$GP <- GP(Y.rep,Y)
 
      out
      
    }else{
      
      for(s in 1:n.samples){
        
        Q <- Y-X%*%beta[s,]-w[,s]
        
        d[s] <- n*log(tau.sq[s])+(t(Q)%*%Q)/tau.sq[s]

        status <- 0
        
        ##GP
        mu <- X%*%beta[s,]+w[,s]
        Y.rep[,s] <- rnorm(n, mu, sqrt(tau.sq[s]))
                
        if(verbose){
          if(status == n.report){
            cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
            status <- 0
          }
          status <- status+1
        }
        
      }
      
      d.bar <- mean(d)
      
      ##Get d.bar.omega
   
      Q <- Y-X%*%beta.mu-w.mu
      
      d.bar.omega <- n*log(tau.sq.mu)+(t(Q)%*%Q)/tau.sq.mu
      
      pd <- d.bar - d.bar.omega
      dic <- d.bar + pd
      
      rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
      colnames(DIC) <- c("value")
      DIC[1,1] <- d.bar
      DIC[2,1] <- d.bar.omega
      DIC[3,1] <- pd
      DIC[4,1] <- dic
      
      out <- list()
      out$DIC <- DIC
      out$GP <- GP(Y.rep,Y)
      
      out
    }

  }else if(class(sp.obj) == "bayesGeostatExact"){
    
    X <- sp.obj$X
    n <- sp.obj$n
    p <- sp.obj$p
    Y <- sp.obj$Y
    coords <- sp.obj$coords
    cov.model <- sp.obj$cov.model
    phi <- sp.obj$phi
    alpha <- sp.obj$alpha
    
    if(cov.model == "matern")
      nu <- sp.obj$nu
    
    ##get samples
    beta <- as.matrix(samples[,1:p])
    tau.sq <- samples[,"tau.sq"]
    sigma.sq <- samples[,"sigma.sq"]    
    
    ##make R
    D <- as.matrix(dist(coords))
    
    if(cov.model != "matern"){
      R <- spCor(D, cov.model, phi)
    }else{
      R <- spCor(D, cov.model, phi, nu)
    }

    R.eigen <- eigen(R)
    R.vals <- R.eigen$values
    R.vecs <- R.eigen$vectors
    R.vects.t <- t(R.vecs)

    d <- rep(0, n.samples)

    Y.rep <- matrix(0, n, n.samples)
    
    DIC <- matrix(0,4,1)
    status <- 0

    if(verbose)
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
    
    for(s in 1:n.samples){
      
      ## Using rmvnorm is slow - it calculates the matrix square-root each time
      ## sp.effects[s,] <- rmvnorm(1, sp..posterior.mean[s,], sigma.sq.posterior.samples[s]*V.sp)
      
      ## Instead use the pre-computed V.sp.root
      z <- rnorm(nrow(V.sp), 0, 1)
      w[,s] <- sp.posterior.mean[s,] + sqrt(sigma.sq[s])*V.sp.root%*%z
      
      if(verbose){
        if(status == n.report){
          cat(paste("Sampled: ",s," of ",n.samples,", ",round(100*s/n.samples,2),"%\n", sep=""))
          status <- 0
        }
        status <- status+1
      }
    }

    if(verbose)
      cat("-------------------------------------------------\n\tCalculating DIC & GP\n-------------------------------------------------\n")

    for(s in 1:n.samples){
      d[s] <- n*log(tau.sq[s])+1/(tau.sq[s])*(t(Y-X%*%beta[s,]-w[,s])%*%(Y-X%*%beta[s,]-w[,s]))

      ##GP
      Y.rep[,s] <- rnorm(n, X%*%beta[s,]+w[,s], sqrt(tau.sq[s]))
    }
        
    d.bar <- mean(d)
    
    sigma.sq.mu <- mean(sigma.sq)
    tau.sq.mu <- mean(tau.sq)
    beta.mu <- as.matrix(colMeans(beta))
    w.mu <- as.matrix(rowMeans(w))
    
    d.bar.omega <- n*log(tau.sq.mu)+1/(tau.sq.mu)*(t(Y-X%*%beta.mu-w.mu)%*%(Y-X%*%beta.mu-w.mu))
    pd <- d.bar - d.bar.omega
    dic <- d.bar + pd
    
    rownames(DIC) <- c("bar.D", "D.bar.Omega", "pD", "DIC")
    colnames(DIC) <- c("value")
    DIC[1,1] <- d.bar
    DIC[2,1] <- d.bar.omega
    DIC[3,1] <- pd
    DIC[4,1] <- dic
      
    out <- list()
    out$DIC <- DIC
    out$sp.effects <- w
    out$GP <- GP(Y.rep,Y)
    
    out

  }else if(class(sp.obj) == "spGGT"){

    ##if spatial effect previously calculated
    if("sp.effects" %in% names(sp.obj)){
      sp.effect <- TRUE
    }else{
      sp.effect <- FALSE
    }
    
    sp.effect <- as.integer(sp.effect)
 
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

    if(!sp.effect){w <- NULL}
    storage.mode(w) <- "double"
    
    if(no.Psi){
      stop("DIC cannot be computed for a no nugget model.")
    }
    
    ##Removing DIC.marg
    DIC.marg <- FALSE
    DIC.unmarg <- TRUE
    
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
    out$DIC <- out$DIC.unmarg
    out$DIC.unmarg <- NULL
    out
    
  }else{
    stop("error: wrong class of input object.\n")
  }
}
