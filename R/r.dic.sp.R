"r.DIC.sp" <- function(sp.obj, start=2, end=numeric(0), thin=1, verbose=TRUE, ...){

  if(missing(sp.obj)){stop("error: r.DIC.sp expects an output object of class ggt.sp\n")}  
  if(class(sp.obj) != "ggt.sp"){stop("error: r.DIC.sp requires an output object of class ggt.sp\n")}
  if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}
  
  #subsample if specified
  samples <- as.matrix(sp.obj$samples)
  n.samples <- nrow(samples)

  if(missing(end))
    end <- n.samples
  
  if (!is.numeric(start) || start >= n.samples)
    stop("error: invalid start")
  if (!is.numeric(end) || end > n.samples) 
    stop("error: invalid end")
  if (!is.numeric(thin) || thin >= n.samples) 
    stop("error: invalid thin")

  samples <- samples[seq(start, end, by=as.integer(thin)),]
  n.samples <- nrow(samples)
  
  #get other stuff out of sp.obj
  cov.model <- sp.obj$cov.model
  no.tausq <- sp.obj$no.tausq
  X <- as.matrix(sp.obj$X)
  Y <- as.matrix(sp.obj$Y)
  coords <- as.matrix(sp.obj$coords)

  #get var component samples
  tausq <- as.numeric(0)
  if(!no.tausq)
    tausq <- as.matrix(samples[,"tausq"])
  
  sigmasq <- as.matrix(samples[,"sigmasq"])
  phi <- as.matrix(samples[,"phi"])

  nu <- as.numeric(0)
  if(cov.model == "matern")
    nu <- as.matrix(samples[,"nu"])
  
  #as this is a ggt.sp obj all theta are the first ncol(X) in the
  #sample matrix
  theta <- as.matrix(samples[,1:ncol(X)])

  cat("\n------------------------------------------------------------------------\n")
  cat("\tCalculating Deviance Information Criterion (DIC)")
  cat("\n------------------------------------------------------------------------\n")
    
  #get DIC
  D <- as.matrix(dist(coords))
  DD <- rep(0,n.samples)
  
  for(s in 1:n.samples){

    sigmasq.s <- sigmasq[s]
    if(!no.tausq)
      tausq.s <- tausq[s]
    phi.s <- phi[s]
    if(cov.model == "matern")
      nu.s <- nu[s]

    R <- exp(-phi.s*D)
    R <- sigmasq.s*(R)
    diag(R) <- diag(R) + tausq.s

    logDetR <- determinant(R, logarithm = TRUE)$modulus[1]
    
    theta.s <- as.matrix(theta[s,])
    DD[s] <- logDetR+t(Y-X%*%theta.s)%*%chol2inv(chol(R))%*%(Y-X%*%theta.s)
    
  }
  DDBar <- mean(DD)
  
  sigmasq.mean <- mean(sigmasq)
  if(!no.tausq)
    tausq.mean <- mean(tausq)
  phi.mean <- mean(phi)
  if(cov.model == "matern")
    nu.mean <- mean(nu)

  theta.mean <- as.matrix(colMeans(theta))
  
  R <- exp(-phi.mean*D)
  R <- sigmasq.mean*(R)
  diag(R) <- diag(R) + tausq.mean

  logDetR <- determinant(R, logarithm = TRUE)$modulus[1]
  
  DDBarOmega <- logDetR+t(Y-X%*%theta.mean)%*%chol2inv(chol(R))%*%(Y-X%*%theta.mean)

  pD <- DDBar - DDBarOmega
  DIC <- DDBar + DDBar - DDBarOmega

  return.DIC <- matrix(0, 4, 1)
  return.DIC[1,] <- DDBar
  return.DIC[2,] <- DDBarOmega
  return.DIC[3,] <- pD
  return.DIC[4,] <- DIC

  colnames(return.DIC) <- "value"
  rownames(return.DIC) <- c("D.bar","D.bar.Omega","pD","DIC")
  return.DIC
  
}
