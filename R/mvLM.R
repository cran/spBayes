mvLM <- function(formula, data = parent.frame(), starting, tuning, priors, n.samples, sub.samples, verbose=TRUE, n.report=100, ...){
  
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
  ##Starting values
  ####################################################
  beta.starting <- 0
  L.starting <- 0

  if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
  
  names(starting) <- tolower(names(starting))   

  if(!"beta" %in% names(starting)){
    beta.starting <- coefficients(lm(Y~X-1))
  }else{
    beta.starting <- starting[["beta"]]
  }
     
  if(!"l" %in% names(starting)){stop("error: starting value for L must be specified")}

  L.starting <- starting[["l"]]
  if(length(L.starting) != m*(m-1)/2+m){stop(paste("error: L must be of length ",m*(m-1)/2+m," in starting value list",sep=""))}           

  storage.mode(beta.starting) <- "double"
  storage.mode(L.starting) <- "double"
  
  ####################################################
  ##Priors
  ####################################################
  Psi.IW <- 0

  if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
  names(priors) <- tolower(names(priors))
  
  if(!"psi.iw" %in% names(priors)){stop("error: Psi.IW must be specified")}
  Psi.IW <- priors[["psi.iw"]]
  
  if(!is.list(Psi.IW) || length(Psi.IW) != 2){stop("error: Psi.IW must be a list of length 2")}
  if(length(Psi.IW[[1]]) != 1 ){stop("error: Psi.IW[[1]] must be of length 1 (i.e., the IW df hyperparameter)")}   
  if(length(Psi.IW[[2]]) != m^2 ){stop(paste("error: Psi.IW[[2]] must be a vector or matrix of length, ",m^2, ", (i.e., the IW scale matrix hyperparameter)",sep=""))}
  
  ####################################################
  ##Tuning values
  ####################################################
  L.tuning <- 0
  
  if(missing(tuning)){stop("error: tuning value vector for the parameters must be specified")}
  
  names(tuning) <- tolower(names(tuning))
  
  if(!"l" %in% names(tuning)){stop("error: L must be specified in tuning value list")}
  L.tuning <- as.vector(tuning[["l"]])
  if(length(L.tuning) != m*(m-1)/2+m){stop(paste("error: L must be of length ",m*(m-1)/2+m," in tuning value list",sep=""))}
      
  storage.mode(L.tuning) <- "double" 
  
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

  out <- .Call("mvLM", Y, X, p, n, m, Psi.IW, L.starting, beta.starting, L.tuning, n.samples, verbose, n.report)
  out$Y <- Y
  out$X <- X
  out$n <- n
  out$m <- m
  
  out$p <- p
  out$sub.samples <- sub.samples
  
  ##subsample 
  out$p.samples <- mcmc(t(out$p.samples[,seq(sub.samples[1], sub.samples[2], by=as.integer(sub.samples[3]))]))
  out$n.samples <- nrow(out$p.samples)##get adjusted n.samples
  
  col.names <- rep("null",ncol(out$p.samples))
  
  nltr <- m*(m-1)/2+m
  
  col.names[1:p] <- x.names
  col.names[(p+1):ncol(out$p.samples)] <- rep("Psi",nltr)
     
  colnames(out$p.samples) <- col.names

  AtA <- function(x, m){
    A <- matrix(0, m, m)
    A[lower.tri(A, diag=TRUE)] <- x
    (A%*%t(A))[lower.tri(A, diag=TRUE)]
  }

  Psi.names <- paste("Psi[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
  colnames(out$p.samples)[colnames(out$p.samples)%in%"Psi"] <- Psi.names
  out$p.samples[,Psi.names] <- t(apply(out$p.samples[,Psi.names], 1, AtA, m)) 

  class(out) <- "mvLM"
  out  
}

