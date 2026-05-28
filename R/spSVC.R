spSVC <- function(formula, data = parent.frame(), svc.cols=1, coords,  
                  priors, starting, tuning, cov.model, center.scale=FALSE,
                  amcmc, n.samples, n.omp.threads = 1, 
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
        
        holder <- parseFormula2(formula, data, na.action=NULL)
        Y <- holder[[1]]
        X <- as.matrix(holder[[2]])
        x.names <- holder[[3]]
        
    }else{
        stop("error: formula is misspecified")
    }

    ##check for na y or x, use drp on coords below
    dropped.obs <- !apply(cbind(Y,X), 1, function(x){any(is.na(x))})
    Y <- Y[dropped.obs]
    X <- X[dropped.obs,,drop=FALSE]
    if(any(!dropped.obs)){
        warning(paste0("Removing ",sum(!dropped.obs), " data rows due to NA values."))
    }

    p <- ncol(X)
    n <- nrow(X)

    X.sc <- list(apply(X, 2, mean), apply(X, 2, sd))
    
    if(center.scale){
        
        if(any(X.sc[[2]] == 0)){
            X.sc[[1]][X.sc[[2]] == 0] <- 0
            X.sc[[2]][X.sc[[2]] == 0] <- 1
        }
        
        X <- sapply(1:ncol(X), function(i){(X[,i]-X.sc[[1]][i])/X.sc[[2]][i]})
    }

    if(is.character(svc.cols)){
        if(!all(svc.cols %in% x.names)){
            missing.cols <- svc.cols[!(svc.cols %in% x.names)]
            stop(paste("error: column name ", paste(missing.cols, collapse=" and "), " not in design matrix columns", sep=""))
        }
        svc.cols <- (1:p)[x.names %in% svc.cols] ##convert desired column names into column index
        
    }else if(is.numeric(svc.cols)){
        if(!all(svc.cols %in% 1:p)){
            missing.cols <- svc.cols[!(svc.cols %in% (1:p))]
            stop(paste("error: column index ", paste(missing.cols, collapse=" "), " not in design matrix columns", sep=""))
        }
    }

    Z <- X[,svc.cols, drop=FALSE]
    m <- ncol(Z)
    n.ltr <- m*(m+1)/2
    
    storage.mode(Y) <- "double"
    storage.mode(X) <- "double"
    storage.mode(p) <- "integer"
    storage.mode(n) <- "integer"
    storage.mode(Z) <- "double"
    storage.mode(m) <- "integer"
    storage.mode(n.ltr) <- "integer"
    
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
    ##Distance matrices
    ####################################################
    if(missing(coords)){stop("error: coords must be specified")}
    
    if(is.vector(coords)){
        if(is.character(coords)){
            if(all(coords %in% colnames(data))){
                coords <- as.matrix(data[,coords])
            }else{
                stop(paste0("error: coords name ", paste(coords[!(coords %in% colnames(data))], collapse=" and "), " not in data"))
            }
        }else if(all(coords %in% (1:ncol(data)))){
            coords <- as.matrix(data[,coords])
        }else{
            stop(paste0("error: coords column index ", paste(coords[!(coords %in% (1:ncol(data)))], collapse=" and "), " not in data"))
        }
    }else{
        if(!any(is.matrix(coords), is.data.frame(coords))){
            stop("error: coords must n-by-m matrix or dataframe of coordinates or vector indicating the column names or integer indexes in data")
        }
        coords <- as.matrix(coords)
    }

    coords <- coords[dropped.obs,,drop=FALSE]
    
    if(nrow(coords) != n){
        stop("error: the number of rows in spatial coordinates is different than data used in the model formula")
    }
    
    if(any(duplicated(coords))){
        stop("error: duplicated coordinates found. Remove duplicates prior to running spSVC")
    }
    
    q <- ncol(coords)

    storage.mode(q) <- "integer"
    storage.mode(coords) <- "double"
    
    ####################################################
    ##Covariance model
    ####################################################
    if(missing(cov.model)){stop("error: cov.model must be specified")}
    if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}
    
    ####################################################
    ##Priors
    ####################################################
    if(missing(priors)) {stop("error: prior list for the parameters must be specified")}
    
    names(priors) <- tolower(names(priors))

    beta.Norm <- 0
    beta.prior <- "flat"
   
    if("beta.norm" %in% names(priors)){
        beta.Norm <- priors[["beta.norm"]]
        if(!is.list(beta.Norm) || length(beta.Norm) != 2){stop("error: beta.Norm must be a list of length 2")}
        if(length(beta.Norm[[1]]) != p ){stop(paste("error: beta.Norm[[1]] must be a vector of length, ",p, "",sep=""))}
        if(length(beta.Norm[[2]]) != p^2 ){stop(paste("error: beta.Norm[[2]] must be a ",p,"x",p," covariance matrix",sep=""))}
        beta.prior <- "normal"
    }

    K.prior <- 0   
    K.diag <- TRUE
    if("sigma.sq.ig" %in% names(priors)){
        
        K.prior <- priors[["sigma.sq.ig"]]
        if(!is.list(K.prior) || length(K.prior) != 2){stop("error: sigma.sq.IG must be a list of length 2")}
        if(length(K.prior[[1]]) != m){stop(paste("error: sigma.sq.IG[[1]] must be a vector of length, ",m, "",sep=""))}
        if(length(K.prior[[2]]) != m){stop(paste("error: sigma.sq.IG[[2]] must be a vector of length, ",m, "",sep=""))}
                 
    }else if("k.iw" %in% names(priors)){
        
        K.prior <- priors[["k.iw"]]
        if(!is.list(K.prior) || length(K.prior) != 2){stop("error: K.IW must be a list of length 2")}
        if(length(K.prior[[1]]) != 1){stop("error: K.IW[[1]] must be of length 1 (i.e., the IW df hyperparameter)")}   
        if(length(K.prior[[2]]) != m^2){stop(paste("error: K.IW[[2]] must be a vector or matrix of length, ",m^2, ", (i.e., the IW scale matrix hyperparameter)",sep=""))}
        K.diag <- FALSE
        
    }else{
        stop("error: sigma.sq.ig or k.iw must be specified")
    }

    if(!"phi.unif" %in% names(priors)){stop("error: phi.Unif must be specified")}
    phi.Unif <- priors[["phi.unif"]]
    if(!is.list(phi.Unif) || length(phi.Unif) != 2){stop("error: phi.Unif must be a list of length 2")}
    if(length(phi.Unif[[1]]) != m){stop(paste("error: phi.Unif[[1]] must be a vector of length, ",m, "",sep=""))}
    if(length(phi.Unif[[2]]) != m){stop(paste("error: phi.Unif[[2]] must be a vector of length, ",m, "",sep=""))}
    if(any(phi.Unif[[2]]-phi.Unif[[1]] <= 0)){stop("error: phi.Unif has zero support")}
    phi.Unif <- as.vector(t(cbind(phi.Unif[[1]],phi.Unif[[2]])))

    nu.Unif <- 0
 
    if(cov.model == "matern"){
        
        if(!"nu.unif" %in% names(priors)){stop("error: nu.Unif must be specified")}
        nu.Unif <- priors[["nu.unif"]]
        if(!is.list(nu.Unif) || length(nu.Unif) != 2){stop("error: nu.Unif must be a list of length 2")}
        if(length(nu.Unif[[1]]) != m){stop(paste("error: nu.Unif[[1]] must be a vector of length, ",m, "",sep=""))}
        if(length(nu.Unif[[2]]) != m){stop(paste("error: nu.Unif[[2]] must be a vector of length, ",m, "",sep=""))}
        if(any(nu.Unif[[2]]-nu.Unif[[1]] <= 0)){stop("error: nu.Unif has zero support")}
        nu.Unif <- as.vector(t(cbind(nu.Unif[[1]],nu.Unif[[2]])))
    }
    
    tau.sq.IG <- 0
    nugget <- FALSE
    
    if("tau.sq.ig" %in% names(priors)){
        tau.sq.IG <- priors[["tau.sq.ig"]]
        
        if(!is.vector(tau.sq.IG) || length(tau.sq.IG) != 2){stop("error: tau.sq.IG must be a vector of length 2")}
        if(any(tau.sq.IG <= 0)){stop("error: tau.sq.IG must be a positive vector of length 2")}
        nugget <- TRUE
    }
    

    storage.mode(K.prior[[1]]) <- "double"; storage.mode(K.prior[[2]]) <- "double"
    storage.mode(tau.sq.IG) <- "double"
    storage.mode(nu.Unif) <- "double"
    storage.mode(phi.Unif) <- "double"
    storage.mode(nugget) <- "integer"
    storage.mode(K.diag) <- "integer"
    ####################################################
    ##Starting values
    ####################################################
    if(missing(starting)){stop("error: starting value list for the parameters must be specified")}
    
    names(starting) <- tolower(names(starting))   
    
    beta.starting <- 0
    if(!"beta" %in% names(starting)){
        beta.starting <- as.vector(coefficients(lm(Y~X-1)))
    }else{
        beta.starting <- starting[["beta"]]
    }
    
    A.starting <- 0
    
    if(K.diag){
        
        if(!"sigma.sq" %in% names(starting)){stop("error: sigma.sq must be specified in starting")}
        
        A.starting <- as.vector(starting[["sigma.sq"]])
        if(length(A.starting) != m){stop(paste("error: sigma.sq must be of length ",m," in starting value list",sep=""))}
        
    }else{
      
        if(!"a" %in% names(starting)){stop("error: A must be specified in starting")}
        A.starting <- as.vector(starting[["a"]])
        if(length(A.starting) != n.ltr){stop(paste("error: A must be of length ",n.ltr," in starting value list",sep=""))}
        
    }
    
    if(!"phi" %in% names(starting)){stop("error: phi must be specified in starting")}
    phi.starting <- starting[["phi"]]
    if(length(phi.starting) != m){stop(paste("error: phi must be of length ",m," in starting value list",sep=""))}
           
    nu.starting <- 0
    if(cov.model == "matern"){
        if(!"nu" %in% names(starting)){stop("error: nu must be specified in starting")}
        nu.starting <- starting[["nu"]]
        if(length(nu.starting) != m){stop(paste("error: nu must be of length ",m," in starting value list",sep=""))}
    }


    if(any(phi.starting <= phi.Unif[seq(1,2*m,2)]) || any(phi.starting >= phi.Unif[seq(2,2*m,2)])){
        stop("one or more phi.starting is not in its prior support")
    }
    
    if(cov.model == "matern"){
        if(any(nu.starting <= nu.Unif[seq(1,2*m,2)]) || any(nu.starting >= nu.Unif[seq(2,2*m,2)])){
            stop("one or more nu.starting is not in its prior support")
        }
    }
 
        
    tau.sq.starting <- 0
    if(nugget){
        if(!"tau.sq" %in% names(starting)){stop("error: a prior was specified for tau.sq therefore tau.sq must be specified in starting value list")}
        tau.sq.starting <- starting[["tau.sq"]][1]
    }
    
    storage.mode(beta.starting) <- "double"
    storage.mode(phi.starting) <- "double"
    storage.mode(A.starting) <- "double"
    storage.mode(tau.sq.starting) <- "double"
    storage.mode(nu.starting) <- "double"
    
    ####################################################
    ##Tuning values
    ####################################################
    if(missing(tuning)){stop("error: tuning value vector for the spatial parameters must be specified")}
    
    names(tuning) <- tolower(names(tuning))
    
    A.tuning <- 0
    
    if(K.diag){
        
        if(!"sigma.sq" %in% names(tuning)){stop("error: sigma.sq must be specified in tuning")}
        
        A.tuning <- as.vector(tuning[["sigma.sq"]])
        if(length(A.tuning) != m){stop(paste("error: sigma.sq must be of length ",m," in tuning value list",sep=""))}
        
    }else{
        
        if(!"a" %in% names(tuning)){stop("error: A must be specified in tuning")}
        
        A.tuning <- as.vector(tuning[["a"]])
        if(length(A.tuning) != n.ltr){stop(paste("error: A must be of length ",n.ltr," in tuning value list",sep=""))}
        
    }
        
    if(!"phi" %in% names(tuning)){stop("error: phi must be specified in tuning value list")}
    phi.tuning <- tuning[["phi"]]
    if(length(phi.tuning) != m){stop(paste("error: phi must be of length ",m," in tuning value list",sep=""))}
    
    nu.tuning <- 0
    if(cov.model == "matern"){
        if(!"nu" %in% names(tuning)){stop("error: nu must be specified in tuning value list")}
        nu.tuning <- tuning[["nu"]]
        if(length(nu.tuning) != m){stop(paste("error: nu must be of length ",m," in tuning value list",sep=""))}
    }    
    
    if(nugget){
        if(!"tau.sq" %in% names(tuning)){stop("error: tau.sq must be specified in tuning value list")}
        tau.sq.tuning <- tuning[["tau.sq"]][1]
    }
    
    storage.mode(phi.tuning) <- "double"
    storage.mode(A.tuning) <- "double"
    storage.mode(tau.sq.tuning) <- "double"
    storage.mode(nu.tuning) <- "double"
  
    ####################################################
    ##Other stuff
    ####################################################
    storage.mode(n.report) <- "integer"
    storage.mode(verbose) <- "integer"
    storage.mode(n.omp.threads) <- "integer"
    
    ####################################################
    ##Pack it up and off it goes
    ####################################################
    ptm <- proc.time()

    out <- .Call("spSVC", Y, X, p, Z, m, n, coords, q,
                 K.diag, beta.prior, beta.Norm, K.prior, tau.sq.IG, nu.Unif, phi.Unif,
                 phi.starting, A.starting, tau.sq.starting, nu.starting,
                 phi.tuning, A.tuning, tau.sq.tuning, nu.tuning,
                 nugget, cov.model, is.amcmc, n.batch, batch.length, accept.rate, verbose, n.report, n.omp.threads)   
    
    
    run.time <- proc.time() - ptm
    
    out$p.theta.samples <- mcmc(t(out$p.theta.samples))

    sigma.sq.names <- paste0("sigma.sq.",x.names[svc.cols])
    phi.names <- paste0("phi.",x.names[svc.cols])
    nu.names <- paste0("nu.",x.names[svc.cols])
    
    if(!nugget && cov.model != "matern"){
        col.names <- c(sigma.sq.names, phi.names)
    }else if(nugget && cov.model != "matern"){
        col.names <- c(sigma.sq.names, "tau.sq", phi.names)
    }else if(!nugget && cov.model == "matern"){
        col.names <- c(sigma.sq.names, phi.names, nu.names)
    }else{
        col.names <- c(sigma.sq.names, "tau.sq", phi.names, nu.names)
    }

    if(!K.diag){
        col.names <- col.names[!(col.names %in% sigma.sq.names)]
        K.names <- paste("K[",matrix(apply(cbind(expand.grid(1:m,1:m)), 1, function(x) paste(x, collapse=",")),m,m)[lower.tri(matrix(0,m,m), diag=TRUE)],"]",sep="")
        col.names <- c(K.names, col.names)
    }

    colnames(out$p.theta.samples) <- col.names
    
    if(!K.diag){
        AtA <- function(x, m){
            A <- matrix(0, m, m)
            A[lower.tri(A, diag=TRUE)] <- x
            (A%*%t(A))[lower.tri(A, diag=TRUE)]
        }
        
        out$p.theta.samples[,K.names] <- t(apply(out$p.theta.samples[,K.names,drop=FALSE], 1, AtA, m))   
    }
         
    out$Y <- Y
    out$X <- X
    out$Z <- Z
    out$center.scale <- center.scale
    out$X.sc <- X.sc
    out$coords <- coords
    out$cov.model <- cov.model
    out$nugget <- nugget
    out$beta.prior <- beta.prior
    out$beta.Norm <- beta.Norm
    out$x.names <- x.names
    out$run.time <- run.time
    out$K.diag <- K.diag
    out$svc.cols <- svc.cols
    out$dropped.obs <- !dropped.obs
    class(out) <- "spSVC"
    out  
}

