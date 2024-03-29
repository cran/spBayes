spPredict <- function(sp.obj, pred.coords, pred.covars, joint=FALSE, start=1, end, thin=1, verbose=TRUE, n.report=100, n.omp.threads=1, ...){
  
    ####################################################
    ##Check for unused args
    ####################################################
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for(i in elip.args){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }
    
    if(missing(sp.obj)){stop("error: spPredict expects sp.obj\n")}
    ##if(!class(sp.obj) %in% c("spLM", "spMvLM", "spGLM", "spMvGLM","bayesLMRef","bayesLMConjugate","bayesGeostatExact","nonSpGLM","nonSpMvGLM","spMisalignLM","spMisalignGLM", "spSVC")){
    if(!inherits(sp.obj, c("spLM", "spMvLM", "spGLM", "spMvGLM","bayesLMRef","bayesLMConjugate","bayesGeostatExact","nonSpGLM","nonSpMvGLM","spMisalignLM","spMisalignGLM", "spSVC"))){    
        stop("error: requires an output object of class spLM, spMvLM, spGLM, spMvGLM, bayesGeostatExact, bayesLMRef, bayesLMConjugate, nonSpGLM, nonSpMvGLM, spMisalignLM, spMisalignGLM, or spSVC \n")
    }
    
    ##obj.class <- class(sp.obj)
    
    ##
    ##spSVC
    ##
    ##if(obj.class == "spSVC"){
    if(inherits(sp.obj, "spSVC")){        

        if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
        if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
                
        ##check that spRecover was previously run
        if(all(c("p.beta.recover.samples", "p.w.recover.samples", "p.theta.recover.samples") %in% names(sp.obj))){
            
            n.samples <- nrow(sp.obj$p.theta.recover.samples)
            
            if(missing(end)){end <- n.samples}
            if(!is.numeric(start) || start >= n.samples){stop("error: invalid start")}
            if(!is.numeric(end) || end > n.samples){stop("error: invalid end")}
            if(!is.numeric(thin) || thin >= n.samples){stop("error: invalid thin")}
            
            s.indx <- seq(as.integer(start), as.integer(end), by=as.integer(thin))
            
            p.beta.samples <- sp.obj$p.beta.recover.samples[s.indx,,drop=FALSE]       
            p.theta.samples <- sp.obj$p.theta.recover.samples[s.indx,,drop=FALSE]
            p.w.samples <- sp.obj$p.w.recover.samples[,s.indx,drop=FALSE]
            
            n.samples <- nrow(p.theta.samples)

            if(verbose){
                message(paste0("Using ", n.samples, " posterior samples from previous spRecover call."))
            }
            
        }else{
            stop("error: run sp.obj <- spRecover(sp.obj, get.beta=T, get.w=T, ...) before calling spDiag")
        }
        
        Y <-sp.obj$Y
        X <- sp.obj$X
        Z <- sp.obj$Z
        p <- ncol(X)
        n <- nrow(X)
        m <- ncol(Z)
        n.ltr <- m*(m+1)/2
        K.diag <- sp.obj$K.diag
        coords <- sp.obj$coords
        n.params <- ncol(p.theta.samples)
        cov.model <- sp.obj$cov.model
        nugget <- sp.obj$nugget
        beta.prior <- sp.obj$beta.prior
        beta.Norm <- sp.obj$beta.Norm
        x.names <- sp.obj$x.names
        svc.cols <- sp.obj$svc.cols
        center.scale <- sp.obj$center.scale
        X.sc <- sp.obj$X.sc
        
        if(ncol(coords) != ncol(pred.coords)){
            stop(paste0("Number of columns in coords does not equal that of pred.coords."))
        }
        
        if(nrow(pred.covars) != nrow(pred.coords)){
            stop(paste0("Number of rows in pred.covars does not equal that of pred.coords."))
        }

        if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
        pred.covars <- as.matrix(pred.covars)
  
        if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
        
        if(center.scale){
            pred.covars <- sapply(1:ncol(X), function(i){(pred.covars[,i]-X.sc[[1]][i])/X.sc[[2]][i]})
            message("spSVC argument center.scale=TRUE, centering and scaling non intercept columns of pred.covars respectively.")
        }
        
        q <- nrow(pred.coords)
        
        obs.D <- iDist(coords)
        pred.obs.D <- iDist(pred.coords, coords)
        if(joint){
            pred.D <- iDist(pred.coords)
        }
        
        ##parameters
        A.indx <- 0; tau.sq.indx <- 0; phi.indx <- 0; nu.indx <- 0
        
        if(K.diag){
            j <- m
        }else{
            j <- n.ltr
        }
        
        if(!nugget && cov.model != "matern"){
            A.indx <- 0; phi.indx <- A.indx+j
        }else if(nugget && cov.model != "matern"){
            A.indx <- 0; tau.sq.indx <- A.indx+j; phi.indx <- tau.sq.indx+1
        }else if(!nugget && cov.model == "matern"){
            A.indx <- 0; phi.indx <- A.indx+j; nu.indx <- phi.indx+m
        }else{
            A.indx <- 0; tau.sq.indx <- A.indx+j; phi.indx <- tau.sq.indx+1; nu.indx <- phi.indx+m
        }
        
        storage.mode(m) <- "integer"
        storage.mode(n) <- "integer"
        storage.mode(K.diag) <- "integer"
        storage.mode(obs.D) <- "double"
        storage.mode(pred.obs.D) <- "double"
        if(joint){
            storage.mode(pred.D) <- "double"
        }
        storage.mode(q) <- "integer" 
        storage.mode(p.theta.samples) <- "double"
        storage.mode(p.w.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(A.indx) <- "integer"
        storage.mode(phi.indx) <- "integer"
        storage.mode(nu.indx) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(n.omp.threads) <- "integer"

        ptm <- proc.time()
        
        if(joint){
            out <- .Call("spSVCPredictJoint", m, n, K.diag, obs.D, pred.obs.D, pred.D, q,
                         p.theta.samples, p.w.samples, n.samples, 
                         A.indx, phi.indx, nu.indx,
                         cov.model, 
                         verbose, n.report, n.omp.threads)
        }else{
            out <- .Call("spSVCPredictMarg", m, n, K.diag, obs.D, pred.obs.D, q,
                         p.theta.samples, p.w.samples, n.samples, 
                         A.indx, phi.indx, nu.indx,
                         cov.model, 
                         verbose, n.report, n.omp.threads)
        }

     
        run.time <- proc.time() - ptm
        out$run.time <- run.time

        out$center.scale.pred.covars <- pred.covars
        
        ##save the samples used for prediction
        colnames(p.beta.samples) <- x.names
        out$p.beta.recover.samples <- mcmc(p.beta.samples)
        out$p.theta.recover.samples <- mcmc(p.theta.samples)
        out$recover.args <- c(start, end, thin)
        
        ##also save as list with elements corresponding to svc.cols
        out$p.w.predictive.samples.list  <- vector("list", length = m)
        names(out$p.w.predictive.samples.list) <- paste0("w.",x.names[svc.cols])
        indx <- rep(1:m, q)
        for(i in 1:m){
            out$p.w.predictive.samples.list[[i]] <- out$p.w.predictive.samples[indx == i,]
        }
       
        ##B.tilde = B + w for each svc.cols
        out$p.tilde.beta.predictive.samples.list  <- vector("list", length = m)
        names(out$p.tilde.beta.predictive.samples.list) <- paste0("tilde.beta.",x.names[svc.cols])
        for(i in 1:m){
            out$p.tilde.beta.predictive.samples.list[[i]] <- out$p.w.predictive.samples.list[[i]] + matrix(rep(as.matrix(out$p.beta.recover.samples)[,svc.cols[i]], each=q), nrow=q)
        }
       
        ##y predicted
        X.tilde <- t(bdiag(as.list(as.data.frame(t(pred.covars[,svc.cols,drop=FALSE])))))        
        out$p.y.predictive.samples <- matrix(0, q, n.samples)
         
        for(i in 1:n.samples){
            out$p.y.predictive.samples[,i] <- rnorm(q, (pred.covars%*%out$p.beta.recover.samples[i,] + X.tilde%*%out$p.w.predictive.samples[,i])[,1], sqrt(out$p.theta.recover.samples[i,"tau.sq"]))
        }
                
        return(out)
    }
    
    
    ##
    ##non spatial model prediction
    ##
    ##if(obj.class %in% c("nonSpGLM", "nonSpMvGLM")){
    if(inherits(sp.obj, c("nonSpGLM", "nonSpMvGLM"))){
            
        X <- sp.obj$X
        Y <- sp.obj$Y
        family <- sp.obj$family
        p.beta.samples <- sp.obj$p.beta.samples
        n.samples <- nrow(p.beta.samples)
        
        ##subsamples
        if(missing(end)){end <- n.samples}
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        s.indx <- seq(start, end, by=as.integer(thin))
        
        p.beta.samples <- p.beta.samples[s.indx,,drop=FALSE]
        n.samples <- nrow(p.beta.samples)
        
        out <- list()
        
        if(family == "binomial"){
            out$p.y.predictive.samples <- apply(p.beta.samples, 1,  function(s){1/(1+exp(-pred.covars%*%s))})          
        }else{##poisson
            out$p.y.predictive.samples <- apply(p.beta.samples, 1,  function(s){exp(pred.covars%*%s)}) 
        }
        
        return(out)
    }
    
    ##
    ##bayesLMRef
    ##
    ##if(obj.class %in% c("bayesLMRef","bayesLMConjugate")){
    if(inherits(sp.obj, c("bayesLMRef","bayesLMConjugate"))){
        
        X <- sp.obj$X
        Y <- sp.obj$Y
        p <- ncol(X)
        p.beta.tauSq.samples <- sp.obj$p.beta.tauSq.samples
        n.samples <- nrow(p.beta.tauSq.samples)
        
        ##subsamples
        if(missing(end)){end <- n.samples}
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        s.indx <- seq(start, end, by=as.integer(thin))
        
        p.beta.tauSq.samples <- p.beta.tauSq.samples[s.indx,]
        n.samples <- nrow(p.beta.tauSq.samples)
        
        ##check covars
        if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
        if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
        if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
        
        out <- list()
        out$p.y.predictive.samples <- apply(p.beta.tauSq.samples, 1, function(s){rnorm(nrow(pred.covars), pred.covars%*%s[1:p], sqrt(s[p+1]))})
        
        return(out)
    }
    
    ##
    ##bayesGeostatExact
    ##
    ##if(obj.class == "bayesGeostatExact"){
    if(inherits(sp.obj, "bayesGeostatExact")){
        
        X <- sp.obj$X
        n <- sp.obj$n
        p <- sp.obj$p
        Y <- sp.obj$Y
        coords <- sp.obj$coords
        cov.model <- sp.obj$cov.model
        samples <- sp.obj$p.samples
        phi <- sp.obj$phi
        n.samples <- sp.obj$n.samples
        
        if(cov.model == "matern")
            nu <- sp.obj$nu
        
        ##check covars
        if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
        
        if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
        
        if(ncol(pred.covars) != ncol(X))
            stop(paste("error: pred.covars must have ",p," columns\n"))
        
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
        
        n.pred <- nrow(pred.coords)
        
        y.pred <- matrix(0, n.pred, n.samples)
        
        R.eigen <- eigen(R)
        R.vals <- R.eigen$values
        R.vecs <- R.eigen$vectors
        R.vects.t <- t(R.vecs)
        
        if(verbose){
            cat("-------------------------------------------------\n")
            cat("\t\tPredicting\n")
            cat("-------------------------------------------------\n")
            bar <- txtProgressBar(min=1, max=n.samples, initial=0, style=3, char="*", width=40)
        }
        
        ##for each pred point by each sample
        for(s in 1:n.samples){
            
            R.inv <- R.vecs%*%diag(1/(R.vals+tau.sq[s]/sigma.sq[s]))%*%t(R.vecs)
            
            for(i in 1:n.pred){
                
                D.pred <- sqrt((pred.coords[i,1]-coords[,1])^2 + (pred.coords[i,2]-coords[,2])^2)
                
                if(cov.model == "exponential"){
                    gamma <- exp(-phi*D.pred)
                }else if(cov.model == "matern"){
                    gamma <- (D.pred*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D.pred*phi, nu=nu)
                }else if(cov.model == "gaussian"){
                    gamma <- exp(-1*((phi*D.pred)^2))
                }else if(cov.model == "spherical"){
                    gamma <- D.pred
                    gamma[TRUE] <- 1
                    gamma[D.pred > 0 & D.pred < 1/phi] <- 1-1.5*phi*D.pred[D.pred > 0 & D.pred <= 1/phi]+0.5*((phi*D.pred[D.pred > 0 & D.pred <= 1/phi])^3)
                    gamma[D.pred >= 1/phi] <- 0   
                }else{
                    stop("error: in spPredict, specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
                }
                
                gamma <- as.matrix(gamma)
                
                mu <- pred.covars[i,]%*%beta[s,]+t(gamma)%*%R.inv%*%(Y-X%*%beta[s,])
                S <- sigma.sq[s]*(1-t(gamma)%*%R.inv%*%gamma)+tau.sq[s]
                
                y.pred[i,s] <- rnorm(1, mu, sqrt(S))
                
            }
            
            if(verbose){
                setTxtProgressBar(bar, s)
            }
        }
        
        out <- list()
        out$p.y.predictive.samples <- y.pred
        
        return(out)
    }
    
    
    ##
    ##spatial model prediction
    ##
    ##if(obj.class %in% c("spGLM", "spMvGLM")){
    if(inherits(sp.obj, c("spGLM", "spMvGLM"))){
        
        if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
        if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
        if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}
        
        family <- sp.obj$family
        X <- sp.obj$X
        Y <- sp.obj$Y
        p <- ncol(X)
        m <- 1 ##for spGLM
        ##if(obj.class == "spMvGLM"){m <- sp.obj$m}
        if(inherits(sp.obj, "spMvGLM")){m <- sp.obj$m}
        obs.coords <- sp.obj$coords
        n <- nrow(obs.coords)
        cov.model <- sp.obj$cov.model
        p.beta.theta.samples <- sp.obj$p.beta.theta.samples
        n.samples <- nrow(p.beta.theta.samples)
        is.pp <- sp.obj$is.pp
        q <- nrow(pred.coords)
        
        ##check covars
        if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
        if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
        if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
        
        ##subsamples
        if(missing(end)){end <- n.samples}
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        s.indx <- seq(start, end, by=as.integer(thin))
        
        p.beta.theta.samples <- t(p.beta.theta.samples[s.indx,,drop=FALSE])
        n.samples <- ncol(p.beta.theta.samples)
        
        knots.obs.D  <- NULL
        knots.D <- NULL
        knots.pred.D <- NULL
        obs.pred.D <- NULL
        obs.D <- NULL
        pred.D <- NULL
        p.w.samples <- NULL
        
        if(is.pp){
            knot.coords <- sp.obj$knot.coords
            g <- nrow(knot.coords)
            p.w.samples <- sp.obj$p.w.knots.samples[,s.indx,drop=FALSE]
            knots.D <- iDist(knot.coords)
            pred.knots.D <- iDist(pred.coords, knot.coords)
        }else{
            p.w.samples <- sp.obj$p.w.samples[,s.indx,drop=FALSE]
            obs.pred.D <- iDist(obs.coords, pred.coords)
            obs.D <- iDist(obs.coords)
        }
        
        storage.mode(X) <- "double"
        storage.mode(Y) <- "double"
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"
        storage.mode(m) <- "integer"
        storage.mode(pred.covars) <- "double"
        storage.mode(q) <- "integer"
        storage.mode(p.beta.theta.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(p.w.samples) <- "double"   
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        
        if(is.pp){
            
            storage.mode(g) <- "integer"
            storage.mode(knots.D) <- "double"
            storage.mode(pred.knots.D) <- "double"
            
            out <- .Call("spPPMvGLMPredict", family, X, Y, n, m, g, p, pred.covars, q, knots.D, pred.knots.D, 
                         p.beta.theta.samples, p.w.samples, n.samples,
                         cov.model, verbose, n.report)
        }else{
            
            storage.mode(obs.pred.D) <- "double"
            storage.mode(obs.D) <- "double"
            
            out <- .Call("spMvGLMPredict", family, X, Y, n, m, p, pred.covars, q, obs.D, obs.pred.D, 
                         p.beta.theta.samples, p.w.samples, n.samples, cov.model,
                         verbose, n.report)
        }
        
        out
        
    ##}else if(obj.class == "spLM"){
    }else if(inherits(sp.obj, "spLM")){
        
        if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
        if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
        if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}
        
        X <- sp.obj$X
        Y <- sp.obj$Y
        p <- ncol(X)
        n <- nrow(X)
        obs.coords <- sp.obj$coords
        cov.model <- sp.obj$cov.model
        p.theta.samples <- sp.obj$p.theta.samples
        n.samples <- nrow(p.theta.samples)
        is.pp <- sp.obj$is.pp
        nugget <- sp.obj$nugget
        beta.prior <- sp.obj$beta.prior
        beta.Norm <- sp.obj$beta.Norm
        ##r.indx <- sp.obj$r.indx
        get.beta <- sp.obj$get.beta
        
        ##check covars
        if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
        if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
        if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
        
        ##subsamples
        if(missing(end)){end <- n.samples}
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        s.indx <- seq(start, end, by=as.integer(thin))
        
        beta <- NULL;
        
        if(is.pp){
            beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
            knot.coords <- sp.obj$knot.coords
            m <- nrow(knot.coords)
            modified.pp <- sp.obj$modified.pp
        }   
        
        p.theta.samples <- p.theta.samples[s.indx,,drop=FALSE]##note, for flat we could use the p.theta.recover.samples
        
        n.samples <- nrow(p.theta.samples)
        
        ##recover beta if needed
        if(!is.pp && beta.prior == "flat"){
            ## if(all(s.indx %in% r.indx)){
            ##   beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
            ## }else{
            beta <- spRecover(sp.obj, get.beta=TRUE, get.w=FALSE, start=start, end=end, thin=thin)$p.beta.recover.samples
            ## }
        }
        
        sigma.sq.indx <- 0; tau.sq.indx <- 0; phi.indx <- 0; nu.indx <- 0
        
        if(!nugget && cov.model != "matern"){
            sigma.sq.indx <- 0; phi.indx <- 1
        }else if(nugget && cov.model != "matern"){
            sigma.sq.indx <- 0; tau.sq.indx <- 1; phi.indx <- 2
        }else if(!nugget && cov.model == "matern"){
            sigma.sq.indx <- 0; phi.indx <- 1; nu.indx <- 2
        }else{
            sigma.sq.indx <- 0; tau.sq.indx <- 1; phi.indx <- 2; nu.indx <- 3
        }
        
        knots.obs.D  <- NULL
        knots.D <- NULL
        knots.pred.D <- NULL
        obs.pred.D <- NULL
        obs.D <- NULL
        
        if(is.pp){
            knots.obs.D <- iDist(knot.coords, obs.coords)
            knots.D <- iDist(knot.coords)
            knots.pred.D <- iDist(knot.coords, pred.coords)
        }else{
            obs.pred.D <- iDist(obs.coords, pred.coords)
            obs.D <- iDist(obs.coords)
        }
        
        n.pred <- nrow(pred.coords)
        
        storage.mode(X) <- "double"
        storage.mode(Y) <- "double"
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"
        storage.mode(pred.covars) <- "double"
        storage.mode(n.pred) <- "integer"
        storage.mode(p.theta.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(beta) <- "double"
        storage.mode(sigma.sq.indx) <- "integer"
        storage.mode(tau.sq.indx) <- "integer"
        storage.mode(phi.indx) <- "integer"
        storage.mode(nu.indx) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        
        if(is.pp){
            
            storage.mode(m) <- "integer"
            storage.mode(knots.D) <- "double"
            storage.mode(knots.obs.D) <- "double"
            storage.mode(knots.pred.D) <- "double"
            storage.mode(modified.pp) <- "integer"
            
            out <- .Call("spPPLMPredict", X, Y, n, p, m, pred.covars, n.pred,
                         p.theta.samples, n.samples,
                         beta, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                         nugget, knots.D, knots.obs.D, knots.pred.D, cov.model, modified.pp, verbose, n.report)
        }else{
            
            storage.mode(obs.pred.D) <- "double"
            storage.mode(obs.D) <- "double"
            
            out <- .Call("spLMPredict", X, Y, n, p, pred.covars, n.pred,
                         p.theta.samples, n.samples,
                         beta.prior, beta.Norm, beta, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                         obs.D, obs.pred.D, cov.model, nugget, verbose, n.report)
        }
        
        out
        
    ##}else if(obj.class == "spMvLM"){
    }else if(inherits(sp.obj, "spMvLM")){
        
        if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
        if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
        if(!ncol(pred.coords) == 2){stop("error: pred.coords must have two columns (assumed to be X, Y)\n")}
        
        X <- sp.obj$X
        Y <- sp.obj$Y
        p <- ncol(X)
        m <- sp.obj$m
        obs.coords <- sp.obj$coords
        n <- nrow(obs.coords)
        cov.model <- sp.obj$cov.model
        p.theta.samples <- sp.obj$p.theta.samples
        n.samples <- nrow(p.theta.samples)
        is.pp <- sp.obj$is.pp
        nugget <- sp.obj$nugget
        Psi.diag <- sp.obj$Psi.diag
        beta.prior <- sp.obj$beta.prior
        beta.Norm <- sp.obj$beta.Norm
        ##r.indx <- sp.obj$r.indx
        get.beta <- sp.obj$get.beta
        
        q <- nrow(pred.coords)
        
        ##check covars
        if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
        if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
        if(ncol(pred.covars) != ncol(X)){ stop(paste("error: pred.covars must have ",p," columns\n"))}
        
        ##subsamples
        if(missing(end)){end <- n.samples}
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        s.indx <- seq(start, end, by=as.integer(thin))
        
        beta <- NULL
        
        if(is.pp){
            beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
            knot.coords <- sp.obj$knot.coords
            g <- nrow(knot.coords)
            modified.pp <- sp.obj$modified.pp
        }  
        
        p.theta.samples <- t(p.theta.samples[s.indx,,drop=FALSE])##note, for flat we could use the p.theta.recover.samples
        n.samples <- ncol(p.theta.samples)
        
        ##recover beta if needed (note, beta samples not needed for beta normal)
        if(!is.pp && beta.prior == "flat"){
            ## if(all(s.indx %in% r.indx)){
            ##   beta <- sp.obj$p.beta.samples[s.indx,,drop=FALSE]
            ## }else{
            beta <- spRecover(sp.obj, get.beta=TRUE, get.w=FALSE, start=start, end=end, thin=thin)$p.beta.recover.samples
            ## }
        }
        
        knots.obs.D  <- NULL
        knots.D <- NULL
        knots.pred.D <- NULL
        obs.pred.D <- NULL
        obs.D <- NULL
        
        if(is.pp){
            knots.obs.D <- iDist(knot.coords, obs.coords)
            knots.D <- iDist(knot.coords)
            knots.pred.D <- iDist(knot.coords, pred.coords)
        }else{
            obs.pred.D <- iDist(obs.coords, pred.coords)
            obs.D <- iDist(obs.coords)
        }
        
        storage.mode(X) <- "double"
        storage.mode(Y) <- "double"
        storage.mode(n) <- "integer"
        storage.mode(p) <- "integer"
        storage.mode(m) <- "integer"
        storage.mode(pred.covars) <- "double"
        storage.mode(q) <- "integer"
        storage.mode(p.theta.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(beta) <- "double"
        storage.mode(nugget) <- "integer"
        storage.mode(Psi.diag) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        
        if(is.pp){
            
            storage.mode(g) <- "integer"
            storage.mode(knots.D) <- "double"
            storage.mode(knots.obs.D) <- "double"
            storage.mode(knots.pred.D) <- "double"
            storage.mode(modified.pp) <- "integer"
            
            out <- .Call("spPPMvLMPredict", X, Y, n, m, g, p, pred.covars, q,
                         knots.D, knots.obs.D, knots.pred.D, 
                         p.theta.samples, beta, n.samples,
                         nugget, Psi.diag, cov.model,
                         modified.pp, verbose, n.report)
        }else{
            
            storage.mode(obs.pred.D) <- "double"
            storage.mode(obs.D) <- "double"
            
            
            Z <- t(pred.covars)
            storage.mode(Z) <- "double"
            
            out <- .Call("spMvLMPredict", X, Y, n, m, p, Z, q, obs.D, obs.pred.D, 
                         p.theta.samples, beta, n.samples,
                         beta.prior, beta.Norm,
                         nugget, Psi.diag, cov.model,
                         verbose, n.report)
        }
        
        out
        
    ##}else if(obj.class == "spMisalignLM"){
    }else if(inherits(sp.obj, "spMisalignLM")){
        
        X <- sp.obj$X
        Y <- sp.obj$Y
        m <- sp.obj$m
        obs.coords <- sp.obj$coords
        misalign.p <- sp.obj$misalign.p
        misalign.n <- sp.obj$misalign.n
        cov.model <- sp.obj$cov.model
        p.theta.samples <- sp.obj$p.theta.samples
        n.samples <- nrow(p.theta.samples)
        nugget <- sp.obj$nugget
        beta.prior <- sp.obj$beta.prior
        beta.Norm <- sp.obj$beta.Norm
        
        ##get prediction covariates and coordinates
        if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
        if(any(!is.list(pred.covars), length(pred.covars) != m)){stop(paste("error: pred.covars must be a list of length have ",m,"\n"))}
        if(!all(unlist(lapply(pred.covars, is.matrix)))){stop("error: pred.covars must be a list of matrices\n")}
        
        misalign.n.pred <- unlist(lapply(pred.covars, nrow))
        misalign.p.pred <- unlist(lapply(pred.covars, ncol))
        X.pred <- do.call(adiag, pred.covars)
        
        if(!identical(misalign.p, misalign.p.pred)){stop("error: the number of columns in the model matrices do not match those listed in pred.covars\n")}
        
        if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
        if(any(!is.list(pred.coords), length(pred.coords) != m)){stop(paste("error: pred.coords must be a list of length have ",m,"\n"))}
        if(!all(unlist(lapply(pred.coords, is.matrix)))){stop("error: pred.coords must be a list of matrices\n")}
        
        if(!identical(unlist(lapply(pred.coords, nrow)), misalign.n.pred)){stop("error: number of rows in the pred.coords matrices do not match those listed in pred.covars\n")}
        if(any(unlist(lapply(pred.coords, ncol)) != 2)){stop("error: all matrices in pred.coords must have two columns\n")}
        
        pred.coords <- as.matrix(do.call(rbind, pred.coords))
        
        ##subsamples
        if(missing(end)){end <- n.samples}
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        s.indx <- seq(start, end, by=as.integer(thin))
        
        p.theta.samples <- t(p.theta.samples[s.indx,,drop=FALSE])##note, for flat we could use the p.theta.recover.samples
        n.samples <- ncol(p.theta.samples)
        
        ##recover beta if needed (note, beta samples not needed for beta normal)
        beta <- NULL
        if(beta.prior == "flat"){
            beta <- t(spRecover(sp.obj, get.beta=TRUE, get.w=FALSE, start=start, end=end, thin=thin)$p.beta.recover.samples)
        }
        
        pred.obs.D <- iDist(pred.coords, obs.coords)
        obs.D <- iDist(obs.coords)
        
        storage.mode(X) <- "double"
        storage.mode(Y) <- "double"
        storage.mode(m) <- "integer"
        storage.mode(misalign.n) <- "integer"
        storage.mode(misalign.p) <- "integer"
        storage.mode(p.theta.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(beta) <- "double"
        storage.mode(nugget) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(misalign.n.pred) <- "integer"
        storage.mode(misalign.p.pred) <- "integer"
        storage.mode(pred.obs.D) <- "double"
        storage.mode(obs.D) <- "double"
        
        Z <- X.pred
        storage.mode(Z) <- "double"
        
        out <- .Call("spMisalignPredict", Y, X, m, misalign.n, misalign.p,
                     Z, misalign.n.pred, misalign.p.pred,
                     obs.D, pred.obs.D, 
                     p.theta.samples, beta, n.samples,
                     beta.prior, beta.Norm,
                     nugget, cov.model,
                     verbose, n.report)
        
        
        out$p.beta.samples.recover <- mcmc(t(beta))
        out
        
    ##}else if(obj.class == "spMisalignGLM"){
    }else if(inherits(sp.obj, "spMisalignGLM")){
        
        family <- sp.obj$family
        X <- sp.obj$X
        Y <- sp.obj$Y
        m <- sp.obj$m
        obs.coords <- sp.obj$coords
        misalign.p <- sp.obj$misalign.p
        misalign.n <- sp.obj$misalign.n
        cov.model <- sp.obj$cov.model
        p.w.samples <- sp.obj$p.w.samples
        p.beta.theta.samples <- sp.obj$p.beta.theta.samples
        n.samples <- nrow(p.beta.theta.samples)
        
        ##get prediction covariates and coordinates
        if(missing(pred.covars)){stop("error: pred.covars must be specified\n")}
        if(any(!is.list(pred.covars), length(pred.covars) != m)){stop(paste("error: pred.covars must be a list of length have ",m,"\n"))}
        if(!all(unlist(lapply(pred.covars, is.matrix)))){stop("error: pred.covars must be a list of matrices\n")}
        
        misalign.n.pred <- unlist(lapply(pred.covars, nrow))
        misalign.p.pred <- unlist(lapply(pred.covars, ncol))
        X.pred <- do.call(adiag, pred.covars)
        
        if(!identical(misalign.p, misalign.p.pred)){stop("error: the number of columns in the model matrices do not match those listed in pred.covars\n")}
        
        if(missing(pred.coords)){stop("error: pred.coords must be specified\n")}
        if(any(!is.list(pred.coords), length(pred.coords) != m)){stop(paste("error: pred.coords must be a list of length have ",m,"\n"))}
        if(!all(unlist(lapply(pred.coords, is.matrix)))){stop("error: pred.coords must be a list of matrices\n")}
        
        if(!identical(unlist(lapply(pred.coords, nrow)), misalign.n.pred)){stop("error: number of rows in the pred.coords matrices do not match those listed in pred.covars\n")}
        if(any(unlist(lapply(pred.coords, ncol)) != 2)){stop("error: all matrices in pred.coords must have two columns\n")}
        
        pred.coords <- as.matrix(do.call(rbind, pred.coords))
        
        ##subsamples
        if(missing(end)){end <- n.samples}
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        s.indx <- seq(start, end, by=as.integer(thin))
        
        p.w.samples <- p.w.samples[,s.indx,drop=FALSE]
        
        p.samples <- t(p.beta.theta.samples[s.indx,,drop=FALSE])
        n.samples <- ncol(p.samples)
        
        pred.obs.D <- iDist(pred.coords, obs.coords)
        obs.D <- iDist(obs.coords)
        
        storage.mode(X) <- "double"
        storage.mode(Y) <- "double"
        storage.mode(m) <- "integer"
        storage.mode(misalign.n) <- "integer"
        storage.mode(misalign.p) <- "integer"
        storage.mode(p.samples) <- "double"
        storage.mode(p.w.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(misalign.n.pred) <- "integer"
        storage.mode(misalign.p.pred) <- "integer"
        storage.mode(pred.obs.D) <- "double"
        storage.mode(obs.D) <- "double"
        
        Z <- X.pred
        storage.mode(Z) <- "double"
        
        out <- .Call("spMisalignGLMPredict", family, Y, X, m, misalign.n, misalign.p,
                     Z, misalign.n.pred, misalign.p.pred,
                     obs.D, pred.obs.D, 
                     p.samples, p.w.samples, n.samples,
                     cov.model,
                     verbose, n.report)
        out
        
    }else{
        stop("error: wrong class\n")
    }
    
}
