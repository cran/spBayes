spRecover <- function(sp.obj, get.beta=TRUE, get.w=TRUE, start=1, end, thin=1, verbose=TRUE, n.report=100, n.omp.threads=1, ...){
  
  ####################################################
  ##Check for unused args
  ####################################################
    formal.args <- names(formals(sys.function(sys.parent())))
    elip.args <- names(list(...))
    for(i in elip.args){
        if(! i %in% formal.args)
            warning("'",i, "' is not an argument")
    }
    
    if(missing(sp.obj)){stop("error: spRecover expects sp.obj\n")}
    ##if(!class(sp.obj) %in% c("spLM", "spMvLM", "spMisalignLM", "spSVC")){
    if(!inherits(sp.obj, c("spLM", "spMvLM", "spMisalignLM", "spSVC"))){
        stop("error: requires an output object of class spLM, spMvLM, spMisalignLM, or spSVC\n")
    }
    
    ##obj.class <- class(sp.obj)
    
    ##need beta to get.w
    if(get.w){
        get.beta <- TRUE
    }
    
    ##if(obj.class == "spSVC"){
    if(inherits(sp.obj, "spSVC")){
            
        Y <- sp.obj$Y
        X <- sp.obj$X
        Z <- sp.obj$Z
        p <- ncol(X)
        n <- nrow(X)
        m <- ncol(Z)
        n.ltr <- m*(m+1)/2
        K.diag <- sp.obj$K.diag
        coords <- sp.obj$coords
        q <- ncol(coords)
        p.theta.samples <- sp.obj$p.theta.samples
        n.samples <- nrow(p.theta.samples)
        cov.model <- sp.obj$cov.model
        nugget <- sp.obj$nugget
        beta.prior <- sp.obj$beta.prior
        beta.Norm <- sp.obj$beta.Norm
        x.names <- sp.obj$x.names
        svc.cols <- sp.obj$svc.cols
        
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
        
        if(missing(end)){
            end <- n.samples
        }
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        p.theta.samples <- p.theta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE]
        n.samples <- nrow(p.theta.samples)
        
        ##for now no w for no nugget model
        if(!nugget && get.w){
            warning("Recovery of w is not available for the no nugget model")
            get.w <- FALSE
        }
        
        storage.mode(Y) <- "double"
        storage.mode(X) <- "double"
        storage.mode(Z) <- "double"
        storage.mode(p) <- "integer"
        storage.mode(m) <- "integer"
        storage.mode(n) <- "integer"
        storage.mode(coords) <- "double"
        storage.mode(q) <- "integer"
        storage.mode(K.diag) <- "integer"
        storage.mode(p.theta.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(A.indx) <- "integer"
        storage.mode(tau.sq.indx) <- "integer"
        storage.mode(phi.indx) <- "integer"
        storage.mode(nu.indx) <- "integer"
        storage.mode(nugget) <- "integer"
        storage.mode(get.beta) <- "integer"
        storage.mode(get.w) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        storage.mode(n.omp.threads) <- "integer"
        
        out <- .Call("spSVCRecover", Y, X, p, Z, m, n, coords, q, K.diag,
                     p.theta.samples, n.samples, 
                     A.indx, tau.sq.indx, phi.indx, nu.indx,
                     beta.prior, beta.Norm, 	   
                     nugget, cov.model,
                     get.beta, get.w,
                     verbose, n.report, n.omp.threads)
        
        rownames(out$p.beta.samples) <- x.names
        sp.obj$p.beta.recover.samples <- mcmc(t(out$p.beta.samples))
        sp.obj$p.theta.recover.samples <- mcmc(p.theta.samples)
        sp.obj$recover.args <- c(start, end, thin)
        
        if(get.w){
            sp.obj$p.w.recover.samples <- out$p.w.samples
            
            ##also save as list with elements corresponding to svc.cols
            sp.obj$p.w.recover.samples.list  <- vector("list", length = m)
            names(sp.obj$p.w.recover.samples.list) <- paste0("w.",x.names[svc.cols])
            indx <- rep(1:m, n)
            for(i in 1:m){
                sp.obj$p.w.recover.samples.list[[i]] <- sp.obj$p.w.recover.samples[indx == i,]
            }
            
            ##B.tilde = B + w for each svc.cols
            sp.obj$p.tilde.beta.recover.samples.list  <- vector("list", length = m)
            names(sp.obj$p.tilde.beta.recover.samples.list) <- paste0("tilde.beta.",x.names[svc.cols])
            for(i in 1:m){
                sp.obj$p.tilde.beta.recover.samples.list[[i]] <- sp.obj$p.w.recover.samples.list[[i]] + matrix(rep(as.matrix(sp.obj$p.beta.recover.samples)[,svc.cols[i]], each=n), nrow=n)
            }
        }
        
        ##y fitted
        X.tilde <- t(bdiag(as.list(as.data.frame(t(X[,svc.cols,drop=FALSE])))))    
        
        sp.obj$p.y.samples <- matrix(0, n, n.samples)
        for(i in 1:n.samples){
            sp.obj$p.y.samples[,i] <- rnorm(n, (X%*%sp.obj$p.beta.recover.samples[i,] + X.tilde%*%sp.obj$p.w.recover.samples[,i])[,1], sqrt(sp.obj$p.theta.recover.samples[i,"tau.sq"]))
        }
        
        class(sp.obj) <- "spSVC"
        
        sp.obj

    ##}else if(obj.class == "spLM"){
    }else if(inherits(sp.obj, "spLM")){        

        Y <-sp.obj$Y
        X <- sp.obj$X
        p <- ncol(X)
        n <- nrow(X)
        coords <- sp.obj$coords
        p.theta.samples <- sp.obj$p.theta.samples
        n.samples <- nrow(p.theta.samples)
        n.params <- ncol(p.theta.samples)
        cov.model <- sp.obj$cov.model
        nugget <- sp.obj$nugget
        beta.prior <- sp.obj$beta.prior
        beta.Norm <- sp.obj$beta.Norm
        x.names <- sp.obj$x.names
        is.pp <- sp.obj$is.pp
        modified.pp <- sp.obj$modified.pp
        
        coords.D <- iDist(coords)
        
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
        
        if(missing(end))
            end <- n.samples
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        p.theta.samples <- p.theta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE]
        n.samples <- nrow(p.theta.samples)
        
        ##if pp
        p.beta.samples <- NULL
        knots.D <- NULL
        knots.obs.D <- NULL
        m <- NULL
        if(is.pp){
            p.beta.samples <- sp.obj$p.beta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE]
            knots.D <- iDist(sp.obj$knot.coords)
            knots.obs.D <- iDist(sp.obj$knot.coords, coords)
            m <- nrow(knots.D)
        }
        
        storage.mode(Y) <- "double"
        storage.mode(X) <- "double"
        storage.mode(p) <- "integer"
        storage.mode(n) <- "integer"
        storage.mode(coords.D) <- "double"
        storage.mode(p.theta.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(n.params) <- "integer"
        storage.mode(sigma.sq.indx) <- "integer"
        storage.mode(tau.sq.indx) <- "integer"
        storage.mode(phi.indx) <- "integer"
        storage.mode(nu.indx) <- "integer"
        storage.mode(nugget) <- "integer"
        storage.mode(get.beta) <- "integer"
        storage.mode(get.w) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        
        storage.mode(m) <- "integer"
        storage.mode(is.pp) <- "integer"
        storage.mode(modified.pp) <- "integer"
        storage.mode(knots.D) <- "double"
        storage.mode(knots.obs.D) <- "double"
        storage.mode(p.beta.samples) <- "double"
        
        if(is.pp){
            out <- .Call("spPPLMRecover", X, Y, n, p, m, 
                         p.theta.samples, n.samples, 
                         p.beta.samples, sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                         nugget, knots.D, knots.obs.D, cov.model, modified.pp, verbose, n.report)
            
            sp.obj$p.beta.recover.samples <- mcmc(p.beta.samples)
            sp.obj$p.theta.recover.samples <- mcmc(p.theta.samples)
            
            sp.obj$p.w.recover.samples <- out$p.w.samples
            sp.obj$p.wStr.recover.samples <- out$p.wStr.samples
            
            
        }else{
            out <- .Call("spLMRecover", Y, X, p, n, coords.D,
                         p.theta.samples, n.samples, 
                         sigma.sq.indx, tau.sq.indx, phi.indx, nu.indx,
                         beta.prior, beta.Norm, 	   
                         nugget, cov.model,
                         get.beta, get.w,
                         verbose, n.report)
            
            rownames(out$p.beta.samples) <- x.names
            sp.obj$p.beta.recover.samples <- mcmc(t(out$p.beta.samples))
            sp.obj$p.theta.recover.samples <- mcmc(p.theta.samples)
            
            if(get.w){
                sp.obj$p.w.recover.samples <- out$p.w.samples
            }
            
        }
        
        class(sp.obj) <- "spLM"
        
        sp.obj
        
    ##}else if(obj.class == "spMvLM"){
    }else if(inherits(sp.obj, "spMvLM")){
        
        Y <-sp.obj$Y
        X <- sp.obj$X
        p <- ncol(X)
        m <- sp.obj$m
        coords <- sp.obj$coords
        n <- nrow(coords)
        p.theta.samples <- sp.obj$p.theta.samples
        n.samples <- nrow(p.theta.samples)
        cov.model <- sp.obj$cov.model
        nugget <- sp.obj$nugget
        Psi.diag <- sp.obj$Psi.diag
        beta.prior <- sp.obj$beta.prior
        beta.Norm <- sp.obj$beta.Norm
        x.names <- sp.obj$x.names
        is.pp <- sp.obj$is.pp
        modified.pp <- sp.obj$modified.pp
        
        coords.D <- iDist(coords)
        
        if(missing(end))
            end <- n.samples
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        p.theta.samples <- t(p.theta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE])
        n.samples <- ncol(p.theta.samples)
        
        ##if pp
        p.beta.samples <- NULL
        knots.D <- NULL
        knots.obs.D <- NULL
        if(is.pp){
            p.beta.samples <- sp.obj$p.beta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE]
            knots.D <- iDist(sp.obj$knot.coords)
            knots.obs.D <- iDist(sp.obj$knot.coords, coords)
            g <- nrow(knots.D)
        }
        
        storage.mode(Y) <- "double"
        storage.mode(X) <- "double"
        storage.mode(p) <- "integer"
        storage.mode(n) <- "integer"
        storage.mode(m) <- "integer"
        storage.mode(coords.D) <- "double"
        storage.mode(p.theta.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(nugget) <- "integer"
        storage.mode(Psi.diag) <- "integer"
        storage.mode(get.beta) <- "integer"
        storage.mode(get.w) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        
        if(is.pp){
            storage.mode(g) <- "integer"
            storage.mode(knots.D) <- "double"
            storage.mode(knots.obs.D) <- "double"
            storage.mode(modified.pp) <- "integer"
            
            out <- .Call("spPPMvLMRecover", X, Y, n, m, g, p,
                         knots.D, knots.obs.D,
                         p.theta.samples, p.beta.samples, n.samples,   
                         nugget, Psi.diag, cov.model,
                         modified.pp, verbose, n.report)
            
            sp.obj$p.beta.recover.samples <- mcmc(p.beta.samples)
            sp.obj$p.theta.recover.samples <- mcmc(t(p.theta.samples))
            
            sp.obj$p.w.recover.samples <- out$p.w.samples
            sp.obj$p.wStr.recover.samples <- out$p.wStr.samples
            
        }else{
            
            out <- .Call("spMvLMRecover", Y, X, p, n, m, coords.D,
                         p.theta.samples, n.samples,
                         beta.prior, beta.Norm, 	   
                         nugget, Psi.diag, cov.model,
                         get.beta, get.w,
                         verbose, n.report)
            
            rownames(out$p.beta.samples) <- x.names
            sp.obj$p.beta.recover.samples <- mcmc(t(out$p.beta.samples))
            sp.obj$p.theta.recover.samples <- mcmc(t(p.theta.samples))
            
            if(get.w){
                sp.obj$p.w.recover.samples <- out$p.w.samples
            }
            
        }
        
        class(sp.obj) <- "spMvLM"
        
        sp.obj
        
    ##}else if(obj.class == "spMisalignLM"){
    }else if(inherits(sp.obj, "spMisalignLM")){
        
        Y <-sp.obj$Y
        X <- sp.obj$X
        m <- sp.obj$m
        coords <- sp.obj$coords
        misalign.p <- sp.obj$misalign.p
        misalign.n <- sp.obj$misalign.n
        p.theta.samples <- sp.obj$p.theta.samples
        n.samples <- nrow(p.theta.samples)
        cov.model <- sp.obj$cov.model
        nugget <- sp.obj$nugget
        beta.prior <- sp.obj$beta.prior
        beta.Norm <- sp.obj$beta.Norm
        x.names <- sp.obj$x.names
        
        coords.D <- iDist(coords)
        
        if(missing(end))
            end <- n.samples
        
        if(!is.numeric(start) || start >= n.samples)
            stop("error: invalid start")
        if(!is.numeric(end) || end > n.samples) 
            stop("error: invalid end")
        if(!is.numeric(thin) || thin >= n.samples) 
            stop("error: invalid thin")
        
        p.theta.samples <- t(p.theta.samples[seq(start, end, by=as.integer(thin)),,drop=FALSE])
        n.samples <- ncol(p.theta.samples)
        
        storage.mode(Y) <- "double"
        storage.mode(X) <- "double"
        storage.mode(misalign.p) <- "integer"
        storage.mode(misalign.n) <- "integer"
        storage.mode(m) <- "integer"
        storage.mode(coords.D) <- "double"
        storage.mode(p.theta.samples) <- "double"
        storage.mode(n.samples) <- "integer"
        storage.mode(nugget) <- "integer"
        storage.mode(get.beta) <- "integer"
        storage.mode(get.w) <- "integer"
        storage.mode(verbose) <- "integer"
        storage.mode(n.report) <- "integer"
        
        out <- .Call("spMisalignRecover", Y, X, misalign.p, misalign.n, m, coords.D,
                     p.theta.samples, n.samples,
                     beta.prior, beta.Norm, 	   
                     nugget, cov.model,
                     get.beta, get.w,
                     verbose, n.report)
        
        rownames(out$p.beta.samples) <- x.names
        sp.obj$p.beta.recover.samples <- mcmc(t(out$p.beta.samples))
        sp.obj$p.theta.recover.samples <- mcmc(t(p.theta.samples))
        
        if(get.w){
            sp.obj$p.w.recover.samples <- out$p.w.samples
        }
        
        class(sp.obj) <- "spMisalignLM"
        
        sp.obj
        
    }else{
        stop("error: wrong class\n")
    }
    
}
