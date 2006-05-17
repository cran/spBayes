"ggt.sp" <- function(formula, data = parent.frame(), coords, run.control,
                     beta.control, cov.model, var.update.control,
                     DIC=TRUE, DIC.start = 1, verbose=TRUE, ...){

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
  ##verbose
  ####################################################
  if(!is.logical(verbose)){stop("error: verbose must be of type logical")}
  
  ##start report
  if(verbose){
    cat("\n------------------------------------------------------------------------\n")
    cat("\t\t\tModel specifications")
    cat("\n------------------------------------------------------------------------\n")
  }

  ####################################################
  ##run
  ####################################################
  if(missing(run.control)){stop("error: run.control must be specified")}
  if(!is.list(run.control)){stop("error: run.control must be a list")}
  if(!"n.samples"%in%names(run.control)){stop("error: n.samples must be specified in run.control list")}
  n.samples <- run.control$n.samples
  if(!is.numeric(n.samples) || n.samples <= 0){stop("error: n.samples must be numeric and greater than zero")}

  ##make the run list
  run <- list("n.samples" = as.integer(n.samples))

  ####################################################
  ##DIC
  ####################################################
  if(!is.logical(DIC)){stop("error: DIC must be of type logical")}
  if(DIC){
    if(!is.numeric(DIC.start)){stop("error: DIC.start must be of type logical")}
    if(DIC.start < 1 || DIC.start >= as.integer(n.samples))
      {stop("error: DIC.start is not reasonable, DIC.start should be > 1 and < n.samples")}
    
    DIC.start <- as.integer(DIC.start)
  }else{
    DIC.start <- NULL
  }
  
  ####################################################
  ##covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}

  ####################################################
  ##Beta
  ####################################################
  if(missing(formula)){stop("error: model formula must be a specified\n")}
   
  ##form response and model matrices
  holder <- parse.formula(formula, data)
  Y <- holder[[1]]
  X <- holder[[2]]
  xnames <- holder[[3]]
  
  ##get the beta.control stuff
  if(missing(beta.control)){stop("error: beta.control must be specified\n")}
  if(!is.list(beta.control)){stop("error: beta.control must be list\n")}

  ##get the method for updating the beta
  if(!"beta.update.method"%in%names(beta.control)){
    if(verbose)
      cat("note: 'beta.update.method' is missing, updates will be made with Gibbs\n\n")
    update.beta.method <- "GIBBS"
  }else if(beta.control$beta.update.method%in%c("mh","MH","metrop","Metrop")){
    update.beta.method <- "MH"
  }else if(beta.control$beta.update.method%in%c("gibbs","Gibbs","GIBBS")){
    update.beta.method <- "GIBBS"
  }else{
    stop("error: beta.update.method in the beta.control list is not specified correctly, it should be set to 'gibbs' or 'metorp'\n")
  }

  ##get the tuning matrix if metrop is the update method
  if(update.beta.method == "MH"){
    ##get the Cholesky factored tuning matrix
    if(!"tuning.matrix"%in%names(beta.control)){
      if(verbose)
        cat("note: Beta tuning matrix 'tuning.matrix' is missing, will use the default tuning matrix chol(vcov(lm(Y~X)))\n\n")
      if(ncol(X) == 1){
        tuning.matrix <- t(chol(vcov(lm(Y ~ 1)))) ##recall I want the lower
      }else{
        tuning.matrix <- t(chol(vcov(lm(Y ~ X)))) ##recall I want the lower
      }
    }else{
      if(is.matrix(beta.control$tuning.matrix)){
        tuning.matrix <- beta.control$tuning.matrix ##just a temp for the next test
        
        if(dim(tuning.matrix)[1] == dim(tuning.matrix)[2] && dim(tuning.matrix)[1] == ncol(X)){
          ##recall I want the lower, assume they give the upper with chol()
          tuning.matrix <- t(beta.control$tuning.matrix)
        }else{
          stop("error: beta tuning matrix is not square or does not have the same dimension as the number of beta in the model formula\n")
        }
      }else{
        stop("error: beta tuning matrix is not a matrix\n")
      }     
    }
  }else{
    ##just set the tuning matrix to NULL so I don't have to send two different beta lists to ggt
    tuning.matrix <- NULL
  }

  ##get the starting values
  if(!"beta.starting"%in%names(beta.control)){
    beta.starting <- rep(0.0,ncol(X))
    if(verbose)
      cat("note: Beta staring values 'beta.starting' is missing, will using default starting values of 0 for each beta\n\n")
  }else{
    beta.starting <- as.double(beta.control$beta.starting)
    if(length(beta.starting) != ncol(X)){
      stop(paste("error: number of beta starting values (", length(beta.starting),") is different than the number of specified covariates (",ncol(X),")\n",sep=""))
    }
  }

  ##get prior if needed
  if("beta.prior"%in%names(beta.control)){
    if(class(beta.control$beta.prior) == "ggt.prior" && beta.control$beta.prior$dist %in% c("NORMAL")){
      beta.prior <- "NORMAL"
      beta.prior.mu <- as.matrix(beta.control$beta.prior$params$mu)
      beta.prior.precision <- as.matrix(beta.control$beta.prior$params$precision)
      
      ##check the hyperparameters
      if(nrow(beta.prior.precision) != ncol(beta.prior.precision))
        stop("error: precision matrix for the beta's normal prior is not square\n")
      
      if(nrow(beta.prior.mu) != ncol(beta.prior.precision))
        stop("error: for the specified normal prior on the beta parameters, the number of rows in the mu matrix does not equal the number of rows in the square precision matrix\n")
      
      if((ncol(beta.prior.mu) != 1) || (nrow(beta.prior.mu) != ncol(X)))
        stop("error: mu matrix for the beta's normal prior should be 1 x p, where p is the number of beta parameters")
            
    }else{
      stop("error: the prior on the beta parameters is not specified correctly or is not an valid prior distribution\n")
    }
  }else{ ##assume flat
    beta.prior <- "FLAT"
    beta.prior.mu <- NULL
    beta.prior.precision <- NULL
  }


  ##get sample.order
  if("sample.order"%in%names(beta.control)){
    if(update.beta.method == "GIBBS"){
      if(verbose)
        cat("note: 'sample.order' is specified in the 'beta.control' list; however, sample order only applies when 'beta.update.method' is Metropolis-Hastings (MH)\n")
    }
    
    beta.sample.order <- beta.control$sample.order
  }else{
    if(update.beta.method == "MH"){
      if(verbose)
        cat("note: 'beta.update.method' is Metropolis-Hastings (MH), but no 'sample.order' is specified in the 'beta.control' list; therefore, relative sample order is set to 0\n")

      beta.sample.order <- 0
    }else{##else update.beta.method == "GIBBS"
      beta.sample.order <- NULL
    }
      
  }

  
  ##make the beta list
  beta <- list("Y"=Y, "X"=as.matrix(as.data.frame(X)),
               "fixed.effects.starting"=beta.starting,
               "update.theta.method"=update.beta.method, "tuning.matrix"=tuning.matrix, "prior"=beta.prior,
               "prior.mu"=beta.prior.mu, "prior.precision"=beta.prior.precision, "sample.order"=beta.sample.order)
 
  ####################################################
  ##coordinate matrix
  ####################################################
  if(missing(coords)){stop("error: coords must be specified")}
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(ncol(coords) != 2 || nrow(coords) != nrow(Y)){
    stop("error: either the coords have more than two columns or then number of rows is different than data used in the model formula")
  }
  
  ####################################################
  ##variance
  ####################################################
  if(missing(var.update.control)){stop("error: var.update.control must be specified")}
  if(!is.list(var.update.control)){stop("error: var.update.control must be list")}

  ##check for the parameters required for the specified cov.model
  ##sigma^2, phi, nu
  var.params <- names(var.update.control)
  
  if(!"sigmasq"%in%var.params)
    stop("error: sigmasq must be included in the var.update.control list.")

  if(!"phi"%in%var.params)
    stop("error: phi must be included in the var.update.control list.")
  
  if(!"nu"%in%var.params && cov.model=="matern")
    stop("error: because the matern cov.model was selected, the nu parameter must be included in the var.update.control list.")

  if("nu"%in%var.params && cov.model!="matern")
    stop("error: the selected cov.model '",cov.model,"' does not use the nu parameter, remove nu from the var.update.control list.")

  no.tausq <- FALSE
  if(!"tausq"%in%var.params)
    no.tausq <- TRUE

  ##
  ##check for fixed parameters and add to fixed.var.params list
  ##

  ##make sure they didn't put 'fixed' as a value not a tag
  for(i in 1:length(var.update.control)){
    tmp <- var.update.control[[i]]
    for(j in 1:length(tmp)){
      if(tmp[j]=="fixed"){
        stop("error: a fixed parameter in the var.update.control list must have 'fixed' as a list tag and associated value (e.g., phi = list(fixed=0.1))\n")
      }
    }
  }
  
  ##
  ##add parameters to their respective lists fixed.var.parmas and var.update.params
  ##
  fixed.var.params <- list()

  ##get just the fixed parameters
  for(i in 1:length(var.update.control)){
    if("fixed"%in%names(var.update.control[[i]])){
      fixed.var.params <- c(fixed.var.params, var.update.control[i])
    }
  }
  
  ##get just the free parameters
  tmp <- list()
  for(i in 1:length(var.update.control)){
    if(!"fixed"%in%names(var.update.control[[i]])){
      tmp <- c(tmp, var.update.control[i])
    }
  }
  var.update.control <- tmp
  
  ##make sure there is at least one free parameter
  if(length(var.update.control) == 0)
    stop("error: at least one parameter in the var.update.control list must be free (i.e., not fixed)\n")

  ##check if all the update information is in the parameters' list
  ##get the MH sampling order
  cnt <- 0
  for(i in 1:length(var.update.control)){
    if("sample.order"%in%names(var.update.control[[i]])){
      cnt <- cnt+1
    }
  }
  
  if(cnt != length(var.update.control)){
      if(verbose){
        cat("note: 'sample.order' is not specified for the parameters in the var.update.control list, or not specified for all parameters; therefore, single block updating will be used in the Metropolis-Hastings\n\n")
      }
      for(i in 1:length(var.update.control)){
        var.update.control[[i]]$sample.order <- 0
      }
  }

  ##check for starting values
  for(i in 1:length(var.update.control)){
    if(!"starting"%in%names(var.update.control[[i]])){
      stop("error: 'starting' value is not specified for at least one of the non-fixed parameters defined in the var.update.control list\n")
    }
  }

  ##check for tuning values
  for(i in 1:length(var.update.control)){
    if(!"tuning"%in%names(var.update.control[[i]])){
      stop("error: 'tuning' value is not specified for at least one of the non-fixed parameters defined in the var.update.control list\n")
    }
    
    if(var.update.control[[i]]$tuning <= 0){
      stop("error: at least one of the tuning values associated with a non-fixed parameter in the var.update.control list is <=0, tuning must be > 0\n")
    }
  }

  ##check for prior
  for(i in 1:length(var.update.control)){
    if(!"prior"%in%names(var.update.control[[i]])){
      stop("error: 'prior' is not specified for at least one of the non-fixed parameters defined in the var.update.control list\n")
    }
    
    if(class(var.update.control[[i]]$prior) != "ggt.prior" || !var.update.control[[i]]$prior$dist %in% c("IG", "UNIF", "LOGUNIF", "HC")){
      stop(" error: one of the 'priors' associated with a non-fixed parameter in the var.update.control is not specified correctly or is not an valid prior distribution for this parameter\n")
    }
  }

  ##
  ##make the random list
  ##
  var <- list("var.update.control" = var.update.control,
              "fixed.var.params" = fixed.var.params,
              "no.tausq" = as.integer(no.tausq), "env"=.GlobalEnv)
  
  ####################################################
  ##report
  ####################################################
  if(verbose){
    cat("---------------------------------------------------------\n")
    ##show the covariates
    cat("Number of observations: ",nrow(X),"\n\n")
    print.model.matrix(X)
    cat("Beta will be updated using: ")

    update.beta.method
    if(update.beta.method=="MH"){
      cat("Metropolis-Hastings\n\n")
    }else{
      cat("Gibbs\n\n")
    }
    cat("Beta prior is: ", beta.prior,"\n\n")

    cat("Beta starting values (same order as formula covariates):\n")
    for(i in beta.starting)
      cat("\t",i)
    cat("\n\n")

    if(update.beta.method=="MH")
      cat("Beta sample order (relative): ", beta.sample.order,"\n\n")
        
    cat("Covariance model: ",cov.model,"\n\n")
    
    cat("Covariance model parameters to be estimated:\n")
      for(i in 1:length(var.update.control)){
        cat("\t",names(var.update.control[i]),":\n")
        
        cat("\t\tsample order (relative):\t\t",var.update.control[[i]]$sample.order,"\n")
        cat("\t\tstarting value:\t\t\t\t",var.update.control[[i]]$starting,"\n")
        cat("\t\tMetropolis-Hastings tuning value:\t",var.update.control[[i]]$tuning,"\n")
        cat("\t\tprior:\t\t\t\t\t",var.update.control[[i]]$prior[["dist"]],"\n")
      }
    
    if(length(fixed.var.params)){
      cat("\nCovariance model fixed parameters and values:\n\t")
      for(i in 1:length(fixed.var.params)){
        cat(names(fixed.var.params[i]),"=",fixed.var.params[[i]]$fixed," ")
      }
    }

    if(no.tausq){
      cat("\n\nNote: tau^2 was not specified as fixed or free so it will not be used in the covariance model\n")
    }

  }
  
  ####################################################
  ##set up args and off it goes
  ####################################################
  args <- list("run" = run, "fixed"=beta, "cov.model" =cov.model,
               "coords"=coords, "var" = var, "DIC"=as.integer(DIC),
               "DIC.start"=as.integer(DIC.start), "verbose"=as.integer(verbose))
  
  out <- .Call("ggtsp", args)

  out1 <- list()
  out1$cov.model <- cov.model
  out1$no.tausq <- as.logical(no.tausq)
  out1$coords <- coords
  out1$X <- X
  out1$Y <- Y
  out1$p.samples <- cbind(t(out$fixedEffectSamples),out$varParameterSamples)
  out1$acceptance <- out$acceptance
  if(DIC)
    out1$DIC <- out$DIC
  class(out1) <- c("ggt.sp")
  out1
}
         
