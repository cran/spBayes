spGGT <- function(formula, data = parent.frame(), coords, run.control, var.update.control, beta.update.control,
                   cov.model, verbose=TRUE, ...){
  
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
  ##run control
  ####################################################
  if(missing(run.control)){stop("error: run.control must be specified")}
  
  if(!"n.samples"%in%names(run.control))
    {stop("error: 'n.samples' must be specified in run.control")}
  
  n.samples <- as.integer(run.control$n.samples)
 
  if(n.samples <= 0)
    stop("error: 'n.samples' must be > 0")

  if(!"linpack"%in%names(run.control))##make default TRUE
    linpack <- as.integer(TRUE)
  else
    linpack <- as.integer(run.control$linpack)

  if(!"sp.effects"%in%names(run.control))
    sp.effects <- as.integer(FALSE)
  else
    sp.effects <- as.integer(run.control$sp.effects)

  ######################
  ##DIC in run control
  ######################
  if(!"DIC"%in%names(run.control))
    DIC <- as.integer(FALSE)
  else
    DIC <- as.integer(run.control$DIC)

  ##just ignore on c++ side if DIC == FALSE
  if(!"DIC.start"%in%names(run.control))
    DIC.start <- 1
  else
    DIC.start <- run.control$DIC.start
  
  if(DIC.start <= 0 || DIC.start > n.samples)
    stop("error: 'DIC.start' must be > 0 and <= n.samples")

  ##change to c++ index
  DIC.start <- DIC.start-1
  
  DIC.start <- as.integer(DIC.start)
  
    
  ####################################################
  ##covariance model
  ####################################################
  if(missing(cov.model)){stop("error: cov.model must be specified")}
  if(!cov.model%in%c("gaussian","exponential","matern","spherical"))
    {stop("error: specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")}
  
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
  
  ##make sure storage mode is correct
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(m) <- "integer"
  n.beta <- ncol(X)

  ####################################################
  ##Distance matrix
  ####################################################
  if(missing(coords)){stop("error: coords must be specified")}
  if(!is.matrix(coords)){stop("error: coords must n-by-2 matrix of xy-coordinate locations")}
  if(ncol(coords) != 2 || nrow(coords) != nrow(Y)/m){
    stop("error: either the coords have more than two columns or then number of rows is different than data used in the model formula")
  }
  
  D <- as.matrix(dist(coords))
  storage.mode(D) <- "double"
  
  ####################################################
  ##variance parameter
  ####################################################
  if(missing(var.update.control)){stop("error: var.update.control must be specified")}
  if(!is.list(var.update.control)){stop("error: var.update.control must be list")}

  ##check for the parameters required for the specified cov.model
  ##K, phi, nu
  var.params <- names(var.update.control)
  
  if(!"K"%in%var.params)
    stop("error: K must be included in the var.update.control list.")

  if(!"phi"%in%var.params)
    stop("error: phi must be included in the var.update.control list.")
  
  if(!"nu"%in%var.params && cov.model=="matern")
    stop("error: because the matern cov.model was selected, the nu parameter must be included in the var.update.control list.")

  if("nu"%in%var.params && cov.model!="matern")
    stop("error: the selected cov.model '",cov.model,"' does not use the nu parameter, remove nu from the var.update.control list.")
  
  no.Psi <- FALSE
  if(!"Psi"%in%var.params)
    no.Psi <- TRUE

  ####################################################
  ##Parameter prep for matrix, vector, or scaler 
  ####################################################
  mkPar <- function(param.name, K.list, m){

    ##determine the case
    ##case 1: if m == 1 or if m > 1 and single non-fixed prior is given
    ##case 2: if m > 1 and prior is a list of length m
    ##case 3: if m > 1 and prior is WISH
    
    ##prior
    if(!"prior"%in%names(K.list))
      stop("error: 'prior' value is not specified for the ",param.name," parameter in the var.update.control list\n")
    
    if(class(K.list$prior) == "list"){
      if(m == 1)
        stop("error: misspecification for ",param.name," parameter in the var.update.control list, a only a single response model is specified but 'prior' is a list\n")

      if(length(K.list$prior) != m)
         stop("error: misspecification for ",param.name," parameter in the var.update.control list, 'prior' must be a list of length ",m,"\n")

      ##check each prior in the list
      for(i in K.list$prior){
        if(class(i) != "prior" || !i$dist %in% c("IG", "UNIF", "HC", "FIXED"))
          stop("error: 'prior' value is misspecified for the ",param.name," parameter in the var.update.control list. It seems you are trying to specify independent priors for each element on the diag(",param.name,"); however, one or more of these priors is misspecified, they must IG, UNIF, HC, or FIXED\n")
      }
      K.list$case <- 2
    }else if(class(K.list$prior) == "prior"){
      if(m == 1){
        if(K.list$prior$dist %in% c("IG", "UNIF", "HC", "FIXED")){
          K.list$case <- 1
        }else{
          stop("error: 'prior' value is misspecifed for the ",param.name," parameter in the var.update.control list.  For the univariate response model the prior must be IG, UNIF, HC, or FIXED\n")
        }
      }else if(m > 1){

        ##two options based on the starting value, if dist is IWISH then starting must be a mxm matrix
        ##if dist is in "IG", "UNIF", "HC" then starting must be a scaler
        ##now if dist is FIXED then the starting data type determines which case we have
        ##if dist is FIXED and starting is matrix then assumed fixed matrix
        ##if dist is Fixed and starting is scaler then assumed single parame for all diag elements

        ##first check for starting in list.
        if(!"starting"%in%names(K.list))
          stop("error: 'starting' value is not specified for the ",param.name," parameter is not in the var.update.control list\n")
              
        if(K.list$prior$dist == "IWISH"){
          K.list$case <- 3
        }else if(K.list$prior$dist%in%c("IG", "UNIF", "HC")){
          K.list$case <- 1
        }else if(K.list$prior$dist == "FIXED"){
          if(is.matrix(K.list$starting)){
            if(dim(K.list$starting)[1] != dim(K.list$starting)[2] || dim(K.list$starting)[1] != m)
              stop("error: 'starting' value must be a ",m,"x",m," matrix as ",param.name,"'s 'prior' is FIXED and 'starting' is a matrix")
            K.list$case <- 3
            
          }else if(is.numeric(K.list$starting) || length(K.list$starting) == 1){
            K.list$case <- 1
          }else{
            stop("error: 'prior' value is misspecified for the ",param.name," parameter in the var.update.control list, with a multivariate response and FIXED 'prior' the 'starting' must be either a mxm matrix or scaler depending on which specification you want\n")
          }
        }
      }else{
        stop("error: 'prior' value is misspecifed for the ",param.name," parameter in the var.update.control list.\n")
      }
    }else{
      stop("error: 'prior' value is misspecified for the ",param.name," parameter in the var.update.control list.\n")
    }

    
    ##make list for K and assign case
    K <- list()
    K$case <- as.integer(K.list$case)
    
    ##cases
    if(K$case == 1){
      
      ##prior
      K$prior <- K.list$prior
      
      ##starting
      if(!"starting"%in%names(K.list))
        stop("error: 'starting' value is not specified for the ",param.name," parameter is not in the var.update.control list\n")
      
      if(length(K.list$starting) != 1)
        stop("error: 'starting' value is the wrong length for diag(",param.name,"), because this is a multiple response model and you specified a single prior it seems you want all diagonal elements of ",param.name," to share a common parameter therefore only specify one 'starting' value")
      
      K$starting <- K.list$starting
      
      ##tuning
      if(K$prior$dist != "FIXED"){ ##disregard if fixed
        if(!"tuning"%in%names(K.list))
          stop("error: 'tuning' value is not specified for the ",param.name," parameter in the var.update.control list\n")
        
        if(length(K.list$tuning) != 1)
                  stop("error: 'tuning' value is the wrong length for diag(",param.name,"), because this is a multiple response model and you specified a single prior it seems you want all diagonal elements of ",param.name," to share a common parameter therefore only specify one 'tuning' value")
        
        K$tuning <- K.list$tuning
        
        ##sample order
        if(!"sample.order"%in%names(K.list)){
          if(verbose){
            cat("note: 'sample.order' is not specified for the ",param.name," parameter in the var.update.control list; therefore, it is set to zero\n\n")
          }
          K$sample.order <- 0
        }else if(length(K.list$sample.order) != 1){
          stop("error: 'sample.order' should be of length 1 or m")
        }else{
          
          if(any(K.list$sample.order < 0))
            stop("error: for the ",param.name," parameter in the var.update.control list 'sample.order' must be >= 0\n")
          
          K$sample.order <- K.list$sample.order
        }
      }else{##fixed
        K$tuning <- -1
        K$sample.order <- -1
      }
            
      ##make sure storage mode is correct
      storage.mode(K$starting) <- "double"
      storage.mode(K$tuning) <- "double"
      storage.mode(K$sample.order) <- "integer"
      
    }else if(K$case == 2){
      
      ##prior
      K$prior <- K.list$prior
      
      ##starting
      if(!"starting"%in%names(K.list))
        stop("error: 'starting' value is not specified for the ",param.name," parameter is not in the var.update.control list\n")
      
      if(length(K.list$starting) == 1){
        ##expand to m
        K.list$starting <- rep(K.list$starting, m)
      }else if(length(K.list$starting) != m){
        stop("error: 'starting' value is the wrong length for diag(",param.name,")")
      }    
      
      K$starting <- K.list$starting
      
      ##tuning
      if(!"tuning"%in%names(K.list))
        stop("error: 'tuning' value is not specified for the ",param.name," parameter in the var.update.control list\n")
      
      if(length(K.list$tuning) == 1){
        ##expand to m
        K.list$tuning <- rep(K.list$tuning, m)
      }else if(length(K.list$tuning) != m){
        stop("error: 'tuning' value is the wrong length for diag(",param.name,")")
      }    
      
      K$tuning <- K.list$tuning
      
      ##sample order
      if(!"sample.order"%in%names(K.list)){
        if(verbose){
          cat("note: 'sample.order' is not specified for the ",param.name," parameter in the var.update.control list; therefore, it is set to zero\n\n")
        }
        K$sample.order <- rep(0, m)
      }else if(length(K.list$sample.order) == 1){
        if(any(K.list$sample.order < 0))
          stop("error: for the ",param.name," parameter in the var.update.control list 'sample.order' must be >= 0\n")
        
        ##expand to m
        K$sample.order <- rep(K.list$sample.order, m)
      }else if(length(K.list$sample.order) != m){
        stop("error: 'sample.order' should be of length 1 or m")
      }else{
        if(any(K.list$sample.order < 0))
          stop("error: for the ",param.name," parameter in the var.update.control list 'sample.order' must be >= 0\n")
        
        K$sample.order <- K.list$sample.order
      }

      
      ##make sure storage mode is correct
      storage.mode(K$starting) <- "double"
      storage.mode(K$tuning) <- "double"
      storage.mode(K$sample.order) <- "integer"
      
    }else if(K$case == 3){
      
      ##must be fixed or IWISH
      if(K.list$prior$dist=="FIXED"){
        
        if(!"starting"%in%names(K.list))
          stop("error: 'starting' value is not specified for the ",param.name," parameter is not in the var.update.control list\n")
        
        if(!is.matrix(K.list$starting))
          stop("error: 'starting' value must be a mxm matrix as ",param.name,"'s")
        
        if(dim(K.list$starting)[1] != dim(K.list$starting)[2] || dim(K.list$starting)[1] != m)
          stop("error: 'starting' value must be a mxm matrix as ",param.name,"'s")
        
        K$starting <- K.list$starting[lower.tri(K.list$starting, diag=TRUE)]

        storage.mode(K$starting) <- "double"
        K$sample.order <- -1
        storage.mode(K$sample.order) <- "integer"
        K$fixed <- as.integer(TRUE)      
        
      }else{##IWISH
        
        ##prior
        if(!K.list$prior$dist=="IWISH")
          stop("error: ",param.name,"'s prior distribution must be IWISH\n")
        
        ##check dim of S, already checked for square in prior formation
        if(dim(K.list$prior$params$S)[1] != m)
          stop("error: the IWISH prior on ",param.name," has the wrong dimension S hyperparameter")
        
        K$prior <- K.list$prior
        
        ##starting
        if(!"starting"%in%names(K.list))
          stop("error: 'starting' value is not specified for the ",param.name," parameter in the var.update.control list\n")
        
        if(!is.matrix(K.list$starting))
          stop("error: 'starting' value must be a mxm matrix as ",param.name,"'s")
        
        if(dim(K.list$starting)[1] != dim(K.list$starting)[2] || dim(K.list$starting)[1] != m)
          stop("error: 'starting' value must be a mxm matrix as ",param.name,"'s")
        
        K$starting <- K.list$starting[lower.tri(K.list$starting, diag=TRUE)] 
        
        ##tuning
        if(!"tuning"%in%names(K.list))
          stop("error: 'tuning' value is not specified for the ",param.name," parameter in the var.update.control list\n")
        
        if(!is.matrix(K.list$tuning))
          stop("error: 'tuning' value must be a matrix as IWISH was specified for ",param.name,"'s prior")
        
        if(dim(K.list$tuning)[1] != dim(K.list$tuning)[2] || dim(K.list$tuning)[1] != (m^2-m)/2+m)
          stop("error: 'tuning' value must be the upper Cholesky of a matrix with dimension (m^2-m)/2+m as IWISH was specified for ",param.name,"'s prior")
        
        K$tuning <- t(K.list$tuning)##I want lower Cholesky
                
        ##sample order
        if(!"sample.order"%in%names(K.list)){
          if(verbose){
            cat("note: 'sample.order' is not specified for the ",param.name," parameter in the var.update.control list; therefore, it is set to zero\n\n")
          }
          K$sample.order <- 0
        }else if(length(K.list$sample.order) != 1){
          stop("error: 'sample.order' should be of length 1")
        }else{
          if(K.list$sample.order < 0)
            stop("error: for the ",param.name," parameter in the var.update.control list 'sample.order' must be >= 0\n")
         
          K$sample.order <- K.list$sample.order
        }
        
        ##make sure storage mode is correct
        storage.mode(K$starting) <- "double"
        storage.mode(K$tuning) <- "double"
        storage.mode(K$sample.order) <- "integer"
        K$fixed <- as.integer(FALSE)
      }
      
    }else{
      stop("error: invalid case number in ",param.name," parameter var.update.control list\n")
    }
    ##return
    K
    
  }## end mkKPsi function
   
  ####################################################
  ##make the variance parameter list
  ####################################################
  ##
  ##mk K
  ##
  K <- mkPar("K", var.update.control$K, m)
  var <- list("K"=K)
  
  ##
  ##mk Psi
  ##
  if(!no.Psi){
    Psi <- mkPar("Psi", var.update.control$Psi, m)
    var$Psi <- Psi
  }
  var$no.Psi <- as.integer(no.Psi)
  
  ##
  ##mk phi
  ##
  phi <- list()
  
  ##prior
  if(!"prior"%in%names(var.update.control$phi))
    stop("error: 'prior' value is not specified for the phi parameter in the var.update.control list\n")
  
  if(class(var.update.control$phi$prior) == "prior"){##for m >= 1 with single prior assume separable
    
    var.update.control$phi$case <- 1
    
  }else if(class(var.update.control$phi$prior) == "list"){
    
    if(length(var.update.control$phi$prior) != m)
      stop("error: it seems you are trying to specify a non-separable model but the number of priors on the phi vector != m\n")      
    
    var.update.control$phi$case <- 2
    
  }else{
    stop("error: 'prior' value is misspecified for the phi parameter in the var.update.control list\n")
  }
  
  phi$case <- var.update.control$phi$case
  phi$prior <- var.update.control$phi$prior
  
  ##starting
  if(!"starting"%in%names(var.update.control$phi))
    stop("error: 'starting' value is not specified for the phi parameter is not in the var.update.control list\n")
  
  if(phi$case == 1){
    if(length(var.update.control$phi$starting) != 1)
      stop("error: 'starting' value is the wrong length for phi in the var.update.control list, should be length 1\n")
    
  }else{ ##case is 2
    if(length(var.update.control$phi$starting) == 1)
      var.update.control$phi$starting <- rep(var.update.control$phi$starting, m)
    
    if(length(var.update.control$phi$starting) != m)
      stop("error: 'starting' value is the wrong length for phi in the var.update.control list, should be length m\n")
  }
  
  phi$starting <- var.update.control$phi$starting
  
  ##tuning
  if(phi$case == 1){
    if(phi$prior$dist != "FIXED"){
      
      if(!"tuning"%in%names(var.update.control$phi))
        stop("error: 'tuning' value is not specified for the phi parameter is not in the var.update.control list\n")
      
      if(length(var.update.control$phi$tuning) != 1)
        stop("error: 'tuning' value is the wrong length for phi in the var.update.control list, should be length 1\n")
      
    }else{## fixed
      var.update.control$phi$tuning <- -1
    }
    
  }else{ ##case is 2
    
    if(!"tuning"%in%names(var.update.control$phi))
      stop("error: 'tuning' value is not specified for the phi parameter is not in the var.update.control list\n")
    
    if(length(var.update.control$phi$tuning) == 1)      ##expand
      var.update.control$phi$tuning <- rep(var.update.control$phi$tuning, m)
    
    if(length(var.update.control$phi$tuning) != m)
      stop("error: 'tuning' value is the wrong length for phi in the var.update.control list, should be length m\n")
  }
  
  phi$tuning <- var.update.control$phi$tuning
  
  ##sample order
  if(!"sample.order"%in%names(var.update.control$phi)){
    if(verbose){
      cat("note: 'sample.order' is not specified for the phi parameter in the var.update.control list; therefore, it is set to zero\n\n")
    }
    var.update.control$phi$sample.order <- 0
  }
  
  if(phi$case == 1){
    if(length(var.update.control$phi$sample.order) != 1)
      stop("error: 'sample.order' value is the wrong length for phi in the var.update.control list, should be length 1\n")
    
  }else{ ##case is 2
    
    if(length(var.update.control$phi$sample.order) == 1)      ##expand
      var.update.control$phi$sample.order <- rep(var.update.control$phi$sample.order, m)
    
    if(length(var.update.control$phi$sample.order) != m)
      stop("error: 'sample.order' value is the wrong length for phi in the var.update.control list, should be length m\n")
  }
  
  phi$sample.order <- var.update.control$phi$sample.order

  if(any(phi$sample.order < 0))
    stop("error: for the phi parameter in the var.update.control list 'sample.order' must be >= 0\n")
  
  ##make sure storage mode is correct
  storage.mode(phi$case) <- "integer"  
  storage.mode(phi$starting) <- "double"
  storage.mode(phi$tuning) <- "double"
  storage.mode(phi$sample.order) <- "integer"
  
  var$phi <- phi
  
  ##
  ##mk nu
  ##
  if(cov.model=="matern"){
    nu <- list()
    
    if(class(var.update.control$nu$prior) == "prior"){##for m >= 1 with single prior assume separable
      
      var.update.control$nu$case <- 1
      
    }else if(class(var.update.control$nu$prior) == "list"){
      
      if(length(var.update.control$nu$prior) != m)
        stop("error: it seems you are trying to specify a non-separable model but the number of priors on the nu vector != m\n")      
      
      var.update.control$nu$case <- 2
      
    }else{
      stop("error: 'prior' value is misspecified for the nu parameter in the var.update.control list\n")
    }
        
    nu$case <- var.update.control$nu$case
    if(nu$case != phi$case)
      stop("error: specifications for nu and phi must match (i.e., they both must either describe a separable or non-separable model)\n")

    
    nu$prior <- var.update.control$nu$prior
    
    ##starting
    if(!"starting"%in%names(var.update.control$nu))
      stop("error: 'starting' value is not specified for the nu parameter is not in the var.update.control list\n")
    
    if(nu$case == 1){
      if(length(var.update.control$nu$starting) != 1)
        stop("error: 'starting' value is the wrong length for nu in the var.update.control list, should be length 1\n")
      
    }else{ ##case is 2
      if(length(var.update.control$nu$starting) == 1)
        var.update.control$nu$starting <- rep(var.update.control$nu$starting, m)
      
      if(length(var.update.control$nu$starting) != m)
        stop("error: 'starting' value is the wrong length for nu in the var.update.control list, should be length m\n")
    }
    
    nu$starting <- var.update.control$nu$starting
    
    ##tuning
    if(nu$case == 1){
      if(nu$prior$dist != "FIXED"){
        
        if(!"tuning"%in%names(var.update.control$nu))
          stop("error: 'tuning' value is not specified for the nu parameter is not in the var.update.control list\n")
        
        if(length(var.update.control$nu$tuning) != 1)
          stop("error: 'tuning' value is the wrong length for nu in the var.update.control list, should be length 1\n")
        
      }else{## fixed
        var.update.control$nu$tuning <- -1
      }
      
    }else{ ##case is 2
      
      if(!"tuning"%in%names(var.update.control$nu))
        stop("error: 'tuning' value is not specified for the nu parameter is not in the var.update.control list\n")
      
      if(length(var.update.control$nu$tuning) == 1)      ##expand
        var.update.control$nu$tuning <- rep(var.update.control$nu$tuning, m)
      
      if(length(var.update.control$nu$tuning) != m)
        stop("error: 'tuning' value is the wrong length for nu in the var.update.control list, should be length m\n")
    }
    
    nu$tuning <- var.update.control$nu$tuning
    
    ##sample order
    if(!"sample.order"%in%names(var.update.control$nu)){
      if(verbose){
        cat("note: 'sample.order' is not specified for the nu parameter in the var.update.control list; therefore, it is set to zero\n\n")
      }
      var.update.control$nu$sample.order <- 0
    }
    
    if(nu$case == 1){
      if(length(var.update.control$nu$sample.order) != 1)
        stop("error: 'sample.order' value is the wrong length for nu in the var.update.control list, should be length 1\n")
      
    }else{ ##case is 2
      
      if(length(var.update.control$nu$sample.order) == 1)      ##expand
        var.update.control$nu$sample.order <- rep(var.update.control$nu$sample.order, m)
      
      if(length(var.update.control$nu$sample.order) != m)
        stop("error: 'sample.order' value is the wrong length for nu in the var.update.control list, should be length m\n")
    }
    
    nu$sample.order <- var.update.control$nu$sample.order
    
    if(any(nu$sample.order < 0))
      stop("error: for the nu parameter in the var.update.control list 'sample.order' must be >= 0\n")
    
    ##make sure storage mode is correct
    storage.mode(nu$case) <- "integer"  
    storage.mode(nu$starting) <- "double"
    storage.mode(nu$tuning) <- "double"
    storage.mode(nu$sample.order) <- "integer"
    
    var$nu <- nu
    
  }

  ####################################################
  ##make the beta parameter list
  ####################################################
  
  if(missing(beta.update.control)){stop("error: beta.update.control must be specified")}
  if(!is.list(beta.update.control)){stop("error: beta.update.control must be list")}

  ##check
  beta.list <- beta.update.control##less typing
  beta.list.names <- names(beta.list)
  beta <- list()

  ##prior
  if(!"prior"%in%beta.list.names){
    if(verbose){
      cat("note: 'prior' is not specified for the Beta parameters in the beta.update.control list; therefore, is set to flat\n\n")
    }
    beta$prior <- prior(dist="FLAT")
  }else{    
    if(class(beta.list$prior) != "prior" || !beta.list$prior$dist %in% c("NORMAL", "FLAT"))##someday allow fixed
      stop("error: for the Beta parameters the prior distribution must be NORMAL or FLAT\n")
    beta$prior <- beta.list$prior
  }

  
  if(beta$prior == "NORMAL" && length(beta$prior$params$mu) != n.beta)
    stop("error: for the Beta parameters the prior distribution is NORMAL but the hyperparameter dimensions do not equal the number of parameters\n")
  

  ##update method
  if(!"update"%in%beta.list.names){
    if(verbose){
      cat("note: 'update' method is not specified for the Beta parameters in the beta.update.control list; therefore, Gibbs will be used\n\n")
    }
    beta$update <- "GIBBS"
  }else if(!beta.list$update %in% c("GIBBS", "MH")){
    stop("error: for the Beta parameters the 'update' method must be GIBBS or MH, in the beta.update.control list\n")
  }else{
    beta$update <- beta.list$update
  }

  ##starting
  if(!"starting"%in%beta.list.names){
    if(verbose){
      cat("note: 'starting' is not specified for the Beta parameters in the beta.update.control list; therefore, starting values will be taken from lm(Y~X)\n\n")
    }
    beta$starting <- unname(coefficients(lm(Y ~ X-1)))
  }else if(length(beta.list$starting) == 1){
    ##expand to n.beta
    beta$starting <- rep(beta.list$starting, n.beta)  
  }else if(length(beta.list$starting) != n.beta){
    stop("error: 'starting' value is the wrong length for Beta, in the beta.update.control list")
  }else{
   beta$starting <- beta.list$starting
  }
  
  ##tuning
  if(beta$update == "MH"){
    if(!"tuning"%in%beta.list.names){
      if(verbose)
        cat("note: Beta tuning matrix 'tuning' is missing, will use the default tuning matrix chol(vcov(lm(Y~X)))\n\n")
      
      beta.list$tuning <- chol(vcov(lm(Y ~ X-1)))
    }

    ##print(dim(beta.list$tuning))

    if(n.beta > 1 && !is.matrix(beta.list$tuning))
      stop("error: 'tuning' value must be a matrix when the number of Beta > 1, in the beta.update.control list")

    if(n.beta == 1)
      beta.list$tuning <- as.matrix(beta.list$tuning[1])
    
    if(dim(beta.list$tuning)[1] != dim(beta.list$tuning)[2] || dim(beta.list$tuning)[1] != n.beta)
      stop("error: 'tuning' value must be the upper Cholesky of a matrix with dimension equal to the number of Beta parameters, in the beta.update.control list")
    
    beta$tuning <- t(beta.list$tuning)##I want lower Cholesky
  }
  
  ##sample order
  if(!"sample.order"%in%beta.list.names){
    if(verbose){
      cat("note: 'sample.order' is not specified for the Beta parameter in the beta.update.control list; therefore, it is set to zero\n\n")
    }
    beta$sample.order <- 0
  }else if(length(beta.list$sample.order) != 1){
    stop("error: 'sample.order' should be of length 1, in the beta.update.control list")
  }else{
    beta$sample.order <- beta.list$sample.order
  }

  if(beta$sample.order < 0)
    stop("error: for the Beta parameter in the var.update.control list 'sample.order' must be >= 0\n")
  
  ##make sure storage mode is correct
  storage.mode(beta$starting) <- "double"
  if(beta$update == "MH")
    storage.mode(beta$tuning) <- "double"
  storage.mode(beta$sample.order) <- "integer"
  beta$n.beta <- as.integer(n.beta)

  ####################################################
  ##Check support on UNIF prior
  ####################################################

  ##This is done in prior.R for now since any parameter that can receive a uniform prior is strictly greater than zero
  ##i.e., phi, K, Psi, nu. This will change later.
  
  
  #####################################################
  ##off it goes
  ####################################################

  args <- list("n.samples"=n.samples, "linpack"=linpack, "sp.effects"=sp.effects, "DIC"=DIC, "DIC.start"=DIC.start,
               "Y"=Y, "X"=X, "m"=m, "D"=D, "var.control"=var, "beta.control"=beta, "cov.model"=cov.model,
               "verbose"=as.integer(verbose))

  out <- .Call("ggtSp", args)
  rownames(out$p.samples)[(nrow(out$p.samples)-n.beta+1):nrow(out$p.samples)] <- x.names
  out$p.samples <- t(out$p.samples)
  
  ##if K$case or Psi$case is 3 then the corresponding p.samples are the t(chol(K)) or t(chol(Psi)), respectively.
  ##So return A%*%t(A) if needed for each.
  AtA <- function(x, m){
    A <- matrix(0, m, m)
    A[lower.tri(A, diag=TRUE)] <- x
    (A%*%t(A))[lower.tri(A, diag=TRUE)]
  }

  ##get the names right too (i.e., column major lower triangle)
  if(K$case == 3){
    K.names <- paste("K_",1:((m^2-m)/2+m),sep="")
    colnames(out$p.samples)[colnames(out$p.samples)%in%"K"] <- K.names
    out$p.samples[,K.names] <- t(apply(out$p.samples[,K.names], 1, AtA, m))
  }
  
  if(!no.Psi){
    if(Psi$case == 3){
      Psi.names <- paste("Psi_",1:((m^2-m)/2+m),sep="")
      colnames(out$p.samples)[colnames(out$p.samples)%in%"Psi"] <- Psi.names
      out$p.samples[,Psi.names] <- t(apply(out$p.samples[,Psi.names], 1, AtA, m))
    }
  }


  ##drop the first row of p.samples (it's the starting values)
  out$p.samples <- mcmc(out$p.samples[-1,])
  out$X <- X
  out$Y <- Y
  out$m <- m
  out$K.case <- K$case
  out$no.Psi <- no.Psi
  if(!no.Psi)
    out$Psi.case <- Psi$case
  out$phi.case <- phi$case
  out$cov.model <- cov.model
  out$coords <- coords
  class(out) <- "spGGT"
  out
}

