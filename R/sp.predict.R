"sp.predict" <- function(ggt.sp.obj, pred.joint = TRUE, pred.coords, pred.covars, use.covar.names = TRUE, start=1, end, thin=1, verbose=TRUE, ...){

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  if(missing(ggt.sp.obj)){stop("error: predict.sp expects an output object of class ggt.sp\n")}  
  if(class(ggt.sp.obj) != "ggt.sp"){stop("error: predict.sp requires an output object of class ggt.sp\n")}
  if(missing(pred.coords)){stop("error: predict.sp requires a data.frame or matrix of points to predict, pred.coords\n")}
  if(!any(is.data.frame(pred.coords), is.matrix(pred.coords))){stop("error: pred.coords must be a data.frame or matrix\n")}
  if(!any(ncol(pred.coords) == 2, all(c("X","Y") %in% colnames(pred.coords)))){stop("error: pred.coords must have two columns (assumed to be X,Y) or include column names 'X' and 'Y'\n")}

  if(!missing(pred.covars)){
    if(!any(is.data.frame(pred.covars), is.matrix(pred.covars))){stop("error: pred.covars must be a data.frame or matrix\n")}
    if(nrow(pred.coords) != nrow(pred.covars)){stop("error: number of rows in pred.covars must be the same as pred.covars, and coordinate rows should correspond to covariate rows, although this is obviously not checked\n")}
  }
  
  if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}
  if(!is.logical(use.covar.names)){stop("error: use.covar.names must be of type logical\n")}
  if(!is.logical(pred.joint)){stop("error: pred.joint must be of type logical\n")}

  
  ##subsample if specified
  samples <- as.matrix(ggt.sp.obj$p.samples)
  n.samples <- nrow(samples)

  if(missing(end))
    end <- n.samples
  
  if(!is.numeric(start) || start >= n.samples)
    stop("error: invalid start")
  if(!is.numeric(end) || end > n.samples) 
    stop("error: invalid end")
  if(!is.numeric(thin) || thin >= n.samples) 
    stop("error: invalid thin")

  ##note that the first row in ggt.sp.obj$p.samples is the starting values
  ##so if start is 1 make it 2, the offset is fixed on the C++ side
  if(start == 1)
    start <- 2
    
  samples <- samples[seq(start, end, by=as.integer(thin)),]
  
  ##get other stuff out of ggt.sp.obj
  cov.model <- ggt.sp.obj$cov.model
  no.tausq <- ggt.sp.obj$no.tausq
  X <- as.matrix(ggt.sp.obj$X)
  Y <- as.matrix(ggt.sp.obj$Y)
  obs.coords <- as.matrix(ggt.sp.obj$coords)
  
  ##get var component samples
  tausq <- as.numeric(0)
  if(!no.tausq)
    tausq <- as.matrix(samples[,"tausq"])
  
  sigmasq <- as.matrix(samples[,"sigmasq"])
  phi <- as.matrix(samples[,"phi"])
  
  nu <- as.numeric(0)
  if(cov.model == "matern")
    nu <- as.matrix(samples[,"nu"])
  
  ##as this is a ggt.sp obj all theta are the first ncol(X) in the
  ##sample matrix, also transpose theta to make BLAS access simple
  theta <- t(as.matrix(samples[,1:ncol(X)]))
  
  ##check and get the covariates associated with each prediction
  if(ncol(X) != 1 && missing(pred.covars)){stop("error: pred.covars needs to be specified, pred.covars should be a matrix or data frame with columns corresponding to the ",ncol(X)," fixed effects samples in the ggt.sp.obj\n")}

  pred.X <- NULL

  if(ncol(X) == 1 && missing(pred.covars)){#this assumes just an intercept

    pred.X <- as.matrix(rep(1, nrow(pred.coords)))
    
  }else if(!use.covar.names){
    
    if(ncol(pred.covars) != ncol(X))
      {stop("error: as use.covar.names is FALSE, the number of columns in the matrix or data frame specified by pred.covars must be equal to the number of covariates specified in ggt.sp.obj, which is ",ncol(X),"\n")}
    
    pred.X <- as.matrix(pred.covars)
    
  }else{##use.covar.names == TRUE
    
    if(any(is.null(colnames(X)), is.null(colnames(pred.covars))))
      stop("error: argument use.covar.names is TURE but pred.covars and/or ggt.sp.obj$p.samples does not have column names\n")

    ##if the '(Intercept)' is in X but not in pred.covars then add it
    if(all("(Intercept)" %in% colnames(X), !"(Intercept)" %in% colnames(pred.covars))){
      intercept <- as.matrix(rep(1.0, nrow(pred.coords)))
      colnames(intercept) <- "(Intercept)"
      pred.covars <- cbind(intercept, pred.covars)
      warning("warning: '(Intercept)' added to 'pred.covars'\n")
    }
          
    if(!all(colnames(X) %in% colnames(pred.covars)))
      stop("error: argument use.covar.names is TURE but the column name(s) ", colnames(X)[!colnames(X)%in%colnames(pred.covars)]," are missing from pred.covars\n")

    pred.X <- as.matrix(pred.covars[,colnames(X)])

  }
  ##just double check dims
  if(ncol(pred.X) != nrow(theta))
    stop("error: this should never happen, the number of prediction covariates != number of sampled theta\n")

  
  ##get distance from each prediction point to each observed point
  ##note, the pred.coords were checked above
  if(is.null(colnames(pred.coords))){##assumed 2 columns X,Y
    pred.X.coords <- matrix(pred.coords[,1], nrow(obs.coords), nrow(pred.coords), byrow=T)
    pred.Y.coords <- matrix(pred.coords[,2], nrow(obs.coords), nrow(pred.coords), byrow=T)
  }else{
    pred.X.coords <- matrix(pred.coords[,"X"], nrow(obs.coords), nrow(pred.coords), byrow=T)
    pred.Y.coords <- matrix(pred.coords[,"Y"], nrow(obs.coords), nrow(pred.coords), byrow=T)
  }

  gamma.D <- sqrt((pred.X.coords-obs.coords[,1])^2+(pred.Y.coords-obs.coords[,2])^2)
  gamma.D <- t(gamma.D) ##for the dsymm

  D <- as.matrix(dist(obs.coords))

  pred.D <- as.matrix(dist(pred.coords))
  
  #other stuff
  if(pred.joint && nrow(pred.coords)==1)
    pred.joint <- FALSE #helps on the C++ side


  ##call
  out <- .Call("predictSp", "cov.model"=cov.model, "no.tausq"=no.tausq, "X"=X, "Y"=Y,
               "coords"=obs.coords, "theta"=theta, "sigmasq"=sigmasq, "tausq"=tausq,
               "phi"=phi, "nu"=nu, "D"=D, "gamma.D"=gamma.D, "pred.D"=pred.D, "pred.X"=pred.X,
               "pred.joint"=as.integer(pred.joint), "verbose"=as.integer(verbose))
  
  out1 <- list()
  #out1$cov.model <- cov.model
  #out1$no.tausq <- as.logical(no.tausq)
  out1$obs.coords <- obs.coords
  out1$pred.coords <- pred.coords
  #out1$gamma.D <- gamma.D
  #out1$pred.D <- pred.D
  #out1$obs.D <- D
  out1$pp.samples <- out

  out1
}
 
