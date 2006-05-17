"sp.DIC" <- function(ggt.sp.obj, start=1, end, thin=1, verbose=TRUE, ...){

  ####################################################
  ##Check for unused args
  ####################################################
  formal.args <- names(formals(sys.function(sys.parent())))
  elip.args <- names(list(...))
  for(i in elip.args){
    if(! i %in% formal.args)
      warning("'",i, "' is not an argument")
  }
  
  if(missing(ggt.sp.obj)){stop("error: DIC.sp expects an output object of class ggt.sp\n")}  
  if(class(ggt.sp.obj) != "ggt.sp"){stop("error:DIC.sp  requires an output object of class ggt.sp\n")}
  if(!is.logical(verbose)){stop("error: verbose must be of type logical\n")}
 
  ##subsample if specified
  samples <- as.matrix(ggt.sp.obj$p.samples)
  n.samples <- nrow(samples)

  if(missing(end))
    end <- n.samples
  
  if (!is.numeric(start) || start >= n.samples || start < 1)
    stop("error: invalid start")
  if (!is.numeric(end) || end > n.samples) 
    stop("error: invalid end")
  if (!is.numeric(thin) || thin >= n.samples) 
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
  coords <- as.matrix(ggt.sp.obj$coords)

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

  ##call
  out <- .Call("dicSp", "cov.model" = cov.model, "no.tausq"=no.tausq, "X"=X, "Y"=Y,
               "coords"=coords, "theta"=theta, "sigmasq"=sigmasq, "tausq"=tausq,
               "phi"=phi, "nu"=nu, "verbose"=as.integer(verbose))

  out1 <- list()
  out1$DIC <- out
  out1
}
