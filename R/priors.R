prior <- function(dist, ...){

  dist.params <- list(...)

  is.sq <- function(mtrx){
    if(dim(mtrx)[1] != dim(mtrx)[2])
      FALSE
    else
      TRUE
  }
  
  if(missing(dist))
    stop("dist type must be specified")
  ##
  ##distribution not found
  ##
  if(!dist %in% c("IG", "ig", "UNIF", "unif", "LOGUNIF", "logunif", "HC", "hc", "NORMAL", "normal", "IWISH", "iwish", "FIXED", "fixed", "flat", "FLAT")){
    stop("dist type ",dist," not found\n")
  }

  if(dist %in% c("LOGUNIF", "logunif")){
    stop("As of spBayes_0.0-6 the there is no log-uniform prior. Please replace with the uniform prior.")
  }

  
  ##
  ##Inverse-Wishart
  ##
  if(dist %in% c("IWISH", "iwish")){
    if("df"%in%names(dist.params) && is.numeric(dist.params$df)){
      df <- as.matrix(dist.params$df)
    }else{
      stop("df not specified (or specified incorrectly) for IWISH prior\n")
    }
    if("S"%in%names(dist.params) && is.matrix(dist.params$S) && is.sq(dist.params$S)){
      S <- dist.params$S
     }else{
       stop("S not specified (or specified incorrectly) for IWISH prior\n")
     }

    ##check if S is PD
    options(show.error.messages = FALSE)
    pd <- try(chol(S))
    options(show.error.messages = TRUE)
    if(class(pd) == "try-error")
      stop("S must be a positive definate matrix for the IWISH prior\n")
    
    ##storage mode check
    storage.mode(df) <- "double"
    storage.mode(S) <- "double"
    
    prior <- list(dist="IWISH",params=list("df"=df,"S"=S))   
  }

  
  ##
  ##Normal
  ##
  if(dist %in% c("NORMAL", "normal")){
    if("mu"%in%names(dist.params) && is.numeric(dist.params$mu)){
      mu <- as.matrix(dist.params$mu)
    }else{
      stop("mu not specified (or specified incorrectly) for NORMAL prior\n")
    }

    if(!"precision"%in%names(dist.params))
      stop("precision not specified for NORMAL prior\n")
    
    if(length(mu) == 1)##assume just one parameter so precision can be a scalar
      dist.params$precision <- as.matrix(dist.params$precision)
    
    if(is.matrix(dist.params$precision) && is.sq(dist.params$precision)){
      precision <- dist.params$precision
    }else{
      stop("precision specified incorrectly for NORMAL prior\n")
    }
    
    if(dim(precision)[1] != length(mu))
      stop("the dimension of the square precision matrix must equal the length of the mu vector for NORMAL prior\n")


    ##check if precision is PD
    options(show.error.messages = FALSE)
    pd <- try(chol(precision))
    options(show.error.messages = TRUE)
    if(class(pd) == "try-error")
      stop("precision must be a positive definate matrix for the NORMAL prior\n")

    
    ##storage mode check
    storage.mode(mu) <- "double"
    storage.mode(precision) <- "double"
    
    prior <- list(dist="NORMAL",params=list("mu"=mu,"precision"=precision))   
  }
  
  ##
  ##flat
  ##
  if(dist %in% c("FLAT", "flat")){
    prior <- list(dist="FLAT", params=-1)
  }
  
  ##
  ##Inverse-gamma
  ##
  if(dist %in% c("IG", "ig")){
    if("shape"%in%names(dist.params) && is.numeric(dist.params$shape) && dist.params$shape > 0){
      shape <- dist.params$shape
    }else{
      warning("shape not specified (or specified incorrectly) for IG prior, set to default 2\n")
      shape <- 2
    }
    if("scale"%in%names(dist.params) && is.numeric(dist.params$scale) && dist.params$scale > 0){
      scale <- dist.params$scale
     }else{
       warning("scale not specified (or specified incorrectly) for IG prior, set to default 0.001\n")
       scale <- 0.001
     }

    ##storage mode check
    storage.mode(shape) <- "double"
    storage.mode(scale) <- "double"
    
    prior <- list(dist="IG",params=list("shape"=shape,"scale"=scale))   
  }

  ##
  ##Uniform
  ## 
  if(dist %in% c("UNIF", "unif")){

    if("a"%in%names(dist.params) && is.numeric(dist.params$a)){
      a <- dist.params$a
    }else{
      stop("'a' not specified (or specified incorrectly) for UNIF prior\n")
    }
    if("b"%in%names(dist.params) && is.numeric(dist.params$b) && dist.params$b > dist.params$a){
      b <- dist.params$b
     }else{
       stop("'b' not specified (or specified incorrectly) for UNIF prior (check if b > a\n")
     }
   
    ##storage mode check
    storage.mode(a) <- "double"
    storage.mode(b) <- "double"

    ##for now any parameter that receives a UNIF in this model ggt.sp must have support > 0
    if(a <= 0)
      stop("For this parameter the UNIF prior support must be greater than zero (i.e., 'a' > 0)\n")
    
    prior <- list(dist="UNIF", params=list("a"=a,"b"=b))
  }


  ##
  ##Log Half-Cauchy
  ## 
  if(dist %in% c("HC", "hc")){
    
    if("a"%in%names(dist.params) && is.numeric(dist.params$a)){
      a <- dist.params$a
    }else{
      stop("'a' not specified (or specified incorrectly) for HC prior\n")
    }

    ##storage mode check
    storage.mode(a) <- "double"
    
    prior <- list(dist="HC", params=list("a"=a))
  }

  ##
  ##fixed
  ## 
  if(dist %in% c("FIXED", "fixed")){
    prior <- list(dist="FIXED", params=-1)
  }
  
  class(prior) <- "prior"
  prior
}
