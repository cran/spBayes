prior <- function(dist, ...){

  dist.params <- list(...)
  
  if(missing(dist))
    stop("dist type must be specified")
  #
  #distribution not found
  #
  if(!dist %in% c("IG", "ig", "UNIF", "unif", "LOGUNIF", "logunif", "HC", "hc", "NORMAL", "normal")){
    stop("dist type ",dist," not found\n")
  }

  #
  #
  #
  if(dist %in% c("NORMAL", "normal")){
    if("mu"%in%names(dist.params) && is.numeric(dist.params$mu)){
      mu <- as.matrix(dist.params$mu)
    }else{
      stop("mu not specified (or specified incorrectly) for NORMAL prior\n")
    }
    if("precision"%in%names(dist.params) && is.numeric(dist.params$precision)){
      precision <- as.matrix(dist.params$precision)
     }else{
       stop("precision not specified (or specified incorrectly) for NORMAL prior\n")
     }
    prior <- list(dist="NORMAL",params=list("mu"=mu,"precision"=precision))   
  }

  
  #
  #Inverse-gamma
  #
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

    prior <- list(dist="IG",params=list("shape"=shape,"scale"=scale))   
  }

  #
  #Uniform
  # 
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

    prior <- list(dist="UNIF", params=list("a"=a,"b"=b))
  }


  #
  #Log Uniform
  # 
  if(dist %in% c("LOGUNIF", "logunif")){

    
    if("a"%in%names(dist.params) && is.numeric(dist.params$a)){
      a <- dist.params$a
    }else{
      stop("'a' not specified (or specified incorrectly) for LOGUNIF prior\n")
    }
    if("b"%in%names(dist.params) && is.numeric(dist.params$b) && dist.params$b > dist.params$a){
      b <- dist.params$b
     }else{
       stop("'b' not specified (or specified incorrectly) for LOGUNIF prior (check if b > a\n")
     }

    prior <- list(dist="LOGUNIF", params=list("a"=a,"b"=b))
  }

    #
  #Log Half-Cauchy
  # 
  if(dist %in% c("HC", "hc")){
    
    if("a"%in%names(dist.params) && is.numeric(dist.params$a)){
      a <- dist.params$a
    }else{
      stop("'a' not specified (or specified incorrectly) for HC prior\n")
    }

    prior <- list(dist="HC", params=list("a"=a))
  }
  
  class(prior) <- "ggt.prior"
  prior
}
