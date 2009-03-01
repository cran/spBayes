covGivens <- function(theta, lambda, p, strict=TRUE)
{

  if(strict){
    if(length(theta) != p*(p-1)/2 || any(theta > pi/2 && theta > -pi/2))
      {stop("error: 'theta' must be a vector of length p*(p-1)/2 and in the interval (pi/2, -pi/2)")}
    
    if(length(lambda) != p || any(lambda < 0))
      {stop("error: 'lambda' must be a positive vector of length p")}
  }else{
    if(length(theta) != p*(p-1)/2)
      {stop("error: 'theta' must be a vector of length p*(p-1)/2")}
    
    if(length(lambda) != p)
      {stop("error: 'lambda' must be a vector of length p")}
  }
    
  cov <- matrix(0, p, p)
  tmp1 <- cov
  P <- cov

  storage.mode(theta) <- "double"
  storage.mode(lambda) <- "double"
  storage.mode(tmp1) <- "double"
  storage.mode(P) <- "double"
  
  junk <- .Call("covGivens", as.integer(p), theta, lambda, cov, tmp1, P)

  list("C"=cov, "P"=P, "Lambda"=lambda)

}
