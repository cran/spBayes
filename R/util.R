
#parse formula and return a list that contains the model response
#matrix as element one, and the model matrix as element two
#lifted from MCMCpack (thank you)

"parse.formula" <-  function(formula, data, intercept=TRUE, justX=FALSE){
    
    # extract Y, X, and variable names for model formula and frame
    mt <- terms(formula, data=data)
    if(missing(data)) data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$intercept <- mf$justX <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    if (!intercept){
      attributes(mt)$intercept <- 0
    }

    # null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    X <- as.matrix(X)         # X matrix
    xvars <- dimnames(X)[[2]] # X variable names
    xobs  <- dimnames(X)[[1]] # X observation names
    if (justX){
      Y <- NULL
    }
    else {
      Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
    }
    return(list(Y, X, xvars, xobs))
  }




"print.model.matrix" <- function(x){
  cat("Regressors (i.e., Beta):\n")
  cat("\t", unlist(dimnames(x)[2]) ,"\n\n")
  if(!is.null(attributes(x)$contrasts)){
    cat("Factor regressors and associated contrast:\n")
    contrastList <- attributes(x)$contrasts
    for(i in 1:length(contrastList))
      cat("\t",names(attributes(x)$contrasts)[i], "\t", contrastList[[i]], "\n")              
  } 
}

#trim whitespace from the ends of a string
"trim" <- function(x){
  sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", x))
}
