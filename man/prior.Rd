\name{prior}
\alias{prior}

\title{Creates prior distribution definitions for use in ggt and ggt.sp}

\description{
  The function \code{prior} creates a valid prior definition for use in
  the \code{\link{ggt}} and \code{\link{ggt.sp}} functions.
}

\usage{
prior(dist, ...)
}

\arguments{
  \item{dist}{a quoted key word that identifies the prior distribution.  The
    choices are inverse gamma \code{"IG"}, uniform \code{"UNIF"},
    log-uniform, \code{"LOGUNIF"}, half-Cauchy \code{"HC"}, and normal \code{"NORMAL"}.  }
  \item{...}{the hyperparameters of the chosen prior distribution,
    passed as quoted key words with associated values. See details below.  }
}

\value{
  \item{prior}{a list object of class \code{ggt.prior} which is used for the value portion of the prior tag in
  the \code{\link{ggt}} and \code{\link{ggt.sp}} functions. }
}


\details{
  Up to a proportionality constant the possible priors are:
  \itemize{
    \item{Uniform: }{\eqn{1/(b-a), b > a}{1/(b-a)}, where \eqn{b >
	a}. Use key words \code{"a"}, and \code{"b"}. }
    
    \item{Log-uniform: }{\eqn{1/(b-a), b > a}{1/(b-a)}, where \eqn{b > a
	> 0}.  When the associated parameter is known to be strictly
      greater than zero, the log-uniform should be used to improve
      sampler efficiency.
      Use key words \code{"a"}, and \code{"b"}.  }
    
    \item{Inverse gamma: }{\eqn{x^{-(a+1)} exp(-b/x)}{x^{-(a+1)}
	exp(-b/x)}, where \eqn{a} is shape \eqn{> 0} and \eqn{b} is scale \eqn{> 0}.
      Use key words \code{"shape"}, and \code{"scale"}. }
    
    \item{Half-Cauchy:}{\eqn{(x^2 + a^2)^{-1}}{(x^2 + a^2)^{-1}}. Use key word \code{"a"}. }
    
    \item{Normal:}{\eqn{exp(-1/2 (\beta - mu)^t V^{-1} (\beta -
	mu))}{exp(-1/2 (B - mu)^t V^{-1} (B - mu))}, where \eqn{mu} is the
      mean vector of length \eqn{p} and \eqn{V^{-1}}{V^{-1}} is the
      \eqn{p \times p}{p x p} precision matrix.  Use key words
      \code{"mu"} and \code{"precision"}. }
  }
  The normal prior can only be used for the regressors specified in \code{\link{ggt}} and \code{\link{ggt.sp}}.
}


\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.

  Gelman, A. (2005). Prior distributions for variance parmaters in
  hierarchical models. 1(2):1-19, Bayesian Analysis.
  
  Further information on the package \pkg{spBayes} can be found
  at: \url{http://blue.fr.umn.edu/spatialBayes}.
}

\author{
  Andrew O. Finley \email{afinley@stat.umn.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
  Bradley P. Carlin \email{brad@biostat.umn.edu}
}

\seealso{\code{\link{ggt}},  \code{\link{ggt.sp}} }

\keyword{misc}
