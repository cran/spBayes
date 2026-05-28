\name{mkMvX}
\alias{mkMvX}          
\title{Make a multivariate design matrix}

\description{
  Given \eqn{q}{q} univariate design matrices, the function
  \code{mkMvX} creates a multivariate design matrix suitable for use in \code{\link{spPredict}}.
}

\usage{
  mkMvX(X)
}

\arguments{
  \item{X}{a list of \eqn{q}{q} univariate design matrices.  The
    matrices must have the same number of rows (i.e., observations) but
    may have different number of columns (i.e., regressors). }
}

\value{
  A multivariate design matrix suitable for use in \code{\link{spPredict}}.
}

\seealso{\code{\link{spPredict}}}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}.
}

\examples{
\dontrun{
##Define some univariate model design matrices
##with intercepts.
X.1 <- cbind(rep(1, 10), matrix(rnorm(50), nrow=10))
X.2 <- cbind(rep(1, 10), matrix(rnorm(20), nrow=10))
X.3 <- cbind(rep(1, 10), matrix(rnorm(30), nrow=10))

##Make a multivariate design matrix suitable
##for use in spPredict.
X.mv <- mkMvX(list(X.1, X.2, X.3))
}
}

\keyword{manip}