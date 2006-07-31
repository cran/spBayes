\name{mk.mv.X}
\alias{mk.mv.X}          
\title{Make a multivariate design matrix}

\description{
  Given \eqn{m}{m} univariate design matrices, the function
  \code{mk.mv.X} creates a multivariate design matrix suitable for use in \code{\link{sp.predict}}.
}

\usage{
  mk.mv.X(X)
}

\arguments{
  \item{X}{a list of \eqn{m}{m} univariate design matrices.  The
    matrices must have the same number of rows (i.e., observations) but
    may have different number of columns (i.e., regressors). }
}

\value{
  A multivariate design matrix suitable for use in \code{\link{sp.predict}}.
}

\seealso{\code{\link{sp.predict}}}

\author{
  Andrew O. Finley \email{afinley@stat.umn.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
  Bradley P. Carlin \email{brad@biostat.umn.edu}.
}

\examples{
##Define some univariate model design matrices
##with intercepts.
X.1 <- cbind(rep(1, 10), matrix(rnorm(50), nrow=10))
X.2 <- cbind(rep(1, 10), matrix(rnorm(20), nrow=10))
X.3 <- cbind(rep(1, 10), matrix(rnorm(30), nrow=10))

##Make a multivariate design matrix suitable
##for use in sp.predict.
X.mv <- mk.mv.X(list(X.1, X.2, X.3))
}

\keyword{misc}
