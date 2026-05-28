\name{iDist}
\alias{iDist}          
\title{Euclidean distance matrix}

\description{
  Computes the inter-site Euclidean distance matrix for one or two sets of points.
}

\usage{
 iDist(coords.1, coords.2, ...)
}

\arguments{
  \item{coords.1}{an \eqn{n\times p}{nxp} matrix with each row
  corresponding to a point in \eqn{p}{p} dimensional space. }
  \item{coords.2}{an \eqn{m\times p}{mxp} matrix with each row
  corresponding to a point in \eqn{p}{p} dimensional space. If this is
  missing then \code{coords.1} is used. }
  \item{...}{currently no additional arguments.}
}

\value{
  The \eqn{n\times n}{nxn} or \eqn{n\times m}{nxm} inter-site Euclidean distance matrix.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudiptob@biostat.umn.edu}, \cr
}

\examples{
\dontrun{
n <- 10
p1 <- cbind(runif(n),runif(n))

m <- 5
p2 <- cbind(runif(m),runif(m))

D <- iDist(p1, p2)
}
}

\keyword{utilities}
