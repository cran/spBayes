\name{rf.n1000.dat}
\docType{data}

\alias{rf.n1000.dat}

\title{Synthetic univariate data used for illustrations}

\description{
   The synthetic dataset was generated from a stationary and
   isotropic spatial process using an exponential correlation function (code below).
}

\usage{
data(rf.n1000.dat)
}

\format{
  The rf.n1000.dat files contain n rows (i.e., points) and 4
  columns. Column 1, Y, is the response variable. Column 2, w, is the
  vector of random spatial effects. Columns 4 and 5 are the x and y coordinates.
  }

\examples{
\dontrun{

##These synthetic data were generated with the following:

n <- 1000

coords <- cbind(runif(n, 1, 10), runif(n, 1, 10))

sigma.sq <- 10
tau.sq <- 1
phi <- 3/5

C <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- mvrnorm(1, rep(0,n), C)

Y <- w+rnorm(n, 0, sqrt(tau.sq))

}

}
  
\keyword{datasets}

 