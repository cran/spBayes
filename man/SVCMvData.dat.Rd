\name{SVCMvData.dat}
\docType{data}

\alias{SVCMvData.dat}

\title{Synthetic data from a space-varying coefficients model}

\description{
  Data simulated from a space-varying coefficients model.
}

\usage{
data(SVCMvData.dat)
}
 
\format{
  The data frame generated from the code in the example section below.
}

\examples{

\dontrun{
##The dataset was generated with the code below.

library(Matrix)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

set.seed(1)
n <- 200

coords <- cbind(runif(n,0,1), runif(n,0,1))
colnames(coords) <- c("x.coords","y.coords")

X <- as.matrix(cbind(1, rnorm(n), rnorm(n)))
colnames(X) <- c("intercept","a","b")

Z <- t(bdiag(as.list(as.data.frame(t(X)))))

beta <- c(1, 10, -10)
p <- length(beta)

q <- 3
A <- matrix(0, q, q)
A[lower.tri(A, T)] <- c(1, -1, 0, 1, 1, 0.1)
K <- A%*%t(A)
K
cov2cor(K)

phi <- c(3/0.75, 3/0.5, 3/0.5)

Psi <- diag(0,q)
C <- mkSpCov(coords, K, Psi, phi, cov.model="exponential")

tau.sq <- 0.1

w <- rmvn(1, rep(0,q*n), C)

y <- rnorm(n, as.vector(X\%*\%beta + Z\%*\%w), sqrt(tau.sq))

w.0 <- w[seq(1, length(w), by=q)]
w.a <- w[seq(2, length(w), by=q)]
w.b <- w[seq(3, length(w), by=q)]

SVCMvData <- data.frame(cbind(coords, y, X[,2:3], w.0, w.a, w.b))

}
}

\keyword{datasets}

