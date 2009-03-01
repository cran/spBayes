\name{mvCovInvLogDet}
\alias{mvCovInvLogDet}
\title{Utility function for constructing multivariate predictive process covariance matrices}

\description{
Utility function for constructing multivariate predictive process covariance matrices
}

\usage{
   mvCovInvLogDet(coords, knots, cov.model, V,
                  Psi, theta, modified.pp=TRUE, SWM=TRUE, ...)
}

\arguments{
  \item{coords}{a \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{knots}{a \eqn{m \times 2}{m x 2} matrix of knot
    coordinates. }
  \item{V}{spatial cross-covariance matrix. }  
  \item{Psi}{residual cross-covariance matrix.} 
  \item{theta}{if \code{cov.model} is \code{matern} then \code{theta} is a
  vector of the multivariate spatial ranges followed by the associated smoothness parameters; otherwise,
  \code{theta} is the spatial range vector.  } 
  \item{modified.pp}{a logical value indicating if the \emph{modified
    predictive process} should be used. }
  \item{cov.model}{a quoted key word that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }
  \item{SWM}{ a logical value indicating if the
  Sherman-Woodbury-Morrison equation should be used for computing the
  inverse (this is ignored if \code{modified.pp} is \code{FALSE}).}  
  \item{...}{currently no additional arguments.  }
}

\value{
  Returns the covariance matrix, log-determinant, and inverse.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
set.seed(1)

n <- 100 ##number of coords
m <- 25  ##number of knots
q <- 5   ##number of response variables

coords <- cbind(runif(n),runif(n))
knots <- cbind(runif(m), runif(m))

theta <- c(rep(6,q),rep(0.5,q)) #phi and nu

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- rnorm(q*(q-1)/2+q)

V <- A\%*\%t(A)

P <- matrix(0,q,q)
P[lower.tri(A,TRUE)] <- rnorm(q*(q-1)/2+q)

Psi <- P\%*\%t(P)

tol <- 1.0e-8


##
##non bias-adjusted predictive process
##
c1 <- mvCovInvLogDet(coords=coords, knots=knots, cov.model="exponential",
                   V=V, Psi=Psi, theta=theta,
                   modified.pp=FALSE, SWM=FALSE)

if(max(abs(chol2inv(chol(c1$C)) - c1$C.inv)) > tol) stop("test-1 failed")
if(max(abs(2*sum(log(diag(chol(c1$C)))) - c1$log.det)) > tol)
stop("test-1 failed")


c2 <- mvCovInvLogDet(coords=coords, knots=knots, cov.model="matern",
                   V=V, Psi=Psi, theta=theta,
                   modified.pp=FALSE, SWM=FALSE)

if(max(abs(chol2inv(chol(c2$C)) - c2$C.inv)) > tol) stop("test-2 failed")
if(max(abs(2*sum(log(diag(chol(c2$C)))) - c2$log.det)) > tol)
stop("test-2 failed")

##
##bias-adjusted predictive process
##
c3 <- mvCovInvLogDet(coords=coords, knots=knots, cov.model="exponential",
                   V=V, Psi=Psi, theta=theta,
                   modified.pp=TRUE, SWM=FALSE)

if(max(abs(chol2inv(chol(c3$C)) - c3$C.inv)) > tol) stop("test-3 failed")
if(max(abs(2*sum(log(diag(chol(c3$C)))) - c3$log.det)) > tol)
stop("test-3 failed")


c4 <- mvCovInvLogDet(coords=coords, knots=knots, cov.model="matern",
                   V=V, Psi=Psi, theta=theta,
                   modified.pp=TRUE, SWM=FALSE)

if(max(abs(chol2inv(chol(c4$C)) - c4$C.inv)) > tol) stop("test-4 failed")
if(max(abs(2*sum(log(diag(chol(c4$C)))) - c4$log.det)) > tol)
stop("test-4 failed")

##
##non bias-adjusted predictive process using SWM
##
c5 <- mvCovInvLogDet(coords=coords, knots=knots, cov.model="exponential",
                   V=V, Psi=Psi, theta=theta,
                   modified.pp=FALSE, SWM=TRUE)

if(max(abs(chol2inv(chol(c5$C)) - c5$C.inv)) > tol) stop("test-5 failed")
if(max(abs(2*sum(log(diag(chol(c5$C)))) - c5$log.det)) > tol)
stop("test-5 failed")

c6 <- mvCovInvLogDet(coords=coords, knots=knots, cov.model="matern",
                   V=V, Psi=Psi, theta=theta,
                   modified.pp=FALSE, SWM=TRUE)

if(max(abs(chol2inv(chol(c6$C)) - c6$C.inv)) > tol) stop("test-6 failed")
if(max(abs(2*sum(log(diag(chol(c6$C)))) - c6$log.det)) > tol)
stop("test-6 failed")

##
##bias-adjusted predictive process using SWM
##
c7 <- mvCovInvLogDet(coords=coords, knots=knots, cov.model="exponential",
                   V=V, Psi=Psi, theta=theta,
                   modified.pp=TRUE, SWM=TRUE)

if(max(abs(chol2inv(chol(c7$C)) - c7$C.inv)) > tol) stop("test-7 failed")
if(max(abs(2*sum(log(diag(chol(c7$C)))) - c7$log.det)) > tol)
stop("test-7 failed")

c8 <- mvCovInvLogDet(coords=coords, knots=knots, cov.model="matern",
                   V=V, Psi=Psi, theta=theta,
                   modified.pp=TRUE, SWM=TRUE)

if(max(abs(chol2inv(chol(c8$C)) - c8$C.inv)) > tol) stop("test-8 failed")
if(max(abs(2*sum(log(diag(chol(c8$C)))) - c8$log.det)) > tol)
stop("test-8 failed")

}

\keyword{misc}
