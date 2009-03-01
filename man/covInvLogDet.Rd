\name{covInvLogDet}
\alias{covInvLogDet}
\title{Utility function for constructing predictive process covariance matrices}

\description{
Utility function for constructing predictive process covariance matrices
}

\usage{
   covInvLogDet(coords, knots, cov.model, sigma.sq,
                  tau.sq, theta, modified.pp=TRUE, SWM=TRUE, ...)
}

\arguments{
  \item{coords}{a \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{knots}{a \eqn{m \times 2}{m x 2} matrix of knot
    coordinates. }
  \item{sigma.sq}{partial sill value. }  
  \item{tau.sq}{nugget value.} 
  \item{theta}{if \code{cov.model} is \code{matern} then \code{theta} is a
  vector of the spatial range and smoothness parameters; otherwise,
  \code{theta} is the spatial range.  } 
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

n <- 100
coords <- cbind(runif(n),runif(n))

m <- 50
knots <- cbind(runif(m), runif(m))

theta <- c(6, 2)
sigma.sq <- 1
tau.sq <- 1
tol <- 1.0e-5

##
##non bias-adjusted predictive process
##
c1 <- covInvLogDet(coords=coords, knots=knots, cov.model="exponential",
                   sigma.sq=sigma.sq, tau.sq=tau.sq, theta=theta,
                   modified.pp=FALSE, SWM=FALSE)

if(max(abs(chol2inv(chol(c1$C)) - c1$C.inv)) > tol) stop("test-1 failed")
if(max(abs(2*sum(log(diag(chol(c1$C)))) - c1$log.det)) > tol)
stop("test-1 failed")


c2 <- covInvLogDet(coords=coords, knots=knots, cov.model="matern",
                   sigma.sq=sigma.sq, tau.sq=tau.sq, theta=theta,
                   modified.pp=FALSE, SWM=FALSE)

if(max(abs(chol2inv(chol(c2$C)) - c2$C.inv)) > tol) stop("test-2 failed")
if(max(abs(2*sum(log(diag(chol(c2$C)))) - c2$log.det)) > tol)
stop("test-2 failed")

##
##bias-adjusted predictive process
##
c3 <- covInvLogDet(coords=coords, knots=knots, cov.model="exponential",
                   sigma.sq=sigma.sq, tau.sq=tau.sq, theta=theta,
                   modified.pp=TRUE, SWM=FALSE)

if(max(abs(chol2inv(chol(c3$C)) - c3$C.inv)) > tol) stop("test-3 failed")
if(max(abs(2*sum(log(diag(chol(c3$C)))) - c3$log.det)) > tol)
stop("test-3 failed")


c4 <- covInvLogDet(coords=coords, knots=knots, cov.model="matern",
                   sigma.sq=sigma.sq, tau.sq=tau.sq, theta=theta,
                   modified.pp=TRUE, SWM=FALSE)

if(max(abs(chol2inv(chol(c4$C)) - c4$C.inv)) > tol) stop("test-4 failed")
if(max(abs(2*sum(log(diag(chol(c4$C)))) - c4$log.det)) > tol)
stop("test-4 failed")

##
##non bias-adjusted predictive process using SWM
##
c5 <- covInvLogDet(coords=coords, knots=knots, cov.model="exponential",
                   sigma.sq=sigma.sq, tau.sq=tau.sq, theta=theta,
                   modified.pp=FALSE, SWM=TRUE)

if(max(abs(chol2inv(chol(c5$C)) - c5$C.inv)) > tol) stop("test-5 failed")
if(max(abs(2*sum(log(diag(chol(c5$C)))) - c5$log.det)) > tol)
stop("test-5 failed")

c6 <- covInvLogDet(coords=coords, knots=knots, cov.model="matern",
                   sigma.sq=sigma.sq, tau.sq=tau.sq, theta=theta,
                   modified.pp=FALSE, SWM=TRUE)

if(max(abs(chol2inv(chol(c6$C)) - c6$C.inv)) > tol) stop("test-6 failed")
if(max(abs(2*sum(log(diag(chol(c6$C)))) - c6$log.det)) > tol)
stop("test-6 failed")

##
##Check SWM
##
##c.str <- sigma.sq*exp(-theta[1]*iDist(knots))
##ct <- sigma.sq*exp(-theta[1]*iDist(coords,knots))
##C.inv.swm <- diag(1/tau.sq,n)-diag(1/tau.sq,n)/%*/%ct/%*/%
##  chol2inv(chol(c.str+t(ct)/%*/%diag(1/tau.sq,n)/%*/%ct))/%*/%
##  t(ct)/%*/%diag(1/tau.sq,n)
##if(max(abs(C.inv.swm - c5$C.inv)) > tol) stop("test-5a failed")

##
##bias-adjusted predictive process using SWM
##
c7 <- covInvLogDet(coords=coords, knots=knots, cov.model="exponential",
                   sigma.sq=sigma.sq, tau.sq=tau.sq, theta=theta,
                   modified.pp=TRUE, SWM=TRUE)

if(max(abs(chol2inv(chol(c7$C)) - c7$C.inv)) > tol) stop("test-7 failed")
if(max(abs(2*sum(log(diag(chol(c7$C)))) - c7$log.det)) > tol)
stop("test-7 failed")

c8 <- covInvLogDet(coords=coords, knots=knots, cov.model="matern",
                   sigma.sq=sigma.sq, tau.sq=tau.sq, theta=theta,
                   modified.pp=TRUE, SWM=TRUE)

if(max(abs(chol2inv(chol(c8$C)) - c8$C.inv)) > tol) stop("test-8 failed")
if(max(abs(2*sum(log(diag(chol(c8$C)))) - c8$log.det)) > tol)
stop("test-8 failed")

}

\keyword{misc}
