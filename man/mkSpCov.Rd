\name{mkSpCov}
\alias{mkSpCov}
\title{Function for calculating univariate and multivariate covariance matrices}

\description{
  The function \code{mkSpCov} calculates a spatial covariance matrix
  given spatial locations and spatial covariance parameters. 
}

\usage{
mkSpCov(coords, K, Psi, theta, cov.model)
}

\arguments{

  \item{coords}{an \eqn{n \times 2}{nx2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  
  \item{K}{the \eqn{q \times q}{q x q} spatial cross-covariance
  matrix. For a univariate model this corresponds to the partial sill,
  \eqn{\sigma^2}{sigma^2}.}

\item{Psi}{the \eqn{q \times q}{qxq} non-spatial covariance
matrix. For a univariate model this corresponds to the nugget,
\eqn{\tau^2}{tau^2}.}

  \item{theta}{a vector of \eqn{q}{q} spatial decay parameters. If
    \code{cov.model} is \code{"matern"} then \code{theta} is a vector of
    length \eqn{2\times q}{2xq} with the spatial decay parameters in the first
    \eqn{q}{q} elements and the spatial smoothness parameters in the last
    \eqn{q}{q} elements.}
  
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }

  \item{...}{currently no additional arguments.  }
}

\value{
\item{C}{the \eqn{nq\times nq}{nqxnq} spatial covariance matrix. }
}

\details{

   Covariance functions return the covariance
  \eqn{C(h)} between a pair locations separated by distance \eqn{h}. The covariance function can be written as a product of a variance parameter \eqn{\sigma^2}{\sigma^2} and a positive definite \emph{correlation function} \eqn{\rho(h)}{\rho(h)}: \eqn{C(h) = \sigma^2 \rho(h)}{C(h) = \sigma^2 * \rho(h)}, see, e.g.,
  Banerjee et al. (2004) p. 27 for more details. The expressions of the correlations functions available in \pkg{spBayes} are given below. More will be added upon request.
  
  For all correlations functions, \eqn{\phi}{phi} is the spatial \emph{decay} parameter.
  Some of the correlation functions will have an extra parameter
  \eqn{\nu}{nu}, the \emph{smoothness} parameter.
  \eqn{K_\nu(x)}{K_\nu(x)} denotes the modified Bessel
  function of the third kind of order \eqn{\nu}{nu}. See
  documentation of the function \code{besselK} for further details.
  The following functions are valid for \eqn{\phi>0}{\phi >
    0} and \eqn{\nu>0}{\nu > 0}, unless stated otherwise.
  
  \bold{gaussian}\cr
  \deqn{\rho(h) = \exp[-(\phi h)^2]}{\rho(h) = exp(-(\phi*h)^2)} 

  \bold{exponential}\cr
  \deqn{\rho(h) = \exp(-\phi h)}{\rho(h) = exp(-\phi*h)}

  \bold{matern}\cr
  \deqn{\rho(h) =
    \frac{1}{2^{\nu-1}\Gamma(\nu)}(\phi h)^\nu
    K_{\nu}(\phi h)}{\rho(h) =
    (1/(2^(\nu-1) * \Gamma(\nu))) * ((\phi*h)^\nu) * K_{\nu}(\phi*h)} 

  \bold{spherical}\cr
  \deqn{\rho(h) = \left\{ \begin{array}{ll}
    1 - 1.5\phi h + 0.5(\phi h)^3
    \mbox{ , if $h$ < $\frac{1}{\phi}$} \cr
    0    \mbox{ ,  otherwise}
    \end{array} \right.}{\rho(h) = 
    1 - 1.5 * (\phi*h) + 0.5*(\phi*h)^3   if h < 1/\phi , 
    0   otherwise}
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
\dontrun{
##A bivariate spatial covariance matrix

n <- 2 ##number of locations
q <- 2 ##number of responses at each location
nltr <- q*(q+1)/2 ##number of triangular elements in the cross-covariance matrix

coords <- cbind(runif(n,0,1), runif(n,0,1))

##spatial decay parameters
theta <- rep(6,q)

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- rnorm(nltr, 5, 1)
K <- A\%*\%t(A)

Psi <- diag(1,q)

C <- mkSpCov(coords, K, Psi, theta, cov.model="exponential")
}
}

\keyword{utilities}