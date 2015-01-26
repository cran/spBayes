\name{spMvLM}
\alias{spMvLM}
\title{Function for fitting multivariate Bayesian spatial regression models}

\description{
  The function \code{spMvLM} fits Gaussian multivariate Bayesian
  spatial regression models. Given a set of knots, \code{spMvLM} will
  also fit a \emph{predictive process} model (see references below).
}

\usage{
spMvLM(formula, data = parent.frame(), coords, knots,
      starting, tuning, priors, cov.model,
      modified.pp = TRUE, amcmc, n.samples, 
      verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a list of \eqn{q}{q} symbolic regression model descriptions to be fit. See example below. }
  
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{spMvLM} is called.  }
  
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates in \eqn{R^2}{R^2} (e.g., easting and northing). }
  
  \item{knots}{either a \eqn{m \times 2}{m x 2} matrix of the \emph{predictive process} knot coordinates in \eqn{R^2}{R^2} (e.g., easting and northing) or a vector of length two or three with the first and second elements recording the  number of columns and rows in the desired knot grid. The third, optional, element sets the offset of the outermost knots from the extent of the \code{coords}. }

  \item{starting}{a list with tags corresponding to \code{beta}, \code{A}, \code{phi}, and \code{nu}. Depending on the specification of the non-spatial residual, tags are \code{L} or \code{Psi}
    for a block diagonal or diagonal covariance matrix, respectively.

The value portion of each tag is a vector
  that holds the parameter's starting values and are of length
  \eqn{p}{p} for \code{beta} (where \eqn{p}{p} is the total number of regression coefficients in the multivariate model), \eqn{\frac{q(q+1)}{2}}{q(q+1)/2} for
  \code{A} and \code{L}, and \eqn{q}{q} for \code{Psi}, \code{phi}, and \code{nu}.  Here,
    \code{A} and \code{L} hold the lower-triangle elements in column major ordering of the Cholesky square root
    of the spatial and non-spatial cross-covariance matrices,
    respectively. }
  
  \item{tuning}{a list with tags \code{A}, \code{phi}, and \code{nu}. Depending on the specification of
  the non-spatial residual, tags are \code{L} or \code{Psi}
  for a block diagonal or diagonal covariance matrix, respectively. The value portion of each tag defines the
    variance of the Metropolis sampler Normal proposal distribution.  For \code{A} and \code{L}
    the vectors are of length \eqn{\frac{q(q+1)}{2}}{q(q+1)/2} and \eqn{q}{q} for \code{Psi}, \code{phi}, and \code{nu}.}
  
  \item{priors}{a list with tags \code{beta.flat},
    \code{beta.norm}, \code{K.iw}, \code{Psi.iw}, \code{Psi.ig}, \code{phi.unif},
    and \code{nu.unif}. If the regression coefficients, i.e., \code{beta} vector, are
    assumed to follow a multivariate Normal distribution then pass the
    hyperparameters as a list of length two with the first
    and second elements corresponding to the mean vector and positive
    definite covariance matrix, respectively. If \code{beta} is assumed flat then no arguments are passed. The default is a flat prior.
    Use \code{Psi.iw} if the non-spatial residual covariance matrix is assumed block diagonal.  Otherwise if the non-spatial residual covariance matrix is assumed diagonal then each of the \eqn{q}{q} diagonal element are assumed to
    follow an inverse-Gamma in which case use \code{Psi.ig}. The hyperparameters of the inverse-Wishart, i.e., for
    cross-covariance matrices \eqn{AA'}{AA'} \code{K.iw}
    and \eqn{LL'}{LL'} \code{Psi.iw},  are
    passed as a list of length two, with the first and second elements corresponding
    to the \eqn{df}{df} and \eqn{q\times q}{qxq} \emph{scale} matrix,
    respectively. If \code{Psi.ig} is specified, the inverse-Gamma
    hyperparameters of the diagonal variance elements are pass using a
    list of length two with the first and second list elements consisting of
    vectors of the \eqn{q}{q} \emph{shape} and \emph{scale} hyperparameters,
    respectively. The hyperparameters of the Uniform \code{phi.unif},
    and \code{nu.unif} are also passed as a list of vectors with the first
    and second list elements corresponding to the lower and upper
    support, respectively. 
    }

  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }
  \item{modified.pp}{a logical value indicating if the \emph{modified
      predictive process} should be used (see references below for
    details). Note, if a predictive process model is not used (i.e., \code{knots} is not specified) then
    this argument is ignored. }
    \item{amcmc}{a list with tags \code{n.batch}, \code{batch.length}, and
  \code{accept.rate}. Specifying this argument invokes an adaptive MCMC
  sampler see Roberts and Rosenthal (2007) for an explanation.}
  \item{n.samples}{the number of MCMC iterations. This argument is
  ignored if \code{amcmc} is specified. }
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report Metropolis acceptance and MCMC progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{spMvLM}, which is a list with the following
  tags:
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix specified by
    \code{coords}.  }
  \item{knot.coords}{the \eqn{m \times 2}{m x 2} matrix as specified by \code{knots}.  }
  \item{p.theta.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  \item{acceptance}{the Metropolis sampling
    acceptance percent. Reported at \code{batch.length} or \code{n.report}
  intervals for \code{amcmc} specified and non-specified, respectively}
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\details{
Model parameters can be fixed at their \code{starting} values by setting their
\code{tuning} values to zero.

The \emph{no nugget} model is specified by removing \code{Psi} and \code{L} from the \code{starting} list.
}
  
\seealso{
  \code{\link{spLM}}
}

\references{
  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian Predictive Process Models for Large Spatial Datasets. \emph{Journal of the Royal Statistical Society Series B}, 70:825--848.
  
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.

  Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software}, 63:1--28. \url{http://www.jstatsoft.org/v63/i13}.
  
  Finley, A.O., H. Sang, S. Banerjee, and A.E. Gelfand. (2009) Improving the performance of predictive process modeling for large datasets. \emph{Computational Statistics and Data Analysis}, 53:2873--2884.
  
  Finley, A.O., S. Banerjee, A.R. Ek, and R.E. McRoberts. (2008) Bayesian multivariate process modeling for prediction of forest attributes. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, 13:60--83.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
\dontrun{
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p)))){stop("Dimension problem!")}
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

set.seed(1)

##Generate some data
n <- 25 ##number of locations
q <- 2 ##number of outcomes at each location
nltr <- q*(q+1)/2 ##number of triangular elements in the cross-covariance matrix

coords <- cbind(runif(n,0,1), runif(n,0,1))

##Parameters for the bivariate spatial random effects
theta <- rep(3/0.5,q)

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- c(1,-1,0.25)
K <- A\%*\%t(A)

Psi <- diag(0,q)

C <- mkSpCov(coords, K, Psi, theta, cov.model="exponential")

w <- rmvn(1, rep(0,nrow(C)), C)

w.1 <- w[seq(1,length(w),q)]
w.2 <- w[seq(2,length(w),q)]

##Covariate portion of the mean
x.1 <- cbind(1, rnorm(n))
x.2 <- cbind(1, rnorm(n))
x <- mkMvX(list(x.1, x.2))

B.1 <- c(1,-1)
B.2 <- c(-1,1)
B <- c(B.1, B.2)

Psi <- diag(c(0.1, 0.5))

y <- rnorm(n*q, x\%*\%B+w, diag(n)\%x\%Psi)

y.1 <- y[seq(1,length(y),q)]
y.2 <- y[seq(2,length(y),q)]

##Call spMvLM
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
n.samples <- 1000

starting <- list("phi"=rep(3/0.5,q), "A"=A.starting, "Psi"=rep(1,q))
tuning <- list("phi"=rep(1,q), "A"=rep(0.01,length(A.starting)), "Psi"=rep(0.01,q))
priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
               "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(c(2,2), c(0.1,0.1)))

m.1 <- spMvLM(list(y.1~x.1-1, y.2~x.2-1), 
               coords=coords, starting=starting, tuning=tuning, priors=priors,
               n.samples=n.samples, cov.model="exponential", n.report=100)

burn.in <- 0.75*n.samples

m.1 <- spRecover(m.1, start=burn.in)

round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)

m.1.w.hat <- summary(mcmc(t(m.1$p.w.recover.samples)))$quantiles[,c(3,1,5)]
m.1.w.1.hat <- m.1.w.hat[seq(1, nrow(m.1.w.hat), q),]
m.1.w.2.hat <- m.1.w.hat[seq(2, nrow(m.1.w.hat), q),]

par(mfrow=c(1,2))
plot(w.1, m.1.w.1.hat[,1], xlab="Observed w.1", ylab="Fitted w.1",
     xlim=range(w), ylim=range(m.1.w.hat), main="Spatial random effects w.1")
arrows(w.1, m.1.w.1.hat[,1], w.1, m.1.w.1.hat[,2], length=0.02, angle=90)
arrows(w.1, m.1.w.1.hat[,1], w.1, m.1.w.1.hat[,3], length=0.02, angle=90)
lines(range(w), range(w))

plot(w.2, m.1.w.2.hat[,1], xlab="Observed w.2", ylab="Fitted w.2",
     xlim=range(w), ylim=range(m.1.w.hat), main="Spatial random effects w.2")
arrows(w.2, m.1.w.2.hat[,1], w.2, m.1.w.2.hat[,2], length=0.02, angle=90)
arrows(w.2, m.1.w.2.hat[,1], w.2, m.1.w.2.hat[,3], length=0.02, angle=90)
lines(range(w), range(w))
}
}

\keyword{model}