\name{spSVC}
\alias{spSVC}
\title{Function for fitting univariate Bayesian spatially-varying
  coefficient regression models}

\description{
  The function \code{spSVC} fits Gaussian univariate Bayesian spatially-varying
  coefficient regression models.
}

\usage{
spSVC(formula, data = parent.frame(), svc.cols=1, coords, 
      priors, starting, tuning, cov.model, center.scale=FALSE,
      amcmc, n.samples, n.omp.threads = 1,
      verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which
    \code{spSVC} is called. }
  
  \item{svc.cols}{a vector indicating which columns of the regression
    design matrix \eqn{X}{X} should be space-varying. \code{svc.cols} can
    be an integer vector with values indicating \eqn{X}{X} columns or a
    character vector with values corresponding to \eqn{X}{X} column
    names. \code{svc.cols} default argument
    of 1 results in a space-varying intercept model (assuming an intercept
    is specified in the first column of the design matrix).}
  
  \item{coords}{an \eqn{n \times m}{n x m} matrix of the observation coordinates
    in \eqn{R^m}{R^m} (e.g., \eqn{R^2}{R^2} might be easting and northing). }

  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{sigma.sq.ig}, \code{k.iw},
    \code{tau.sq.ig}, \code{phi.unif}, \code{nu.unif},
    \code{beta.norm}, and \code{beta.flat}. Scalar variance parameters \code{simga.sq} and
    \code{tau.sq} are assumed to follow an
    inverse-Gamma distribution. Cross-covariance matrix parameter \code{K} is
    assumed to follow an inverse-Wishart. The spatial decay \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform
    distributions. The regression coefficient priors can be either flat
    or multivariate Normal. 
    
    There are two specification for the Gaussian Process (GP) on the \code{svc.cols} columns: 1)
    univariate GPs on each column; 2) multivariate GP on the \eqn{r}{r}
    columns (i.e., where \eqn{r}{r} equals \code{length(svc.cols)}). If
    univariate GPs are desired, specify \code{sigma.sq.ig} as a
    list of length two with the first and second elements corresponding
    to the length \eqn{r}{r} \emph{shape} and \emph{scale}
    hyperparameter vectors,
    respectively. If a multivariate GP is desired, specify \code{k.iw}
    as a list of length two with the first and second elements corresponding
    to the degrees-of-freedom \eqn{df}{df} and \eqn{r\times r}{rxr} \emph{scale} matrix,
    respectively. This inverse-Wishart prior is on the \eqn{r\times r}{rxr} multivariate GP
    cross-covariance matrix defined as \eqn{K=AA'}{K=AA'} where \eqn{A}{A} is the
    lower-triangle Cholesky square root of \eqn{K}{K}.
    
    If the regression coefficients, i.e., \code{beta}
    vector, are assumed to follow a multivariate Normal distribution then pass the
    hyperparameters as a list of length two with the first
    and second elements corresponding to the mean vector and positive
    definite covariance matrix, respectively. If
    \code{beta} is assumed flat then no arguments are passed. The default
    is a flat prior. Similarly, \code{phi}
    and \code{nu} are specified as lists of length two with
    the first and second elements holding vectors of length \eqn{r}{r} lower
    and upper bounds of the Uniforms' support, respectively.}
  
  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{beta}, \code{sigma.sq}, \code{A},
    \code{tau.sq}, \code{phi}, and \code{nu}.  The value portion of each
    tag is the parameter's starting value(s). Starting values must be
    set for the \eqn{r}{r} univariate or multivariate GP \code{phi} and
    \code{nu}. For univariate GPs \code{sigma.sq.ig} is
    specified as a vector of length \eqn{r}{r} and for a multivariate GP \code{A} is specified as a vector of 
    \eqn{r\times(r+1)/2}{r*(r+1)/2} that gives the lower-triangle elements in column major ordering of the Cholesky square root
    of the cross-covaraince matrix \eqn{K=AA'}{K=AA'}. \code{tau.sq} is
    a single value. See Finley and Banerjee (2019) for more details.}
  
  \item{tuning}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{sigma.sq}, \code{A},
    \code{tau.sq}, \code{phi}, and \code{nu}. The value portion of each
    tag defines the variance of the Metropolis sampler Normal proposal
    distribution. For \code{sigma.sq}, \code{phi}, and \code{nu} the
    tuning value vectors are of length \eqn{r}{r} and \code{A} is of
    length \eqn{r\times(r+1)/2}{r*(r+1)/2}. \code{tuning} vector elements correspond to \code{starting} vector elements. \code{tau.sq} is
    a single value.}
  
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.}

  \item{center.scale}{if \code{TRUE}, non-constant columns of \eqn{X}{X}
    are centered on zero and scaled to have variance one. If
    \code{spPredict} is subsequently called this centering and scaling
    is applied to \code{pred.covars}. }
  
  \item{amcmc}{a list with tags \code{n.batch}, \code{batch.length}, and
    \code{accept.rate}. Specifying this argument invokes an adaptive MCMC
    sampler, see Roberts and Rosenthal (2007) for an explanation.}
  
  \item{n.samples}{the number of MCMC iterations. This argument is
    ignored if \code{amcmc} is specified.}
  
  \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we
    recommend setting \code{n.omp.threads} up to the number of
    hyperthreaded cores.}
  
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}
  
  \item{n.report}{the interval to report Metropolis sampler acceptance
    and MCMC progress.}  
  
  \item{...}{currently no additional arguments.}
}

\value{
  An object of class \code{spSVC}, which is a list comprising:
  \item{coords}{the \eqn{n \times m}{n x m} matrix specified by
    \code{coords}.  }
  \item{p.theta.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  \item{acceptance}{the Metropolis sampling
    acceptance percent. Reported at \code{batch.length} or \code{n.report}
  intervals for \code{amcmc} specified and non-specified, respectively.}

    The return object will include additional objects used for subsequent
  parameter recovery, prediction, and model fit evaluation using
  \code{\link{spRecover}}, \code{\link{spPredict}},
  \code{\link{spDiag}}, respectively. 
}

\details{
Model parameters can be fixed at their \code{starting} values by setting their
\code{tuning} values to zero.

The \emph{no nugget} model is specified by removing \code{tau.sq} from the \code{starting} list.
}
  
\seealso{
  \code{\link{spRecover}, \link{spDiag}, \link{spPredict}}
}

\references{
  Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software}, 63:1--28. \url{http://www.jstatsoft.org/v63/i13}.
  
  Roberts G.O. and Rosenthal J.S. (2006). Examples of Adaptive MCMC. \url{http://probability.ca/jeff/ftpdir/adaptex.pdf}.

  Finley, A.O. and S. Banerjee (2019) Bayesian spatially varying
  coefficient models in the spBayes R package. \url{http://arxiv.org}.
  
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudipto@ucla.edu}
}

\examples{
\dontrun{

library(Matrix)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}


##Assume both columns of X are space-varying and the two GPs don't covary
set.seed(1)
n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))

X <- as.matrix(cbind(1, rnorm(n)))
colnames(X) <- c("x.1", "x.2")

Z <- t(bdiag(as.list(as.data.frame(t(X)))))

B <- as.matrix(c(1,5))
p <- length(B)

sigma.sq <- c(1,5)
tau.sq <- 1
phi <- 3/0.5

D <- as.matrix(dist(coords))

C <- exp(-phi*D)\%x\%diag(sigma.sq)

w <- rmvn(1, rep(0,p*n), C)

mu <- as.vector(X\%*\%B + Z\%*\%w)

y <- rnorm(n, mu, sqrt(tau.sq))

##fit a model to the simulated dat
starting <- list("phi"=rep(3/0.5, p), "sigma.sq"=rep(1, p), "tau.sq"=1)

tuning <- list("phi"=rep(0.1, p), "sigma.sq"=rep(0.1, p), "tau.sq"=0.1)

cov.model <- "exponential"

priors <- list("phi.Unif"=list(rep(3/2, p), rep(3/0.0001, p)),
               "sigma.sq.IG"=list(rep(2, p), rep(2, p)),
               "tau.sq.IG"=c(2, 1))

##fit model
n.samples <- 2000

m.1 <- spSVC(y~X-1, coords=coords, starting=starting, svc.cols=c(1,2),
             tuning=tuning, priors=priors, cov.model=cov.model,
             n.samples=n.samples, n.omp.threads=4)

plot(m.1$p.theta.samples, density=FALSE)

##recover posterior samples
m.1 <- spRecover(m.1, start=floor(0.75*n.samples), thin=2, n.omp.threads=4)

summary(m.1$p.beta.recover.samples)
summary(m.1$p.theta.recover.samples)

##check fitted values
quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}

##fitted y
y.hat <- apply(m.1$p.y.samples, 1, quant)

rng <- c(-15, 20)
plot(y, y.hat[2,], pch=19, cex=0.5, xlab="Fitted y", ylab="Observed y",
     xlim=rng, ylim=rng)
arrows(y, y.hat[2,], y, y.hat[1,], angle=90, length=0.05)
arrows(y, y.hat[2,], y, y.hat[3,], angle=90, length=0.05)
lines(rng, rng, col="blue")

##recovered w
w.hat <- apply(m.1$p.w.recover.samples, 1, quant)

w.1.indx <- seq(1, p*n, p)
w.2.indx <- seq(2, p*n, p)

par(mfrow=c(1,2))

rng <- c(-5,5)
plot(w[w.1.indx], w.hat[2,w.1.indx], pch=19, cex=0.5, xlab="Fitted w.1", ylab="Observed w.1",
     xlim=rng, ylim=rng)
arrows(w[w.1.indx], w.hat[2,w.1.indx], w[w.1.indx], w.hat[1,w.1.indx], angle=90, length=0.05)
arrows(w[w.1.indx], w.hat[2,w.1.indx], w[w.1.indx], w.hat[3,w.1.indx], angle=90, length=0.05)
lines(rng, rng, col="blue")

rng <- c(-10,10)
plot(w[w.2.indx], w.hat[2,w.2.indx], pch=19, cex=0.5, xlab="Fitted w.2", ylab="Observed w.2",
     xlim=rng, ylim=rng)
arrows(w[w.2.indx], w.hat[2,w.2.indx], w[w.2.indx], w.hat[1,w.2.indx], angle=90, length=0.05)
arrows(w[w.2.indx], w.hat[2,w.2.indx], w[w.2.indx], w.hat[3,w.2.indx], angle=90, length=0.05)
lines(rng, rng, col="blue")
}
}

\keyword{model}