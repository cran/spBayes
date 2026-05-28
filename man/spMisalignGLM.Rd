\name{spMisalignGLM}
\alias{spMisalignGLM}
\title{Function for fitting multivariate generalized linear Bayesian spatial regression models to misaligned data}

\description{
  The function \code{spMisalignGLM} fits Gaussian multivariate Bayesian
  generalized linear spatial regression models to misaligned data. 
}

\usage{
spMisalignGLM(formula, family="binomial", weights, data = parent.frame(), coords, 
      starting, tuning, priors, cov.model,
      amcmc, n.samples, verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a list of \eqn{q}{q} symbolic regression models to be fit. See example below. }

  \item{family}{currently only supports \code{binomial} and
    \code{poisson} data using the logit and log link functions,
    respectively.}
  
  \item{weights}{an optional list of weight vectors associated with each model
  in the formula list. Weights correspond to number of trials and \emph{offset} for
  each location for the \code{binomial} and \code{poisson} family,
  respectively.}
  
  \item{data}{an optional data frame containing the variables in the
    model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{spMisalignGLM} is called.  }
  
  \item{coords}{a list of \eqn{q}{q} \eqn{n_i \times 2}{n_i x 2}
  matrices of the observation coordinates in \eqn{R^2}{R^2} (e.g.,
  easting and northing) where \eqn{i=(1,2,\ldots,q)}{i=(1,2,\ldots,q)}. }
  
  \item{starting}{a list with tags corresponding to \code{A},
    \code{phi}, and \code{nu}. The value portion of each tag is a vector
    that holds the parameter's starting values. \code{A} is of length \eqn{\frac{q(q+1)}{2}}{q(q+1)/2} and holds the lower-triangle elements in column major ordering of the Cholesky square root
    of the spatial cross-covariance matrix \eqn{K=AA'}{K=AA'}. \code{phi} and \code{nu} are
    of length \eqn{q}{q}.}
  
  \item{tuning}{a list with tags \code{A}, \code{phi}, and \code{nu}. The value portion of each tag defines the
    variance of the Metropolis sampler Normal proposal distribution.
   \code{A} is of length \eqn{\frac{q(q+1)}{2}}{q(q+1)/2} and \code{phi}
    and \code{nu} are of length \eqn{q}{q}.}
  
  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{beta.flat},
    \code{beta.norm}, \code{K.iw}, \code{phi.unif}, and
    \code{nu.unif}. If the regression coefficients are each assumed to follow a Normal distribution, i.e., \code{beta.norm}, then mean and variance hyperparameters are passed as the first and second list elements, respectively. If
  \code{beta} is assumed flat then no arguments are passed. The default
  is a flat prior. The spatial cross-covariance matrix \eqn{K=AA'}{K=AA'} is assumed to follow an
    inverse-Wishart distribution, whereas the spatial decay \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform distributions. The
    hyperparameters of the inverse-Wishart are
    passed as a list of length two, with the first and second elements corresponding
    to the \eqn{df}{df} and \eqn{q\times q}{qxq} \emph{scale} matrix,
    respectively. The hyperparameters of the Uniform are also passed as a list of vectors with the first and second list elements corresponding to the lower and upper
    support, respectively. }

  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }

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
 An object of class \code{spMisalignGLM}, which is a list with the following
  tags:

  \item{p.beta.theta.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  
  \item{acceptance}{the Metropolis sampler
    acceptance rate. If \code{amcmc} is used then this will be a matrix of
    each parameter's acceptance rate at the end of each
    batch. Otherwise, the sampler is a Metropolis with a joint proposal
    of all parameters.}
  
  \item{acceptance.w}{if \code{amcmc} is used then this will be a matrix of the Metropolis sampler acceptance rate for each location's spatial random effect. }
  
    \item{p.w.samples}{a matrix that holds samples from the posterior
    distribution of the locations' spatial random effects. Posterior samples are organized with the first response variable
  \eqn{n_1}{n_1} locations held in rows \eqn{1\ldots,n_1}{1,\ldots,n_1} rows, then the
  next response variable samples held in the
  \eqn{(n_1+1),\ldots,(n_1+n_2)}{(n_1+1),\ldots,(n_1+n_2)}, etc.}
  
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\details{
  If a \code{binomial} model is specified the response vector is the
  number of successful trials at each location and \code{weights} is the
  total number of trials at each location.

  For a \code{poisson} specification, the \code{weights} vector is the
  count offset, e.g., population, at each location. This differs from
  the \code{\link{glm}} \code{offset} argument which is passed as the
  log of this value.  
}
  
\seealso{
  \code{\link{spMvGLM}} \code{\link{spMisalignLM}}
}

\references{
  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian Predictive Process Models for Large Spatial Datasets. \emph{Journal of the Royal Statistical Society Series B}, 70:825--848.
  
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.

Finley, A.O., S. Banerjee, and B.D. Cook. (2014) Bayesian hierarchical models for spatially misaligned data. \emph{Methods in Ecology and Evolution}, 5:514--523.
  
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

##generate some data
n <- 100 ##number of locations
q <- 3 ##number of outcomes at each location
nltr <- q*(q+1)/2 ##number of triangular elements in the cross-covariance matrix

coords <- cbind(runif(n,0,1), runif(n,0,1))

##parameters for generating a multivariate spatial GP covariance matrix
theta <- rep(3/0.5,q) ##spatial decay

A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- c(1,1,-1,1,0.5,0.25)
K <- A\%*\%t(A)
K ##spatial cross-covariance
cov2cor(K) ##spatial cross-correlation

C <- mkSpCov(coords, K, diag(0,q), theta, cov.model="exponential")

w <- rmvn(1, rep(0,nrow(C)), C) ##spatial random effects

w.a <- w[seq(1,length(w),q)]
w.b <- w[seq(2,length(w),q)]
w.c <- w[seq(3,length(w),q)]

##covariate portion of the mean
x.a <- cbind(1, rnorm(n))
x.b <- cbind(1, rnorm(n))
x.c <- cbind(1, rnorm(n))
x <- mkMvX(list(x.a, x.b, x.c))

B.1 <- c(1,-1)
B.2 <- c(-1,1)
B.3 <- c(-1,-1)
B <- c(B.1, B.2, B.3)

y <- rpois(nrow(C), exp(x\%*\%B+w))

y.a <- y[seq(1,length(y),q)]
y.b <- y[seq(2,length(y),q)]
y.c <- y[seq(3,length(y),q)]

##subsample to make spatially misaligned data
sub.1 <- 1:50
y.1 <- y.a[sub.1]
w.1 <- w.a[sub.1]
coords.1 <- coords[sub.1,]
x.1 <- x.a[sub.1,]

sub.2 <- 25:75
y.2 <- y.b[sub.2]
w.2 <- w.b[sub.2]
coords.2 <- coords[sub.2,]
x.2 <- x.b[sub.2,]

sub.3 <- 50:100
y.3 <- y.c[sub.3]
w.3 <- w.c[sub.3]
coords.3 <- coords[sub.3,]
x.3 <- x.c[sub.3,]

##call spMisalignGLM
q <- 3
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

n.batch <- 200
batch.length <- 25
n.samples <- n.batch*batch.length

starting <- list("beta"=rep(0,length(B)), "phi"=rep(3/0.5,q), "A"=A.starting, "w"=0)
                 
tuning <- list("beta"=rep(0.1,length(B)), "phi"=rep(1,q), "A"=rep(0.1,length(A.starting)), "w"=1)

priors <- list("phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
               "K.IW"=list(q+1, diag(0.1,q)),  rep(0.1,q))

m.1 <- spMisalignGLM(list(y.1~x.1-1, y.2~x.2-1, y.3~x.3-1), family="poisson",
                     coords=list(coords.1, coords.2, coords.3),
                     starting=starting, tuning=tuning, priors=priors,
                     amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
                     cov.model="exponential", n.report=10)

burn.in <- floor(0.75*n.samples)

plot(m.1$p.beta.theta.samples, density=FALSE)

##predict for all locations, i.e., observed and not observed
out <- spPredict(m.1, start=burn.in, thin=10, pred.covars=list(x.a, x.b, x.c), 
                 pred.coords=list(coords, coords, coords))

##summary and check
quants <- function(x){quantile(x, prob=c(0.5,0.025,0.975))}

y.hat <- apply(out$p.y.predictive.samples, 1, quants)

##unstack and plot
y.a.hat <- y.hat[,1:n]
y.b.hat <- y.hat[,(n+1):(2*n)]
y.c.hat <- y.hat[,(2*n+1):(3*n)]

par(mfrow=c(1,3))
plot(y.a ,y.a.hat[1,], xlab="Observed y.a", ylab="Fitted & predicted y.a")
plot(y.b, y.b.hat[1,], xlab="Observed y.b", ylab="Fitted & predicted y.b")
plot(y.c, y.c.hat[1,], xlab="Observed y.c", ylab="Fitted & predicted y.c")

  }
}

\keyword{model}