\name{spMisalignLM}
\alias{spMisalignLM}
\title{Function for fitting multivariate Bayesian spatial regression
  models to misaligned data}

\description{
  The function \code{spMisalignLM} fits Gaussian multivariate Bayesian
  spatial regression models to misaligned data. 
}

\usage{
spMisalignLM(formula, data = parent.frame(), coords, 
      starting, tuning, priors, cov.model,
      amcmc, n.samples, verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a list of \eqn{q}{q} symbolic regression models to be fit. See example below. }
  
  \item{data}{an optional data frame containing the variables in the
    model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{spMisalignLM} is called.  }
  
  \item{coords}{a list of \eqn{q}{q} \eqn{n_i \times 2}{n_i x 2}
  matrices of the observation coordinates in \eqn{R^2}{R^2} (e.g.,
  easting and northing) where \eqn{i=(1,2,\ldots,q)}{i=(1,2,\ldots,q)}. . }
  
  \item{starting}{a list with tags corresponding to \code{A},
    \code{phi}, \code{nu}, and \code{Psi}. The value portion of each tag is a vector
    that holds the parameter's starting values.

    \code{A} is of length \eqn{\frac{q(q+1)}{2}}{q(q+1)/2} and holds the lower-triangle elements in column major ordering of the Cholesky square root
    of the spatial cross-covariance matrix.

    \code{phi} and \code{nu} are
    of length \eqn{q}{q}. The vector of residual variances \code{Psi} is also of length \eqn{q}{q}.}
  
  \item{tuning}{a list with tags \code{A}, \code{phi}, \code{nu}, and \code{Psi}. The value portion of each tag defines the
    variance of the Metropolis sampler Normal proposal distribution.
  \code{A} is of length \eqn{\frac{q(q+1)}{2}}{q(q+1)/2} and \code{Psi},
    \code{phi}, and \code{nu} are of length \eqn{q}{q}.}
    
  \item{priors}{a list with tags \code{beta.flat}, \code{K.iw}, \code{Psi.ig}, \code{phi.unif}
    and \code{nu.unif}. The hyperparameters of the inverse-Wishart for
    the cross-covariance matrix \eqn{K=AA'}{K=AA'} are
    passed as a list of length two, with the first and second elements corresponding
    to the \eqn{df}{df} and \eqn{q\times q}{qxq} \emph{scale} matrix,
    respectively.  The inverse-Gamma hyperparameters for the non-spatial
    residual variances are specified as a list \code{Psi.ig} of length two with the first and second list elements consisting of
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
  An object of class \code{spMisalignLM}, which is a list with the following
  tags:
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
}
  
\seealso{
  \code{\link{spMvLM}}\code{\link{spMisalignGLM}}
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

Psi <- c(0.1, 0.1, 0.1) ##non-spatial residual variance, i.e., nugget

y <- rnorm(n*q, x\%*\%B+w, rep(sqrt(Psi),n))

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

##call spMisalignLM
q <- 3
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
n.samples <- 5000

starting <- list("phi"=rep(3/0.5,q), "A"=A.starting, "Psi"=rep(1,q))
tuning <- list("phi"=rep(0.5,q), "A"=rep(0.01,length(A.starting)), "Psi"=rep(0.1,q))
priors <- list("phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
               "K.IW"=list(q+1, diag(0.1,q)), "Psi.ig"=list(rep(2,q), rep(0.1,q)))

m.1 <- spMisalignLM(list(y.1~x.1-1, y.2~x.2-1, y.3~x.3-1), 
                    coords=list(coords.1, coords.2, coords.3),
                    starting=starting, tuning=tuning, priors=priors, 
                    n.samples=n.samples, cov.model="exponential", n.report=100)

burn.in <- floor(0.75*n.samples)

plot(m.1$p.theta.samples, density=FALSE)

##recover regression coefficients and random effects
m.1 <- spRecover(m.1, start=burn.in)

round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)

##predict for all locations, i.e., observed and not observed
out <- spPredict(m.1, start=burn.in, thin=10, pred.covars=list(x.a, x.b,
x.c), 
                 pred.coords=list(coords, coords, coords))

##summary and check
quants <- function(x){quantile(x, prob=c(0.5,0.025,0.975))}

y.hat <- apply(out$p.y.predictive.samples, 1, quants)

##unstack and plot
y.a.hat <- y.hat[,1:n]
y.b.hat <- y.hat[,(n+1):(2*n)]
y.c.hat <- y.hat[,(2*n+1):(3*n)]

par(mfrow=c(1,3))
plot(y.a, y.a.hat[1,], xlab="Observed y.a", ylab="Fitted & predicted y.a",
     xlim=range(y), ylim=range(y.hat), main="")
arrows(y.a[-sub.1], y.a.hat[1,-sub.1], y.a[-sub.1], y.a.hat[2,-sub.1], length=0.02, angle=90)
arrows(y.a[-sub.1], y.a.hat[1,-sub.1], y.a[-sub.1], y.a.hat[3,-sub.1], length=0.02, angle=90)
lines(range(y.a), range(y.a))

plot(y.b, y.b.hat[1,], xlab="Observed y.b", ylab="Fitted & predicted y.b",
     xlim=range(y), ylim=range(y.hat), main="")
arrows(y.b[-sub.2], y.b.hat[1,-sub.2], y.b[-sub.2], y.b.hat[2,-sub.2], length=0.02, angle=90)
arrows(y.b[-sub.2], y.b.hat[1,-sub.2], y.b[-sub.2], y.b.hat[3,-sub.2], length=0.02, angle=90)
lines(range(y.b), range(y.b))

plot(y.c, y.c.hat[1,], xlab="Observed y.c", ylab="Fitted & predicted y.c",
     xlim=range(y), ylim=range(y.hat), main="")
arrows(y.c[-sub.3], y.c.hat[1,-sub.3], y.c[-sub.3], y.c.hat[2,-sub.3], length=0.02, angle=90)
arrows(y.c[-sub.3], y.c.hat[1,-sub.3], y.c[-sub.3], y.c.hat[3,-sub.3], length=0.02, angle=90)
lines(range(y.c), range(y.c))
}
}


\keyword{model}