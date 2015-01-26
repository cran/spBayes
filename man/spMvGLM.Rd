\name{spMvGLM}
\alias{spMvGLM}
\title{Function for fitting multivariate Bayesian generalized linear spatial regression models}

\description{
  The function \code{spMvGLM} fits multivariate Bayesian
  generalized linear spatial regression models. Given a set of knots,
  \code{spMvGLM} will also fit a \emph{predictive process} model (see references below).
}

\usage{
spMvGLM(formula, family="binomial", weights, data = parent.frame(), coords, knots,
      starting, tuning, priors, cov.model,
      amcmc, n.samples, 
      verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a list of \eqn{q}{q} symbolic regression model descriptions to be fit. See example below. }
  \item{family}{currently only supports \code{binomial} and
    \code{poisson} data using the logit and log link functions,
    respectively.}
  \item{weights}{an optional \eqn{n \times q}{n x q} matrix of weights
  to be used in the fitting process. The order of the
  columns correspond to the univariate models in the formula list. Weights correspond to number of trials and \emph{offset} for
  each location for the \code{binomial} and \code{poisson} family,
  respectively.}
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spMvGLM} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{knots}{either a \eqn{m \times 2}{m x 2} matrix of the
  \emph{predictive process} knot coordinates in \eqn{R^2}{R^2} (e.g.,
  easting and northing) or a vector of length two or three with the
  first and second elements recording the  number of columns and rows in
  the desired knot grid. The third, optional, element sets the offset of
  the outermost knots from the extent of the \code{coords}. }

\item{starting}{a list with each tag corresponding to a parameter name. Valid tags are \code{beta}, \code{A}, \code{phi},
  \code{nu}, and \code{w}. The value portion of each tag is a vector
  that holds the parameter's starting values and are of length
  \eqn{p}{p} for \code{beta} (where \eqn{p}{p} is the total number of regression coefficients in the multivariate model), \eqn{\frac{q(q+1)}{2}}{q(q+1)/2} for
  \code{A}, and \eqn{q}{q} for \code{phi}, and \code{nu}. Here, \code{A} holds the the lower-triangle elements in column major ordering of the Cholesky square root
  of the spatial cross-covariance matrix. If the
    \emph{predictive process} is used then \code{w} 
    must be of length \eqn{qm}{qm}; otherwise,
    it must be of length \eqn{qn}{qn}. Alternatively, \code{w} can be set
    as a scalar, in which case the value is repeated. }

  \item{tuning}{a list with tags \code{beta}, \code{A}, \code{phi},
    \code{nu}, and \code{w}. The value portion of each tag defines the
    variance of the Metropolis sampler Normal proposal distribution. The value portion of these tags is of length \eqn{p}{p} for \code{beta}, \eqn{\frac{q(q+1)}{2}}{q(q+1)/2} for
  \code{A}, and \eqn{q}{q} for \code{phi}, and \code{nu}. Here, \code{A}
  holds the tuning values corresponding to the lower-triangle elements in column major ordering of the Cholesky square root
  of the spatial cross-covariance matrix.   If the
    \emph{predictive process} is used then \code{w} 
    must be of length \eqn{qm}{qm}; otherwise,
    it must be of length \eqn{qn}{qn}. Alternatively, \code{w} can be set
    as a scalar, in which case the value is repeated.  The tuning value for \code{beta} can be a vector of length \eqn{p}{p} or, if an adaptive MCMC is not used,
    i.e., \code{amcmc} is not specified, the lower-triangle of the
  \eqn{p\times p}{pxp} Cholesky square-root of the desired proposal
  covariance matrix.

  }

  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{beta.flat},
    \code{beta.norm}, \code{K.iw}, \code{phi.unif}, and
    \code{nu.unif}. If the regression coefficients are each assumed to follow a Normal distribution, i.e., \code{beta.norm}, then mean and variance hyperparameters are passed as the first and second list elements, respectively. If
  \code{beta} is assumed flat then no arguments are passed. The default
  is a flat prior. The spatial cross-covariance matrix \code{K} is assumed to follow an
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
  \item{n.report}{the interval to report Metropolis sampler acceptance and MCMC progress. }  
  \item{...}{currently no additional arguments.  }
}

\details{
  If a \code{binomial} model is specified the response vector is the
  number of successful trials at each location and \code{weights} is the
  total number of trials at each location.

  For a \code{poisson} specification, the \code{weights} vector is the
  count offset, e.g., population, at each location. This differs from
  the \code{\link{glm}} \code{offset} argument which is passed as the
  log of this value.  
  
  A non-spatial model is fit when \code{coords} is not specified. See
  example below.  
}

\value{
 An object of class \code{spMvGLM}, which is a list with the following
  tags:
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix specified by
    \code{coords}.  }
  \item{knot.coords}{the \eqn{m \times 2}{m x 2} matrix as specified by \code{knots}.  }
  \item{p.beta.theta.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  \item{acceptance}{the Metropolis sampler
    acceptance rate. If \code{amcmc} is used then this will be a matrix of
    each parameter's acceptance rate at the end of each
    batch. Otherwise, the sampler is a Metropolis with a joint proposal
    of all parameters.}
  
  \item{acceptance.w}{if this is a non-predictive process model and
  \code{amcmc} is used then this will be a matrix of the Metropolis sampler acceptance rate for each location's spatial random effect. }
  
  \item{acceptance.w.knots}{if this is a \emph{predictive process} model and \code{amcmc} is used then this will be a matrix of
    the Metropolis sampler acceptance rate for each knot's spatial random effect. }
  
  \item{p.w.knots.samples}{a matrix that holds samples from the posterior
    distribution of the knots' spatial random effects. The rows of this matrix
    correspond to the \eqn{q\times m}{q x m} knot locations and the columns are the
    posterior samples. This is only returned if a \emph{predictive process} model
    is used.}
  
    \item{p.w.samples}{a matrix that holds samples from the posterior
    distribution of the locations' spatial random effects. The rows of this matrix
    correspond to the \eqn{q\times n}{q x n} point observations and the columns are the
    posterior samples. }
  
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\seealso{
  \code{\link{spGLM}}
}

\references{
    Finley, A.O., S. Banerjee, and R.E. McRoberts. (2008) A Bayesian
    approach to quantifying uncertainty in multi-source forest area
    estimates. \emph{Environmental and Ecological Statistics},
15:241--258.

  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian
  Predictive Process Models for Large Spatial Datasets. \emph{Journal of
  the Royal Statistical Society Series B}, 70:825--848.

Finley, A.O., H. Sang, S. Banerjee, and A.E. Gelfand. (2009) Improving the performance of predictive process modeling for large datasets. \emph{Computational Statistics and Data Analysis}, 53:2873-2884.

    Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software}, 63:1--28. \url{http://www.jstatsoft.org/v63/i13}.

  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.
  
   Roberts G.O. and Rosenthal J.S. (2006) Examples of Adaptive MCMC. \url{http://probability.ca/jeff/ftpdir/adaptex.pdf} Preprint. 
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
\dontrun{
library(MBA)

##Some useful functions
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

weight <- 10 ##i.e., trials 
p <- 1/(1+exp(-(x\%*\%B+w)))
y <- rbinom(n*q, size=rep(weight,n*q), prob=p)

y.1 <- y[seq(1,length(y),q)]
y.2 <- y[seq(2,length(y),q)]

##Call spMvLM
fit <- glm((y/weight)~x-1, weights=rep(weight, n*q), family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]

n.batch <- 100
batch.length <- 50
n.samples <- n.batch*batch.length

starting <- list("beta"=beta.starting, "phi"=rep(3/0.5,q), "A"=A.starting, "w"=0)
tuning <- list("beta"=beta.tuning, "phi"=rep(1,q), "A"=rep(0.1,length(A.starting)),
               "w"=0.5)
priors <- list("beta.Flat", "phi.Unif"=list(rep(3/0.75,q), rep(3/0.25,q)),
               "K.IW"=list(q+1, diag(0.1,q)))

m.1 <- spMvGLM(list(y.1~x.1-1, y.2~x.2-1),
               coords=coords, weights=matrix(weight,n,q),
               starting=starting, tuning=tuning, priors=priors,
               amcmc=list("n.batch"=n.batch,"batch.length"=batch.length,"accept.rate"=0.43),
               cov.model="exponential", n.report=25)

burn.in <- 0.75*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in))$quantiles[,c(3,1,5)])

beta.hat <- t(m.1$p.beta.theta.samples[sub.samps,1:length(B)])
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x\%*\%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n*q, size=rep(weight, n*q), prob=p)})

y.hat.mu <- apply(y.hat, 1, mean)

##Unstack to get each response variable fitted values
y.hat.mu.1 <- y.hat.mu[seq(1,length(y.hat.mu),q)]
y.hat.mu.2 <- y.hat.mu[seq(2,length(y.hat.mu),q)]

##Take a look
par(mfrow=c(2,2))
surf <- mba.surf(cbind(coords,y.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed y.1 positive trials")
contour(surf, add=TRUE)
points(coords)
zlim <- range(surf[["z"]], na.rm=TRUE)

surf <- mba.surf(cbind(coords,y.hat.mu.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, zlim=zlim, main="Fitted y.1 positive trials")
contour(surf, add=TRUE)
points(coords)

surf <- mba.surf(cbind(coords,y.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed y.2 positive trials")
contour(surf, add=TRUE)
points(coords)
zlim <- range(surf[["z"]], na.rm=TRUE)

surf <- mba.surf(cbind(coords,y.hat.mu.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, zlim=zlim, main="Fitted y.2 positive trials")
contour(surf, add=TRUE)
points(coords)
}
}

\keyword{model}