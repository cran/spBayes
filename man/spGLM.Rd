\name{spGLM}
\alias{spGLM}
\title{Function for fitting univariate Bayesian generalized linear spatial regression models}

\description{
  The function \code{spGLM} fits univariate Bayesian
  generalized linear spatial regression models. Given a set of knots,
  \code{spGLM} will also fit a \emph{predictive process} model (see references below).
}

\usage{
spGLM(formula, family="binomial", weights, data = parent.frame(),
      coords, knots, starting, tuning, priors, cov.model,
      amcmc, n.samples, verbose=TRUE,
      n.report=100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  \item{family}{currently only supports \code{binomial} and
  \code{poisson} data using the logit and log link functions,
  respectively.}
  \item{weights}{an optional vector of weights to be used in the fitting
  process. Weights correspond to number of trials and \emph{offset} for
  each location for the \code{binomial} and \code{poisson} family,
  respectively. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spGLM} is called.  }
 \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{knots}{either a \eqn{m \times 2}{m x 2} matrix of the
  \emph{predictive process} knot coordinates in \eqn{R^2}{R^2} (e.g., easting and northing)
    or a vector of length two or three with the first and second elements recording the
    number of columns and rows in the desired knot grid. The third,
    optional, element sets the offset of the outermost knots from the
    extent of the \code{coords}. }
  \item{starting}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{beta}, \code{sigma.sq},
    \code{phi}, \code{nu}, and \code{w}.  The value portion of each tag is the parameter's starting value. If the predictive
    process is used then \code{w} must be of length \eqn{m}{m}; otherwise,
    it must be of length \eqn{n}{n}. Alternatively, \code{w} can be set
    as a scalar, in which case the value is repeated.}
  
  \item{tuning}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{beta}, \code{sigma.sq},
    \code{phi}, \code{nu}, and \code{w}. The value portion of each tag
    defines the variance of the Metropolis sampler Normal proposal
    distribution. 

    The tuning value for \code{beta} can be a
    vector of length \eqn{p}{p} (where \eqn{p}{p} is the number of
    regression coefficients) or, if an adaptive MCMC is not used,
    i.e., \code{amcmc} is not specified, the lower-triangle of the
  \eqn{p\times p}{pxp} Cholesky square-root of the desired proposal
  covariance matrix. If the \emph{predictive process} is used then \code{w} must be of length \eqn{m}{m}; otherwise,
    it must be of length \eqn{n}{n}. Alternatively, \code{w} can be set
    as a scalar, in which case the value is repeated.  }
   \item{priors}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{sigma.sq.ig},
    \code{phi.unif}, \code{nu.unif},
  \code{beta.norm},  and \code{beta.flat}. Variance parameter
  \code{simga.sq} is assumed to follow an
    inverse-Gamma distribution, whereas the spatial decay \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform distributions. The
    hyperparameters of the inverse-Gamma are
    passed as a vector of length two, with the first and second elements corresponding
    to the \emph{shape} and \emph{scale}, respectively. The hyperparameters
    of the Uniform are also passed as a vector of length two with the first
    and second elements corresponding to the lower and upper support,
    respectively. If the regression coefficients are each assumed to follow a Normal distribution, i.e., \code{beta.norm}, then mean and variance hyperparameters are passed as the first and second list elements, respectively. If
  \code{beta} is assumed flat then no arguments are passed. The default
  is a flat prior.
 }  
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.  }
  \item{amcmc}{a list with tags \code{n.batch}, \code{batch.length}, and
  \code{accept.rate}. Specifying this argument invokes an adaptive MCMC
  sampler, see Roberts and Rosenthal (2007) for an explanation.}

\item{n.samples}{the number of MCMC iterations. This argument is
  ignored if \code{amcmc} is specified.}

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
  An object of class \code{spGLM}, which is a list with the following
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
    correspond to the \eqn{m}{m} knot locations and the columns are the
    posterior samples. This is only returned if a \emph{predictive process} model
    is used.}
  
    \item{p.w.samples}{a matrix that holds samples from the posterior
    distribution of the locations' spatial random effects. The rows of this matrix
    correspond to the \eqn{n}{n} point observations and the columns are the
    posterior samples. }
  
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\seealso{
  \code{\link{spMvGLM}}
}

\references{

  
  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian Predictive Process Models for Large Spatial Datasets. \emph{Journal of the Royal Statistical Society Series B}, 70:825--848.

  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004) Hierarchical modeling and analysis for spatial data. Chapman and Hall/CRC Press, Boca Raton, Fla.

    Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software}, 63:1--28. \url{http://www.jstatsoft.org/v63/i13}.

    Finley, A.O., S. Banerjee, and R.E. McRoberts. (2008) A Bayesian approach to quantifying uncertainty in multi-source forest area estimates. \emph{Environmental and Ecological Statistics}, 15:241--258.
  
  Roberts G.O. and Rosenthal J.S. (2006) Examples of Adaptive MCMC. \url{http://probability.ca/jeff/ftpdir/adaptex.pdf} Preprint.

}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
\dontrun{
library(MBA)

set.seed(1)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p) \%*\% D + rep(mu,rep(n,p)))
}

################################
##Spatial binomial
################################

##Generate binary data
coords <- as.matrix(expand.grid(seq(0,100,length.out=8), seq(0,100,length.out=8)))
n <- nrow(coords)

phi <- 3/50
sigma.sq <- 2

R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- rmvn(1, rep(0,n), R)

x <- as.matrix(rep(1,n))
beta <- 0.1
p <- 1/(1+exp(-(x\%*\%beta+w)))

weights <- rep(1, n)
weights[coords[,1]>mean(coords[,1])] <- 10

y <- rbinom(n, size=weights, prob=p)

##Collect samples
fit <- glm((y/weights)~x-1, weights=weights, family="binomial")
beta.starting <- coefficients(fit)
beta.tuning <- t(chol(vcov(fit)))

n.batch <- 200
batch.length <- 50
n.samples <- n.batch*batch.length

m.1 <- spGLM(y~1, family="binomial", coords=coords, weights=weights, 
             starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=beta.tuning, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Normal"=list(0,10), "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

p.hat <- 1/(1+exp(-(x\%*\%beta.hat+w.hat)))

y.hat <- apply(p.hat, 2, function(x){rbinom(n, size=weights, prob=p)})

y.hat.mu <- apply(y.hat, 1, mean)
y.hat.var <- apply(y.hat, 1, var)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Interpolated mean of posterior rate\n(observed rate)")
contour(surf, add=TRUE)
text(coords, label=paste("(",y,")",sep=""))

surf <- mba.surf(cbind(coords,y.hat.var),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Interpolated variance of posterior rate\n(observed #
of trials)")
contour(surf, add=TRUE)
text(coords, label=paste("(",weights,")",sep=""))

###########################
##Spatial poisson
###########################
##Generate count data
set.seed(1)

n <- 100

coords <- cbind(runif(n,1,100),runif(n,1,100))

phi <- 3/50
sigma.sq <- 2

R <- sigma.sq*exp(-phi*as.matrix(dist(coords)))
w <- rmvn(1, rep(0,n), R)

x <- as.matrix(rep(1,n))
beta <- 0.1
y <- rpois(n, exp(x\%*\%beta+w))

##Collect samples
beta.starting <- coefficients(glm(y~x-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~x-1, family="poisson"))))

n.batch <- 500
batch.length <- 50
n.samples <- n.batch*batch.length

##Note tuning list is now optional

m.1 <- spGLM(y~1, family="poisson", coords=coords,
             starting=list("beta"=beta.starting, "phi"=0.06,"sigma.sq"=1, "w"=0),
             tuning=list("beta"=0.1, "phi"=0.5, "sigma.sq"=0.5, "w"=0.5),
             priors=list("beta.Flat", "phi.Unif"=c(0.03, 0.3), "sigma.sq.IG"=c(2, 1)),
             amcmc=list("n.batch"=n.batch, "batch.length"=batch.length, "accept.rate"=0.43),
             cov.model="exponential", verbose=TRUE, n.report=10)

##Just for fun check out the progression of the acceptance
##as it moves to 43\% (same can be seen for the random spatial effects).
plot(mcmc(t(m.1$acceptance)), density=FALSE, smooth=FALSE)

##Now parameter summaries, etc.
burn.in <- 0.9*n.samples
sub.samps <- burn.in:n.samples

m.1$p.samples[,"phi"] <- 3/m.1$p.samples[,"phi"]

plot(m.1$p.beta.theta.samples)
print(summary(window(m.1$p.beta.theta.samples, start=burn.in)))

beta.hat <- m.1$p.beta.theta.samples[sub.samps,"(Intercept)"]
w.hat <- m.1$p.w.samples[,sub.samps]

y.hat <- apply(exp(x\%*\%beta.hat+w.hat), 2, function(x){rpois(n, x)})

y.hat.mu <- apply(y.hat, 1, mean)

##Take a look
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,y),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed counts")
contour(surf, add=TRUE)
text(coords, labels=y, cex=1)

surf <- mba.surf(cbind(coords,y.hat.mu),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Fitted counts")
contour(surf, add=TRUE)
text(coords, labels=round(y.hat.mu,0), cex=1)
}
}

\keyword{model}