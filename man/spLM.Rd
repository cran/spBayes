\name{spLM}
\alias{spLM}
\title{Function for fitting univariate Bayesian spatial regression models}

\description{
  The function \code{spLM} fits Gaussian univariate Bayesian
  spatial regression models. Given a set of knots, \code{spLM} will also
  fit a \emph{predictive process} model (see references below).
}

\usage{
spLM(formula, data = parent.frame(), coords, knots,
      starting, tuning, priors, cov.model,
      modified.pp = TRUE, amcmc, n.samples, 
      verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spLM} is called.  }
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
    \code{tau.sq}, \code{phi}, and \code{nu}.  The value portion of each tag is the parameter's starting value. }

  \item{tuning}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{sigma.sq},
    \code{tau.sq}, \code{phi}, and \code{nu}. The value portion of each tag defines the variance of the Metropolis sampler Normal proposal distribution.}

  \item{priors}{a list with each tag corresponding to a
    parameter name. Valid tags are \code{sigma.sq.ig},
    \code{tau.sq.ig}, \code{phi.unif}, \code{nu.unif},
  \code{beta.norm}, and \code{beta.flat}. Variance parameters, \code{simga.sq} and
    \code{tau.sq}, are assumed to follow an
    inverse-Gamma distribution, whereas the spatial decay \code{phi}
    and smoothness \code{nu} parameters are assumed to follow Uniform distributions. The
    hyperparameters of the inverse-Gamma are
    passed as a vector of length two, with the first and second elements corresponding
    to the \emph{shape} and \emph{scale}, respectively. The hyperparameters
    of the Uniform are also passed as a vector of length two with the first
    and second elements corresponding to the lower and upper support,
    respectively. If the regression coefficients, i.e., \code{beta}
    vector, are assumed to follow a multivariate Normal distribution then pass the
  hyperparameters as a list of length two with the first
  and second elements corresponding to the mean vector and positive
  definite covariance matrix, respectively. If
  \code{beta} is assumed flat then no arguments are passed. The default
  is a flat prior.
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
  sampler, see Roberts and Rosenthal (2007) for an explanation.}
  \item{n.samples}{the number of MCMC iterations. This argument is
  ignored if \code{amcmc} is specified.}
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report Metropolis sampler acceptance and MCMC progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{spLM}, which is a list with the following
  tags:
  \item{coords}{the \eqn{n \times 2}{n x 2} matrix specified by
    \code{coords}.  }
  \item{knot.coords}{the \eqn{m \times 2}{m x 2} matrix as specified by \code{knots}.  }
  \item{p.theta.samples}{a \code{coda} object of posterior samples for the defined
    parameters. }
  \item{acceptance}{the Metropolis sampling
    acceptance percent. Reported at \code{batch.length} or \code{n.report}
  intervals for \code{amcmc} specified and non-specified, respectively.}
  The return object might include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\details{
Model parameters can be fixed at their \code{starting} values by setting their
\code{tuning} values to zero.

The \emph{no nugget} model is specified by removing \code{tau.sq} from the \code{starting} list.
}
  
\seealso{
  \code{\link{spMvLM}}
}

\references{
  Banerjee, S., A.E. Gelfand, A.O. Finley, and H. Sang. (2008) Gaussian Predictive Process Models for Large Spatial Datasets. \emph{Journal of the Royal Statistical Society Series B}, 70:825--848.

  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical
  modeling and analysis for spatial data. Chapman and Hall/CRC Press,
  Boca Raton, FL.

  Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software}, 63:1--28. \url{http://www.jstatsoft.org/v63/i13}.
  
  Finley, A.O., H. Sang, S. Banerjee, and A.E. Gelfand. (2009) Improving the performance of predictive process modeling for large datasets. \emph{Computational Statistics and Data Analysis}, 53:2873--2884.
  
  Roberts G.O. and Rosenthal J.S. (2006). Examples of Adaptive MCMC. \url{http://probability.ca/jeff/ftpdir/adaptex.pdf}.

}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{baner009@umn.edu}
}

\examples{
\dontrun{
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)\%*\%D + rep(mu,rep(n,p)))
}

set.seed(1)

n <- 100
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))
p <- length(B)

sigma.sq <- 2
tau.sq <- 0.1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X\%*\%B + w, sqrt(tau.sq))

n.samples <- 2000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                 "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))

priors.2 <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
                 "sigma.sq.IG"=c(2, 2), "tau.sq.IG"=c(2, 0.1))

cov.model <- "exponential"

n.report <- 500
verbose <- TRUE

m.1 <- spLM(y~X-1, coords=coords, starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

m.2 <- spLM(y~X-1, coords=coords, starting=starting,
            tuning=tuning, priors=priors.2, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

burn.in <- 0.5*n.samples

##recover beta and spatial random effects
m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)
m.2 <- spRecover(m.2, start=burn.in, verbose=FALSE)

round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.2$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.2$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)

m.1.w.summary <- summary(mcmc(t(m.1$p.w.recover.samples)))$quantiles[,c(3,1,5)]
m.2.w.summary <- summary(mcmc(t(m.2$p.w.recover.samples)))$quantiles[,c(3,1,5)]

plot(w, m.1.w.summary[,1], xlab="Observed w", ylab="Fitted w",
     xlim=range(w), ylim=range(m.1.w.summary), main="Spatial random effects")
arrows(w, m.1.w.summary[,1], w, m.1.w.summary[,2], length=0.02, angle=90)
arrows(w, m.1.w.summary[,1], w, m.1.w.summary[,3], length=0.02, angle=90)
lines(range(w), range(w))

points(w, m.2.w.summary[,1], col="blue", pch=19, cex=0.5)
arrows(w, m.2.w.summary[,1], w, col="blue", m.2.w.summary[,2], length=0.02, angle=90)
arrows(w, m.2.w.summary[,1], w, col="blue", m.2.w.summary[,3], length=0.02, angle=90)

###########################
##Predictive process model
###########################
m.1 <- spLM(y~X-1, coords=coords, knots=c(6,6,0.1), starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

m.2 <- spLM(y~X-1, coords=coords, knots=c(6,6,0.1), starting=starting,
            tuning=tuning, priors=priors.2, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

burn.in <- 0.5*n.samples

round(summary(window(m.1$p.beta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.2$p.beta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)

round(summary(window(m.1$p.theta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
round(summary(window(m.2$p.theta.samples, start=burn.in))$quantiles[,c(3,1,5)],2)
}
}

\keyword{model}