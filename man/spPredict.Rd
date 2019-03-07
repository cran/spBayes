\name{spPredict}
\alias{spPredict}
\title{Function for new locations given a model object
}

\description{
  The function \code{spPredict} collects posterior predictive samples
  for a set of new locations given a \code{\link{spLM}}, \code{\link{spMvLM}},
  \code{\link{spGLM}}, \code{\link{spMvGLM}},
  \code{\link{spMisalignLM}}, \code{\link{spMisalignGLM}},
  \code{\link{bayesGeostatExact}}, \code{\link{bayesLMConjugate}} \code{\link{bayesLMRef}} or \code{\link{spSVC}} object.
}

\usage{
spPredict(sp.obj, pred.coords, pred.covars, joint=FALSE, start=1, end, thin=1,
          verbose=TRUE, n.report=100, n.omp.threads=1, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{spLM}}, \code{\link{spMvLM}},
    \code{\link{spGLM}}, \code{\link{spMvGLM}},
    \code{\link{spMisalignLM}}, \code{\link{spMisalignGLM}},
    \code{\link{bayesGeostatExact}}, \code{\link{bayesLMConjugate}} or
    \code{\link{bayesLMRef}}. For \code{\link{spSVC}}, \code{sp.obj} is an object from \code{\link{spRecover}}. }
  
  \item{pred.coords}{for \code{\link{spLM}}, \code{\link{spMvLM}},
    \code{\link{spGLM}}, \code{\link{spMvGLM}}, and \code{\link{bayesGeostatExact}} \code{pred.coords} is a \eqn{n^{\ast} \times 2}{n* x 2} matrix of \eqn{n^{\ast}}{n*} prediction
    location coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
    For \code{\link{spMisalignLM}} and
    \code{\link{spMisalignGLM}} \code{pred.coords} is a list of \eqn{q}{q} \eqn{n^{\ast}_i \times 2}{n*_i x 2}
    matrices of prediction location coordinates where
    \eqn{i=(1,2,\ldots,q)}{i=(1,2,\ldots,q)}. For \code{\link{spSVC}}
    \code{pred.coords} is an \eqn{n^{\ast} \times m}{n* x m} matrix of \eqn{n^{\ast}}{n*} prediction
    location coordinates in \eqn{R^m}{R^m}.}
  
  \item{pred.covars}{for \code{\link{spLM}}, \code{\link{spMvLM}},
    \code{\link{spGLM}}, \code{\link{spMvGLM}},
    \code{\link{bayesGeostatExact}}, \code{\link{bayesLMConjugate}},
    \code{\link{bayesLMRef}}, and \code{\link{spSVC}} \code{pred.covars} is a \eqn{n^{\ast} \times p}{n* x p} design matrix associated
    with the new locations (including the intercept if one is specified
    in \code{sp.obj}'s formula argument).  If this is a multivariate prediction defined
    by \eqn{q}{q} models, i.e., for \code{\link{spMvLM}} or \code{\link{spMvGLM}}, the multivariate design matrix can be created
    by passing a list of the \eqn{q}{q} univariate design matrices to
    the \code{\link{mkMvX}} function.  For \code{\link{spMisalignLM}} and
    \code{\link{spMisalignGLM}} \code{pred.covars} is a list of \eqn{q}{q} \eqn{n^{\ast}_i \times p_i}{n*_i x p_i}
    design matrices where \eqn{i=(1,2,\ldots,q)}{i=(1,2,\ldots,q)}}
  
  \item{joint}{specifies whether posterior samples should be drawn
    from the joint or point-wise predictive distribution. This argument is only
    implemented for \code{\link{spSVC}}. Prediction for all other models
    uses the point-wise posterior predictive distribution.}

  \item{start}{specifies the first sample included in the composition sampling.}

  \item{end}{specifies the last sample included in the composition.
    The default is to use all posterior samples in \code{sp.obj}. }
  
  \item{thin}{a sample thinning factor.  The default of 1 considers all
    samples between \code{start} and \code{end}.  For example, if \code{thin = 10} then 1 in 10 samples are considered between \code{start} and
    \code{end}. }

  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }

  \item{n.report}{the interval to report sampling progress. }  

  \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we
    recommend setting \code{n.omp.threads} up to the number of
    hyperthreaded cores. This argument is only
    implemented for \code{\link{spSVC}}.}
  
  \item{...}{currently no additional arguments.  }
}

\value{
  \item{p.y.predictive.samples}{a matrix that holds the response variable(s) posterior
    predictive samples. For multivariate models \code{\link{spMvLM}} or
    \code{\link{spMvGLM}} the rows of this matrix
    correspond to the predicted locations and the columns are the posterior
    predictive samples.  If prediction is for \eqn{q}{q} response
    variables the \code{p.y.predictive.samples} matrix has
    \eqn{qn^{\ast}}{qn*} rows, where \eqn{n^{\ast}}{n*} is the number of
  prediction locations.  The predictions for locations are held in rows
    \eqn{1:q, (q+1):2q, \ldots, ((n^{\ast}-1)q+1):qn^{\ast}}{1:q, (q+1):2q, \ldots,
      ((n*-1)q+1):qn*} (i.e., the samples for the first location's \eqn{q}{q}
    response variables are in rows \eqn{1:q}{1:q}, second location in rows \eqn{(q+1):2q}{(q+1):2q},
    etc.).

    For \code{\link{spMisalignLM}} and \code{\link{spMisalignGLM}}
    the posterior predictive samples are organized differently in
    \code{p.y.predictive.samples} with the first response variable
    \eqn{n^{\ast}_1}{n_*1} locations held in rows \eqn{1\ldots,n^{\ast}_1}{1,\ldots,n*_1} rows, then the
    next response variable samples held in the
    \eqn{(n^{\ast}_1+1),\ldots,(n^{\ast}_1+n^{\ast}_2)}{(n*_1+1),\ldots,(n*_1+n*_2)}, etc.
    
    For \code{\link{spSVC}} given the \eqn{r}{r} space-varying
    coefficients, \code{p.y.predictive.samples} has
    \eqn{rn^{\ast}}{rn*} rows and the columns are the posterior
    predictive samples. The predictions for coefficient are held in rows
    \eqn{1:r, (r+1):2r, \ldots, ((n^{\ast}-1)r+1):rn^{\ast}}{1:r, (r+1):2r, \ldots,
      ((n-1)r+1):rn*} (i.e., the samples for the first location's
    \eqn{r}{r} regression coefficients are in rows 1:r, second location in rows \eqn{(r+1):2r}{(r+1):2r},
    etc.).
    
    For \code{\link{spGLM}} and \code{\link{spMisalignGLM}} the \code{p.y.predictive.samples} matrix holds
    posterior predictive samples \eqn{\frac{1}{1+\exp(-x(s)'\beta-w(s))}}{1/1+(exp(-x(s)'B-w(s)))} and
    \eqn{\exp(x(s)'\beta+w(s))}{exp(x(s)'B+w(s))} for
    \code{family} binomial and poisson, respectively. Here \eqn{s}{s} indexes
    the prediction location, \eqn{\beta}{B} is the vector of regression
    coefficients, and \eqn{w}{w} is the associated spatial random
    spatial effect. These values can be fed directly into \code{rbinom}
    or \code{rpois} to generate the realization from the respective
    distribution. }
  
  \item{p.w.predictive.samples}{a matrix organized the same as
    \code{p.y.predictive.samples}, that holds the spatial random effects posterior
    predictive samples. }
  
  \item{p.w.predictive.samples.list}{only returned for
    \code{\link{spSVC}}. This provides \code{p.w.predictive.samples} in a
    different (more convenient form). Elements in this list hold
    samples for each of the \eqn{r}{r} coefficients. List element names
    indicate either the coefficient index or name specified in
    \code{\link{spSVC}}'s \code{svc.cols} argument. The sample matrices
    are \eqn{n^{\ast}}{n*} rows and predictive samples along the columns. }
    
  \item{p.tilde.beta.predictive.samples.list}{only returned for
    \code{\link{spSVC}}. Like \code{p.w.predictive.samples.list} but with
    the addition of the corresponding \eqn{\beta}{beta} posterior
    samples (i.e., \eqn{\beta+w(s)}{beta+w(s)}). }

  \item{center.scale.pred.covars}{only returned for the
    \code{\link{spSVC}} when its \code{center.scale} argument is \code{TRUE}. This is the
    prediction design matrix centered and scaled with respect to column means and variances of the design matrix used to estimate model
    parameters, i.e., the one defined in \code{\link{spSVC}}'s formula argument.}
}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical
  modeling and analysis for spatial data. Chapman and Hall/CRC Press,
  Boca Raton, FL.
  
  Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software},
  63:1--28. \url{http://www.jstatsoft.org/v63/i13}.

  Finley, A.O. and S. Banerjee (2019) Bayesian spatially varying
  coefficient models in the spBayes R package. \url{http://arxiv.org}.
      
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Sudipto Banerjee \email{sudipto@ucla.edu}
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

n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))
p <- length(B)
sigma.sq <- 10
tau.sq <- 0.01
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X\%*\%B + w, sqrt(tau.sq))

##partition the data for out of sample prediction
mod <- 1:100
y.mod <- y[mod]
X.mod <- X[mod,]
coords.mod <- coords[mod,]

n.samples <- 1000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
               "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 0.01))
cov.model <- "exponential"

m.1 <- spLM(y.mod~X.mod-1, coords=coords.mod, starting=starting, tuning=tuning,
priors=priors, cov.model=cov.model, n.samples=n.samples)

m.1.pred <- spPredict(m.1, pred.covars=X, pred.coords=coords,
start=0.5*n.samples)

y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, mean)

quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}

y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, quant)

plot(y, y.hat[2,], pch=19, cex=0.5, xlab="observed y", ylab="predicted y")
arrows(y[-mod], y.hat[2,-mod], y[-mod], y.hat[1,-mod], angle=90, length=0.05)
arrows(y[-mod], y.hat[2,-mod], y[-mod], y.hat[3,-mod], angle=90, length=0.05)
}
}

\keyword{utilities}