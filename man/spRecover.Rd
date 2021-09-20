\name{spRecover}
\alias{spRecover}
\title{Function for recovering regression coefficients and spatial
  random effects for \code{\link{spLM}}, \code{\link{spMvLM}},
  \code{\link{spMisalignLM}}, \code{\link{spSVC}} using composition sampling
}

\description{
 Function for recovering regression coefficients and spatial random
 effects for \code{\link{spLM}}, \code{\link{spMvLM}}, and
  \code{\link{spMisalignLM}} using composition sampling.
}

\usage{
spRecover(sp.obj, get.beta=TRUE, get.w=TRUE, start=1, end, thin=1,
          verbose=TRUE, n.report=100, n.omp.threads=1, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{spLM}},
  \code{\link{spMvLM}}, \code{\link{spMisalignLM}}, or \code{\link{spSVC}}.}
  \item{get.beta}{if \code{TRUE}, regression coefficients will be recovered.}
  \item{get.w}{if \code{TRUE}, spatial random effects will be recovered.}
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
  The input \code{sp.obj} with posterior samples of regression coefficients and/or spatial random effects appended. 
  tags:
    \item{p.theta.recover.samples}{those \code{p.theta.samples} used in the composition sampling. }
    \item{p.beta.recover.samples}{a \code{coda} object of regression coefficients posterior samples. }
    \item{p.w.recover.samples}{a \code{coda} object of spatial random
      effects posterior samples. Rows correspond to locations' random
      effects and columns are posterior samples. Given \eqn{q}{q} responses, the \code{p.w.recover.samples} matrix
    for \code{\link{spMvLM}} has \eqn{qn}{qn} rows. The recovered
    random effects for locations are held in rows
    \eqn{1:q, (q+1):2q, \ldots, ((n-1)q+1):qn}{1:q, (q+1):2q, \ldots,
      ((n*-1)q+1):qn} (i.e., the samples for the first location's \eqn{q}{q}
    response variables are in rows 1:q, second location in rows \eqn{(q+1):2q}{(q+1):2q},
    etc.).

    For \code{\link{spSVC}} given the \eqn{r}{r} space-varying
    coefficients, \code{p.w.recover.samples} has
    \eqn{rn}{rn} rows. The recovered
    random effects for locations are held in rows
    \eqn{1:r, (r+1):2r, \ldots, ((n-1)r+1):rn}{1:r, (r+1):2r, \ldots,
      ((n-1)r+1):rn} (i.e., the samples for the first location's
    \eqn{r}{r} regression coefficients are in rows 1:r, second location in rows \eqn{(r+1):2r}{(r+1):2r},
    etc.).}

  \item{p.w.recover.samples.list}{only returned for
    \code{\link{spSVC}}. This provides \code{p.w.recover.samples} in a
    different (more convenient form). Elements in this list hold
    samples for each of the \eqn{r}{r} coefficients. List element names
    indicate either the coefficient index or name specified in
    \code{\link{spSVC}}'s \code{svc.cols} argument. The sample matrices
    are \eqn{n}{n} rows and predictive samples along the columns. }
    
  \item{p.tilde.beta.recover.samples.list}{only returned for
    \code{\link{spSVC}}. Like \code{p.w.recover.samples.list} but with
    the addition of the corresponding \eqn{\beta}{beta} posterior
    samples (i.e., \eqn{\beta+w(s)}{beta+w(s)}). }

    \item{p.y.samples}{only returned for
    \code{\link{spSVC}}. These posterior are the fitted values with locations on
    the rows and samples on the columns. For a given sample the fitted
    value for the \eqn{i^{th}}{i-th} location is \eqn{N(x(s_i)\beta +
    z(s_i)w(s_i), \tau^2)}{N(x(s_i)beta +
    z(s_i)w(s_i), tau^2)}. }
  
}

\references{
  Banerjee, S., Carlin, B.P., and Gelfand, A.E. (2004). Hierarchical
  modeling and analysis for spatial data. Chapman and Hall/CRC Press,
  Boca Raton, FL.

    Finley, A.O., S. Banerjee, and A.E. Gelfand. (2015) spBayes for large
  univariate and multivariate point-referenced spatio-temporal data
  models. \emph{Journal of Statistical Software},
  63:1--28. \url{https://www.jstatsoft.org/article/view/v063i13}.

  Finley, A.O. and S. Banerjee (2019) Bayesian spatially varying
  coefficient models in the spBayes R package. \url{https://arxiv.org/abs/1903.03028}.

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

n <- 50
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

n.samples <- 1000

starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
               "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 0.01))
cov.model <- "exponential"

m.1 <- spLM(y~X-1, coords=coords, starting=starting, tuning=tuning,
            priors=priors, cov.model=cov.model, n.samples=n.samples)

m.1 <- spRecover(m.1, start=0.5*n.samples, thin=2)

summary(window(m.1$p.beta.recover.samples))

w.hat <- apply(m.1$p.w.recover.samples, 1, mean)
plot(w, w.hat, xlab="Observed w", ylab="Fitted w")
}
}

\keyword{utilities}